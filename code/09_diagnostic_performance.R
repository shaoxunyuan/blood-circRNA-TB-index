#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Purpose:      Calculate AUROC, DeLong 95% CI, AUPRC, and operating
#               characteristics at the two locked IDC-derived TB-index
#               thresholds (sensitivity, specificity, PPV, NPV, TP, FN, TN, FP).
# Inputs:       results/locked_TB_index_scores.tsv   (output of script 08)
#               config/locked_thresholds.tsv
# Outputs:      results/diagnostic_performance.tsv   (threshold-free metrics)
#               results/diagnostic_thresholds.tsv    (operating characteristics)
# Key parameters:
#               estimator     pROC::roc, DeLong 95% CI
#               case/control  case = 1, control = 0; direction = "<" so that a
#                             higher TB-index calls a case
#               thresholds    read from config/locked_thresholds.tsv
#                             (Youden = 1.407, high-sensitivity = 0.261)
#               positive call TB-index >= threshold
#               AUPRC         trapezoidal integration of the precision-recall
#                             curve after ranking samples by decreasing TB-index
# Dependencies: R >= 4.4.3; data.table; pROC
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({library(data.table); library(pROC)})
source("code/_paths.R")

TOL <- 1e-3

## ---- 1. Thresholds (single source of truth) --------------------------------
thr <- fread(file.path(CONFIG, "locked_thresholds.tsv"))
youden_cut <- thr[Threshold_name == "Youden", TB_index_cutoff]
hs_cut     <- thr[Threshold_name == "High_sensitivity", TB_index_cutoff]
stopifnot(length(youden_cut) == 1, length(hs_cut) == 1)

## ---- 2. Scores -------------------------------------------------------------
s_f <- file.path(RESULTS, "locked_TB_index_scores.tsv")
paths_require(s_f, "locked TB-index scores (output of script 08)",
              "run code/08_TB_index_scoring.R, or NCBI GEO GSE348100")
s <- fread(s_f)
stopifnot(sum(is.na(s$TB_index)) == 0)

## ---- 3. AUPRC: trapezoidal integration of the PR curve ---------------------
## AUPRC: trapezoidal integration of the precision-recall curve after ranking
## samples by decreasing TB-index (manuscript Methods). This is the original
## implementation used for Table 2: the (recall = 0, precision = 1) anchor is
## prepended, points at which recall does not increase are dropped, and the
## remaining points are integrated by the trapezoidal rule.
auprc <- function(labels, scores) {
  ord <- order(scores, decreasing = TRUE)
  y <- labels[ord]
  n_pos <- sum(y); n <- length(y)
  if (n_pos == 0 || n_pos == n) return(NA_real_)
  tp <- cumsum(y); fp <- cumsum(1 - y)
  prec <- tp / (tp + fp); rec <- tp / n_pos
  prec <- c(1, prec); rec <- c(0, rec)
  k <- c(TRUE, diff(rec) > 0)
  prec <- prec[k]; rec <- rec[k]
  sum(diff(rec) * (prec[-1] + prec[-length(prec)]) / 2)
}

## ---- 4. Operating characteristics at a threshold ---------------------------
op_at <- function(y, score, cut, name) {
  pred <- as.integer(score >= cut)             # positive call: TB-index >= threshold
  TP <- sum(pred == 1 & y == 1); FP <- sum(pred == 1 & y == 0)
  FN <- sum(pred == 0 & y == 1); TN <- sum(pred == 0 & y == 0)
  den_sen <- TP + FN; den_spe <- TN + FP
  data.table(
    Threshold_name = name, Threshold = cut,
    TP = TP, FN = FN, TN = TN, FP = FP,
    Sensitivity = TP / den_sen, Specificity = TN / den_spe,
    PPV = if (TP + FP > 0) TP / (TP + FP) else NA_real_,
    NPV = if (TN + FN > 0) TN / (TN + FN) else NA_real_,
    Sens_num = TP, Sens_den = den_sen, Spec_num = TN, Spec_den = den_spe
  )
}

## ---- 5. Comparisons --------------------------------------------------------
compute <- function(coh, case, ctrl, baseline_only = FALSE) {
  d <- s[Cohort == coh & Group %in% c(case, ctrl)]
  if (baseline_only) {
    stopifnot("Baseline_diagnostic_analysis" %in% names(d))
    d <- d[Baseline_diagnostic_analysis == "Yes"]
  }
  y <- as.integer(d$Group == case)             # case = 1, control = 0
  score <- d$TB_index
  r <- roc(response = y, predictor = score, levels = c(0, 1), direction = "<", quiet = TRUE)
  ci <- as.numeric(ci.auc(r, method = "delong"))
  tf <- data.table(
    Cohort = coh, Comparison = paste(case, "vs", ctrl),
    Cases_n = sum(y == 1), Controls_n = sum(y == 0),
    AUROC = as.numeric(auc(r)),
    AUROC_CI_lower = ci[1], AUROC_CI_upper = ci[3],
    AUPRC = auprc(y, score))
  ops <- rbindlist(list(op_at(y, score, youden_cut, "Youden"),
                        op_at(y, score, hs_cut, "High_sensitivity")))
  list(tf = tf, ops = cbind(tf, ops))
}

## Included comparisons: IDC (four contrasts + pooled non-APTB), ITC (HC, LC),
## PVC (PTB/EPTB vs HC). PVC PTB/EPTB versus incident LTBI are added below.
##
## PVC diagnostic comparisons use BASELINE libraries only
## (Included_in_baseline_diagnostic_analysis == "Yes"), i.e. PTB n = 8 and
## EPTB n = 8 against HC n = 14, as reported in Table 2. Later PVC time points
## are reserved for the longitudinal analysis (script 11).
BASE <- list(c("IDC", "APTB", "HC"), c("IDC", "APTB", "LTBI"), c("IDC", "APTB", "PNEU"),
             c("IDC", "APTB", "LC"), c("ITC", "APTB", "HC"), c("ITC", "APTB", "LC"),
             c("PVC", "PTB", "HC"), c("PVC", "EPTB", "HC"))

res <- lapply(BASE, function(k) compute(k[1], k[2], k[3], baseline_only = (k[1] == "PVC")))
tf  <- rbindlist(lapply(res, `[[`, "tf"))
ops <- rbindlist(lapply(res, `[[`, "ops"))

## IDC APTB versus pooled non-APTB
d <- s[Cohort == "IDC"]; y <- as.integer(d$Group == "APTB"); score <- d$TB_index
r <- roc(y, score, quiet = TRUE, direction = "<"); ci <- as.numeric(ci.auc(r, method = "delong"))
tf <- rbind(tf, data.table(Cohort = "IDC", Comparison = "APTB vs pooled non-APTB",
                           Cases_n = sum(y == 1), Controls_n = sum(y == 0),
                           AUROC = as.numeric(auc(r)), AUROC_CI_lower = ci[1],
                           AUROC_CI_upper = ci[3], AUPRC = auprc(y, score)), fill = TRUE)
ops <- rbind(ops, cbind(tf[Comparison == "APTB vs pooled non-APTB"],
                        rbindlist(list(op_at(y, score, youden_cut, "Youden"),
                                       op_at(y, score, hs_cut, "High_sensitivity")))),
             fill = TRUE)

## PVC baseline cases versus incident (post-conversion) LTBI
pvc <- s[Cohort == "PVC"]; ltbi <- pvc[Group == "LTBI"]
for (g in c("PTB", "EPTB")) {
  cases <- pvc[Group == g & TimePoint == "BL"]
  if (nrow(cases) < 2 || nrow(ltbi) < 2) next
  d <- rbind(cases, ltbi); y <- as.integer(d$Group == g); score <- d$TB_index
  r <- roc(y, score, quiet = TRUE, direction = "<"); ci <- as.numeric(ci.auc(r, method = "delong"))
  row_tf <- data.table(Cohort = "PVC", Comparison = paste(g, "vs incident LTBI"),
                       Cases_n = nrow(cases), Controls_n = nrow(ltbi),
                       AUROC = as.numeric(auc(r)), AUROC_CI_lower = ci[1],
                       AUROC_CI_upper = ci[3], AUPRC = auprc(y, score))
  tf <- rbind(tf, row_tf, fill = TRUE)
  ops <- rbind(ops, cbind(row_tf, rbindlist(list(op_at(y, score, youden_cut, "Youden"),
                                                 op_at(y, score, hs_cut, "High_sensitivity")))),
               fill = TRUE)
}

fwrite(tf,  file.path(RESULTS, "diagnostic_performance.tsv"), sep = "\t")
fwrite(ops, file.path(RESULTS, "diagnostic_thresholds.tsv"),  sep = "\t")

## ---- 6. QA against the manuscript (Supplementary Tables S3/S6, Table 2) ----
EXPECTED <- data.table(
  Comparison = c("APTB vs HC", "APTB vs LTBI", "APTB vs PNEU", "APTB vs LC",
                 "APTB vs pooled non-APTB", "APTB vs HC", "APTB vs LC",
                 "PTB vs HC", "EPTB vs HC", "PTB vs incident LTBI", "EPTB vs incident LTBI"),
  Cohort = c(rep("IDC", 5), rep("ITC", 2), rep("PVC", 4)),
  AUROC = c(0.917, 0.880, 0.969, 0.989, 0.949, 0.933, 0.885, 0.795, 0.357, 0.825, 0.312),
  AUPRC = c(0.937, 0.914, 0.988, 0.992, NA,      0.940, 0.878, 0.778, 0.377, NA,    NA),
  Youden_sens = c(0.783, 0.783, 0.783, 0.783, NA, 0.083, 0.083, 0.125, 0.000, NA, NA),
  Youden_spec = c(0.958, 0.846, 1.000, 1.000, NA, 1.000, 1.000, 0.929, 0.929, NA, NA),
  HiSens_sens = c(0.870, 0.870, 0.870, 0.870, NA, 0.333, 0.333, 0.625, 0.000, NA, NA),
  HiSens_spec = c(0.792, 0.846, 0.929, 1.000, NA, 1.000, 1.000, 0.929, 0.929, NA, NA)
)
fail <- character(0)
for (i in seq_len(nrow(EXPECTED))) {
  e <- EXPECTED[i]
  cur <- tf[Cohort == e$Cohort & Comparison == e$Comparison]
  if (nrow(cur) != 1) { fail <- c(fail, paste(e$Cohort, e$Comparison, "(missing)")); next }
  if (abs(cur$AUROC - e$AUROC) > TOL) fail <- c(fail, sprintf("%s %s AUROC %.4f vs %.4f", e$Cohort, e$Comparison, cur$AUROC, e$AUROC))
  if (!is.na(e$AUPRC) && abs(cur$AUPRC - e$AUPRC) > TOL) fail <- c(fail, sprintf("%s %s AUPRC %.4f vs %.4f", e$Cohort, e$Comparison, cur$AUPRC, e$AUPRC))
  oy <- ops[Cohort == e$Cohort & Comparison == e$Comparison & Threshold_name == "Youden"]
  oh <- ops[Cohort == e$Cohort & Comparison == e$Comparison & Threshold_name == "High_sensitivity"]
  if (!is.na(e$Youden_sens) && abs(oy$Sensitivity - e$Youden_sens) > TOL) fail <- c(fail, sprintf("%s %s Youden sens", e$Cohort, e$Comparison))
  if (!is.na(e$Youden_spec) && abs(oy$Specificity - e$Youden_spec) > TOL) fail <- c(fail, sprintf("%s %s Youden spec", e$Cohort, e$Comparison))
  if (!is.na(e$HiSens_sens) && abs(oh$Sensitivity - e$HiSens_sens) > TOL) fail <- c(fail, sprintf("%s %s HiSens sens", e$Cohort, e$Comparison))
  if (!is.na(e$HiSens_spec) && abs(oh$Specificity - e$HiSens_spec) > TOL) fail <- c(fail, sprintf("%s %s HiSens spec", e$Cohort, e$Comparison))
}
if (length(fail)) {
  stop("Table 2 / Supplementary Table QA failed for:\n  ", paste(fail, collapse = "\n  "),
       call. = FALSE)
}
message("QA: ", nrow(tf), " comparisons and ", nrow(ops),
        " threshold rows reproduce Table 2 / Supplementary Tables S3 and S6")
