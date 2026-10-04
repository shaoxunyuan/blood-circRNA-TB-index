#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Purpose:      Paired DeLong comparisons between the circRNA TB-index and the
#               Sweeney3 / RISK6 benchmark signatures.
# Inputs:       data/benchmark/SF6_Sweeney3_RISK6_per_sample.csv
#                 per-sample TB-index and benchmark signature scores
#                 (Supplementary Table S7 / NCBI GEO GSE348100)
# Outputs:      results/DeLong_comparisons.tsv
# Key parameters:
#               test        pROC::roc.test, paired design (same participants)
#               family      all 16 paired signature comparisons form a SINGLE
#                           family; Benjamini-Hochberg is applied across those 16
#               reporting   both the unadjusted and the BH-adjusted P value;
#                           no significance stars are produced
#               interpretation  secondary and exploratory benchmark analyses,
#                           not a test of signature superiority
# Dependencies: R >= 4.4.3; data.table; pROC
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({library(data.table); library(pROC)})
source("code/_paths.R")

score_f <- file.path(DATA_BENCH, "SF6_Sweeney3_RISK6_per_sample.csv")
paths_require(score_f, "per-sample benchmark scores (TB-index, Sweeney3, RISK6)",
              "Supplementary Table S7 / NCBI GEO GSE348100")
sc <- fread(score_f)

CONTRASTS <- list(
  c("IDC", "APTB", "HC"), c("IDC", "APTB", "LTBI"), c("IDC", "APTB", "PNEU"),
  c("IDC", "APTB", "LC"), c("ITC", "APTB", "HC"), c("ITC", "APTB", "LC"),
  c("PVC", "PTB", "HC"), c("PVC", "EPTB", "HC"))

rows <- rbindlist(lapply(CONTRASTS, function(k) {
  d <- sc[Cohort == k[1] & Group %in% k[2:3]]
  y <- as.integer(d$Group == k[2])
  rbindlist(lapply(c("Sweeney3", "RISK6"), function(bm) {
    rt <- roc.test(roc(y, d$TB_index, quiet = TRUE, direction = "<"),
                   roc(y, d[[bm]],   quiet = TRUE, direction = "<"),
                   paired = TRUE, method = "delong")
    data.table(Cohort = k[1], Comparison = paste(k[2], "vs", k[3]), Benchmark = bm,
               AUROC_TB_index = as.numeric(rt$roc1$auc),
               AUROC_benchmark = as.numeric(rt$roc2$auc),
               P_value = rt$p.value)
  }))
}), fill = TRUE)

rows[, BH_adjusted_P_value := p.adjust(P_value, method = "BH")]   # single family of 16
setcolorder(rows, c("Cohort", "Comparison", "Benchmark", "AUROC_TB_index",
                    "AUROC_benchmark", "P_value", "BH_adjusted_P_value"))
fwrite(rows, file.path(RESULTS, "DeLong_comparisons.tsv"), sep = "\t")
message("Paired DeLong comparisons: ", nrow(rows), " (one BH family)")
