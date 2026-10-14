#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Purpose:      Assess the stability of the complete IDC discovery and
#               greedy-selection procedure using 500 stratified bootstrap
#               resamples.
# Inputs:       data/qn/IDC_QN_matrix.tsv      (genes x samples, QN log2(TPM+1))
#               data/metadata/IDC_metadata.tsv (group labels)
# Outputs:      results/bootstrap_panel_size.tsv            panel-size distribution
#               results/bootstrap_selection_frequency.tsv   per-gene selection frequency
#               results/bootstrap_auc_summary.tsv            in-bag / OOB AUROC
#               results/bootstrap_summary.tsv                headline summary
# Key parameters:
#               iterations  500
#               seed        2026
#               resampling  stratified within each clinical group, at the original
#                           group sizes (APTB 23 / HC 24 / LTBI 13 / PNEU 14 / LC 35)
#               per-iteration procedure
#                           1. stratified bootstrap resample
#                           2. four APTB-centred limma comparisons, each fitted on
#                              the two relevant groups
#                           3. screen |log2FC| > 0.5 and nominal P < 0.05
#                           4. four-way intersection
#                           5. direction concordance
#                           6. greedy forward selection, APTB vs pooled non-APTB
#                           7. stopping rule: incremental AUROC <= 0.02
#               The panel is re-derived in every iteration; no fixed path is used.
# Dependencies: R >= 4.4.3; data.table; limma; pROC
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({library(data.table); library(limma); library(pROC)})
source("code/_paths.R")

ITERS <- 500L
SEED  <- 2026L
## Optional overrides for a quick smoke test (defaults are the manuscript values):
##   Rscript code/07_bootstrap_stability.R [iterations] [seed]
.args <- commandArgs(trailingOnly = TRUE)
if (length(.args) >= 1) ITERS <- as.integer(.args[1])
if (length(.args) >= 2) SEED  <- as.integer(.args[2])
COMPARATORS <- c("HC", "LTBI", "PNEU", "LC")
DELTA_AUC_CUTOFF <- 0.02

## ---- 1. Inputs -------------------------------------------------------------
mat_f  <- file.path(DATA_QN, "IDC_QN_matrix.tsv")
meta_f <- file.path(DATA_META, "IDC_metadata.tsv")
paths_require(mat_f, "IDC QN log2(TPM+1) host-gene matrix (genes x samples)",
              "output of code/03_frozen_IDC_reference_normalization.R")
paths_require(meta_f, "IDC per-library metadata", "the de-identified cohort sample metadata")

mat  <- fread(mat_f)
meta <- fread(meta_f)
samples <- setdiff(names(mat), "Gene")
if (anyDuplicated(meta$Sample_ID)) stop("Duplicate Sample_ID in IDC metadata")
missing <- setdiff(samples, meta$Sample_ID); extra <- setdiff(meta$Sample_ID, samples)
if (length(missing) || length(extra)) stop("Sample ID mismatch (IDC)", call. = FALSE)
meta <- meta[match(samples, meta$Sample_ID)]

X <- as.matrix(mat[, ..samples]); rownames(X) <- mat$Gene   # genes x samples
grp <- setNames(as.character(meta$Group), meta$Sample_ID)

## Locked panel used for the in-bag / out-of-bag stability summaries.
lock <- paths_read_tsv(file.path(CONFIG, "locked_panel.tsv"))
LOCK_G <- lock$Gene; LOCK_W <- setNames(lock$Weight, lock$Gene)
stopifnot(all(LOCK_G %in% rownames(X)))
group_sizes <- table(grp)
message("Groups: ", paste(names(group_sizes), group_sizes, sep = "=", collapse = " "))

## ---- 2. One full discovery run on a sample set -----------------------------
screen <- function(Xs, gs) {
  m <- data.table(Sample_ID = colnames(Xs), Group = gs)
  de <- rbindlist(lapply(COMPARATORS, function(cmp) {
    keep <- m[Group %in% c("APTB", cmp), Sample_ID]
    g <- factor(m[match(keep, Sample_ID), Group], levels = c(cmp, "APTB"))
    fit <- eBayes(lmFit(Xs[, keep, drop = FALSE], model.matrix(~ g)))
    tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
    data.table(Gene = rownames(tt), Comparison = cmp, log2FC = tt$logFC, P = tt$P.Value)
  }))
  de[, hit := abs(log2FC) > 0.5 & P < 0.05]
  w <- dcast(de, Gene ~ Comparison, value.var = c("log2FC", "hit"))
  ## Restore the feature-matrix row order. The original procedure builds its
  ## candidate set with Reduce(intersect, ...), which preserves the order of the
  ## first argument (the expression-matrix row order). dcast sorts instead, and
  ## that reordering changes which gene wins an exact AUROC tie in the greedy
  ## step, so the order is restored here.
  w <- w[match(rownames(Xs), w$Gene)]
  hit_cols <- paste0("hit_", COMPARATORS); lfc_cols <- paste0("log2FC_", COMPARATORS)
  cand <- w[rowSums(w[, ..hit_cols]) == length(COMPARATORS)]
  if (!nrow(cand)) return(NULL)
  concordant <- apply(cand[, ..lfc_cols], 1, function(v) length(unique(sign(v))) == 1L)
  cand <- cand[concordant]
  if (!nrow(cand)) return(NULL)
  data.table(Gene = cand$Gene,
             Weight = ifelse(cand[[lfc_cols[1]]] > 0, 1, -1))
}

greedy <- function(Xs, y, cand) {
  genes <- cand$Gene; w <- setNames(cand$Weight, cand$Gene)
  auc_of <- function(sel) {
    s <- if (length(sel) == 1L) Xs[sel, ] * w[[sel]] else
      colSums(Xs[sel, , drop = FALSE] * w[sel])
    as.numeric(auc(roc(y, as.numeric(s), quiet = TRUE, direction = "<")))
  }
  sel <- character(0); path <- numeric(0)
  repeat {
    remaining <- setdiff(genes, sel)
    if (!length(remaining)) break
    gains <- vapply(remaining, function(g) auc_of(c(sel, g)), numeric(1))
    best <- names(which.max(gains))
    prev <- if (length(path)) tail(path, 1) else 0
    if (length(sel) >= 1L && (gains[[best]] - prev) <= DELTA_AUC_CUTOFF) break
    sel <- c(sel, best); path <- c(path, gains[[best]])
  }
  list(genes = sel, weights = w[sel], in_bag = if (length(path)) tail(path, 1) else NA_real_)
}

## ---- 3. Stratified bootstrap ----------------------------------------------
set.seed(SEED)
idx_by_group <- split(seq_along(grp), grp)
rows <- vector("list", ITERS); freq <- integer(0); sizes <- integer(ITERS)
in_bag <- numeric(ITERS); oob <- numeric(ITERS)

for (b in seq_len(ITERS)) {
  take <- unlist(lapply(idx_by_group, function(ii) sample(ii, length(ii), replace = TRUE)))
  Xs <- X[, take, drop = FALSE]; gs <- grp[take]
  cand <- screen(Xs, gs)
  sizes[b] <- if (is.null(cand)) 0L else {
    y <- as.integer(gs == "APTB")
    g <- greedy(Xs, y, cand)
    freq <- c(freq, setNames(rep(1L, length(g$genes)), g$genes))

    ## In-bag / out-of-bag AUROC of the LOCKED panel (config/locked_panel.tsv),
    ## which is what the bootstrap out-of-bag AUROC summary summarises; the per-iteration greedy panel is
    ## used only for the selection-frequency and panel-size summaries.
    score_of <- function(M) as.numeric(colSums(M * LOCK_W))
    in_bag[b] <- as.numeric(auc(roc(y, score_of(Xs[LOCK_G, , drop = FALSE]), quiet = TRUE, direction = "<")))
    oob_idx <- setdiff(seq_along(grp), unique(take))
    if (length(oob_idx)) {
      yo <- as.integer(grp[oob_idx] == "APTB")
      if (length(unique(yo)) == 2L) {
        oob[b] <- as.numeric(auc(roc(yo, score_of(X[LOCK_G, oob_idx, drop = FALSE]),
                                    quiet = TRUE, direction = "<")))
      }
    }
    rows[[b]] <- data.table(iteration = b, Gene = g$genes,
                            panel_size = length(g$genes), in_bag_auc = g$in_bag)
    length(g$genes)
  }
}

## ---- 4. Outputs ------------------------------------------------------------
sel_tab <- rbindlist(rows, fill = TRUE)
## Count selections per gene. Note: character indexing into a named vector with
## duplicated names returns only the first match, so the counts are built with
## table() rather than freq[Gene].
freq_tab <- as.data.table(table(names(freq)))
setnames(freq_tab, c("Gene", "n_selected"))
freq_tab[, n_selected := as.integer(n_selected)]
freq_tab[, pct := round(100 * n_selected / ITERS, 1)]
setorder(freq_tab, -n_selected, Gene)

size_tab <- data.table(panel_size = as.integer(names(table(sizes))),
                       n = as.integer(table(sizes)))
size_tab[, pct := round(100 * n / ITERS, 1)]

auc_tab <- data.table(
  metric = c("median_in_bag_AUROC", "median_OOB_AUROC",
             "median_panel_size", "panel_size_IQR_low", "panel_size_IQR_high",
             "panel_size_le3_pct", "panel_size_3_pct"),
  value = c(round(median(in_bag, na.rm = TRUE), 3), round(median(oob, na.rm = TRUE), 3),
            median(sizes), quantile(sizes, .25), quantile(sizes, .75),
            round(100 * mean(sizes <= 3), 1), round(100 * mean(sizes == 3), 1))
)

consensus <- head(freq_tab$Gene, 3)
summary_tab <- rbindlist(list(
  data.table(Item = "iterations", Value = ITERS),
  data.table(Item = "seed", Value = SEED),
  data.table(Item = "median_panel_size", Value = median(sizes)),
  data.table(Item = "panel_size_IQR", Value = sprintf("%d-%d", quantile(sizes, .25), quantile(sizes, .75))),
  data.table(Item = "panel_size_le3_pct", Value = round(100 * mean(sizes <= 3), 1)),
  rbindlist(lapply(c("SOD2", "ASXL1", "GOLPH3", "PCNX2"), function(g)
    data.table(Item = paste0("selection_frequency_", g),
               Value = freq_tab[Gene == g, pct][1]))),
  data.table(Item = "median_in_bag_AUROC", Value = round(median(in_bag, na.rm = TRUE), 3)),
  data.table(Item = "median_OOB_AUROC", Value = round(median(oob, na.rm = TRUE), 3)),
  data.table(Item = "consensus_panel", Value = paste(consensus, collapse = "+"))
))

fwrite(sel_tab,   file.path(RESULTS, "bootstrap_panel_size.tsv"), sep = "\t")
fwrite(freq_tab,  file.path(RESULTS, "bootstrap_selection_frequency.tsv"), sep = "\t")
fwrite(auc_tab,   file.path(RESULTS, "bootstrap_auc_summary.tsv"), sep = "\t")
fwrite(summary_tab, file.path(RESULTS, "bootstrap_summary.tsv"), sep = "\t")

## ---- 5. QA against the manuscript (the bootstrap-stability summaries) ---------------------------
EXP <- list(median_panel_size = 3, panel_size_le3_pct = 91.4, panel_size_3_pct = 51.8,
            SOD2 = 49.8, ASXL1 = 28.4, GOLPH3 = 17.8, PCNX2 = 5.0,
            median_in_bag = 0.952, median_OOB = 0.953)
fail <- character(0)
chk <- function(label, got, exp, tol = 0.05) {
  if (is.na(got) || abs(got - exp) > tol) fail <<- c(fail, sprintf("%s: %.3f vs %.3f", label, got, exp))
}
chk("median panel size", median(sizes), EXP$median_panel_size, 0.5)
chk("panel size <= 3 (%)", round(100 * mean(sizes <= 3), 1), EXP$panel_size_le3_pct, 0.1)
chk("panel size 3 (%)", round(100 * mean(sizes == 3), 1), EXP$panel_size_3_pct, 0.1)
## Selection frequencies are checked to within one bootstrap selection (0.2 pp at
## 500 iterations). The greedy step breaks exact AUROC ties by candidate order, so
## a one-selection difference is expected; the observed difference is reported by
## the message below rather than absorbed silently.
for (g in c("SOD2", "ASXL1", "GOLPH3", "PCNX2")) {
  got <- freq_tab[Gene == g, pct][1]
  chk(paste("freq", g), got, EXP[[g]], 0.25)
  if (abs(got - EXP[[g]]) > 0.05)
    message(sprintf("  note: %s selection frequency %.1f%% vs %.1f%% in the manuscript (%.0f vs %.0f of %d iterations)",
                    g, got, EXP[[g]], got * ITERS / 100, EXP[[g]] * ITERS / 100, ITERS))
}
chk("median in-bag AUROC", round(median(in_bag, na.rm = TRUE), 3), EXP$median_in_bag, 0.001)
chk("median OOB AUROC", round(median(oob, na.rm = TRUE), 3), EXP$median_OOB, 0.001)
if (length(fail)) stop("Bootstrap summary does not reproduce the bootstrap-stability summaries:\n  ", paste(fail, collapse = "\n  "), call. = FALSE)
message("QA: bootstrap summary reproduces the bootstrap-stability summaries; consensus panel = ", paste(consensus, collapse = "+"))
