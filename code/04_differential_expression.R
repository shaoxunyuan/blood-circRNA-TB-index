#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Purpose:      Four IDC APTB-centred limma comparisons; complete gene-level
#               statistics including nominal and BH-adjusted P values.
# Inputs:       data/qn/IDC_QN_matrix.tsv        (genes x samples, QN log2(TPM+1))
#               data/metadata/IDC_metadata.tsv
# Outputs:      results/DE_complete.tsv
# Key parameters:
#               contrasts      APTB vs HC / LTBI / PNEU / LC, each fitted on the
#                              two relevant clinical groups only (manuscript Methods)
#               screen         |log2FC| > 0.5 and nominal P < 0.05
#               multiplicity   Benjamini-Hochberg applied within each contrast
#                              across the 2,754 genes (one family per contrast);
#                              retained as secondary sensitivity information
# Dependencies: R >= 4.4.3; limma; data.table
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({library(data.table); library(limma)})
source("code/_paths.R")

EXPECTED_GENES <- 2754L
## Up/down counts of the prespecified screen, as reported in the manuscript
## (the differential-expression screen counts). Checked after the fit; a mismatch stops the run.
EXPECTED_DE_COUNTS <- list(
  "APTB vs HC"   = c(up = 117L, down = 146L),
  "APTB vs LTBI" = c(up = 141L, down = 151L),
  "APTB vs PNEU" = c(up = 208L, down = 235L),
  "APTB vs LC"   = c(up = 199L, down = 236L)
)
CONTRASTS <- c("HC", "LTBI", "PNEU", "LC")

## ---- 1. Expression matrix (genes x samples) -------------------------------
mat_f  <- file.path(DATA_QN, "IDC_QN_matrix.tsv")
meta_f <- file.path(DATA_META, "IDC_metadata.tsv")
paths_require(mat_f, "IDC QN log2(TPM+1) host-gene matrix (genes x samples)",
              "output of code/03_frozen_IDC_reference_normalization.R")
paths_require(meta_f, "IDC per-library metadata", "the de-identified cohort sample metadata")

mat  <- fread(mat_f)
meta <- fread(meta_f)
stopifnot(!anyDuplicated(meta$Sample_ID))

genes <- mat$Gene
sample_cols <- setdiff(names(mat), "Gene")
expr_matrix <- as.matrix(mat[, ..sample_cols])   # genes x samples
rownames(expr_matrix) <- genes
stopifnot(ncol(expr_matrix) == length(sample_cols))
stopifnot(nrow(expr_matrix) == EXPECTED_GENES)

## ---- 2. Sample-ID alignment ----------------------------------------------
if (anyDuplicated(meta$Sample_ID)) stop("Duplicate Sample_ID in IDC metadata")
if (anyDuplicated(colnames(expr_matrix))) stop("Duplicate sample columns in the matrix")
missing_in_expr <- setdiff(meta$Sample_ID, colnames(expr_matrix))
extra_in_expr   <- setdiff(colnames(expr_matrix), meta$Sample_ID)
if (length(missing_in_expr) || length(extra_in_expr)) {
  stop("Sample ID mismatch: ", length(missing_in_expr), " sample(s) missing from the ",
       "expression matrix, ", length(extra_in_expr), " extra.", call. = FALSE)
}
expr_matrix <- expr_matrix[, meta$Sample_ID, drop = FALSE]
stopifnot(identical(colnames(expr_matrix), meta$Sample_ID))

## ---- 3. Contrasts: each restricted to its two relevant groups -------------
## Manuscript Methods: "APTB was compared separately with HC, LTBI, PNEU, and LC,
## with each contrast restricted to samples belonging to the two relevant
## clinical groups." A single multi-group design does not reproduce the reported
## P values, because the residual variance would be pooled across all five
## groups; the per-contrast two-group fit is the analysis that was performed.
meta[, Group := factor(Group, levels = c("APTB", "HC", "LTBI", "PNEU", "LC"))]

fit_contrast <- function(cmp) {
  m <- meta[Group %in% c("APTB", cmp)]
  X <- expr_matrix[, m$Sample_ID, drop = FALSE]
  g <- factor(m$Group, levels = c(cmp, "APTB"))
  design <- model.matrix(~ g)
  fit <- eBayes(lmFit(X, design))
  topTable(fit, coef = 2, number = Inf, adjust.method = "BH", sort.by = "none")
}

## ---- 4. Complete statistics per contrast ---------------------------------
res <- rbindlist(lapply(CONTRASTS, function(cmp) {
  tt <- fit_contrast(cmp)
  data.table(
    Gene          = rownames(tt),
    Comparison    = paste("APTB vs", cmp),
    log2FC        = tt$logFC,
    AveExpr       = tt$AveExpr,
    t             = tt$t,
    P_value       = tt$P.Value,
    BH_adjusted_P = tt$adj.P.Val,
    B             = tt$B,
    Direction     = ifelse(tt$logFC > 0, "Up", "Down"),
    Screened      = abs(tt$logFC) > 0.5 & tt$P.Value < 0.05
  )
}), fill = TRUE)

setcolorder(res, c("Gene", "Comparison", "log2FC", "AveExpr", "t",
                   "P_value", "BH_adjusted_P", "B", "Direction", "Screened"))
fwrite(res, file.path(RESULTS, "DE_complete.tsv"), sep = "\t")

## ---- 5. QA: the prespecified screen still reproduces the reported counts ---
stopifnot(nrow(res) == EXPECTED_GENES * length(CONTRASTS))
for (cmp in names(EXPECTED_DE_COUNTS)) {
  d <- res[Comparison == cmp & Screened == TRUE]
  n_up <- sum(d$Direction == "Up"); n_down <- sum(d$Direction == "Down")
  exp <- EXPECTED_DE_COUNTS[[cmp]]
  if (n_up != exp[["up"]] || n_down != exp[["down"]]) {
    stop(sprintf("DE count mismatch for %s: observed up=%d down=%d, expected up=%d down=%d",
                 cmp, n_up, n_down, exp[["up"]], exp[["down"]]), call. = FALSE)
  }
  message(sprintf("  %-14s screened up=%3d down=%3d  (matches the differential-expression screen counts)",
                  cmp, n_up, n_down))
}
message("QA: ", EXPECTED_GENES, " genes x ", length(CONTRASTS),
        " contrasts; BH within each contrast; screen = |log2FC| > 0.5 & nominal P < 0.05")
