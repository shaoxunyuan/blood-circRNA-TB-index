#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Purpose:      Greedy forward selection over the 18 direction-concordant
#               candidates; locks the three-gene panel.
# Inputs:       data/qn/IDC_QN_matrix.tsv        (genes x samples, QN log2(TPM+1))
#               results/core18_candidates.tsv     (output of script 05)
#               data/metadata/IDC_metadata.tsv
# Outputs:      results/greedy_selection_path.tsv
#               config/locked_panel.tsv           (SOD2, GOLPH3, PCNX2)
# Key parameters:
#               candidate set  18 four-way direction-concordant genes
#               endpoint       APTB vs pooled non-APTB AUROC
#               weights        equal magnitude: +1 for genes upregulated in APTB,
#                              -1 for genes downregulated in APTB
#               step 1         the single gene with the highest AUROC
#               later steps    the remaining gene giving the largest increase in
#                              cumulative AUROC
#               stopping       the best achievable incremental AUROC <= 0.02; the
#                              subthreshold candidate is NOT added to the panel
# Dependencies: R >= 4.4.3; data.table; pROC
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({library(data.table); library(pROC)})
source("code/_paths.R")

EXPECTED_CANDIDATES <- 18L
DELTA_AUC_CUTOFF <- 0.02
## Manuscript Figure 2D reference path (tolerance-checked, never hard-coded as output).
EXPECTED_PATH <- data.table(
  gene = c("SOD2", "GOLPH3", "PCNX2"),
  cumulative_auc = c(0.863, 0.924, 0.949)
)
TOL <- 1e-3

## ---- 1. Inputs -------------------------------------------------------------
mat_f  <- file.path(DATA_QN, "IDC_QN_matrix.tsv")
meta_f <- file.path(DATA_META, "IDC_metadata.tsv")
cand_f <- file.path(RESULTS, "core18_candidates.tsv")
paths_require(mat_f,  "IDC QN log2(TPM+1) host-gene matrix (genes x samples)",
              "output of code/03_frozen_IDC_reference_normalization.R")
paths_require(meta_f, "IDC per-library metadata", "Supplementary Table S1")
paths_require(cand_f, "direction-concordant candidates (output of script 05)",
              "run code/05_four_way_candidate_intersection.R")

mat  <- fread(mat_f)
meta <- fread(meta_f)
cand <- fread(cand_f)

genes <- cand$Gene
stopifnot(length(genes) == EXPECTED_CANDIDATES)
stopifnot(all(genes %in% mat$Gene))

## ---- 2. Sample-ID alignment ------------------------------------------------
stopifnot(!anyDuplicated(meta$Sample_ID))
samples <- setdiff(names(mat), "Gene")
missing <- setdiff(meta$Sample_ID, samples); extra <- setdiff(samples, meta$Sample_ID)
if (length(missing) || length(extra)) stop("Sample ID mismatch (IDC)", call. = FALSE)
meta <- meta[match(samples, meta$Sample_ID)]

## ---- 3. Candidate submatrix: rows = IDC samples, columns = 18 genes --------
sub <- mat[Gene %in% genes]
sub <- sub[match(genes, sub$Gene)]           # lock the candidate order
stopifnot(nrow(sub) == EXPECTED_CANDIDATES)
sample_cols <- setdiff(names(sub), "Gene")
Xg <- as.matrix(sub[, ..sample_cols]); rownames(Xg) <- sub$Gene   # genes x samples
X <- t(Xg)                                   # samples x genes
stopifnot(nrow(X) == length(samples), ncol(X) == EXPECTED_CANDIDATES)
stopifnot(identical(colnames(X), genes))

## ---- 4. Direction weights (equal magnitude, signed by APTB direction) -----
w <- setNames(ifelse(cand$Direction == "Up", 1, -1), cand$Gene)
y <- as.integer(meta$Group == "APTB")         # APTB vs pooled non-APTB

score_of <- function(sel) as.numeric(X[, sel, drop = FALSE] %*% w[sel])
auc_of   <- function(sel) as.numeric(auc(roc(y, score_of(sel), quiet = TRUE, direction = "<")))

## ---- 5. Greedy forward selection ------------------------------------------
selected <- character(0); path <- list()
repeat {
  remaining <- setdiff(genes, selected)
  gain <- vapply(remaining, function(g) auc_of(c(selected, g)), numeric(1))
  best <- names(which.max(gain))
  prev <- if (length(path)) path[[length(path)]]$cumulative_auc else 0
  if (length(selected) >= 1 && (gain[[best]] - prev) <= DELTA_AUC_CUTOFF) break
  selected <- c(selected, best)
  path[[length(path) + 1]] <- data.table(step = length(selected), gene = best,
                                         cumulative_auc = gain[[best]])
}

sel_path <- rbindlist(path)
fwrite(sel_path, file.path(RESULTS, "greedy_selection_path.tsv"), sep = "\t")

lock <- data.table(Gene = selected,
                   Direction = ifelse(w[selected] > 0, "Up", "Down"),
                   Weight = unname(w[selected]))
fwrite(lock, file.path(CONFIG, "locked_panel.tsv"), sep = "\t", quote = FALSE)

## ---- 6. QA: the locked panel and its path reproduce Figure 2D -------------
if (!identical(selected, EXPECTED_PATH$gene)) {
  stop("Locked panel mismatch: obtained [", paste(selected, collapse = ", "),
       "], expected [", paste(EXPECTED_PATH$gene, collapse = ", "), "]", call. = FALSE)
}
for (i in seq_len(nrow(EXPECTED_PATH))) {
  got <- sel_path$cumulative_auc[i]; exp <- EXPECTED_PATH$cumulative_auc[i]
  if (abs(got - exp) > TOL) {
    stop(sprintf("Greedy AUROC mismatch at step %d (%s): observed %.4f, expected %.4f",
                 i, EXPECTED_PATH$gene[i], got, exp), call. = FALSE)
  }
  message(sprintf("  step %d: %-7s cumulative AUROC = %.4f (expected %.3f)",
                  i, sel_path$gene[i], got, exp))
}
message("QA: locked panel = ", paste(selected, collapse = " / "))
