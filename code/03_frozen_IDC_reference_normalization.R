#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Purpose:      Fit the frozen IDC-reference quantile normalization (QN) on the
#               IDC discovery cohort and apply it unchanged to IDC, ITC and PVC.
# Inputs:       data/qn/{IDC,ITC,PVC}_common_2754.tsv   (output of script 02;
#                 WIDE matrix, genes x samples, raw TPM)
#               data/metadata/{IDC,ITC,PVC}_metadata.tsv
# Outputs:      data/qn/QN_reference_IDC.tsv      (rank, qn_rank_mean)
#               data/qn/{IDC,ITC,PVC}_QN_matrix.tsv  (genes x samples,
#                 log2(TPM+1) mapped onto the frozen IDC reference)
# Key params:   Reference fitted on the IDC only; no ITC or PVC sample enters
#               the fit. Within each sample genes are ranked with
#               rank(..., ties.method = "first"); the cross-sample mean
#               log2(TPM+1) at each rank position is the frozen reference.
#               Every sample of every cohort is mapped onto that reference by
#               replacing its rank positions with the reference values.
# Dependencies: R >= 4.4.3; data.table
#
# Input format contract (see docs/reproducibility.md): script 03 consumes the
# common-gene matrices written by script 02 in wide format.
# ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(data.table))
source("code/_paths.R")

COHORTS <- c("IDC", "ITC", "PVC")
EXPECTED_LIBRARIES <- c(IDC = 109L, ITC = 35L, PVC = 81L)
EXPECTED_GENES <- 2754L

## ---- 1. Common-gene matrices from script 02 -------------------------------
common_f <- file.path(RESULTS, "common_feature_genes.tsv")
paths_require(common_f, "common feature space gene list (output of script 02)",
              "run code/02_build_common_feature_matrix.R")
common <- fread(common_f)$Gene
stopifnot(length(common) == EXPECTED_GENES)

read_common <- function(coh) {
  f <- file.path(DATA_QN, sprintf("%s_common_%d.tsv", coh, EXPECTED_GENES))
  paths_require(f, paste0(coh, " common-gene matrix (output of script 02)"),
                "run code/02_build_common_feature_matrix.R")
  d <- fread(f)
  stopifnot(identical(d$Gene, common))
  samples <- setdiff(names(d), "Gene")
  stopifnot(length(samples) == EXPECTED_LIBRARIES[[coh]])
  X <- as.matrix(d[, ..samples]); rownames(X) <- d$Gene
  X
}

## ---- 2. Sample-ID alignment against the metadata --------------------------
samples_of <- function(coh) {
  f <- file.path(DATA_META, paste0(coh, "_metadata.tsv"))
  paths_require(f, paste0(coh, " metadata"), "the de-identified cohort sample metadata")
  m <- fread(f)
  stopifnot(!anyDuplicated(m$Sample_ID))
  m
}

tpm <- list(); meta <- list()
for (coh in COHORTS) {
  meta[[coh]] <- samples_of(coh)
  X <- read_common(coh)
  missing_in_meta <- setdiff(colnames(X), meta[[coh]]$Sample_ID)
  extra_in_meta   <- setdiff(meta[[coh]]$Sample_ID, colnames(X))
  if (length(missing_in_meta) || length(extra_in_meta)) {
    stop("Sample ID mismatch for ", coh, ": ", length(missing_in_meta),
         " not in metadata, ", length(extra_in_meta), " extra in metadata.",
         call. = FALSE)
  }
  ## Reorder explicitly onto the metadata order.
  tpm[[coh]] <- X[, meta[[coh]]$Sample_ID, drop = FALSE]
  message(sprintf("  %s: %d genes x %d samples", coh, nrow(tpm[[coh]]), ncol(tpm[[coh]])))
}

## ---- 3. log2(TPM + 1) -----------------------------------------------------
logtpm <- lapply(tpm, function(X) log2(X + 1))

## ---- 4. Frozen IDC reference ---------------------------------------------
## Rank genes within each IDC sample (ties.method = "first"), then average the
## log2(TPM+1) value at each rank position across IDC samples.
idc <- logtpm$IDC
## An ascending stable sort gives the same ordering as rank(..., ties.method =
## "first") applied to the untransformed values, so the reference is the mean
## log2(TPM+1) at each rank position across the IDC samples.
ranked <- vapply(seq_len(ncol(idc)), function(j) idc[order(idc[, j]), j], numeric(nrow(idc)))
qn_ref <- rowMeans(ranked)
stopifnot(length(qn_ref) == EXPECTED_GENES)

fwrite(data.table(rank = seq_along(qn_ref), qn_rank_mean = qn_ref),
       file.path(DATA_QN, "QN_reference_IDC.tsv"), sep = "\t")

## ---- 5. Map every sample of every cohort onto the reference ---------------
qn_map <- function(X) {
  out <- vapply(seq_len(ncol(X)), function(j) {
    r <- rank(X[, j], ties.method = "first")
    qn_ref[r]
  }, numeric(nrow(X)))
  dimnames(out) <- dimnames(X)
  out
}

for (coh in COHORTS) {
  Q <- qn_map(logtpm[[coh]])
  out <- data.table(Gene = rownames(Q), Q)
  fwrite(out, file.path(DATA_QN, paste0(coh, "_QN_matrix.tsv")), sep = "\t")
  message(sprintf("  wrote data/qn/%s_QN_matrix.tsv (%d genes x %d samples)",
                  coh, nrow(Q), ncol(Q)))
}

## ---- 6. QA ----------------------------------------------------------------
for (coh in COHORTS) {
  X <- fread(file.path(DATA_QN, paste0(coh, "_QN_matrix.tsv")))
  stopifnot(nrow(X) == 2754L)
  stopifnot(ncol(X) - 1L == EXPECTED_LIBRARIES[[coh]])
}
message("QA: QN reference = ", length(qn_ref), " values; matrices ",
        paste(sprintf("%s=%d", COHORTS, EXPECTED_LIBRARIES), collapse = " "))
