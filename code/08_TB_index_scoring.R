#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Purpose:      Compute the locked three-gene TB-index for every library.
# Inputs:       data/qn/{IDC,ITC,PVC}_QN_matrix.tsv  (genes x samples)
#               data/metadata/{IDC,ITC,PVC}_metadata.tsv
#               config/locked_panel.tsv
# Outputs:      results/locked_TB_index_scores.tsv
# Key parameters:
#               TB-index = SOD2_QN - GOLPH3_QN - PCNX2_QN
#               The panel and its weights are locked in the IDC and applied
#               unchanged to ITC and PVC (no cohort-specific re-optimisation).
# Dependencies: R >= 4.4.3; data.table
# ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(data.table))
source("code/_paths.R")

EXPECTED_LIBRARIES <- c(IDC = 109L, ITC = 35L, PVC = 81L)
PANEL <- c("SOD2", "GOLPH3", "PCNX2")
## Locked weights, applied explicitly by gene name.
WEIGHTS <- c(SOD2 = 1, GOLPH3 = -1, PCNX2 = -1)

panel <- paths_read_tsv(file.path(CONFIG, "locked_panel.tsv"))
panel <- panel[Gene %in% PANEL]
stopifnot(nrow(panel) == length(PANEL))
panel <- panel[match(PANEL, panel$Gene)]
if (!identical(panel$Gene, PANEL) ||
    !isTRUE(all.equal(as.numeric(panel$Weight), as.numeric(WEIGHTS)))) {
  stop("config/locked_panel.tsv does not define the expected panel ",
       paste(PANEL, collapse = " / "), " with weights ",
       paste(WEIGHTS, collapse = " / "), call. = FALSE)
}

scores <- rbindlist(lapply(c("IDC", "ITC", "PVC"), function(coh) {
  mat_f  <- file.path(DATA_QN, paste0(coh, "_QN_matrix.tsv"))
  meta_f <- file.path(DATA_META, paste0(coh, "_metadata.tsv"))
  paths_require(mat_f,  paste0(coh, " QN log2(TPM+1) host-gene matrix (genes x samples)"),
                "output of code/03_frozen_IDC_reference_normalization.R")
  paths_require(meta_f, paste0(coh, " per-library metadata"), "the de-identified cohort sample metadata")

  mat  <- fread(mat_f)
  meta <- fread(meta_f)
  stopifnot(!anyDuplicated(meta$Sample_ID))
  stopifnot(all(PANEL %in% mat$Gene))

  samples <- setdiff(names(mat), "Gene")
  missing <- setdiff(meta$Sample_ID, samples); extra <- setdiff(samples, meta$Sample_ID)
  if (length(missing) || length(extra)) {
    stop("Sample ID mismatch for ", coh, ": ", length(missing),
         " missing from the matrix, ", length(extra), " extra.", call. = FALSE)
  }
  if (length(samples) != EXPECTED_LIBRARIES[[coh]]) {
    stop(sprintf("%s: expected %d libraries, found %d",
                 coh, EXPECTED_LIBRARIES[[coh]], length(samples)), call. = FALSE)
  }

  ## Explicit gene-by-gene extraction from a genes x samples matrix:
  ## no reliance on column order or on dynamic column selection.
  sub <- mat[Gene %in% PANEL]
  sub <- sub[match(PANEL, sub$Gene)]
  stopifnot(identical(sub$Gene, PANEL))
  M <- as.matrix(sub[, ..samples]); rownames(M) <- sub$Gene   # 3 genes x samples

  meta <- meta[match(samples, meta$Sample_ID)]
  stopifnot(identical(meta$Sample_ID, samples))

  data.table(
    Sample_ID = samples,
    Participant_ID =
      if ("Participant_ID" %in% names(meta)) as.character(meta$Participant_ID) else samples,
    Cohort    = coh,
    Group     = meta$Group,
    TimePoint = if ("TimePoint" %in% names(meta)) meta$TimePoint else NA_character_,
    Baseline_diagnostic_analysis =
      if ("Baseline_diagnostic_analysis" %in% names(meta)) meta$Baseline_diagnostic_analysis else NA_character_,
    SOD2_QN   = as.numeric(M["SOD2", ]),
    GOLPH3_QN = as.numeric(M["GOLPH3", ]),
    PCNX2_QN  = as.numeric(M["PCNX2", ]),
    TB_index  = as.numeric(M["SOD2", ] - M["GOLPH3", ] - M["PCNX2", ])
  )
}), fill = TRUE)

fwrite(scores, file.path(RESULTS, "locked_TB_index_scores.tsv"), sep = "\t")

## ---- QA --------------------------------------------------------------------
stopifnot(sum(is.na(scores$TB_index)) == 0)
for (coh in names(EXPECTED_LIBRARIES)) {
  n <- sum(scores$Cohort == coh)
  if (n != EXPECTED_LIBRARIES[[coh]]) {
    stop(sprintf("%s: expected %d scores, found %d", coh, EXPECTED_LIBRARIES[[coh]], n),
         call. = FALSE)
  }
  message(sprintf("  %s: %d scores", coh, n))
}
stopifnot(nrow(scores) == 225L)
message("QA: 225 libraries scored; TB-index = SOD2_QN - GOLPH3_QN - PCNX2_QN")
