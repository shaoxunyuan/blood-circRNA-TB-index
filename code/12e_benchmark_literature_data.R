#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Purpose:      Literature-formula scores for the benchmark (the literature-benchmark results).
# Inputs:       Benchmark workbook from 12a (per_sample_scores with lin_* QN columns).
# Outputs:      output/benchmark_literature.xlsx with Sweeney3 / RISK6 scores.
# Key params:   Sweeney3 = (GBP5+DUSP3)/2 - KLF2;
#               RISK6 = GM(GBP2,FCGR1B,SERPING1) - GM(TUBGCP6,TRMT2A,SDR39U1)
#               (QN log2(TPM+1); GM on linear scale). Column name RISK6.
# Dependencies: R >= 4.4.3; openxlsx; data.table.
# ---------------------------------------------------------------------------
# =============================================================================
# build_benchmark_literature_data.R — 文献公式版：Sweeney3 + RISK6
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
})


## Package-relative paths
.args_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.code_dir <- if (length(.args_file)) dirname(normalizePath(sub("^--file=", "", .args_file[1]))) else getwd()
source(file.path(.code_dir, "_paths.R"))

S8       <- file.path(DATA_BENCH, "benchmark_tbindex.xlsx")
OUT_FILE <- file.path(OUT, "benchmark_literature.xlsx")

stopifnot(file.exists(S8))

gm3_log2 <- function(a, b, c) 2^((as.numeric(a) + as.numeric(b) + as.numeric(c)) / 3)
risk6geo <- function(gbp2, fcgr1b, serping1, tubgcp6, trmt2a, sdr39u1) {
  gm3_log2(gbp2, fcgr1b, serping1) - gm3_log2(tubgcp6, trmt2a, sdr39u1)
}

cat("================================================================\n")
cat(" Literature-formula benchmark: Sweeney3 + RISK6\n")
cat("================================================================\n")

file.copy(S8, OUT_FILE, overwrite = TRUE)
wb <- loadWorkbook(OUT_FILE)
pss  <- as.data.table(read.xlsx(OUT_FILE, sheet = "per_sample_scores"))
ss   <- as.data.table(read.xlsx(OUT_FILE, sheet = "Sample_scores"))

pss[, Sweeney3_lit := (lin_GBP5 + lin_DUSP3)/2 - lin_KLF2]
pss[, RISK6_lit := risk6geo(lin_GBP2, lin_FCGR1B, lin_SERPING1,
                            lin_TUBGCP6, lin_TRMT2A, lin_SDR39U1)]
pss[, Sweeney3 := round(Sweeney3_lit, 4)]
pss[, RISK6 := round(RISK6_lit, 4)]
pss[, `:=`(Sweeney3_lit = NULL, RISK6_lit = NULL)]

m <- merge(ss[, .(Sample_ID, Cohort)], pss[, .(Sample_ID, Cohort, Sweeney3, RISK6)],
           by = c("Sample_ID", "Cohort"), all.x = TRUE)
ss[, Sweeney3 := m$Sweeney3]
ss[, RISK6 := m$RISK6]

removeWorksheet(wb, "per_sample_scores"); addWorksheet(wb, "per_sample_scores")
writeData(wb, "per_sample_scores", pss)
setColWidths(wb, "per_sample_scores", cols = 1:ncol(pss), widths = "auto")
writeComment(wb, "per_sample_scores", col = 1, row = 1,
  comment = createComment(paste0(
    "Literature-formula scores:\n",
    "Sweeney3 = (GBP5+DUSP3)/2 - KLF2; ",
    "RISK6 = GM(GBP2,FCGR1B,SERPING1) - GM(TUBGCP6,TRMT2A,SDR39U1) ",
    "(QN log2(TPM+1) -> GM on linear scale); circRNA TB-index unchanged.")))

removeWorksheet(wb, "Sample_scores"); addWorksheet(wb, "Sample_scores")
writeData(wb, "Sample_scores", ss)
setColWidths(wb, "Sample_scores", cols = 1:ncol(ss), widths = "auto")

sm <- as.data.table(read.xlsx(OUT_FILE, sheet = "Signature_metadata"))
sm[Signature == "Sweeney3", Formula := "(GBP5 + DUSP3)/2 - KLF2 (linear log2TPM QN)"]
sm[Signature %in% c("RISK6", "adapted RISK6", "RISK6geo"), Signature := "RISK6"]
sm[Signature == "RISK6", Formula :=
     "GM(GBP2,FCGR1B,SERPING1) - GM(TUBGCP6,TRMT2A,SDR39U1) (QN log2(TPM+1) -> GM on linear scale)"]
sm[Signature == "RISK6", Reference :=
     "Penn-Nicholson et al. Sci Rep 2020 (RNA-seq/microarray geometric-mean formulation)"]
sm[Signature != "circRNA TB-index", Reference := paste0(Reference, " (literature formula)")]
removeWorksheet(wb, "Signature_metadata"); addWorksheet(wb, "Signature_metadata")
writeData(wb, "Signature_metadata", sm)
setColWidths(wb, "Signature_metadata", cols = 1:ncol(sm), widths = "auto")

for (s in c("AUC_benchmark", "Delong_test", "Joint_model", "Joint_model_DeLong")) {
  if (s %in% names(wb)) removeWorksheet(wb, s)
}

nms <- names(wb)
wanted <- c("per_sample_scores", "Sample_scores", "Gene_direction_circular", "Gene_direction_linear",
            "Sweeney3_genes_expression", "Sweeney3_DE_ATB",
            "RISK6_genes_expression", "RISK6_DE_ATB",
            "Expression_matrix_3circ_9lin",
            "Signature_metadata")
ord <- unlist(lapply(wanted, function(s) which(nms == s)))
worksheetOrder(wb) <- ord

saveWorkbook(wb, OUT_FILE, overwrite = TRUE)
cat(sprintf("\nDone → %s (%d sheets)\n", OUT_FILE, length(names(wb))))
cat("  Sweeney3 = (GBP5+DUSP3)/2 - KLF2\n")
cat("  RISK6 = GM(GBP2,FCGR1B,SERPING1) - GM(TUBGCP6,TRMT2A,SDR39U1)\n")
cat("================================================================\n")
