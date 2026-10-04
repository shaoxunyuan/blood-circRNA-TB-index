#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Purpose:      Derive the 18 four-way direction-concordant candidate genes and
#               the per-contrast BH-adjusted P values reported for them.
# Inputs:       results/DE_complete.tsv
# Outputs:      results/core18_candidates.tsv
# Key parameters:
#               screen      |log2FC| > 0.5 and nominal P < 0.05 in ALL four
#                           APTB-centred contrasts
#               concordance the sign of log2FC must agree across the four contrasts
#               BH          adjusted within each contrast across the 2,754 genes
#               Discordant genes are excluded from downstream feature selection.
# Dependencies: R >= 4.4.3; data.table
# ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(data.table))
source("code/_paths.R")

de_f <- file.path(RESULTS, "DE_complete.tsv")
paths_require(de_f, "complete gene-level DE statistics (step 04 output)", "run code/04_differential_expression.R first")
de <- fread(de_f)
de[, screened := abs(log2FC) > 0.5 & P_value < 0.05]

sm <- dcast(de, Gene ~ Comparison, value.var = c("log2FC", "P_value", "BH_adjusted_P", "screened"))
lf  <- grep("^log2FC_", names(sm), value = TRUE)
sc  <- grep("^screened_", names(sm), value = TRUE)
bh  <- grep("^BH_adjusted_P_", names(sm), value = TRUE)

keep <- sm[rowSums(sm[, ..sc]) == length(sc)]
concordant <- apply(keep[, ..lf], 1, function(x) length(unique(sign(x))) == 1L)
candidates <- keep[concordant]

## Report the contrasts in the manuscript order: HC, LTBI, PNEU, LC.
ORD <- c("APTB vs HC", "APTB vs LTBI", "APTB vs PNEU", "APTB vs LC")
lf <- paste0("log2FC_", ORD); sc <- paste0("screened_", ORD); bh <- paste0("BH_adjusted_P_", ORD)
short <- sub("^log2FC_APTB vs ", "", lf)
out <- data.table(
  Gene = candidates$Gene,
  Direction = ifelse(candidates[[lf[1]]] > 0, "Up", "Down"))
for (i in seq_along(lf)) {
  out[[paste0("log2FC_APTB_vs_", short[i])]] <- candidates[[lf[i]]]
}
for (i in seq_along(bh)) {
  out[[paste0("BH_adjusted_P_APTB_vs_", short[i])]] <- candidates[[bh[i]]]
}
out[, BH_significant_comparisons_n :=
      rowSums(.SD < 0.05), .SDcols = grep("^BH_adjusted_P_", names(out), value = TRUE)]

fwrite(out, file.path(RESULTS, "core18_candidates.tsv"), sep = "\t")
message("Direction-concordant candidates: ", nrow(out))
