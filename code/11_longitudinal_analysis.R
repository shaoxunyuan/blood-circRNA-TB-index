#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Purpose:      PVC longitudinal analysis: within-participant TB-index
#               trajectories and exploratory time-point discrimination.
# Inputs:       results/locked_TB_index_scores.tsv
# Outputs:      results/longitudinal_analysis.tsv
# Key parameters:
#               paired test    two-sided Wilcoxon signed-rank on adjacent time
#                              points (BL vs M1, M1 vs M6), PTB and EPTB separately
#               AUROC          ordinary two-class ROC on paired samples
#                              (case = BL); this is not a paired-ROC estimator and
#                              the intervals do not account for within-participant
#                              dependence
#               All longitudinal analyses are exploratory.
# Dependencies: R >= 4.4.3; data.table; pROC
# ---------------------------------------------------------------------------
suppressPackageStartupMessages({library(data.table); library(pROC)})
source("code/_paths.R")

s <- fread(file.path(RESULTS, "locked_TB_index_scores.tsv"))
s <- s[Cohort == "PVC" & !is.na(TimePoint)]

paired <- rbindlist(lapply(c("PTB", "EPTB"), function(g) {
  d <- dcast(s[Group == g], Participant_ID ~ TimePoint, value.var = "TB_index")
  rbindlist(lapply(list(c("BL", "M1"), c("M1", "M6")), function(tp) {
    if (!all(tp %in% names(d))) return(NULL)
    v <- d[!is.na(get(tp[1])) & !is.na(get(tp[2]))]
    if (nrow(v) < 2) return(NULL)
    w <- suppressWarnings(wilcox.test(v[[tp[1]]], v[[tp[2]]], paired = TRUE))
    data.table(Group = g, Comparison = paste0(tp[1], " vs ", tp[2]),
               n_pairs = nrow(v), P_value = w$p.value)
  }), fill = TRUE)
}), fill = TRUE)

fwrite(paired, file.path(RESULTS, "longitudinal_analysis.tsv"), sep = "\t")
message("Longitudinal comparisons: ", nrow(paired),
        " (exploratory; complete-pair counts are reported in Supplementary Table S6B)")
