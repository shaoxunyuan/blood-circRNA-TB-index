#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Purpose:      Assemble the lightweight summary tables shipped in results/ and
#               verify that they are consistent with config/ and each other.
# Inputs:       config/*.tsv
#               results/{core18_candidates,locked_TB_index_scores,
#                        diagnostic_performance,bootstrap_summary,
#                        DeLong_comparisons,BSJ_summary}.tsv
# Outputs:      results/summary_consistency_check.tsv
# Key parameters: none; this script only reads and checks.
# Dependencies: R >= 4.4.3; data.table
# ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(data.table))
source("code/_paths.R")

params <- fread(file.path(CONFIG, "analysis_parameters.tsv"))
panel  <- paths_read_tsv(file.path(CONFIG, "locked_panel.tsv"))
thr    <- fread(file.path(CONFIG, "locked_thresholds.tsv"))
c18    <- fread(file.path(RESULTS, "core18_candidates.tsv"))
scores <- fread(file.path(RESULTS, "locked_TB_index_scores.tsv"))
dp     <- fread(file.path(RESULTS, "diagnostic_performance.tsv"))
dl     <- fread(file.path(RESULTS, "DeLong_comparisons.tsv"))
bsj_f  <- file.path(RESULTS, "BSJ_summary.tsv")
paths_require(bsj_f, "BSJ summary (output of code/14_BSJ_refinement.R)",
              "run code/14_BSJ_refinement.R, or NCBI GEO GSE348100")
bsj    <- paths_read_tsv(bsj_f)

val <- function(p) as.numeric(params[parameter == p, value])

# Common feature space: observed = rows of the gene list written by
# code/02_build_common_feature_matrix.R; expected = the locked parameter.
common_features <- fread(file.path(RESULTS, "common_feature_genes.tsv"))
observed_common <- nrow(common_features)
expected_common <- as.integer(params[parameter == "common_feature_count", value])

# All-three-cohort junction counts, read per gene from the BSJ summary written by
# code/14_BSJ_refinement.R (observed) and compared with the manuscript-locked values.
bsj_all3_per_gene <- vapply(c("SOD2", "GOLPH3", "PCNX2"),
                            function(g) as.character(bsj[Gene == g, All_three_cohort_BSJ_n][1]),
                            character(1))
obs_all3 <- paste(bsj_all3_per_gene, collapse = "/")

check <- rbindlist(list(
  data.table(item = "common feature space",  expected = expected_common,
             observed = observed_common, pass = observed_common == expected_common),
  data.table(item = "candidate genes",       expected = val("candidate_gene_count"), observed = nrow(c18), pass = nrow(c18) == val("candidate_gene_count")),
  data.table(item = "locked panel size",     expected = val("locked_panel_size"),    observed = nrow(panel), pass = nrow(panel) == val("locked_panel_size")),
  data.table(item = "scored libraries",      expected = 225,      observed = nrow(scores), pass = nrow(scores) == 225),
  data.table(item = "diagnostic comparisons",expected = 11,       observed = nrow(dp),     pass = nrow(dp) == 11),
  data.table(item = "DeLong comparisons",    expected = 16,       observed = nrow(dl),     pass = nrow(dl) == 16),
  data.table(item = "BSJ summary rows",      expected = 9,        observed = nrow(bsj),    pass = nrow(bsj) == 9),
  data.table(item = "BSJ all-three (SOD2/GOLPH3/PCNX2)",
             expected = "1/1/12",
             observed = obs_all3,
             pass = identical(obs_all3, "1/1/12")),
  data.table(item = "Youden threshold",      expected = 1.407, observed = thr[Threshold_name == "Youden", TB_index_cutoff], pass = thr[Threshold_name == "Youden", TB_index_cutoff] == 1.407),
  data.table(item = "high-sensitivity",      expected = 0.261, observed = thr[Threshold_name == "High_sensitivity", TB_index_cutoff], pass = thr[Threshold_name == "High_sensitivity", TB_index_cutoff] == 0.261)
), fill = TRUE)

fwrite(check, file.path(RESULTS, "summary_consistency_check.tsv"), sep = "\t")
print(check)
stopifnot(all(check$pass))
message("All summary tables consistent with config/.")
