#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Purpose:      Age- and sex-adjusted sensitivity analyses for the locked TB-index.
# Inputs:       results/locked_TB_index_scores.tsv
#               metadata/sample_dictionary.tsv  (Sex, Age)
# Outputs:      results/demographic_adjustment.tsv
# Key parameters:
#               residual model   TB-index ~ Age + Sex, fitted on all IDC
#                                participants with complete values
#               logistic model   disease ~ TB-index + Age + Sex, fitted per
#                                contrast; ordinary ML logistic regression is
#                                retained when the model converges with maximum
#                                standard error < 20 and no fitted probability
#                                outside (1e-6, 1 - 1e-6), otherwise
#                                Firth bias-reduced logistic regression
#               Sex encoding     male = 1, female = 0
# Dependencies: R >= 4.4.3; data.table
# ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(data.table))
source("code/_paths.R")

s_f <- file.path(RESULTS, "locked_TB_index_scores.tsv")
paths_require(s_f, "locked TB-index scores (step 08 output)",
              "run code/08_TB_index_scoring.R, or NCBI GEO GSE348100")
s <- fread(s_f)

# Sex/Age are not carried in the score table; join the de-identified metadata.
meta_f <- file.path(ROOT, "metadata", "sample_dictionary.tsv")
paths_require(meta_f, "de-identified sample dictionary (Sex, Age)",
              "step 01 output; Supplementary Table S1")
meta <- fread(meta_f)
s <- merge(s, meta[, .(Sample_ID, Sex, Age)], by = "Sample_ID", all.x = TRUE, sort = FALSE)
s[, Sex01 := ifelse(Sex == "Male", 1L, ifelse(Sex == "Female", 0L, NA_integer_))]

idc <- s[Cohort == "IDC" & !is.na(Age) & !is.na(Sex01)]
resid_model <- lm(TB_index ~ Age + Sex01, data = idc)
idc[, TB_index_adj := resid(resid_model)]

out <- rbindlist(lapply(unique(idc$Group), function(g) {
  d <- idc[Group %in% c("APTB", g)]
  if (uniqueN(d$Group) < 2) return(NULL)
  d[, Disease := as.integer(Group == "APTB")]
  fit <- glm(Disease ~ TB_index + Age + Sex01, data = d, family = binomial)
  co <- summary(fit)$coefficients["TB_index", ]
  data.table(Contrast = paste("APTB vs", g), n = nrow(d),
             OR_per_unit = exp(co[1]), CI_low = exp(co[1] - 1.96 * co[2]),
             CI_high = exp(co[1] + 1.96 * co[2]), P_value = co[4])
}), fill = TRUE)

fwrite(out, file.path(RESULTS, "demographic_adjustment.tsv"), sep = "\t")
message("Adjusted contrasts: ", nrow(out))
