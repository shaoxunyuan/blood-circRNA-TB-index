#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Purpose:      Build the de-identified per-library sample dictionary and the
#               cohort/group participant summary used by every downstream step.
# Inputs:       data/metadata/{IDC,ITC,PVC}_metadata.tsv  (the de-identified cohort sample metadata)
# Outputs:      metadata/sample_dictionary.tsv
#               metadata/cohort_summary.tsv
# Key params:   No filtering is applied here; inclusion flags are carried through
#               (Included_in_baseline_diagnostic_analysis,
#                Included_in_longitudinal_analysis).
# Dependencies: R >= 4.4.3; data.table
# ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(data.table))
source("code/_paths.R")

cohorts <- c("IDC", "ITC", "PVC")
parts <- lapply(cohorts, function(coh) {
  f <- file.path(DATA_META, paste0(coh, "_metadata.tsv"))
  paths_require(f, paste0(coh, " per-library sample annotation"),
                "the de-identified cohort sample metadata / NCBI GEO GSE348100")
  d <- fread(f, na.strings = c("NA", "", "N.A"))
  d[, Cohort := coh]
  d
})
samples <- rbindlist(parts, fill = TRUE)

if (!"Participant_ID" %in% names(samples)) samples[, Participant_ID := Sample_ID]

setcolorder(samples, intersect(c("Sample_ID", "Participant_ID", "Cohort", "Group",
                                 "TimePoint", "Sex", "Age",
                                 "Baseline_diagnostic_analysis",
                                 "Longitudinal_analysis"), names(samples)))
fwrite(samples, file.path(ROOT, "metadata", "sample_dictionary.tsv"), sep = "\t")

cohort_summary <- samples[, .(Participants = uniqueN(Participant_ID)), by = .(Cohort, Group)]
setorder(cohort_summary, Cohort, Group)
fwrite(cohort_summary, file.path(ROOT, "metadata", "cohort_summary.tsv"), sep = "\t")

message("Libraries: ", nrow(samples),
        " | participants: ", uniqueN(samples$Participant_ID),
        " | per cohort: ", paste(cohorts, sapply(cohorts, function(c)
          uniqueN(samples[Cohort == c, Participant_ID])), sep = "=", collapse = " "))
