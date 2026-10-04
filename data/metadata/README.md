# data/metadata/

**Required file**

| File | Content |
|------|---------|
| `IDC_metadata.tsv`, `ITC_metadata.tsv`, `PVC_metadata.tsv` | Per-library sample annotation used as analysis input: sample identifier, cohort, group, time point, inclusion flags, sex, age |

Every script that needs cohort or sample information reads these three files; no
second metadata input format is used anywhere in the repository.

**Where to obtain it**

- Same content as **Supplementary Table S1** of the manuscript.
- A de-identified, repository-ready version is already shipped at
  `metadata/sample_dictionary.tsv` (all 225 libraries). Split it by cohort if the
  scripts are run per cohort.

**Main columns**

`Sample_ID`, `Participant_ID`, `Cohort`, `Group`, `TimePoint`, `Sex`, `Age`,
`Baseline_diagnostic_analysis`, `Longitudinal_analysis`

**Read by**

- `code/01_prepare_metadata.R`
- `code/03_frozen_IDC_reference_normalization.R`
- `code/09_diagnostic_performance.R`
- `code/10_demographic_adjustment.R`
- `code/11_longitudinal_analysis.R`
