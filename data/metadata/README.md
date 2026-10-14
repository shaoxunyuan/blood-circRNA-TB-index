# data/metadata/

**Required file**

| File | Content |
|------|---------|
| `IDC_metadata.tsv`, `ITC_metadata.tsv`, `PVC_metadata.tsv` | Per-library sample annotation used as analysis input: sample identifier, cohort, group, time point, inclusion flags, sex, age |

Every script that needs cohort or sample information reads these three files; no
second metadata input format is used anywhere in the repository.

**Where to obtain it**

- A de-identified, repository-ready version is shipped at
  `metadata/sample_dictionary.tsv` (all 225 libraries). Split that file by `Cohort` to
  obtain the three analysis input tables listed above, or download the matching
  de-identified sample metadata from NCBI GEO **GSE348100**.

  `Sample_ID` must be the sample identifier used by the expression and quantification
  inputs of the same cohort (`data/qn/...` matrices, `data/bsj/...` tables); the scripts
  join on that string. When the deposited matrices use a different identifier set from
  the published sample dictionary, map the identifiers before splitting the dictionary
  by cohort — do not edit the identifier inside a script.

**Main columns**

`Sample_ID`, `Participant_ID`, `Cohort`, `Group`, `TimePoint`, `Sex`, `Age`,
`Baseline_diagnostic_analysis`, `Longitudinal_analysis`

**Read by**

These cohort metadata files are used throughout the analysis pipeline wherever sample
membership, cohort, clinical group, age, sex, inclusion status or time point information
is required: `code/01`, `code/02`, `code/03`, `code/04`, `code/06`, `code/07`,
`code/08`, `code/12a` and `code/14`.

(`code/09`, `code/10` and `code/11` read the shipped
`metadata/sample_dictionary.tsv` instead, because they work from the already-scored
`results/locked_TB_index_scores.tsv`.)
