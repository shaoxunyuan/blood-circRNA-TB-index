# code/

Analysis scripts numbered in execution order, matching the manuscript Methods.
Every script starts with a standard header:

```
Purpose:      what the script does
Inputs:       what it reads and where the inputs come from
Outputs:      what it writes
Key parameters: parameters that fix the result
Dependencies: R version and packages
```

## Scripts

| Script | Step | Role |
|--------|------|------|
| `_paths.R` | — | Resolves repository-relative paths (`ROOT`, `CONFIG`, `DATA*`, `RESULTS`, …) and provides `paths_require()` for missing-input errors. Sourced by every script. No absolute paths. |
| `palette_manuscript.R` | — | Single source of truth for the locked final colours and shapes (group / cohort / feature / marker / direction). |
| `01_prepare_metadata.R` | 1 | De-identified sample dictionary and cohort/group summary |
| `02_build_common_feature_matrix.R` | 2 | Common feature space (2,754 circRNA host genes) |
| `03_frozen_IDC_reference_normalization.R` | 3 | Frozen IDC-reference quantile normalization |
| `04_differential_expression.R` | 4 | Four APTB-centred limma comparisons with BH-adjusted P values |
| `05_four_way_candidate_intersection.R` | 5 | 18 four-way direction-concordant candidates |
| `06_greedy_forward_selection.R` | 6 | Greedy forward selection; locks the panel |
| `07_bootstrap_stability.R` | 7 | 500× stratified bootstrap stability |
| `08_TB_index_scoring.R` | 8 | TB-index = SOD2_QN − GOLPH3_QN − PCNX2_QN |
| `09_diagnostic_performance.R` | 9 | AUROC with 95% CI and locked-threshold operating points |
| `10_demographic_adjustment.R` | 10 | Age- and sex-adjusted sensitivity analyses |
| `11_longitudinal_analysis.R` | 11 | PVC longitudinal analysis (exploratory) |
| `12a`–`12f_benchmark_*.R` | 12 | Sweeney3 / RISK6 benchmark pipeline (run in alphabetical order) |
| `13_DeLong_comparisons.R` | 13 | 16 paired DeLong comparisons; one BH family |
| `14_BSJ_refinement.R` | 14 | BSJ contribution analysis for the three panel genes |
| `15_generate_summary_tables.R` | 15 | Assembles and cross-checks the summary tables |

## Usage

Run from the repository root:

```bash
Rscript code/01_prepare_metadata.R
```

`code/_paths.R` resolves the repository root from the working directory (or the script
location), and can be overridden with `TBCIRCRNA_ROOT`.

## Inputs

Input matrices have been deposited in NCBI GEO (accession **GSE348100**). The GEO record
is currently private and is planned to be made publicly available upon manuscript
acceptance; they are not redistributed here. `data/` documents exactly where each file
must be placed; see `data/README.md` and `docs/reproducibility.md`.

When an input is missing, scripts stop with an explicit message:

```
Required input file not found.
  File : data/qn/IDC_QN_matrix.tsv
  What : IDC QN log2(TPM+1) host-gene matrix
  Where: NCBI GEO GSE348100
```
