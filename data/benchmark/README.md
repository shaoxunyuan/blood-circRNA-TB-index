# data/benchmark/

**Required files**

| File | Content |
|------|---------|
| `circ_qn_{IDC,ITC,PVC}.rds` | circRNA QN matrices used for the benchmark |
| `linear_qn_{IDC,ITC,PVC}.rds` | Linear-RNA QN matrices (13,989 genes) used for Sweeney3 / RISK6 |
| `{IDC,ITC,PVC}_TBindex_per_sample.csv` | Per-sample TB-index used in the benchmark |
| `benchmark_literature.xlsx` | Curated literature comparison for the signature benchmark |
| `benchmark_Sweeney3_RISK6_per_sample.csv` | Per-sample TB-index and Sweeney3 / RISK6 scores for paired DeLong comparisons |

**Where to obtain it**

- circRNA matrices: NCBI GEO **GSE348100** (currently private; planned public release
  upon manuscript acceptance).
- Linear-RNA matrices: derived from the same RNA-seq quantification; raw sequencing data
  for the public validation cohort are publicly available from NCBI BioProject
  **PRJNA588242**. In-house raw sequencing data are not publicly available.
- Literature comparison workbook: provide as `benchmark_literature.xlsx` when re-running
  the literature-comparison step.
- Per-sample benchmark scores (`benchmark_Sweeney3_RISK6_per_sample.csv`): shipped in this
  directory for the DeLong step. Full expression matrices remain in NCBI GEO
  **GSE348100**.

**Main columns**

`Sample_ID`, `Cohort`, `Group`, `TimePoint`, `TBindex`, and per-gene QN expression columns.

**Read by**

Benchmarking is performed using the scripts

- `code/12a_benchmark_tbindex_data.R`
- `code/12b_benchmark_tbindex_matrix.R`
- `code/12c_benchmark_tbindex_scores.R`
- `code/12d_benchmark_risk6_genes.R`
- `code/12e_benchmark_literature_data.R`
- `code/12f_benchmark_add_genes.R`

run in this order. Paired ROC comparisons are then performed using

- `code/13_DeLong_comparisons.R` (reads `benchmark_Sweeney3_RISK6_per_sample.csv`)

The reported benchmark results are shipped as `results/DeLong_comparisons.tsv`.


## RISK6 note

The RISK6 score in this package follows the published RNA-seq/microarray formula
(Penn-Nicholson et al., 2020):
`GM(GBP2, FCGR1B, SERPING1) − GM(TUBGCP6, TRMT2A, SDR39U1)`,
with geometric means evaluated on the linear (TPM+1) scale from frozen IDC-reference
QN-normalized log2(TPM+1) values (Penn-Nicholson et al., 2020 RNA-seq/microarray
formula). The internal column name remains `RISK6`.
