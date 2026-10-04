# data/benchmark/

**Required files**

| File | Content |
|------|---------|
| `circ_qn_{IDC,ITC,PVC}.rds` | circRNA QN matrices used for the benchmark |
| `linear_qn_{IDC,ITC,PVC}.rds` | Linear-RNA QN matrices (13,989 genes) used for Sweeney3 / RISK6 |
| `{IDC,ITC,PVC}_TBindex_per_sample.csv` | Per-sample TB-index used in the benchmark |
| `benchmark_literature.xlsx` | Curated literature comparison (**Supplementary Table S7A**) |

**Where to obtain it**

- circRNA matrices: NCBI GEO **GSE348100**.
- Linear-RNA matrices: derived from the same RNA-seq quantification; the public
  validation cohort is available from NCBI BioProject **PRJNA588242**.
- Literature comparison: **Supplementary Table S7A** of the manuscript.

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

- `code/13_DeLong_comparisons.R`

The reported benchmark results are shipped as `results/DeLong_comparisons.tsv`.
