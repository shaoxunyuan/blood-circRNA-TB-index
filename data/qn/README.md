# data/qn/

**Required files**

| File | Content |
|------|---------|
| `IDC_QN_matrix.tsv`, `ITC_QN_matrix.tsv`, `PVC_QN_matrix.tsv` | Quantile-normalised log2(TPM+1) circRNA host-gene expression, 2,754 genes × libraries |
| `IDC_QN_matrix_before.tsv`, … | The same matrices before normalisation (log2(TPM+1)) |
| `QN_reference_IDC.tsv` | Frozen IDC reference: `rank`, `qn_rank_mean` |

**Where to obtain it**

- Processed expression matrices and associated sample metadata have been deposited in
  NCBI GEO under accession **GSE348100**. The GEO record is currently private and is
  planned to be made publicly available upon manuscript acceptance.
- Place the matrices in this directory using the filenames listed above before re-running
  the scripts.

**Matrix orientation**

Each matrix is stored with **one row per host gene**. The first column is `Gene`,
followed by one column for each RNA-seq library:

| File | Dimensions |
|------|------------|
| `IDC_QN_matrix.tsv` | 2,754 genes × 109 libraries |
| `ITC_QN_matrix.tsv` | 2,754 genes × 35 libraries |
| `PVC_QN_matrix.tsv` | 2,754 genes × 81 libraries |

`QN_reference_IDC.tsv` has the columns `rank` and `qn_rank_mean` (2,754 rows).

**Read by**

- `code/03_frozen_IDC_reference_normalization.R`
- `code/04_differential_expression.R`
- `code/06_greedy_forward_selection.R`
- `code/07_bootstrap_stability.R`
- `code/08_TB_index_scoring.R`
