# Benchmark pipeline

Data flow of the Sweeney3 / RISK6 benchmark (steps 12a–12f and 13). Every path is
relative to the repository root and every input is declared in
`data/benchmark/README.md` or `data/metadata/README.md`.

| Script | Input | Output | Next |
|--------|-------|--------|------|
| `12a_benchmark_tbindex_data.R` | `data/benchmark/{IDC,ITC,PVC}_TBindex_per_sample.csv`; `data/qn/{cohort}_QN_matrix.tsv` (linear-RNA QN values); `data/metadata/{IDC,ITC,PVC}_metadata.tsv` | per-sample benchmark scores: `Sample_ID, Cohort, Group, TimePoint, circRNA_TB_index, Sweeney3, RISK6` (written to the benchmark workbook under `output/`) | 12b, 12c, 12f |
| `12b_benchmark_tbindex_matrix.R` | `data/benchmark/circ_qn_{IDC,ITC,PVC}.rds`; `data/benchmark/linear_qn_{IDC,ITC,PVC}.rds`; the 12a per-sample table | expression matrix sheet: 3 circRNA genes + 9 linear-RNA genes, samples aligned to 12a | 12c |
| `12c_benchmark_tbindex_scores.R` | the 12b expression matrix | per-sample Sweeney3 / RISK6 scores as signed weighted sums of QN values | 12f |
| `12d_benchmark_risk6_genes.R` | `data/benchmark/linear_qn_{IDC,ITC,PVC}.rds` | per-gene expression and ATB differential-expression sheets for the Sweeney3 and RISK6 gene sets | 12e, 12f |
| `12e_benchmark_literature_data.R` | curated literature summary (shipped as Supplementary Table S7A) | literature comparison sheet | — |
| `12f_benchmark_add_genes.R` | the 12a per-sample table; linear-RNA QN values for GBP5, DUSP3, KLF2 | per-sample table extended with the Sweeney3 constituent genes | 13 |
| `13_DeLong_comparisons.R` | `data/benchmark/SF6_Sweeney3_RISK6_per_sample.csv` | `results/DeLong_comparisons.tsv` (16 paired comparisons, one BH family) | — |

## Inputs and dependencies

Every benchmark input is either shipped in `results/` or documented in
`data/*/README.md`; no undeclared input is read. `12a` reads the same per-cohort metadata
TSVs that steps 2–11 use, so no second metadata file is needed.

## Metadata fields used by the benchmark

`Sample_ID` (matches the TB-index per-sample tables), `Group` (contrast definition) and
`TimePoint` (PVC baseline selection). These come from the same three TSVs that steps 2–11
use; no second metadata file is maintained.
