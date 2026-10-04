# A three-gene peripheral blood circRNA-derived TB-index for tuberculosis

## Overview

This repository provides the code and documented inputs required to reproduce the
principal downstream analyses reported in the accompanying manuscript for the
circRNA-derived three-gene TB-index (SOD2, GOLPH3, PCNX2). Cohort-specific pre-filter
matrices used to derive the common 2,754-gene feature space are documented separately in
the data-availability information. It is intended for
reproducing the analysis, not for restating the manuscript results.

## Study cohorts

| Cohort | Participants |
|--------|--------------|
| IDC (in-house discovery cohort) | 109 |
| ITC (in-house test cohort) | 35 |
| PVC (public pediatric validation cohort) | 42 |
| **Total** | **186 participants, 225 RNA-seq libraries** |

IDC and ITC comprise adult/adolescent participants; the PVC comprises pediatric
participants (public validation dataset, NCBI BioProject PRJNA588242).

## Analytical framework

| Item | Value |
|------|-------|
| Common feature space | 2,754 circRNA-derived host genes |
| Direction-concordant candidates | 18 genes |
| Locked panel | SOD2, GOLPH3, PCNX2 |
| TB-index | `SOD2_QN − GOLPH3_QN − PCNX2_QN` |
| Locked Youden threshold | 1.407 |
| Locked high-sensitivity threshold | 0.261 |

The prespecified discovery screen was |log2FC| > 0.5 **and** nominal *P* < 0.05 in all
four APTB-centred comparisons, combined with four-way direction concordance.
Benjamini–Hochberg-adjusted *P* values are reported as secondary multiplicity
sensitivity information and do not redefine the screen.

## Repository structure

| Path | Contents |
|------|----------|
| `code/` | Numbered analysis scripts (01–15) matching the manuscript Methods; each starts with Purpose / Inputs / Outputs / Key parameters / Dependencies, plus the shared `palette_manuscript.R` and `_paths.R`. |
| `config/` | `analysis_parameters.tsv`, `locked_panel.tsv`, `locked_thresholds.tsv` — the single source of truth for the final parameters. |
| `data/` | Input placeholders with per-directory READMEs: where each deposited input file must be placed, its source, its columns, and which scripts read it. |
| `metadata/` | `sample_dictionary.tsv` (de-identified per-library annotation, 225 libraries) and `cohort_summary.tsv`. |
| `results/` | Lightweight derived results: common feature gene list, candidate genes, locked TB-index scores, diagnostic performance (threshold-free metrics plus operating characteristics at the two locked thresholds), bootstrap panel-size / selection-frequency / AUROC summaries, paired DeLong comparisons, BSJ overlap / contribution / summary tables. |
| `env/` | `sessionInfo.txt` and `package_versions.tsv` captured from the analysis environment. |
| `docs/` | `analysis_workflow.md`, `data_availability.md`, `reproducibility.md`. |
| `output/` | Default write directory when scripts are re-run. |

## Reproducibility

Run from the repository root; every script sources `code/_paths.R`, which resolves the
repository root without any absolute path.

```r
# 01  de-identified sample dictionary and cohort summary
Rscript code/01_prepare_metadata.R
# 02  common feature space (2,754 host genes)
Rscript code/02_build_common_feature_matrix.R
# 03  frozen IDC-reference quantile normalization
Rscript code/03_frozen_IDC_reference_normalization.R
# 04  four APTB-centred limma comparisons
Rscript code/04_differential_expression.R
# 05  18 four-way direction-concordant candidates
Rscript code/05_four_way_candidate_intersection.R
# 06  greedy forward selection -> locked panel
Rscript code/06_greedy_forward_selection.R
# 07  bootstrap stability
Rscript code/07_bootstrap_stability.R
# 08  TB-index scoring
Rscript code/08_TB_index_scoring.R
# 09  diagnostic performance
Rscript code/09_diagnostic_performance.R
# 10  age/sex-adjusted sensitivity analyses
Rscript code/10_demographic_adjustment.R
# 11  longitudinal analysis (PVC)
Rscript code/11_longitudinal_analysis.R
# 12a-12f  Sweeney3 / RISK6 benchmark
Rscript code/12a_benchmark_tbindex_data.R
Rscript code/12b_benchmark_tbindex_matrix.R
Rscript code/12c_benchmark_tbindex_scores.R
Rscript code/12d_benchmark_risk6_genes.R
Rscript code/12e_benchmark_literature_data.R
Rscript code/12f_benchmark_add_genes.R
# 13  paired DeLong comparisons
Rscript code/13_DeLong_comparisons.R
# 14  BSJ refinement
Rscript code/14_BSJ_refinement.R
# 15  summary tables and consistency check
Rscript code/15_generate_summary_tables.R
```

The benchmark steps require the deposited processed benchmark inputs documented in
`data/benchmark/README.md`.

Details, including which steps require the deposited matrices and what each script
expects, are in `docs/analysis_workflow.md` and `docs/reproducibility.md`.

## Data availability

In-house raw sequencing data are not publicly released because the original consent did
not authorize public release of individual-level sequencing data.

Processed circRNA-derived expression matrices have been deposited in NCBI GEO under
accession **GSE348100** and will be released upon publication.

The public validation RNA-seq dataset is available from NCBI BioProject **PRJNA588242**.

This repository contains analysis code, configuration files, lightweight derived
results, and reproducibility documentation.

## Software environment

See `env/sessionInfo.txt` and `env/package_versions.tsv` (captured from the real
analysis environment).

## Citation

See `CITATION.cff`.

## License

Code is released under the MIT License (`LICENSE`).
