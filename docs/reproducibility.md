# Reproducibility

## Environment

- R >= 4.4.3 (analysis performed with R 4.4.3).
- Package versions as actually used are recorded in `env/package_versions.tsv`, with the
  full `env/sessionInfo.txt` captured from the analysis environment. Versions were
  extracted from the real environment; none were entered by hand.
- No absolute paths are used anywhere in `code/`. The repository root is resolved at run
  time by `code/_paths.R` and can be overridden with the environment variable
  `TBCIRCRNA_ROOT`.

## Execution order

Run every script from the repository root. All 15 steps are listed below in the order in
which they must be run; each entry names the script, its inputs, its outputs and the
corresponding manuscript section.

| # | Step | Script | Input | Output | Manuscript section |
|---|------|--------|-------|--------|--------------------|
| 1 | Metadata preparation | `code/01_prepare_metadata.R` | `data/metadata/{IDC,ITC,PVC}_metadata.tsv` | `metadata/sample_dictionary.tsv`, `metadata/cohort_summary.tsv` | Methods — Study cohorts and sample annotation |
| 2 | Common feature matrix | `code/02_build_common_feature_matrix.R` | `data/qn/raw/{IDC,ITC,PVC}_genelevel.txt` (LONG: gene_name / type / BioSample / TPM); `data/metadata/*_metadata.tsv` | `results/common_feature_genes.tsv` (2,754 genes); `results/feature_filter_counts.tsv`; `data/qn/{IDC,ITC,PVC}_common_2754.tsv` (genes x samples) | Methods — circRNA quantification; common feature space |
| 3 | Frozen IDC-reference normalization | `code/03_frozen_IDC_reference_normalization.R` | `data/qn/{cohort}_common_2754.tsv` (output of step 2); `data/metadata/*_metadata.tsv` | `data/qn/QN_reference_IDC.tsv` (rank, qn_rank_mean); `data/qn/{cohort}_QN_matrix.tsv` (genes x samples) | Methods — IDC-reference quantile normalization |
| 4 | Differential expression | `code/04_differential_expression.R` | `data/qn/IDC_QN_matrix.tsv`; `data/metadata/IDC_metadata.tsv` | `results/DE_complete.tsv` (logFC, AveExpr, t, P, BH-adjusted P, B per contrast) | Methods — Multi-comparator differential expression (each contrast restricted to its two relevant groups) |
| 5 | Four-way candidate intersection | `code/05_four_way_candidate_intersection.R` | `results/DE_complete.tsv` | `results/core18_candidates.tsv` (18 genes) | Methods — Candidate definition |
| 6 | Greedy forward selection | `code/06_greedy_forward_selection.R` | `data/qn/IDC_QN_matrix.tsv`, `results/core18_candidates.tsv` | `results/greedy_selection_path.tsv`, `config/locked_panel.tsv` | Methods — Greedy forward selection and panel locking |
| 7 | Bootstrap stability | `code/07_bootstrap_stability.R` | `data/qn/IDC_QN_matrix.tsv`; `data/metadata/IDC_metadata.tsv`; `config/locked_panel.tsv` | `results/bootstrap_panel_size.tsv`, `results/bootstrap_selection_frequency.tsv`, `results/bootstrap_auc_summary.tsv`, `results/bootstrap_summary.tsv` | Methods — Bootstrap stability analysis; Figure 2E–2G |
| 8 | TB-index scoring | `code/08_TB_index_scoring.R` | `data/qn/{cohort}_QN_matrix.tsv`, `config/locked_panel.tsv` | `results/locked_TB_index_scores.tsv` | Methods — Construction of the TB-index |
| 9 | Diagnostic performance | `code/09_diagnostic_performance.R` | `results/locked_TB_index_scores.tsv`; `config/locked_thresholds.tsv` | `results/diagnostic_performance.tsv` (AUROC, DeLong CI, AUPRC); `results/diagnostic_thresholds.tsv` (Sens, Spec, PPV, NPV, TP, FN, TN, FP at both locked thresholds) | Methods — Diagnostic performance; Table 2; Figure 3 |
| 10 | Demographic adjustment | `code/10_demographic_adjustment.R` | `results/locked_TB_index_scores.tsv`, `metadata/sample_dictionary.tsv` | `results/demographic_adjustment.tsv` | Methods — Age- and sex-adjusted analyses; Supplementary Table S5 |
| 11 | Longitudinal analysis | `code/11_longitudinal_analysis.R` | `results/locked_TB_index_scores.tsv` | `results/longitudinal_analysis.tsv` | Methods — Longitudinal analyses; Supplementary Table S6B |
| 12 | Benchmarking against Sweeney3 and RISK6 | `code/12a_benchmark_tbindex_data.R`, `code/12b_benchmark_tbindex_matrix.R`, `code/12c_benchmark_tbindex_scores.R`, `code/12d_benchmark_risk6_genes.R`, `code/12e_benchmark_literature_data.R`, `code/12f_benchmark_add_genes.R` (run in this order) | `data/benchmark/circ_qn_*.rds`, `data/benchmark/linear_qn_*.rds`, `data/benchmark/{cohort}_TBindex_per_sample.csv` | per-sample benchmark scores; Supplementary Table S7A | Methods — Benchmarking against published signatures |
| 13 | Paired DeLong comparisons | `code/13_DeLong_comparisons.R` | `data/benchmark/SF6_Sweeney3_RISK6_per_sample.csv` | `results/DeLong_comparisons.tsv` (16 comparisons) | Methods — Statistical analysis; Supplementary Table S7D |
| 14 | BSJ refinement | `code/14_BSJ_refinement.R` | `data/bsj/{IDC,ITC,PVC}_GOLPH3_SOD2_PCNX2_circular.tsv` (circular isoforms of the three panel genes, restricted to the final cohort libraries: IDC 109 / ITC 35 / PVC 81) | `results/BSJ_overlap_summary.tsv`, `results/BSJ_contribution_detail.tsv`, `results/BSJ_summary.tsv` | Methods — BSJ-level refinement; Supplementary Table S8; Figure 4 |
| 15 | Summary tables and final outputs | `code/15_generate_summary_tables.R` | `config/*.tsv`, `results/*.tsv` | `results/summary_consistency_check.tsv` | internal consistency check |

`code/07_bootstrap_stability.R` accepts two optional command-line arguments (iterations,
seed) and defaults to the manuscript values 500 and 2026; see its header for usage.

## Input format contract

The pre-filter quantification is **LONG** — one row per host gene x library with
`gene_name`, `type`, `BioSample` and `TPM`. Everything downstream of step 2 is
**WIDE**: features in rows, samples in columns, with the first column holding the gene
symbol (`gene_name` / `Gene`).

- Script **02** is the only place that reads the raw matrices. It performs the
  cohort-specific detectability filter and the three-cohort intersection, and
  writes one common-gene matrix per cohort.
- Script **03** consumes those matrices and never re-reads the raw files, so the
  wide format is assumed in exactly one place.
- Scripts **04, 06 and 08** reorder every matrix explicitly onto the sample
  order of the metadata after checking for duplicate, missing and extra sample
  IDs; no implicit ordering is assumed.

## What can be re-run from this repository alone

| Script | Re-runnable from this repository? | Required external input |
|--------|-----------------------------------|-------------------------|
| `01_prepare_metadata.R` | Yes, with Supplementary Table S1 saved as `data/metadata/*_metadata.tsv` | Supplementary Table S1 |
| `02_build_common_feature_matrix.R` | No | Pre-filter gene-level circRNA host-gene table per cohort (long format) |
| `03_frozen_IDC_reference_normalization.R` | No | Gene-level circular host-gene TPM (GEO GSE348100) |
| `04_differential_expression.R` | No | IDC QN matrix (step 03 output) |
| `05_four_way_candidate_intersection.R` | No | `results/DE_complete.tsv` (step 04 output) |
| `06_greedy_forward_selection.R` | No | IDC QN matrix, `results/core18_candidates.tsv` |
| `07_bootstrap_stability.R` | No | IDC QN matrix and metadata (GEO GSE348100 / Supplementary Table S1) |
| `08_TB_index_scoring.R` | No | Per-cohort QN matrices (GEO GSE348100) |
| `09_diagnostic_performance.R` | Yes, using the shipped `results/locked_TB_index_scores.tsv` | — |
| `10_demographic_adjustment.R` | Yes, using the shipped score table and `metadata/sample_dictionary.tsv` | — |
| `11_longitudinal_analysis.R` | Yes, using the shipped score table | — |
| `12a`–`12f_benchmark_*.R` | No | circRNA and linear-RNA QN matrices |
| `13_DeLong_comparisons.R` | No | Benchmark per-sample scores |
| `14_BSJ_refinement.R` | No | Circular isoform tables for the three panel genes, per cohort |
| `15_generate_summary_tables.R` | Yes | — |

When an input is missing, the script stops with:

```
Required input file not found.
  File : data/qn/IDC_QN_matrix.tsv
  What : IDC QN log2(TPM+1) host-gene matrix
  Where: NCBI GEO GSE348100
```

## What is fixed rather than re-computed

The reported numbers are frozen in `config/` and `results/` and match the manuscript and
Supplementary Tables S1–S8:

- cohort sizes and library counts (`metadata/cohort_summary.tsv`, `metadata/sample_dictionary.tsv`);
- common feature space, candidate count, panel size, bootstrap settings and locked
  thresholds (`config/*.tsv`);
- the 18 direction-concordant candidates with per-contrast log2FC and BH-adjusted
  *P* values (`results/core18_candidates.tsv`);
- per-library TB-index scores (`results/locked_TB_index_scores.tsv`);
- AUROC with 95% CI for every reported comparison (`results/diagnostic_performance.tsv`);
- bootstrap stability summary (`results/bootstrap_summary.tsv`);
- all 16 paired DeLong comparisons with unadjusted and BH-adjusted *P* values
  (`results/DeLong_comparisons.tsv`);
- BSJ overlap counts and largest-shared-BSJ contributions (`results/BSJ_summary.tsv`,
  `results/BSJ_overlap_summary.tsv`, `results/BSJ_contribution_detail.tsv`), all produced
  by `code/14_BSJ_refinement.R` and checked against Figure 4B inside the script.

Selected analysis steps were re-run for reproducibility verification using the
documented processed inputs; the shipped result tables correspond to the
manuscript-reported final results. Not every upstream data-processing step is rerun
automatically from raw sequencing data, because in-house raw sequencing files are not
publicly released.

## Interpretation notes

- Bootstrap out-of-bag AUROC is a resampling summary and is **not** an independent
  validation estimate.
- The PVC PTB-versus-LTBI and EPTB-versus-LTBI comparisons are cross-time-point
  (baseline cases versus post-conversion follow-up controls) and are reported as
  exploratory secondary analyses.
- BH-adjusted *P* values are a secondary multiplicity sensitivity summary; they do not
  redefine the prespecified discovery screen.
