# Analysis workflow

Fifteen analysis steps matching the manuscript Methods. Scripts are numbered in
execution order; every script sources `code/_paths.R`, and all filenames below are the
ones the scripts actually use.

---

### Step 1 — Sample annotation

- **Script:** `code/01_prepare_metadata.R`
- **Input:** `data/metadata/{IDC,ITC,PVC}_metadata.tsv`
- **Output:** `metadata/sample_dictionary.tsv`, `metadata/cohort_summary.tsv`
- **Purpose:** Build the de-identified per-library dictionary and the cohort/group
  participant summary used by every downstream step.
- **Manuscript:** Methods — Study cohorts and sample annotation

### Step 2 — Common feature matrix

- **Script:** `code/02_build_common_feature_matrix.R`
- **Input:** `data/qn/raw/{IDC,ITC,PVC}_genelevel.txt` (long: gene_name / type / BioSample / TPM), `data/metadata/{cohort}_metadata.tsv`
- **Output:** `data/qn/{IDC,ITC,PVC}_common_2754.tsv`, `results/common_feature_genes.tsv`, `results/feature_filter_counts.tsv`
- **Purpose:** Apply the cohort-specific detectability filter and take the three-cohort
  intersection of circRNA host genes.
- **Manuscript:** Methods — circRNA quantification; definition of the common feature space

### Step 3 — Frozen IDC-reference normalization

- **Script:** `code/03_frozen_IDC_reference_normalization.R`
- **Input:** `data/qn/{IDC,ITC,PVC}_common_2754.tsv` (Step 2 output), `data/metadata/{cohort}_metadata.tsv`
- **Output:** `data/qn/QN_reference_IDC.tsv`, `data/qn/{IDC,ITC,PVC}_QN_matrix.tsv`
- **Purpose:** Fit the frozen IDC-reference rank-mean distribution on the IDC only and map
  every sample of every cohort onto it.
- **Manuscript:** Methods — IDC-reference quantile normalization

### Step 4 — Differential expression

- **Script:** `code/04_differential_expression.R`
- **Input:** `data/qn/IDC_QN_matrix.tsv`, `data/metadata/IDC_metadata.tsv`
- **Output:** `results/DE_complete.tsv`
- **Purpose:** Four APTB-centred limma comparisons, each fitted on its two relevant
  clinical groups, with nominal and BH-adjusted P values for all 2,754 genes.
- **Manuscript:** Methods — Multi-comparator differential expression and candidate definition

### Step 5 — Four-way direction-concordant candidates

- **Script:** `code/05_four_way_candidate_intersection.R`
- **Input:** `results/DE_complete.tsv`
- **Output:** `results/core18_candidates.tsv`
- **Purpose:** Apply the prespecified screen (|log2FC| > 0.5, nominal P < 0.05) in all four
  contrasts, require direction concordance, and report the per-contrast BH-adjusted P values.
- **Manuscript:** Methods — Candidate definition; Results — 18 direction-concordant genes

### Step 6 — Greedy forward panel selection

- **Script:** `code/06_greedy_forward_selection.R`
- **Input:** `data/qn/IDC_QN_matrix.tsv`, `results/core18_candidates.tsv`, `data/metadata/IDC_metadata.tsv`
- **Output:** `results/greedy_selection_path.tsv`, `config/locked_panel.tsv`
- **Purpose:** Greedy forward selection over the 18 candidates against APTB versus pooled
  non-APTB, stopping when the incremental AUROC is <= 0.02.
- **Manuscript:** Methods — Greedy forward selection and panel locking; Figure 2D

### Step 7 — Bootstrap stability

- **Script:** `code/07_bootstrap_stability.R` (optional arguments: iterations, seed)
- **Input:** `data/qn/IDC_QN_matrix.tsv`, `data/metadata/IDC_metadata.tsv`, `config/locked_panel.tsv`
- **Output:** `results/bootstrap_panel_size.tsv`, `results/bootstrap_selection_frequency.tsv`, `results/bootstrap_auc_summary.tsv`, `results/bootstrap_summary.tsv`
- **Purpose:** Re-run the complete discovery and greedy-selection procedure in 500
  stratified bootstrap resamples and summarise panel size, selection frequency and the
  in-bag / out-of-bag AUROC of the locked panel.
- **Manuscript:** Methods — Bootstrap stability analysis; Figure 2E–2G

### Step 8 — TB-index scoring

- **Script:** `code/08_TB_index_scoring.R`
- **Input:** `data/qn/{IDC,ITC,PVC}_QN_matrix.tsv`, `data/metadata/{cohort}_metadata.tsv`, `config/locked_panel.tsv`
- **Output:** `results/locked_TB_index_scores.tsv`
- **Purpose:** Compute TB-index = SOD2_QN − GOLPH3_QN − PCNX2_QN for all 225 libraries,
  applying the IDC-locked panel unchanged to ITC and PVC.
- **Manuscript:** Methods — Construction of the TB-index

### Step 9 — Diagnostic performance

- **Script:** `code/09_diagnostic_performance.R`
- **Input:** `results/locked_TB_index_scores.tsv`, `config/locked_thresholds.tsv`
- **Output:** `results/diagnostic_performance.tsv`, `results/diagnostic_thresholds.tsv`
- **Purpose:** AUROC with DeLong 95% CI and AUPRC for every reported comparison, plus
  sensitivity, specificity, PPV, NPV and the TP/FN/TN/FP counts at both locked thresholds.
- **Manuscript:** Methods — Diagnostic performance; Table 2; Figure 3

### Step 10 — Age- and sex-adjusted analyses

- **Script:** `code/10_demographic_adjustment.R`
- **Input:** `results/locked_TB_index_scores.tsv`, `metadata/sample_dictionary.tsv`
- **Output:** `results/demographic_adjustment.tsv`
- **Purpose:** Fit the TB-index ~ Age + Sex residual model and the per-contrast
  disease ~ TB-index + Age + Sex logistic models.
- **Manuscript:** Methods — Age- and sex-adjusted analyses; Supplementary Table S5

### Step 11 — Longitudinal analysis

- **Script:** `code/11_longitudinal_analysis.R`
- **Input:** `results/locked_TB_index_scores.tsv`
- **Output:** `results/longitudinal_analysis.tsv`
- **Purpose:** Paired within-participant comparisons of adjacent PVC time points, reported
  as exploratory.
- **Manuscript:** Methods — Longitudinal analyses; Supplementary Table S6B

### Step 12 — Benchmark against Sweeney3 and RISK6

- **Scripts:** `code/12a_benchmark_tbindex_data.R`, `code/12b_benchmark_tbindex_matrix.R`,
  `code/12c_benchmark_tbindex_scores.R`, `code/12d_benchmark_risk6_genes.R`,
  `code/12e_benchmark_literature_data.R`, `code/12f_benchmark_add_genes.R` (run in this order)
- **Input:** `data/benchmark/{circ,linear}_qn_*.rds`, `data/benchmark/{cohort}_TBindex_per_sample.csv`, `data/metadata/{cohort}_metadata.tsv`
- **Output:** per-sample benchmark scores and the literature comparison (Supplementary Table S7A)
- **Purpose:** Apply the published Sweeney3 and RISK6 formulas to the linear-RNA QN values
  without re-estimating coefficients or directions.
- **Manuscript:** Methods — Benchmarking against published signatures

### Step 13 — Paired DeLong comparisons

- **Script:** `code/13_DeLong_comparisons.R`
- **Input:** `data/benchmark/SF6_Sweeney3_RISK6_per_sample.csv` (per-sample TB-index and benchmark signature scores)
- **Output:** `results/DeLong_comparisons.tsv`
- **Purpose:** Paired DeLong tests between the TB-index and each comparator signature for
  all 16 comparisons, with unadjusted and BH-adjusted P values from one family.
- **Manuscript:** Methods — Statistical analysis; Supplementary Table S7D

### Step 14 — BSJ refinement

- **Script:** `code/14_BSJ_refinement.R`
- **Input:** `data/bsj/{IDC,ITC,PVC}_GOLPH3_SOD2_PCNX2_circular.tsv` (BSJ abundance and annotation for the three panel genes; long format with `gene_name`, `BSJ_ID`, `BioSample`, `TPM`) **+** `data/metadata/{IDC,ITC,PVC}_metadata.tsv` (final cohort membership). The script restricts each BSJ table internally to the final `Sample_ID` list of the cohort metadata, reports the libraries it drops, and expects the final cohort sizes IDC 109 / ITC 35 / PVC 81. `data/bsj/README.md` documents how the input tables are built.
- **Output:** `results/BSJ_overlap_summary.tsv`, `results/BSJ_contribution_detail.tsv`, `results/BSJ_summary.tsv`
- **Purpose:** Cross-cohort BSJ overlap and contribution analysis for the locked genes. A BSJ is matched across cohorts by its genomic junction (`BSJ_ID` = chromosome|donor|acceptor|strand), never by row index. For each gene and cohort the analysis reports the all-three-cohort count, the shared and cohort-specific junction counts, and the contribution of the largest shared BSJ, of the other shared BSJs and of the cohort-specific BSJs to the cohort-pooled circular TPM of that gene.
- **Manuscript:** Methods — BSJ-level refinement; Results — BSJ-level refinement; **Figure 4**; **Supplementary Table S8**

**Mapping to the manuscript deliverables**

| Script 14 output | Supplementary Table S8 sheet | Figure 4 |
|------------------|------------------------------|----------|
| `BSJ_overlap_summary.tsv` | S8A (per-gene cohort overlap counts) | **Figure 4A** (cross-cohort BSJ overlap) |
| `BSJ_contribution_detail.tsv` | S8B / S8C / S8D (per-BSJ annotation, pooled TPM and contribution) | **Figure 4B** (largest shared, other shared, cohort-specific) |
| `BSJ_summary.tsv` | S8D summary and S8E–S8G (per-sample BSJ and gene-level detail) | **Figure 4B** |

### Step 15 — Summary tables and consistency check

- **Script:** `code/15_generate_summary_tables.R`
- **Input:** `config/*.tsv`, `results/*.tsv`
- **Output:** `results/summary_consistency_check.tsv`
- **Purpose:** Assemble the lightweight summary tables and verify each shipped result table against the locked parameters. Every check compares an *expected* value from `config/` with an *observed* value read from the corresponding result file (common feature space, candidate genes, locked panel, scored libraries, diagnostic comparisons, DeLong comparisons, BSJ summary rows and the all-three-cohort BSJ counts).
- **Manuscript:** — (internal consistency check; no manuscript section)

---

## Shared conventions

- **Input format:** the pre-filter quantification is **long** (`gene_name`, `type`,
  `BioSample`, `TPM`); everything downstream is **wide** (features in rows, samples in
  columns, first column holding the gene symbol). Script 02 is the only script that reads
  the pre-filter tables; script 03 consumes its output.
- **Sample-ID alignment:** steps 02, 03, 04, 06 and 08 check for duplicated, missing and
  extra sample IDs and then reorder explicitly onto the metadata order.
- **Feature space:** 2,754 circRNA-derived host genes — the intersection of the
  detectability-filtered gene sets of IDC, ITC and PVC (cohort-level retention:
  3,535 / 3,058 / 4,457).
- **Normalization:** log2(TPM + 1) mapped onto a frozen IDC-reference rank-mean
  distribution; no validation-cohort label, group or outcome enters the normalization.
- **Discovery screen:** |log2FC| > 0.5 and nominal *P* < 0.05, four-comparator
  intersection, direction concordance. BH-adjusted *P* values are secondary multiplicity
  sensitivity information.
- **Locked panel:** SOD2 (+1), GOLPH3 (−1), PCNX2 (−1);
  TB-index = SOD2_QN − GOLPH3_QN − PCNX2_QN.
- **Locked thresholds** (IDC APTB vs HC, applied unchanged elsewhere): Youden = 1.407;
  high-sensitivity (10th percentile of the IDC APTB TB-index) = 0.261.

## Working from the shipped files only

`results/` already contains the final derived values for every step, matching the
manuscript and Supplementary Tables S1–S8. Each script stops with an explicit
`Required input file not found` message naming the file and its source when an input is
missing.
