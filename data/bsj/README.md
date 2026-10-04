# data/bsj/

**Required files**

| File | Content |
|------|---------|
| `IDC_GOLPH3_SOD2_PCNX2_circular.tsv`, `ITC_…`, `PVC_…` | Isoform-level quant table filtered to circular isoforms of GOLPH3, SOD2 and PCNX2 (one row per isoform × library) |

**Where to obtain it**

- Derived from the isoform-level quantification used in the manuscript; the
  corresponding summary tables are published as **Supplementary Table S8**.
- If the isoform-level tables are not available, `code/14_BSJ_refinement.R` cannot be
  re-run, but `results/BSJ_summary.tsv` already contains the reported values.

**How each file must be built**

Take the circular isoforms of GOLPH3, SOD2 and PCNX2 from the isoform-level
quantification (`type == "circular"`), and keep **only the libraries of the final
cohort**:

| File | Final cohort libraries |
|------|------------------------|
| `IDC_…_circular.tsv` | 109 |
| `ITC_…_circular.tsv` | 35 |
| `PVC_…_circular.tsv` | 81 |

The cohort restriction is part of the analysis definition, not a formatting detail:
the source tables contain libraries that are not part of the final cohort (the ITC
pre-filter table contains 36 libraries, and the PVC isoform table contains follow-up and
post-treatment libraries). Leaving them in changes the junction sets and the
overlap/contribution values, so a BSJ that is detected only in an excluded library must
not enter the analysis. A BSJ is identified by `BSJ_ID`
(`chromosome|donor|acceptor|strand`, hg38) and is matched across cohorts by these
genomic coordinates only — never by row order or row number.

**Main columns**

`Name`, `TPM`, `NumReads`, `BioSample`, `type`, `gene_id`, `gene_name`, `BSJ_ID`,
`chr`, `circRNAStart`, `circRNAEnd`, `strand`

**Read by**

- `code/14_BSJ_refinement.R`
