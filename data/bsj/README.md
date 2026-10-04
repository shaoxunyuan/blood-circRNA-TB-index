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

The cohort restriction is part of the analysis definition, not a formatting detail.
Source isoform-level tables may contain libraries outside the final analysis set, so the
BSJ input files must be restricted by the final sample identifiers documented in the
cohort metadata. The final BSJ analysis uses **109 IDC, 35 ITC and 81 PVC RNA-seq
libraries**.

For each cohort, retain only the `BioSample` identifiers present in the final cohort
metadata used by the manuscript analysis. Do not select libraries by time point, baseline
/ follow-up status or treatment status: the 81 PVC libraries, for example, already
include the longitudinal time points of the final cohort, and the only reliable
selection key is the final sample identifier list.

A BSJ is identified by `BSJ_ID` (`chromosome|donor|acceptor|strand`, hg38) and is matched
across cohorts by these genomic coordinates only — never by row order or row number.

**Main columns**

`Name`, `TPM`, `NumReads`, `BioSample`, `type`, `gene_id`, `gene_name`, `BSJ_ID`,
`chr`, `circRNAStart`, `circRNAEnd`, `strand`

**Read by**

- `code/14_BSJ_refinement.R`
