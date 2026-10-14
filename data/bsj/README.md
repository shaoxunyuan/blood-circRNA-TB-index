# data/bsj/

**Required files**

| File | Content |
|------|---------|
| `IDC_GOLPH3_SOD2_PCNX2_circular.tsv`, `ITC_…`, `PVC_…` | Isoform-level quant table filtered to circular isoforms of GOLPH3, SOD2 and PCNX2 (one row per isoform × library) |

**Where to obtain it**

- Derived from isoform-level circular quantification for the three locked panel genes.
- If the isoform-level tables are not available, `code/14_BSJ_refinement.R` cannot be
  re-run, but `results/BSJ_summary.tsv` already contains the reported values.

**How each file is built**

Take the circular isoforms of GOLPH3, SOD2 and PCNX2 from the isoform-level
quantification (`type == "circular"`).

The source BSJ tables may contain additional libraries. `code/14_BSJ_refinement.R`
restricts each table to the final `Sample_ID` list in the corresponding cohort metadata
(`data/metadata/{IDC,ITC,PVC}_metadata.tsv`) and reports any excluded libraries.

The final BSJ analysis uses **109 IDC, 35 ITC and 81 PVC RNA-seq libraries**. Cohort
membership is taken from the final `Sample_ID` list in the cohort metadata; the 81 PVC
libraries already include the longitudinal time points of the final cohort.

A BSJ is identified by `BSJ_ID` (`chromosome|donor|acceptor|strand`, hg38) and is matched
across cohorts by these genomic coordinates.

**Main columns**

`Name`, `TPM`, `NumReads`, `BioSample`, `type`, `gene_id`, `gene_name`, `BSJ_ID`,
`chr`, `circRNAStart`, `circRNAEnd`, `strand`

**Read by**

- `code/14_BSJ_refinement.R` (also reads `data/metadata/{IDC,ITC,PVC}_metadata.tsv` to
  restrict each table to the final cohort `Sample_ID` list)
