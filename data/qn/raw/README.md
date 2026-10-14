# data/qn/raw/

> **Note on the directory name.** Here, "raw" means pre-normalization gene-level
> processed matrices used as analysis inputs (circular host-gene TPM tables).

**Required input files**

| File | Content |
|------|---------|
| `IDC_genelevel.txt`, `ITC_genelevel.txt`, `PVC_genelevel.txt` | Pre-filter gene-level circRNA host-gene quantification, **long** format: `gene_id`, `gene_name`, `gene_biotype`, `type`, `BioSample`, `TPM`, `NumReads`, `n_isoforms` (one row per host gene x library; only `type == "circular"` rows are used) |

**Where to obtain it**

- These processed gene-level matrices accompany the deposited data: processed expression
  matrices and associated sample metadata have been deposited in NCBI GEO under accession
  **GSE348100**. The GEO record is currently private and is planned to be made publicly
  available upon manuscript acceptance.
- Raw sequencing data (FASTQ) from the in-house cohorts are not publicly available because
  the original informed consent did not authorize public release of individual-level
  sequencing data.

**Main columns**

`gene_name`, `type` (`circular` / `linear`), `BioSample`, `TPM`. The library labels in
`BioSample` must match `Sample_ID` in the corresponding metadata file.

**Cohort completeness**

Each pre-filter table must contain all libraries required for the final cohort analysis.
Libraries present in the pre-filter table but absent from the final metadata are excluded
explicitly by sample ID and reported by the processing script; a metadata library missing
from the table stops the run with an explicit message. The final cohort sizes are
IDC 109, ITC 35 and PVC 81 libraries.

**Read by**

- `code/02_build_common_feature_matrix.R` (reads `{cohort}_genelevel.txt`)
- `code/03_frozen_IDC_reference_normalization.R`

These files are large and are excluded by `.gitignore`; keep them local only.
