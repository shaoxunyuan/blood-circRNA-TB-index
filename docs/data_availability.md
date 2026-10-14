# Data availability

## In-house cohorts

Processed expression matrices and associated sample metadata from the in-house discovery
cohort (IDC) and in-house test cohort (ITC) have been deposited in NCBI GEO under
accession **GSE348100**. The GEO record is currently private and is planned to be made
publicly available upon manuscript acceptance. During peer review, the private GEO
submission can be accessed using a reviewer token provided through the journal.

Raw sequencing data (FASTQ) from the in-house cohorts are not publicly available because
the original informed consent did not authorize public release of individual-level
sequencing data. Raw in-house sequencing data are not deposited in GEO and are not
planned for public release.

## Public validation cohort

Raw sequencing data for the public validation cohort (PVC) are publicly available from
NCBI BioProject **PRJNA588242**.

## Code and derived outputs

Analysis code, configuration files, reproducibility documentation, and supporting
processed data that are permitted for public sharing are provided in this GitHub
repository.

## What is where

| Material | Location | Status |
|----------|----------|--------|
| Analyses code and workflow | this repository (`code/`, `docs/`) | Public (this repository) |
| Fixed parameters, locked panel, locked thresholds | this repository (`config/`) | Public (this repository) |
| De-identified per-library annotation | this repository (`metadata/sample_dictionary.tsv`) | Public (this repository) |
| Lightweight derived results | this repository (`results/`) | Public (this repository) |
| Processed expression matrices and associated sample metadata (IDC / ITC; PVC processed matrices as deposited) | NCBI GEO **GSE348100** | Deposited; currently private; planned public release upon manuscript acceptance |
| Public pediatric validation cohort (raw RNA-seq) | NCBI BioProject **PRJNA588242** | Publicly available |
| In-house raw sequencing (IDC / ITC FASTQ) | — | Restricted; not publicly available |

## What is deliberately not in this repository

- FASTQ / BAM / CRAM / SAM and any raw or intermediate alignment files.
- In-house individual-level sequencing data. The original informed consent did not
  authorize public release of individual-level sequencing data.
- Full-scale expression matrices that are already deposited in GEO.
- Clinical source records (participant names, medical-record numbers, contact details,
  addresses or other direct identifiers). Only de-identified, published annotation is
  provided.

## Notes

- `metadata/sample_dictionary.tsv` uses the public sample identifiers
  (`IDC_001…IDC_109`, `ITC_001…ITC_035`, and the public PVC sample accessions) and the
  same group labels as the de-identified cohort sample metadata. No internal identifiers are used.
- Sex and age are included because they are required by the age- and sex-adjusted
  sensitivity analyses and are already reported in the de-identified cohort sample metadata.
- No GEO reviewer token is stored in this repository; tokens for peer review are provided
  through the journal editorial system.
