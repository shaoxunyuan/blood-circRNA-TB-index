# Data availability

## In-house cohorts

Processed expression matrices and associated sample metadata have been deposited in NCBI
GEO under accession **GSE348100** and will be released upon publication.

Raw sequencing data from the in-house cohorts are not publicly available because the
original informed consent did not authorize public release of individual-level sequencing
data.

## Public validation cohort

Raw sequencing data for the public validation cohort are publicly available from NCBI
BioProject **PRJNA588242**.

## Code and derived outputs

Analysis code, configuration files, reproducibility documentation, and selected
lightweight derived outputs are provided in this GitHub repository.

## What is where

| Material | Location |
|----------|----------|
| Analyses code and workflow | this repository (`code/`, `docs/`) |
| Fixed parameters, locked panel, locked thresholds | this repository (`config/`) |
| De-identified per-library annotation | this repository (`metadata/sample_dictionary.tsv`) |
| Lightweight derived results | this repository (`results/`) |
| Processed circRNA-derived expression matrices (IDC / ITC / PVC) | NCBI GEO **GSE348100** |
| Supplementary Tables S1–S8 (sample annotation, expression matrices, differential expression, TB-index scores, diagnostics, benchmarks, BSJ) | manuscript submission package |
| Public pediatric validation cohort (raw RNA-seq) | NCBI BioProject **PRJNA588242** |

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
  same group labels as Supplementary Table S1. No internal identifiers are used.
- Sex and age are included because they are required by the age- and sex-adjusted
  sensitivity analyses and are already reported in Supplementary Table S1.
