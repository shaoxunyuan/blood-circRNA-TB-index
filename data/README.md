# data/ — analysis inputs (placeholders)

No expression matrix is redistributed in this repository. This directory defines
**where** the deposited inputs belong so that the scripts in `code/` can be re-run.

## Data sources

| Data | Source | Status |
|------|--------|--------|
| In-house raw sequencing (IDC, ITC) | — | **Not publicly available.** The original informed consent did not authorize public release of individual-level sequencing data. Raw FASTQ are not deposited in GEO and are not planned for public release. |
| Processed expression matrices and associated sample metadata (IDC, ITC, PVC) | NCBI GEO **GSE348100** | Deposited; currently private; planned public release upon manuscript acceptance (peer-review access via GEO reviewer token provided through the journal) |
| Public pediatric validation cohort (raw sequencing) | NCBI BioProject **PRJNA588242** | Publicly available |

> Processed expression matrices and associated sample metadata from IDC and ITC have been
> deposited in NCBI GEO under accession GSE348100. The GEO record is currently private and
> is planned to be made publicly available upon manuscript acceptance.

## Layout

```
data/
├── metadata/    per-cohort sample annotation used as the analysis input
├── qn/          quantile-normalised expression matrices
│   └── raw/     gene-level TPM matrices before normalisation
├── bsj/         isoform-level circular quant subsets for the three panel genes
└── benchmark/   circRNA and linear-RNA QN matrices for the Sweeney3 / RISK6 benchmark
```

Each subdirectory has its own `README.md` stating:

1. which files are required;
2. where to obtain them;
3. the corresponding GEO accession / supporting table / public repository;
4. the main column names;
5. which scripts read the file.

## Contract for the scripts

Every script resolves this directory through `code/_paths.R` and stops with an
explicit message when a file is missing:

```
Required input file not found.
  File : data/qn/IDC_QN_matrix.tsv
  What : IDC quantile-normalised log2(TPM+1) host-gene matrix
  Where: NCBI GEO GSE348100
```

Bulk data are excluded from version control by `.gitignore` (`raw/`, `fastq/`, `bam/`,
`cram/`), so the placeholders stay in the repository while the matrices do not.
