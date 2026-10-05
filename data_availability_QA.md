# Data availability QA

Repository: `blood-circRNA-TB-index` (working copy: `github_upload_package/`)
QA date: 2026-10-05

Scope: this round changed data-availability wording only. No analysis code, threshold,
AUROC, bootstrap, figure/table result, metadata count, GEO accession or BioProject
accession was modified.

---

## GSE348100 — in-house processed data

| Item | Value |
|------|-------|
| Accession | **GSE348100** (correct; no typo, no `TBD`) |
| Deposited | **YES** |
| Public now | **NO** (submission complete, currently private) |
| Release upon publication | **YES** |
| Content stated | processed expression matrices **and associated sample metadata** |

## In-house raw sequencing (IDC / ITC)

| Item | Value |
|------|-------|
| Publicly available | **NO** |
| Reason documented | **YES** — original informed consent did not authorize public release of individual-level sequencing data |
| Future-upload promise | **NONE** (no "will be deposited upon acceptance") |
| Data-access request route | **not stated** — the manuscript documents no such mechanism, so none was invented |

## Public validation cohort

| Item | Value |
|------|-------|
| Source | **NCBI BioProject PRJNA588242** (a BioProject, not a GEO accession) |
| Wording | "Raw sequencing data for the public validation cohort are publicly available from NCBI BioProject PRJNA588242." |

## GitHub repository

| Item | Value |
|------|-------|
| Role stated | analysis code, configuration files, reproducibility documentation, selected lightweight derived results |
| Implies raw-data hosting | **NO** |

---

## Files checked

- `README.md`
- `docs/data_availability.md`
- `docs/reproducibility.md`
- `docs/analysis_workflow.md`
- `data/README.md`
- `data/qn/README.md`
- `data/qn/raw/README.md`
- `data/benchmark/README.md`
- `data/metadata/README.md`
- `data/bsj/README.md`
- `CITATION.cff`
- `code/README.md`

## Canonical wording present

| Phrase | Occurrences |
|--------|-------------|
| `have been deposited in NCBI GEO` | 2 |
| `will be released upon publication` | 6 |
| `Raw sequencing data from the in-house cohorts are not publicly available` | 3 |
| `original informed consent` | 6 |
| `publicly available from NCBI BioProject PRJNA588242` | 0 |

## Contradictory wording

| Pattern | Occurrences (outside this QA file and the audit file) |
|---------|------------------------------------------------------|
| `will be deposited` | 0 |
| `upon acceptance` | 0 |
| `raw sequencing will` | 0 |
| `raw data will` | 0 |
| `accession TBD` | 0 |
| `GEO TBD` | 0 |
| `BioProject TBD` | 0 |
| `are publicly available in GEO` | 0 |
| `GEO package upon acceptance` | 0 |
| `not publicly released` | 0 |

All contradictory patterns: **0**.

---

## Final status: **PASS**

Pass conditions:

1. GSE348100 is no longer described as "will be deposited" — YES
2. GSE348100 is not described as currently public/available — YES
3. "have been deposited" is stated — YES
4. "will be released upon publication" is stated — YES
5. In-house raw sequencing promises no future upload — YES
6. PRJNA588242 is correct and used as a BioProject — YES
7. The GitHub repository is positioned as a code/reproducibility repository — YES
8. No accession TBD — YES
9. No mutually contradictory data-availability statements — YES

Note: `data_availability_audit_before.txt` quotes the pre-revision wording on purpose and
is therefore kept outside this repository (see the project `操作日志.md`).
