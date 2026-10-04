#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Purpose:      Build the common circRNA host-gene feature space shared by the
#               three cohorts (IDC, ITC, PVC) and write one common-gene matrix
#               per cohort.
# Inputs:       data/qn/raw/{IDC,ITC,PVC}_genelevel.txt
#                 LONG table, one row per host gene x library:
#                 gene_id, gene_name, gene_biotype, type, BioSample, TPM,
#                 NumReads, n_isoforms. Only type == "circular" rows are used.
#               data/metadata/{IDC,ITC,PVC}_metadata.tsv   (Sample_ID, Group)
# Outputs:      results/common_feature_genes.tsv        (2,754 gene symbols)
#               results/feature_filter_counts.tsv        (per-cohort counts)
#               data/qn/{IDC,ITC,PVC}_common_2754.tsv    (genes x samples, TPM)
# Key params:   Detectability filter: within a cohort a host gene is retained
#               when at least three libraries of any single clinical group have
#               TPM > 0. The common feature space is the intersection of the
#               three retained sets. Every library of a cohort must appear in
#               its pre-filter table; libraries present in the table but absent
#               from the final metadata are dropped.
# Dependencies: R >= 4.4.3; data.table
# ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(data.table))
source("code/_paths.R")

COHORTS <- c("IDC", "ITC", "PVC")
EXPECTED_LIBRARIES <- c(IDC = 109L, ITC = 35L, PVC = 81L)
EXPECTED_RETAINED  <- c(IDC = 3535L, ITC = 3058L, PVC = 4457L)
EXPECTED_GENES <- 2754L
MIN_DETECT_SAMPLES <- 3L

## ---- 1. Cohort metadata ----------------------------------------------------
read_meta <- function(coh) {
  f <- file.path(DATA_META, paste0(coh, "_metadata.tsv"))
  paths_require(f, paste0(coh, " per-library metadata (Sample_ID, Group)"),
                "Supplementary Table S1; see data/metadata/README.md")
  m <- fread(f)
  if (anyDuplicated(m$Sample_ID)) stop("Duplicate Sample_ID in ", coh, " metadata", call. = FALSE)
  m
}

## ---- 2. Long pre-filter table ---------------------------------------------
read_long <- function(coh) {
  f <- file.path(DATA_QN_RAW, paste0(coh, "_genelevel.txt"))
  paths_require(f, paste0(coh, " pre-filter gene-level host-gene table (long format)"),
                "NCBI GEO GSE348100; see data/qn/raw/README.md")
  d <- fread(f)
  need <- c("gene_name", "type", "BioSample", "TPM")
  miss <- setdiff(need, names(d))
  if (length(miss)) {
    stop("Missing column(s) in ", basename(f), ": ", paste(miss, collapse = ", "), call. = FALSE)
  }
  d[type == "circular" & !is.na(gene_name) & gene_name != ""]
}

## ---- 3. Per-cohort detectability filter -----------------------------------
retained <- list(); counts <- list()
for (coh in COHORTS) {
  meta <- read_meta(coh)
  samples <- as.character(meta$Sample_ID)
  if (length(samples) != EXPECTED_LIBRARIES[[coh]]) {
    stop(sprintf("%s: %d libraries in the metadata but %d are expected",
                 coh, length(samples), EXPECTED_LIBRARIES[[coh]]), call. = FALSE)
  }

  d <- read_long(coh)
  absent <- setdiff(samples, unique(d$BioSample))
  if (length(absent)) {
    stop(sprintf("%s: %d metadata sample(s) are absent from the pre-filter table (first: %s). The pre-filter table must contain all libraries listed in the final cohort metadata; extra libraries are allowed and are dropped explicitly; see data/qn/raw/README.md.",
                 coh, length(absent), absent[1]), call. = FALSE)
  }
  extra <- setdiff(unique(d$BioSample), samples)
  if (length(extra)) {
    message(sprintf("  %s: %d library(ies) in the pre-filter table are not part of the final cohort and were dropped",
                    coh, length(extra)))
  }

  d <- d[BioSample %in% samples]
  grp <- setNames(as.character(meta$Group), samples)
  d[, Group := grp[BioSample]]

  det <- d[TPM > 0, .(n_detected = uniqueN(BioSample)), by = .(gene_name, Group)]
  keep <- det[, .(max_group_detected = max(n_detected)), by = gene_name][
    max_group_detected >= MIN_DETECT_SAMPLES, gene_name]

  retained[[coh]] <- keep
  counts[[coh]] <- data.table(Cohort = coh, Libraries = length(samples),
                              Genes_in_prefilter = uniqueN(d$gene_name),
                              Genes_retained = length(keep))
  message(sprintf("  %s: %d libraries, %d of %d genes retained (expected %d)",
                  coh, length(samples), length(keep), uniqueN(d$gene_name),
                  EXPECTED_RETAINED[[coh]]))
}

common <- Reduce(intersect, retained)
counts_tab <- rbind(rbindlist(counts),
                    data.table(Cohort = "intersection", Libraries = NA_integer_,
                               Genes_in_prefilter = NA_integer_,
                               Genes_retained = length(common)), fill = TRUE)
fwrite(counts_tab, file.path(RESULTS, "feature_filter_counts.tsv"), sep = "\t")

## ---- 4. Common-gene matrices (genes x samples, TPM) ------------------------
for (coh in COHORTS) {
  meta <- read_meta(coh); samples <- as.character(meta$Sample_ID)
  d <- read_long(coh)
  d <- d[BioSample %in% samples & gene_name %in% common]
  M <- dcast(d, gene_name ~ BioSample, value.var = "TPM", fun.aggregate = sum, fill = 0)
  M <- M[match(common, M$gene_name)]
  setnames(M, "gene_name", "Gene")
  M <- M[, c("Gene", samples), with = FALSE]
  stopifnot(nrow(M) == EXPECTED_GENES, ncol(M) - 1L == EXPECTED_LIBRARIES[[coh]])
  fwrite(M, file.path(DATA_QN, sprintf("%s_common_%d.tsv", coh, EXPECTED_GENES)), sep = "\t")
  message(sprintf("  wrote data/qn/%s_common_%d.tsv (%d genes x %d samples)",
                  coh, EXPECTED_GENES, nrow(M), ncol(M) - 1L))
}

fwrite(data.table(Gene = common), file.path(RESULTS, "common_feature_genes.tsv"), sep = "\t")

## ---- 5. QA -----------------------------------------------------------------
for (coh in COHORTS) {
  n <- counts_tab[Cohort == coh, Genes_retained]
  if (n != EXPECTED_RETAINED[[coh]]) {
    stop(sprintf("Feature-filter mismatch for %s: %d genes retained, %d expected",
                 coh, n, EXPECTED_RETAINED[[coh]]), call. = FALSE)
  }
}
if (length(common) != EXPECTED_GENES) {
  stop(sprintf("Common feature space: %d genes, %d expected",
               length(common), EXPECTED_GENES), call. = FALSE)
}
message("QA: retained ", paste(sprintf("%s=%d", COHORTS, EXPECTED_RETAINED), collapse = " "),
        "; intersection = ", length(common))
