#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Purpose:      Perform cross-cohort BSJ (back-splice junction) overlap and
#               contribution analysis for the locked TB-index genes SOD2,
#               GOLPH3 and PCNX2, in IDC, ITC and PVC.
# Inputs:       data/bsj/{IDC,ITC,PVC}_GOLPH3_SOD2_PCNX2_circular.tsv
#                 Long table, one row per circular isoform x library, for the
#                 three panel genes. Required columns:
#                 gene_name, BSJ_ID, BioSample, TPM.
#                 BSJ_ID identifies the junction as
#                 chromosome|donor|acceptor|strand (hg38), so a BSJ is matched
#                 across cohorts by genomic coordinates and never by row index.
#                 The BSJ input tables are restricted internally to the final
#                 cohort sample identifiers: libraries present in a table but
#                 not listed in the cohort metadata are reported and dropped.
#               data/metadata/{IDC,ITC,PVC}_metadata.tsv
#                 Final per-library annotation; only Sample_ID is used, to
#                 define final cohort membership (IDC 109 / ITC 35 / PVC 81
#                 libraries). Membership is never inferred from time point,
#                 baseline / follow-up status or treatment status.
# Outputs:      results/BSJ_overlap_summary.tsv      per-gene cohort overlap counts
#               results/BSJ_contribution_detail.tsv  per gene x cohort x BSJ detail
#               results/BSJ_summary.tsv              per gene x cohort summary
# Definitions:
#               all-three-cohort BSJ  detected in IDC, ITC and PVC
#               shared BSJ            detected in at least two cohorts
#               cohort-specific BSJ   detected in exactly one cohort
#               largest shared BSJ    within a gene and cohort, the shared
#                                     junction with the highest cohort-pooled TPM
#               cohort-pooled TPM     sum of TPM over the libraries of that cohort
#               contribution_pct      pooled_TPM_BSJ / total_gene_circular_TPM * 100
# Dependencies: R >= 4.4.3; data.table
# ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(data.table))
source("code/_paths.R")

GENES <- c("SOD2", "GOLPH3", "PCNX2")
COHORTS <- c("IDC", "ITC", "PVC")
N_LIB <- c(IDC = 109L, ITC = 35L, PVC = 81L)   # locked final cohort sizes

## ---- 1. Final cohort membership --------------------------------------------
read_final_ids <- function(coh) {
  f <- file.path(DATA_META, paste0(coh, "_metadata.tsv"))
  paths_require(f, paste0(coh, " final per-library metadata (Sample_ID)"),
                "Supplementary Table S1; see data/metadata/README.md")
  m <- fread(f)
  if (!"Sample_ID" %in% names(m)) {
    stop("Missing column 'Sample_ID' in ", basename(f), call. = FALSE)
  }
  ids <- unique(as.character(m$Sample_ID))
  if (anyDuplicated(ids)) stop("Duplicate Sample_ID in ", basename(f), call. = FALSE)
  if (length(ids) != N_LIB[[coh]]) {
    stop(sprintf("%s metadata lists %d libraries but the final cohort has %d; see data/metadata/README.md",
                 coh, length(ids), N_LIB[[coh]]), call. = FALSE)
  }
  ids
}

## ---- 2. Per-cohort BSJ tables, restricted to the final cohort --------------
read_bsj <- function(coh, final_ids) {
  f <- file.path(DATA_BSJ, sprintf("%s_GOLPH3_SOD2_PCNX2_circular.tsv", coh))
  paths_require(f, paste0(coh, " circular isoform table for the three panel genes"),
                "NCBI GEO GSE348100; see data/bsj/README.md")
  d <- fread(f, sep = "\t")
  need <- c("gene_name", "BSJ_ID", "BioSample", "TPM")
  miss <- setdiff(need, names(d))
  if (length(miss)) {
    stop("Missing column(s) in ", basename(f), ": ", paste(miss, collapse = ", "),
         ". Required: ", paste(need, collapse = ", "), call. = FALSE)
  }
  d <- d[gene_name %in% GENES & !is.na(BSJ_ID) & BSJ_ID != ""]
  d[, BioSample := as.character(BioSample)]

  ## Final cohort membership is enforced here from Sample_ID; the input table is
  ## never trusted to be pre-restricted.
  extra_ids   <- setdiff(unique(d$BioSample), final_ids)
  missing_ids <- setdiff(final_ids, unique(d$BioSample))
  if (length(missing_ids) == length(final_ids)) {
    stop(sprintf("%s: none of the %d final cohort Sample_IDs occur in %s. The metadata and the BSJ table appear to use different sample identifiers; see data/metadata/README.md.",
                 coh, length(final_ids), basename(f)), call. = FALSE)
  }
  if (length(extra_ids)) {
    message(sprintf("  %s: dropping %d BSJ-input library(ies) not present in the final cohort metadata (%s%s)",
                    coh, length(extra_ids), paste(head(extra_ids, 5), collapse = ", "),
                    if (length(extra_ids) > 5) ", ..." else ""))
  }
  if (length(missing_ids)) {
    message(sprintf("  %s: %d final-cohort library(ies) have no quantified circular isoform row for the panel genes and contribute zero (%s%s)",
                    coh, length(missing_ids), paste(head(missing_ids, 5), collapse = ", "),
                    if (length(missing_ids) > 5) ", ..." else ""))
  }
  d <- d[BioSample %in% final_ids]
  if (!nrow(d)) {
    stop(sprintf("%s: no BSJ rows remain after restricting to the final cohort", coh), call. = FALSE)
  }
  if (uniqueN(d$BioSample) != length(final_ids) - length(missing_ids)) {
    stop(sprintf("%s: restriction to the final cohort did not converge", coh), call. = FALSE)
  }
  message(sprintf("  %s: %d of %d final-cohort libraries contribute quantified panel-gene junction rows",
                  coh, uniqueN(d$BioSample), length(final_ids)))
  d[, .(gene_name, BSJ_ID, BioSample, TPM = as.numeric(TPM))]
}

final_ids <- lapply(COHORTS, read_final_ids); names(final_ids) <- COHORTS
raw <- Map(read_bsj, COHORTS, final_ids); names(raw) <- COHORTS

## Cohort-pooled TPM per gene x BSJ (one row per detected junction).
pooled <- lapply(raw, function(d) d[, .(pooled_TPM = sum(TPM)), by = .(gene_name, BSJ_ID)])

## ---- 2. Detection sets and overlap categories ------------------------------
detected <- lapply(pooled, function(d) d[pooled_TPM > 0])
sets <- lapply(detected, function(d) split(d$BSJ_ID, d$gene_name))

overlap_rows <- list(); detail_rows <- list(); summary_rows <- list()
for (g in GENES) {
  sets_g <- lapply(sets, function(s) if (is.null(s[[g]])) character(0) else s[[g]])
  in_all3 <- Reduce(intersect, sets_g)
  ## shared = present in at least two cohorts
  all_ids <- unique(unlist(sets_g))
  n_cohorts <- vapply(all_ids, function(id) sum(vapply(sets_g, function(s) id %in% s, logical(1))), integer(1))
  shared <- all_ids[n_cohorts >= 2]
  specific <- all_ids[n_cohorts == 1]

  overlap_rows[[g]] <- data.table(
    Gene = g,
    IDC_only = sum(n_cohorts == 1 & vapply(all_ids, function(id) id %in% sets_g$IDC, logical(1))),
    ITC_only = sum(n_cohorts == 1 & vapply(all_ids, function(id) id %in% sets_g$ITC, logical(1))),
    PVC_only = sum(n_cohorts == 1 & vapply(all_ids, function(id) id %in% sets_g$PVC, logical(1))),
    IDC_ITC_only = sum(n_cohorts == 2 & vapply(all_ids, function(id) id %in% sets_g$IDC & id %in% sets_g$ITC, logical(1))),
    IDC_PVC_only = sum(n_cohorts == 2 & vapply(all_ids, function(id) id %in% sets_g$IDC & id %in% sets_g$PVC, logical(1))),
    ITC_PVC_only = sum(n_cohorts == 2 & vapply(all_ids, function(id) id %in% sets_g$ITC & id %in% sets_g$PVC, logical(1))),
    All_three = length(in_all3))

  for (coh in COHORTS) {
    d <- detected[[coh]][gene_name == g]
    total <- sum(d$pooled_TPM)
    if (total <= 0) stop("No circular TPM for ", g, " in ", coh, call. = FALSE)

    d[, category := fifelse(BSJ_ID %in% specific, "Cohort-specific BSJs",
                     fifelse(BSJ_ID %in% shared,   "Shared BSJ", "Unclassified"))]
    sh <- d[category == "Shared BSJ"]
    largest_id <- if (nrow(sh)) sh[which.max(pooled_TPM), BSJ_ID] else NA_character_
    d[, category := fifelse(BSJ_ID %in% largest_id, "Largest shared BSJ",
                     fifelse(BSJ_ID %in% shared, "Other shared BSJs", category))]
    d[, `:=`(Cohort = coh, Gene = g, Contribution_pct = 100 * pooled_TPM / total,
             Is_largest_shared = BSJ_ID %in% largest_id)]
    detail_rows[[paste(g, coh)]] <- d[, .(Gene, Cohort, BSJ_ID, Detection_category = category,
                                          Is_largest_shared, Cohort_pooled_TPM = pooled_TPM,
                                          Contribution_pct)]

    pct <- function(cat) sum(d[category == cat, pooled_TPM]) / total * 100
    summary_rows[[paste(g, coh)]] <- data.table(
      Gene = g, Cohort = coh,
      All_three_cohort_BSJ_n = length(in_all3),
      Shared_BSJ_n = nrow(sh),
      Cohort_specific_BSJ_n = nrow(d[category == "Cohort-specific BSJs"]),
      Largest_shared_BSJ = largest_id,
      Largest_shared_contribution_pct = pct("Largest shared BSJ"),
      Other_shared_contribution_pct = pct("Other shared BSJs"),
      Cohort_specific_contribution_pct = pct("Cohort-specific BSJs"),
      Shared_total_contribution_pct = pct("Largest shared BSJ") + pct("Other shared BSJs"))
  }
}

overlap_tab <- rbindlist(overlap_rows)
detail_tab  <- rbindlist(detail_rows)
summary_tab <- rbindlist(summary_rows)

fwrite(overlap_tab, file.path(RESULTS, "BSJ_overlap_summary.tsv"), sep = "\t")
fwrite(detail_tab,  file.path(RESULTS, "BSJ_contribution_detail.tsv"), sep = "\t")
fwrite(summary_tab, file.path(RESULTS, "BSJ_summary.tsv"), sep = "\t")

## ---- 3. QA -----------------------------------------------------------------
EXPECTED_ALL3 <- c(SOD2 = 1L, GOLPH3 = 1L, PCNX2 = 12L)
EXPECTED_LARGEST <- data.table(
  Gene = rep(GENES, each = 3), Cohort = rep(COHORTS, times = 3),
  pct = c(95.5, 89.7, 73.5,    # SOD2
          99.1, 99.7, 100.0,   # GOLPH3
          26.4, 21.8, 24.9))   # PCNX2
EXPECTED_SHARED_TOTAL <- data.table(
  Gene = c("SOD2","SOD2","SOD2"), Cohort = COHORTS, pct = c(95.5, 100.0, 77.0))

fail <- character(0)
for (g in GENES) {
  got <- overlap_tab[Gene == g, All_three]
  if (got != EXPECTED_ALL3[[g]]) fail <- c(fail, sprintf("%s all-three-cohort BSJs: %d vs %d", g, got, EXPECTED_ALL3[[g]]))
}
for (i in seq_len(nrow(EXPECTED_LARGEST))) {
  e <- EXPECTED_LARGEST[i]
  got <- summary_tab[Gene == e$Gene & Cohort == e$Cohort, Largest_shared_contribution_pct]
  if (abs(got - e$pct) > 0.05) fail <- c(fail, sprintf("%s %s largest shared: %.1f vs %.1f", e$Gene, e$Cohort, got, e$pct))
}
for (i in seq_len(nrow(EXPECTED_SHARED_TOTAL))) {
  e <- EXPECTED_SHARED_TOTAL[i]
  got <- summary_tab[Gene == e$Gene & Cohort == e$Cohort, Shared_total_contribution_pct]
  if (abs(got - e$pct) > 0.05) fail <- c(fail, sprintf("%s %s shared total: %.1f vs %.1f", e$Gene, e$Cohort, got, e$pct))
}
## Each gene x cohort decomposition must sum to 100%.
summary_tab[, total_pct := Largest_shared_contribution_pct +
                         Other_shared_contribution_pct +
                         Cohort_specific_contribution_pct]
bad <- summary_tab[abs(total_pct - 100) > 1e-6]
if (nrow(bad)) fail <- c(fail, sprintf("contribution does not sum to 100%% for %s",
                                       paste(bad$Gene, bad$Cohort, sep = "/", collapse = ", ")))
if (nrow(summary_tab) != 9L) fail <- c(fail, "expected 9 gene x cohort rows")
if (length(fail)) stop("BSJ QA failed:\n  ", paste(fail, collapse = "\n  "), call. = FALSE)

message("QA: all-three-cohort BSJs ", paste(sprintf("%s=%d", GENES, EXPECTED_ALL3), collapse = " "),
        "; largest-shared contributions and shared totals match Figure 4B")
print(summary_tab[, .(Gene, Cohort, All_three_cohort_BSJ_n, Shared_BSJ_n,
                      Largest_shared_contribution_pct, Other_shared_contribution_pct,
                      Cohort_specific_contribution_pct, Shared_total_contribution_pct)])
