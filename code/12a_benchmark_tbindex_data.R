#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Purpose:      Assemble per-sample TB-index and Sweeney3 / RISK6 scores.
# Inputs:       circRNA TB-index per sample; linear-RNA QN matrices; sample metadata.
# Outputs:      Benchmark workbook sheets underlying results/DeLong_comparisons.tsv
#               (the literature-benchmark results).
# Key params:   Locked panel SOD2 (+1) / GOLPH3 (-1) / PCNX2 (-1); published Sweeney3;
#               RISK6 = GM(GBP2,FCGR1B,SERPING1) - GM(TUBGCP6,TRMT2A,SDR39U1)
#               (Penn-Nicholson RNA-seq/microarray geometric-mean formulation; column RISK6).
# Note:         QN log2(TPM+1) with GM on linear (TPM+1) scale. Not original qRT-PCR RISK6.
# Dependencies: R >= 4.4.3; data.table; openxlsx.
# ---------------------------------------------------------------------------
# =============================================================================
# build_benchmark_tbindex_data.R — 整合 circRNA TB-index + Sweeney3/RISK6 到 benchmark workbook
#
# 数据来源:
#   circRNA TB-index: data/benchmark/{IDC,ITC,PVC}_TBindex_per_sample.csv
#   Sweeney3/RISK6:   data/benchmark/linear_qn_*.rds (linear QN)
#   Sample metadata:  data/metadata/{IDC,ITC,PVC}_metadata.tsv  (the de-identified cohort sample metadata)
#
# 输出:
#   在 Figure_benchmark_TBindex.xlsx 写入 "per_sample_scores"
#   (Sample_ID, Cohort, Group, circRNA_TB_index,
#    GBP5, DUSP3, KLF2, Sweeney3, RISK6, TimePoint)
#   GBP5/DUSP3/KLF2 = Sweeney3 三基因的 linearRNA QN 表达;
#   Sweeney3 = (GBP5 + DUSP3)/2 - KLF2
#

suppressPackageStartupMessages({
  library(data.table); library(openxlsx)
})

## Package-relative paths (do not hard-code absolute workstation paths)
.args_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.code_dir <- if (length(.args_file)) dirname(normalizePath(sub("^--file=", "", .args_file[1]))) else getwd()
source(file.path(.code_dir, "_paths.R"))
## Per-cohort metadata TSVs are the single metadata input used across the
## repository; no workbook form is maintained.
META_FILES <- file.path(DATA_META, sprintf("%s_metadata.tsv", c("IDC", "ITC", "PVC")))
S8   <- file.path(DATA_BENCH, "benchmark_tbindex.xlsx")
# shipped per-sample TB-index and linear QN
V4_IDC <- file.path(DATA_BENCH, "IDC_TBindex_per_sample.csv")
V4_ITC <- file.path(DATA_BENCH, "ITC_TBindex_per_sample.csv")
V4_PVC <- file.path(DATA_BENCH, "PVC_TBindex_per_sample.csv")
LIN_IDC <- file.path(DATA_BENCH, "linear_qn_IDC.rds")
LIN_ITC <- file.path(DATA_BENCH, "linear_qn_ITC.rds")
LIN_PVC <- file.path(DATA_BENCH, "linear_qn_PVC.rds")


cat("================================================================\n")
cat(" build_benchmark_tbindex_data.R\n")
cat(" 整合 TBindex + Sweeney3/RISK6 → benchmark per_sample_scores\n")
cat("================================================================\n")

# ---- 1. circRNA TB-index ----
cat("\n[1/4] 读取 circRNA TB-index ...\n")
load_tb <- function(path) {
  d <- fread(path)
  bs <- intersect(c("BioSample", "Run"), names(d))[1]
  setnames(d, bs, "BioSample")
  d[, .(BioSample, TBindex)]
}
idc_tb <- load_tb(V4_IDC)
itc_tb <- load_tb(V4_ITC)
pvc_tb <- load_tb(V4_PVC)
idc_tb[, Cohort := "IDC"]; itc_tb[, Cohort := "ITC"]; pvc_tb[, Cohort := "PVC"]
# PVC: the LTBI group is made of incident converters (post-conversion samples);
# no baseline LTBI library exists, so the LTBI TB-index values are used as supplied.
circ_tb <- rbindlist(list(
  idc_tb[, .(TBindex = mean(TBindex)), by = .(BioSample, Cohort)],
  itc_tb[, .(TBindex = mean(TBindex)), by = .(BioSample, Cohort)],
  pvc_tb[, .(TBindex = mean(TBindex)), by = .(BioSample, Cohort)]
), fill = TRUE)
cat(sprintf("  TBindex: %d 样本\n", nrow(circ_tb)))

# ---- 2. Sweeney3 / RISK6 (linear QN; Penn-Nicholson RNA-seq/microarray GM formula) ----
# RISK6 = GM(GBP2,FCGR1B,SERPING1) - GM(TUBGCP6,TRMT2A,SDR39U1);
# GM on linear (TPM+1) scale from QN log2(TPM+1) = 2^mean(log2). Column name remains RISK6.

cat("\n[2/4] 计算 Sweeney3 / RISK6 ...\n")
load_qn <- function(path) { x <- readRDS(path); rownames(x$logTPM_QN) <- x$genes; x$logTPM_QN }
gm3_log2 <- function(a, b, c) 2^((as.numeric(a) + as.numeric(b) + as.numeric(c)) / 3)
score_from_mat <- function(mat, samples) {
  m <- match(sub("^LTBI_", "LTB_", samples), colnames(mat))
  ok <- !is.na(m)
  data.table(BioSample = samples[ok],
    # Sweeney3 三基因的 linearRNA QN 表达 (log2(TPM+1))
    GBP5  = as.numeric(mat["GBP5",  m[ok]]),
    DUSP3 = as.numeric(mat["DUSP3", m[ok]]),
    KLF2  = as.numeric(mat["KLF2",  m[ok]]),
    # Sweeney3 = (GBP5 + DUSP3)/2 - KLF2
    Sweeney3 = as.numeric((mat["GBP5", m[ok]] + mat["DUSP3", m[ok]])/2 - mat["KLF2", m[ok]]),
    RISK6 = as.numeric(
      gm3_log2(mat["GBP2", m[ok]], mat["FCGR1B", m[ok]], mat["SERPING1", m[ok]]) -
      gm3_log2(mat["TUBGCP6", m[ok]], mat["TRMT2A", m[ok]], mat["SDR39U1", m[ok]])))
}
qn_list <- list(
  IDC = load_qn(LIN_IDC),
  ITC = load_qn(LIN_ITC),
  PVC = load_qn(LIN_PVC))
sc_all <- rbindlist(lapply(c("IDC", "ITC", "PVC"), function(c) {
  bios <- unique(circ_tb[Cohort == c, BioSample])
  qn <- qn_list[[c]]
  if (mean(bios %in% colnames(qn)) > 0.5) {
    s <- score_from_mat(qn, bios); s[, Cohort := c]; return(s)
  }
  NULL
}), fill = TRUE)
cat(sprintf("  Scores: %d 样本\n", nrow(sc_all)))

# ---- 3. Per-cohort metadata (TSV) ----
cat("\n[3/4] 读取各队列 metadata (TSV) ...\n")
clin <- rbindlist(lapply(seq_along(c("IDC", "ITC", "PVC")), function(i) {
  coh <- c("IDC", "ITC", "PVC")[i]
  f <- META_FILES[i]
  paths_require(f, paste0(coh, " per-library metadata"),
                "the de-identified cohort sample metadata; see data/metadata/README.md")
  d <- fread(f)
  d <- d[!is.na(Group) & Group != ""]
  data.table(BioSample = as.character(d$Sample_ID), Cohort = coh,
             Group = as.character(d$Group),
             TimePoint = if ("TimePoint" %in% names(d)) as.character(d$TimePoint) else NA_character_)
}))
cat(sprintf("  Clinical: %d 样本\n", nrow(clin)))

# ---- 4. 合并 → per_sample_scores ----
cat("\n[4/4] 合并写入 benchmark workbook ...\n")
ps <- merge(circ_tb, sc_all, by = c("BioSample", "Cohort"), all.x = TRUE)
ps <- merge(ps, clin, by = c("BioSample", "Cohort"), all.x = TRUE)
ps <- ps[!is.na(Group) & !is.na(Sweeney3) & !is.na(RISK6)]
setorder(ps, Cohort, Group, BioSample)
setnames(ps, c("TBindex", "Sweeney3", "RISK6"), c("circRNA_TB_index", "Sweeney3", "RISK6"))
cat(sprintf("  per_sample_scores: %d 样本\n", nrow(ps)))

wb <- loadWorkbook(S8)
if ("per_sample_scores" %in% names(wb)) removeWorksheet(wb, "per_sample_scores")
addWorksheet(wb, "per_sample_scores")
writeData(wb, "per_sample_scores", ps)
setColWidths(wb, "per_sample_scores", cols = 1:ncol(ps), widths = "auto")
writeComment(wb, "per_sample_scores", col = 1, row = 1,
  comment = createComment("Input data for the benchmark table (the literature-benchmark results). Sources: circRNA TB-index and linear-RNA QN scores."))
# 放 README 后
nms <- names(wb); rd <- which(nms == "README")
ord <- c(rd, which(nms == "per_sample_scores"), setdiff(seq_along(nms), c(rd, which(nms == "per_sample_scores"))))
worksheetOrder(wb) <- ord

saveWorkbook(wb, S8, overwrite = TRUE)
cat(sprintf("\nDone → %s (%d sheets)\n", S8, length(names(wb))))
cat("  per_sample_scores 已写入\n")
cat("================================================================\n")
