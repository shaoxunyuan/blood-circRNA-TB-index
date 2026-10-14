#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Purpose:      Build the 3 circRNA + 9 linear-RNA expression matrix for the benchmark.
# Inputs:       circRNA QN matrix; linear-RNA QN matrix; benchmark sample set.
# Outputs:      Expression matrix sheet used by the benchmark figures/tables.
# Key params:   circRNA genes GOLPH3 / SOD2 / PCNX2; Sweeney3 genes GBP5, DUSP3, KLF2;
#               RISK6 genes GBP2, SDR39U1, FCGR1B, TUBGCP6, SERPING1, TRMT2A.
# Dependencies: R >= 4.4.3; data.table; openxlsx.
# ---------------------------------------------------------------------------
# =============================================================================
# build_benchmark_expression_matrix.R — 给 benchmark workbook 新增"3 circRNA + 9 linearRNA"表达矩阵
#
# 需求（2026-08-08）：Figure_benchmark_TBindex.xlsx 应包含一个矩阵 =
#   本项目 3 个 circRNA 基因（GOLPH3/SOD2/PCNX2）的表达
#   + Sweeney3（GBP5/DUSP3/KLF2）与 RISK6（GBP2/FCGR1B/SERPING1/TUBGCP6/TRMT2A/SDR39U1）
#   共 9 个 linearRNA 基因的表达。确认现有 xlsx 无此合并矩阵后新增。
#
# 数据来源:
#   circRNA QN: data/benchmark/circ_qn_{IDC,ITC,PVC}.rds
#   linear  QN: data/benchmark/linear_qn_{IDC,ITC,PVC}.rds
#   样本集:    基准分析的共同样本（the literature-benchmark results per-sample score table）
#
# 输出: 新增 sheet "Expression_matrix_3circ_9lin"
#   (Sample_ID, Cohort, Group, TimePoint,
#    circ_GOLPH3, circ_SOD2, circ_PCNX2,
#    lin_GBP5, lin_DUSP3, lin_KLF2, lin_GBP2, lin_FCGR1B, lin_SERPING1,
#    lin_TUBGCP6, lin_TRMT2A, lin_SDR39U1)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
})


## Package-relative paths
.args_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.code_dir <- if (length(.args_file)) dirname(normalizePath(sub("^--file=", "", .args_file[1]))) else getwd()
source(file.path(.code_dir, "_paths.R"))

S8   <- file.path(DATA_BENCH, "benchmark_tbindex.xlsx")
CIRC_QN <- c(
  IDC = file.path(DATA_BENCH, "circ_qn_IDC.rds"),
  ITC = file.path(DATA_BENCH, "circ_qn_ITC.rds"),
  PVC = file.path(DATA_BENCH, "circ_qn_PVC.rds")
)
LIN_QN <- c(
  IDC = file.path(DATA_BENCH, "linear_qn_IDC.rds"),
  ITC = file.path(DATA_BENCH, "linear_qn_ITC.rds"),
  PVC = file.path(DATA_BENCH, "linear_qn_PVC.rds")
)

CIRC3 <- c("GOLPH3", "SOD2", "PCNX2")
LIN9  <- c("GBP5", "DUSP3", "KLF2", "GBP2", "FCGR1B", "SERPING1",
           "TUBGCP6", "TRMT2A", "SDR39U1")
COHORTS <- c("IDC", "ITC", "PVC")
QN_FILE <- c(IDC = "qn_discovery.rds", ITC = "qn_TBnovogene2.rds", PVC = "qn_PRJNA588242.rds")

cat("================================================================\n")
cat(" 给 benchmark workbook 新增 3circ+9lin 表达矩阵\n")
cat("================================================================\n")
stopifnot(file.exists(S8))

# 读取 QN 矩阵: rownames = 基因, colnames = 样本
load_qn <- function(path) { x <- readRDS(path); rownames(x$logTPM_QN) <- x$genes; x$logTPM_QN }

# 样本集与 per_sample_scores 对齐
ps <- as.data.table(read.xlsx(S8, sheet = "per_sample_scores"))
info <- ps[, .(Sample_ID, Cohort, Group, TimePoint)]
cat(sprintf("per_sample_scores 样本: %d\n", nrow(info)))

# 逐队列提取 circ 3 基因 + linear 9 基因表达
rows <- list()
for (c in COHORTS) {
  sid <- ps[Cohort == c, Sample_ID]
  qc <- load_qn(CIRC_QN[[c]])
  ql <- load_qn(LIN_QN[[c]])
  # 样本匹配（含 PVC LTBI_→LTB_ 前缀惯例）
  mc <- match(sid, colnames(qc))
  if (any(is.na(mc))) mc[is.na(mc)] <- match(sub("^LTBI_", "LTB_", sid), colnames(qc))[is.na(mc)]
  ml <- match(sid, colnames(ql))
  if (any(is.na(ml))) ml[is.na(ml)] <- match(sub("^LTBI_", "LTB_", sid), colnames(ql))[is.na(ml)]
  ok <- !is.na(mc) & !is.na(ml)

  dt <- data.table(Sample_ID = sid, Cohort = c)
  for (g in CIRC3) dt[[paste0("circ_", g)]] <- NA_real_
  for (g in LIN9)  dt[[paste0("lin_",  g)]] <- NA_real_
  if (any(ok)) {
    for (g in CIRC3) dt[ok, (paste0("circ_", g)) := as.numeric(qc[g, mc[ok]])]
    for (g in LIN9)  dt[ok, (paste0("lin_",  g)) := as.numeric(ql[g, ml[ok]])]
  }
  cat(sprintf("  %s: %d 样本 (匹配 circ=%d, lin=%d)\n", c, length(sid), sum(!is.na(mc)), sum(!is.na(ml))))
  rows[[c]] <- dt
}
mx <- rbindlist(rows)
mx <- merge(mx, info, by = c("Sample_ID", "Cohort"), all.x = TRUE)
setcolorder(mx, c("Sample_ID", "Cohort", "Group", "TimePoint",
                  paste0("circ_", CIRC3), paste0("lin_", LIN9)))
setorder(mx, Cohort, Group, Sample_ID)
cat(sprintf("合并矩阵: %d 样本 × %d 列（12 个基因）\n", nrow(mx), ncol(mx) - 4))

# 写入 benchmark workbook
wb <- loadWorkbook(S8)
if ("Expression_matrix_3circ_9lin" %in% names(wb)) removeWorksheet(wb, "Expression_matrix_3circ_9lin")
addWorksheet(wb, "Expression_matrix_3circ_9lin")
writeData(wb, "Expression_matrix_3circ_9lin", mx)
setColWidths(wb, "Expression_matrix_3circ_9lin", cols = 1:ncol(mx), widths = "auto")
writeComment(wb, "Expression_matrix_3circ_9lin", col = 1, row = 1,
  comment = createComment(sprintf(
    "Combined expression matrix: 3 circRNA genes (GOLPH3/SOD2/PCNX2) + 9 linearRNA genes (Sweeney3 GBP5/DUSP3/KLF2 + RISK6 GBP2/FCGR1B/SERPING1/TUBGCP6/TRMT2A/SDR39U1). Values = QN-normalized log2(TPM+1); circRNA from data/benchmark/circ_qn_*.rds, linear-RNA from data/benchmark/linear_qn_*.rds. Samples aligned to the benchmark per-sample score table (B=%d).", nrow(mx))))

saveWorkbook(wb, S8, overwrite = TRUE)
cat(sprintf("\nDone → %s (%d sheets)\n", S8, length(names(wb))))
cat("  新增: Expression_matrix_3circ_9lin\n")
cat("================================================================\n")
