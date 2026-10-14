#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Purpose:      Append the Sweeney3 constituent gene expression columns to the
#               per-sample benchmark table.
# Inputs:       Linear-RNA QN expression values for GBP5, DUSP3 and KLF2.
# Outputs:      Extended per-sample benchmark scores sheet.
# Key params:   Sweeney3 = (GBP5 + DUSP3)/2 - KLF2.
# Dependencies: R >= 4.4.3; data.table; openxlsx.
# ---------------------------------------------------------------------------
# =============================================================================
# build_benchmark_add_genes.R — 在 benchmark workbook 的 per_sample_scores 中补充 Sweeney3 三基因
#   (GBP5 / DUSP3 / KLF2) 的 mRNA 表达列, 与已计算的 Sweeney3 score 同表存放
#
# 说明:
#   Sweeney3 score 已由 build_benchmark_tbindex_data.R 计算并存在于 per_sample_scores;
#   本脚本仅把三个基因的表达列补进去, 使该表自包含 (基因表达 → 计算得分)。
#   基因表达来自 linearRNA QN 矩阵 (log2(TPM+1) QN 归一化)。
#
# Input: linear-RNA QN matrices at data/benchmark/linear_qn_{IDC,ITC,PVC}.rds
#   (build_benchmark_tbindex_data.R 中的 results/流程linearRNA_v1 路径已失效)
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(openxlsx); library(readxl)})


## Package-relative paths (do not hard-code absolute workstation paths)
.args_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.code_dir <- if (length(.args_file)) dirname(normalizePath(sub("^--file=", "", .args_file[1]))) else getwd()
source(file.path(.code_dir, "_paths.R"))

S8_IN  <- file.path(DATA_BENCH, "benchmark_tbindex.xlsx")
S8     <- file.path(OUT, "benchmark_tbindex.xlsx")
file.copy(S8_IN, S8, overwrite = TRUE)
LIN_QN <- c(
  IDC = file.path(DATA_BENCH, "linear_qn_IDC.rds"),
  ITC = file.path(DATA_BENCH, "linear_qn_ITC.rds"),
  PVC = file.path(DATA_BENCH, "linear_qn_PVC.rds")
)

# 读取 QN 矩阵: rownames = 基因, colnames = 样本
load_qn <- function(path) { x <- readRDS(path); rownames(x$logTPM_QN) <- x$genes; x$logTPM_QN }
qn_list <- list(
  IDC = load_qn(LIN_QN[["IDC"]]),
  ITC = load_qn(LIN_QN[["ITC"]]),
  PVC = load_qn(LIN_QN[["PVC"]]))

ps <- as.data.table(read_excel(S8, sheet = "per_sample_scores"))

# 逐队列提取三个基因表达; 合并键为 (Cohort, Sample_ID) —— 因为 IDC/ITC 都有 TB_xx 命名
gene_dt <- rbindlist(lapply(c("IDC", "ITC", "PVC"), function(c) {
  sid <- ps[Cohort == c, Sample_ID]
  col <- sub("^LTBI_", "LTB_", sid)               # PVC 的 LTBI_ 前缀按 build 惯例映射为 LTB_
  m   <- match(col, colnames(qn_list[[c]]))
  ok  <- !is.na(m)
  dt  <- data.table(Sample_ID = sid, Cohort = c,
                    GBP5 = NA_real_, DUSP3 = NA_real_, KLF2 = NA_real_)
  if (any(ok)) dt[ok, `:=`(
    GBP5  = as.numeric(qn_list[[c]]["GBP5",  m[ok]]),
    DUSP3 = as.numeric(qn_list[[c]]["DUSP3", m[ok]]),
    KLF2  = as.numeric(qn_list[[c]]["KLF2",  m[ok]]))]
  dt
}))

ps2 <- merge(ps, gene_dt, by = c("Cohort", "Sample_ID"), all.x = TRUE, sort = FALSE)

# 校验: Sweeney3 是否等于 (GBP5 + DUSP3)/2 - KLF2 (用未舍入值比对)
ps2[, s3_check := (GBP5 + DUSP3) / 2 - KLF2]
d <- max(abs(ps2$s3_check - ps2$Sweeney3), na.rm = TRUE)
cat(sprintf("Sweeney3 vs (GBP5+DUSP3)/2-KLF2 最大差异: %.3e\n", d))
stopifnot(d < 0.001)   # 全量精度下应当几乎为 0 (差异仅来自 xlsx 存 4 位小数)

# 基因表达保留 4 位小数 (与 Sweeney3 精度风格一致)
for (g in c("GBP5", "DUSP3", "KLF2")) ps2[, (g) := round(get(g), 4)]
ps2[, s3_check := NULL]

# 列顺序与排序
setcolorder(ps2, c("Sample_ID", "Cohort", "Group", "circRNA_TB_index",
                   "GBP5", "DUSP3", "KLF2", "Sweeney3", "RISK6", "TimePoint"))
setorder(ps2, Cohort, Group, Sample_ID)

# 写回 xlsx (仅重建 per_sample_scores, 保留其他 sheets)
wb <- loadWorkbook(S8)
removeWorksheet(wb, "per_sample_scores")
addWorksheet(wb, "per_sample_scores")
writeData(wb, "per_sample_scores", ps2)
setColWidths(wb, "per_sample_scores", cols = 1:ncol(ps2), widths = "auto")
writeComment(wb, "per_sample_scores", col = 1, row = 1,
  comment = createComment("Per-sample benchmark scores (the literature-benchmark results). GBP5/DUSP3/KLF2 = linear-RNA QN log2(TPM+1); Sweeney3 = (GBP5+DUSP3)/2 - KLF2; RISK6 = GM(GBP2,FCGR1B,SERPING1)-GM(TUBGCP6,TRMT2A,SDR39U1) (not original qRT-PCR RISK6)."))
nms <- names(wb); rd <- which(nms == "README")
ord <- c(rd, which(nms == "per_sample_scores"), setdiff(seq_along(nms), c(rd, which(nms == "per_sample_scores"))))
worksheetOrder(wb) <- ord
saveWorkbook(wb, S8, overwrite = TRUE)
cat(sprintf("Done. per_sample_scores 已更新 (含三基因表达列), 共 %d 样本\n", nrow(ps2)))
