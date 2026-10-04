#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Purpose:      Literature-comparison data for the benchmark (Supplementary Table S7A).
# Inputs:       Curated summary of previously published circRNA-based TB biomarker studies.
# Outputs:      Literature comparison sheet.
# Key params:   Descriptive comparison only; no re-analysis of published datasets.
# Dependencies: R >= 4.4.3; openxlsx.
# ---------------------------------------------------------------------------
# =============================================================================
# build_S8_Benchmark_literature_data.R — S8 文献公式版数据构建
#
# 背景: Figure_benchmark_TBindex.* 主版本自 2026-08-08 起将
#   Sweeney3/RISK6 改为 TB-index 式加权和 (权重 = IDC sign(logFC), 弃用文献
#   方向)。本脚本生成"文献公式版" Figure_benchmark_literature.xlsx,
#   在同一批 184 样本上按文献原始公式重算 Sweeney3/RISK6:
#
#   Sweeney3 (Sweeney et al. Lancet Respir Med 2016):
#       = (lin_GBP5 + lin_DUSP3)/2 − lin_KLF2
#   RISK6 (Penn-Nicholson et al. Sci Rep 2020):
#       = ((lin_GBP2 − lin_SDR39U1) + (lin_FCGR1B − lin_TUBGCP6)
#          + (lin_SERPING1 − lin_TRMT2A))/3
#
# 样本集与主版本完全一致 (per_sample_scores 的 184 样本, 含 9 个 linearRNA
#   QN log2(TPM+1) 表达列), circRNA TB-index 不变。仅替换 Sweeney3/RISK6 两列
#   及 Signature_metadata 的元数据; 基因表达 / 方向 sheet 不变。
#
# 输出: output/Figure_benchmark_literature.xlsx
#   (数据 sheet; AUC/Joint/DeLong 结果由 Figure_benchmark_literature.R
#    生成)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
})


## Package-relative paths (do not hard-code absolute workstation paths)
.args_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.code_dir <- if (length(.args_file)) dirname(normalizePath(sub("^--file=", "", .args_file[1]))) else getwd()
source(file.path(.code_dir, "_paths.R"))

S8       <- file.path(DATA_BENCH, "benchmark_tbindex.xlsx")
OUT_FILE <- file.path(OUT, "benchmark_literature.xlsx")

stopifnot(file.exists(S8))

cat("================================================================\n")
cat(" 文献公式版 S8: Sweeney3/RISK6 按文献原式重算\n")
cat("================================================================\n")

# ---- 0. 复制主版本工作簿作为底稿 ----
file.copy(S8, OUT_FILE, overwrite = TRUE)
wb <- loadWorkbook(OUT_FILE)
pss  <- as.data.table(read.xlsx(OUT, sheet = "per_sample_scores"))
ss   <- as.data.table(read.xlsx(OUT, sheet = "Sample_scores"))

# ---- 1. 文献公式重算 Sweeney3 / RISK6 (同一批 184 样本) ----
pss[, Sweeney3_lit := (lin_GBP5 + lin_DUSP3)/2 - lin_KLF2]
pss[, RISK6_lit := ((lin_GBP2 - lin_SDR39U1) +
                    (lin_FCGR1B - lin_TUBGCP6) +
                    (lin_SERPING1 - lin_TRMT2A))/3]
pss[, Sweeney3 := round(Sweeney3_lit, 4)]
pss[, RISK6 := round(RISK6_lit, 4)]
pss[, `:=`(Sweeney3_lit = NULL, RISK6_lit = NULL)]

# ---- 2. Sample_scores 同步 (circRNA_TB_index 不变, 仅换 Sweeney3/RISK6) ----
m <- merge(ss[, .(Sample_ID, Cohort)], pss[, .(Sample_ID, Cohort, Sweeney3, RISK6)],
           by = c("Sample_ID", "Cohort"), all.x = TRUE)
ss[, Sweeney3 := m$Sweeney3]
ss[, RISK6 := m$RISK6]

# ---- 3. 重写 per_sample_scores 与 Sample_scores ----
removeWorksheet(wb, "per_sample_scores"); addWorksheet(wb, "per_sample_scores")
writeData(wb, "per_sample_scores", pss)
setColWidths(wb, "per_sample_scores", cols = 1:ncol(pss), widths = "auto")
writeComment(wb, "per_sample_scores", col = 1, row = 1,
  comment = createComment(
    "Literature-formula scores (same 184 samples as benchmark version, 2026-08-10):\nSweeney3 = (GBP5+DUSP3)/2 - KLF2; RISK6 = ((GBP2-SDR39U1)+(FCGR1B-TUBGCP6)+(SERPING1-TRMT2A))/3\n(linearRNA QN log2(TPM+1)); circRNA TB-index unchanged."))

removeWorksheet(wb, "Sample_scores"); addWorksheet(wb, "Sample_scores")
writeData(wb, "Sample_scores", ss)
setColWidths(wb, "Sample_scores", cols = 1:ncol(ss), widths = "auto")

# ---- 4. Signature_metadata: 明确文献公式版 + 修正 RISK6 引用 (Penn-Nicholson) ----
sm <- as.data.table(read.xlsx(OUT, sheet = "Signature_metadata"))
sm[Signature == "Sweeney3", Formula := "(GBP5 + DUSP3)/2 - KLF2 (linear log2TPM QN)"]
sm[Signature == "RISK6", Formula := "((GBP2-SDR39U1)+(FCGR1B-TUBGCP6)+(SERPING1-TRMT2A))/3 (linear log2TPM QN)"]
sm[Signature == "RISK6", Reference := "Penn-Nicholson et al. Sci Rep 2020"]
sm[Signature != "circRNA TB-index", Reference := paste0(Reference, " (literature formula)")]
removeWorksheet(wb, "Signature_metadata"); addWorksheet(wb, "Signature_metadata")
writeData(wb, "Signature_metadata", sm)
setColWidths(wb, "Signature_metadata", cols = 1:ncol(sm), widths = "auto")

# ---- 5. 移除旧的 AUC/Joint/DeLong 结果 (由图脚本重新生成) ----
for (s in c("AUC_benchmark", "Delong_test", "Joint_model", "Joint_model_DeLong")) {
  if (s %in% names(wb)) removeWorksheet(wb, s)
}

# ---- 6. sheet 顺序对齐主版本 (去掉结果 sheet) ----
nms <- names(wb)
wanted <- c("per_sample_scores", "Sample_scores", "Gene_direction_circular", "Gene_direction_linear",
            "Sweeney3_genes_expression", "Sweeney3_DE_ATB",
            "RISK6_genes_expression", "RISK6_DE_ATB",
            "Expression_matrix_3circ_9lin",
            "Signature_metadata")
ord <- unlist(lapply(wanted, function(s) which(nms == s)))
worksheetOrder(wb) <- ord

saveWorkbook(wb, OUT_FILE, overwrite = TRUE)
cat(sprintf("\nDone → %s (%d sheets)\n", OUT, length(names(wb))))
cat("  Sweeney3 = (GBP5+DUSP3)/2 - KLF2 (literature)\n")
cat("  RISK6    = ((GBP2-SDR39U1)+(FCGR1B-TUBGCP6)+(SERPING1-TRMT2A))/3 (literature)\n")
cat("  circRNA TB-index 不变, 样本集与主版本一致 (184 samples)\n")
cat("================================================================\n")
