#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Purpose:      Recompute the benchmark signature scores used in the manuscript.
# Inputs:       Benchmark expression matrix (circRNA and linear-RNA QN values).
# Outputs:      Per-sample benchmark scores (the literature-benchmark results).
# Key params:   Scores computed as signed weighted sums of QN expression values.
# Dependencies: R >= 4.4.3; data.table; openxlsx.
# ---------------------------------------------------------------------------
# =============================================================================
# build_benchmark_tbindex_scores.R — 把 benchmark 的 Sweeney3/RISK6 评分改为 TB-index 式加权和
#
# 需求（2026-08-08）：Panel C 中 Sweeney3/RISK6 不再用文献公式，改用 TB-index
#   的公式 = Σ(权重 × linearRNA 表达)，权重 = 由 IDC 数据 APTB vs non-APTB 的
#   limma 方向 (sign logFC) 决定，完全弃用文献方向。circRNA TB-index 不变。
#
# 数据来源: benchmark workbook 的 Expression_matrix_3circ_9lin（184 样本，12 基因 QN 表达）
# 输出: 更新 per_sample_scores sheet —— 替换 Sweeney3/RISK6 两列为新值，
#       并补充 9 个 linearRNA 基因表达列（自包含：表达 → 评分）
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(openxlsx)
})


## Package-relative paths
.args_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.code_dir <- if (length(.args_file)) dirname(normalizePath(sub("^--file=", "", .args_file[1]))) else getwd()
source(file.path(.code_dir, "_paths.R"))
S8 <- file.path(DATA_BENCH, "benchmark_tbindex.xlsx")

cat("================================================================\n")
cat(" 把 Sweeney3/RISK6 改为 TB-index 式加权和\n")
cat("================================================================\n")
stopifnot(file.exists(S8))

pss <- as.data.table(read.xlsx(S8, sheet = "per_sample_scores"))
ex  <- as.data.table(read.xlsx(S8, sheet = "Expression_matrix_3circ_9lin"))

# ---- 保留原 circRNA TB-index 与临床列 ----
base <- pss[, .(Sample_ID, Cohort, Group, circRNA_TB_index, TimePoint)]

# ---- 权重: IDC APTB vs non-APTB 的 limma 方向 (sign logFC)，完全弃用文献方向 ----
idcE <- ex[Cohort == "IDC"]
y <- fifelse(idcE$Group == "APTB", 1L, 0L)
LIN9 <- c("lin_GBP5", "lin_DUSP3", "lin_KLF2",
          "lin_GBP2", "lin_FCGR1B", "lin_SERPING1",
          "lin_TUBGCP6", "lin_TRMT2A", "lin_SDR39U1")
suppressPackageStartupMessages(library(limma))
M <- t(as.matrix(idcE[, LIN9, with = FALSE]))
design <- model.matrix(~factor(y, levels = c(0, 1)))
fit <- eBayes(lmFit(M, design))
tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
w <- setNames(sign(tt$logFC), rownames(tt))            # 权重 = sign(logFC)
cat("数据方向权重 (sign logFC, IDC APTB vs non-APTB):\n")
for (g in LIN9) cat(sprintf("  %-10s %+d\n", sub("lin_", "", g), w[g]))

# ---- 按 TB-index 公式重算 Sweeney3 / RISK6（Σ 权重×表达）----
ex[, Sweeney3 := Reduce(`+`, Map(function(g) w[g] * get(g),
                                 c("lin_GBP5", "lin_DUSP3", "lin_KLF2")))]
ex[, RISK6 := Reduce(`+`, Map(function(g) w[g] * get(g),
                              c("lin_GBP2", "lin_FCGR1B", "lin_SERPING1",
                                "lin_TUBGCP6", "lin_TRMT2A", "lin_SDR39U1")))]
ex[, Sweeney3 := round(Sweeney3, 4)]
ex[, RISK6 := round(RISK6, 4)]

# ---- 合并: 保留 9 个 linearRNA 基因表达列 + 新评分 ----
keep_lin <- c("lin_GBP5", "lin_DUSP3", "lin_KLF2",
              "lin_GBP2", "lin_FCGR1B", "lin_SERPING1",
              "lin_TUBGCP6", "lin_TRMT2A", "lin_SDR39U1")
ps2 <- merge(base, ex[, c("Sample_ID", "Cohort", keep_lin, "Sweeney3", "RISK6"), with = FALSE],
             by = c("Sample_ID", "Cohort"), all.x = TRUE)
setcolorder(ps2, c("Sample_ID", "Cohort", "Group", "circRNA_TB_index",
                   keep_lin, "Sweeney3", "RISK6", "TimePoint"))
setorder(ps2, Cohort, Group, Sample_ID)

cat(sprintf("per_sample_scores 更新: %d 样本\n", nrow(ps2)))
cat("  Sweeney3 = Σ(w×表达), w=sign(logFC) IDC APTB vs non-APTB\n")
cat("  RISK6    = Σ(w×表达), w=sign(logFC) IDC APTB vs non-APTB\n")
cat("  circRNA TB-index 不变\n")

wb <- loadWorkbook(S8)
removeWorksheet(wb, "per_sample_scores")
addWorksheet(wb, "per_sample_scores")
writeData(wb, "per_sample_scores", ps2)
setColWidths(wb, "per_sample_scores", cols = 1:ncol(ps2), widths = "auto")
writeComment(wb, "per_sample_scores", col = 1, row = 1,
  comment = createComment(
    "Scores recomputed 2026-08-08 as TB-index-style weighted sums:\nSweeney3 = +GBP5 + DUSP3 - KLF2; RISK6 = +GBP2 + FCGR1B + SERPING1 - SDR39U1 - TUBGCP6 - TRMT2A\n(w = up/down direction; linearRNA QN log2(TPM+1)); circRNA TB-index unchanged."))

# per_sample_scores 放 README 后（无 README 则放最前）
nms <- names(wb)
rd  <- which(nms == "per_sample_scores")
worksheetOrder(wb) <- c(rd, setdiff(seq_along(nms), rd))

saveWorkbook(wb, S8, overwrite = TRUE)
cat(sprintf("\nDone → %s (%d sheets)\n", S8, length(names(wb))))
cat("================================================================\n")
