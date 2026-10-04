#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# Purpose:      Expression matrices and differential expression for the published
#               signature gene sets Sweeney3 and RISK6.
# Inputs:       Linear-RNA QN expression matrix.
# Outputs:      Per-gene expression sheets for the benchmark workbook.
# Key params:   Sweeney3: GBP5, DUSP3 (up), KLF2 (down); RISK6: GBP2, FCGR1B, SERPING1
#               (up), SDR39U1, TUBGCP6, TRMT2A (down).
# Dependencies: R >= 4.4.3; data.table; limma; openxlsx.
# ---------------------------------------------------------------------------
# =============================================================================
# build_S8_RISK6_genes.R — 文献评分基因(Sweeney3 / RISK6)的表达矩阵 + ATB 差异表达
#
# 对两个已发表评分基因集分别生成两个 sheet:
#   Sweeney3: GBP5, DUSP3 上调; KLF2 下调
#   RISK6   : GBP2, FCGR1B, SERPING1 上调; TUBGCP6, TRMT2A, SDR39U1 下调
#
# 每个基因集输出:
#   Sheet 1: {Signature}_genes_expression — 每样本基因 mRNA QN 表达 (log2(TPM+1))
#   Sheet 2: {Signature}_DE_ATB — IDC 中 ATB vs 其他四组 (HC/LTBI/PNEU/LC) 的差异
#            (log2FC, Welch t 检验 P 值, 上下调方向), 与文献方向对比,
#            不一致行标浅红 (#FFC7CE)
#
# 数据: linearRNA QN 矩阵 (data/benchmark/linear_qn_{IDC,ITC,PVC}.rds)
# 样本: 沿用 per_sample_scores 的样本 (IDC 109 / ITC 35 / PVC 81)
# =============================================================================
suppressPackageStartupMessages({library(data.table); library(openxlsx); library(readxl)})


## Package-relative paths (do not hard-code absolute workstation paths)
.args_file <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
.code_dir <- if (length(.args_file)) dirname(normalizePath(sub("^--file=", "", .args_file[1]))) else getwd()
source(file.path(.code_dir, "_paths.R"))

S8_IN  <- file.path(DATA_BENCH, "benchmark_tbindex.xlsx")
S8     <- file.path(OUT, "benchmark_tbindex.xlsx")  # write outputs here
file.copy(S8_IN, S8, overwrite = TRUE)
LIN_QN <- c(
  IDC = file.path(DATA_BENCH, "linear_qn_IDC.rds"),
  ITC = file.path(DATA_BENCH, "linear_qn_ITC.rds"),
  PVC = file.path(DATA_BENCH, "linear_qn_PVC.rds")
)
LIN_VS_CIRC <- file.path(DATA_BENCH, "linear_vs_circRNA.xlsx")

# 两个评分基因集的定义 (文献上下调方向) + RNA 类型/方向来源
GENE_SETS <- list(
  Sweeney3 = list(up = c("GBP5", "DUSP3"), down = c("KLF2"),
                  rna = "linear", source = "literature"),
  RISK6    = list(up = c("GBP2", "FCGR1B", "SERPING1"),
                  down = c("TUBGCP6", "TRMT2A", "SDR39U1"),
                  rna = "linear", source = "literature")
)

# circRNA TB-index 三基因 (dirY4 推导方向, 均为 circular)
CIRC_GENES <- c("GOLPH3", "SOD2", "PCNX2")
CIRC_DIR   <- c(GOLPH3 = "down", SOD2 = "up", PCNX2 = "down")   # 对应权重 -1/+1/-1

cat("===== 读取 linearRNA QN 矩阵与 per_sample_scores =====\n")
load_qn <- function(path) { x <- readRDS(path); rownames(x$logTPM_QN) <- x$genes; x$logTPM_QN }
qn_list <- list(
  IDC = load_qn(LIN_QN[["IDC"]]),
  ITC = load_qn(LIN_QN[["ITC"]]),
  PVC = load_qn(LIN_QN[["PVC"]]))
pss <- as.data.table(read_excel(S8, sheet = "per_sample_scores"))

# ---------------------------------------------------------------------------
# 提取某基因集的每样本表达 (双键合并, 因 IDC/ITC 都有 TB_xx 命名)
# ---------------------------------------------------------------------------
extract_expr <- function(c, sid, genes) {
  qn  <- qn_list[[c]]
  col <- sub("^LTBI_", "LTB_", sid)               # PVC 的 LTBI_ → LTB_ 映射
  m   <- match(col, colnames(qn)); ok <- !is.na(m)
  dt  <- data.table(Sample_ID = sid, Cohort = c)
  for (g in genes) dt[[g]] <- NA_real_
  if (any(ok)) for (g in genes) dt[ok, (g) := as.numeric(qn[g, m[ok]])]
  dt
}

wb <- loadWorkbook(S8)
de_summary <- list()   # 收集各基因集的 DE, 用于汇总 sheet

for (sig in names(GENE_SETS)) {
  up   <- GENE_SETS[[sig]]$up
  down <- GENE_SETS[[sig]]$down
  genes <- c(up, down)

  cat(sprintf("处理 %s ...\n", sig))

  # ---- Sheet 1: 表达矩阵 ----
  expr_dt <- rbindlist(lapply(c("IDC","ITC","PVC"), function(c) extract_expr(c, pss[Cohort==c, Sample_ID], genes)))
  expr_dt <- merge(pss[, .(Sample_ID, Cohort, Group, TimePoint)], expr_dt, by = c("Sample_ID","Cohort"))
  setcolorder(expr_dt, c("Sample_ID","Cohort","Group","TimePoint", genes))
  for (g in genes) expr_dt[, (g) := round(get(g), 4)]
  setorder(expr_dt, Cohort, Group, Sample_ID)
  sheet_expr <- paste0(sig, "_genes_expression")
  if (sheet_expr %in% names(wb)) removeWorksheet(wb, sheet_expr)
  addWorksheet(wb, sheet_expr)
  writeData(wb, sheet_expr, expr_dt)
  setColWidths(wb, sheet_expr, cols = 1:ncol(expr_dt), widths = "auto")
  writeComment(wb, sheet_expr, col = 1, row = 1,
    comment = createComment(sprintf("%s 基因 mRNA QN 表达 (log2(TPM+1)). Genes: %s. Source: linearRNA QN matrix.",
                                    sig, paste(genes, collapse=", "))))

  # ---- Sheet 2: 差异表达 (IDC, ATB vs 其他四组) + 文献方向对比 ----
  idc_expr <- expr_dt[Cohort == "IDC"]
  de_rows <- list()
  for (g in genes) {
    for (oth in c("HC", "LTBI", "PNEU", "LC")) {
      a <- idc_expr[Group == "APTB", ..g][[1]]
      b <- idc_expr[Group == oth,    ..g][[1]]
      if (length(a) < 2 || length(b) < 2) next
      l2fc <- mean(a) - mean(b)                        # log2FC = ATB 均值 − 其他组均值
      p    <- t.test(a, b)$p.value                     # Welch t 检验
      dir  <- fifelse(l2fc > 0, "up", "down")
      lit  <- fifelse(g %in% up, "up", "down")
      # P 值完整显示: < 0.0001 用科学计数法 (避免舍入成 0.0000)
      p_txt <- fifelse(p < 0.0001, formatC(p, format = "e", digits = 2),
                       formatC(round(p, 4), format = "f", digits = 4))
      de_rows[[length(de_rows) + 1]] <- data.table(
        Gene = g, Comparison = paste0("APTB vs ", oth),
        log2FC = round(l2fc, 4), P_value = p_txt,
        Direction = dir, Literature = lit,
        Consistent = fifelse(dir == lit, "Yes", "No"))
    }
  }
  de_dt <- rbindlist(de_rows)
  de_dt[, Comparison := factor(Comparison, c("APTB vs HC","APTB vs LTBI","APTB vs PNEU","APTB vs LC"))]
  de_dt[, Gene := factor(Gene, genes)]
  setorder(de_dt, Comparison, Gene)
  # 标注签名与 RNA 类型, 收集用于汇总 sheet
  de_dt[, `:=`(Signature = sig, RNA_type = GENE_SETS[[sig]]$rna,
               Direction_source = GENE_SETS[[sig]]$source)]
  de_summary[[sig]] <- copy(de_dt)

  sheet_de <- paste0(sig, "_DE_ATB")
  if (sheet_de %in% names(wb)) removeWorksheet(wb, sheet_de)
  addWorksheet(wb, sheet_de)
  writeData(wb, sheet_de, de_dt)
  setColWidths(wb, sheet_de, cols = 1:ncol(de_dt), widths = "auto")
  # 颜色标记: 方向与文献不一致的行整行标浅红
  no_rows <- which(de_dt$Consistent == "No") + 1L
  if (length(no_rows) > 0) {
    st_bad <- createStyle(fgFill = "#FFC7CE", fontColour = "#9C0006")
    for (r in no_rows) addStyle(wb, sheet_de, st_bad, rows = r, cols = 1:ncol(de_dt), gridExpand = TRUE)
  }
  writeComment(wb, sheet_de, col = 1, row = 1,
    comment = createComment(sprintf(
      "IDC cohort. log2FC = mean(ATB) - mean(other) on log2 QN expr; P = Welch t-test (full precision). Direction = our data up/down; Literature = up for %s, down for %s; Consistent = match. Red fill = mismatch.",
      paste(up, collapse="/"), paste(down, collapse="/"))))
  cat(sprintf("  %s: 表达 %d 样本, DE %d 行, 不一致 %d 行标红\n",
              sig, nrow(expr_dt), nrow(de_dt), length(no_rows)))
}

# =============================================================================
# 汇总: Gene_direction_summary — 全部 9 基因 (3 circ + 3 Sween + 6 RISK6)
# =============================================================================
cat("构建 Gene_direction_summary ...\n")

# circRNA 三基因的 IDC DE (来自 S7 xlsx 的 circRNA QN 表达), 方向来源 = dirY4
s7 <- as.data.table(read_excel(LIN_VS_CIRC,
                               sheet = "Per_sample_data"))
idc_circ <- s7[Cohort == "IDC"]
circ_rows <- list()
for (g in CIRC_GENES) {
  col <- paste0("circ_", g)
  for (oth in c("HC", "LTBI", "PNEU", "LC")) {
    a <- idc_circ[Group == "APTB", ..col][[1]]
    b <- idc_circ[Group == oth,    ..col][[1]]
    if (length(a) < 2 || length(b) < 2) next
    l2fc <- mean(a) - mean(b)
    p    <- t.test(a, b)$p.value
    dir  <- fifelse(l2fc > 0, "up", "down")
    exp  <- CIRC_DIR[[g]]
    p_txt <- fifelse(p < 0.0001, formatC(p, format = "e", digits = 2),
                     formatC(round(p, 4), format = "f", digits = 4))
    circ_rows[[length(circ_rows) + 1]] <- data.table(
      Signature = "circRNA TB-index", Gene = g, RNA_type = "circular",
      Direction_source = "dirY4", Expected_direction = exp,
      Comparison = paste0("APTB vs ", oth), log2FC = round(l2fc, 4), P_value = p_txt,
      Data_direction = dir, Consistent = fifelse(dir == exp, "Yes", "No"))
  }
}
circ_sum <- rbindlist(circ_rows)

# 统一线性基因集 (Sween/RISK6) 为同构列
lin_sum <- rbindlist(lapply(de_summary, function(d) {
  d[, .(Signature, Gene, RNA_type, Direction_source,
        Expected_direction = Literature,
        Comparison = as.character(Comparison),
        log2FC = as.numeric(log2FC), P_value = as.character(P_value),
        Data_direction = Direction, Consistent)]
}))

sum_dt <- rbind(lin_sum, circ_sum)
# 固定基因与比较的顺序
gene_order <- c(CIRC_GENES, GENE_SETS$Sweeney3$up, GENE_SETS$Sweeney3$down,
                GENE_SETS$RISK6$up, GENE_SETS$RISK6$down)
sum_dt[, Gene := factor(Gene, gene_order)]
sum_dt[, Comparison := factor(Comparison, c("APTB vs HC","APTB vs LTBI","APTB vs PNEU","APTB vs LC"))]
setorder(sum_dt, Comparison, Gene)
sum_dt[, Gene := as.character(Gene)]

# 拆成两个表: circular (circRNA TB-index 基因) 与 linear (Sweeney3 + RISK6 基因)
# 清理旧的合并表 (早期版本残留)
if ("Gene_direction_summary" %in% names(wb)) removeWorksheet(wb, "Gene_direction_summary")
# drop_consistent = TRUE 时写出前去掉 Consistent 列 (标红逻辑仍先用其计算行号)
write_gene_tab <- function(dt, sheet_name, comment, drop_consistent = FALSE) {
  if (sheet_name %in% names(wb)) removeWorksheet(wb, sheet_name)
  addWorksheet(wb, sheet_name)
  no_rows <- which(dt$Consistent == "No") + 1L
  dt_out <- if (drop_consistent) dt[, -"Consistent"] else dt
  writeData(wb, sheet_name, dt_out)
  setColWidths(wb, sheet_name, cols = 1:ncol(dt_out), widths = "auto")
  if (length(no_rows) > 0) {
    st_bad <- createStyle(fgFill = "#FFC7CE", fontColour = "#9C0006")
    for (r in no_rows) addStyle(wb, sheet_name, st_bad, rows = r,
                                cols = 1:ncol(dt_out), gridExpand = TRUE)
  }
  writeComment(wb, sheet_name, col = 1, row = 1, comment = createComment(comment))
  cat(sprintf("  %s: %d 行, 不一致 %d 行标红\n", sheet_name, nrow(dt_out), length(no_rows)))
}

# circular 表: circRNA TB-index 三基因 (dirY4 方向), 去掉 Consistent 列
circ_dt <- sum_dt[RNA_type == "circular"]
write_gene_tab(circ_dt, "Gene_direction_circular",
  paste("circRNA TB-index 三基因 (GOLPH3/SOD2/PCNX2), circular RNA 表达.",
        "Direction_source = dirY4 (数据推导, 权重 -1/+1/-1).",
        "log2FC = mean(ATB) - mean(other) on circular QN expr. P = Welch t-test."),
  drop_consistent = TRUE)

# linear 表: Sweeney3 + RISK6 九基因 (文献方向), 按 Signature 排序
lin_dt <- sum_dt[RNA_type == "linear"]
lin_dt[, Signature := factor(Signature, c("Sweeney3", "RISK6"))]
setorder(lin_dt, Signature, Gene, Comparison)
lin_dt[, Signature := as.character(Signature)]
write_gene_tab(lin_dt, "Gene_direction_linear",
  paste("Sweeney3 (GBP5/DUSP3/KLF2) 与 RISK6 (6 genes), linear mRNA 表达.",
        "Direction_source = literature (Sweeney3: GBP5/DUSP3 up, KLF2 down; RISK6: GBP2/FCGR1B/SERPING1 up, TUBGCP6/TRMT2A/SDR39U1 down).",
        "log2FC = mean(ATB) - mean(other) on linear QN expr. P = Welch t-test. Red fill = mismatch."))

# 新 sheet 放到 per_sample_scores 之后
nms <- names(wb); base <- which(nms == "per_sample_scores")
new_sheets <- c("Gene_direction_circular", "Gene_direction_linear",
                unlist(lapply(names(GENE_SETS), function(s) c(paste0(s, "_genes_expression"), paste0(s, "_DE_ATB")))))
ord <- c(base, which(nms %in% new_sheets), setdiff(seq_along(nms), c(base, which(nms %in% new_sheets))))
worksheetOrder(wb) <- ord

saveWorkbook(wb, S8, overwrite = TRUE)
cat(sprintf("Done → %s (%d sheets)\n", S8, length(names(wb))))
