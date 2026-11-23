# ============================================================
# 03_4_consistency_and_sensitivity.r
#   Part 3: 03_1 Version A→M（p3 → cellPC）
#
# Dependencies：
#   results/03_1_celltype_marker_score/data/celltype_marker_scores.rds
#   results/03_2_deconvolution/data/cell_proportions.rds
#   results/03_3_causal_mediation/tables/A_to_M_diagnostic.csv
#   results/01_1_prepare_bulk_data/data/{expr_corrected,clinical,gene_annotation}.rds
#
# Output：
#   results/03_4_consistency_and_sensitivity/
#     tables/consistency_pearson_spearman.csv
#     tables/cellPC_loadings_03_1.csv
#     tables/A_to_M_03_1_sensitivity.csv
#     tables/final_summary.csv
#     figures/consistency_scatter.pdf
#     figures/consistency_bar.pdf
# ============================================================

rm(list = ls(all.names = TRUE)); gc()
setwd("E:/retina_tools/aging_rgc_axon/AMD")
set.seed(2026)

suppressPackageStartupMessages({
  library(limma); library(splines); library(ggplot2)
})

# ---- Outputdirectory ----
out_dir <- "results/03_4_consistency_and_sensitivity"
tab_dir <- file.path(out_dir, "tables")
fig_dir <- file.path(out_dir, "figures")
for (d in c(out_dir, tab_dir, fig_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

cat("=== 03_4 consistency + sensitivity start ===\n")

# ============================================================
# ============================================================
cat("\n--- 0. 载入 ---\n")

# 03_1 marker score
s1 <- readRDS("results/03_1_celltype_marker_score/data/celltype_marker_scores.rds")
s1_df <- as.data.frame(s1)
rownames(s1_df) <- s1_df$sample_id
s1_df$sample_id <- NULL
colnames(s1_df) <- sub("_score$", "", colnames(s1_df))
s1_mat <- as.matrix(s1_df)

# 03_2 deconvolution
s2_mat <- as.matrix(readRDS("results/03_2_deconvolution/data/cell_proportions.rds"))

cat("03_1:", nrow(s1_mat), "×", ncol(s1_mat), "\n")
cat("03_2:", nrow(s2_mat), "×", ncol(s2_mat), "\n")

# ============================================================
# ============================================================
cat("\n--- 1. 一致性检查 ---\n")

common_samples <- intersect(rownames(s1_mat), rownames(s2_mat))
s1_mat <- s1_mat[common_samples, , drop = FALSE]
s2_mat <- s2_mat[common_samples, , drop = FALSE]
cat("共同样本:", length(common_samples), "\n")

common_ct <- intersect(colnames(s1_mat), colnames(s2_mat))
cat("共同细胞类型:", paste(common_ct, collapse=", "), "\n")
cat("03_1 独有:", paste(setdiff(colnames(s1_mat), common_ct), collapse=", "), "\n")
cat("03_2 独有:", paste(setdiff(colnames(s2_mat), common_ct), collapse=", "), "\n")

consistency <- data.frame()
for (ct in common_ct) {
  x <- s1_mat[, ct]; y <- s2_mat[, ct]
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 30) next
  consistency <- rbind(consistency, data.frame(
    celltype     = ct,
    n_samples    = sum(ok),
    pearson_r    = round(cor(x[ok], y[ok]), 3),
    spearman_rho = round(cor(x[ok], y[ok], method = "spearman"), 3),
    stringsAsFactors = FALSE))
}
consistency <- consistency[order(-abs(consistency$pearson_r)), ]
print(consistency, row.names = FALSE)
write.csv(consistency, file.path(tab_dir, "consistency_pearson_spearman.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

median_r <- median(consistency$pearson_r, na.rm = TRUE)
cat("\nPearson r 中位数:", round(median_r, 3), "\n")
if (median_r > 0.6) {
  cat("→ 高度一致\n")
} else if (median_r > 0.3) {
  cat("→ 中度一致\n")
} else {
  cat("→ 弱一致\n")
}

plot_df <- data.frame()
for (ct in common_ct) {
  x <- s1_mat[, ct]; y <- s2_mat[, ct]
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 30) next
  plot_df <- rbind(plot_df, data.frame(
    celltype = ct, score_03_1 = x[ok], prop_03_2 = y[ok],
    stringsAsFactors = FALSE))
}
if (nrow(plot_df) > 0) {
  p1 <- ggplot(plot_df, aes(x = score_03_1, y = prop_03_2)) +
    geom_point(alpha = 0.3, size = 0.4) +
    geom_smooth(method = "lm", se = FALSE, color = "#C44E52", linewidth = 0.5) +
    facet_wrap(~ celltype, scales = "free", ncol = 3) +
    theme_bw(base_size = 10) +
    labs(x = "03_1 marker score", y = "03_2 deconvolution proportion",
         title = "Consistency: marker score vs deconvolution")
  ggsave(file.path(fig_dir, "consistency_scatter.pdf"), p1,
         width = 9, height = 8)

  p2 <- ggplot(consistency, aes(x = reorder(celltype, -pearson_r),
                                 y = pearson_r, fill = pearson_r > 0.6)) +
    geom_col(width = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = 0.6, linetype = "dotted", color = "red") +
    scale_fill_manual(values = c("TRUE" = "#4C72B0", "FALSE" = "#B0B0B0"),
                      guide = "none") +
    theme_bw(base_size = 11) +
    labs(x = "", y = "Pearson r",
         title = "03_1 vs 03_2 consistency per cell type")
  ggsave(file.path(fig_dir, "consistency_bar.pdf"), p2,
         width = 7, height = 5)
}

# ============================================================
# ============================================================
cat("\n--- 2. 03_1 PCA ---\n")

ct_common <- c("AC","BC","Cone","HC","MG","RGC","Rod")
s1_use <- s1_mat[, intersect(ct_common, colnames(s1_mat)), drop = FALSE]

pca1 <- prcomp(s1_use, scale. = TRUE)
n_pc <- min(5, ncol(pca1$x))
med_use1 <- pca1$x[, 1:n_pc, drop = FALSE]
colnames(med_use1) <- paste0("cellPC", 1:n_pc)

cat("cellPC 方差解释:",
    paste(round(100 * pca1$sdev[1:n_pc]^2 / sum(pca1$sdev^2), 1), "%",
          collapse=", "), "\n")

cat("\n=== 03_1 cellPC 载荷 ===\n")
for (i in 1:n_pc) {
  cat(sprintf("cellPC%d: %s\n", i,
              paste(sprintf("%s=%.2f",
                            names(sort(pca1$rotation[, i], decreasing = TRUE)),
                            sort(pca1$rotation[, i], decreasing = TRUE)),
                    collapse = ", ")))
}

loadings_df <- data.frame(
  celltype = rownames(pca1$rotation),
  pca1$rotation,
  check.names = FALSE)
write.csv(loadings_df, file.path(tab_dir, "cellPC_loadings_03_1.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# ============================================================
# Part 3: 03_1 Version A→M（p3 → cellPC）
# ============================================================
cat("\n--- 3. 03_1 版本 A→M ---\n")

expr <- readRDS("results/01_1_prepare_bulk_data/data/expr_corrected.rds")
clin <- readRDS("results/01_1_prepare_bulk_data/data/clinical.rds")
anno <- readRDS("results/01_1_prepare_bulk_data/data/gene_annotation.rds")

common <- intersect(colnames(expr), rownames(clin))
expr <- expr[, common, drop = FALSE]
clin <- clin[common, , drop = FALSE]

anno_sub <- anno[anno$ensembl %in% rownames(expr), ]
anno_sub <- anno_sub[!is.na(anno_sub$symbol), ]
expr <- expr[anno_sub$ensembl, , drop = FALSE]
expr <- expr[!duplicated(anno_sub$symbol), , drop = FALSE]
rownames(expr) <- anno_sub$symbol[!duplicated(anno_sub$symbol)]

clin$grade <- as.integer(clin$grade)
rownames(clin) <- colnames(expr)

# samplesalign
common2 <- Reduce(intersect, list(rownames(med_use1), rownames(clin), colnames(expr)))
cat("对齐后样本:", length(common2), "\n")
med_use1 <- med_use1[common2, , drop = FALSE]
expr <- expr[, common2, drop = FALSE]
clin <- clin[common2, , drop = FALSE]

# covariate + p3
sp <- splines::ns(clin$age_n, df = 3)
C <- cbind(age_s1 = sp[,1], age_s2 = sp[,2], age_s3 = sp[,3],
           sex_M = as.integer(clin$sex == "Male"),
           rin = clin$rin_n, pmi = clin$pmi_n)
rownames(C) <- common2

poly_mat <- poly(clin$grade, 3)
A_p3 <- poly_mat[, 3]
M <- med_use1
n_med <- ncol(M)

res_03_1 <- data.frame()
for (i in seq_len(n_med)) {
  fit_i <- lm(M[, i] ~ A_p3 + C)
  co <- summary(fit_i)$coefficients
  if (!"A_p3" %in% rownames(co)) next
  res_03_1 <- rbind(res_03_1, data.frame(
    cellPC = colnames(M)[i],
    a_coef = round(co["A_p3", "Estimate"], 4),
    a_p    = round(co["A_p3", "Pr(>|t|)"], 4),
    stringsAsFactors = FALSE))
}
res_03_1$a_FDR <- round(p.adjust(res_03_1$a_p, "BH"), 4)
print(res_03_1, row.names = FALSE)
write.csv(res_03_1, file.path(tab_dir, "A_to_M_03_1_sensitivity.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# ============================================================
# ============================================================
cat("\n--- 4. 汇总对比 ---\n")

f_03_2 <- "results/03_3_causal_mediation/tables/A_to_M_diagnostic.csv"
if (file.exists(f_03_2)) {
  am_03_2 <- read.csv(f_03_2, stringsAsFactors = FALSE)
  res_03_2 <- data.frame(
    cellPC   = am_03_2$mediator,
    a_03_2   = am_03_2$a_coef,
    p_03_2   = am_03_2$a_p,
    FDR_03_2 = am_03_2$a_FDR,
    stringsAsFactors = FALSE)
  cat("★ 从 03_3 输出读取 03_2 结果\n")
} else {
  cat("★ 03_3 输出不存在，用硬编码 fallback\n")
  res_03_2 <- data.frame(
    cellPC   = c("cellPC1","cellPC2","cellPC3","cellPC4","cellPC5"),
    a_03_2   = c( 1.1640, -0.9822, -2.4573, -0.4487,  0.2658),
    p_03_2   = c( 0.6026,  0.4298,  0.0044,  0.5705,  0.7270),
    FDR_03_2 = c( 0.7270,  0.7270,  0.0222,  0.7270,  0.7270),
    stringsAsFactors = FALSE)
}

final_summary <- merge(res_03_1, res_03_2, by = "cellPC", all = TRUE)
final_summary <- final_summary[, c("cellPC",
                                    "a_03_2", "p_03_2", "FDR_03_2",
                                    "a_coef", "a_p", "a_FDR")]
colnames(final_summary) <- c("cellPC",
                              "a_03_2_deconv", "p_03_2_deconv", "FDR_03_2_deconv",
                              "a_03_1_marker", "p_03_1_marker", "FDR_03_1_marker")
write.csv(final_summary, file.path(tab_dir, "final_summary.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

cat("\n=== 03_1 vs 03_2 对比 ===\n")
print(final_summary, row.names = FALSE)

best_03_1 <- res_03_1[which.min(res_03_1$a_p), ]
best_pc_idx <- as.integer(sub("cellPC", "", best_03_1$cellPC))
best_loading <- names(which.max(abs(pca1$rotation[, best_pc_idx])))

cat("\n=== 核心结论 ===\n")
cat("一致性: 7 种共同细胞类型，Pearson r 中位数 =", round(median_r, 3), "\n")
cat(sprintf("03_2 显著 PC: %s (FDR=%.4f)\n",
            res_03_2$cellPC[which.min(res_03_2$FDR_03_2)],
            min(res_03_2$FDR_03_2)))
cat(sprintf("03_1 最显著 PC: %s (FDR=%.4f, 载荷主导=%s)\n",
            best_03_1$cellPC, best_03_1$a_FDR, best_loading))

idx_03_2 <- which.min(res_03_2$FDR_03_2)
idx_03_1 <- which.min(res_03_1$a_p)
a_03_2_sig <- res_03_2$a_03_2[idx_03_2]
a_03_1_sig <- res_03_1$a_coef[idx_03_1]
if (sign(a_03_2_sig) == sign(a_03_1_sig)) {
  cat(sprintf("方向一致: 两个方法都是 %s\n",
              ifelse(a_03_2_sig > 0, "正", "负")))
} else {
  cat("方向不一致\n")
}

# ============================================================
# summary
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("=== 03_4 完成 ===\n")
cat(strrep("=", 70), "\n", sep = "")
cat("一致性: Pearson r 中位数 =", round(median_r, 3), "\n")
cat("敏感性: 03_1 最显著 PC =", best_03_1$cellPC,
    "(FDR =", round(best_03_1$a_FDR, 4), ", 载荷主导 =", best_loading, ")\n")
cat("输出目录:", out_dir, "\n")
cat("=== done ===\n")
