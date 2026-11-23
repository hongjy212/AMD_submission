# ============================================================
# 02_5_mashr_shared_effects.r
# Input：
#   results/01_1_prepare_bulk_data/data/{expr_corrected,clinical,gene_annotation}.rds
#   results/02_1_differential_expression_main/tables/sig16_genes.csv
#   tables/mash_PM_all.csv
#   tables/mash_LFSR_all.csv
#   tables/mash_pi_est.csv
#   tables/mash_16_long.csv
#   tables/mash_16_logFC_wide.csv
#   tables/mash_16_lfsr_wide.csv
#   tables/mash_only_significant.csv
#   tables/mash_gene_pattern.csv
#   tables/mash_contrast_pattern.csv
#   figures/mash_lfsr_heatmap.pdf/.png
#   figures/mash_logFC_heatmap.pdf/.png
# ============================================================

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(limma); library(splines); library(mashr)
  library(tidyr); library(pheatmap)
})

DIR_THIS <- setup_output_dir("02_5_mashr_shared_effects")
log_info("=== 02_5_mashr_shared_effects.r start ===")

# ============================================================
# ============================================================
log_info("--- 1. 载入数据 ---")

expr <- read_upstream("01_1_prepare_bulk_data", "data", "expr_corrected.rds")
clin <- read_upstream("01_1_prepare_bulk_data", "data", "clinical.rds")
anno <- read_upstream("01_1_prepare_bulk_data", "data", "gene_annotation.rds")

sig16_file <- file.path(DIR_RESULT,
                         "02_1_differential_expression_main",
                         "tables", "sig16_genes.csv")
if (!file.exists(sig16_file)) {
  stop("请先运行 02_1_differential_expression_main.r", call. = FALSE)
}
sig16 <- read.csv(sig16_file, stringsAsFactors = FALSE)$gene
log_info("读入 16 基因: {length(sig16)} 个")

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

log_info("表达矩阵: {nrow(expr)} × {ncol(expr)}")

# ============================================================
# ============================================================
log_info("--- 2. limma multi (4 对比) ---")

sp <- splines::ns(clin$age_n, df = 3)
covar_df <- data.frame(
  age_s1 = sp[,1], age_s2 = sp[,2], age_s3 = sp[,3],
  sex_M = as.integer(clin$sex == "Male"),
  rin = clin$rin_n, pmi = clin$pmi_n
)
rownames(covar_df) <- rownames(clin)

design_df <- data.frame(grade = factor(clin$grade), covar_df)
rownames(design_df) <- rownames(clin)
design_df <- design_df[complete.cases(design_df), , drop = FALSE]
expr_sub <- expr[, rownames(design_df), drop = FALSE]

design_multi <- model.matrix(
  ~ 0 + grade + age_s1 + age_s2 + age_s3 + sex_M + rin + pmi,
  data = design_df)
colnames(design_multi)[grep("^grade", colnames(design_multi))] <- paste0("g", 1:4)

cm <- limma::makeContrasts(
  g2_vs_g1 = g2 - g1,
  g3_vs_g2 = g3 - g2,
  g4_vs_g3 = g4 - g3,
  g4_vs_g1 = g4 - g1,
  levels = design_multi)
fit_multi <- eBayes(contrasts.fit(lmFit(expr_sub, design_multi), cm))

# ============================================================
# 3. build mashr Input (Bhat / Shat)
# ============================================================
log_info("--- 3. 构建 Bhat / Shat ---")

mash_contrasts <- c("g2_vs_g1", "g3_vs_g2", "g4_vs_g3", "g4_vs_g1")
Bhat_m <- fit_multi$coefficients[, mash_contrasts, drop = FALSE]
sigma_post <- sqrt(fit_multi$s2.post)
Shat_m <- fit_multi$stdev.unscaled[, mash_contrasts, drop = FALSE] * sigma_post

ok <- complete.cases(Bhat_m, Shat_m) &
      apply(Shat_m, 1, function(x) all(x > 0 & is.finite(x)))
Bhat_m <- Bhat_m[ok, , drop = FALSE]
Shat_m <- Shat_m[ok, , drop = FALSE]
cn <- colnames(Bhat_m)
log_info("mash 输入: {nrow(Bhat_m)} × {ncol(Bhat_m)}")

# ============================================================
# ============================================================
log_info("--- 4. mashr 模型 ---")

mash_data <- mashr::mash_set_data(Bhat_m, Shat_m)

# Canonical（exclude equal_effects）
U.c_full <- mashr::cov_canonical(mash_data)
U.c <- U.c_full[!grepl("equal_effects", names(U.c_full))]
log_info("Canonical: {length(U.c)}")

# PCA (npc=2)
U.pca <- tryCatch(mashr::cov_pca(mash_data, npc = 2),
                   error = function(e) { log_warn("cov_pca 失败"); NULL })
U.ed <- if (!is.null(U.pca)) {
  tryCatch(mashr::cov_ed(mash_data, U.pca), error = function(e) NULL)
} else NULL

U.global <- diag(length(cn))
dimnames(U.global) <- list(cn, cn)
U.global["g4_vs_g1", "g4_vs_g1"] <- 1
U.global <- list(global_only = U.global)

U.adj <- diag(length(cn))
dimnames(U.adj) <- list(cn, cn)
U.adj["g2_vs_g1", "g2_vs_g1"] <- 1.0
U.adj["g3_vs_g2", "g3_vs_g2"] <- 1.0
U.adj["g4_vs_g3", "g4_vs_g3"] <- 1.0
U.adj["g2_vs_g1", "g3_vs_g2"] <- 0.3; U.adj["g3_vs_g2", "g2_vs_g1"] <- 0.3
U.adj["g3_vs_g2", "g4_vs_g3"] <- 0.3; U.adj["g4_vs_g3", "g3_vs_g2"] <- 0.3
U.adj <- list(adjacent_chain = U.adj)

Ulist <- c()
if (!is.null(U.ed)) Ulist <- c(Ulist, U.ed)
Ulist <- c(Ulist, U.global, U.adj, U.c)
log_info("Ulist: {length(Ulist)}")

log_info("训练 mash 模型（outputlevel=2）...")
t0 <- Sys.time()
m <- mashr::mash(mash_data, Ulist = Ulist, outputlevel = 2)
log_info("训练用时: {round(as.numeric(difftime(Sys.time(), t0, units='secs')), 1)}s")

# ============================================================
# 5. extractResults
# ============================================================
log_info("--- 5. 提取结果 ---")
log_info("m$result 字段: {paste(names(m[['result']]), collapse=', ')}")

PM <- m[["result"]][["PosteriorMean"]]
PSD <- m[["result"]][["PosteriorSD"]]
LFSR <- m[["result"]][["lfsr"]]
rownames(PM) <- rownames(Bhat_m); colnames(PM) <- cn
rownames(PSD) <- rownames(Bhat_m); colnames(PSD) <- cn
rownames(LFSR) <- rownames(Bhat_m); colnames(LFSR) <- cn

save_table(as.data.frame(PM), DIR_THIS, "mash_PM_all.csv")
save_table(as.data.frame(LFSR), DIR_THIS, "mash_LFSR_all.csv")

pi_est <- sort(mashr::get_estimated_pi(m), decreasing = TRUE)
save_table(data.frame(component = names(pi_est), pi = as.numeric(pi_est)),
           DIR_THIS, "mash_pi_est.csv")
cat("\n=== Pi 估计 top 5 ===\n")
print(head(pi_est, 5))

# ============================================================
# 6. extract 16 genes
# ============================================================
log_info("--- 6. 提取 16 基因 ---")

sig16_in <- intersect(sig16, rownames(LFSR))
log_info("16 基因匹配 mash: {length(sig16_in)}")

pm_16 <- PM[sig16_in, , drop = FALSE]
lfsr_16 <- LFSR[sig16_in, , drop = FALSE]

# limma FDR
limma_fdr <- data.frame()
for (ct in cn) {
  tt <- topTable(fit_multi, coef = ct, number = Inf, sort.by = "none")
  limma_fdr <- rbind(limma_fdr, data.frame(
    gene = rownames(tt), contrast = ct,
    limma_FDR = tt$adj.P.Val, stringsAsFactors = FALSE))
}

mash_long <- data.frame()
for (g in sig16_in) {
  for (ct in cn) {
    mash_long <- rbind(mash_long, data.frame(
      gene = g, contrast = ct,
      limma_logFC = fit_multi$coefficients[g, ct],
      limma_FDR = limma_fdr$limma_FDR[limma_fdr$gene == g &
                                       limma_fdr$contrast == ct],
      mash_logFC = pm_16[g, ct],
      mash_SD = PSD[g, ct],
      lfsr = lfsr_16[g, ct],
      stringsAsFactors = FALSE))
  }
}
mash_long$limma_sig <- mash_long$limma_FDR < 0.05
mash_long$mash_sig <- mash_long$lfsr < 0.05
save_table(mash_long, DIR_THIS, "mash_16_long.csv")

mash_wide_lfc <- tidyr::pivot_wider(
  mash_long[, c("gene","contrast","mash_logFC")],
  names_from = contrast, values_from = mash_logFC)
save_table(as.data.frame(mash_wide_lfc), DIR_THIS, "mash_16_logFC_wide.csv")

mash_wide_lfsr <- tidyr::pivot_wider(
  mash_long[, c("gene","contrast","lfsr")],
  names_from = contrast, values_from = lfsr)
save_table(as.data.frame(mash_wide_lfsr), DIR_THIS, "mash_16_lfsr_wide.csv")

# ============================================================
# ============================================================
log_info("--- 7. 对比统计 ---")

cross <- table(limma = mash_long$limma_sig, mash = mash_long$mash_sig)
cat("\n=== limma vs mash 交叉表 ===\n")
print(cross)
n_both <- sum(mash_long$limma_sig & mash_long$mash_sig)
n_mash_only <- sum(!mash_long$limma_sig & mash_long$mash_sig)
n_limma_only <- sum(mash_long$limma_sig & !mash_long$mash_sig)
log_info("两者都显著: {n_both}")
log_info("仅 mash 显著: {n_mash_only}")
log_info("仅 limma 显著: {n_limma_only}")

if (n_mash_only > 0) {
  mash_only <- mash_long[!mash_long$limma_sig & mash_long$mash_sig, ]
  save_table(mash_only, DIR_THIS, "mash_only_significant.csv")
  cat("\n=== mash 独有显著 top 20 ===\n")
  print(head(mash_only[order(mash_only$lfsr),
                        c("gene","contrast","limma_logFC","mash_logFC",
                          "lfsr","limma_FDR")], 20), row.names = FALSE)
}

gene_pattern <- aggregate(mash_sig ~ gene, data = mash_long, FUN = sum)
colnames(gene_pattern) <- c("gene", "n_sig_contrasts")
gene_pattern <- gene_pattern[order(-gene_pattern$n_sig_contrasts), ]
save_table(gene_pattern, DIR_THIS, "mash_gene_pattern.csv")

contrast_pattern <- aggregate(mash_sig ~ contrast, data = mash_long, FUN = sum)
colnames(contrast_pattern) <- c("contrast", "n_sig_genes")
contrast_pattern <- contrast_pattern[order(-contrast_pattern$n_sig_genes), ]
save_table(contrast_pattern, DIR_THIS, "mash_contrast_pattern.csv")

cat("\n=== 每基因显著对比数 ===\n"); print(gene_pattern, row.names = FALSE)
cat("\n=== 每对比显著基因数 ===\n"); print(contrast_pattern, row.names = FALSE)

# ============================================================
# ============================================================
log_info("--- 8. 图 ---")

lfsr_mat <- as.matrix(mash_wide_lfsr[, -1])
rownames(lfsr_mat) <- mash_wide_lfsr$gene
pheatmap(lfsr_mat, cluster_rows = TRUE, cluster_cols = TRUE,
         display_numbers = TRUE, number_format = "%.3f",
         color = colorRampPalette(c("#C44E52","#FFFFFF","#4C72B0"))(50),
         main = "mashr LFSR (16 genes x 4 contrasts)",
         filename = file.path(DIR_THIS, "figures", "mash_lfsr_heatmap.png"),
         width = 9, height = 7)

lfc_mat <- as.matrix(mash_wide_lfc[, -1])
rownames(lfc_mat) <- mash_wide_lfc$gene
lim <- max(abs(lfc_mat), na.rm = TRUE)
pheatmap(lfc_mat, cluster_rows = TRUE, cluster_cols = TRUE,
         display_numbers = TRUE, number_format = "%.2f",
         color = colorRampPalette(c("#4C72B0","#FFFFFF","#C44E52"))(50),
         breaks = seq(-lim, lim, length.out = 51),
         main = "mashr posterior logFC (16 genes x 4 contrasts)",
         filename = file.path(DIR_THIS, "figures", "mash_logFC_heatmap.png"),
         width = 9, height = 7)

# ============================================================
# 9. summary
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("=== 02_5 mashr 汇总 ===\n")
cat(strrep("=", 70), "\n", sep = "")
cat("16 基因 mash 分析:\n")
cat(sprintf("  两者都显著: %d / %d\n", n_both, nrow(mash_long)))
cat(sprintf("  仅 mash 显著: %d (%.1f%%)\n", n_mash_only,
            100 * n_mash_only / nrow(mash_long)))

cat("\n--- 每基因显著对比数 ---\n")
print(gene_pattern, row.names = FALSE)
cat("\n--- 每对比显著基因数 ---\n")
print(contrast_pattern, row.names = FALSE)

if (n_mash_only > 0) {
  cat("\n--- mash 独有显著 top 10 ---\n")
  print(head(mash_only[order(mash_only$lfsr),
                        c("gene","contrast","limma_logFC","mash_logFC","lfsr")], 10),
        row.names = FALSE)
}

cat("\n--- 输出位置 ---\n")
cat("目录:", DIR_THIS, "\n")
cat("  tables/mash_PM_all.csv\n")
cat("  tables/mash_LFSR_all.csv\n")
cat("  tables/mash_pi_est.csv\n")
cat("  tables/mash_16_long.csv\n")
cat("  tables/mash_16_logFC_wide.csv\n")
cat("  tables/mash_16_lfsr_wide.csv\n")
cat("  tables/mash_only_significant.csv\n")
cat("  tables/mash_gene_pattern.csv\n")
cat("  tables/mash_contrast_pattern.csv\n")
cat("  figures/mash_lfsr_heatmap.pdf/.png\n")
cat("  figures/mash_logFC_heatmap.pdf/.png\n")

log_info("=== 02_5_mashr_shared_effects.r done ===")
