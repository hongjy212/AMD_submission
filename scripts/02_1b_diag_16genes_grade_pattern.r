#!/usr/bin/env Rscript
# ============================================================
# diag_16genes_grade_pattern.r
#
# Background：
#
#
# Output（results/diag/diag_16genes_grade_pattern/）：
#   figures/genes_poly_terms_heatmap.pdf/.png
#   figures/genes_pattern_breakdown.pdf/.png
# ============================================================

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(limma); library(splines); library(ggplot2)
  library(pheatmap); library(tidyr)
})

DIR_THIS <- setup_output_dir("diag_16genes_grade_pattern", subdir = "diag")
log_info("=== diag_16genes_grade_pattern.r start ===")

# ============================================================
# ============================================================
log_info("--- 1. 载入 ---")

expr <- read_upstream("01_1_prepare_bulk_data", "data", "expr_corrected.rds")
clin <- read_upstream("01_1_prepare_bulk_data", "data", "clinical.rds")
anno <- read_upstream("01_1_prepare_bulk_data", "data", "gene_annotation.rds")

common <- intersect(colnames(expr), rownames(clin))
expr <- expr[, common, drop = FALSE]
clin <- clin[common, , drop = FALSE]

anno_sub <- anno[anno$ensembl %in% rownames(expr), ]
anno_sub <- anno_sub[!is.na(anno_sub$symbol), ]
expr <- expr[anno_sub$ensembl, , drop = FALSE]
expr <- expr[!duplicated(anno_sub$symbol), , drop = FALSE]
rownames(expr) <- toupper(anno_sub$symbol[!duplicated(anno_sub$symbol)])

clin$grade <- as.integer(clin$grade)
clin$age_v <- as.numeric(clin$age_n)
clin$pmi_v <- as.numeric(clin$pmi_n)
clin$rin_v <- as.numeric(clin$rin_n)
clin$batch_f <- factor(clin$batch)
rownames(clin) <- colnames(expr)

GENES_16 <- c("FCF1","UGGT1","CCDC3","S100B","BPNT1","BRI3BP","BTBD19",
              "DPH7","AKR1B1","ZNF623","ZNF587","RSPH3","TAF3","NDUFB8",
              "ZNF417","C6orf120")
genes_in <- intersect(GENES_16, rownames(expr))
log_info("16 基因命中: {length(genes_in)}")

# hUSI
husi_file <- file.path(DIR_RESULT, "04_1_aging_score", "data",
                       "aging_scores.rds")
aging_scores <- if (file.exists(husi_file)) readRDS(husi_file) else NULL

# ============================================================
# ============================================================
log_info("--- 2. limma 单对比 ---")

run_limma_16 <- function(expr_sub, clin, use_batch = FALSE) {
  covar <- data.frame(
    age_s1 = splines::ns(clin$age_v, df = 3)[, 1],
    age_s2 = splines::ns(clin$age_v, df = 3)[, 2],
    age_s3 = splines::ns(clin$age_v, df = 3)[, 3],
    sex_M = as.integer(clin$sex == "Male"),
    rin = clin$rin_v, pmi = clin$pmi_v,
    stringsAsFactors = FALSE)
  if (use_batch) covar$batch_f <- clin$batch_f
  rownames(covar) <- rownames(clin)
  design_df <- data.frame(grade = factor(clin$grade), covar)
  rownames(design_df) <- rownames(clin)
  design_df <- design_df[complete.cases(design_df), , drop = FALSE]
  expr_use <- expr_sub[, rownames(design_df), drop = FALSE]

  fml <- as.formula(paste("~ 0 + grade +",
                          paste(colnames(covar), collapse = " + ")))
  design <- model.matrix(fml, data = design_df)
  colnames(design)[grep("^grade", colnames(design))] <- paste0("g", 1:4)

  cm <- limma::makeContrasts(
    g2_vs_g1 = g2 - g1,
    g3_vs_g2 = g3 - g2,
    g4_vs_g3 = g4 - g3,
    g4_vs_g1 = g4 - g1,
    levels = design)
  fit <- eBayes(contrasts.fit(lmFit(expr_use, design), cm))

  out <- data.frame()
  for (ct in colnames(cm)) {
    tt <- topTable(fit, coef = ct, number = Inf, sort.by = "none")
    out <- rbind(out, data.frame(
      gene = rownames(tt), contrast = ct,
      logFC = tt$logFC, p = tt$P.Value, FDR = tt$adj.P.Val,
      stringsAsFactors = FALSE))
  }
  out
}

expr_16 <- expr[genes_in, , drop = FALSE]
common_16 <- intersect(colnames(expr_16), rownames(clin))
expr_16 <- expr_16[, common_16, drop = FALSE]
clin_16 <- clin[common_16, ]

limma_nobatch <- run_limma_16(expr_16, clin_16, use_batch = FALSE)
limma_batch   <- run_limma_16(expr_16, clin_16, use_batch = TRUE)

limma_16 <- merge(limma_nobatch, limma_batch,
                   by = c("gene","contrast"),
                   suffixes = c("_nobatch","_batch"))
limma_16 <- limma_16[order(limma_16$gene, limma_16$contrast), ]
save_table(limma_16, DIR_THIS, "genes_limma_per_contrast.csv")

cat("\n=== 16 基因 limma 各对比（不控 batch）===\n")
print(limma_16[, c("gene","contrast","logFC_nobatch",
                    "p_nobatch","FDR_nobatch")],
      row.names = FALSE, digits = 4)

cat("\n=== 各对比 FDR<0.05 基因数（不控/控 batch）===\n")
for (ct in c("g2_vs_g1","g3_vs_g2","g4_vs_g3","g4_vs_g1")) {
  sub <- limma_16[limma_16$contrast == ct, ]
  n1 <- sum(sub$FDR_nobatch < 0.05, na.rm = TRUE)
  n2 <- sum(sub$FDR_batch < 0.05, na.rm = TRUE)
  cat(sprintf("  %-12s  nobatch=%d  batch=%d\n", ct, n1, n2))
}

# ============================================================
# ============================================================
log_info("--- 3. poly 项系数 ---")

poly_mat <- poly(clin_16$grade, 3)
colnames(poly_mat) <- c("grade_p1","grade_p2","grade_p3")
rownames(poly_mat) <- rownames(clin_16)

covar <- data.frame(
  poly_mat,
  age_s1 = splines::ns(clin_16$age_v, df = 3)[, 1],
  age_s2 = splines::ns(clin_16$age_v, df = 3)[, 2],
  age_s3 = splines::ns(clin_16$age_v, df = 3)[, 3],
  sex_M = as.integer(clin_16$sex == "Male"),
  rin = clin_16$rin_v, pmi = clin_16$pmi_v,
  batch_f = clin_16$batch_f,
  stringsAsFactors = FALSE)
rownames(covar) <- rownames(clin_16)
covar <- covar[complete.cases(covar), , drop = FALSE]
expr_use <- expr_16[, rownames(covar), drop = FALSE]

design <- model.matrix(~ grade_p1 + grade_p2 + grade_p3 +
                        age_s1 + age_s2 + age_s3 +
                        sex_M + rin + pmi + batch_f,
                        data = covar)
fit <- eBayes(lmFit(expr_use, design))

coef_poly <- c("grade_p1","grade_p2","grade_p3")
tt_all <- topTable(fit, coef = coef_poly, number = Inf, sort.by = "none")
poly_terms <- data.frame(
  gene = rownames(tt_all),
  p1_logFC = tt_all$grade_p1,
  p2_logFC = tt_all$grade_p2,
  p3_logFC = tt_all$grade_p3,
  stringsAsFactors = FALSE)

poly_p <- data.frame(gene = rownames(expr_use))
for (cf in coef_poly) {
  tt <- topTable(fit, coef = cf, number = Inf, sort.by = "none")
  poly_p[[paste0(cf, "_p")]] <- tt$P.Value
  poly_p[[paste0(cf, "_FDR")]] <- tt$adj.P.Val
}
poly_terms <- merge(poly_terms, poly_p, by = "gene")
save_table(poly_terms, DIR_THIS, "genes_poly_terms.csv")

cat("\n=== 16 基因 poly 项系数（控 batch）===\n")
print(poly_terms, row.names = FALSE, digits = 4)

# ============================================================
# ============================================================
log_info("--- 4. grade 均值 ---")

expr_z <- t(scale(t(expr_16)))
expr_z[is.na(expr_z)] <- 0

grade_means <- data.frame()
for (g in genes_in) {
  for (gr in 1:4) {
    idx <- clin_16$grade == gr
    vals <- expr_z[g, idx]
    grade_means <- rbind(grade_means, data.frame(
      gene = g, grade = gr,
      mean_z = mean(vals, na.rm = TRUE),
      se = sd(vals, na.rm = TRUE) / sqrt(sum(idx)),
      stringsAsFactors = FALSE))
  }
}
save_table(grade_means, DIR_THIS, "genes_grade_means.csv")

# ============================================================
# ============================================================
log_info("--- 5. 模式分类 ---")

pattern_class <- data.frame()
for (g in genes_in) {
  means <- grade_means$mean_z[grade_means$gene == g]
  mono_up <- all(diff(means) > 0)
  mono_down <- all(diff(means) < 0)
  zigzag_1234 <- means[2] > means[1] && means[2] > means[3] && means[4] > means[3]
  zigzag_3214 <- means[1] > means[2] && means[3] > means[2] && means[3] > means[4]

  pattern <- "其他"
  if (mono_up) pattern <- "单调上升"
  if (mono_down) pattern <- "单调下降"
  if (zigzag_1234) pattern <- "1<2>3<4 锯齿"
  if (zigzag_3214) pattern <- "1>2<3>4 倒锯齿"

  peak_grade <- which.max(means)
  trough_grade <- which.min(means)

  pattern_class <- rbind(pattern_class, data.frame(
    gene = g,
    pattern = pattern,
    peak_grade = peak_grade,
    trough_grade = trough_grade,
    range = max(means) - min(means),
    p1 = means[1], p2 = means[2], p3 = means[3], p4 = means[4],
    stringsAsFactors = FALSE))
}
save_table(pattern_class, DIR_THIS, "genes_pattern_class.csv")

cat("\n=== 16 基因模式分类 ===\n")
print(pattern_class[, c("gene","pattern","peak_grade","trough_grade","range")],
      row.names = FALSE, digits = 3)

cat("\n模式分布:\n")
print(table(pattern_class$pattern))

# ============================================================
# ============================================================
log_info("--- 6. 与 hUSI / batch 相关 ---")

cor_out <- data.frame()
if (!is.null(aging_scores) && "hUSI" %in% rownames(aging_scores)) {
  common_h <- intersect(colnames(aging_scores), rownames(clin_16))
  hUSI_v <- aging_scores["hUSI", common_h]
  for (g in genes_in) {
    r_husi <- cor(expr_16[g, common_h], hUSI_v,
                   method = "spearman", use = "complete.obs")
    r_batch <- cor(expr_16[g, ],
                    as.integer(clin_16$batch_f),
                    method = "spearman", use = "complete.obs")
    r_grade <- cor(expr_16[g, ], clin_16$grade,
                    method = "spearman", use = "complete.obs")
    cor_out <- rbind(cor_out, data.frame(
      gene = g,
      r_hUSI = round(r_husi, 3),
      r_batch = round(r_batch, 3),
      r_grade = round(r_grade, 3),
      stringsAsFactors = FALSE))
  }
}
save_table(cor_out, DIR_THIS, "genes_hUSI_correlation.csv")

cat("\n=== 16 基因与 hUSI / batch / grade 相关 ===\n")
print(cor_out, row.names = FALSE)

# ============================================================
# 7. figure
# ============================================================
log_info("--- 7. 图 ---")

grade_means$gene_f <- factor(grade_means$gene, levels = genes_in)
p1 <- ggplot(grade_means, aes(x = grade, y = mean_z,
                                group = gene)) +
  geom_line(linewidth = 0.8, color = "#4C72B0") +
  geom_point(size = 2, color = "#4C72B0") +
  geom_errorbar(aes(ymin = mean_z - se, ymax = mean_z + se),
                 width = 0.15, alpha = 0.5, color = "#4C72B0") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray70") +
  facet_wrap(~ gene_f, scales = "free_y", ncol = 4) +
  scale_x_continuous(breaks = 1:4, labels = paste0("G", 1:4)) +
  theme_bw(base_size = 10) +
  labs(x = "AMD grade", y = "z-score (mean ± SE)",
       title = "16 poly3_F genes: grade means")
save_fig(p1, DIR_THIS, "genes_grade_curves", width = 12, height = 10)

poly_hm <- as.matrix(poly_terms[, c("p1_logFC","p2_logFC","p3_logFC")])
rownames(poly_hm) <- poly_terms$gene
pheatmap(poly_hm,
         cluster_rows = TRUE, cluster_cols = FALSE,
         display_numbers = TRUE, number_format = "%.2f",
         color = colorRampPalette(c("#4C72B0","#FFFFFF","#C44E52"))(50),
         main = "Poly terms (p1/p2/p3 coefficients)",
         filename = file.path(DIR_THIS, "figures",
                              "genes_poly_terms_heatmap.png"),
         width = 5, height = 7)

pattern_count <- as.data.frame(table(pattern_class$pattern))
p3 <- ggplot(pattern_count, aes(x = reorder(Var1, Freq), y = Freq,
                                  fill = Var1)) +
  geom_col(alpha = 0.8) + coord_flip() + theme_bw(base_size = 11) +
  labs(x = "", y = "Count", fill = "Pattern",
       title = "16 genes: pattern classification") +
  theme(legend.position = "none")
save_fig(p3, DIR_THIS, "genes_pattern_breakdown", width = 6, height = 4)

# ============================================================
# 8. summary
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("=== diag_16genes_grade_pattern.r 汇总 ===\n")
cat(strrep("=", 70), "\n", sep = "")

cat("\n各对比 FDR<0.05 基因数:\n")
for (ct in c("g2_vs_g1","g3_vs_g2","g4_vs_g3","g4_vs_g1")) {
  sub <- limma_16[limma_16$contrast == ct, ]
  n1 <- sum(sub$FDR_nobatch < 0.05, na.rm = TRUE)
  n2 <- sum(sub$FDR_batch < 0.05, na.rm = TRUE)
  cat(sprintf("  %-12s  nobatch=%d  batch=%d\n", ct, n1, n2))
}

cat("\n模式分布:\n")
print(table(pattern_class$pattern))

cat("\npoly 项显著的基因数（控 batch）:\n")
cat("  p1 FDR<0.05:", sum(poly_terms$grade_p1_FDR < 0.05, na.rm=TRUE), "\n")
cat("  p2 FDR<0.05:", sum(poly_terms$grade_p2_FDR < 0.05, na.rm=TRUE), "\n")
cat("  p3 FDR<0.05:", sum(poly_terms$grade_p3_FDR < 0.05, na.rm=TRUE), "\n")

cat("\n--- 输出位置 ---\n")
cat("目录:", DIR_THIS, "\n")

log_info("=== diag_16genes_grade_pattern.r done ===")
