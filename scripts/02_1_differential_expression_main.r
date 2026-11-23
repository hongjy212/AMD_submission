# ============================================================
# 02_1_differential_expression_main.r (v4)
#
#
# Input：
#   results/01_1b_batch_enhance/data/expr_wgcna.rds        ★ 419
#   results/01_1b_batch_enhance/data/clinical_wgcna.rds    ★ 419
#   results/01_1_prepare_bulk_data/data/gene_annotation.rds
# Output：results/02_1_differential_expression_main/
# ============================================================

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(limma); library(splines)
})

DIR_THIS <- setup_output_dir("02_1_differential_expression_main")
log_info("=== 02_1_differential_expression_main.r v4 (419, sig16 from 453) start ===")

# ============================================================
# ============================================================
log_info("--- 1. 载入数据 ---")

expr <- read_upstream("01_1b_batch_enhance", "data", "expr_wgcna.rds")
clin <- read_upstream("01_1b_batch_enhance", "data", "clinical_wgcna.rds")
anno <- read_upstream("01_1_prepare_bulk_data", "data", "gene_annotation.rds")

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
log_info("grade: {paste(names(table(clin$grade)), table(clin$grade), sep='=', collapse=' ')}")

# ============================================================
# ============================================================
log_info("--- 2. 协变量 ---")

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

log_info("样本数: {nrow(design_df)}")

# ============================================================
# ============================================================
log_info("--- 3. limma multi (7 对比) ---")

design_multi <- model.matrix(
  ~ 0 + grade + age_s1 + age_s2 + age_s3 + sex_M + rin + pmi,
  data = design_df)
grade_idx <- grep("^grade", colnames(design_multi))
colnames(design_multi)[grade_idx] <- paste0("g", 1:4)

cm <- limma::makeContrasts(
  g2_vs_g1      = g2 - g1,
  g3_vs_g2      = g3 - g2,
  g4_vs_g3      = g4 - g3,
  g4_vs_g1      = g4 - g1,
  anyAMD_vs_g1  = (g2 + g3 + g4)/3 - g1,
  late_vs_early = (g3 + g4)/2 - (g1 + g2)/2,
  g4_vs_rest    = g4 - (g1 + g2 + g3)/3,
  levels = design_multi
)
contrast_names <- colnames(cm)
log_info("7 对比: {paste(contrast_names, collapse=', ')}")

fit_multi <- eBayes(contrasts.fit(lmFit(expr_sub, design_multi), cm))

multi_res <- data.frame()
for (ct in contrast_names) {
  tt <- topTable(fit_multi, coef = ct, number = Inf, sort.by = "none")
  multi_res <- rbind(multi_res, data.frame(
    gene = rownames(tt), contrast = ct,
    logFC = tt$logFC, t = tt$t, p = tt$P.Value, FDR = tt$adj.P.Val,
    stringsAsFactors = FALSE))
}
save_table(multi_res, DIR_THIS, "limma_multi_full.csv")
log_info("limma multi: {nrow(multi_res)} 行")

cat("\n=== limma multi 各对比 FDR<0.05 ===\n")
for (ct in contrast_names) {
  n_sig <- sum(multi_res$FDR[multi_res$contrast == ct] < 0.05, na.rm = TRUE)
  cat(sprintf("  %-15s  %d\n", ct, n_sig))
}

# ============================================================
# ============================================================
log_info("--- 4. poly3_F (419) ---")

poly_mat <- poly(clin$grade, 3)
colnames(poly_mat) <- c("grade_p1", "grade_p2", "grade_p3")
rownames(poly_mat) <- rownames(clin)
design_poly_df <- data.frame(poly_mat, covar_df)
rownames(design_poly_df) <- rownames(clin)
design_poly_df <- design_poly_df[complete.cases(design_poly_df), , drop = FALSE]
expr_poly <- expr[, rownames(design_poly_df), drop = FALSE]

design_poly <- model.matrix(
  ~ grade_p1 + grade_p2 + grade_p3 + age_s1 + age_s2 + age_s3 + sex_M + rin + pmi,
  data = design_poly_df)
fit_poly <- eBayes(lmFit(expr_poly, design_poly))

coef_poly <- grep("^grade_p", colnames(design_poly), value = TRUE)
tt_F <- topTable(fit_poly, coef = coef_poly, number = Inf, sort.by = "none")
polyF <- data.frame(
  gene = rownames(tt_F),
  t = tt_F$F,
  p = tt_F$P.Value,
  FDR = tt_F$adj.P.Val,
  stringsAsFactors = FALSE)
save_table(polyF, DIR_THIS, "poly3_F_full.csv")

log_info("poly3_F (419): FDR<0.05 = {sum(polyF$FDR < 0.05, na.rm=TRUE)}")
log_info("poly3_F (419): p<0.001   = {sum(polyF$p < 0.001, na.rm=TRUE)}")
log_info("poly3_F (419): p<0.01    = {sum(polyF$p < 0.01, na.rm=TRUE)}")

# ============================================================
# ============================================================
log_info("--- 5. 16 基因（从 453 继承） ---")

GENES_16_FALLBACK <- c(
  "AKR1B1", "FCF1", "RSPH3", "UGGT1", "DPH7", "CCDC3", "S100B", "BPNT1",
  "TAF3", "NDUFB8", "ZNF417", "ZNF623", "BRI3BP", "C6orf120", "ZNF587", "BTBD19"
)

sig16_file <- file.path(DIR_THIS, "tables", "sig16_genes.csv")
if (file.exists(sig16_file)) {
  sig16_prev <- read.csv(sig16_file, stringsAsFactors = FALSE)$gene
  if (length(sig16_prev) >= 16) {
    sig16 <- sig16_prev
    log_info("sig16 loaded from sig16_genes.csv: {length(sig16)} genes")
  } else {
    log_warn("sig16_genes.csv has {length(sig16_prev)} genes; using fallback 16-gene list")
    sig16 <- GENES_16_FALLBACK
  }
} else {
  log_warn("sig16_genes.csv not found; using fallback 16-gene list")
  sig16 <- GENES_16_FALLBACK
}

save_table(data.frame(gene = sig16), DIR_THIS, "sig16_genes.csv")

sig16_detail <- polyF[polyF$gene %in% sig16, ]
sig16_detail <- sig16_detail[order(sig16_detail$FDR), ]
save_table(sig16_detail, DIR_THIS, "sig16_poly3_detail.csv")

log_info("sig16 详情（419 上的 poly3_F）: {nrow(sig16_detail)} 基因")
cat("\n=== 16 基因在 419 上的 poly3_F ===\n")
print(sig16_detail, row.names = FALSE, digits = 4)

sig16_detail$FDR_internal <- p.adjust(sig16_detail$p, method = "BH")
cat("\n候选集内部 FDR<0.05:", sum(sig16_detail$FDR_internal < 0.05, na.rm = TRUE),
    "/", nrow(sig16_detail), "\n")

save_table(sig16_detail, DIR_THIS, "sig16_poly3_detail_internalFDR.csv")

# ============================================================
# ============================================================
log_info("--- 5b. 16 基因年龄双模型分类 ---")

age_class_tbl <- data.frame()
for (g in sig16) {
  if (!g %in% rownames(expr)) {
    log_warn("gene not in expr: {g}")
    next
  }
  df <- data.frame(y = expr[g, ], grade = clin$grade,
                   age = clin$age_n, sex = clin$sex,
                   rin = clin$rin_n, pmi = clin$pmi_n)
  df <- df[complete.cases(df), ]

  sp_g <- splines::ns(df$age, df = 3)

  fit1 <- lm(y ~ poly(grade, 3) + sp_g + sex + rin + pmi, data = df)
  p_g <- anova(fit1)["poly(grade, 3)", "Pr(>F)"]

  fit2 <- lm(y ~ sp_g + poly(grade, 3) + sex + rin + pmi, data = df)
  p_a <- anova(fit2)["sp_g", "Pr(>F)"]

  r_age <- cor(df$y, df$age, method = "pearson")
  r_age_s <- cor(df$y, df$age, method = "spearman")

  resid_y <- residuals(lm(y ~ factor(grade) + sex + rin + pmi, data = df))
  resid_age <- residuals(lm(age ~ factor(grade) + sex + rin + pmi, data = df))
  r_age_partial <- cor(resid_y, resid_age)

  age_class_tbl <- rbind(age_class_tbl, data.frame(
    gene = g,
    p_grade = p_g,
    p_age = p_a,
    r_age = round(r_age, 3),
    r_age_spearman = round(r_age_s, 3),
    r_age_partial_grade = round(r_age_partial, 3),
    stringsAsFactors = FALSE))
}

if (nrow(age_class_tbl) > 0) {
  age_class_tbl$type <- "无关"
  age_class_tbl$type[age_class_tbl$p_grade < 0.05 &
                     age_class_tbl$p_age > 0.10] <- "纯疾病"
  age_class_tbl$type[age_class_tbl$p_grade < 0.05 &
                     age_class_tbl$p_age < 0.05] <- "年龄+疾病"
  age_class_tbl$type[age_class_tbl$p_grade > 0.10 &
                     age_class_tbl$p_age < 0.05] <- "纯年龄"

  age_class_tbl <- age_class_tbl[order(age_class_tbl$type,
                                        -abs(age_class_tbl$r_age)), ]
  save_table(age_class_tbl, DIR_THIS, "sig16_age_classification.csv")

  log_info("16 基因年龄分类：")
  print(table(age_class_tbl$type))
  cat("\n=== 16 基因年龄分类 ===\n")
  print(age_class_tbl, row.names = FALSE)
}

# ============================================================
# 6. figure
# ============================================================
log_info("--- 6. 图 ---")

if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)

  if (nrow(sig16_detail) > 0) {
    sig16_detail$gene <- factor(
      sig16_detail$gene,
      levels = rev(sig16_detail$gene[order(sig16_detail$t)]))
    p1 <- ggplot(sig16_detail, aes(x = gene, y = t, fill = FDR < 0.05)) +
      geom_col() +
      coord_flip() +
      theme_bw(base_size = 11) +
      scale_fill_manual(values = c("TRUE" = "#C44E52", "FALSE" = "#4C72B0"),
                         guide = "none") +
      labs(x = "", y = "poly3_F t statistic",
           title = "16 gene poly3_F (age=spline, 419)",
           subtitle = sprintf("FDR<0.05: %d genes",
                              sum(sig16_detail$FDR < 0.05, na.rm = TRUE)))
    save_fig(p1, DIR_THIS, "sig16_poly3_t", width = 6, height = 5)
  }

  if (nrow(age_class_tbl) > 0) {
    plot_df <- age_class_tbl
    plot_df$neglogp_grade <- -log10(pmax(plot_df$p_grade, 1e-300))

    if (requireNamespace("ggrepel", quietly = TRUE)) {
      library(ggrepel)
      p2 <- ggplot(plot_df,
                   aes(x = r_age, y = neglogp_grade,
                       color = type, label = gene)) +
        geom_point(size = 3) +
        geom_text_repel(size = 3, max.overlaps = 20) +
        geom_vline(xintercept = 0, color = "gray70") +
        geom_hline(yintercept = -log10(0.05),
                   linetype = "dashed", color = "gray50") +
        theme_bw(base_size = 12) +
        scale_color_manual(values = c(
          "纯疾病" = "#4C72B0",
          "年龄+疾病" = "#C44E52",
          "无关" = "#B0B0B0",
          "纯年龄" = "#FDB462")) +
        labs(x = "Pearson r with age",
             y = "-log10(p_grade)",
             color = "类型",
             title = "16 基因：年龄相关 vs 疾病显著性（419）")
      save_fig(p2, DIR_THIS, "sig16_age_vs_grade", width = 7, height = 5)
    }
  }
}

# ============================================================
# 7. summary
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("=== 02_1 汇总 (v4, 419) ===\n")
cat(strrep("=", 70), "\n", sep = "")

cat("\n--- 1. 数据 ---\n")
cat("表达矩阵:", nrow(expr), "×", ncol(expr), "\n")
cat("样本数:", nrow(design_df), "\n")
cat("grade 分布:", paste(names(table(clin$grade)),
                         table(clin$grade), sep = "=", collapse = " "), "\n")

cat("\n--- 2. limma multi (7 对比) ---\n")
for (ct in contrast_names) {
  n_sig <- sum(multi_res$FDR[multi_res$contrast == ct] < 0.05, na.rm = TRUE)
  cat(sprintf("  %-15s FDR<0.05: %d\n", ct, n_sig))
}

cat("\n--- 3. poly3_F (419) ---\n")
cat("FDR<0.05 (全转录组):", sum(polyF$FDR < 0.05, na.rm = TRUE), "\n")
cat("p<0.001            :", sum(polyF$p < 0.001, na.rm = TRUE), "\n")

cat("\n--- 4. 16 基因（从 453 继承）---\n")
cat("sig16 长度:", length(sig16), "\n")
cat("在 419 上的 poly3_F 候选集内部 FDR<0.05:",
    sum(sig16_detail$FDR_internal < 0.05, na.rm = TRUE), "/",
    nrow(sig16_detail), "\n")
cat("\n")
print(sig16_detail[, c("gene", "t", "p", "FDR", "FDR_internal")],
      row.names = FALSE, digits = 4)

if (nrow(age_class_tbl) > 0) {
  cat("\n--- 5. 16 基因年龄分类 ---\n")
  print(table(age_class_tbl$type))
}

log_info("=== 02_1_differential_expression_main.r v4 done ===")

