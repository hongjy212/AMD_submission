#!/usr/bin/env Rscript
# ============================================================
# 04_1b_aging_score_layer_batch.r (v2: 419)
#
#
# Input：
#   results/01_1b_batch_enhance/data/expr_wgcna.rds        ★ 419
#   results/01_1b_batch_enhance/data/clinical_wgcna.rds    ★ 419
#   results/01_1_prepare_bulk_data/data/gene_annotation.rds
#
#   tables/hUSI_limma_multi_adjusted.csv
#   tables/hUSI_poly3F_adjusted.csv
#   tables/hUSI_age_comparison_adjusted.csv
#   tables/hUSI_age_classification_adjusted.csv
#   tables/hUSI_batch_within.csv
#   tables/hUSI_pathway_correlation.csv
#   tables/hUSI_16genes_correlation.csv
#   tables/hUSI_compare_models.csv
#   figures/hUSI_grade_curves_adjusted.pdf/.png
#   figures/hUSI_batch_within.pdf/.png
#   figures/hUSI_pathway_correlation.png
# ============================================================

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(limma); library(splines); library(ggplot2)
  library(pheatmap); library(tidyr)
})

DIR_THIS <- setup_output_dir("04_1_aging_score")
log_info("=== 04_1b_aging_score_layer_batch.r v2 (419) start ===")

# ============================================================
# ============================================================
log_info("--- 1. 载入 ---")

clin <- read_upstream("01_1b_batch_enhance", "data", "clinical_wgcna.rds")
expr <- read_upstream("01_1b_batch_enhance", "data", "expr_wgcna.rds")
anno <- read_upstream("01_1_prepare_bulk_data", "data", "gene_annotation.rds")
aging_scores <- read_upstream("04_1_aging_score", "data", "aging_scores.rds")
score_mat <- read_upstream("04_2_pathway_score", "data",
                            "pathway_scores_ssgsea.rds")

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

sid <- as.character(clin$sample_id)
clin$layer <- sub(".*_(L[0-9]+)$", "\\1", sid)
clin$layer <- factor(clin$layer, levels = paste0("L", 1:8))

rownames(clin) <- colnames(expr)

common_samp <- Reduce(intersect, list(
  colnames(aging_scores),
  rownames(clin)))
aging_sub <- aging_scores[, common_samp, drop = FALSE]
clin_sub <- clin[common_samp, ]

log_info("aging scores: {nrow(aging_sub)} × {ncol(aging_sub)}")
log_info("评分: {paste(rownames(aging_sub), collapse=', ')}")
log_info("样本: {ncol(aging_sub)}")

# ============================================================
# ============================================================
make_covar <- function(clin, use_layer = TRUE, use_batch = TRUE) {
  sp <- splines::ns(clin$age_v, df = 3)
  covar <- data.frame(
    age_s1 = sp[, 1], age_s2 = sp[, 2], age_s3 = sp[, 3],
    sex_M = as.integer(clin$sex == "Male"),
    rin = clin$rin_v, pmi = clin$pmi_v,
    stringsAsFactors = FALSE)
  if (use_layer) covar$layer <- clin$layer
  if (use_batch) covar$batch <- clin$batch_f
  rownames(covar) <- rownames(clin)
  covar
}

# ============================================================
# ============================================================
log_info("--- 3. limma multi (3 种模型) ---")

run_limma <- function(score_mat, clin, use_layer, use_batch, model_name) {
  covar <- make_covar(clin, use_layer, use_batch)
  design_df <- data.frame(grade = factor(clin$grade), covar)
  rownames(design_df) <- rownames(clin)
  design_df <- design_df[complete.cases(design_df), , drop = FALSE]
  score_sub <- score_mat[, rownames(design_df), drop = FALSE]

  rhs <- c("grade", colnames(covar))
  fml <- as.formula(paste("~ 0 +", paste(rhs, collapse = " + ")))
  design <- model.matrix(fml, data = design_df)
  colnames(design)[grep("^grade", colnames(design))] <- paste0("g", 1:4)

  cm <- limma::makeContrasts(
    g2_vs_g1 = g2 - g1,
    g3_vs_g2 = g3 - g2,
    g4_vs_g3 = g4 - g3,
    g4_vs_g1 = g4 - g1,
    levels = design)
  fit <- eBayes(contrasts.fit(lmFit(score_sub, design), cm))

  out <- data.frame()
  for (ct in colnames(cm)) {
    tt <- topTable(fit, coef = ct, number = Inf, sort.by = "none")
    out <- rbind(out, data.frame(
      feature = rownames(tt), model = model_name, contrast = ct,
      logFC = tt$logFC, p = tt$P.Value, FDR = tt$adj.P.Val,
      stringsAsFactors = FALSE))
  }
  out
}

models_limma <- rbind(
  run_limma(aging_sub, clin_sub, FALSE, FALSE, "M1_age_PMI_RIN"),
  run_limma(aging_sub, clin_sub, FALSE, TRUE,  "M1_plus_batch"),
  run_limma(aging_sub, clin_sub, TRUE,  TRUE,  "M1_plus_layer_batch")
)
save_table(models_limma, DIR_THIS, "hUSI_limma_multi_adjusted.csv")

cat("\n=== 4 种衰老评分 limma (3 模型对比) ===\n")
print(models_limma[, c("feature","model","contrast",
                        "logFC","p","FDR")],
      row.names = FALSE, digits = 4)

cmp <- reshape(models_limma[, c("feature","model","contrast","logFC","p","FDR")],
                idvar = c("feature","contrast"),
                timevar = "model", direction = "wide")
save_table(cmp, DIR_THIS, "hUSI_compare_models.csv")

# ============================================================
# ============================================================
log_info("--- 4. poly3_F (3 种模型) ---")

run_poly3F <- function(score_mat, clin, use_layer, use_batch, model_name) {
  poly_mat <- poly(clin$grade, 3)
  colnames(poly_mat) <- c("grade_p1", "grade_p2", "grade_p3")
  rownames(poly_mat) <- rownames(clin)

  covar <- make_covar(clin, use_layer, use_batch)
  design_df <- data.frame(poly_mat, covar)
  rownames(design_df) <- rownames(clin)
  design_df <- design_df[complete.cases(design_df), , drop = FALSE]
  score_sub <- score_mat[, rownames(design_df), drop = FALSE]

  rhs <- c("grade_p1", "grade_p2", "grade_p3", colnames(covar))
  fml <- as.formula(paste("~", paste(rhs, collapse = " + ")))
  design <- model.matrix(fml, data = design_df)
  fit <- eBayes(lmFit(score_sub, design))

  coef_poly <- grep("^grade_p", colnames(design), value = TRUE)
  tt_F <- topTable(fit, coef = coef_poly, number = Inf, sort.by = "none")

  out <- data.frame(
    feature = rownames(tt_F), model = model_name,
    F_stat = tt_F$F, p = tt_F$P.Value, FDR = tt_F$adj.P.Val,
    stringsAsFactors = FALSE)
  out
}

poly_models <- rbind(
  run_poly3F(aging_sub, clin_sub, FALSE, FALSE, "M1_age_PMI_RIN"),
  run_poly3F(aging_sub, clin_sub, FALSE, TRUE,  "M1_plus_batch"),
  run_poly3F(aging_sub, clin_sub, TRUE,  TRUE,  "M1_plus_layer_batch")
)
save_table(poly_models, DIR_THIS, "hUSI_poly3F_adjusted.csv")

cat("\n=== 4 种衰老评分 poly3_F (3 模型) ===\n")
print(poly_models, row.names = FALSE, digits = 4)

# ============================================================
# 5. agediagnosis
# ============================================================
log_info("--- 5. 年龄诊断 ---")

run_age_mode <- function(score_mat, clin, age_mode) {
  covar <- data.frame(
    sex_M = as.integer(clin$sex == "Male"),
    rin = clin$rin_v, pmi = clin$pmi_v,
    layer = clin$layer, batch = clin$batch_f,
    stringsAsFactors = FALSE)
  if (age_mode == "linear") covar$age <- clin$age_v
  else if (age_mode == "spline") {
    sp <- splines::ns(clin$age_v, df = 3)
    colnames(sp) <- c("age_s1", "age_s2", "age_s3")
    covar <- cbind(covar, sp)
  }
  rownames(covar) <- rownames(clin)
  design_df <- data.frame(grade = factor(clin$grade), covar)
  rownames(design_df) <- rownames(clin)
  design_df <- design_df[complete.cases(design_df), , drop = FALSE]
  score_sub <- score_mat[, rownames(design_df), drop = FALSE]

  fml <- as.formula(paste("~ 0 + grade +",
                          paste(colnames(covar), collapse = " + ")))
  design <- model.matrix(fml, data = design_df)
  colnames(design)[grep("^grade", colnames(design))] <- paste0("g", 1:4)
  cm <- limma::makeContrasts(g3_vs_g2 = g3 - g2, g4_vs_g1 = g4 - g1,
                              levels = design)
  fit <- eBayes(contrasts.fit(lmFit(score_sub, design), cm))

  out <- data.frame()
  for (ct in colnames(cm)) {
    tt <- topTable(fit, coef = ct, number = Inf, sort.by = "none")
    out <- rbind(out, data.frame(
      feature = rownames(tt), age_mode = age_mode, contrast = ct,
      logFC = tt$logFC, p = tt$P.Value, FDR = tt$adj.P.Val,
      stringsAsFactors = FALSE))
  }
  out
}

age_cmp <- rbind(
  run_age_mode(aging_sub, clin_sub, "none"),
  run_age_mode(aging_sub, clin_sub, "linear"),
  run_age_mode(aging_sub, clin_sub, "spline"))
save_table(age_cmp, DIR_THIS, "hUSI_age_comparison_adjusted.csv")

cat("\n=== g3_vs_g2 年龄诊断 ===\n")
print(age_cmp[age_cmp$contrast == "g3_vs_g2", ],
      row.names = FALSE, digits = 4)

# ============================================================
# ============================================================
log_info("--- 6. 年龄双模型分类 ---")

age_class <- data.frame()
for (i in seq_len(nrow(aging_sub))) {
  sc_name <- rownames(aging_sub)[i]
  y <- aging_sub[i, ]

  df <- data.frame(y = y, grade = clin_sub$grade, age = clin_sub$age_v,
                    sex = clin_sub$sex, rin = clin_sub$rin_v,
                    pmi = clin_sub$pmi_v,
                    layer = clin_sub$layer, batch = clin_sub$batch_f)
  df <- df[complete.cases(df), ]

  sp_g <- splines::ns(df$age, df = 3)
  fit1 <- lm(y ~ poly(grade, 3) + sp_g + sex + rin + pmi + layer + batch,
              data = df)
  p_g <- anova(fit1)["poly(grade, 3)", "Pr(>F)"]
  fit2 <- lm(y ~ sp_g + poly(grade, 3) + sex + rin + pmi + layer + batch,
              data = df)
  p_a <- anova(fit2)["sp_g", "Pr(>F)"]

  age_class <- rbind(age_class, data.frame(
    feature = sc_name, p_grade = p_g, p_age = p_a,
    stringsAsFactors = FALSE))
}

age_class$type <- "无关"
age_class$type[age_class$p_grade < 0.05 & age_class$p_age > 0.10] <- "纯疾病"
age_class$type[age_class$p_grade < 0.05 & age_class$p_age < 0.05] <- "年龄+疾病"
age_class$type[age_class$p_grade > 0.10 & age_class$p_age < 0.05] <- "纯年龄"
save_table(age_class, DIR_THIS, "hUSI_age_classification_adjusted.csv")

cat("\n=== 年龄双模型分类 ===\n")
print(age_class, row.names = FALSE, digits = 4)

# ============================================================
# ============================================================
log_info("--- 7. batch 内验证 ---")

run_within <- function(score_mat, clin, batch_val) {
  idx <- clin$batch_f == batch_val
  sub_clin <- clin[idx, ]
  sub_score <- score_mat[, rownames(sub_clin), drop = FALSE]

  sp <- splines::ns(sub_clin$age_v, df = 3)
  covar <- data.frame(
    age_s1 = sp[, 1], age_s2 = sp[, 2], age_s3 = sp[, 3],
    sex_M = as.integer(sub_clin$sex == "Male"),
    rin = sub_clin$rin_v, pmi = sub_clin$pmi_v,
    layer = sub_clin$layer,
    stringsAsFactors = FALSE)
  rownames(covar) <- rownames(sub_clin)
  design_df <- data.frame(grade = factor(sub_clin$grade), covar)
  rownames(design_df) <- rownames(sub_clin)
  design_df <- design_df[complete.cases(design_df), , drop = FALSE]
  score_use <- sub_score[, rownames(design_df), drop = FALSE]

  design <- model.matrix(
    ~ 0 + grade + age_s1 + age_s2 + age_s3 + sex_M + rin + pmi + layer,
    data = design_df)
  colnames(design)[grep("^grade", colnames(design))] <- paste0("g", 1:4)
  cm <- limma::makeContrasts(g3_vs_g2 = g3 - g2, levels = design)
  fit <- eBayes(contrasts.fit(lmFit(score_use, design), cm))
  tt <- topTable(fit, coef = "g3_vs_g2", number = Inf, sort.by = "none")
  data.frame(feature = rownames(tt), batch = batch_val,
             logFC = tt$logFC, p = tt$P.Value, FDR = tt$adj.P.Val,
             stringsAsFactors = FALSE)
}

batch_within <- do.call(rbind, lapply(levels(clin_sub$batch_f), function(b) {
  run_within(aging_sub, clin_sub, b)
}))
save_table(batch_within, DIR_THIS, "hUSI_batch_within.csv")

cat("\n=== batch 内 g3_vs_g2 ===\n")
print(batch_within, row.names = FALSE, digits = 4)

# ============================================================
# ============================================================
log_info("--- 8. 与 04_2 通路相关 ---")

common_pw <- intersect(colnames(score_mat), common_samp)
hUSI_v <- aging_sub["hUSI", common_pw]
cor_tbl <- data.frame()
for (gs in rownames(score_mat)) {
  r <- cor(score_mat[gs, common_pw], hUSI_v, use = "complete.obs")
  cor_tbl <- rbind(cor_tbl, data.frame(
    pathway = gs, r_hUSI = round(r, 3),
    stringsAsFactors = FALSE))
}
gs_used_path <- file.path(DIR_RESULT, "04_2_pathway_score", "tables",
                          "gene_sets_used.csv")
if (file.exists(gs_used_path)) {
  gs_used <- read.csv(gs_used_path, stringsAsFactors = FALSE)
  cor_tbl$group <- gs_used$group[match(cor_tbl$pathway, gs_used$geneset)]
}
cor_tbl <- cor_tbl[order(-abs(cor_tbl$r_hUSI)), ]
save_table(cor_tbl, DIR_THIS, "hUSI_pathway_correlation.csv")

cat("\n=== hUSI 与通路相关 top 20 ===\n")
print(head(cor_tbl, 20), row.names = FALSE)

# ============================================================
# ============================================================
log_info("--- 9. 与 16 基因相关 ---")

GENES_16 <- c("FCF1","UGGT1","CCDC3","S100B","BPNT1","BRI3BP","BTBD19",
              "DPH7","AKR1B1","ZNF623","ZNF587","RSPH3","TAF3","NDUFB8",
              "ZNF417","C6orf120")
genes_in <- intersect(GENES_16, rownames(expr))
common_16 <- intersect(common_samp, colnames(expr))
hUSI_16 <- aging_sub["hUSI", common_16]

cor_16 <- data.frame()
for (g in genes_in) {
  r <- cor(expr[g, common_16], hUSI_16, use = "complete.obs")
  cor_16 <- rbind(cor_16, data.frame(
    gene = g, r_hUSI = round(r, 3),
    stringsAsFactors = FALSE))
}
cor_16 <- cor_16[order(-abs(cor_16$r_hUSI)), ]
save_table(cor_16, DIR_THIS, "hUSI_16genes_correlation.csv")

cat("\n=== 16 基因与 hUSI 相关 ===\n")
print(cor_16, row.names = FALSE)

# ============================================================
# 10. figure
# ============================================================
log_info("--- 10. 图 ---")

grade_means <- data.frame()
for (sc in rownames(aging_sub)) {
  for (g in 1:4) {
    idx <- clin_sub$grade == g
    vals <- aging_sub[sc, idx]
    grade_means <- rbind(grade_means, data.frame(
      score = sc, grade = g,
      mean_z = mean(vals, na.rm = TRUE),
      se = sd(vals, na.rm = TRUE) / sqrt(sum(idx)),
      stringsAsFactors = FALSE))
  }
}
p1 <- ggplot(grade_means, aes(x = grade, y = mean_z,
                                color = score, group = score)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean_z - se, ymax = mean_z + se),
                 width = 0.15, alpha = 0.5) +
  theme_bw(base_size = 11) +
  scale_x_continuous(breaks = 1:4, labels = paste0("G", 1:4)) +
  labs(x = "AMD grade", y = "Score (mean ± SE)",
       title = "4 aging scores by grade (layer+batch adjusted, 419)",
       color = "Score")
save_fig(p1, DIR_THIS, "hUSI_grade_curves_adjusted", width = 8, height = 5)

p2 <- ggplot(batch_within, aes(x = reorder(feature, logFC),
                                 y = logFC, fill = batch)) +
  geom_col(position = "dodge", alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  coord_flip() + theme_bw(base_size = 11) +
  labs(x = "", y = "logFC (g3_vs_g2)", fill = "Batch",
       title = "hUSI & aging scores: batch within-batch g3_vs_g2 (419)")
save_fig(p2, DIR_THIS, "hUSI_batch_within", width = 7, height = 4)

top_n <- 30
top_cor <- head(cor_tbl[order(-abs(cor_tbl$r_hUSI)), ], top_n)
hm_mat <- matrix(top_cor$r_hUSI, ncol = 1,
                  dimnames = list(top_cor$pathway, "r_hUSI"))
pheatmap(hm_mat,
         cluster_rows = TRUE, cluster_cols = FALSE,
         display_numbers = TRUE, number_format = "%.2f",
         color = colorRampPalette(c("#4C72B0","#FFFFFF","#C44E52"))(50),
         main = "hUSI correlation with pathways (top 30, 419)",
         filename = file.path(DIR_THIS, "figures",
                              "hUSI_pathway_correlation.png"),
         width = 5, height = 8)

# ============================================================
# 11. summary
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("=== 04_1b 汇总 (v2, 419) ===\n")
cat(strrep("=", 70), "\n", sep = "")

cat("\n【核心问题 1】hUSI 的 g3_vs_g2 效应\n")
hUSI_g3 <- models_limma[models_limma$feature == "hUSI" &
                          models_limma$contrast == "g3_vs_g2", ]
for (m in unique(hUSI_g3$model)) {
  r <- hUSI_g3[hUSI_g3$model == m, ]
  cat(sprintf("  %-22s  logFC=%+.4f  p=%.4f  FDR=%.4f\n",
              m, r$logFC, r$p, r$FDR))
}

cat("\n【核心问题 2】hUSI poly3_F\n")
hUSI_poly <- poly_models[poly_models$feature == "hUSI", ]
for (m in unique(hUSI_poly$model)) {
  r <- hUSI_poly[hUSI_poly$model == m, ]
  cat(sprintf("  %-22s  F=%.3f  p=%.4f  FDR=%.4f\n",
              m, r$F_stat, r$p, r$FDR))
}

cat("\n【核心问题 3】hUSI batch 内效应\n")
print(batch_within[batch_within$feature == "hUSI", ],
      row.names = FALSE, digits = 4)

cat("\n【核心问题 4】hUSI 与线粒体相关\n")
mito_pw <- c("Mito_CI_Subunits","Mito_CIII_Subunits","Mito_CIV_Subunits",
             "Mito_CV_Subunits","Mito_Complex_V","Mito_TCA","OXPHOS")
mito_cor <- cor_tbl[cor_tbl$pathway %in% mito_pw, ]
print(mito_cor, row.names = FALSE)

cat("\n--- 输出位置 ---\n")
cat("目录:", DIR_THIS, "\n")

log_info("=== 04_1b v2 (419) done ===")

