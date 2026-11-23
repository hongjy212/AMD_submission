#!/usr/bin/env Rscript
# ============================================================
# ============================================================

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(limma); library(splines); library(ggplot2)
  library(pheatmap); library(tidyr)
})

DIR_THIS <- setup_output_dir("diag_rin_sensitivity", subdir = "diag")
log_info("=== diag_rin_sensitivity.r (v2) start ===")

# ============================================================
# ============================================================
log_info("--- 1. 载入 ---")

expr <- read_upstream("01_1_prepare_bulk_data", "data", "expr_corrected.rds")
clin <- read_upstream("01_1_prepare_bulk_data", "data", "clinical.rds")
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
clin$layer <- factor(sub(".*_(L[0-9]+)$", "\\1", sid),
                      levels = paste0("L", 1:8))
rownames(clin) <- colnames(expr)

common_samp <- Reduce(intersect, list(
  colnames(aging_scores),
  colnames(score_mat),
  rownames(clin)))
aging_sub <- aging_scores[, common_samp, drop = FALSE]
score_sub <- score_mat[, common_samp, drop = FALSE]
expr_sub  <- expr[, common_samp, drop = FALSE]
clin_sub  <- clin[common_samp, ]

log_info("样本: {length(common_samp)}")
log_info("RIN 范围: {round(min(clin_sub$rin_v),1)} - {round(max(clin_sub$rin_v),1)}")
log_info("RIN 中位数: {round(median(clin_sub$rin_v), 1)}")

# ============================================================
# ============================================================
log_info("--- 2. RIN 在 grade 间分布 ---")

rin_dist <- do.call(rbind, lapply(sort(unique(clin_sub$grade)), function(g) {
  idx <- clin_sub$grade == g
  x <- clin_sub$rin_v[idx]
  data.frame(
    grade = g, n = sum(idx),
    mean = round(mean(x, na.rm = TRUE), 2),
    sd = round(sd(x, na.rm = TRUE), 2),
    median = round(median(x, na.rm = TRUE), 2),
    q25 = round(quantile(x, 0.25, na.rm = TRUE), 2),
    q75 = round(quantile(x, 0.75, na.rm = TRUE), 2),
    stringsAsFactors = FALSE)
}))
save_table(rin_dist, DIR_THIS, "rin_grade_distribution.csv")

cat("\n=== RIN 在 grade 间分布 ===\n")
print(rin_dist, row.names = FALSE)

kw_rin <- kruskal.test(clin_sub$rin_v ~ clin_sub$grade)
rho_rin_grade <- cor(clin_sub$rin_v, clin_sub$grade,
                      method = "spearman", use = "complete.obs")
cat(sprintf("\nKW p = %.4g\n", kw_rin$p.value))
cat(sprintf("Spearman rho(RIN, grade) = %.3f\n", rho_rin_grade))

# ============================================================
# ============================================================
log_info("--- 3. RIN 与协变量相关 ---")

covars <- list(
  grade = clin_sub$grade,
  age = clin_sub$age_v,
  pmi = clin_sub$pmi_v,
  batch = as.integer(clin_sub$batch_f),
  layer = as.integer(clin_sub$layer)
)
cross_out <- data.frame()
for (cv_name in names(covars)) {
  cv <- covars[[cv_name]]
  r <- cor(clin_sub$rin_v, cv, method = "spearman", use = "complete.obs")
  n <- sum(!is.na(clin_sub$rin_v) & !is.na(cv))
  t_stat <- r * sqrt((n - 2) / (1 - r^2))
  p_two <- 2 * pt(-abs(t_stat), df = n - 2)
  cross_out <- rbind(cross_out, data.frame(
    variable = cv_name,
    spearman_r = round(r, 3),
    n = n,
    p_two_sided = signif(p_two, 3),
    stringsAsFactors = FALSE))
}
if ("sex" %in% colnames(clin_sub)) {
  wt <- wilcox.test(clin_sub$rin_v ~ clin_sub$sex)
  cross_out <- rbind(cross_out, data.frame(
    variable = "sex(Wilcox)",
    spearman_r = NA,
    n = nrow(clin_sub),
    p_two_sided = signif(wt$p.value, 3),
    stringsAsFactors = FALSE))
}
save_table(cross_out, DIR_THIS, "rin_covariate_correlation.csv")

cat("\n=== RIN 与协变量相关 ===\n")
print(cross_out, row.names = FALSE)

# ============================================================
# ============================================================
log_info("--- 4. RIN 与通路评分相关 ---")

rin_cor <- data.frame()
for (gs in rownames(score_sub)) {
  r <- cor(score_sub[gs, ], clin_sub$rin_v,
           method = "spearman", use = "complete.obs")
  rin_cor <- rbind(rin_cor, data.frame(
    pathway = gs, r_RIN = round(r, 3), stringsAsFactors = FALSE))
}
gs_used_path <- file.path(DIR_RESULT, "04_2_pathway_score", "tables",
                          "gene_sets_used.csv")
if (file.exists(gs_used_path)) {
  gs_used <- read.csv(gs_used_path, stringsAsFactors = FALSE)
  rin_cor$group <- gs_used$group[match(rin_cor$pathway, gs_used$geneset)]
}
rin_cor <- rin_cor[order(-abs(rin_cor$r_RIN)), ]
save_table(rin_cor, DIR_THIS, "rin_pathway_correlation.csv")

cat("\n=== RIN 与通路相关 top 20 ===\n")
print(head(rin_cor[, c("pathway","group","r_RIN")], 20),
      row.names = FALSE)

cat("\n=== |r_RIN| 按 group 中位数 ===\n")
grp_med <- aggregate(abs(r_RIN) ~ group, data = rin_cor, FUN = median)
grp_med <- grp_med[order(-grp_med$`abs(r_RIN)`), ]
print(grp_med, row.names = FALSE)

r_hUSI_RIN <- cor(aging_sub["hUSI", ], clin_sub$rin_v,
                   method = "spearman", use = "complete.obs")
cat(sprintf("\nhUSI vs RIN: Spearman r = %.3f\n", r_hUSI_RIN))

# ============================================================
# ============================================================
make_covar <- function(clin, rin_mode = "linear", use_layer = TRUE,
                        use_batch = TRUE) {
  sp <- splines::ns(clin$age_v, df = 3)
  covar <- data.frame(
    age_s1 = sp[, 1], age_s2 = sp[, 2], age_s3 = sp[, 3],
    sex_M = as.integer(clin$sex == "Male"),
    pmi = clin$pmi_v,
    stringsAsFactors = FALSE)
  if (rin_mode == "linear") {
    covar$rin <- clin$rin_v
  } else if (rin_mode == "spline") {
    rs <- splines::ns(clin$rin_v, df = 3)
    colnames(rs) <- c("rin_s1", "rin_s2", "rin_s3")
    covar <- cbind(covar, rs)
  }
  if (use_layer) covar$layer <- clin$layer
  if (use_batch) covar$batch <- clin$batch_f
  rownames(covar) <- rownames(clin)
  covar
}

run_limma <- function(score_mat, clin, rin_mode, model_name,
                       use_layer = TRUE, use_batch = TRUE) {
  covar <- make_covar(clin, rin_mode, use_layer, use_batch)
  design_df <- data.frame(grade = factor(clin$grade), covar)
  rownames(design_df) <- rownames(clin)
  design_df <- design_df[complete.cases(design_df), , drop = FALSE]
  score_use <- score_mat[, rownames(design_df), drop = FALSE]

  rhs <- c("grade", colnames(covar))
  fml <- as.formula(paste("~ 0 +", paste(rhs, collapse = " + ")))
  design <- model.matrix(fml, data = design_df)
  colnames(design)[grep("^grade", colnames(design))] <- paste0("g", 1:4)
  cm <- limma::makeContrasts(
    g3_vs_g2 = g3 - g2,
    g4_vs_g3 = g4 - g3,
    levels = design)
  fit <- eBayes(contrasts.fit(lmFit(score_use, design), cm))

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

# ============================================================
# ============================================================
log_info("--- 6. 控 RIN 前后效应变化 ---")

rin_modes <- c("none", "linear", "spline")

all_effects <- data.frame()
for (rm in rin_modes) {
  model_name <- paste0("RIN_", rm)
  e1 <- run_limma(aging_sub, clin_sub, rm, model_name)
  e2 <- run_limma(score_sub, clin_sub, rm, model_name)
  all_effects <- rbind(all_effects, e1, e2)
}

all_effects$type <- ifelse(all_effects$feature %in% rownames(aging_sub),
                            "aging_score", "pathway")
save_table(all_effects, DIR_THIS, "rin_adjusted_effects.csv")

cat("\n=== 各模型下 g3_vs_g2 FDR<0.05 数量 ===\n")
for (m in unique(all_effects$model)) {
  n <- sum(all_effects$FDR[all_effects$model == m &
                            all_effects$contrast == "g3_vs_g2"] < 0.05,
           na.rm = TRUE)
  cat(sprintf("  %-15s  %d\n", m, n))
}

cat("\n=== hUSI 在 3 种 RIN 处理下 ===\n")
hUSI_eff <- all_effects[all_effects$feature == "hUSI" &
                          all_effects$contrast == "g3_vs_g2", ]
print(hUSI_eff[, c("model","logFC","p","FDR")],
      row.names = FALSE, digits = 4)

cor_all <- data.frame()
for (g in rownames(expr_sub)) {
  r <- cor(expr_sub[g, ], aging_sub["hUSI", ], method = "spearman",
           use = "complete.obs")
  cor_all <- rbind(cor_all, data.frame(
    gene = g, r_hUSI = r, stringsAsFactors = FALSE))
}
cor_all$abs_r <- abs(cor_all$r_hUSI)
cor_all <- cor_all[order(-cor_all$abs_r), ]
axis_genes <- head(cor_all, 100)$gene

cat("\n=== hUSI 轴 top 100 在 3 种 RIN 处理下 ===\n")
axis_effects <- data.frame()
for (rm in rin_modes) {
  model_name <- paste0("RIN_", rm)
  e <- run_limma(expr_sub[axis_genes, , drop = FALSE],
                  clin_sub, rm, model_name)
  e <- e[e$contrast == "g3_vs_g2", ]
  e$type <- "hUSI_axis_gene"
  axis_effects <- rbind(axis_effects, e)
}
axis_effects$direction_neg <- axis_effects$logFC < 0
save_table(axis_effects, DIR_THIS, "rin_hUSI_axis_correlation.csv")

cat("hUSI 轴 top 100 中 g3_vs_g2 FDR<0.05 数量:\n")
for (m in unique(axis_effects$model)) {
  sub <- axis_effects[axis_effects$model == m, ]
  n <- sum(sub$FDR < 0.05, na.rm = TRUE)
  cat(sprintf("  %-15s  %d / 100\n", m, n))
}

# ============================================================
# ============================================================
log_info("--- 7. RIN 分层内部分析 ---")

rin_q <- cut(clin_sub$rin_v,
              breaks = quantile(clin_sub$rin_v, probs = 0:3/3),
              labels = c("RIN_low", "RIN_mid", "RIN_high"),
              include.lowest = TRUE)
clin_sub$rin_stratum <- rin_q

cat("\n=== RIN 分层 × grade 样本数 ===\n")
print(table(clin_sub$rin_stratum, clin_sub$grade))

run_stratified <- function(score_mat, clin, stratum, model_name) {
  idx <- clin$rin_stratum == stratum & !is.na(clin$rin_stratum)
  sub_clin <- clin[idx, ]
  sub_score <- score_mat[, rownames(sub_clin), drop = FALSE]

  sp <- splines::ns(sub_clin$age_v, df = 3)
  covar <- data.frame(
    age_s1 = sp[, 1], age_s2 = sp[, 2], age_s3 = sp[, 3],
    sex_M = as.integer(sub_clin$sex == "Male"),
    pmi = sub_clin$pmi_v,
    layer = sub_clin$layer,
    batch = sub_clin$batch_f,
    stringsAsFactors = FALSE)
  rownames(covar) <- rownames(sub_clin)
  design_df <- data.frame(grade = factor(sub_clin$grade), covar)
  rownames(design_df) <- rownames(sub_clin)
  design_df <- design_df[complete.cases(design_df), , drop = FALSE]
  if (length(unique(design_df$grade)) < 3) return(NULL)
  score_use <- sub_score[, rownames(design_df), drop = FALSE]

  design <- model.matrix(
    ~ 0 + grade + age_s1 + age_s2 + age_s3 + sex_M + pmi + layer + batch,
    data = design_df)
  colnames(design)[grep("^grade", colnames(design))] <- paste0("g", 1:4)
  cm <- limma::makeContrasts(g3_vs_g2 = g3 - g2, levels = design)
  fit <- eBayes(contrasts.fit(lmFit(score_use, design), cm))
  tt <- topTable(fit, coef = "g3_vs_g2", number = Inf, sort.by = "none")
  data.frame(
    feature = rownames(tt), stratum = stratum, model = model_name,
    logFC = tt$logFC, p = tt$P.Value, FDR = tt$adj.P.Val,
    n = nrow(design_df),
    stringsAsFactors = FALSE)
}

aging_features <- rownames(aging_sub)
pathway_features <- intersect(
  c("Mito_CI_Subunits","Mito_CIII_Subunits","Mito_CIV_Subunits",
    "Mito_CV_Subunits","Mito_TCA","OXPHOS","Mito_Complex_V",
    "Mito_Complex_I","Mito_Complex_III","Mito_Complex_IV"),
  rownames(score_sub))
score_focus <- rbind(
  aging_sub[aging_features, , drop = FALSE],
  score_sub[pathway_features, , drop = FALSE])
score_focus <- score_focus[!duplicated(rownames(score_focus)), ]
log_info("focus 特征: {nrow(score_focus)}")
log_info("  aging: {paste(aging_features, collapse=', ')}")
log_info("  pathway: {paste(pathway_features, collapse=', ')}")

stratified_out <- data.frame()
for (st in c("RIN_low", "RIN_mid", "RIN_high")) {
  r <- tryCatch(run_stratified(score_focus, clin_sub, st, st),
                 error = function(e) {
                   log_warn("{st}: {e$message}")
                   NULL
                 })
  if (!is.null(r)) stratified_out <- rbind(stratified_out, r)
}

if (nrow(stratified_out) > 0) {
  save_table(stratified_out, DIR_THIS, "rin_stratified_effects.csv")

  cat("\n=== RIN 分层内 g3_vs_g2 ===\n")
  print(stratified_out[, c("feature","stratum","logFC","p","n")],
        row.names = FALSE, digits = 4)

  cat("\nhUSI 在 RIN 三层内的效应:\n")
  print(stratified_out[stratified_out$feature == "hUSI",
                        c("stratum","logFC","p","n")],
        row.names = FALSE, digits = 4)
} else {
  log_warn("分层分析全失败，检查 grade 分布")
}

# ============================================================
# ============================================================
log_info("--- 8. 长基因 vs 短基因 RIN 敏感性 ---")

find_col <- function(nms, pat) {
  h <- grep(pat, nms, ignore.case = TRUE, value = TRUE)
  if (length(h) == 0) return(NA_character_)
  h[1]
}

gene_length_col <- find_col(colnames(anno), "length|gene_length|bp")
has_length <- !is.na(gene_length_col) && gene_length_col %in% colnames(anno)

if (has_length) {
  log_info("anno 有 gene_length，执行长度分析")
} else {
  log_warn("anno 里没有 gene_length 列，跳过长度分析")
}

# ============================================================
# ============================================================
log_info("--- 9. 稳健性总结 ---")

robust_summary <- data.frame()

for (sc in rownames(aging_sub)) {
  sc_eff <- all_effects[all_effects$feature == sc &
                          all_effects$contrast == "g3_vs_g2", ]
  robust_summary <- rbind(robust_summary, data.frame(
    feature = sc,
    model = sc_eff$model,
    logFC = sc_eff$logFC,
    p = sc_eff$p,
    robust_p05 = sc_eff$p < 0.05,
    stringsAsFactors = FALSE))
}

save_table(robust_summary, DIR_THIS, "rin_robust_summary.csv")

cat("\n=== 稳健性总结 ===\n")
print(robust_summary, row.names = FALSE, digits = 4)

cat("\n每个 feature 在 3 种 RIN 处理下 p<0.05 的模型数:\n")
robust_count <- aggregate(robust_p05 ~ feature, data = robust_summary,
                            FUN = sum, na.rm = TRUE)
print(robust_count, row.names = FALSE)

# ============================================================
# 10. figure
# ============================================================
log_info("--- 10. 图 ---")

clin_sub$grade_f <- factor(clin_sub$grade)
p1 <- ggplot(clin_sub, aes(x = grade_f, y = rin_v, fill = grade_f)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 0.3, alpha = 0.3) +
  theme_bw(base_size = 11) +
  labs(x = "AMD grade", y = "RIN",
       title = sprintf("RIN across grades (KW p = %.3g)",
                       kw_rin$p.value)) +
  theme(legend.position = "none")
save_fig(p1, DIR_THIS, "rin_grade_boxplot", width = 7, height = 5)

top_n <- 40
top_rin <- head(rin_cor[order(-abs(rin_cor$r_RIN)), ], top_n)
hm_mat <- matrix(top_rin$r_RIN, ncol = 1,
                  dimnames = list(top_rin$pathway, "r_RIN"))
pheatmap(hm_mat,
         cluster_rows = TRUE, cluster_cols = FALSE,
         display_numbers = TRUE, number_format = "%.2f",
         color = colorRampPalette(c("#4C72B0","#FFFFFF","#C44E52"))(50),
         main = "RIN correlation with pathways (top 40)",
         filename = file.path(DIR_THIS, "figures",
                              "rin_pathway_correlation.png"),
         width = 5, height = 9)

if (nrow(stratified_out) > 0) {
  focus_show <- intersect(c("hUSI","ML_Aging","SenMayo","CoreScence",
                              "Mito_CI_Subunits","OXPHOS","Mito_TCA"),
                            unique(stratified_out$feature))
  plot_df <- stratified_out[stratified_out$feature %in% focus_show, ]
  plot_df$feature <- factor(plot_df$feature, levels = focus_show)

  p3 <- ggplot(plot_df, aes(x = stratum, y = logFC, fill = stratum)) +
    geom_col(alpha = 0.8) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    facet_wrap(~ feature, scales = "free_y", ncol = 4) +
    theme_bw(base_size = 10) +
    labs(x = "RIN stratum", y = "logFC (g3_vs_g2)",
         title = "g3_vs_g2 effects within RIN strata") +
    theme(legend.position = "none")
  save_fig(p3, DIR_THIS, "rin_stratified_effects", width = 10, height = 6)
}

# ============================================================
# 11. summary
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("=== diag_rin_sensitivity.r 汇总 (v2) ===\n")
cat(strrep("=", 70), "\n", sep = "")

cat("\n【Q1】RIN 在 grade 间分布\n")
cat(sprintf("  KW p = %.4g, Spearman rho = %.3f\n",
            kw_rin$p.value, rho_rin_grade))

cat("\n【Q2】RIN 与协变量相关\n")
print(cross_out, row.names = FALSE)

cat("\n【Q3】hUSI 在 3 种 RIN 处理下\n")
print(hUSI_eff[, c("model","logFC","p","FDR")],
      row.names = FALSE, digits = 4)

cat("\n【Q4】4 衰老评分稳健性\n")
print(robust_count, row.names = FALSE)

cat("\n【Q5】hUSI 轴 top 100 中 g3_vs_g2 FDR<0.05 数量\n")
for (m in unique(axis_effects$model)) {
  sub <- axis_effects[axis_effects$model == m, ]
  n <- sum(sub$FDR < 0.05, na.rm = TRUE)
  cat(sprintf("  %-15s  %d / 100\n", m, n))
}

cat("\n【Q6】RIN 分层内 hUSI 效应\n")
if (nrow(stratified_out) > 0) {
  print(stratified_out[stratified_out$feature == "hUSI",
                        c("stratum","logFC","p","FDR","n")],
        row.names = FALSE, digits = 4)
} else {
  cat("  分层分析无结果\n")
}

cat("\n【Q7】hUSI vs RIN 直接相关\n")
cat(sprintf("  hUSI vs RIN: r = %.3f\n", r_hUSI_RIN))

cat("\n--- 输出位置 ---\n")
cat("目录:", DIR_THIS, "\n")

log_info("=== diag_rin_sensitivity.r (v2) done ===")

