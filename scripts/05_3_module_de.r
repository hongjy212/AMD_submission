#!/usr/bin/env Rscript
# ============================================================
#
#     "age+disease" → "Age + disease"
#
# v3 Fix：
#
#
# ============================================================

if (exists(".AMD_SETUP_LOADED", envir = .GlobalEnv)) {
  rm(".AMD_SETUP_LOADED", envir = .GlobalEnv)
}
source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(limma); library(splines); library(ggplot2)
  library(GSVA); library(pheatmap); library(tidyr)
})

DIR_THIS <- setup_output_dir("05_3_module_de")
log_info("=== 05_3_module_de.r (v4) start ===")

MIN_GENES <- 10
MAX_GENES <- 5000

# ============================================================
# ============================================================
log_info("--- 1. 载入 ---")

expr <- read_upstream("01_1b_batch_enhance", "data", "expr_wgcna.rds")
clin <- read_upstream("01_1b_batch_enhance", "data", "clinical_wgcna.rds")
module_df <- read_upstream("05_1_consensus_wgcna", "tables",
                            "consensus_module_genes_annotated.csv")
MEs_global <- read_upstream("05_2_module_trait", "data", "MEs_global.rds")
aging_scores <- read_upstream("04_1_aging_score", "data", "aging_scores.rds")
score_pathway <- read_upstream("04_2_pathway_score", "data",
                                "pathway_scores_ssgsea.rds")
anno <- read_upstream("01_1_prepare_bulk_data", "data", "gene_annotation.rds")

common <- intersect(colnames(expr), rownames(clin))
expr <- expr[, common, drop = FALSE]
clin <- clin[common, , drop = FALSE]

log_info("expr: {nrow(expr)} × {ncol(expr)}")
log_info("MEs_global: {nrow(MEs_global)} 样本 × {ncol(MEs_global)} 模块")

# ============================================================
# ============================================================
log_info("--- 2. 整理变量 ---")

clin$grade <- as.integer(clin$grade)
clin$age_v <- as.numeric(clin$age_n)
clin$rin_v <- as.numeric(clin$rin_n)
clin$pmi_v <- as.numeric(clin$pmi_n)
clin$batch_f <- factor(clin$batch)
sid <- as.character(clin$sample_id)
clin$layer <- factor(sub("^.*_(L[0-9]+)$", "\\1", sid),
                      levels = paste0("L", 1:8))
rownames(clin) <- colnames(expr)

log_info("grade: {paste(names(table(clin$grade)), table(clin$grade), sep='=', collapse=' ')}")
log_info("layer NA: {sum(is.na(clin$layer))}")

# ============================================================
# ============================================================
log_info("--- 3. 构建模块基因集 ---")

strip_version <- function(x) sub("\\..*$", "", x)
expr_ens <- expr
rownames(expr_ens) <- strip_version(rownames(expr))

all_modules <- setdiff(unique(module_df$Module), "grey")
log_info("真实模块: {length(all_modules)} 个")

gs_keep <- list()
gs_meta <- data.frame()
for (mod in all_modules) {
  mod_ens <- strip_version(module_df$ENSEMBL[module_df$Module == mod])
  mod_in_expr <- intersect(mod_ens, rownames(expr_ens))
  n <- length(mod_in_expr)
  if (n >= MIN_GENES && n <= MAX_GENES) {
    gs_keep[[mod]] <- mod_in_expr
    gs_meta <- rbind(gs_meta, data.frame(
      module = mod, n_genes = n, stringsAsFactors = FALSE))
  }
}
log_info("保留模块: {length(gs_keep)} / {length(all_modules)}")
print(gs_meta, row.names = FALSE)

# ============================================================
# ============================================================
log_info("--- 4. ssGSEA 打分 ---")

set.seed(SEED)
ssgsea_param <- gsvaParam(
  exprData = as.matrix(expr_ens),
  geneSets = gs_keep,
  kcdf = "Gaussian",
  minSize = MIN_GENES,
  maxSize = MAX_GENES
)
ssgsea_scores <- gsva(ssgsea_param)
log_info("ssGSEA: {nrow(ssgsea_scores)} × {ncol(ssgsea_scores)}")
save_rds(ssgsea_scores, DIR_THIS, "module_scores_ssgsea.rds")

# z-score crossvalidation
expr_z <- t(scale(t(expr_ens)))
expr_z[is.na(expr_z)] <- 0
zscore_mat <- do.call(rbind, lapply(gs_keep, function(genes) {
  colMeans(expr_z[genes, , drop = FALSE])
}))
save_rds(zscore_mat, DIR_THIS, "module_scores_zscore.rds")

common_gs <- intersect(rownames(ssgsea_scores), rownames(zscore_mat))
if (length(common_gs) >= 5) {
  r_cmp <- sapply(common_gs, function(gs) {
    cor(ssgsea_scores[gs, ], zscore_mat[gs, ])
  })
  save_table(data.frame(module = common_gs,
                          r_ssgsea_zscore = round(r_cmp, 3)),
              DIR_THIS, "ssgsea_vs_zscore_correlation.csv")
  log_info("ssGSEA vs zscore 中位相关: {round(median(r_cmp, na.rm=TRUE), 3)}")
}

# ============================================================
# ============================================================
log_info("--- 5. 设计矩阵 ---")

sp <- splines::ns(clin$age_v, df = 3)
covar_df <- data.frame(
  age_s1 = sp[, 1], age_s2 = sp[, 2], age_s3 = sp[, 3],
  sex_M = as.integer(clin$sex == "Male"),
  rin = clin$rin_v, pmi = clin$pmi_v,
  layer = clin$layer, batch = clin$batch_f,
  stringsAsFactors = FALSE)
rownames(covar_df) <- rownames(clin)

design_df <- data.frame(grade = factor(clin$grade), covar_df)
rownames(design_df) <- rownames(clin)
design_df <- design_df[complete.cases(design_df), , drop = FALSE]

design_multi <- model.matrix(
  ~ 0 + grade + age_s1 + age_s2 + age_s3 + sex_M + rin + pmi + layer + batch,
  data = design_df)
colnames(design_multi)[grep("^grade", colnames(design_multi))] <- paste0("g", 1:4)

log_info("design: {nrow(design_multi)} × {ncol(design_multi)}")

cm <- limma::makeContrasts(
  g2_vs_g1 = g2 - g1, g3_vs_g2 = g3 - g2, g4_vs_g3 = g4 - g3,
  g4_vs_g1 = g4 - g1,
  anyAMD_vs_g1 = (g2+g3+g4)/3 - g1,
  late_vs_early = (g3+g4)/2 - (g1+g2)/2,
  g4_vs_rest = g4 - (g1+g2+g3)/3,
  levels = design_multi)

# ============================================================
# ============================================================
run_limma_multi <- function(score_mat, method_name) {
  score_sub <- score_mat[, rownames(design_df), drop = FALSE]
  fit <- eBayes(contrasts.fit(lmFit(score_sub, design_multi), cm))
  out <- data.frame()
  for (ct in colnames(cm)) {
    tt <- topTable(fit, coef = ct, number = Inf, sort.by = "none")
    out <- rbind(out, data.frame(
      module = gsub("^ME", "", rownames(tt)),
      contrast = ct, method = method_name,
      logFC = tt$logFC, t = tt$t, p = tt$P.Value, FDR = tt$adj.P.Val,
      stringsAsFactors = FALSE))
  }
  out
}

# ============================================================
# ============================================================
run_poly3F <- function(score_mat, method_name) {
  poly_mat <- poly(clin$grade, 3)
  colnames(poly_mat) <- c("grade_p1", "grade_p2", "grade_p3")
  rownames(poly_mat) <- rownames(clin)
  design_poly_df <- data.frame(poly_mat, covar_df)
  rownames(design_poly_df) <- rownames(clin)
  design_poly_df <- design_poly_df[complete.cases(design_poly_df), , drop = FALSE]
  score_poly <- score_mat[, rownames(design_poly_df), drop = FALSE]

  design_poly <- model.matrix(
    ~ grade_p1 + grade_p2 + grade_p3 + age_s1 + age_s2 + age_s3 +
      sex_M + rin + pmi + layer + batch,
    data = design_poly_df)
  fit_poly <- eBayes(lmFit(score_poly, design_poly))

  coef_poly <- grep("^grade_p", colnames(design_poly), value = TRUE)
  tt_F <- topTable(fit_poly, coef = coef_poly, number = Inf, sort.by = "none")
  tt_p1 <- topTable(fit_poly, coef = "grade_p1", number = Inf, sort.by = "none")
  tt_p2 <- topTable(fit_poly, coef = "grade_p2", number = Inf, sort.by = "none")
  tt_p3 <- topTable(fit_poly, coef = "grade_p3", number = Inf, sort.by = "none")

  out <- data.frame(
    module = gsub("^ME", "", rownames(tt_F)),
    method = method_name,
    F_stat = tt_F$F, p = tt_F$P.Value, FDR = tt_F$adj.P.Val,
    p1_logFC = tt_p1$logFC, p1_p = tt_p1$P.Value, p1_FDR = tt_p1$adj.P.Val,
    p2_logFC = tt_p2$logFC, p2_p = tt_p2$P.Value, p2_FDR = tt_p2$adj.P.Val,
    p3_logFC = tt_p3$logFC, p3_p = tt_p3$P.Value, p3_FDR = tt_p3$adj.P.Val,
    stringsAsFactors = FALSE)
  out[order(out$p), ]
}

# ============================================================
# ============================================================
run_age_class <- function(score_mat, method_name) {
  out <- data.frame()
  for (i in seq_len(nrow(score_mat))) {
    sc_name <- gsub("^ME", "", rownames(score_mat)[i])
    y <- score_mat[i, ]
    df <- data.frame(y = y, grade = clin$grade, age = clin$age_v,
                     sex = clin$sex, rin = clin$rin_v, pmi = clin$pmi_v,
                     layer = clin$layer, batch = clin$batch_f)
    df <- df[complete.cases(df), ]
    sp_g <- splines::ns(df$age, df = 3)
    fit1 <- lm(y ~ poly(grade, 3) + sp_g + sex + rin + pmi + layer + batch, data = df)
    p_g <- anova(fit1)["poly(grade, 3)", "Pr(>F)"]
    fit2 <- lm(y ~ sp_g + poly(grade, 3) + sex + rin + pmi + layer + batch, data = df)
    p_a <- anova(fit2)["sp_g", "Pr(>F)"]
    out <- rbind(out, data.frame(
      module = sc_name, method = method_name, p_grade = p_g, p_age = p_a,
      stringsAsFactors = FALSE))
  }
  out$type <- "Unrelated"
  out$type[out$p_grade < 0.05 & out$p_age > 0.10] <- "Pure disease"
  out$type[out$p_grade < 0.05 & out$p_age < 0.05] <- "Age + disease"
  out$type[out$p_grade > 0.10 & out$p_age < 0.05] <- "Pure age"
  out
}

# ============================================================
# ============================================================
log_info("--- 9. limma multi（ssGSEA + ME）---")

multi_ssgsea <- run_limma_multi(ssgsea_scores, "ssGSEA")
multi_ME <- run_limma_multi(t(MEs_global), "ME")
save_table(multi_ssgsea, DIR_THIS, "module_limma_multi_ssgsea.csv")
save_table(multi_ME, DIR_THIS, "module_limma_multi_ME.csv")

n_mod_ss <- nrow(ssgsea_scores)
n_mod_me <- ncol(MEs_global)

cat("\n=== ssGSEA limma multi FDR<0.05 ===\n")
for (ct in colnames(cm)) {
  n <- sum(multi_ssgsea$FDR[multi_ssgsea$contrast == ct] < 0.05, na.rm = TRUE)
  cat(sprintf("  %-15s  %d / %d\n", ct, n, n_mod_ss))
}

cat("\n=== ME limma multi FDR<0.05 ===\n")
for (ct in colnames(cm)) {
  n <- sum(multi_ME$FDR[multi_ME$contrast == ct] < 0.05, na.rm = TRUE)
  cat(sprintf("  %-15s  %d / %d\n", ct, n, n_mod_me))
}

g3_ss <- multi_ssgsea[multi_ssgsea$contrast == "g3_vs_g2", ]
g3_me <- multi_ME[multi_ME$contrast == "g3_vs_g2", ]
cmp <- merge(g3_ss[, c("module","logFC","p","FDR")],
              g3_me[, c("module","logFC","p","FDR")],
              by = "module", suffixes = c("_ssGSEA","_ME"))
cmp$consistent <- (cmp$FDR_ssGSEA < 0.05) & (cmp$FDR_ME < 0.05)
cmp <- cmp[order(cmp$FDR_ssGSEA), ]
save_table(cmp, DIR_THIS, "module_limma_multi_compare.csv")

cat("\n=== g3_vs_g2 两方法对比（按 ssGSEA FDR 排序）===\n")
print(cmp, row.names = FALSE, digits = 4)

cat("\n两方法都 FDR<0.05 的模块:\n")
consistent_mods <- cmp$module[cmp$consistent]
if (length(consistent_mods) > 0) {
  print(consistent_mods)
} else {
  cat("  无\n")
}

# ============================================================
# 10. poly3_F
# ============================================================
log_info("--- 10. poly3_F ---")

polyF_ssgsea <- run_poly3F(ssgsea_scores, "ssGSEA")
polyF_ME <- run_poly3F(t(MEs_global), "ME")
save_table(polyF_ssgsea, DIR_THIS, "module_poly3F_ssgsea.csv")
save_table(polyF_ME, DIR_THIS, "module_poly3F_ME.csv")

cat("\n=== ssGSEA poly3_F ===\n")
print(polyF_ssgsea[, c("module","F_stat","p","FDR",
                         "p1_FDR","p2_FDR","p3_FDR")],
      row.names = FALSE, digits = 4)

cat("\n=== ME poly3_F ===\n")
print(polyF_ME[, c("module","F_stat","p","FDR",
                     "p1_FDR","p2_FDR","p3_FDR")],
      row.names = FALSE, digits = 4)

cat("\n=== p1/p2/p3 FDR<0.05 数量（两方法）===\n")
for (m_name in c("ssGSEA","ME")) {
  pf <- if (m_name == "ssGSEA") polyF_ssgsea else polyF_ME
  cat(sprintf("  %-8s p1=%d  p2=%d  p3=%d\n", m_name,
              sum(pf$p1_FDR < 0.05, na.rm = TRUE),
              sum(pf$p2_FDR < 0.05, na.rm = TRUE),
              sum(pf$p3_FDR < 0.05, na.rm = TRUE)))
}

pf_cmp <- merge(
  polyF_ssgsea[, c("module","p3_logFC","p3_FDR")],
  polyF_ME[, c("module","p3_logFC","p3_FDR")],
  by = "module", suffixes = c("_ssGSEA","_ME"))
pf_cmp$p3_consistent <- (pf_cmp$p3_FDR_ssGSEA < 0.05) & (pf_cmp$p3_FDR_ME < 0.05)
pf_cmp$p3_direction_consistent <- sign(pf_cmp$p3_logFC_ssGSEA) ==
                                    sign(pf_cmp$p3_logFC_ME)
pf_cmp <- pf_cmp[order(pf_cmp$p3_FDR_ssGSEA), ]
save_table(pf_cmp, DIR_THIS, "module_poly3F_compare.csv")

cat("\n=== poly3_F p3 两方法对比 ===\n")
print(pf_cmp, row.names = FALSE, digits = 4)

# ============================================================
# ============================================================
log_info("--- 11. 年龄双模型分类 ---")

age_class_ssgsea <- run_age_class(ssgsea_scores, "ssGSEA")
age_class_ME <- run_age_class(t(MEs_global), "ME")
save_table(age_class_ssgsea, DIR_THIS, "module_age_classification_ssgsea.csv")
save_table(age_class_ME, DIR_THIS, "module_age_classification_ME.csv")

cat("\n=== ssGSEA 年龄双模型分类 ===\n")
print(table(age_class_ssgsea$type))
pure_ss <- age_class_ssgsea[age_class_ssgsea$type == "Pure disease", ]
if (nrow(pure_ss) > 0) {
  cat("\n纯疾病模块（ssGSEA）:\n")
  print(pure_ss, row.names = FALSE, digits = 4)
}

cat("\n=== ME 年龄双模型分类 ===\n")
print(table(age_class_ME$type))
pure_me <- age_class_ME[age_class_ME$type == "Pure disease", ]
if (nrow(pure_me) > 0) {
  cat("\n纯疾病模块（ME）:\n")
  print(pure_me, row.names = FALSE, digits = 4)
}

age_cmp <- merge(
  age_class_ssgsea[, c("module","type")],
  age_class_ME[, c("module","type")],
  by = "module", suffixes = c("_ssGSEA","_ME"))
age_cmp$type_consistent <- age_cmp$type_ssGSEA == age_cmp$type_ME
save_table(age_cmp, DIR_THIS, "module_age_classification_compare.csv")

cat("\n年龄分类一致的模块数:",
    sum(age_cmp$type_consistent), "/", nrow(age_cmp), "\n")

# ============================================================
# ============================================================
log_info("--- 12. 与生物学轴相关 ---")

common_aging <- intersect(colnames(ssgsea_scores), colnames(aging_scores))
hUSI_v <- aging_scores["hUSI", common_aging]
ML_v <- aging_scores["ML_Aging", common_aging]

GENES_16 <- c("FCF1","UGGT1","CCDC3","S100B","BPNT1","BRI3BP","BTBD19",
              "DPH7","AKR1B1","ZNF623","ZNF587","RSPH3","TAF3","NDUFB8",
              "ZNF417","C6orf120")
sym_map <- setNames(toupper(anno$symbol), strip_version(anno$ensembl))
expr_row_sym <- sym_map[strip_version(rownames(expr))]
genes_16_in <- rownames(expr)[expr_row_sym %in% GENES_16]
expr_z_16 <- t(scale(t(expr[genes_16_in, , drop = FALSE])))
expr_z_16[is.na(expr_z_16)] <- 0
score_16 <- colMeans(expr_z_16)

mito_pw <- intersect(c("Mito_CI_Subunits","Mito_CIII_Subunits",
                        "Mito_CIV_Subunits","Mito_CV_Subunits",
                        "Mito_Complex_V","Mito_TCA","OXPHOS"),
                      rownames(score_pathway))
mito_score <- colMeans(score_pathway[mito_pw, , drop = FALSE])

cor_axis <- data.frame()
for (mod in rownames(ssgsea_scores)) {
  common_h <- intersect(colnames(ssgsea_scores), common_aging)
  r_h <- cor(ssgsea_scores[mod, common_h], hUSI_v, use = "complete.obs")
  r_m <- cor(ssgsea_scores[mod, common_h], ML_v, use = "complete.obs")

  common_16 <- intersect(colnames(ssgsea_scores), names(score_16))
  r_16 <- cor(ssgsea_scores[mod, common_16], score_16[common_16], use = "complete.obs")

  common_mt <- intersect(colnames(ssgsea_scores), names(mito_score))
  r_mt <- cor(ssgsea_scores[mod, common_mt], mito_score[common_mt], use = "complete.obs")

  cor_axis <- rbind(cor_axis, data.frame(
    module = mod,
    r_hUSI = round(r_h, 3),
    r_ML_Aging = round(r_m, 3),
    r_genes16 = round(r_16, 3),
    r_mito = round(r_mt, 3),
    stringsAsFactors = FALSE))
}
cor_axis <- cor_axis[order(-abs(cor_axis$r_mito)), ]
save_table(cor_axis, DIR_THIS, "module_hUSI_correlation.csv")

cat("\n=== 模块 vs 生物学轴 ===\n")
print(cor_axis, row.names = FALSE)

# ============================================================
# ============================================================
log_info("--- 13. 模块汇总表 ---")

summary_df <- data.frame()
for (mod in rownames(ssgsea_scores)) {
  g3_ss <- multi_ssgsea[multi_ssgsea$module == mod &
                         multi_ssgsea$contrast == "g3_vs_g2", ]
  g3_me <- multi_ME[multi_ME$module == mod &
                     multi_ME$contrast == "g3_vs_g2", ]
  pf_ss <- polyF_ssgsea[polyF_ssgsea$module == mod, ]
  pf_me <- polyF_ME[polyF_ME$module == mod, ]
  ac_ss <- age_class_ssgsea[age_class_ssgsea$module == mod, ]
  ax_row <- cor_axis[cor_axis$module == mod, ]

  summary_df <- rbind(summary_df, data.frame(
    module = mod,
    n_genes = sum(module_df$Module == mod),

    g3v2_logFC_ssGSEA = if (nrow(g3_ss) > 0) as.numeric(g3_ss$logFC) else NA_real_,
    g3v2_FDR_ssGSEA = if (nrow(g3_ss) > 0) as.numeric(g3_ss$FDR) else NA_real_,
    g3v2_logFC_ME = if (nrow(g3_me) > 0) as.numeric(g3_me$logFC) else NA_real_,
    g3v2_FDR_ME = if (nrow(g3_me) > 0) as.numeric(g3_me$FDR) else NA_real_,

    poly3F_p3_FDR_ssGSEA = if (nrow(pf_ss) > 0) as.numeric(pf_ss$p3_FDR) else NA_real_,
    poly3F_p3_FDR_ME = if (nrow(pf_me) > 0) as.numeric(pf_me$p3_FDR) else NA_real_,

    age_class_ssGSEA = if (nrow(ac_ss) > 0) as.character(ac_ss$type) else NA_character_,

    r_hUSI = if (nrow(ax_row) > 0) as.numeric(ax_row$r_hUSI) else NA_real_,
    r_ML_Aging = if (nrow(ax_row) > 0) as.numeric(ax_row$r_ML_Aging) else NA_real_,
    r_genes16 = if (nrow(ax_row) > 0) as.numeric(ax_row$r_genes16) else NA_real_,
    r_mito = if (nrow(ax_row) > 0) as.numeric(ax_row$r_mito) else NA_real_,
    stringsAsFactors = FALSE))
}
summary_df <- summary_df[order(summary_df$g3v2_FDR_ssGSEA), ]
save_table(summary_df, DIR_THIS, "module_summary_all.csv")

cat("\n=== 模块完整汇总 ===\n")
print(summary_df, row.names = FALSE, digits = 4)

# ============================================================
# 14. figure
# ============================================================
log_info("--- 14. 图 ---")

top16 <- polyF_ssgsea
top16$module <- factor(top16$module,
                        levels = top16$module[order(top16$F_stat)])
top16$sig <- ifelse(top16$FDR < 0.05, "FDR<0.05",
             ifelse(top16$p < 0.05, "p<0.05", "ns"))

p1 <- ggplot(top16, aes(x = module, y = F_stat, fill = sig)) +
  geom_col() + coord_flip() + theme_bw(base_size = 10) +
  scale_fill_manual(values = c("FDR<0.05" = "#C44E52",
                                "p<0.05" = "#FDB462",
                                "ns" = "#B0B0B0")) +
  labs(x = "", y = "poly3_F F statistic", fill = "",
       title = "Module poly3_F (ssGSEA, non-monotonic)")
save_fig(p1, DIR_THIS, "module_poly3F_bar", width = 7, height = 5)

top8 <- head(polyF_ssgsea$module, 8)
curve_data <- data.frame()
for (mod in top8) {
  for (g in 1:4) {
    idx <- clin$grade == g
    curve_data <- rbind(curve_data, data.frame(
      module = mod, grade = g,
      mean_score = mean(ssgsea_scores[mod, idx], na.rm = TRUE),
      se = sd(ssgsea_scores[mod, idx], na.rm = TRUE) / sqrt(sum(idx)),
      stringsAsFactors = FALSE))
  }
}

p2 <- ggplot(curve_data, aes(x = grade, y = mean_score,
                               color = module, group = module)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_score - se, ymax = mean_score + se),
                 width = 0.15, alpha = 0.5) +
  theme_bw(base_size = 11) +
  scale_x_continuous(breaks = 1:4, labels = paste0("G", 1:4)) +
  labs(x = "AMD grade", y = "ssGSEA score (mean ± SE)",
       title = "Top 8 modules by poly3_F", color = "Module")
save_fig(p2, DIR_THIS, "module_grade_curves", width = 10, height = 6)

heat_mat <- matrix(NA_real_, nrow = nrow(ssgsea_scores), ncol = 4,
                    dimnames = list(rownames(ssgsea_scores), paste0("G", 1:4)))
for (mod in rownames(ssgsea_scores)) {
  for (g in 1:4) {
    heat_mat[mod, g] <- mean(ssgsea_scores[mod, clin$grade == g], na.rm = TRUE)
  }
}
pheatmap(heat_mat, cluster_rows = TRUE, cluster_cols = FALSE,
         display_numbers = TRUE, number_format = "%.2f",
         color = colorRampPalette(c("#4C72B0","#FFFFFF","#C44E52"))(50),
         main = "Modules × grade (ssGSEA mean)",
         filename = file.path(DIR_THIS, "figures", "module_heatmap.png"),
         width = 7, height = 7)

# ============================================================
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("=== 05_3 汇总 (v4) ===\n")
cat(strrep("=", 70), "\n", sep = "")

cat("\n【Q1】limma multi g3_vs_g2 FDR<0.05\n")
cat(sprintf("  ssGSEA: %d / %d\n",
            sum(multi_ssgsea$FDR[multi_ssgsea$contrast == "g3_vs_g2"] < 0.05,
                na.rm = TRUE), n_mod_ss))
cat(sprintf("  ME:     %d / %d\n",
            sum(multi_ME$FDR[multi_ME$contrast == "g3_vs_g2"] < 0.05,
                na.rm = TRUE), n_mod_me))

cat("\n  两方法都显著的模块:\n")
if (length(consistent_mods) > 0) {
  print(consistent_mods)
} else {
  cat("  无\n")
}

cat("\n  g3_vs_g2 对比明细（前 6 行）:\n")
print(head(cmp, 6), row.names = FALSE, digits = 4)

cat("\n【Q2】poly3_F p3_FDR<0.05\n")
cat(sprintf("  ssGSEA: %d\n", sum(polyF_ssgsea$p3_FDR < 0.05, na.rm = TRUE)))
cat(sprintf("  ME:     %d\n", sum(polyF_ME$p3_FDR < 0.05, na.rm = TRUE)))

cat("\n  两方法 p3 都显著的模块:\n")
p3_consistent <- pf_cmp$module[pf_cmp$p3_consistent]
if (length(p3_consistent) > 0) {
  print(p3_consistent)
} else {
  cat("  无\n")
}

cat("\n  p3 明细（前 6 行）:\n")
print(head(pf_cmp, 6), row.names = FALSE, digits = 4)

cat("\n【Q3】年龄双模型分类\n")
cat(sprintf("  ssGSEA pure disease: %d\n",
            sum(age_class_ssgsea$type == "Pure disease")))
cat(sprintf("  ME     pure disease: %d\n",
            sum(age_class_ME$type == "Pure disease")))
cat(sprintf("  分类一致:     %d / %d\n",
            sum(age_cmp$type_consistent), nrow(age_cmp)))

cat("\n【Q4】与生物学轴（|r|>0.4）\n")
strong <- cor_axis[abs(cor_axis$r_hUSI) > 0.4 |
                    abs(cor_axis$r_ML_Aging) > 0.4 |
                    abs(cor_axis$r_genes16) > 0.4 |
                    abs(cor_axis$r_mito) > 0.4, ]
print(strong, row.names = FALSE)

cat("\n--- 输出位置 ---\n")
cat("目录:", DIR_THIS, "\n")

log_info("=== 05_3 (v4) done ===")

