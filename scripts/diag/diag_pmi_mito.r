#!/usr/bin/env Rscript
# ============================================================
# diag_pmi_mito.r
#
# Input：
#   results/01_1_prepare_bulk_data/data/{expr_corrected,clinical,gene_annotation}.rds
#   results/01_2_prepare_gene_sets/data/gene_sets.rds
#   results/04_2_pathway_score/data/pathway_scores_ssgsea.rds
#
#   tables/pmi_grade_distribution.csv
#   tables/pmi_covariate_correlation.csv
#   tables/pmi_pathway_correlation.csv
#   tables/pmi_adjusted_effects.csv
#   tables/mtdna_vs_nuc_effects.csv
#   tables/pmi_robust_pathways.csv
#   figures/pmi_grade_boxplot.pdf/.png
#   figures/pmi_pathway_correlation.pdf/.png
#   figures/mtdna_vs_nuc_curves.pdf/.png
# ============================================================

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(limma); library(splines); library(ggplot2)
  library(tidyr); library(pheatmap)
})

DIR_THIS <- setup_output_dir("diag_pmi_mito")
log_info("=== diag_pmi_mito.r start ===")

# ============================================================
# ============================================================
log_info("--- 1. 载入数据 ---")

expr  <- read_upstream("01_1_prepare_bulk_data", "data", "expr_corrected.rds")
clin  <- read_upstream("01_1_prepare_bulk_data", "data", "clinical.rds")
anno  <- read_upstream("01_1_prepare_bulk_data", "data", "gene_annotation.rds")
GENE_SETS <- read_upstream("01_2_prepare_gene_sets", "data", "gene_sets.rds")
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
rownames(clin) <- colnames(expr)

log_info("expr: {nrow(expr)} × {ncol(expr)}")
log_info("score_mat: {nrow(score_mat)} × {ncol(score_mat)}")

find_col <- function(nms, pat) {
  h <- grep(pat, nms, ignore.case = TRUE, value = TRUE)
  if (length(h) == 0) return(NA_character_)
  h[1]
}
col_pmi <- find_col(colnames(clin), "^pmi$|pmi_n|post.?mortem")
col_rin <- find_col(colnames(clin), "^rin$|rin_n|RIN")
col_age <- find_col(colnames(clin), "age_n|^age$")
log_info("PMI 列: {col_pmi}")
log_info("RIN 列: {col_rin}")
log_info("Age 列: {col_age}")

if (is.na(col_pmi)) {
  stop("未找到 PMI 列，无法继续")
}

clin$pmi_v <- as.numeric(clin[[col_pmi]])
if (!is.na(col_rin)) clin$rin_v <- as.numeric(clin[[col_rin]])
clin$age_v <- as.numeric(clin[[col_age]])

# ============================================================
# ============================================================
log_info("--- 2. PMI 在 grade 间分布 ---")

pmi_dist <- do.call(rbind, lapply(sort(unique(clin$grade)), function(g) {
  idx <- clin$grade == g
  x <- clin$pmi_v[idx]
  data.frame(
    grade = g, n = sum(idx),
    mean = mean(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    median = median(x, na.rm = TRUE),
    q25 = quantile(x, 0.25, na.rm = TRUE),
    q75 = quantile(x, 0.75, na.rm = TRUE),
    stringsAsFactors = FALSE)
}))
save_table(pmi_dist, DIR_THIS, "pmi_grade_distribution.csv")
cat("\n=== PMI 在 grade 间的分布 ===\n")
print(pmi_dist, row.names = FALSE, digits = 3)

if (length(unique(clin$grade)) >= 3) {
  kw <- kruskal.test(clin$pmi_v ~ clin$grade)
  sp_r <- cor(clin$pmi_v, clin$grade, method = "spearman",
              use = "complete.obs")
  cat(sprintf("\nKruskal-Wallis p = %.4g\n", kw$p.value))
  cat(sprintf("Spearman rho(PMI, grade) = %.3f\n", sp_r))
}

covars_to_check <- c("age_v", if (!is.na(col_rin)) "rin_v")
cor_out <- data.frame()
for (cv in covars_to_check) {
  if (!cv %in% colnames(clin)) next
  r <- cor(clin$pmi_v, clin[[cv]], method = "spearman",
           use = "complete.obs")
  cor_out <- rbind(cor_out, data.frame(
    var1 = "PMI", var2 = cv, rho = round(r, 3)))
}
# PMI vs sex
if ("sex" %in% colnames(clin)) {
  wt <- wilcox.test(clin$pmi_v ~ clin$sex)
  cor_out <- rbind(cor_out, data.frame(
    var1 = "PMI", var2 = "sex(Wilcox p)",
    rho = round(wt$p.value, 4)))
}
# PMI vs batch
if ("batch" %in% colnames(clin)) {
  kw2 <- kruskal.test(clin$pmi_v ~ clin$batch)
  cor_out <- rbind(cor_out, data.frame(
    var1 = "PMI", var2 = "batch(KW p)",
    rho = round(kw2$p.value, 4)))
}
save_table(cor_out, DIR_THIS, "pmi_covariate_correlation.csv")
cat("\n=== PMI 与协变量相关 ===\n")
print(cor_out, row.names = FALSE)

# ============================================================
# ============================================================
log_info("--- 3. PMI vs 通路评分 ---")

common_samp <- intersect(colnames(score_mat), rownames(clin))
pmi_v <- clin[common_samp, "pmi_v"]
rin_v <- if (!is.na(col_rin)) clin[common_samp, "rin_v"] else NULL
age_v <- clin[common_samp, "age_v"]

pmi_cor <- data.frame()
for (gs in rownames(score_mat)) {
  y <- score_mat[gs, common_samp]
  r_pmi <- cor(y, pmi_v, method = "spearman", use = "complete.obs")
  r_age <- cor(y, age_v, method = "spearman", use = "complete.obs")
  r_rin <- if (!is.null(rin_v)) {
    cor(y, rin_v, method = "spearman", use = "complete.obs")
  } else NA_real_
  pmi_cor <- rbind(pmi_cor, data.frame(
    pathway = gs,
    r_pmi = round(r_pmi, 3),
    r_age = round(r_age, 3),
    r_rin = round(r_rin, 3),
    stringsAsFactors = FALSE))
}

gs_info_path <- file.path(DIR_RESULT, "04_2_pathway_score", "tables",
                          "gene_sets_used.csv")
if (file.exists(gs_info_path)) {
  gs_used <- read.csv(gs_info_path, stringsAsFactors = FALSE)
  pmi_cor$group <- gs_used$group[match(pmi_cor$pathway, gs_used$geneset)]
}
pmi_cor <- pmi_cor[order(-abs(pmi_cor$r_pmi)), ]
save_table(pmi_cor, DIR_THIS, "pmi_pathway_correlation.csv")

cat("\n=== PMI 相关 top 25 ===\n")
print(head(pmi_cor[, c("pathway", "group", "r_pmi", "r_age", "r_rin")], 25),
      row.names = FALSE)

if ("group" %in% colnames(pmi_cor)) {
  cat("\n=== |r_pmi| 按 group 中位数 ===\n")
  grp_med <- aggregate(abs(r_pmi) ~ group, data = pmi_cor, FUN = median)
  grp_med <- grp_med[order(-grp_med$`abs(r_pmi)`), ]
  print(grp_med, row.names = FALSE)
}

# ============================================================
# ============================================================
log_info("--- 4. 控 PMI 前后效应变化 ---")

run_model <- function(score_mat, clin, model_name,
                       use_age = TRUE, use_pmi = FALSE,
                       use_rin = FALSE, pmi_spline = FALSE) {
  covar <- data.frame(
    sex_M = as.integer(clin$sex == "Male"),
    stringsAsFactors = FALSE)
  # PMI
  if (use_pmi) {
    if (pmi_spline) {
      sp <- splines::ns(clin$pmi_v, df = 3)
      colnames(sp) <- c("pmi_s1", "pmi_s2", "pmi_s3")
      covar <- cbind(covar, sp)
    } else {
      covar$pmi <- clin$pmi_v
    }
  }
  # RIN
  if (use_rin && "rin_v" %in% colnames(clin)) {
    covar$rin <- clin$rin_v
  }
  # Age
  if (use_age) {
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
                          paste(setdiff(colnames(covar), ""),
                                collapse = " + ")))
  design <- model.matrix(fml, data = design_df)
  colnames(design)[grep("^grade", colnames(design))] <- paste0("g", 1:4)

  cm <- limma::makeContrasts(
    g3_vs_g2 = g3 - g2,
    g4_vs_g3 = g4 - g3,
    levels = design)
  fit <- eBayes(contrasts.fit(lmFit(score_sub, design), cm))

  out <- data.frame()
  for (ct in colnames(cm)) {
    tt <- topTable(fit, coef = ct, number = Inf, sort.by = "none")
    out <- rbind(out, data.frame(
      pathway = rownames(tt), model = model_name, contrast = ct,
      logFC = tt$logFC, p = tt$P.Value, FDR = tt$adj.P.Val,
      stringsAsFactors = FALSE))
  }
  out
}

all_models <- rbind(
  run_model(score_mat, clin, "M1_age",
            use_age = TRUE, use_pmi = FALSE, use_rin = FALSE),
  run_model(score_mat, clin, "M2_age_PMI_linear",
            use_age = TRUE, use_pmi = TRUE, use_rin = FALSE),
  run_model(score_mat, clin, "M3_age_PMI_spline",
            use_age = TRUE, use_pmi = TRUE, use_rin = FALSE,
            pmi_spline = TRUE),
  run_model(score_mat, clin, "M4_age_PMI_spline_RIN",
            use_age = TRUE, use_pmi = TRUE, use_rin = TRUE,
            pmi_spline = TRUE)
)

if (exists("gs_used")) {
  all_models$group <- gs_used$group[match(all_models$pathway, gs_used$geneset)]
}
save_table(all_models, DIR_THIS, "pmi_adjusted_effects.csv")

cat("\n=== g3_vs_g2：控 PMI 前后 FDR<0.05 数量 ===\n")
for (m in unique(all_models$model)) {
  n <- sum(all_models$FDR[all_models$model == m &
                          all_models$contrast == "g3_vs_g2"] < 0.05,
           na.rm = TRUE)
  cat(sprintf("  %-25s  %d\n", m, n))
}

m1 <- all_models[all_models$model == "M1_age" &
                 all_models$contrast == "g3_vs_g2", c("pathway","logFC","FDR")]
m4 <- all_models[all_models$model == "M4_age_PMI_spline_RIN" &
                 all_models$contrast == "g3_vs_g2", c("pathway","logFC","FDR")]
cmp <- merge(m1, m4, by = "pathway", suffixes = c("_M1","_M4"))
cmp$pmi_driven <- cmp$FDR_M1 < 0.05 & cmp$FDR_M4 >= 0.05
cmp$pmi_robust <- cmp$FDR_M1 < 0.05 & cmp$FDR_M4 < 0.05
cmp <- cmp[order(cmp$FDR_M1), ]

cat("\n=== g3_vs_g2：PMI 驱动 vs PMI 稳健 ===\n")
cat("  PMI 驱动（M1 显著、M4 不显著）:", sum(cmp$pmi_driven), "\n")
cat("  PMI 稳健（两者都显著）:", sum(cmp$pmi_robust), "\n")

if (sum(cmp$pmi_driven) > 0) {
  cat("\n  -- PMI 驱动 pathway --\n")
  print(cmp[cmp$pmi_driven, ], row.names = FALSE, digits = 4)
}
if (sum(cmp$pmi_robust) > 0) {
  cat("\n  -- PMI 稳健 pathway --\n")
  print(cmp[cmp$pmi_robust, ], row.names = FALSE, digits = 4)
}

# ============================================================
# ============================================================
log_info("--- 5. mtDNA vs 核编码线粒体 ---")

MTDNA_GENES <- c("MT-ND1","MT-ND2","MT-ND3","MT-ND4","MT-ND4L",
                 "MT-ND5","MT-ND6","MT-CYB","MT-CO1","MT-CO2",
                 "MT-CO3","MT-ATP6","MT-ATP8")
MTDNA_ALT <- c("ND1","ND2","ND3","ND4","ND4L","ND5","ND6",
               "CYTB","CO1","CO2","CO3","ATP6","ATP8")

mt_in <- intersect(MTDNA_GENES, rownames(expr))
mt_alt_in <- intersect(MTDNA_ALT, rownames(expr))
log_info("mtDNA 编码基因命中: {length(mt_in)}")
log_info("mtDNA 简写命中: {length(mt_alt_in)}")
mt_all <- unique(c(mt_in, mt_alt_in))
cat("  mtDNA 基因:", paste(mt_all, collapse = ", "), "\n")

mito_genes <- unique(unlist(lapply(GENE_SETS[
  grepl("^Mito_", names(GENE_SETS))], function(x) x$genes)))
nuc_mito <- setdiff(mito_genes, mt_all)
nuc_mito_in <- intersect(nuc_mito, rownames(expr))
log_info("核编码线粒体基因命中: {length(nuc_mito_in)}")

non_mito <- setdiff(rownames(expr), c(mt_all, mito_genes))
set.seed(SEED)
non_mito_sample <- sample(non_mito, min(500, length(non_mito)))

calc_score <- function(genes, expr, common_samp) {
  genes_in <- intersect(genes, rownames(expr))
  if (length(genes_in) < 3) return(NULL)
  expr_z <- t(scale(t(expr[genes_in, common_samp, drop = FALSE])))
  expr_z[is.na(expr_z)] <- 0
  colMeans(expr_z)
}

s_mt   <- calc_score(mt_all, expr, common_samp)
s_nuc  <- calc_score(nuc_mito_in, expr, common_samp)
s_rand <- calc_score(non_mito_sample, expr, common_samp)
s_all  <- list(mtDNA_encoded = s_mt,
               nuclear_mito  = s_nuc,
               non_mito_ctrl = s_rand)

mt_cor <- data.frame()
for (nm in names(s_all)) {
  s <- s_all[[nm]]
  if (is.null(s)) next
  mt_cor <- rbind(mt_cor, data.frame(
    geneset = nm,
    n_genes = length(intersect(
      switch(nm,
        mtDNA_encoded = mt_all,
        nuclear_mito = nuc_mito_in,
        non_mito_ctrl = non_mito_sample),
      rownames(expr))),
    r_pmi = round(cor(s, pmi_v, method = "spearman",
                       use = "complete.obs"), 3),
    r_rin = if (!is.null(rin_v)) round(cor(s, rin_v, method = "spearman",
                                            use = "complete.obs"), 3) else NA,
    r_age = round(cor(s, age_v, method = "spearman",
                       use = "complete.obs"), 3),
    stringsAsFactors = FALSE))
}
save_table(mt_cor, DIR_THIS, "mtdna_vs_nuc_effects.csv")
cat("\n=== mtDNA vs 核编码 vs 非线粒体对照 ===\n")
print(mt_cor, row.names = FALSE)

run_one <- function(score_v, clin, common_samp) {
  df <- data.frame(
    y = score_v[common_samp],
    grade = factor(clin[common_samp, "grade"]),
    age = clin[common_samp, "age_v"],
    sex = clin[common_samp, "sex"],
    pmi = clin[common_samp, "pmi_v"])
  if ("rin_v" %in% colnames(clin)) df$rin <- clin[common_samp, "rin_v"]
  df <- df[complete.cases(df), ]
  sp <- splines::ns(df$age, df = 3)
  if ("rin" %in% colnames(df)) {
    fit <- lm(y ~ grade + sp + sex + pmi + rin, data = df)
    fit2 <- lm(y ~ grade + sp + sex + rin, data = df)  # 不控 PMI
  } else {
    fit <- lm(y ~ grade + sp + sex + pmi, data = df)
    fit2 <- lm(y ~ grade + sp + sex, data = df)
  }
  # g3_vs_g2
  L <- rep(0, length(coef(fit)))
  names(L) <- names(coef(fit))
  if ("grade3" %in% names(L) && "grade2" %in% names(L)) {
    L["grade3"] <- 1; L["grade2"] <- -1
  }
  est <- sum(L * coef(fit))
  se <- sqrt(L %*% vcov(fit) %*% L)
  p_ctrl_pmi <- 2 * pnorm(-abs(est / se))

  L2 <- rep(0, length(coef(fit2)))
  names(L2) <- names(coef(fit2))
  if ("grade3" %in% names(L2) && "grade2" %in% names(L2)) {
    L2["grade3"] <- 1; L2["grade2"] <- -1
  }
  est2 <- sum(L2 * coef(fit2))
  se2 <- sqrt(L2 %*% vcov(fit2) %*% L2)
  p_no_pmi <- 2 * pnorm(-abs(est2 / se2))

  data.frame(
    logFC_ctrl_PMI = round(est, 4),
    p_ctrl_PMI = signif(p_ctrl_pmi, 3),
    logFC_no_PMI = round(est2, 4),
    p_no_PMI = signif(p_no_pmi, 3),
    stringsAsFactors = FALSE)
}

mt_res <- data.frame()
for (nm in names(s_all)) {
  if (is.null(s_all[[nm]])) next
  r <- run_one(s_all[[nm]], clin, common_samp)
  mt_res <- rbind(mt_res, data.frame(geneset = nm, r, stringsAsFactors = FALSE))
}
save_table(mt_res, DIR_THIS, "mtdna_vs_nuc_g3v2.csv")
cat("\n=== g3_vs_g2 效应：控/不控 PMI ===\n")
print(mt_res, row.names = FALSE)

# ============================================================
# ============================================================
log_info("--- 6. PMI 稳健通路清单 ---")

robust <- all_models[all_models$model == "M4_age_PMI_spline_RIN" &
                     all_models$contrast == "g3_vs_g2" &
                     all_models$FDR < 0.05 &
                     abs(all_models$logFC) > 0.05, ]
robust <- merge(robust, pmi_cor[, c("pathway","r_pmi")], by = "pathway")
robust <- robust[abs(robust$r_pmi) < 0.3, ]
robust <- robust[order(robust$FDR), ]
save_table(robust, DIR_THIS, "pmi_robust_pathways.csv")

cat("\n=== PMI 稳健通路（M4 FDR<0.05, |logFC|>0.05, |r_PMI|<0.3）===\n")
if (nrow(robust) > 0) {
  print(robust[, c("pathway","group","logFC","p","FDR","r_pmi")],
        row.names = FALSE, digits = 4)
} else {
  cat("  无——所有 g3_vs_g2 显著通路都与 PMI 相关或控 PMI 后消失。\n")
}

# ============================================================
# 7. figure
# ============================================================
log_info("--- 7. 图 ---")

clin$grade_f <- factor(clin$grade)
p1 <- ggplot(clin, aes(x = grade_f, y = pmi_v, fill = grade_f)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.7) +
  geom_jitter(width = 0.15, size = 0.5, alpha = 0.4) +
  theme_bw(base_size = 11) +
  labs(x = "AMD grade", y = "PMI",
       title = sprintf("PMI across grades (KW p = %.3g)",
                       kruskal.test(clin$pmi_v ~ clin$grade)$p.value)) +
  theme(legend.position = "none")
save_fig(p1, DIR_THIS, "pmi_grade_boxplot", width = 7, height = 5)

if ("group" %in% colnames(pmi_cor)) {
  top_n <- 40
  top_pmi <- head(pmi_cor[order(-abs(pmi_cor$r_pmi)), ], top_n)
  hm_mat <- as.matrix(top_pmi[, c("r_pmi","r_age","r_rin")])
  rownames(hm_mat) <- top_pmi$pathway
  hm_mat[is.na(hm_mat)] <- 0
  pheatmap(hm_mat,
           cluster_rows = TRUE, cluster_cols = FALSE,
           display_numbers = TRUE, number_format = "%.2f",
           color = colorRampPalette(c("#4C72B0","#FFFFFF","#C44E52"))(50),
           main = sprintf("Top %d pathways by |rho(PMI)|", top_n),
           filename = file.path(DIR_THIS, "figures",
                                "pmi_pathway_correlation.png"),
           width = 6, height = 9)
}

curve_mt <- data.frame()
for (nm in names(s_all)) {
  if (is.null(s_all[[nm]])) next
  for (g in 1:4) {
    idx <- clin[common_samp, "grade"] == g
    x <- s_all[[nm]][idx]
    curve_mt <- rbind(curve_mt, data.frame(
      geneset = nm, grade = g,
      mean_score = mean(x, na.rm = TRUE),
      se = sd(x, na.rm = TRUE) / sqrt(sum(idx)),
      stringsAsFactors = FALSE))
  }
}
p3 <- ggplot(curve_mt, aes(x = grade, y = mean_score,
                             color = geneset, group = geneset)) +
  geom_line(linewidth = 1) + geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean_score - se, ymax = mean_score + se),
                 width = 0.15, alpha = 0.5) +
  theme_bw(base_size = 11) +
  scale_x_continuous(breaks = 1:4, labels = paste0("G", 1:4)) +
  labs(x = "AMD grade", y = "z-score (mean ± SE)",
       title = "mtDNA-encoded vs nuclear-mito vs control",
       color = "")
save_fig(p3, DIR_THIS, "mtdna_vs_nuc_curves", width = 8, height = 5)

# ============================================================
# 8. summary
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("=== diag_pmi_mito.r 汇总 ===\n")
cat(strrep("=", 70), "\n", sep = "")

cat("PMI 在 grade 间 KW p:", signif(
  kruskal.test(clin$pmi_v ~ clin$grade)$p.value, 3), "\n")
cat("Spearman rho(PMI, grade):",
    round(cor(clin$pmi_v, clin$grade, method = "spearman",
              use = "complete.obs"), 3), "\n\n")

cat("g3_vs_g2 控 PMI 后 FDR<0.05 数量:\n")
for (m in unique(all_models$model)) {
  n <- sum(all_models$FDR[all_models$model == m &
                          all_models$contrast == "g3_vs_g2"] < 0.05,
           na.rm = TRUE)
  cat(sprintf("  %-30s  %d\n", m, n))
}

cat("\nPMI 稳健 pathway:", nrow(robust), "\n")
cat("PMI 驱动 pathway:", sum(cmp$pmi_driven, na.rm = TRUE), "\n")

cat("\n--- 输出位置 ---\n")
cat("目录:", DIR_THIS, "\n")

log_info("=== diag_pmi_mito.r done ===")
