#!/usr/bin/env Rscript
# ============================================================
# diag_oxphos_pca.r
#
#
# Input：
#   results/01_1_prepare_bulk_data/data/{expr_corrected,clinical,gene_annotation}.rds
#   results/01_2_prepare_gene_sets/data/gene_sets.rds
#
#   figures/pc1_variance_bar.pdf/.png
#   figures/pc1_grade_scatter.pdf/.png
#   figures/pc1_grade_boxplot.pdf/.png
#   figures/top_loading_genes.pdf/.png
#   figures/random_null.pdf/.png
# ============================================================

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(ggplot2)
  library(pheatmap)
  library(tidyr)
})

DIR_THIS <- setup_output_dir("diag_oxphos_pca")
log_info("=== diag_oxphos_pca.r start ===")

# ============================================================
# ============================================================
log_info("--- 1. 载入数据 ---")

expr <- read_upstream("01_1_prepare_bulk_data", "data", "expr_corrected.rds")
clin <- read_upstream("01_1_prepare_bulk_data", "data", "clinical.rds")
anno <- read_upstream("01_1_prepare_bulk_data", "data", "gene_annotation.rds")
GENE_SETS <- read_upstream("01_2_prepare_gene_sets", "data", "gene_sets.rds")

# optional：hUSI
husi_file <- file.path(DIR_RESULT, "04_1_aging_score", "data", "aging_scores.rds")
aging_scores <- if (file.exists(husi_file)) readRDS(husi_file) else NULL

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
rownames(clin) <- colnames(expr)

log_info("expr: {nrow(expr)} × {ncol(expr)}")
log_info("gene sets: {length(GENE_SETS)}")

# ============================================================
# ============================================================
TARGET_GS <- c(
  "Mito_CI_Subunits", "Mito_CII_Subunits", "Mito_CIII_Subunits",
  "Mito_CIV_Subunits", "Mito_CV_Subunits",
  "Mito_Complex_I", "Mito_Complex_II", "Mito_Complex_III",
  "Mito_Complex_IV", "Mito_Complex_V",
  "Mito_CI_Assembly", "Mito_CIV_Assembly", "Mito_CV_Assembly",
  "Mito_TCA", "Mito_TCA_associated",
  "Mito_Translation", "Mito_Ribosome", "Mito_Proteostasis",
  "Mito_Import", "Mito_Carriers", "Mito_FattyAcidOxidation",
  "Mito_ROS_Glutathione", "Mito_Respirasome", "Mito_Cristae",
  "OXPHOS", "Oxidative_stress"
)
TARGET_GS <- intersect(TARGET_GS, names(GENE_SETS))
log_info("待检查基因集: {length(TARGET_GS)} 个")

# ============================================================
# ============================================================

run_gs_pca <- function(genes, expr, clin, n_pc = 5) {
  genes_in <- intersect(genes, rownames(expr))
  if (length(genes_in) < 3) return(NULL)
  sub <- expr[genes_in, , drop = FALSE]
  sub_z <- t(scale(t(sub)))
  sub_z[is.na(sub_z)] <- 0
  pc <- prcomp(t(sub_z), center = FALSE, scale. = FALSE)
  var_explained <- (pc$sdev^2) / sum(pc$sdev^2)
  n_use <- min(n_pc, length(var_explained))
  list(
    genes = genes_in,
    pca = pc,
    var_explained = var_explained[1:n_use],
    scores = pc$x[, 1:n_use, drop = FALSE],
    loadings = pc$rotation[, 1:n_use, drop = FALSE],
    n_genes = length(genes_in)
  )
}

cor_pc_covar <- function(pca_res, clin, n_pc = 5) {
  out <- data.frame()
  covars <- list(
    grade = clin$grade,
    age = clin$age_v,
    pmi = clin$pmi_v,
    rin = clin$rin_v,
    sex_M = as.integer(clin$sex == "Male")
  )
  for (pc_name in colnames(pca_res$scores)) {
    pc_v <- pca_res$scores[, pc_name]
    for (cv_name in names(covars)) {
      cv <- covars[[cv_name]]
      if (all(is.na(cv))) next
      sp_r <- suppressWarnings(cor(pc_v, cv, method = "spearman",
                                    use = "complete.obs"))
      out <- rbind(out, data.frame(
        PC = pc_name, covar = cv_name,
        spearman_r = round(sp_r, 3),
        stringsAsFactors = FALSE))
    }
    if (!is.null(aging_scores) &&
        "hUSI" %in% rownames(aging_scores)) {
      common_samp <- intersect(colnames(aging_scores), names(pc_v))
      if (length(common_samp) > 30) {
        hUSI_v <- aging_scores["hUSI", common_samp]
        r <- cor(pc_v[common_samp], hUSI_v, method = "spearman",
                 use = "complete.obs")
        out <- rbind(out, data.frame(
          PC = pc_name, covar = "hUSI",
          spearman_r = round(r, 3),
          stringsAsFactors = FALSE))
      }
    }
  }
  out
}

analyze_gs <- function(gs_name, GENE_SETS, expr, clin) {
  gs <- GENE_SETS[[gs_name]]
  if (is.null(gs)) return(NULL)
  pca_res <- run_gs_pca(gs$genes, expr, clin)
  if (is.null(pca_res)) return(NULL)

  # summary
  ve <- pca_res$var_explained
  summary_row <- data.frame(
    geneset = gs_name,
    group = if (!is.null(gs$group)) gs$group else NA,
    n_genes = pca_res$n_genes,
    PC1_var = round(ve[1], 3),
    PC2_var = if (length(ve) >= 2) round(ve[2], 3) else NA,
    PC3_var = if (length(ve) >= 3) round(ve[3], 3) else NA,
    PC1_PC2_cum = if (length(ve) >= 2) round(ve[1] + ve[2], 3) else NA,
    PC1_PC3_cum = if (length(ve) >= 3) round(sum(ve[1:3]), 3) else NA,
    stringsAsFactors = FALSE)

  # PC vs covar
  pc_cor <- cor_pc_covar(pca_res, clin)
  pc_cor$geneset <- gs_name

  # gene loadings PC1
  load_df <- data.frame(
    geneset = gs_name,
    gene = rownames(pca_res$loadings),
    PC1_loading = round(pca_res$loadings[, 1], 4),
    stringsAsFactors = FALSE)
  load_df$abs_loading <- abs(load_df$PC1_loading)
  load_df <- load_df[order(-load_df$abs_loading), ]

  list(summary = summary_row,
       pc_cor = pc_cor,
       loadings = load_df,
       pca = pca_res)
}

# ============================================================
# ============================================================
log_info("--- 4. 跑所有目标基因集 PCA ---")

all_res <- lapply(TARGET_GS, function(gs) {
  tryCatch(analyze_gs(gs, GENE_SETS, expr, clin),
           error = function(e) { log_warn("{gs}: {e$message}"); NULL })
})
names(all_res) <- TARGET_GS
all_res <- all_res[!sapply(all_res, is.null)]

# summary
gs_summary <- do.call(rbind, lapply(all_res, function(x) x$summary))
gs_summary <- gs_summary[order(-gs_summary$PC1_var), ]
save_table(gs_summary, DIR_THIS, "geneset_summary.csv")

cat("\n=== 基因集 PC1 解释方差 ===\n")
print(gs_summary, row.names = FALSE)

# PC vs covar
pc_cor_all <- do.call(rbind, lapply(all_res, function(x) x$pc_cor))
save_table(pc_cor_all, DIR_THIS, "pc_covar_correlation.csv")

pc1_grade <- pc_cor_all[pc_cor_all$PC == "PC1" &
                          pc_cor_all$covar == "grade", ]
pc1_grade <- pc1_grade[order(-abs(pc1_grade$spearman_r)), ]
cat("\n=== PC1 与 grade 的相关（按 |r| 排序）===\n")
print(pc1_grade, row.names = FALSE)

# gene loadings
all_loadings <- do.call(rbind, lapply(all_res, function(x) x$loadings))
save_table(all_loadings, DIR_THIS, "gene_loadings.csv")

# ============================================================
# ============================================================
log_info("--- 5. 随机对照 null 分布 ---")

set.seed(SEED)
N_RANDOM <- 200

random_null <- data.frame()
for (gs_name in names(all_res)) {
  n_g <- all_res[[gs_name]]$summary$n_genes
  real_pc1_var <- all_res[[gs_name]]$summary$PC1_var
  real_pc1_grade <- pc1_grade$spearman_r[pc1_grade$geneset == gs_name]
  if (length(real_pc1_grade) == 0) real_pc1_grade <- NA

  mito_genes <- unique(unlist(lapply(GENE_SETS[
    grepl("^Mito_", names(GENE_SETS))], function(x) x$genes)))
  pool <- setdiff(rownames(expr), mito_genes)

  null_pc1_var <- numeric(N_RANDOM)
  null_pc1_grade <- numeric(N_RANDOM)
  for (i in seq_len(N_RANDOM)) {
    rnd_genes <- sample(pool, n_g)
    pca_r <- run_gs_pca(rnd_genes, expr, clin)
    if (is.null(pca_r)) next
    null_pc1_var[i] <- pca_r$var_explained[1]
    r <- cor(pca_r$scores[, 1], clin$grade, method = "spearman",
             use = "complete.obs")
    null_pc1_grade[i] <- r
  }

  random_null <- rbind(random_null, data.frame(
    geneset = gs_name,
    n_genes = n_g,
    real_PC1_var = real_pc1_var,
    null_PC1_var_mean = round(mean(null_pc1_var, na.rm = TRUE), 3),
    null_PC1_var_sd = round(sd(null_pc1_var, na.rm = TRUE), 3),
    null_PC1_var_p = round(mean(null_pc1_var >= real_pc1_var, na.rm = TRUE), 4),
    real_PC1_grade_r = round(real_pc1_grade, 3),
    null_PC1_grade_abs_mean = round(mean(abs(null_pc1_grade), na.rm = TRUE), 3),
    null_PC1_grade_abs_p = round(mean(abs(null_pc1_grade) >= abs(real_pc1_grade),
                                       na.rm = TRUE), 4),
    stringsAsFactors = FALSE))
}
save_table(random_null, DIR_THIS, "random_null.csv")

cat("\n=== 随机对照 null ===\n")
print(random_null, row.names = FALSE)

# ============================================================
# ============================================================
log_info("--- 6. 核编码 vs mtDNA 分层 ---")

MTDNA_GENES <- c("MT-ND1","MT-ND2","MT-ND3","MT-ND4","MT-ND4L",
                 "MT-ND5","MT-ND6","MT-CYB","MT-CO1","MT-CO2",
                 "MT-CO3","MT-ATP6","MT-ATP8")

COMPLEX_LIST <- c("Mito_Complex_I", "Mito_Complex_III",
                   "Mito_Complex_IV", "Mito_Complex_V")

nuc_vs_mt <- data.frame()
for (cp in COMPLEX_LIST) {
  if (!cp %in% names(GENE_SETS)) next
  genes <- GENE_SETS[[cp]]$genes
  genes_in <- intersect(genes, rownames(expr))
  nuc_genes <- setdiff(genes_in, MTDNA_GENES)
  mt_genes  <- intersect(genes_in, MTDNA_GENES)

  for (sub_name in c("nuclear", "mtDNA", "all")) {
    sub_genes <- switch(sub_name,
      nuclear = nuc_genes,
      mtDNA   = mt_genes,
      all     = genes_in)
    if (length(sub_genes) < 3) next
    pca_res <- run_gs_pca(sub_genes, expr, clin)
    if (is.null(pca_res)) next
    r_grade <- cor(pca_res$scores[, 1], clin$grade,
                   method = "spearman", use = "complete.obs")
    r_pmi <- cor(pca_res$scores[, 1], clin$pmi_v,
                 method = "spearman", use = "complete.obs")
    r_rin <- cor(pca_res$scores[, 1], clin$rin_v,
                 method = "spearman", use = "complete.obs")
    nuc_vs_mt <- rbind(nuc_vs_mt, data.frame(
      complex = cp, subset = sub_name,
      n_genes = length(sub_genes),
      PC1_var = round(pca_res$var_explained[1], 3),
      PC1_grade_r = round(r_grade, 3),
      PC1_pmi_r = round(r_pmi, 3),
      PC1_rin_r = round(r_rin, 3),
      stringsAsFactors = FALSE))
  }
}
save_table(nuc_vs_mt, DIR_THIS, "nuc_vs_mtdna.csv")

cat("\n=== 核编码 vs mtDNA 分层 ===\n")
print(nuc_vs_mt, row.names = FALSE)

# ============================================================
# ============================================================
log_info("--- 7. 跨基因集 PC1 一致性 ---")

PC1_MAT_GS <- intersect(c("Mito_Complex_I", "Mito_Complex_III",
                            "Mito_Complex_IV", "Mito_Complex_V",
                            "Mito_CI_Subunits", "Mito_CIV_Subunits",
                            "Mito_CV_Subunits", "Mito_CIII_Subunits"),
                          names(all_res))
if (length(PC1_MAT_GS) >= 3) {
  pc1_mat <- do.call(cbind, lapply(PC1_MAT_GS, function(gs) {
    pc1 <- all_res[[gs]]$pca$scores[, 1]
    r <- cor(pc1, clin$grade, method = "spearman", use = "complete.obs")
    if (r > 0) pc1 <- -pc1
    pc1
  }))
  colnames(pc1_mat) <- PC1_MAT_GS
  # correlationmatrix
  pc1_cor <- cor(pc1_mat, method = "spearman")
  pc1_cor_long <- data.frame()
  for (i in seq_len(ncol(pc1_cor))) {
    for (j in seq_len(ncol(pc1_cor))) {
      if (i >= j) next
      pc1_cor_long <- rbind(pc1_cor_long, data.frame(
        gs1 = colnames(pc1_cor)[i],
        gs2 = colnames(pc1_cor)[j],
        spearman_r = round(pc1_cor[i, j], 3),
        stringsAsFactors = FALSE))
    }
  }
  save_table(pc1_cor_long, DIR_THIS, "pc_concordance.csv")
  cat("\n=== 跨基因集 PC1 一致性 ===\n")
  print(pc1_cor_long, row.names = FALSE)
}

# ============================================================
# 8. figure
# ============================================================
log_info("--- 8. 图 ---")

gs_summary_plot <- gs_summary
gs_summary_plot$null_mean <- random_null$null_PC1_var_mean[
  match(gs_summary_plot$geneset, random_null$geneset)]
gs_summary_plot$geneset_f <- factor(gs_summary_plot$geneset,
  levels = gs_summary_plot$geneset[order(gs_summary_plot$PC1_var)])

p1 <- ggplot(gs_summary_plot, aes(x = geneset_f, y = PC1_var)) +
  geom_col(fill = "#4C72B0", alpha = 0.7) +
  geom_point(aes(y = null_mean), color = "red", size = 1.5) +
  geom_hline(yintercept = 0.4, linetype = "dashed", color = "gray40") +
  coord_flip() + theme_bw(base_size = 10) +
  labs(x = "", y = "PC1 variance explained",
       title = "PC1 variance: real (blue) vs random null mean (red)",
       subtitle = "Dashed line = 0.4 threshold for coordinated module")
save_fig(p1, DIR_THIS, "pc1_variance_bar", width = 9, height = 7)

pc1_df <- data.frame()
for (gs in names(all_res)) {
  pc1_v <- all_res[[gs]]$pca$scores[, 1]
  r <- cor(pc1_v, clin$grade, method = "spearman", use = "complete.obs")
  if (r > 0) pc1_v <- -pc1_v
  pc1_df <- rbind(pc1_df, data.frame(
    geneset = gs,
    sample = names(pc1_v),
    PC1 = pc1_v,
    grade = clin[names(pc1_v), "grade"],
    stringsAsFactors = FALSE))
}
focus_gs <- intersect(c("Mito_CI_Subunits", "Mito_CIV_Subunits",
                          "Mito_CV_Subunits", "Mito_CIII_Subunits",
                          "Mito_Complex_I", "Mito_Complex_IV",
                          "Mito_Complex_V", "OXPHOS",
                          "Oxidative_stress", "Mito_TCA",
                          "Mito_Translation", "Mito_Proteostasis"),
                        unique(pc1_df$geneset))
pc1_df_focus <- pc1_df[pc1_df$geneset %in% focus_gs, ]

p2 <- ggplot(pc1_df_focus, aes(x = factor(grade), y = PC1)) +
  geom_boxplot(aes(fill = factor(grade)), alpha = 0.6,
                outlier.size = 0.5) +
  geom_jitter(width = 0.15, size = 0.4, alpha = 0.3) +
  facet_wrap(~ geneset, scales = "free_y", ncol = 3) +
  theme_bw(base_size = 10) +
  labs(x = "AMD grade", y = "PC1 score (oriented so ↓ with grade)",
       title = "PC1 vs grade across key pathways") +
  theme(legend.position = "none")
save_fig(p2, DIR_THIS, "pc1_grade_boxplot", width = 12, height = 10)

# 8c. Top loading genes
loadings_top <- do.call(rbind, lapply(focus_gs, function(gs) {
  l <- all_loadings[all_loadings$geneset == gs, ]
  head(l, 10)
}))
loadings_top$geneset_f <- factor(loadings_top$geneset, levels = focus_gs)

p3 <- ggplot(loadings_top, aes(x = reorder(gene, PC1_loading),
                                 y = PC1_loading)) +
  geom_col(aes(fill = PC1_loading > 0), alpha = 0.8) +
  facet_wrap(~ geneset_f, scales = "free", ncol = 3) +
  coord_flip() + theme_bw(base_size = 9) +
  scale_fill_manual(values = c("TRUE" = "#C44E52", "FALSE" = "#4C72B0")) +
  labs(x = "", y = "PC1 loading",
       title = "Top 10 PC1 loading genes per pathway") +
  theme(legend.position = "none")
save_fig(p3, DIR_THIS, "top_loading_genes", width = 12, height = 12)

p4 <- ggplot(random_null, aes(x = reorder(geneset, real_PC1_var),
                                y = real_PC1_var)) +
  geom_point(color = "red", size = 2) +
  geom_errorbar(aes(ymin = null_PC1_var_mean - null_PC1_var_sd,
                     ymax = null_PC1_var_mean + null_PC1_var_sd),
                 color = "gray40", width = 0.3) +
  geom_point(aes(y = null_PC1_var_mean), color = "black", size = 1.5) +
  coord_flip() + theme_bw(base_size = 10) +
  labs(x = "", y = "PC1 variance explained",
       title = "Real (red) vs null (black, mean ± SD)",
       subtitle = "If red > black + 2SD: coordinated module")
save_fig(p4, DIR_THIS, "random_null", width = 9, height = 7)

# ============================================================
# 9. summary
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("=== diag_oxphos_pca.r 汇总 ===\n")
cat(strrep("=", 70), "\n", sep = "")

cat("\n=== 关键判据 ===\n")
key_gs <- c("Mito_CI_Subunits", "Mito_CIII_Subunits",
            "Mito_CIV_Subunits", "Mito_CV_Subunits",
            "Mito_Complex_I", "Mito_Complex_IV", "Mito_Complex_V",
            "OXPHOS", "Oxidative_stress")
for (gs in key_gs) {
  if (!gs %in% names(all_res)) next
  pc1_var <- all_res[[gs]]$summary$PC1_var
  null_var <- random_null$null_PC1_var_mean[random_null$geneset == gs]
  null_sd <- random_null$null_PC1_var_sd[random_null$geneset == gs]
  r_grade <- pc1_grade$spearman_r[pc1_grade$geneset == gs]
  r_hUSI <- pc_cor_all$spearman_r[pc_cor_all$geneset == gs &
                                    pc_cor_all$PC == "PC1" &
                                    pc_cor_all$covar == "hUSI"]
  if (length(r_hUSI) == 0) r_hUSI <- NA
  cat(sprintf("  %-25s PC1_var=%.3f  null=%.3f±%.3f  r_grade=%+.3f  r_hUSI=%+.3f\n",
              gs, pc1_var, null_var, null_sd,
              if (length(r_grade) > 0) r_grade else NA,
              r_hUSI))
}

cat("\n--- 输出位置 ---\n")
cat("目录:", DIR_THIS, "\n")
cat("  tables/geneset_summary.csv\n")
cat("  tables/pc_covar_correlation.csv\n")
cat("  tables/gene_loadings.csv\n")
cat("  tables/random_null.csv\n")
cat("  tables/nuc_vs_mtdna.csv\n")
cat("  tables/pc_concordance.csv\n")

log_info("=== diag_oxphos_pca.r done ===")
