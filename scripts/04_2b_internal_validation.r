#!/usr/bin/env Rscript
# ============================================================
# ============================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2)
})

DIR_ROOT   <- "E:/retina_tools/aging_rgc_axon/AMD"
DIR_RESULT <- file.path(DIR_ROOT, "results")
SEED       <- 20260925

log_info <- function(msg) cat(sprintf("[%s] [INFO] %s
", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
log_warn <- function(msg) cat(sprintf("[%s] [WARN] %s
", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))

read_upstream <- function(module, subdir, filename) {
  path <- file.path(DIR_RESULT, module, subdir, filename)
  if (!file.exists(path)) stop("missing: ", path)
  readRDS(path)
}
save_table <- function(df, dir, name) {
  write.csv(df, file.path(dir, "tables", name), row.names = FALSE)
}
save_fig <- function(p, dir, name, width = 6, height = 4) {
  ggsave(file.path(dir, "figures", paste0(name, ".pdf")), p, width = width, height = height)
  ggsave(file.path(dir, "figures", paste0(name, ".png")), p, width = width, height = height, dpi = 300)
}
setup_output_dir <- function(name, subdir = NULL) {
  d <- if (is.null(subdir)) file.path(DIR_RESULT, name) else file.path(DIR_RESULT, subdir, name)
  for (sub in c("data", "tables", "figures"))
    dir.create(file.path(d, sub), recursive = TRUE, showWarnings = FALSE)
  d
}

DIR_THIS <- setup_output_dir("04_2_pathway_score")
log_info("=== 04_2b v3 start ===")

N_BOOT <- 1000
N_RAND <- 1000

score_mat <- readRDS(file.path(DIR_THIS, "data", "pathway_scores_ssgsea.rds"))
clin <- read_upstream("01_1b_batch_enhance", "data", "clinical_wgcna.rds")
expr <- read_upstream("01_1b_batch_enhance", "data", "expr_wgcna.rds")
anno <- read_upstream("01_1_prepare_bulk_data", "data", "gene_annotation.rds")
GENE_SETS <- read_upstream("01_2_prepare_gene_sets", "data", "gene_sets.rds")

# Ensembl -> Symbol
anno_sub <- anno[anno$ensembl %in% rownames(expr), ]
anno_sub <- anno_sub[!is.na(anno_sub$symbol), ]
expr <- expr[anno_sub$ensembl, , drop = FALSE]
expr <- expr[!duplicated(anno_sub$symbol), , drop = FALSE]
rownames(expr) <- toupper(anno_sub$symbol[!duplicated(anno_sub$symbol)])
log_info(sprintf("expr symbol: %d x %d", nrow(expr), ncol(expr)))

common <- intersect(colnames(score_mat), rownames(clin))
common <- intersect(common, colnames(expr))
score_mat <- score_mat[, common, drop = FALSE]
expr <- expr[, common, drop = FALSE]
grade <- as.integer(clin[common, "grade"])
log_info(sprintf("样本: %d | grade: %s", length(grade),
                 paste(names(table(grade)), table(grade), sep = "=", collapse = " ")))

target_pathways <- c("Mito_CI_Subunits","Mito_CIV_Subunits","Mito_CV_Subunits",
                     "OXPHOS","Complement","Complement_Activation_GOBP",
                     "Fatty_acid_oxidation")
target_pathways <- intersect(target_pathways, rownames(score_mat))
log_info(sprintf("目标通路: %s", paste(target_pathways, collapse = ", ")))

fit_poly3_p3 <- function(y, g) {
  df <- data.frame(y = y, g = g)
  df <- df[complete.cases(df), ]
  if (nrow(df) < 20 || length(unique(df$g)) < 4) return(NA_real_)
  summary(lm(y ~ poly(g, 3), data = df))$coefficients[4, 4]
}

# ---- Bootstrap ----
set.seed(SEED)
boot_res <- data.frame()
for (pw in target_pathways) {
  y <- score_mat[pw, ]
  p3_vals <- sapply(1:N_BOOT, function(i) {
    idx <- sample(length(y), floor(0.9 * length(y)))
    fit_poly3_p3(y[idx], grade[idx])
  })
  obs_p3 <- fit_poly3_p3(y, grade)
  boot_res <- rbind(boot_res, data.frame(
    pathway = pw, obs_p3 = round(obs_p3, 5),
    boot_median = round(median(p3_vals, na.rm = TRUE), 5),
    boot_q025 = round(quantile(p3_vals, 0.025, na.rm = TRUE), 5),
    boot_q975 = round(quantile(p3_vals, 0.975, na.rm = TRUE), 5),
    pct_sig = round(mean(p3_vals < 0.05, na.rm = TRUE), 3),
    stringsAsFactors = FALSE))
  log_info(sprintf("%-30s obs_p3=%.4f pct_sig=%.3f", pw, obs_p3,
                   mean(p3_vals < 0.05, na.rm = TRUE)))
}
save_table(boot_res, DIR_THIS, "bootstrap_p3_results.csv")
print(boot_res, row.names = FALSE)

set.seed(SEED)
expr_z <- t(scale(t(expr)))
expr_z[is.na(expr_z)] <- 0

random_p3 <- data.frame()
for (pw in target_pathways) {
  genes_real <- toupper(GENE_SETS[[pw]]$genes)
  genes_real <- intersect(genes_real, rownames(expr))
  n_genes <- length(genes_real)
  log_info(sprintf("%-30s GENE_SETS n=%d matched n=%d", pw,
                   length(GENE_SETS[[pw]]$genes), n_genes))
  if (n_genes < 10) { log_warn(sprintf("%s skip", pw)); next }

  y_real <- colMeans(expr_z[genes_real, common, drop = FALSE])
  obs_p3 <- fit_poly3_p3(y_real, grade)
  p3_null <- sapply(1:N_RAND, function(i) {
    genes_rand <- sample(rownames(expr), n_genes)
    y_rand <- colMeans(expr_z[genes_rand, common, drop = FALSE])
    fit_poly3_p3(y_rand, grade)
  })
  random_p3 <- rbind(random_p3, data.frame(
    pathway = pw, n_genes = n_genes, obs_p3 = round(obs_p3, 5),
    null_median = round(median(p3_null, na.rm = TRUE), 5),
    null_q05 = round(quantile(p3_null, 0.05, na.rm = TRUE), 5),
    empirical_p = round(mean(p3_null <= obs_p3, na.rm = TRUE), 4),
    stringsAsFactors = FALSE))
  log_info(sprintf("%-30s obs_p3=%.4f empirical_p=%.4f", pw, obs_p3,
                   mean(p3_null <= obs_p3, na.rm = TRUE)))
}
save_table(random_p3, DIR_THIS, "random_geneset_p3.csv")
print(random_p3, row.names = FALSE)

log_info("=== 04_2b v3 done ===")

