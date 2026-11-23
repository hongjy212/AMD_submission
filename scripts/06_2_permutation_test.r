#!/usr/bin/env Rscript
# ============================================================
# 06_2_permutation_test.r
#
#
#   tables/permutation_summary.csv
#   tables/permutation_dist_16genes.csv
#   tables/permutation_dist_modules.csv
#   figures/permutation_histogram.pdf/.png
#   rds_objects/permutation_results.rds
# ============================================================

if (exists(".AMD_SETUP_LOADED", envir = .GlobalEnv)) {
  rm(".AMD_SETUP_LOADED", envir = .GlobalEnv)
}
source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")

suppressPackageStartupMessages({
  library(limma); library(splines); library(ggplot2)
})

DIR_THIS <- setup_output_dir("06_2_permutation_test")
log_info("=== 06_2_permutation_test.r start ===")

N_PERM <- 1000   # 先跑 100 试探；确认后再改 1000

# ============================================================
# ============================================================
log_info("--- 1. 载入 ---")

expr <- read_upstream("01_1b_batch_enhance", "data", "expr_wgcna.rds")
clin <- read_upstream("01_1b_batch_enhance", "data", "clinical_wgcna.rds")
anno <- read_upstream("01_1_prepare_bulk_data", "data", "gene_annotation.rds")
module_df <- read_upstream("05_1_consensus_wgcna", "tables",
                            "consensus_module_genes_annotated.csv")

common <- intersect(colnames(expr), rownames(clin))
expr <- expr[, common, drop = FALSE]
clin <- clin[common, , drop = FALSE]

clin$grade <- as.integer(clin$grade)
clin$age_v <- as.numeric(clin$age_n)
clin$rin_v <- as.numeric(clin$rin_n)
clin$pmi_v <- as.numeric(clin$pmi_n)
clin$batch_f <- factor(clin$batch)
sid <- as.character(clin$sample_id)
clin$layer_f <- factor(sub("^.*_(L[0-9]+)$", "\\1", sid),
                        levels = paste0("L", 1:8))
clin$sex_f <- factor(clin$sex)
rownames(clin) <- colnames(expr)

log_info("expr: {nrow(expr)} × {ncol(expr)}")
log_info("grade: {paste(names(table(clin$grade)), table(clin$grade), sep='=', collapse=' ')}")

# ============================================================
# ============================================================
log_info("--- 2. 基因集准备 ---")

strip_version <- function(x) sub("\\..*$", "", x)
sym_map <- setNames(toupper(anno$symbol), strip_version(anno$ensembl))
expr_row_sym <- sym_map[strip_version(rownames(expr))]

# 16 genes（02_1b）
GENES_16 <- c("FCF1","UGGT1","CCDC3","S100B","BPNT1","BRI3BP","BTBD19",
              "DPH7","AKR1B1","ZNF623","ZNF587","RSPH3","TAF3","NDUFB8",
              "ZNF417","C6orf120")
genes_16_in <- rownames(expr)[expr_row_sym %in% GENES_16]
log_info("16 基因命中: {length(genes_16_in)}")

# WGCNA module（05_3）
mod_names <- c("greenyellow", "black", "green")
mod_genes <- list()
for (mod in mod_names) {
  mod_ens <- strip_version(module_df$ENSEMBL[module_df$Module == mod])
  mod_in_expr <- intersect(mod_ens, strip_version(rownames(expr)))
  mod_genes[[mod]] <- rownames(expr)[strip_version(rownames(expr)) %in% mod_in_expr]
  log_info("模块 {mod}: {length(mod_genes[[mod]])} 基因")
}

# ============================================================
# ============================================================
run_poly3F <- function(expr_sub, grade_vec, clin_sub) {
  poly_mat <- poly(grade_vec, 3)
  colnames(poly_mat) <- c("p1", "p2", "p3")
  rownames(poly_mat) <- rownames(clin_sub)

  sp_age <- splines::ns(clin_sub$age_v, df = 3)
  design_df <- data.frame(
    poly_mat,
    age_s1 = sp_age[, 1],
    age_s2 = sp_age[, 2],
    age_s3 = sp_age[, 3],
    sex_M  = as.integer(clin_sub$sex_f == "Male"),
    rin    = clin_sub$rin_v,
    pmi    = clin_sub$pmi_v,
    batch_f = clin_sub$batch_f,
    layer_f = clin_sub$layer_f,
    stringsAsFactors = FALSE)
  rownames(design_df) <- rownames(clin_sub)

  design <- model.matrix(
    ~ p1 + p2 + p3 + age_s1 + age_s2 + age_s3 +
      sex_M + rin + pmi + batch_f + layer_f,
    data = design_df)

  fit <- eBayes(lmFit(expr_sub, design))
  tt_F  <- topTable(fit, coef = c("p1","p2","p3"),
                     number = Inf, sort.by = "none")
  tt_p3 <- topTable(fit, coef = "p3", number = Inf, sort.by = "none")

  data.frame(
    gene = rownames(tt_F),
    F_stat = tt_F$F,
    F_FDR  = tt_F$adj.P.Val,
    p3_logFC = tt_p3$logFC,
    p3_FDR = tt_p3$adj.P.Val,
    stringsAsFactors = FALSE)
}

# ============================================================
# ============================================================
log_info("--- 4. 真实值 ---")

expr_16 <- expr[genes_16_in, , drop = FALSE]
res_16_real <- run_poly3F(expr_16, clin$grade, clin)
n_p3_real_16 <- sum(res_16_real$p3_FDR < 0.05, na.rm = TRUE)
n_F_real_16  <- sum(res_16_real$F_FDR < 0.05, na.rm = TRUE)

log_info("16 基因真实: p3_FDR<0.05 = {n_p3_real_16}, F_FDR<0.05 = {n_F_real_16}")

res_mod_real <- list()
n_p3_real_mod <- numeric(length(mod_names))
for (i in seq_along(mod_names)) {
  mod <- mod_names[i]
  expr_mod <- expr[mod_genes[[mod]], , drop = FALSE]
  res <- run_poly3F(expr_mod, clin$grade, clin)
  res_mod_real[[mod]] <- res
  n_p3_real_mod[i] <- sum(res$p3_FDR < 0.05, na.rm = TRUE)
  log_info("模块 {mod} 真实 p3_FDR<0.05: {n_p3_real_mod[i]} / {length(mod_genes[[mod]])}")
}

# ============================================================
# 5. Permutation
# ============================================================
log_info("--- 5. Permutation (N={N_PERM}) ---")
log_info("每次 permutation: 打乱 grade → 重跑 poly3F → 数 p3_FDR<0.05")

set.seed(SEED)
n_perm_16 <- numeric(N_PERM)
n_perm_mod <- matrix(NA_integer_, N_PERM, length(mod_names))
colnames(n_perm_mod) <- mod_names

t0 <- Sys.time()
for (i in 1:N_PERM) {
  grade_perm <- sample(clin$grade)

  # 16 genes
  res_16 <- run_poly3F(expr_16, grade_perm, clin)
  n_perm_16[i] <- sum(res_16$p3_FDR < 0.05, na.rm = TRUE)

  # WGCNA module
  for (j in seq_along(mod_names)) {
    mod <- mod_names[j]
    expr_mod <- expr[mod_genes[[mod]], , drop = FALSE]
    res_mod <- run_poly3F(expr_mod, grade_perm, clin)
    n_perm_mod[i, j] <- sum(res_mod$p3_FDR < 0.05, na.rm = TRUE)
  }

  if (i %% 50 == 0) {
    elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
    log_info("  perm {i}/{N_PERM}, 耗时 {elapsed} min")
  }
}

elapsed_total <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
log_info("Permutation 完成，总耗时 {elapsed_total} min")

# ============================================================
# ============================================================
log_info("--- 6. 统计 ---")

perm_summary <- data.frame()

# 16 genes
p_perm_16 <- mean(n_perm_16 >= n_p3_real_16)
perm_summary <- rbind(perm_summary, data.frame(
  object        = "16_genes",
  n_genes       = length(genes_16_in),
  real_p3_sig   = n_p3_real_16,
  real_F_sig    = n_F_real_16,
  perm_mean     = round(mean(n_perm_16), 2),
  perm_sd       = round(sd(n_perm_16), 2),
  perm_max      = max(n_perm_16),
  perm_p95      = quantile(n_perm_16, 0.95),
  p_perm        = round(p_perm_16, 4),
  stringsAsFactors = FALSE))

# module
for (j in seq_along(mod_names)) {
  mod <- mod_names[j]
  real_val <- n_p3_real_mod[j]
  perm_dist <- n_perm_mod[, j]
  p_perm <- mean(perm_dist >= real_val)
  perm_summary <- rbind(perm_summary, data.frame(
    object        = mod,
    n_genes       = length(mod_genes[[mod]]),
    real_p3_sig   = real_val,
    real_F_sig    = NA_integer_,
    perm_mean     = round(mean(perm_dist), 2),
    perm_sd       = round(sd(perm_dist), 2),
    perm_max      = max(perm_dist),
    perm_p95      = quantile(perm_dist, 0.95),
    p_perm        = round(p_perm, 4),
    stringsAsFactors = FALSE))
}

save_table(perm_summary, DIR_THIS, "permutation_summary.csv")

cat("\n=== Permutation 汇总 ===\n")
print(perm_summary, row.names = FALSE)

dist_16 <- data.frame(perm_id = 1:N_PERM, n_p3_sig = n_perm_16)
save_table(dist_16, DIR_THIS, "permutation_dist_16genes.csv")

dist_mod <- data.frame(perm_id = 1:N_PERM, n_perm_mod)
save_table(dist_mod, DIR_THIS, "permutation_dist_modules.csv")

# ============================================================
# 7. figure
# ============================================================
log_info("--- 7. 图 ---")

pdf(file.path(DIR_THIS, "figures", "permutation_histogram.pdf"),
    width = 10, height = 5)
par(mfrow = c(1, 2))

hist(n_perm_16, breaks = seq(-0.5, max(n_perm_16) + 0.5, by = 1),
     col = "#B0B0B0", border = "white",
     main = "16 genes: permuted p3_FDR<0.05",
     xlab = "n genes with p3_FDR<0.05")
abline(v = n_p3_real_16, col = "red", lwd = 3)
legend("topright",
       legend = sprintf("real = %d\np_perm = %.4f",
                         n_p3_real_16, p_perm_16),
       col = "red", lwd = 3, bty = "n", cex = 1.1)

mod_colors <- c("greenyellow" = "#27ae60",
                "black" = "#2c3e50",
                "green" = "#16a085")
plot(NULL,
     xlim = c(0, max(n_perm_mod, na.rm = TRUE) + 5),
     ylim = c(0, 1.2),
     xlab = "n genes with p3_FDR<0.05",
     ylab = "density",
     main = "WGCNA modules: permuted p3_FDR<0.05")
for (j in seq_along(mod_names)) {
  mod <- mod_names[j]
  if (length(unique(n_perm_mod[, j])) > 1) {
    lines(density(n_perm_mod[, j], adjust = 2),
          col = mod_colors[mod], lwd = 2)
  }
  abline(v = n_p3_real_mod[j], col = mod_colors[mod], lwd = 3, lty = 2)
}
legend("topright",
       legend = paste0(mod_names,
                       " real=", n_p3_real_mod,
                       " (p=", perm_summary$p_perm[perm_summary$object %in% mod_names], ")"),
       col = mod_colors, lwd = 2, bty = "n", cex = 0.9)

dev.off()

png(file.path(DIR_THIS, "figures", "permutation_histogram.png"),
    width = 10, height = 5, units = "in", res = 200)
par(mfrow = c(1, 2))
hist(n_perm_16, breaks = seq(-0.5, max(n_perm_16) + 0.5, by = 1),
     col = "#B0B0B0", border = "white",
     main = "16 genes: permuted p3_FDR<0.05",
     xlab = "n genes with p3_FDR<0.05")
abline(v = n_p3_real_16, col = "red", lwd = 3)
legend("topright",
       legend = sprintf("real = %d\np_perm = %.4f",
                         n_p3_real_16, p_perm_16),
       col = "red", lwd = 3, bty = "n", cex = 1.1)

plot(NULL,
     xlim = c(0, max(n_perm_mod, na.rm = TRUE) + 5),
     ylim = c(0, 1.2),
     xlab = "n genes with p3_FDR<0.05",
     ylab = "density",
     main = "WGCNA modules: permuted p3_FDR<0.05")
for (j in seq_along(mod_names)) {
  mod <- mod_names[j]
  if (length(unique(n_perm_mod[, j])) > 1) {
    lines(density(n_perm_mod[, j], adjust = 2),
          col = mod_colors[mod], lwd = 2)
  }
  abline(v = n_p3_real_mod[j], col = mod_colors[mod], lwd = 3, lty = 2)
}
legend("topright",
       legend = paste0(mod_names,
                       " real=", n_p3_real_mod),
       col = mod_colors, lwd = 2, bty = "n", cex = 0.9)

dev.off()

log_info("图已保存")

# ============================================================
# 8. save RDS
# ============================================================
save_rds(list(
  perm_summary   = perm_summary,
  n_perm_16      = n_perm_16,
  n_perm_mod     = n_perm_mod,
  n_p3_real_16   = n_p3_real_16,
  n_p3_real_mod  = n_p3_real_mod,
  res_16_real    = res_16_real,
  res_mod_real   = res_mod_real,
  params = list(
    N_PERM = N_PERM,
    seed   = SEED,
    objects = c("16_genes", mod_names)
  )
), DIR_THIS, "permutation_results.rds")

# ============================================================
# 9. summary
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("=== 06_2 permutation test 汇总 ===\n")
cat(strrep("=", 70), "\n", sep = "")

cat("\n【Q1】Permutation 汇总\n")
print(perm_summary, row.names = FALSE)

cat("\n【Q2】关键判据\n")
for (i in seq_len(nrow(perm_summary))) {
  obj <- perm_summary$object[i]
  real_val <- perm_summary$real_p3_sig[i]
  p_val <- perm_summary$p_perm[i]
  verdict <- if (p_val < 0.01) "✓ 信号真实"
             else if (p_val < 0.05) "⚠ 边缘"
             else "✗ artifact / 过拟合"
  cat(sprintf("  %-15s  real=%-3d  p_perm=%.4f  %s\n",
              obj, real_val, p_val, verdict))
}

cat("\n【Q3】Permutation 分布特征\n")
cat(sprintf("  16 基因:  perm_mean=%.2f, perm_max=%d, real=%d\n",
            mean(n_perm_16), max(n_perm_16), n_p3_real_16))
for (j in seq_along(mod_names)) {
  cat(sprintf("  %-15s  perm_mean=%.2f, perm_max=%d, real=%d\n",
              mod_names[j], mean(n_perm_mod[, j]),
              max(n_perm_mod[, j]), n_p3_real_mod[j]))
}

cat("\n--- 输出位置 ---\n")
cat("目录:", DIR_THIS, "\n")
cat("  tables/permutation_summary.csv\n")
cat("  tables/permutation_dist_16genes.csv\n")
cat("  tables/permutation_dist_modules.csv\n")
cat("  figures/permutation_histogram.pdf/.png\n")

log_info("=== 06_2_permutation_test.r done ===")
