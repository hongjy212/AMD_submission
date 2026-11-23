#!/usr/bin/env Rscript
# =============================================================================
# 07_5_GSE231536_methylation.r (v5)
#
#
#
#
# Input：
#   E:/retina_tools/raw/GSE231536/
#     - GSE231536_M_values_retina.txt.gz
#     - eQTM_significant_pairs_retina.txt
#   results/01_1_prepare_bulk_data/data/gene_annotation.rds
#   results/01_2_prepare_gene_sets/data/gene_sets.rds
#   results/02_1_differential_expression_main/tables/sig16_genes.csv
#   results/04_2_pathway_score/tables/pathway_poly3F_adjusted.csv
#   results/05_1_consensus_wgcna/data/consensusNet.rds
#
# Output：results/07_5_GSE231536_methylation/
# =============================================================================

rm(list = ls()); gc()
setwd("E:/retina_tools/aging_rgc_axon/AMD")
set.seed(2026)

suppressPackageStartupMessages({
  library(data.table); library(readxl); library(dplyr)
  library(tidyr); library(ggplot2); library(logger)
})

out_dir <- "results/07_5_GSE231536_methylation"
tab_dir <- file.path(out_dir, "tables")
fig_dir <- file.path(out_dir, "figures")
for (d in c(out_dir, tab_dir, fig_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(out_dir, "07_5_GSE231536_methylation_v5.log")
log_appender(appender_tee(log_file))
log_layout(layout_glue_generator(format = '{time} [{level}] {msg}'))
log_info("=== 07_5 GSE231536 methylation v5 ===")

save_csv <- function(df, name) {
  write.csv(df, file.path(tab_dir, name),
            row.names = FALSE, fileEncoding = "UTF-8")
}

strip_v <- function(x) sub("\\..*$", "", x)

# ============================================================
# ============================================================
log_info("--- 1. 元数据 ---")
d_raw <- "E:/retina_tools/raw/GSE231536"

meta <- read_excel(file.path(d_raw, "41467_2024_46063_MOESM4_ESM.xlsx"),
                    sheet = "Sheet1")
log_info("原始样本: {nrow(meta)}")
print(table(meta$MGS_level, useNA = "ifany"))
print(table(meta$QC))

meta <- meta[meta$QC == "Passed", ]
log_info("QC Passed: {nrow(meta)}")
print(table(meta$MGS_level))

m_f <- file.path(d_raw, "GSE231536_M_values_retina.txt.gz")
con <- gzfile(m_f, "r"); header <- readLines(con, 1); close(con)
m_samples <- strsplit(header, "\t")[[1]][-1]

meta_matched <- meta[meta$ID %in% m_samples, ]
sample_order <- match(m_samples, meta_matched$ID)
grade_vec <- meta_matched$MGS_level[sample_order]

log_info("最终分级分布:")
print(table(grade_vec, useNA = "ifany"))

# ============================================================
# ============================================================
log_info("--- 2. M 值矩阵 ---")
m_dt <- fread(m_f, sep = "\t", header = TRUE, data.table = FALSE)
rownames(m_dt) <- m_dt[[1]]
m_dt <- m_dt[, -1]

m_mat <- as.matrix(m_dt)
log_info("M 值矩阵: {nrow(m_mat)} CpG × {ncol(m_mat)} 样本")
log_info("矩阵类型: {class(m_mat)[1]}")

# ============================================================
# ============================================================
log_info("--- 3. 基因集 ---")
GENE_SETS <- readRDS("results/01_2_prepare_gene_sets/data/gene_sets.rds")
polyF <- read.csv("results/04_2_pathway_score/tables/pathway_poly3F_adjusted.csv")
genes16 <- toupper(read.csv("results/02_1_differential_expression_main/tables/sig16_genes.csv")$gene)
p3_only_paths <- intersect(polyF$pathway[polyF$pattern == "p3_only"],
                            names(GENE_SETS))

pathway_gene_list <- list()
for (pw in p3_only_paths)
  pathway_gene_list[[pw]] <- toupper(GENE_SETS[[pw]]$genes)

target_genes <- unique(c(genes16, unlist(pathway_gene_list)))
log_info("16 基因: {length(genes16)}")
log_info("19 通路总基因: {length(unique(unlist(pathway_gene_list)))}")
log_info("目标基因总数: {length(target_genes)}")

anno <- readRDS("results/01_1_prepare_bulk_data/data/gene_annotation.rds")
ens2sym <- setNames(toupper(anno$symbol), strip_v(anno$ensembl))
sym2ens <- setNames(strip_v(anno$ensembl), toupper(anno$symbol))
target_ens <- unique(na.omit(sym2ens[target_genes]))
genes16_ens <- na.omit(sym2ens[genes16])
log_info("目标 ENSEMBL: {length(target_ens)}")
log_info("16 基因 ENSEMBL: {length(genes16_ens)}")

# ============================================================
# 4. eQTM match
# ============================================================
log_info("--- 4. eQTM ---")
eqtm <- read.delim(file.path(d_raw, "eQTM_significant_pairs_retina.txt"),
                    stringsAsFactors = FALSE)
log_info("eQTM 总对数: {nrow(eqtm)}")
log_info("beta 中位数: {round(median(eqtm$beta, na.rm=TRUE), 3)}")
log_info("★ eQTM 是 Beta regression：因变量是甲基化，beta 符号 = sign(E-M)")

eqtm_target <- eqtm[eqtm$trait %in% target_ens, ]
target_cpgs <- unique(eqtm_target$geno)
cpgs_in_matrix <- intersect(target_cpgs, rownames(m_mat))
log_info("目标 CpG 在矩阵: {length(cpgs_in_matrix)} / {length(target_cpgs)}")

if (length(cpgs_in_matrix) < 10) stop("CpG 太少")

# ============================================================
# 5. poly2 fit（MGS 1/3/4）
# ============================================================
log_info("--- 5. poly2 拟合 ---")

valid_idx <- which(!is.na(grade_vec) & grade_vec %in% c(1, 3, 4))
m_target <- m_mat[cpgs_in_matrix, valid_idx, drop = FALSE]
grade_valid <- as.numeric(grade_vec[valid_idx])
log_info("有效样本: {length(valid_idx)}")

fit_poly2_cpg <- function(y, g) {
  df <- data.frame(y = y, g = g)
  df <- df[complete.cases(df), ]
  if (nrow(df) < 20 || length(unique(df$g)) < 3)
    return(data.frame(p1 = NA, p2 = NA, R2 = NA, pattern = NA))
  fit <- lm(y ~ poly(g, 2), data = df)
  sm <- summary(fit)$coefficients
  if (nrow(sm) < 3) return(data.frame(p1 = NA, p2 = NA, R2 = NA, pattern = NA))
  p1 <- sm[2, 4]; p2 <- sm[3, 4]
  pattern <- ifelse(p2 < 0.05 & p1 > 0.05, "quadratic_only",
             ifelse(p1 < 0.05, "linear", "none"))
  data.frame(p1 = p1, p2 = p2, R2 = summary(fit)$r.squared, pattern = pattern)
}

poly2_cpg <- do.call(rbind, lapply(seq_len(nrow(m_target)), function(i) {
  cbind(cpg = rownames(m_target)[i], fit_poly2_cpg(m_target[i, ], grade_valid))
}))
poly2_cpg$p2_FDR <- p.adjust(poly2_cpg$p2, method = "BH")
poly2_cpg$p1_FDR <- p.adjust(poly2_cpg$p1, method = "BH")
poly2_cpg <- poly2_cpg[order(poly2_cpg$p2), ]

save_csv(poly2_cpg, "cpg_poly2_results.csv")
log_info("quadratic_only: {sum(poly2_cpg$pattern=='quadratic_only', na.rm=TRUE)}")
log_info("linear: {sum(poly2_cpg$pattern=='linear', na.rm=TRUE)}")
log_info("p2 FDR<0.05: {sum(poly2_cpg$p2_FDR<0.05, na.rm=TRUE)}")

# ============================================================
# 6. MGS1 vs MGS3 / MGS4 differentialmethylation
# ============================================================
log_info("--- 6. 差异甲基化 ---")

compare_delta <- function(group_a, group_b, label) {
  idx_a <- which(grade_vec == group_a)
  idx_b <- which(grade_vec == group_b)
  log_info("{label}: n_a={length(idx_a)}, n_b={length(idx_b)}")

  m_a <- m_mat[cpgs_in_matrix, idx_a, drop = FALSE]
  m_b <- m_mat[cpgs_in_matrix, idx_b, drop = FALSE]

  delta <- rowMeans(m_b, na.rm = TRUE) - rowMeans(m_a, na.rm = TRUE)

  pvals <- sapply(seq_len(nrow(m_a)), function(i) {
    tryCatch(wilcox.test(m_b[i, ], m_a[i, ])$p.value,
             error = function(e) NA_real_)
  })

  data.frame(cpg = rownames(m_a), delta_M = delta,
             p = pvals, p_FDR = p.adjust(pvals, method = "BH"),
             comparison = label, stringsAsFactors = FALSE)
}

diff_m13 <- compare_delta(1, 3, "MGS1_vs_MGS3")
diff_m14 <- compare_delta(1, 4, "MGS1_vs_MGS4")

save_csv(rbind(diff_m13, diff_m14), "deltaM_MGS1_vs_MGS3_MGS4.csv")
log_info("MGS1 vs MGS3 FDR<0.05: {sum(diff_m13$p_FDR<0.05, na.rm=TRUE)}")
log_info("MGS1 vs MGS4 FDR<0.05: {sum(diff_m14$p_FDR<0.05, na.rm=TRUE)}")

# CpG -> genesmapping
for (obj in c("diff_m13", "diff_m14")) {
  df <- get(obj)
  df$gene_ens <- eqtm_target$trait[match(df$cpg, eqtm_target$geno)]
  df$gene_sym <- toupper(names(sym2ens)[match(df$gene_ens, sym2ens)])
  assign(obj, df)
}

gene_diff13 <- diff_m13[!is.na(diff_m13$gene_sym), ] %>%
  group_by(gene_sym) %>%
  summarise(n_cpg = n(), mean_delta_M = mean(delta_M, na.rm = TRUE),
            n_FDR_sig = sum(p_FDR < 0.05, na.rm = TRUE),
            best_p = min(p, na.rm = TRUE), .groups = "drop") %>%
  arrange(best_p)
save_csv(gene_diff13, "gene_MGS1_vs_MGS3_deltaM.csv")

gene_diff14 <- diff_m14[!is.na(diff_m14$gene_sym), ] %>%
  group_by(gene_sym) %>%
  summarise(n_cpg = n(), mean_delta_M = mean(delta_M, na.rm = TRUE),
            n_FDR_sig = sum(p_FDR < 0.05, na.rm = TRUE),
            best_p = min(p, na.rm = TRUE), .groups = "drop") %>%
  arrange(best_p)
save_csv(gene_diff14, "gene_MGS1_vs_MGS4_deltaM.csv")

# ============================================================
# ============================================================
log_info("--- 7. 基因级 eQTM 方向验证 ---")

validate_direction <- function(diff_df, gene_ens_set, label) {
  cpgs_sel <- unique(eqtm_target$geno[eqtm_target$trait %in% gene_ens_set])
  sub <- diff_df[diff_df$cpg %in% cpgs_sel, ]

  merged <- merge(sub, eqtm_target[, c("geno", "beta", "trait")],
                   by.x = "cpg", by.y = "geno", all.x = TRUE)
  merged <- merged[!is.na(merged$beta) & merged$beta != 0, ]
  if (nrow(merged) == 0) { log_info("{label}: 无有效 CpG"); return(NULL) }

  merged$consistent <- sign(merged$delta_M) == sign(merged$beta)
  merged$gene_sym <- toupper(names(sym2ens)[match(merged$trait, sym2ens)])

  n_total <- nrow(merged); n_cons <- sum(merged$consistent, na.rm = TRUE)
  log_info("{label}: {n_cons} / {n_total} = {round(100*n_cons/n_total, 1)}%")
  merged
}

val_g16_m13 <- validate_direction(diff_m13, genes16_ens, "16genes MGS1_vs_MGS3")
val_g16_m14 <- validate_direction(diff_m14, genes16_ens, "16genes MGS1_vs_MGS4")
val_path_m13 <- validate_direction(diff_m13, target_ens, "19pathways MGS1_vs_MGS3")
val_path_m14 <- validate_direction(diff_m14, target_ens, "19pathways MGS1_vs_MGS4")

if (!is.null(val_g16_m13))  save_csv(val_g16_m13,  "gene16_MGS1_vs_MGS3_eQTM.csv")
if (!is.null(val_g16_m14))  save_csv(val_g16_m14,  "gene16_MGS1_vs_MGS4_eQTM.csv")
if (!is.null(val_path_m13)) save_csv(val_path_m13, "pathway_MGS1_vs_MGS3_eQTM.csv")
if (!is.null(val_path_m14)) save_csv(val_path_m14, "pathway_MGS1_vs_MGS4_eQTM.csv")

# ============================================================
# ============================================================
log_info("--- 8. 模块级甲基化验证 ---")

net_1 <- readRDS("results/05_1_consensus_wgcna/data/consensusNet.rds")

module_genes <- list(
  greenyellow = names(net_1$colors)[net_1$colors == "greenyellow"],
  black       = names(net_1$colors)[net_1$colors == "black"],
  green       = names(net_1$colors)[net_1$colors == "green"]
)
for (mod in names(module_genes))
  log_info("{mod}: {length(module_genes[[mod]])} 基因")

# 8a. module CpG overwrite
module_cpgs <- list()
for (mod in names(module_genes)) {
  ens <- strip_v(module_genes[[mod]])
  ens <- intersect(ens, eqtm$trait)
  cpgs <- unique(eqtm$geno[eqtm$trait %in% ens])
  cpgs <- intersect(cpgs, rownames(m_mat))
  module_cpgs[[mod]] <- cpgs
  log_info("{mod}: {length(ens)} 基因有 eQTM, {length(cpgs)} CpG")
}

# 8b. modulemethylationscore
module_scores <- data.frame(
  sample_id = colnames(m_mat), grade = grade_vec,
  stringsAsFactors = FALSE)

for (mod in names(module_cpgs)) {
  cpgs <- module_cpgs[[mod]]
  if (length(cpgs) < 5) next
  score <- colMeans(m_mat[cpgs, , drop = FALSE], na.rm = TRUE)
  module_scores[[paste0("M_", mod)]] <- as.numeric(score)
  log_info("{mod} 评分: mean={round(mean(score, na.rm=TRUE), 3)}")
}

# 8c. moduledifferential（MGS1 vs MGS3 / MGS4）
compare_module <- function(score, grade, a, b) {
  idx_a <- which(grade == a); idx_b <- which(grade == b)
  if (length(idx_a) < 5 || length(idx_b) < 5) return(NULL)
  s_a <- score[idx_a]; s_b <- score[idx_b]
  wt <- wilcox.test(s_b, s_a)
  data.frame(mean_a = mean(s_a, na.rm=TRUE), mean_b = mean(s_b, na.rm=TRUE),
             delta_M = mean(s_b, na.rm=TRUE) - mean(s_a, na.rm=TRUE),
             p = wt$p.value, stringsAsFactors = FALSE)
}

module_diff <- data.frame()
for (mod in names(module_cpgs)) {
  col <- paste0("M_", mod)
  if (!col %in% colnames(module_scores)) next
  r13 <- compare_module(module_scores[[col]], grade_vec, 1, 3)
  r14 <- compare_module(module_scores[[col]], grade_vec, 1, 4)
  module_diff <- rbind(module_diff, data.frame(
    module = mod, n_cpg = length(module_cpgs[[mod]]),
    MGS1_mean = r13$mean_a, MGS3_mean = r13$mean_b,
    delta_M_13 = r13$delta_M, p_13 = r13$p,
    MGS4_mean = r14$mean_b, delta_M_14 = r14$delta_M, p_14 = r14$p,
    stringsAsFactors = FALSE))
}
save_csv(module_diff, "module_methylation_diff.csv")
log_info("模块甲基化差异:")
print(module_diff, row.names = FALSE, digits = 4)

# ============================================================
# ============================================================
log_info("--- 8d. 模块方向验证 ---")

module_direction <- data.frame()
idx_a <- which(grade_vec == 1); idx_b <- which(grade_vec == 3)

for (mod in names(module_cpgs)) {
  cpgs <- module_cpgs[[mod]]
  if (length(cpgs) < 5) next

  sub <- eqtm[eqtm$geno %in% cpgs, c("geno", "beta", "trait")]
  ens <- strip_v(module_genes[[mod]])
  sub <- sub[sub$trait %in% ens, ]
  sub <- sub[!duplicated(sub$geno), ]
  if (nrow(sub) < 5) next

  m_sub <- m_mat[sub$geno, , drop = FALSE]
  delta <- rowMeans(m_sub[, idx_b, drop = FALSE], na.rm = TRUE) -
           rowMeans(m_sub[, idx_a, drop = FALSE], na.rm = TRUE)

  sub$delta_M <- as.numeric(delta)

  sub$consistent <- sign(sub$delta_M) == sign(sub$beta)

  n_na <- sum(is.na(sub$consistent))
  n_total <- nrow(sub); n_cons <- sum(sub$consistent, na.rm = TRUE)
  pct <- 100 * n_cons / max(n_total - n_na, 1)
  binom_p <- binom.test(n_cons, n_total - n_na, 0.5)$p.value

  n_delta_pos <- sum(sub$delta_M > 0, na.rm=TRUE)
  n_delta_neg <- sum(sub$delta_M < 0, na.rm=TRUE)
  n_beta_pos <- sum(sub$beta > 0, na.rm=TRUE)
  n_beta_neg <- sum(sub$beta < 0, na.rm=TRUE)

  log_info("{mod}: {n_cons}/{n_total} = {round(pct,1)}% (binom p={format(binom_p, scientific=TRUE)})")
  log_info("  delta_M: +{n_delta_pos} / -{n_delta_neg}; beta: +{n_beta_pos} / -{n_beta_neg}")

  module_direction <- rbind(module_direction, data.frame(
    module = mod, n_cpg = n_total, n_consistent = n_cons,
    pct = round(pct, 1), binom_p = binom_p,
    n_delta_pos = n_delta_pos, n_delta_neg = n_delta_neg,
    n_beta_pos = n_beta_pos, n_beta_neg = n_beta_neg,
    stringsAsFactors = FALSE))
}
save_csv(module_direction, "module_direction_validation.csv")

bulk_dir <- c(greenyellow = "down", black = "up", green = "down")
module_direction$bulk_dir <- bulk_dir[module_direction$module]
module_direction <- merge(module_direction,
                           module_diff[, c("module", "delta_M_13", "p_13")],
                           by = "module", all.x = TRUE)
module_direction$methyl_dir <- ifelse(module_direction$delta_M_13 > 0, "up", "down")

#          bulk up   → methylation up   → match
module_direction$match_bulk <- with(module_direction,
  (bulk_dir == "down" & methyl_dir == "down") |
  (bulk_dir == "up"   & methyl_dir == "up"))

save_csv(module_direction, "module_direction_with_bulk.csv")
log_info("与 bulk p3 对照:")
print(module_direction[, c("module", "n_cpg", "n_consistent", "pct",
                            "delta_M_13", "p_13", "bulk_dir",
                            "methyl_dir", "match_bulk")],
      row.names = FALSE, digits = 4)

# ============================================================
# 9. visualization
# ============================================================
log_info("--- 9. 可视化 ---")

p1 <- ggplot(diff_m13, aes(x = delta_M)) +
  geom_histogram(bins = 50, fill = "steelblue", color = "white") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "MGS1 vs MGS3 delta M", x = "delta M", y = "Count") +
  theme_minimal(base_size = 12)
ggsave(file.path(fig_dir, "deltaM_MGS1_vs_MGS3.pdf"), p1, width = 8, height = 5)

if (!is.null(val_g16_m13) && nrow(val_g16_m13) > 0) {
  p3 <- ggplot(val_g16_m13, aes(x = beta, y = delta_M, color = consistent)) +
    geom_point(size = 2, alpha = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_color_manual(values = c("TRUE"="#4C72B0", "FALSE"="#C44E52")) +
    labs(title = "16-gene CpGs: eQTM beta vs delta M",
         x = "eQTM beta", y = "delta M", color = "Consistent") +
    theme_minimal(base_size = 12)
  ggsave(file.path(fig_dir, "gene16_direction.pdf"), p3, width = 8, height = 6)
}

plot_df <- module_scores %>%
  pivot_longer(cols = starts_with("M_"),
               names_to = "module", values_to = "score") %>%
  mutate(module = sub("^M_", "", module)) %>%
  filter(!is.na(grade), grade %in% c(1, 3, 4))

p4 <- ggplot(plot_df, aes(x = factor(grade), y = score, fill = factor(grade))) +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(~ module, scales = "free_y", ncol = 3) +
  labs(title = "Module methylation score by MGS grade",
       x = "MGS grade", y = "Mean M value", fill = "MGS") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave(file.path(fig_dir, "module_methylation_score.pdf"), p4,
       width = 12, height = 5)

if (nrow(module_direction) > 0) {
  p5 <- ggplot(module_direction, aes(x = module, y = pct, fill = pct > 50)) +
    geom_col() +
    geom_hline(yintercept = 50, linetype = "dashed", color = "red") +
    geom_text(aes(label = sprintf("%d/%d", n_consistent, n_cpg)),
              vjust = -0.5, size = 3.5) +
    labs(title = "Module eQTM direction consistency (MGS1 vs MGS3)",
         x = "Module", y = "% consistent", fill = ">50%") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  ggsave(file.path(fig_dir, "module_direction.pdf"), p5, width = 8, height = 5)
}

# ============================================================
# 10. summary
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("=== 07_5 v5 汇总 ===\n")
cat(strrep("=", 70), "\n", sep = "")

cat("\n--- 数据 ---\n")
cat("样本数:", ncol(m_mat), " CpG 数:", nrow(m_mat), "\n")
cat("MGS 分级:\n"); print(table(grade_vec, useNA = "ifany"))

cat("\n--- poly2 拟合 ---\n")
cat("quadratic_only:", sum(poly2_cpg$pattern=="quadratic_only", na.rm=TRUE), "\n")
cat("p2 FDR<0.05:", sum(poly2_cpg$p2_FDR<0.05, na.rm=TRUE), "\n")

cat("\n--- 差异甲基化 ---\n")
cat("MGS1 vs MGS3 FDR<0.05:", sum(diff_m13$p_FDR<0.05, na.rm=TRUE), "/", nrow(diff_m13), "\n")
cat("MGS1 vs MGS4 FDR<0.05:", sum(diff_m14$p_FDR<0.05, na.rm=TRUE), "/", nrow(diff_m14), "\n")

cat("\n--- 基因级方向验证 ---\n")
for (obj_name in c("val_g16_m13", "val_g16_m14", "val_path_m13", "val_path_m14")) {
  obj <- get(obj_name)
  if (!is.null(obj) && nrow(obj) > 0) {
    n_cons <- sum(obj$consistent, na.rm=TRUE)
    cat(sprintf("  %-28s  %d / %d = %.1f%%\n",
                obj_name, n_cons, nrow(obj), 100*n_cons/nrow(obj)))
  }
}

cat("\n--- 模块级验证 ---\n")
for (mod in names(module_cpgs))
  cat(sprintf("  %-12s  %d CpG\n", mod, length(module_cpgs[[mod]])))

cat("\n模块甲基化差异:\n")
print(module_diff, row.names = FALSE, digits = 4)

cat("\n模块方向一致率:\n")
print(module_direction[, c("module", "n_cpg", "n_consistent", "pct",
                            "binom_p", "delta_M_13", "p_13",
                            "bulk_dir", "methyl_dir", "match_bulk")],
      row.names = FALSE, digits = 4)

cat("\n方向分布诊断:\n")
print(module_direction[, c("module", "n_delta_pos", "n_delta_neg",
                            "n_beta_pos", "n_beta_neg")],
      row.names = FALSE)

cat("\n目录:", out_dir, "\n")
log_info("=== 07_5 v5 done ===")


