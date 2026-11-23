#!/usr/bin/env Rscript
# ============================================================
# 07_2_gse29801_replication_v5.r
#
#
#   4. gene set overwritediagnosis
#
#   normal  -> G1
#   dry AMD -> G3
#   GA      -> G4
#   CNV     -> G4
#   GA/CNV  -> G4
#
# Input：
#   E:/retina_tools/raw/GSE29801/GSE29801_series_matrix.txt.gz
#   E:/retina_tools/raw/GSE29801/GPL4133_old_annotations.txt.gz
#   results/02_1_differential_expression_main/tables/sig16_genes.csv
#   results/05_1_consensus_wgcna/tables/module_genes_greenyellow.csv
#   results/05_1_consensus_wgcna/tables/module_genes_black.csv
#   results/02_1_differential_expression_main/tables/limma_multi_full.csv
#
# Output：results/07_2_gse29801_replication/
#   data/expr_log2.rds
#   data/meta.rds
#   tables/limma_multi.csv
#   tables/poly3_F.csv
#   tables/gene_set_direction.csv
#   tables/permutation_results.csv
#   tables/consistency_with_GSE115828.csv
#   figures/poly3_F_volcano.pdf
#   figures/gene_set_distribution.pdf
# ============================================================

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(ggplot2)
  library(limma)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

DIR_THIS <- setup_output_dir("07_2_gse29801_replication")
log_info("=== 07_2_gse29801_replication.r v5 start ===")

# ============================================================
# 0. Parameters
# ============================================================
GSE_DIR      <- "E:/retina_tools/raw/GSE29801"
GPL_SKIP     <- 21
AGE_MIN      <- 60
N_PERM       <- 1000
SEED         <- 42

GRADE_MAP <- c(
  "normal"  = "G1",
  "MD1"     = "G1",      # 硬性玻璃膜疣，属正常老化
  "MD2"     = "G2",      # 软性玻璃膜疣，早期 AMD
  "dry AMD" = "G3",
  "GA"      = "G4",
  "CNV"     = "G4",
  "GA/CNV"  = "G4"
)

# ============================================================
# ============================================================
log_info("--- 1. 读 GPL4133 注释 ---")

gpl <- fread(file.path(GSE_DIR, "GPL4133_old_annotations.txt.gz"),
             sep = "\t", header = TRUE, skip = GPL_SKIP,
             quote = "", showProgress = FALSE)

log_info(sprintf("GPL: %d probes x %d cols", nrow(gpl), ncol(gpl)))

gpl_clean <- gpl[!is.na(GENE_SYMBOL) & GENE_SYMBOL != "",
                  .(ID, NAME, GENE_SYMBOL, GENE_NAME, ENSEMBL_ID)]
gpl_clean[, ID := as.integer(ID)]
gpl_clean[, GENE_SYMBOL := toupper(trimws(GENE_SYMBOL))]

log_info(sprintf("有 SYMBOL 的 probe: %d, unique genes: %d",
                 nrow(gpl_clean), uniqueN(gpl_clean$GENE_SYMBOL)))

for (g in c("NDUFB8", "AQR", "UGGT1")) {
  n <- sum(gpl_clean$GENE_SYMBOL == g)
  log_info(sprintf("  %s: %d probes", g, n))
}

# ============================================================
# ============================================================
log_info("--- 2. 读 series matrix ---")

sm_file  <- file.path(GSE_DIR, "GSE29801_series_matrix.txt.gz")
sm_lines <- readLines(gzfile(sm_file))

meta_lines <- grep("^!Sample_characteristics_ch1\t", sm_lines, value = TRUE)

get_field <- function(pattern) {
  idx <- grep(pattern, meta_lines)
  if (length(idx) == 0) return(NULL)
  p <- strsplit(meta_lines[idx[1]], "\t", fixed = TRUE)[[1]]
  vals <- gsub('^"|"$', "", p[-1])
  sub(paste0("^", pattern), "", vals)
}

meta_df <- data.frame(
  sample    = gsub('^"|"$', "",
                   strsplit(grep("^!Sample_geo_accession", sm_lines,
                                 value = TRUE)[1],
                            "\t", fixed = TRUE)[[1]][-1]),
  tissue    = get_field("tissue: "),
  gender    = get_field("gender: "),
  age       = as.integer(get_field("age \\(years\\): ")),
  ocular    = get_field("ocular disease: "),
  grade_raw = get_field("amd classification: "),
  patient   = get_field("patient identifier: "),
  rin       = get_field("rna integrity number \\(rin\\): "),
  stringsAsFactors = FALSE
)

log_info(sprintf("样本数: %d", nrow(meta_df)))

begin_idx <- grep("^!series_matrix_table_begin", sm_lines)[1]
end_idx   <- grep("^!series_matrix_table_end",   sm_lines)[1]

expr_dt <- fread(text = sm_lines[(begin_idx + 1):(end_idx - 1)],
                 sep = "\t", header = TRUE, quote = "",
                 showProgress = FALSE)
colnames(expr_dt) <- gsub('"', '', colnames(expr_dt))

log_info(sprintf("表达矩阵: %d probes x %d samples",
                 nrow(expr_dt), ncol(expr_dt) - 1))

# ============================================================
# ============================================================
log_info("--- 3. 样本过滤 ---")

meta_df$grade <- GRADE_MAP[meta_df$grade_raw]
meta_df$is_retina <- grepl("Retina", meta_df$tissue) &
                     !grepl("RPE-choroid", meta_df$tissue)
meta_df$region <- ifelse(grepl("extramacular", meta_df$tissue, ignore.case = TRUE),
                         "Extramacular",
                  ifelse(grepl("macular",       meta_df$tissue, ignore.case = TRUE),
                         "Macular", NA))

keep <- !is.na(meta_df$grade) &
        meta_df$is_retina &
        !is.na(meta_df$age) &
        meta_df$age >= AGE_MIN
meta_filt <- meta_df[keep, ]

log_info(sprintf("过滤后样本数: %d", nrow(meta_filt)))
log_info("--- 新映射后 grade 分布 ---")
print(table(meta_filt$grade))
print(table(meta_filt$grade, meta_filt$region))
log_info("--- 原始 grade_raw 分布 ---")
print(table(meta_filt$grade_raw))

n_match <- length(intersect(meta_filt$sample, colnames(expr_dt)))
log_info(sprintf("样本与表达矩阵匹配: %d / %d", n_match, nrow(meta_filt)))
stopifnot(n_match == nrow(meta_filt))

save_rds(meta_filt, DIR_THIS, "meta.rds")

# ============================================================
# 4. Probe -> Gene
# ============================================================
log_info("--- 4. Probe -> Gene 映射 ---")

sample_cols <- c("ID_REF", meta_filt$sample)
expr_sub <- expr_dt[, ..sample_cols]
expr_sub[, ID_REF := as.integer(ID_REF)]

setnames(gpl_clean, "ID", "ID_REF")
expr_anno <- merge(expr_sub, gpl_clean[, .(ID_REF, GENE_SYMBOL)],
                   by = "ID_REF", all.x = FALSE)
log_info(sprintf("合并后 probe: %d", nrow(expr_anno)))

if (nrow(expr_anno) == 0) stop("merge 后 0 行")

expr_mat <- as.matrix(expr_anno[, -c("ID_REF", "GENE_SYMBOL")])
qs <- quantile(expr_mat, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE)
log_info(sprintf("表达值分位数: %s", paste(round(qs, 2), collapse = " ")))

if (qs["100%"] > 50) {
  log_info("最大值 > 50，进行 log2(x+1) 变换")
  expr_mat <- log2(expr_mat + 1)
}

expr_dt2 <- data.table(GENE_SYMBOL = expr_anno$GENE_SYMBOL, expr_mat)
expr_gene <- expr_dt2[, lapply(.SD, mean, na.rm = TRUE), by = GENE_SYMBOL]
expr_final <- as.matrix(expr_gene[, -1])
rownames(expr_final) <- expr_gene$GENE_SYMBOL

log_info(sprintf("基因表达矩阵: %d genes x %d samples",
                 nrow(expr_final), ncol(expr_final)))
save_rds(expr_final, DIR_THIS, "expr_log2.rds")

# ============================================================
# 5. limma multi
# ============================================================
log_info("--- 5. limma multi ---")

meta_filt$grade <- factor(meta_filt$grade, levels = c("G1", "G2", "G3", "G4"))

design_df <- data.frame(
  grade  = meta_filt$grade,
  age    = meta_filt$age,
  sex_M  = as.integer(meta_filt$gender == "male"),
  region = factor(meta_filt$region),
  row.names = meta_filt$sample
)

has_region <- length(unique(design_df$region)) > 1
if (has_region) {
  design_m <- model.matrix(~ 0 + grade + age + sex_M + region, data = design_df)
} else {
  log_warn("region 只有一个水平，从 design 中移除")
  design_m <- model.matrix(~ 0 + grade + age + sex_M, data = design_df)
}
colnames(design_m)[1:4] <- c("g1", "g2", "g3", "g4")

cm <- makeContrasts(
  g2_vs_g1 = g2 - g1,
  g3_vs_g2 = g3 - g2,
  g4_vs_g3 = g4 - g3,
  g4_vs_g1 = g4 - g1,
  levels = design_m
)

fit <- eBayes(contrasts.fit(lmFit(expr_final, design_m), cm))

multi_res <- data.frame()
for (ct in colnames(cm)) {
  tt <- topTable(fit, coef = ct, number = Inf, sort.by = "none")
  multi_res <- rbind(multi_res, data.frame(
    gene = rownames(tt), contrast = ct,
    logFC = tt$logFC, p = tt$P.Value, FDR = tt$adj.P.Val,
    stringsAsFactors = FALSE))
}
save_table(multi_res, DIR_THIS, "limma_multi.csv")

cat("\n=== limma multi FDR<0.05 ===\n")
for (ct in colnames(cm)) {
  n <- sum(multi_res$FDR[multi_res$contrast == ct] < 0.05, na.rm = TRUE)
  cat(sprintf("  %-10s  %d\n", ct, n))
}

# ============================================================
# 6. poly3_F
# ============================================================
log_info("--- 6. poly3_F ---")

poly_mat <- poly(as.integer(meta_filt$grade), 3)
colnames(poly_mat) <- c("p1", "p2", "p3")
rownames(poly_mat) <- meta_filt$sample

design_p <- cbind(
  "(Intercept)" = 1,
  p1 = poly_mat[, 1],
  p2 = poly_mat[, 2],
  p3 = poly_mat[, 3],
  age = design_df$age,
  sex_M = design_df$sex_M)
if (has_region) {
  design_p <- cbind(design_p,
                    regionExtramacular = as.integer(design_df$region == "Extramacular"))
}

fit_p <- eBayes(lmFit(expr_final, design_p))
coef_p <- c("p1", "p2", "p3")
tt_F <- topTable(fit_p, coef = coef_p, number = Inf, sort.by = "none")
polyF <- data.frame(
  gene = rownames(tt_F),
  F    = tt_F$F, p = tt_F$P.Value, FDR = tt_F$adj.P.Val,
  stringsAsFactors = FALSE)
save_table(polyF, DIR_THIS, "poly3_F.csv")

log_info(sprintf("poly3_F: FDR<0.05 = %d", sum(polyF$FDR < 0.05, na.rm = TRUE)))
log_info(sprintf("poly3_F: p<0.001   = %d", sum(polyF$p < 0.001, na.rm = TRUE)))

# ============================================================
# ============================================================
log_info("--- 7. 基因集方向一致性 ---")

read_set <- function(orig_path, fix_path) {
  if (file.exists(fix_path)) {
    df <- read.csv(fix_path)
    if ("SYMBOL" %in% colnames(df)) {
      log_info(sprintf("读取 %s", basename(fix_path)))
      sym <- toupper(df$SYMBOL)
      return(sym[!is.na(sym) & sym != ""])
    }
  }
  if (!file.exists(orig_path)) {
    log_warn(sprintf("未找到 %s", basename(orig_path)))
    return(NULL)
  }
  df <- read.csv(orig_path)
  if (!"ENSEMBL" %in% colnames(df)) {
    log_warn(sprintf("%s 无 ENSEMBL 列", basename(orig_path)))
    return(NULL)
  }
  ens <- sub("\\..*$", "", df$ENSEMBL)
  sym <- mapIds(org.Hs.eg.db, keys = ens, column = "SYMBOL",
                keytype = "ENSEMBL", multiVals = "first")
  sym <- toupper(sym)
  log_info(sprintf("现场反查 %s: %d -> %d",
                   basename(orig_path), length(ens), sum(!is.na(sym))))
  sym[!is.na(sym) & sym != ""]
}

genes_16 <- toupper(read_upstream("02_1_differential_expression_main", "tables",
                                  "sig16_genes.csv")$gene)

tbl_dir <- file.path(DIR_RESULT, "05_1_consensus_wgcna", "tables")
genes_gy <- read_set(file.path(tbl_dir, "module_genes_greenyellow.csv"),
                     file.path(tbl_dir, "module_genes_greenyellow_symbolfix.csv"))
genes_bl <- read_set(file.path(tbl_dir, "module_genes_black.csv"),
                     file.path(tbl_dir, "module_genes_black_symbolfix.csv"))

log_info(sprintf("gene set 大小: 16genes=%d, greenyellow=%d, black=%d",
                 length(genes_16),
                 ifelse(is.null(genes_gy), 0, length(genes_gy)),
                 ifelse(is.null(genes_bl), 0, length(genes_bl))))

gene_sets <- list(
  "16genes"     = genes_16,
  "greenyellow" = genes_gy,
  "black"       = genes_bl)
gene_sets <- gene_sets[!sapply(gene_sets, is.null)]

for (sname in names(gene_sets)) {
  n_in <- sum(gene_sets[[sname]] %in% rownames(expr_final))
  log_info(sprintf("%s: %d/%d 在表达矩阵中",
                   sname, n_in, length(gene_sets[[sname]])))
}

check_set <- function(gset, set_name) {
  sub <- multi_res[multi_res$gene %in% gset, ]
  if (nrow(sub) == 0) return(NULL)

  out <- data.frame()
  for (ct in c("g2_vs_g1", "g3_vs_g2", "g4_vs_g3")) {
    s <- sub[sub$contrast == ct, ]
    if (nrow(s) == 0) next
    out <- rbind(out, data.frame(
      gene_set = set_name,
      contrast = ct,
      n        = nrow(s),
      mean_logFC = round(mean(s$logFC, na.rm = TRUE), 4),
      median_logFC = round(median(s$logFC, na.rm = TRUE), 4),
      wilcox_p = if (nrow(s) >= 5) wilcox.test(s$logFC, mu = 0)$p.value else NA_real_,
      n_sig = sum(s$FDR < 0.05, na.rm = TRUE),
      stringsAsFactors = FALSE))
  }
  out
}

direction_tbl <- data.frame()
for (sname in names(gene_sets)) {
  r <- check_set(gene_sets[[sname]], sname)
  if (!is.null(r)) direction_tbl <- rbind(direction_tbl, r)
}
save_table(direction_tbl, DIR_THIS, "gene_set_direction.csv")

cat("\n=== 基因集方向一致性 ===\n")
print(direction_tbl, row.names = FALSE)

# ============================================================
# 8. Permutation
# ============================================================
log_info("--- 8. Permutation ---")

set.seed(SEED)

perm_test <- function(gset, set_name, contrast_name) {
  res_ct <- multi_res[multi_res$contrast == contrast_name, ]
  sub <- res_ct[res_ct$gene %in% gset, ]
  if (nrow(sub) < 5) return(NULL)

  obs_mean <- mean(sub$logFC, na.rm = TRUE)
  all_genes <- res_ct$gene

  perm_means <- replicate(N_PERM, {
    g <- sample(all_genes, length(gset), replace = FALSE)
    mean(res_ct[res_ct$gene %in% g, ]$logFC, na.rm = TRUE)
  })

  data.frame(
    gene_set     = set_name,
    contrast     = contrast_name,
    n            = nrow(sub),
    obs_mean_logFC = round(obs_mean, 4),
    perm_med     = round(median(perm_means), 4),
    perm_q025    = round(quantile(perm_means, 0.025), 4),
    perm_q975    = round(quantile(perm_means, 0.975), 4),
    perm_p       = mean(abs(perm_means) >= abs(obs_mean)),
    stringsAsFactors = FALSE)
}

perm_tbl <- data.frame()
for (sname in names(gene_sets)) {
  for (ct in c("g2_vs_g1", "g3_vs_g2", "g4_vs_g3")) {
    r <- perm_test(gene_sets[[sname]], sname, ct)
    if (!is.null(r)) perm_tbl <- rbind(perm_tbl, r)
  }
}
save_table(perm_tbl, DIR_THIS, "permutation_results.csv")

cat("\n=== Permutation ===\n")
print(perm_tbl, row.names = FALSE)

# ============================================================
# ============================================================
log_info("--- 9. 与 GSE115828 一致性 ---")

gse115_multi <- tryCatch(
  read.csv(file.path(DIR_RESULT, "02_1_differential_expression_main",
                     "tables", "limma_multi_full.csv")),
  error = function(e) NULL)

consistency_tbl <- data.frame()
if (!is.null(gse115_multi)) {
  for (ct in c("g2_vs_g1", "g3_vs_g2", "g4_vs_g3")) {
    gse115_ct <- gse115_multi[gse115_multi$contrast == ct, ]
    gse298_ct <- multi_res[multi_res$contrast == ct, ]

    shared <- intersect(gse115_ct$gene, gse298_ct$gene)
    if (length(shared) < 10) next

    merged <- merge(
      gse115_ct[gse115_ct$gene %in% shared, c("gene", "logFC")],
      gse298_ct[gse298_ct$gene %in% shared, c("gene", "logFC")],
      by = "gene", suffixes = c("_GSE115828", "_GSE29801"))

    r_all <- cor(merged$logFC_GSE115828, merged$logFC_GSE29801,
                 method = "spearman", use = "complete.obs")

    log_info(sprintf("[%s] 共享基因 %d, Spearman r = %.3f",
                     ct, nrow(merged), r_all))

    consistency_tbl <- rbind(consistency_tbl, data.frame(
      contrast    = ct,
      gene_set    = "ALL_shared",
      n           = nrow(merged),
      spearman_r  = round(r_all, 3),
      n_same_dir  = sum(sign(merged$logFC_GSE115828) == sign(merged$logFC_GSE29801)),
      pct_same_dir = round(100 * mean(sign(merged$logFC_GSE115828) ==
                                       sign(merged$logFC_GSE29801)), 1),
      stringsAsFactors = FALSE))

    for (sname in names(gene_sets)) {
      sub <- merged[merged$gene %in% gene_sets[[sname]], ]
      if (nrow(sub) < 5) next
      r <- cor(sub$logFC_GSE115828, sub$logFC_GSE29801,
               method = "spearman", use = "complete.obs")
      n_same <- sum(sign(sub$logFC_GSE115828) == sign(sub$logFC_GSE29801))
      consistency_tbl <- rbind(consistency_tbl, data.frame(
        contrast    = ct,
        gene_set    = sname,
        n           = nrow(sub),
        spearman_r  = round(r, 3),
        n_same_dir  = n_same,
        pct_same_dir = round(100 * n_same / nrow(sub), 1),
        stringsAsFactors = FALSE))
    }
  }

  if (nrow(consistency_tbl) > 0) {
    save_table(consistency_tbl, DIR_THIS, "consistency_with_GSE115828.csv")
    cat("\n=== 与 GSE115828 一致性 ===\n")
    print(consistency_tbl, row.names = FALSE)
  }
}

# ============================================================
# 10. figure
# ============================================================
log_info("--- 10. 图 ---")

plot_df <- polyF
plot_df$neglog10p <- -log10(pmax(plot_df$p, 1e-300))
plot_df$highlight <- "other"
plot_df$highlight[plot_df$gene %in% genes_16] <- "16genes"
if (!is.null(genes_gy)) plot_df$highlight[plot_df$gene %in% genes_gy] <- "greenyellow"
if (!is.null(genes_bl)) plot_df$highlight[plot_df$gene %in% genes_bl] <- "black"

p_vol <- ggplot(plot_df, aes(x = F, y = neglog10p, color = highlight)) +
  geom_point(alpha = 0.4, size = 0.8) +
  scale_color_manual(values = c(
    "other"       = "grey70",
    "16genes"     = "#C44E52",
    "greenyellow" = "#4C72B0",
    "black"       = "#55A868")) +
  theme_bw(base_size = 12) +
  labs(title = "GSE29801 poly3_F (Retina, age>=60)",
       x = "F statistic (poly3)", y = "-log10(p)", color = "gene set")
save_fig(p_vol, DIR_THIS, "poly3_F_volcano", width = 8, height = 5)

if (nrow(direction_tbl) > 0) {
  plot_dir <- direction_tbl
  plot_dir$contrast <- factor(plot_dir$contrast,
                              levels = c("g2_vs_g1", "g3_vs_g2", "g4_vs_g3"))
  p_dir <- ggplot(plot_dir, aes(x = contrast, y = mean_logFC,
                                fill = gene_set)) +
    geom_col(position = position_dodge(0.8)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_bw(base_size = 12) +
    labs(title = "GSE29801 gene set mean logFC",
         x = "", y = "Mean logFC", fill = "Gene set")
  save_fig(p_dir, DIR_THIS, "gene_set_distribution", width = 8, height = 4)
}

# ============================================================
# 11. summary
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("=== 07_2 GSE29801 复现汇总 (v5) ===\n")
cat(strrep("=", 70), "\n", sep = "")

cat("\n--- 1. 样本（v5 修正映射）---\n")
cat(sprintf("总样本: %d (Retina, age>=%d)\n", nrow(meta_filt), AGE_MIN))
print(table(meta_filt$grade))
print(table(meta_filt$grade_raw))
print(table(meta_filt$grade, meta_filt$region))

cat("\n--- 2. 基因 ---\n")
cat(sprintf("表达基因: %d\n", nrow(expr_final)))

cat("\n--- 3. limma multi ---\n")
for (ct in colnames(cm)) {
  n <- sum(multi_res$FDR[multi_res$contrast == ct] < 0.05, na.rm = TRUE)
  cat(sprintf("  %-10s FDR<0.05: %d\n", ct, n))
}

cat("\n--- 4. poly3_F ---\n")
cat(sprintf("FDR<0.05: %d\n", sum(polyF$FDR < 0.05, na.rm = TRUE)))
cat(sprintf("p<0.001:  %d\n", sum(polyF$p < 0.001, na.rm = TRUE)))

cat("\n--- 5. 基因集方向 ---\n")
print(direction_tbl, row.names = FALSE)

cat("\n--- 6. Permutation ---\n")
print(perm_tbl, row.names = FALSE)

if (nrow(consistency_tbl) > 0) {
  cat("\n--- 7. 与 GSE115828 一致性 ---\n")
  print(consistency_tbl, row.names = FALSE)
}

cat("\n--- 输出位置 ---\n"); cat("目录:", DIR_THIS, "\n")
log_info("=== 07_2_gse29801_replication.r v5 done ===")

