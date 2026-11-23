#!/usr/bin/env Rscript
# ============================================================
# 07_1_atac_validation.r  (v4)
#
#   4. report Cohen's d effect size
#
# Input：
#   E:/retina_tools/raw/GSE99287/GSE99287_RPE_ATACSeq_peak_counts.txt.gz
#   E:/retina_tools/raw/GSE99287/GSE99287_Retina_ATACSeq_peak_counts.txt.gz
#   results/02_1_differential_expression_main/tables/sig16_genes.csv
#   results/05_1_consensus_wgcna/tables/module_genes_greenyellow_symbolfix.csv
#   results/05_1_consensus_wgcna/tables/module_genes_black_symbolfix.csv
#
# Output：results/07_1_atac_validation/
#   data/promo_long_rpe.rds
#   data/promo_long_retina.rds
#   tables/gene_set_ATAC_summary.csv
#   tables/gene_level_<set>_<tissue>.csv
#   tables/key_genes_ATAC.csv
#   figures/ATAC_coef_distribution_<tissue>.pdf
#   figures/ATAC_permutation_null_<tissue>.pdf
#   figures/ATAC_key_genes_<tissue>.pdf
# ============================================================

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(tidyr)
  library(tibble); library(ggplot2)
})

DIR_THIS <- setup_output_dir("07_1_atac_validation")
log_info("=== 07_1_atac_validation.r v4 start ===")

# ============================================================
# 0. Parameters
# ============================================================
ATAC_RAW_DIR  <- "E:/retina_tools/raw/GSE99287"
N_PERM        <- 1000
N_RANDOM_SET  <- 50        # 随机基因集负对照
SEED          <- 42
PROMOTER_PAT  <- "promoter"
KEY_GENES     <- c("NDUFB8", "AQR", "UGGT1")

ANNO_COLS <- c("Chr", "Hg19_start", "Hg19_end", "Type", "Gene",
               "Coefficient_of_stage", "Fdr_of_stage_coefficient",
               "Differentially_accessible_region")


# ============================================================
# ============================================================
log_info("--- 1. 读 GSE99287 peak counts ---")

f_rpe    <- file.path(ATAC_RAW_DIR, "GSE99287_RPE_ATACSeq_peak_counts.txt.gz")
f_retina <- file.path(ATAC_RAW_DIR, "GSE99287_Retina_ATACSeq_peak_counts.txt.gz")

stopifnot(file.exists(f_rpe))
rpe <- fread(f_rpe)
log_info(sprintf("RPE: %d peaks x %d cols", nrow(rpe), ncol(rpe)))

retina <- NULL
if (file.exists(f_retina)) {
  retina <- fread(f_retina)
  log_info(sprintf("Retina: %d peaks x %d cols", nrow(retina), ncol(retina)))
} else {
  log_warn("Retina 文件未找到，仅分析 RPE")
}

stopifnot(all(ANNO_COLS %in% colnames(rpe)))
if (!is.null(retina)) stopifnot(all(ANNO_COLS %in% colnames(retina)))

# ============================================================
# ============================================================
log_info("--- 2. 构建 promoter peak -> gene 长表 ---")

build_promo_long <- function(dt, tissue_name) {
  dt <- as.data.table(dt)
  dt[, is_promoter := grepl(PROMOTER_PAT, Type, fixed = TRUE)]
  promo <- dt[is_promoter == TRUE]
  log_info(sprintf("%s: promoter peak = %d / %d",
                   tissue_name, nrow(promo), nrow(dt)))

  promo[, peak_idx := .I]
  long <- promo[, .(gene = trimws(unlist(strsplit(Gene, "[,;]"))),
                    Chr, Hg19_start, Hg19_end, Type,
                    coef = Coefficient_of_stage,
                    fdr  = Fdr_of_stage_coefficient,
                    dar  = Differentially_accessible_region),
                by = peak_idx]
  long <- long[gene != "" & !is.na(gene)]
  long[, gene := toupper(gene)]

  log_info(sprintf("%s: 长表 %d 行, %d 基因",
                   tissue_name, nrow(long), uniqueN(long$gene)))
  long
}

promo_rpe    <- build_promo_long(rpe, "RPE")
promo_retina <- if (!is.null(retina)) build_promo_long(retina, "Retina") else NULL

save_rds(promo_rpe, DIR_THIS, "promo_long_rpe.rds")
if (!is.null(promo_retina)) save_rds(promo_retina, DIR_THIS, "promo_long_retina.rds")

# ============================================================
# ============================================================
log_info("--- 3. 载入 p3 信号基因集 ---")

# 3.1 16 genes
genes_16 <- toupper(read_upstream("02_1_differential_expression_main", "tables",
                                  "sig16_genes.csv")$gene)
log_info(sprintf("16 基因: %d", length(genes_16)))

# 3.2 WGCNA module
read_module_symbols <- function(orig_path, fix_path) {
  if (file.exists(fix_path)) {
    df <- read.csv(fix_path); sym <- df$SYMBOL
    log_info(sprintf("读取 %s (fix)", basename(fix_path)))
  } else {
    if (!file.exists(orig_path)) return(NULL)
    df <- read.csv(orig_path)
    suppressPackageStartupMessages({ library(org.Hs.eg.db); library(AnnotationDbi) })
    ens <- sub("\\..*$", "", df$ENSEMBL)
    sym <- mapIds(org.Hs.eg.db, keys = ens, column = "SYMBOL",
                  keytype = "ENSEMBL", multiVals = "first")
  }
  sym <- toupper(sym); sym[!is.na(sym) & sym != ""]
}

tbl_dir <- file.path(DIR_RESULT, "05_1_consensus_wgcna", "tables")
genes_gy <- read_module_symbols(
  file.path(tbl_dir, "module_genes_greenyellow.csv"),
  file.path(tbl_dir, "module_genes_greenyellow_symbolfix.csv"))
genes_bl <- read_module_symbols(
  file.path(tbl_dir, "module_genes_black.csv"),
  file.path(tbl_dir, "module_genes_black_symbolfix.csv"))

polyF_path <- file.path(DIR_RESULT, "04_2_pathway_score", "tables",
                        "pathway_poly3F_adjusted.csv")
GENE_SETS <- read_upstream("01_2_prepare_gene_sets", "data", "gene_sets.rds")

p3_only_pathways <- character(0)
if (file.exists(polyF_path)) {
  polyF <- read.csv(polyF_path)
  p3_only_pathways <- polyF$pathway[polyF$pattern == "p3_only"]
  log_info(sprintf("从 04_2 读入 p3_only 通路: %d 个", length(p3_only_pathways)))
} else {
  log_warn("04_2 pathway_poly3F_adjusted.csv 不存在，跳过 19 通路")
}

gene_sets <- list(
  "16genes"     = genes_16,
  "greenyellow" = genes_gy,
  "black"       = genes_bl)

for (pw in p3_only_pathways) {
  if (!pw %in% names(GENE_SETS)) next
  gs_sym <- toupper(GENE_SETS[[pw]]$genes)
  gs_sym <- gs_sym[!is.na(gs_sym) & gs_sym != ""]
  if (length(gs_sym) >= 10) gene_sets[[pw]] <- gs_sym
}

gene_sets <- gene_sets[!sapply(gene_sets, is.null)]
log_info(sprintf("总 gene set 数: %d", length(gene_sets)))

# overwritecheck
for (sname in names(gene_sets)) {
  n_in <- sum(gene_sets[[sname]] %in% unique(promo_rpe$gene))
  log_info(sprintf("%s: %d/%d 在 RPE promoter 表中",
                   sname, n_in, length(gene_sets[[sname]])))
}

# ============================================================
# ============================================================
log_info("--- 4. 随机基因集负对照 ---")

set.seed(SEED)
all_genes_in_promo <- unique(promo_rpe$gene)

random_sets <- list()
for (i in seq_len(N_RANDOM_SET)) {
  n <- sample(20:150, 1)  # 随机大小，模拟真实基因集
  random_sets[[paste0("random_", i)]] <- sample(all_genes_in_promo, n)
}

# ============================================================
# ============================================================
log_info("--- 5. 开始 tissue 级分析 ---")

process_tissue <- function(promo_long, tissue_name, gene_sets, key_genes) {

  gene_level <- promo_long[, .(
    mean_coef = mean(coef, na.rm = TRUE),
    n_peaks   = .N,
    n_dar     = sum(dar == 1, na.rm = TRUE)
  ), by = gene]

  log_info(sprintf("[%s] gene-level: %d 基因, mean_coef 中位数=%.3f",
                   tissue_name, nrow(gene_level), median(gene_level$mean_coef)))

  # ---- 5b. gene set summary ----
  summarize_gene_set <- function(gene_set, set_name) {
    sub <- gene_level[gene %in% gene_set]
    if (nrow(sub) < 5) return(NULL)

    # Background：allgenes
    bg <- gene_level[!gene %in% unlist(gene_sets)]
    bg_sd <- sd(bg$mean_coef, na.rm = TRUE)
    bg_mean <- mean(bg$mean_coef, na.rm = TRUE)

    # Cohen's d
    cohen_d <- (mean(sub$mean_coef, na.rm = TRUE) - bg_mean) / bg_sd

    data.frame(
      gene_set      = set_name,
      tissue        = tissue_name,
      n_genes       = nrow(sub),
      n_peaks_total = sum(sub$n_peaks),
      mean_coef     = round(mean(sub$mean_coef, na.rm = TRUE), 4),
      median_coef   = round(median(sub$mean_coef, na.rm = TRUE), 4),
      bg_mean_coef  = round(bg_mean, 4),
      cohen_d       = round(cohen_d, 3),
      n_dar_genes   = sum(sub$n_dar > 0),
      pct_dar_gene  = round(100 * mean(sub$n_dar > 0), 2),
      wilcox_p_zero = wilcox.test(sub$mean_coef, mu = 0)$p.value,
      stringsAsFactors = FALSE)
  }

  summary_tbl <- data.frame()
  for (sname in names(gene_sets)) {
    r <- summarize_gene_set(gene_sets[[sname]], sname)
    if (!is.null(r)) summary_tbl <- rbind(summary_tbl, r)
  }

  set.seed(SEED)

  perm_test_peak <- function(gene_set, set_name) {
    sub <- promo_long[gene %in% gene_set]
    if (nrow(sub) < 5) return(NULL)
    obs_mean <- mean(sub$coef, na.rm = TRUE)
    obs_dar  <- mean(sub$dar == 1, na.rm = TRUE)
    n_target <- nrow(sub)

    perm_res <- replicate(N_PERM, {
      idx <- sample(nrow(promo_long), n_target, replace = FALSE)
      s <- promo_long[idx]
      c(mean(s$coef, na.rm = TRUE), mean(s$dar == 1, na.rm = TRUE))
    })
    perm_coef <- perm_res[1, ]
    perm_dar  <- perm_res[2, ]

    data.frame(
      gene_set = set_name, tissue = tissue_name, level = "peak",
      n_units = n_target,
      obs_mean_coef  = round(obs_mean, 4),
      perm_coef_med  = round(median(perm_coef), 4),
      perm_coef_q025 = round(quantile(perm_coef, 0.025), 4),
      perm_coef_q975 = round(quantile(perm_coef, 0.975), 4),
      perm_p_coef_neg = mean(perm_coef <= obs_mean),
      perm_p_coef_pos = mean(perm_coef >= obs_mean),
      perm_p_coef_two = mean(abs(perm_coef) >= abs(obs_mean)),
      obs_pct_dar    = round(100 * obs_dar, 2),
      perm_dar_med   = round(100 * median(perm_dar), 2),
      perm_p_dar     = mean(perm_dar >= obs_dar),
      stringsAsFactors = FALSE)
  }

  perm_test_gene <- function(gene_set, set_name) {
    sub <- gene_level[gene %in% gene_set]
    if (nrow(sub) < 5) return(NULL)
    obs_mean <- mean(sub$mean_coef, na.rm = TRUE)
    obs_dar  <- mean(sub$n_dar > 0)
    n_target <- nrow(sub)

    all_genes <- gene_level$gene
    perm_res <- replicate(N_PERM, {
      g <- sample(all_genes, n_target, replace = FALSE)
      s <- gene_level[gene %in% g]
      c(mean(s$mean_coef, na.rm = TRUE), mean(s$n_dar > 0))
    })
    perm_coef <- perm_res[1, ]
    perm_dar  <- perm_res[2, ]

    data.frame(
      gene_set = set_name, tissue = tissue_name, level = "gene",
      n_units = n_target,
      obs_mean_coef  = round(obs_mean, 4),
      perm_coef_med  = round(median(perm_coef), 4),
      perm_coef_q025 = round(quantile(perm_coef, 0.025), 4),
      perm_coef_q975 = round(quantile(perm_coef, 0.975), 4),
      perm_p_coef_neg = mean(perm_coef <= obs_mean),
      perm_p_coef_pos = mean(perm_coef >= obs_mean),
      perm_p_coef_two = mean(abs(perm_coef) >= abs(obs_mean)),
      obs_pct_dar    = round(100 * obs_dar, 2),
      perm_dar_med   = round(100 * median(perm_dar), 2),
      perm_p_dar     = mean(perm_dar >= obs_dar),
      stringsAsFactors = FALSE)
  }

  perm_tbl <- data.frame()
  for (sname in names(gene_sets)) {
    r1 <- perm_test_peak(gene_sets[[sname]], sname)
    r2 <- perm_test_gene(gene_sets[[sname]], sname)
    if (!is.null(r1)) perm_tbl <- rbind(perm_tbl, r1)
    if (!is.null(r2)) perm_tbl <- rbind(perm_tbl, r2)
  }

  random_coef <- sapply(names(random_sets), function(rs) {
    sub <- gene_level[gene %in% random_sets[[rs]]]
    if (nrow(sub) < 5) return(NA_real_)
    mean(sub$mean_coef, na.rm = TRUE)
  })
  random_coef <- random_coef[!is.na(random_coef)]

  random_baseline <- data.frame(
    tissue = tissue_name,
    n_random = length(random_coef),
    rand_mean = round(mean(random_coef), 4),
    rand_sd   = round(sd(random_coef), 4),
    rand_q025 = round(quantile(random_coef, 0.025), 4),
    rand_q975 = round(quantile(random_coef, 0.975), 4),
    stringsAsFactors = FALSE)

  key_rpe <- promo_long[gene %in% key_genes]
  key_tbl <- key_rpe[, .(
    n_peaks   = .N,
    mean_coef = round(mean(coef, na.rm = TRUE), 4),
    n_dar     = sum(dar == 1, na.rm = TRUE),
    types     = paste(unique(Type), collapse = "|")
  ), by = gene]
  key_tbl$tissue <- tissue_name

  # ---- 5f. figure ----
  plot_list <- list()
  for (sname in names(gene_sets)) {
    sub <- promo_long[gene %in% gene_sets[[sname]]]
    if (nrow(sub) > 0)
      plot_list[[sname]] <- data.frame(coef = sub$coef, gene_set = sname)
  }
  plot_list[["all_promoter"]] <- data.frame(coef = promo_long$coef,
                                            gene_set = "all_promoter")
  plot_df <- bind_rows(plot_list)
  p_dist <- ggplot(plot_df, aes(x = coef, fill = gene_set)) +
    geom_density(alpha = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme_bw(base_size = 12) +
    labs(title = sprintf("%s promoter accessibility stage coefficient", tissue_name),
         x = "Coefficient_of_stage", y = "Density", fill = "Gene set")
  save_fig(p_dist, DIR_THIS,
           sprintf("ATAC_coef_distribution_%s", tissue_name),
           width = 10, height = 6)

  list(summary = summary_tbl, perm = perm_tbl,
       random_baseline = random_baseline,
       key = key_tbl, gene_level = gene_level)
}

# ============================================================
# ============================================================
log_info("--- 6. RPE 分析 ---")
res_rpe <- process_tissue(promo_rpe, "RPE", gene_sets, KEY_GENES)

res_retina <- NULL
if (!is.null(promo_retina)) {
  log_info("--- 6b. Retina 分析 ---")
  res_retina <- process_tissue(promo_retina, "Retina", gene_sets, KEY_GENES)
}

# ============================================================
# ============================================================
summary_all <- rbind(res_rpe$summary,
                     if (!is.null(res_retina)) res_retina$summary else NULL)
perm_all    <- rbind(res_rpe$perm,
                     if (!is.null(res_retina)) res_retina$perm else NULL)
key_all     <- rbind(res_rpe$key,
                     if (!is.null(res_retina)) res_retina$key else NULL)
random_base <- rbind(res_rpe$random_baseline,
                     if (!is.null(res_retina)) res_retina$random_baseline else NULL)

perm_all$FDR_coef_neg <- p.adjust(perm_all$perm_p_coef_neg, method = "BH")
perm_all$FDR_coef_pos <- p.adjust(perm_all$perm_p_coef_pos, method = "BH")
perm_all$FDR_coef_two <- p.adjust(perm_all$perm_p_coef_two, method = "BH")
perm_all$FDR_dar      <- p.adjust(perm_all$perm_p_dar, method = "BH")

perm_all$direction <- with(perm_all,
  ifelse(perm_p_coef_neg < perm_p_coef_pos, "negative", "positive"))
perm_all$significant <- perm_all$FDR_coef_neg < 0.05 & perm_all$direction == "negative"

save_table(summary_all, DIR_THIS, "gene_set_ATAC_summary.csv")
save_table(perm_all,    DIR_THIS, "permutation_results.csv")
save_table(as.data.frame(key_all), DIR_THIS, "key_genes_ATAC.csv")
save_table(random_base, DIR_THIS, "random_null_baseline.csv")

# ============================================================
# ============================================================
log_info("--- 8. ATAC vs RNA p3 一致性 ---")

if (file.exists(polyF_path)) {
  consistency_tbl <- data.frame()
  for (tissue in unique(perm_all$tissue)) {
    perm_t <- perm_all[perm_all$tissue == tissue & perm_all$level == "gene", ]

    perm_t <- perm_t[perm_t$gene_set %in% p3_only_pathways, ]
    if (nrow(perm_t) == 0) next

    merged <- merge(perm_t,
                    polyF[, c("pathway", "p3_coef", "p3", "pattern")],
                    by.x = "gene_set", by.y = "pathway")
    if (nrow(merged) < 3) next

    r <- cor(merged$obs_mean_coef, merged$p3_coef,
             method = "spearman", use = "complete.obs")

    consistency_tbl <- rbind(consistency_tbl, data.frame(
      tissue       = tissue,
      n_pathways   = nrow(merged),
      spearman_r   = round(r, 3),
      n_same_dir   = sum(sign(merged$obs_mean_coef) == sign(merged$p3_coef)),
      pct_same_dir = round(100 * mean(sign(merged$obs_mean_coef) ==
                                       sign(merged$p3_coef)), 1),
      stringsAsFactors = FALSE))

    # scatter
    p_sc <- ggplot(merged, aes(x = p3_coef, y = obs_mean_coef)) +
      geom_point(alpha = 0.6) +
      geom_smooth(method = "lm", se = FALSE, color = "red") +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_vline(xintercept = 0, linetype = "dashed") +
      theme_bw(base_size = 11) +
      labs(title = sprintf("ATAC vs RNA p3: %s (r = %.3f)", tissue, r),
           x = "RNA p3_coef", y = "ATAC mean_coef") +
      geom_text(aes(label = gene_set), size = 2.5, hjust = -0.1, vjust = -0.3)
    save_fig(p_sc, DIR_THIS,
             sprintf("ATAC_vs_RNA_p3_scatter_%s", tissue),
             width = 7, height = 6)
  }

  if (nrow(consistency_tbl) > 0) {
    save_table(consistency_tbl, DIR_THIS, "ATAC_vs_RNA_p3_consistency.csv")
  }
}

# ============================================================
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("=== 07_1 ATAC 验证汇总 (v4) ===\n")
cat(strrep("=", 70), "\n", sep = "")

cat("\n--- 1. 数据 ---\n")
cat(sprintf("RPE:    %d peaks, promoter %d\n", nrow(rpe), nrow(promo_rpe)))
if (!is.null(retina))
  cat(sprintf("Retina: %d peaks, promoter %d\n", nrow(retina), nrow(promo_retina)))

cat("\n--- 2. Random null baseline ---\n")
print(random_base, row.names = FALSE)

cat("\n--- 3. Permutation (含单侧 + FDR) ---\n")
print(perm_all, row.names = FALSE, digits = 4)

cat("\n--- 4. 显著的通路（FDR_coef_neg < 0.05）---\n")
sig <- perm_all[perm_all$significant == TRUE, ]
if (nrow(sig) > 0) {
  print(sig[, c("gene_set", "tissue", "level", "n_units",
                 "obs_mean_coef", "perm_p_coef_neg", "FDR_coef_neg")],
        row.names = FALSE, digits = 4)
} else {
  cat("无\n")
}

cat("\n--- 5. 关键基因 ---\n")
print(as.data.frame(key_all), row.names = FALSE)

if (file.exists(file.path(DIR_THIS, "tables", "ATAC_vs_RNA_p3_consistency.csv"))) {
  cat("\n--- 6. ATAC vs RNA p3 一致性 ---\n")
  print(read.csv(file.path(DIR_THIS, "tables", "ATAC_vs_RNA_p3_consistency.csv")),
        row.names = FALSE)
}

cat("\n--- 输出位置 ---\n"); cat("目录:", DIR_THIS, "\n")
log_info("=== 07_1_atac_validation.r v4 done ===")

