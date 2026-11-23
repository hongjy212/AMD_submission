#!/usr/bin/env Rscript
# =============================================================================
# 07_4_UGGT1_coexpression.r (v5)
#
#
#
#
#
#       - V3a donor-level
#       - V3b pseudobulk sensitivity
#
# Output：results/07_4_UGGT1_coexpression/
# =============================================================================

rm(list = ls()); gc()
setwd("E:/retina_tools/aging_rgc_axon/AMD")
set.seed(2026)

suppressPackageStartupMessages({
  library(Matrix); library(data.table); library(dplyr)
  library(ggplot2); library(logger); library(tidyr)
})

out_dir <- "results/07_4_UGGT1_coexpression"
tab_dir <- file.path(out_dir, "tables")
fig_dir <- file.path(out_dir, "figures")
for (d in c(out_dir, tab_dir, fig_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

log_file <- file.path(out_dir, "07_4_UGGT1_coexpression.log")
log_appender(appender_tee(log_file))
log_layout(layout_glue_generator(format = '{time} [{level}] {msg}'))
log_info("==================== 07_4 UGGT1 co-expression analysis v5 ====================")
log_info("R version: {R.version.string}")

# ============================================================
# ============================================================
log_info("========== 1. 读入数据 ==========")

genes <- fread("E:/retina_tools/raw/GSE221042/GSE221042_genes.csv.gz",
               header = FALSE, col.names = c("idx", "gene_full"))
genes$symbol <- sub(" \\(.*$", "", genes$gene_full)

anno <- fread("E:/retina_tools/raw/GSE221042/GSE221042_sample_annotations.csv.gz",
              header = TRUE)
colnames(anno) <- c("cell_id", "donor", "condition", "n_genes", "total_counts")

counts <- readMM("E:/retina_tools/raw/GSE221042/GSE221042_data.mtx.gz")
rownames(counts) <- anno$cell_id
colnames(counts) <- genes$symbol
counts <- as(counts, "CsparseMatrix")

total_counts <- Matrix::rowSums(counts)
n_genes <- Matrix::rowSums(counts > 0)
mt_genes <- grep("^MT-", colnames(counts), value = TRUE)
mt_counts <- Matrix::rowSums(counts[, mt_genes, drop = FALSE])
mt_pct <- mt_counts / total_counts * 100

keep <- n_genes >= 200 & n_genes <= 7500 & mt_pct <= 10
counts_qc <- counts[keep, ]
anno_qc <- as.data.frame(anno[keep, ])
log_info("QC 后: {nrow(counts_qc)} 细胞")

celltype_markers <- list(
  Rod = c("RHO", "NRL", "NR2E3", "GNAT1", "PDE6B"),
  Cone = c("ARR3", "OPN1SW", "OPN1MW", "PDE6H", "GNB3"),
  HC = c("LHX1", "CALB1", "ONECUT1", "ONECUT2"),
  BC = c("VSX2", "GRM6", "TRPM1", "ISL1"),
  AC = c("GAD1", "GAD2", "SLC6A9", "TFAP2A"),
  RGC = c("RBPMS", "SLC17A6", "NEFL", "SNCG"),
  MG = c("RLBP1", "SLC1A3", "GLUL", "CLU"),
  RPE = c("RPE65", "BEST1", "TTR", "RGR"),
  Astrocyte = c("GFAP", "S100B", "AQP4"),
  Microglia = c("P2RY12", "TMEM119", "CX3CR1", "AIF1")
)

for (ct in names(celltype_markers)) {
  markers <- intersect(celltype_markers[[ct]], colnames(counts_qc))
  if (length(markers) >= 2) {
    anno_qc[[paste0(ct, "_score")]] <-
      Matrix::rowSums(counts_qc[, markers, drop = FALSE] > 0)
  }
}
score_cols <- grep("_score$", colnames(anno_qc), value = TRUE)
score_mat <- as.matrix(anno_qc[, score_cols])
anno_qc$celltype <- gsub("_score$", "", score_cols[apply(score_mat, 1, which.max)])

log_info("细胞类型分布:")
print(table(anno_qc$celltype))

anno_qc$UGGT1 <- as.numeric(counts_qc[, "UGGT1"])
anno_qc$total_umi <- Matrix::rowSums(counts_qc)

# ============================================================
# ============================================================
GENE_SETS <- readRDS("results/01_2_prepare_gene_sets/data/gene_sets.rds")

er_genes     <- toupper(GENE_SETS[["ER_stress"]]$genes)
ci_genes     <- toupper(GENE_SETS[["Mito_CI_Subunits"]]$genes)
oxphos_genes <- toupper(GENE_SETS[["OXPHOS"]]$genes)
comp_genes   <- toupper(GENE_SETS[["Complement"]]$genes)

er_in     <- intersect(er_genes, colnames(counts_qc))
ci_in     <- intersect(ci_genes, colnames(counts_qc))
oxphos_in <- intersect(oxphos_genes, colnames(counts_qc))
comp_in   <- intersect(comp_genes, colnames(counts_qc))

genes16 <- toupper(read.csv("results/02_1_differential_expression_main/tables/sig16_genes.csv")$gene)
genes16_in <- intersect(genes16, colnames(counts_qc))

genes16_no_uggt1 <- setdiff(genes16_in, "UGGT1")
log_info("16 基因签名去 UGGT1 后: {length(genes16_no_uggt1)} genes")

log_info("ER_stress: {length(er_in)} / {length(er_genes)}")
log_info("Complex I: {length(ci_in)} / {length(ci_genes)}")
log_info("OXPHOS: {length(oxphos_in)} / {length(oxphos_genes)}")
log_info("Complement: {length(comp_in)} / {length(comp_genes)}")
log_info("16 基因签名: {length(genes16_in)} / {length(genes16)}")

# ============================================================
# ============================================================
log_info("========== V1: UGGT1 各细胞类型表达 ==========")

ct_summary <- anno_qc %>%
  group_by(celltype) %>%
  summarise(
    mean_expr = mean(UGGT1),
    pct_expressing = mean(UGGT1 > 0) * 100,
    n_cells = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_expr))

log_info("UGGT1 各细胞类型表达:")
print(ct_summary, row.names = FALSE)
write.csv(ct_summary, file.path(tab_dir, "UGGT1_by_celltype.csv"),
          row.names = FALSE)

overall_pct <- mean(anno_qc$UGGT1 > 0) * 100
log_info("UGGT1 全局表达率: {round(overall_pct, 2)}%")

p1 <- ggplot(ct_summary, aes(x = reorder(celltype, -mean_expr),
                               y = mean_expr, fill = pct_expressing)) +
  geom_col() +
  scale_fill_viridis_c(name = "% expressing") +
  labs(title = "UGGT1 expression across retinal cell types",
       subtitle = sprintf("Overall detection rate: %.2f%% of cells", overall_pct),
       x = "Cell type", y = "Mean expression") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10))
ggsave(file.path(fig_dir, "UGGT1_by_celltype.pdf"), p1, width = 10, height = 6)

# ============================================================
# ============================================================
log_info("========== V2: UGGT1 表达分层下下游基因集差异 ==========")

er_score     <- Matrix::rowMeans(counts_qc[, er_in, drop = FALSE])
ci_score     <- Matrix::rowMeans(counts_qc[, ci_in, drop = FALSE])
oxphos_score <- Matrix::rowMeans(counts_qc[, oxphos_in, drop = FALSE])
comp_score   <- Matrix::rowMeans(counts_qc[, comp_in, drop = FALSE])
sig16_score        <- Matrix::rowMeans(counts_qc[, genes16_in, drop = FALSE])
sig16_no_uggt1_score <- Matrix::rowMeans(counts_qc[, genes16_no_uggt1, drop = FALSE])

gene_set_scores <- list(
  ER_stress             = er_score,
  Complex_I             = ci_score,
  OXPHOS                = oxphos_score,
  Complement            = comp_score,
  Signature16           = sig16_score,
  Signature16_no_UGGT1  = sig16_no_uggt1_score
)

uggt1_expr <- anno_qc$UGGT1
total_umi  <- anno_qc$total_umi

has_expr <- uggt1_expr > 0
log_info("UGGT1 表达细胞数: {sum(has_expr)} / {length(has_expr)}")

uggt1_has <- uggt1_expr[has_expr]
uggt1_high <- uggt1_has > median(uggt1_has)

compare_raw <- function(score_vec, name) {
  s_has <- score_vec[has_expr]
  wt <- wilcox.test(s_has[uggt1_high], s_has[!uggt1_high])
  data.frame(
    gene_set = name,
    mean_high = mean(s_has[uggt1_high]),
    mean_low  = mean(s_has[!uggt1_high]),
    log2_ratio = log2(mean(s_has[uggt1_high]) / mean(s_has[!uggt1_high])),
    wilcox_p = wt$p.value,
    stringsAsFactors = FALSE
  )
}

compare_raw_tbl <- do.call(rbind, lapply(names(gene_set_scores), function(nm) {
  compare_raw(gene_set_scores[[nm]], nm)
}))
compare_raw_tbl <- compare_raw_tbl[order(-abs(compare_raw_tbl$log2_ratio)), ]
write.csv(compare_raw_tbl,
          file.path(tab_dir, "UGGT1_highlow_downstream_raw.csv"),
          row.names = FALSE)

log_info("V2a 原始比较（未校正）:")
print(compare_raw_tbl, row.names = FALSE, digits = 4)

log_info("V2b UMI 分层校正（raw mean）...")

umi_bins <- cut(total_umi, breaks = quantile(total_umi,
                                              probs = seq(0, 1, 0.1)),
                labels = FALSE, include.lowest = TRUE)
anno_qc$umi_bin <- umi_bins

compare_umi_adj <- function(score_vec, name) {
  results <- data.frame()
  for (b in unique(umi_bins[has_expr])) {
    idx <- has_expr & umi_bins == b
    if (sum(idx) < 50) next
    uggt1_b <- uggt1_expr[idx]
    high_b <- uggt1_b > median(uggt1_b)
    if (sum(high_b) < 20 || sum(!high_b) < 20) next

    s_b <- score_vec[idx]
    if (mean(s_b[!high_b]) < 1e-9) next

    wt <- wilcox.test(s_b[high_b], s_b[!high_b])
    results <- rbind(results, data.frame(
      umi_bin = b,
      n_cells = sum(idx),
      log2_ratio = log2(mean(s_b[high_b]) / mean(s_b[!high_b])),
      p = wt$p.value,
      stringsAsFactors = FALSE))
  }
  if (nrow(results) == 0) {
    return(data.frame(gene_set = name, n_bins = 0,
                       mean_log2 = NA, median_log2 = NA,
                       n_sig_bins = 0, mean_p = NA))
  }
  data.frame(
    gene_set = name,
    n_bins = nrow(results),
    mean_log2 = mean(results$log2_ratio, na.rm = TRUE),
    median_log2 = median(results$log2_ratio, na.rm = TRUE),
    n_sig_bins = sum(results$p < 0.05 / nrow(results), na.rm = TRUE),
    mean_p = mean(results$p, na.rm = TRUE)
  )
}

compare_umi_tbl <- do.call(rbind, lapply(names(gene_set_scores), function(nm) {
  compare_umi_adj(gene_set_scores[[nm]], nm)
}))
compare_umi_tbl <- compare_umi_tbl[order(-abs(compare_umi_tbl$mean_log2)), ]
write.csv(compare_umi_tbl,
          file.path(tab_dir, "UGGT1_highlow_downstream_umi_adj.csv"),
          row.names = FALSE)

log_info("V2b UMI 分层校正后（raw mean）:")
print(compare_umi_tbl, row.names = FALSE, digits = 4)

log_info("V2c rank-based 分数...")

rank_signature_score <- function(mat, genes) {
  genes <- intersect(genes, colnames(mat))
  if (length(genes) < 3) return(rep(NA_real_, nrow(mat)))
  sub <- as.matrix(mat[, genes, drop = FALSE])
  ranks_mat <- apply(sub, 2, function(x) rank(x, ties.method = "average"))
  if (is.null(dim(ranks_mat))) ranks_mat <- matrix(ranks_mat, ncol = 1)
  score <- 1 - (rowMeans(ranks_mat) / nrow(sub))
  score
}

er_rank     <- rank_signature_score(counts_qc, er_in)
ci_rank     <- rank_signature_score(counts_qc, ci_in)
oxphos_rank <- rank_signature_score(counts_qc, oxphos_in)
comp_rank   <- rank_signature_score(counts_qc, comp_in)
sig16_rank        <- rank_signature_score(counts_qc, genes16_in)
sig16_no_uggt1_rank <- rank_signature_score(counts_qc, genes16_no_uggt1)

umi_cor_rank <- cor(total_umi, er_rank, method = "spearman")
log_info("rank score vs total UMI: Spearman r = {round(umi_cor_rank, 4)}")
if (abs(umi_cor_rank) > 0.3) {
  log_warn("rank score 与 UMI 强相关（|r| > 0.3），在稀疏数据里并非 UMI-independent")
  log_warn("→ V2c 结果需要谨慎解读，以 V2d（配对匹配）为准")
}

gene_set_scores_rank <- list(
  ER_stress             = er_rank,
  Complex_I             = ci_rank,
  OXPHOS                = oxphos_rank,
  Complement            = comp_rank,
  Signature16           = sig16_rank,
  Signature16_no_UGGT1  = sig16_no_uggt1_rank
)

compare_rank <- function(score_vec, name) {
  s_has <- score_vec[has_expr]
  wt <- wilcox.test(s_has[uggt1_high], s_has[!uggt1_high])
  data.frame(
    gene_set = name,
    mean_high = mean(s_has[uggt1_high], na.rm = TRUE),
    mean_low  = mean(s_has[!uggt1_high], na.rm = TRUE),
    log2_ratio = log2(mean(s_has[uggt1_high], na.rm = TRUE) /
                     mean(s_has[!uggt1_high], na.rm = TRUE)),
    wilcox_p = wt$p.value,
    stringsAsFactors = FALSE
  )
}

compare_rank_tbl <- do.call(rbind, lapply(names(gene_set_scores_rank), function(nm) {
  compare_rank(gene_set_scores_rank[[nm]], nm)
}))
compare_rank_tbl <- compare_rank_tbl[order(-abs(compare_rank_tbl$log2_ratio)), ]
write.csv(compare_rank_tbl,
          file.path(tab_dir, "UGGT1_highlow_rankscore.csv"),
          row.names = FALSE)

log_info("V2c rank-based 分数（含 UMI 混杂）:")
print(compare_rank_tbl, row.names = FALSE, digits = 4)

log_info("V2d UMI-matched pairing...")


n_umi_layers <- 20
umi_layer <- cut(total_umi,
                  breaks = quantile(total_umi, probs = seq(0, 1, length.out = n_umi_layers + 1)),
                  labels = FALSE, include.lowest = TRUE)

pair_df <- data.frame(high_idx = integer(), low_idx = integer())

set.seed(2026)
for (b in unique(umi_layer[has_expr])) {
  idx_b <- which(has_expr & umi_layer == b)
  if (length(idx_b) < 40) next

  uggt1_b <- uggt1_expr[idx_b]
  med_b <- median(uggt1_b)

  high_b <- idx_b[uggt1_b >  med_b]
  low_b  <- idx_b[uggt1_b <= med_b]

  if (length(high_b) < 20 || length(low_b) < 20) next

  high_b <- high_b[order(total_umi[high_b])]
  low_b  <- low_b[order(total_umi[low_b])]

  n_pair <- min(length(high_b), length(low_b))

  if (length(high_b) > n_pair) high_b <- sample(high_b, n_pair)
  if (length(low_b)  > n_pair) low_b  <- sample(low_b,  n_pair)

  high_b <- high_b[order(total_umi[high_b])]
  low_b  <- low_b[order(total_umi[low_b])]

  pair_df <- rbind(pair_df, data.frame(
    high_idx = high_b,
    low_idx  = low_b,
    umi_layer = b,
    umi_high = total_umi[high_b],
    umi_low  = total_umi[low_b]
  ))
}

log_info("UMI-matched pairs: {nrow(pair_df)} 对")

if (nrow(pair_df) >= 100) {
  pair_df$umi_diff <- abs(pair_df$umi_high - pair_df$umi_low)
  pair_df$umi_diff_pct <- pair_df$umi_diff / ((pair_df$umi_high + pair_df$umi_low) / 2) * 100

  log_info("配对后 UMI 差异（绝对值）: median={round(median(pair_df$umi_diff),1)}, IQR=[{round(quantile(pair_df$umi_diff,0.25),1)}, {round(quantile(pair_df$umi_diff,0.75),1)}]")
  log_info("配对后 UMI 差异（相对%）: median={round(median(pair_df$umi_diff_pct),2)}%")
  log_info("配对后 UMI 匹配对照 Spearman: r={round(cor(pair_df$umi_high, pair_df$umi_low, method='spearman'),4)}（应接近 1）")

  write.csv(pair_df,
            file.path(tab_dir, "UGGT1_UMI_matched_pairs.csv"),
            row.names = FALSE)

  compare_matched <- function(score_vec, name) {
    s_high <- score_vec[pair_df$high_idx]
    s_low  <- score_vec[pair_df$low_idx]
    d <- s_high - s_low

    # paired Wilcoxon signed-rank
    wt <- wilcox.test(s_high, s_low, paired = TRUE)

    mean_high <- mean(s_high, na.rm = TRUE)
    mean_low  <- mean(s_low,  na.rm = TRUE)
    mean_diff <- mean(d, na.rm = TRUE)

    data.frame(
      gene_set = name,
      n_pairs  = nrow(pair_df),
      mean_high = mean_high,
      mean_low  = mean_low,
      mean_diff = mean_diff,
      log2_ratio = if (mean_low > 0) log2(mean_high / mean_low) else NA_real_,
      wilcox_p  = wt$p.value,
      stringsAsFactors = FALSE
    )
  }

  compare_matched_tbl <- do.call(rbind,
    lapply(names(gene_set_scores), function(nm) {
      compare_matched(gene_set_scores[[nm]], nm)
    }))
  compare_matched_tbl <- compare_matched_tbl[order(-abs(compare_matched_tbl$log2_ratio)), ]
  write.csv(compare_matched_tbl,
            file.path(tab_dir, "UGGT1_highlow_UMI_matched.csv"),
            row.names = FALSE)

  log_info("V2d UMI-matched pairing（主结果）:")
  print(compare_matched_tbl, row.names = FALSE, digits = 4)

  compare_all <- bind_rows(
    data.frame(method = "V2a raw mean",          compare_raw_tbl[, c("gene_set","log2_ratio")]),
    data.frame(method = "V2b raw + UMI strat",   compare_umi_tbl[, c("gene_set","mean_log2")] %>% rename(log2_ratio = mean_log2)),
    data.frame(method = "V2c rank score",        compare_rank_tbl[, c("gene_set","log2_ratio")]),
    data.frame(method = "V2d UMI-matched",       compare_matched_tbl[, c("gene_set","log2_ratio")])
  )
  write.csv(compare_all,
            file.path(tab_dir, "UGGT1_method_comparison.csv"),
            row.names = FALSE)

  p_compare <- ggplot(compare_all,
                       aes(x = gene_set, y = log2_ratio, fill = method)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    coord_flip() +
    labs(title = "UGGT1 high/low: 4 methods, 4 directions?",
         subtitle = "Directional consistency across methods = robustness",
         x = "", y = "log2 ratio (high / low)") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 9),
          legend.position = "bottom")
  ggsave(file.path(fig_dir, "UGGT1_method_comparison.pdf"), p_compare,
         width = 10, height = 6)

  p_umi_match <- ggplot(pair_df, aes(x = umi_high, y = umi_low)) +
    geom_point(alpha = 0.1, size = 0.5) +
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +
    labs(title = "UMI matching quality",
         subtitle = sprintf("Spearman r = %.4f; median |delta| = %.1f UMI",
                             cor(pair_df$umi_high, pair_df$umi_low, method = "spearman"),
                             median(pair_df$umi_diff)),
         x = "UGGT1-high cell total UMI",
         y = "Matched UGGT1-low cell total UMI") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 9))
  ggsave(file.path(fig_dir, "UGGT1_UMI_match_quality.pdf"), p_umi_match,
         width = 6, height = 6)
} else {
  log_warn("UMI-matched pairs 数量不足（{nrow(pair_df)} < 100），跳过 V2d")
}

layer_results <- data.frame()
for (ct in unique(anno_qc$celltype)) {
  idx <- anno_qc$celltype == ct & has_expr
  if (sum(idx) < 50) next
  uggt1_ct <- uggt1_expr[idx]
  high_ct <- uggt1_ct > median(uggt1_ct)
  if (sum(high_ct) < 20 || sum(!high_ct) < 20) next

  er_ct <- er_score[idx]
  ci_ct <- ci_score[idx]
  ox_ct <- oxphos_score[idx]

  wt_er <- wilcox.test(er_ct[high_ct], er_ct[!high_ct])
  wt_ci <- wilcox.test(ci_ct[high_ct], ci_ct[!high_ct])
  wt_ox <- wilcox.test(ox_ct[high_ct], ox_ct[!high_ct])

  layer_results <- rbind(layer_results, data.frame(
    celltype = ct, n_cells = sum(idx),
    log2_ER = log2(mean(er_ct[high_ct]) / mean(er_ct[!high_ct])),
    log2_CI = log2(mean(ci_ct[high_ct]) / mean(ci_ct[!high_ct])),
    log2_OX = log2(mean(ox_ct[high_ct]) / mean(ox_ct[!high_ct])),
    p_ER = wt_er$p.value, p_CI = wt_ci$p.value, p_OX = wt_ox$p.value,
    stringsAsFactors = FALSE
  ))
}
if (nrow(layer_results) > 0) {
  layer_results <- layer_results[order(-abs(layer_results$log2_ER)), ]
  write.csv(layer_results,
            file.path(tab_dir, "UGGT1_highlow_by_celltype.csv"),
            row.names = FALSE)
  log_info("按细胞类型分层（raw mean）:")
  print(layer_results, row.names = FALSE, digits = 4)
}

# ============================================================
# ============================================================
log_info("========== V3: rank-based signature scores ==========")

anno_qc$UCell_ER  <- er_rank
anno_qc$UCell_CI  <- ci_rank
anno_qc$UCell_OX  <- oxphos_rank
anno_qc$UCell_G16 <- sig16_rank
anno_qc$UCell_G16_noUGGT1 <- sig16_no_uggt1_rank

log_info("评分摘要:")
print(summary(anno_qc[, c("UCell_ER", "UCell_CI", "UCell_OX", "UCell_G16")]))

donor_agg <- anno_qc %>%
  filter(condition %in% c("Healthy", "dAMD")) %>%
  group_by(donor, condition) %>%
  summarise(
    mean_UGGT1     = mean(UGGT1),
    mean_UCell_ER  = mean(UCell_ER, na.rm = TRUE),
    mean_UCell_CI  = mean(UCell_CI, na.rm = TRUE),
    mean_UCell_OX  = mean(UCell_OX, na.rm = TRUE),
    mean_UCell_G16 = mean(UCell_G16, na.rm = TRUE),
    mean_UCell_G16_noUGGT1 = mean(UCell_G16_noUGGT1, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop")

write.csv(donor_agg, file.path(tab_dir, "donor_UGGT1_scores.csv"),
          row.names = FALSE)
log_info("donor-level 评分:")
print(as.data.frame(donor_agg), row.names = FALSE, digits = 4)

# ---- V3b: dAMD vs Healthy ----
log_info("dAMD vs Healthy 对比（donor-level）:")
for (metric in c("mean_UGGT1", "mean_UCell_ER", "mean_UCell_CI",
                  "mean_UCell_OX", "mean_UCell_G16", "mean_UCell_G16_noUGGT1")) {
  h <- donor_agg[[metric]][donor_agg$condition == "Healthy"]
  d <- donor_agg[[metric]][donor_agg$condition == "dAMD"]
  if (length(h) >= 3 && length(d) >= 3 && mean(h) > 0) {
    tt <- t.test(d, h)
    log2fc <- log2(mean(d) / mean(h))
    log_info("  {metric}: Healthy={round(mean(h),4)}, dAMD={round(mean(d),4)}, log2FC={round(log2fc,3)}, p={format(tt$p.value, scientific=TRUE)}")
  }
}

log_info("--- V3c: pseudobulk 敏感性分析 ---")

donor_dh_agg <- donor_agg %>%
  filter(condition %in% c("Healthy", "dAMD"))

donor_dh_agg$weighted_ER  <- donor_dh_agg$mean_UCell_ER  * donor_dh_agg$n_cells
donor_dh_agg$weighted_CI  <- donor_dh_agg$mean_UCell_CI  * donor_dh_agg$n_cells
donor_dh_agg$weighted_OX  <- donor_dh_agg$mean_UCell_OX  * donor_dh_agg$n_cells
donor_dh_agg$weighted_G16 <- donor_dh_agg$mean_UCell_G16 * donor_dh_agg$n_cells

write.csv(donor_dh_agg,
          file.path(tab_dir, "pseudobulk_donor_scores.csv"),
          row.names = FALSE)

pb_results <- data.frame()
for (metric_base in c("mean_UCell_ER", "mean_UCell_CI",
                       "mean_UCell_OX", "mean_UCell_G16")) {
  h_uw <- donor_dh_agg[[metric_base]][donor_dh_agg$condition == "Healthy"]
  d_uw <- donor_dh_agg[[metric_base]][donor_dh_agg$condition == "dAMD"]
  if (length(h_uw) >= 3 && length(d_uw) >= 3) {
    tt_uw <- t.test(d_uw, h_uw)
    pb_results <- rbind(pb_results, data.frame(
      metric = metric_base,
      method = "donor_mean_unweighted",
      mean_healthy = mean(h_uw),
      mean_damd = mean(d_uw),
      log2FC = log2(mean(d_uw) / mean(h_uw)),
      p = tt_uw$p.value,
      stringsAsFactors = FALSE))
  }

  col_name <- sub("mean_UCell_", "UCell_", metric_base)
  pooled_h <- mean(anno_qc[[col_name]][anno_qc$condition == "Healthy"], na.rm = TRUE)
  pooled_d <- mean(anno_qc[[col_name]][anno_qc$condition == "dAMD"], na.rm = TRUE)
  pb_results <- rbind(pb_results, data.frame(
    metric = metric_base,
    method = "pooled_all_cells",
    mean_healthy = pooled_h,
    mean_damd = pooled_d,
    log2FC = log2(pooled_d / pooled_h),
    p = NA_real_,
    stringsAsFactors = FALSE))
}

write.csv(pb_results,
          file.path(tab_dir, "pseudobulk_sensitivity.csv"),
          row.names = FALSE)

log_info("pseudobulk 敏感性对比:")
print(pb_results, row.names = FALSE, digits = 4)

donor_ct <- anno_qc %>%
  filter(condition %in% c("Healthy", "dAMD")) %>%
  group_by(donor, condition, celltype) %>%
  summarise(
    mean_UGGT1 = mean(UGGT1),
    mean_UCell_ER = mean(UCell_ER, na.rm = TRUE),
    mean_UCell_CI = mean(UCell_CI, na.rm = TRUE),
    mean_UCell_OX = mean(UCell_OX, na.rm = TRUE),
    mean_UCell_G16 = mean(UCell_G16, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop") %>%
  filter(n_cells >= 20)

ct_results <- data.frame()
for (ct in unique(donor_ct$celltype)) {
  for (metric in c("mean_UGGT1", "mean_UCell_ER", "mean_UCell_CI",
                    "mean_UCell_OX", "mean_UCell_G16")) {
    sub <- donor_ct[donor_ct$celltype == ct, ]
    h <- sub[[metric]][sub$condition == "Healthy"]
    d <- sub[[metric]][sub$condition == "dAMD"]
    if (length(h) >= 3 && length(d) >= 3 && mean(h) > 0) {
      tt <- t.test(d, h)
      ct_results <- rbind(ct_results, data.frame(
        celltype = ct, metric = metric,
        mean_healthy = mean(h), mean_damd = mean(d),
        log2FC = log2(mean(d) / mean(h)),
        p = tt$p.value, stringsAsFactors = FALSE))
    }
  }
}
if (nrow(ct_results) > 0) {
  ct_results <- ct_results[order(-abs(ct_results$log2FC)), ]
  write.csv(ct_results,
            file.path(tab_dir, "UGGT1_disease_by_celltype.csv"),
            row.names = FALSE)
  log_info("按细胞类型分层（top 15 效应量）:")
  print(head(ct_results, 15), row.names = FALSE, digits = 4)
}

# ---- V3e: visualization ----
plot_df <- anno_qc %>%
  filter(condition %in% c("Healthy", "dAMD")) %>%
  select(donor, condition, celltype, UCell_ER, UCell_CI, UCell_OX, UCell_G16) %>%
  pivot_longer(cols = c("UCell_ER", "UCell_CI", "UCell_OX", "UCell_G16"),
               names_to = "metric", values_to = "score") %>%
  filter(!is.na(score))

p2 <- ggplot(plot_df, aes(x = condition, y = score, fill = condition)) +
  geom_violin(alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.2, alpha = 0.5, outlier.size = 0.2) +
  facet_grid(metric ~ celltype, scales = "free_y") +
  labs(title = "Signature scores by condition and cell type",
       x = "Condition", y = "Score") +
  theme_minimal(base_size = 9) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")
ggsave(file.path(fig_dir, "UCell_by_condition_celltype.pdf"), p2,
       width = 16, height = 10)

# ============================================================
# ============================================================
log_info("========== V4: UGGT1 与 16 基因签名共表达 ==========")

other16 <- setdiff(genes16_in, "UGGT1")
cor_results <- data.frame()
uggt1_has_vec <- uggt1_expr[has_expr]
for (g in other16) {
  g_vec <- as.numeric(counts_qc[has_expr, g])
  if (sd(g_vec) == 0 || sd(uggt1_has_vec) == 0) next
  r <- cor(uggt1_has_vec, g_vec, method = "spearman")
  cor_results <- rbind(cor_results, data.frame(
    gene = g, spearman_r = r, stringsAsFactors = FALSE))
}
if (nrow(cor_results) > 0) {
  cor_results <- cor_results[order(-abs(cor_results$spearman_r)), ]
  write.csv(cor_results,
            file.path(tab_dir, "UGGT1_16gene_correlation.csv"),
            row.names = FALSE)
  log_info("UGGT1 与 16 基因签名其他成员的 Spearman 相关:")
  print(head(cor_results, 10), row.names = FALSE, digits = 4)
}

# ============================================================
# 7. summary
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("=== 07_4 UGGT1 co-expression v5 汇总 ===\n")
cat(strrep("=", 70), "\n", sep = "")

cat("\nV1: UGGT1 表达率 =", round(overall_pct, 2), "%\n")

cat("\nV2a: raw mean（受 UMI 混杂）\n")
print(compare_raw_tbl, row.names = FALSE, digits = 4)

cat("\nV2b: raw + UMI 分层校正（可能过度校正）\n")
print(compare_umi_tbl, row.names = FALSE, digits = 4)

cat("\nV2c: rank score（rank 也 UMI 相关，见日志诊断）\n")
print(compare_rank_tbl, row.names = FALSE, digits = 4)

if (exists("compare_matched_tbl")) {
  cat("\nV2d: UMI-matched pairing（★ v5 主结果）\n")
  print(compare_matched_tbl, row.names = FALSE, digits = 4)
}

cat("\nV3: dAMD vs Healthy（donor-level）\n")
print(as.data.frame(donor_agg), row.names = FALSE, digits = 4)

cat("\nV3c: pseudobulk 敏感性\n")
print(pb_results, row.names = FALSE, digits = 4)

cat("\nV4: UGGT1 与 16 基因签名相关 top 5\n")
if (nrow(cor_results) > 0) print(head(cor_results, 5), row.names = FALSE, digits = 4)

cat("\n目录:", out_dir, "\n")
log_info("=== 07_4 UGGT1 v5 done ===")

