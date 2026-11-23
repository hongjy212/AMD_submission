#!/usr/bin/env Rscript
# =============================================================================
# 07_3_NDUFB8_coexpression.r (v10)
#
#
#
#
#
# Output：results/07_3_NDUFB8_coexpression/
# =============================================================================

rm(list = ls()); gc()
setwd("E:/retina_tools/aging_rgc_axon/AMD")
set.seed(2026)

suppressPackageStartupMessages({
  library(Matrix); library(data.table); library(dplyr)
  library(ggplot2); library(logger); library(tidyr); library(UCell)
})

out_dir <- "results/07_3_NDUFB8_coexpression"
tab_dir <- file.path(out_dir, "tables")
fig_dir <- file.path(out_dir, "figures")
for (d in c(out_dir, tab_dir, fig_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

log_file <- file.path(out_dir, "07_3_NDUFB8_coexpression.log")
log_appender(appender_tee(log_file))
log_layout(layout_glue_generator(format = '{time} [{level}] {msg}'))
log_info("==================== 07_3 NDUFB8 co-expression analysis v10 ====================")
log_info("R version: {R.version.string}")
log_info("UCell version: {as.character(packageVersion('UCell'))}")

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

# ============================================================
# ============================================================
log_info("========== 2. 细胞类型识别 ==========")

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
    ct_score <- Matrix::rowSums(counts_qc[, markers, drop = FALSE] > 0)
    anno_qc[[paste0(ct, "_score")]] <- ct_score
  }
}

score_cols <- grep("_score$", colnames(anno_qc), value = TRUE)
score_mat <- as.matrix(anno_qc[, score_cols])
anno_qc$celltype <- gsub("_score$", "", score_cols[apply(score_mat, 1, which.max)])

log_info("细胞类型分布:")
print(table(anno_qc$celltype))

# ============================================================
# ============================================================
log_info("========== 3. 载入基因集 ==========")

GENE_SETS <- readRDS("results/01_2_prepare_gene_sets/data/gene_sets.rds")

ci_genes <- toupper(GENE_SETS[["Mito_CI_Subunits"]]$genes)
oxphos_genes <- toupper(GENE_SETS[["OXPHOS"]]$genes)
ci_in <- intersect(ci_genes, colnames(counts_qc))
oxphos_in <- intersect(oxphos_genes, colnames(counts_qc))

log_info("Complex I 亚基在数据里: {length(ci_in)} / {length(ci_genes)}")
log_info("OXPHOS 基因在数据里: {length(oxphos_in)} / {length(oxphos_genes)}")

genes16 <- toupper(read.csv("results/02_1_differential_expression_main/tables/sig16_genes.csv")$gene)
genes16_in <- intersect(genes16, colnames(counts_qc))
log_info("16 基因签名在数据里: {length(genes16_in)} / {length(genes16)}")

complex_sets <- list(
  CI   = "Mito_CI_Subunits",
  CII  = "Mito_CII_Subunits",
  CIII = "Mito_CIII_Subunits",
  CIV  = "Mito_CIV_Subunits",
  CV   = "Mito_CV_Subunits"
)
complex_genes <- list()
for (cp in names(complex_sets)) {
  sn <- complex_sets[[cp]]
  if (sn %in% names(GENE_SETS)) {
    g <- toupper(GENE_SETS[[sn]]$genes)
    complex_genes[[cp]] <- intersect(g, colnames(counts_qc))
    log_info("{cp} ({sn}): {length(complex_genes[[cp]])} / {length(g)} genes in data")
  } else {
    log_warn("{cp}: gene set {sn} not found, skip")
  }
}

pick_col <- function(df, base) {
  if (base %in% colnames(df)) return(base)
  alt <- paste0(base, "_UCell")
  if (alt %in% colnames(df)) return(alt)
  stop(sprintf("返回结果里找不到列：%s 或 %s", base, alt))
}

# ============================================================
# ============================================================
log_info("========== V1: NDUFB8 各细胞类型表达 ==========")

if (!"NDUFB8" %in% colnames(counts_qc)) {
  log_warn("NDUFB8 不在数据里")
} else {
  ndufb8_df <- data.frame(
    cell_id = anno_qc$cell_id,
    celltype = anno_qc$celltype,
    condition = anno_qc$condition,
    donor = anno_qc$donor,
    NDUFB8 = as.numeric(counts_qc[, "NDUFB8"])
  )

  ct_summary <- ndufb8_df %>%
    group_by(celltype) %>%
    summarise(
      mean_expr = mean(NDUFB8),
      pct_expressing = mean(NDUFB8 > 0) * 100,
      n_cells = n(),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_expr))

  log_info("NDUFB8 各细胞类型表达:")
  print(ct_summary, row.names = FALSE)
  write.csv(ct_summary, file.path(tab_dir, "NDUFB8_by_celltype.csv"),
            row.names = FALSE)

  overall_pct <- mean(ndufb8_df$NDUFB8 > 0) * 100
  log_info("NDUFB8 全局表达比例: {round(overall_pct, 3)}%")

  sparsity_note <- data.frame(
    n_total = nrow(ndufb8_df),
    n_NDUFB8_pos = sum(ndufb8_df$NDUFB8 > 0),
    pct_pos = round(overall_pct, 4),
    note = ifelse(overall_pct < 1,
                  "Extremely sparse; single-gene high/low split not feasible",
                  "Adequate for high/low split")
  )
  write.csv(sparsity_note, file.path(tab_dir, "NDUFB8_sparsity_note.csv"),
            row.names = FALSE)

  if (overall_pct < 1) {
    log_warn("NDUFB8 单细胞层面过于稀疏，单基因分层分析受限（10X 3' 捕获偏好）")
  }

  p1 <- ggplot(ct_summary, aes(x = reorder(celltype, -mean_expr),
                                 y = mean_expr, fill = pct_expressing)) +
    geom_col() +
    scale_fill_viridis_c(name = "% expressing") +
    labs(title = "NDUFB8 expression across retinal cell types",
         subtitle = sprintf("Overall detection rate: %.3f%% of cells", overall_pct),
         x = "Cell type", y = "Mean expression") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 10))
  ggsave(file.path(fig_dir, "NDUFB8_by_celltype.pdf"), p1,
         width = 10, height = 6)
}

# ============================================================
# ============================================================
log_info("========== V2: NDUFB8 表达分层下的 Complex I / OXPHOS ==========")

if ("NDUFB8" %in% colnames(counts_qc) &&
    length(ci_in) >= 5 && length(oxphos_in) >= 10) {

  ndufb8_expr <- as.numeric(counts_qc[, "NDUFB8"])
  ci_score <- Matrix::rowMeans(counts_qc[, ci_in, drop = FALSE])
  oxphos_score <- Matrix::rowMeans(counts_qc[, oxphos_in, drop = FALSE])

  has_expr <- ndufb8_expr > 0
  log_info("NDUFB8 表达的细胞数: {sum(has_expr)} / {length(has_expr)}")

  if (sum(has_expr) >= 50) {
    ndufb8_has <- ndufb8_expr[has_expr]
    ci_has <- ci_score[has_expr]
    ox_has <- oxphos_score[has_expr]

    ndufb8_high <- ndufb8_has > median(ndufb8_has)

    wt_ci <- wilcox.test(ci_has[ndufb8_high], ci_has[!ndufb8_high])
    wt_ox <- wilcox.test(ox_has[ndufb8_high], ox_has[!ndufb8_high])

    log_info("NDUFB8 高 vs 低（表达细胞内）:")
    log_info("  Complex I: 高={round(mean(ci_has[ndufb8_high]),4)}, 低={round(mean(ci_has[!ndufb8_high]),4)}, p={format(wt_ci$p.value, scientific=TRUE)}")
    log_info("  OXPHOS:    高={round(mean(ox_has[ndufb8_high]),4)}, 低={round(mean(ox_has[!ndufb8_high]),4)}, p={format(wt_ox$p.value, scientific=TRUE)}")

    layer_results <- data.frame()
    for (ct in unique(anno_qc$celltype)) {
      idx <- anno_qc$celltype == ct & has_expr
      if (sum(idx) < 30) next

      ci_ct <- ci_score[idx]
      ndufb8_ct <- ndufb8_expr[idx]
      high_ct <- ndufb8_ct > median(ndufb8_ct)
      if (sum(high_ct) < 10 || sum(!high_ct) < 10) next

      wt <- wilcox.test(ci_ct[high_ct], ci_ct[!high_ct])
      layer_results <- rbind(layer_results, data.frame(
        celltype = ct, n_cells = sum(idx),
        mean_ci_high = mean(ci_ct[high_ct]),
        mean_ci_low  = mean(ci_ct[!high_ct]),
        log2_ratio   = log2(mean(ci_ct[high_ct]) /
                            max(mean(ci_ct[!high_ct]), 1e-6)),
        wilcox_p     = wt$p.value,
        stringsAsFactors = FALSE))
    }

    if (nrow(layer_results) > 0) {
      layer_results <- layer_results[order(-layer_results$log2_ratio), ]
      write.csv(layer_results,
                file.path(tab_dir, "NDUFB8_highlow_by_celltype.csv"),
                row.names = FALSE)
      log_info("按细胞类型分层:")
      print(layer_results, row.names = FALSE)
    } else {
      log_info("按细胞类型分层: 无足够细胞类型满足最小样本数")
    }

    plot_df <- data.frame(
      Complex_I = ci_has,
      group = ifelse(ndufb8_high, "NDUFB8_high", "NDUFB8_low"))

    p2 <- ggplot(plot_df, aes(x = group, y = Complex_I, fill = group)) +
      geom_boxplot(alpha = 0.7, outlier.size = 0.3) +
      labs(title = "Complex I expression in NDUFB8-high vs low cells",
           subtitle = sprintf("Among NDUFB8-expressing cells (n=%d)", sum(has_expr)),
           x = "Group", y = "Complex I mean expression") +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = 10),
            legend.position = "none")
    ggsave(file.path(fig_dir, "NDUFB8_highlow_ComplexI.pdf"), p2,
           width = 7, height = 5)
  } else {
    log_warn("NDUFB8 表达细胞数不足（{sum(has_expr)}），跳过 V2 高/低对比")
  }
}

# ============================================================
# ============================================================
log_info("========== V3: UCell 评分 dAMD vs Healthy ==========")

gene_sets_ucell <- list(
  Complex_I   = ci_in,
  OXPHOS      = oxphos_in,
  Signature16 = genes16_in
)

log_info("计算 UCell 评分（基因集：Complex_I / OXPHOS / Signature16）...")

genes_needed <- unique(c(ci_in, oxphos_in, genes16_in))
log_info("UCell 输入基因数: {length(genes_needed)}")

counts_sub <- as.matrix(counts_qc[, genes_needed, drop = FALSE])

ucell_scores <- as.data.frame(ScoreSignatures_UCell(
  t(counts_sub),
  features = gene_sets_ucell,
  ncores   = 1
))

stopifnot(nrow(ucell_scores) == nrow(anno_qc))
log_info("UCell 返回类型: {class(ucell_scores)[1]}")
log_info("UCell 返回维度: {nrow(ucell_scores)} cells × {ncol(ucell_scores)} signatures")
log_info("UCell 返回列名: {paste(colnames(ucell_scores), collapse=', ')}")

anno_qc$UCell_CI  <- ucell_scores[[pick_col(ucell_scores, "Complex_I")]]
anno_qc$UCell_OX  <- ucell_scores[[pick_col(ucell_scores, "OXPHOS")]]
anno_qc$UCell_G16 <- ucell_scores[[pick_col(ucell_scores, "Signature16")]]

stopifnot(all(c("UCell_CI", "UCell_OX", "UCell_G16") %in% colnames(anno_qc)))
log_info("anno_qc UCell 列已加入: UCell_CI, UCell_OX, UCell_G16")

log_info("UCell 评分摘要:")
print(summary(anno_qc[, c("UCell_CI", "UCell_OX", "UCell_G16")]))

donor_agg <- anno_qc %>%
  group_by(donor, condition) %>%
  summarise(
    mean_UCell_CI  = mean(UCell_CI),
    mean_UCell_OX  = mean(UCell_OX),
    mean_UCell_G16 = mean(UCell_G16),
    n_cells = n(),
    .groups = "drop")

write.csv(donor_agg, file.path(tab_dir, "donor_UCell_scores.csv"),
          row.names = FALSE)
log_info("donor-level UCell 评分（细胞级平均）:")
print(as.data.frame(donor_agg), row.names = FALSE)

donor_dh <- donor_agg[donor_agg$condition %in% c("Healthy", "dAMD"), ]

for (metric in c("mean_UCell_CI", "mean_UCell_OX", "mean_UCell_G16")) {
  h <- donor_dh[[metric]][donor_dh$condition == "Healthy"]
  d <- donor_dh[[metric]][donor_dh$condition == "dAMD"]
  if (length(h) >= 3 && length(d) >= 3) {
    tt <- t.test(d, h)
    log2fc <- log2(mean(d) / mean(h))
    log_info("{metric}: Healthy={round(mean(h),3)}, dAMD={round(mean(d),3)}, log2FC={round(log2fc,3)}, p={format(tt$p.value, scientific=TRUE)}")
  }
}

log_info("--- V3b: pseudobulk 敏感性分析 ---")

donor_dh_agg <- donor_agg %>%
  filter(condition %in% c("Healthy", "dAMD"))

donor_dh_agg$weighted_CI  <- donor_dh_agg$mean_UCell_CI  * donor_dh_agg$n_cells
donor_dh_agg$weighted_OX  <- donor_dh_agg$mean_UCell_OX  * donor_dh_agg$n_cells
donor_dh_agg$weighted_G16 <- donor_dh_agg$mean_UCell_G16 * donor_dh_agg$n_cells

write.csv(donor_dh_agg,
          file.path(tab_dir, "pseudobulk_donor_scores.csv"),
          row.names = FALSE)

log_info("pseudobulk donor 级摘要:")
print(as.data.frame(donor_dh_agg), row.names = FALSE)

pb_results <- data.frame()
for (metric_base in c("mean_UCell_CI", "mean_UCell_OX", "mean_UCell_G16")) {
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

  wcol <- sub("mean_UCell_", "weighted_", metric_base)
  if (wcol %in% colnames(donor_dh_agg)) {
    h_w <- donor_dh_agg[[wcol]][donor_dh_agg$condition == "Healthy"] /
           donor_dh_agg$n_cells[donor_dh_agg$condition == "Healthy"]
    d_w <- donor_dh_agg[[wcol]][donor_dh_agg$condition == "dAMD"] /
           donor_dh_agg$n_cells[donor_dh_agg$condition == "dAMD"]
    if (length(h_w) >= 3 && length(d_w) >= 3) {
      tt_w <- t.test(d_w, h_w)
      pb_results <- rbind(pb_results, data.frame(
        metric = metric_base,
        method = "donor_mean_cellweighted",
        mean_healthy = mean(h_w),
        mean_damd = mean(d_w),
        log2FC = log2(mean(d_w) / mean(h_w)),
        p = tt_w$p.value,
        stringsAsFactors = FALSE))
    }
  }

  pooled_h <- mean(anno_qc[[sub("mean_UCell_", "UCell_", metric_base)]][
    anno_qc$condition == "Healthy"])
  pooled_d <- mean(anno_qc[[sub("mean_UCell_", "UCell_", metric_base)]][
    anno_qc$condition == "dAMD"])
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

log_info("--- V3c: RGC n_cells 检查 + 单次亚采样 ---")

rgc_cells <- anno_qc %>%
  filter(condition %in% c("Healthy", "dAMD"), celltype == "RGC") %>%
  group_by(donor, condition) %>%
  summarise(n_RGC = n(), .groups = "drop")

write.csv(rgc_cells,
          file.path(tab_dir, "RGC_cells_per_donor.csv"),
          row.names = FALSE)

log_info("RGC 细胞数按 condition 分布:")
rgc_by_cond <- rgc_cells %>%
  group_by(condition) %>%
  summarise(
    n_donors = n(),
    mean_n = mean(n_RGC),
    sd_n = sd(n_RGC),
    median_n = median(n_RGC),
    min_n = min(n_RGC),
    max_n = max(n_RGC),
    .groups = "drop")
print(as.data.frame(rgc_by_cond), row.names = FALSE)

rgc_tt <- tryCatch(
  t.test(n_RGC ~ condition, data = rgc_cells),
  error = function(e) NULL)
if (!is.null(rgc_tt)) {
  log_info("RGC n_cells: dAMD vs Healthy 的 t-test p = {format(rgc_tt$p.value, scientific=TRUE)}")
}

min_n <- min(rgc_cells$n_RGC)
log_info("RGC 最小细胞数（用于亚采样）: {min_n}")

rgc_all <- anno_qc %>%
  filter(condition %in% c("Healthy", "dAMD"), celltype == "RGC")

if (min_n >= 10) {
  set.seed(2026)
  rgc_subsampled <- rgc_all %>%
    group_by(donor, condition) %>%
    sample_n(size = min(min_n, n()), replace = FALSE) %>%
    ungroup()

  rgc_sub_summary <- rgc_subsampled %>%
    group_by(donor, condition) %>%
    summarise(
      mean_UCell_CI  = mean(UCell_CI),
      mean_UCell_OX  = mean(UCell_OX),
      mean_UCell_G16 = mean(UCell_G16),
      n_used = n(),
      .groups = "drop")

  write.csv(rgc_sub_summary,
            file.path(tab_dir, "RGC_subsampled_scores.csv"),
            row.names = FALSE)

  log_info("RGC 单次亚采样（每 donor {min_n} 个细胞）后的 dAMD vs Healthy:")
  for (metric in c("mean_UCell_CI", "mean_UCell_OX", "mean_UCell_G16")) {
    h <- rgc_sub_summary[[metric]][rgc_sub_summary$condition == "Healthy"]
    d <- rgc_sub_summary[[metric]][rgc_sub_summary$condition == "dAMD"]
    if (length(h) >= 3 && length(d) >= 3) {
      tt <- t.test(d, h)
      log2fc <- log2(mean(d) / mean(h))
      log_info("  {metric}: Healthy={round(mean(h),3)}, dAMD={round(mean(d),3)}, log2FC={round(log2fc,3)}, p={format(tt$p.value, scientific=TRUE)}")
    }
  }

  rgc_plot <- bind_rows(
    data.frame(donor = rgc_sub_summary$donor,
               condition = rgc_sub_summary$condition,
               score = rgc_sub_summary$mean_UCell_OX,
               method = sprintf("subsampled (%d cells/donor)", min_n)),
    data.frame(donor = donor_agg$donor,
               condition = donor_agg$condition,
               score = donor_agg$mean_UCell_OX,
               method = "all cells")
  ) %>%
    filter(condition %in% c("Healthy", "dAMD"))

  p_rgc <- ggplot(rgc_plot, aes(x = condition, y = score, fill = condition)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.15, alpha = 0.6, size = 1.5) +
    facet_wrap(~ method) +
    labs(title = "RGC OXPHOS score: before vs after subsampling",
         x = "Condition", y = "UCell OXPHOS") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "none")
  ggsave(file.path(fig_dir, "RGC_subsampling_check.pdf"), p_rgc,
         width = 8, height = 4)
} else {
  log_warn("RGC 最小细胞数 < 10，跳过亚采样")
}

log_info("--- V3d: RGC OXPHOS 100 次亚采样稳定性检验 ---")

n_iter <- 100
rgc_sub_iter <- data.frame()

if (min_n >= 10 && nrow(rgc_all) > 0) {
  for (i in seq_len(n_iter)) {
    set.seed(2026 + i)

    rgc_sub_i <- rgc_all %>%
      group_by(donor, condition) %>%
      sample_n(size = min(min_n, n()), replace = FALSE) %>%
      ungroup() %>%
      group_by(donor, condition) %>%
      summarise(
        mean_OX  = mean(UCell_OX),
        mean_CI  = mean(UCell_CI),
        mean_G16 = mean(UCell_G16),
        .groups = "drop")

    h_ox <- rgc_sub_i$mean_OX[rgc_sub_i$condition == "Healthy"]
    d_ox <- rgc_sub_i$mean_OX[rgc_sub_i$condition == "dAMD"]
    h_ci <- rgc_sub_i$mean_CI[rgc_sub_i$condition == "Healthy"]
    d_ci <- rgc_sub_i$mean_CI[rgc_sub_i$condition == "dAMD"]
    h_g16 <- rgc_sub_i$mean_G16[rgc_sub_i$condition == "Healthy"]
    d_g16 <- rgc_sub_i$mean_G16[rgc_sub_i$condition == "dAMD"]

    if (length(h_ox) >= 3 && length(d_ox) >= 3) {
      tt_ox  <- t.test(d_ox,  h_ox)
      tt_ci  <- t.test(d_ci,  h_ci)
      tt_g16 <- t.test(d_g16, h_g16)

      rgc_sub_iter <- rbind(rgc_sub_iter, data.frame(
        iter = i,
        log2FC_OX  = log2(mean(d_ox)  / mean(h_ox)),
        log2FC_CI  = log2(mean(d_ci)  / mean(h_ci)),
        log2FC_G16 = log2(mean(d_g16) / mean(h_g16)),
        p_OX  = tt_ox$p.value,
        p_CI  = tt_ci$p.value,
        p_G16 = tt_g16$p.value,
        stringsAsFactors = FALSE))
    }
  }

  write.csv(rgc_sub_iter,
            file.path(tab_dir, "RGC_subsampling_100iter.csv"),
            row.names = FALSE)

  iter_summary <- data.frame(
    metric = c("OXPHOS", "Complex_I", "Signature16"),
    log2FC_median = c(median(rgc_sub_iter$log2FC_OX,  na.rm = TRUE),
                      median(rgc_sub_iter$log2FC_CI,  na.rm = TRUE),
                      median(rgc_sub_iter$log2FC_G16, na.rm = TRUE)),
    log2FC_IQR_lo = c(quantile(rgc_sub_iter$log2FC_OX,  0.25, na.rm = TRUE),
                      quantile(rgc_sub_iter$log2FC_CI,  0.25, na.rm = TRUE),
                      quantile(rgc_sub_iter$log2FC_G16, 0.25, na.rm = TRUE)),
    log2FC_IQR_hi = c(quantile(rgc_sub_iter$log2FC_OX,  0.75, na.rm = TRUE),
                      quantile(rgc_sub_iter$log2FC_CI,  0.75, na.rm = TRUE),
                      quantile(rgc_sub_iter$log2FC_G16, 0.75, na.rm = TRUE)),
    pct_p_lt_005  = c(mean(rgc_sub_iter$p_OX  < 0.05, na.rm = TRUE) * 100,
                      mean(rgc_sub_iter$p_CI  < 0.05, na.rm = TRUE) * 100,
                      mean(rgc_sub_iter$p_G16 < 0.05, na.rm = TRUE) * 100),
    pct_log2FC_pos = c(mean(rgc_sub_iter$log2FC_OX  > 0, na.rm = TRUE) * 100,
                       mean(rgc_sub_iter$log2FC_CI  > 0, na.rm = TRUE) * 100,
                       mean(rgc_sub_iter$log2FC_G16 > 0, na.rm = TRUE) * 100),
    stringsAsFactors = FALSE)

  write.csv(iter_summary,
            file.path(tab_dir, "RGC_subsampling_100iter_summary.csv"),
            row.names = FALSE)

  log_info("RGC {n_iter} 次亚采样稳定性汇总:")
  print(iter_summary, row.names = FALSE, digits = 4)

  iter_long <- rgc_sub_iter %>%
    select(iter, log2FC_OX, log2FC_CI, log2FC_G16) %>%
    pivot_longer(cols = starts_with("log2FC_"),
                 names_to = "metric", values_to = "log2FC") %>%
    mutate(metric = sub("^log2FC_", "", metric))

  p_iter <- ggplot(iter_long, aes(x = metric, y = log2FC, fill = metric)) +
    geom_violin(alpha = 0.5, scale = "width") +
    geom_boxplot(width = 0.15, alpha = 0.8, outlier.size = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    labs(title = sprintf("RGC UCell scores: %d subsampling iterations", n_iter),
         subtitle = sprintf("Each iteration: %d cells/donor (min n_RGC)", min_n),
         x = "", y = "log2FC (dAMD / Healthy)") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 10),
          legend.position = "none")
  ggsave(file.path(fig_dir, "RGC_100iter_log2FC.pdf"), p_iter,
         width = 7, height = 5)
} else {
  log_warn("RGC 最小细胞数 < 10 或 RGC 数据为空，跳过 100 次亚采样")
}

donor_ct <- anno_qc %>%
  filter(condition %in% c("Healthy", "dAMD")) %>%
  group_by(donor, condition, celltype) %>%
  summarise(
    mean_UCell_CI  = mean(UCell_CI),
    mean_UCell_OX  = mean(UCell_OX),
    mean_UCell_G16 = mean(UCell_G16),
    n_cells = n(),
    .groups = "drop") %>%
  filter(n_cells >= 20)

ct_results <- data.frame()
for (ct in unique(donor_ct$celltype)) {
  for (metric in c("mean_UCell_CI", "mean_UCell_OX", "mean_UCell_G16")) {
    sub <- donor_ct[donor_ct$celltype == ct, ]
    h <- sub[[metric]][sub$condition == "Healthy"]
    d <- sub[[metric]][sub$condition == "dAMD"]
    if (length(h) >= 3 && length(d) >= 3) {
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
            file.path(tab_dir, "UCell_disease_by_celltype.csv"),
            row.names = FALSE)
  log_info("UCell 按细胞类型分层（top 10 效应量）:")
  print(head(ct_results, 10), row.names = FALSE)
}

# ---- 6f. visualization ----
plot_df <- anno_qc %>%
  filter(condition %in% c("Healthy", "dAMD")) %>%
  select(donor, condition, celltype, UCell_CI, UCell_OX, UCell_G16) %>%
  pivot_longer(cols = c("UCell_CI", "UCell_OX", "UCell_G16"),
               names_to = "metric", values_to = "score")

p3 <- ggplot(plot_df, aes(x = condition, y = score, fill = condition)) +
  geom_violin(alpha = 0.7, scale = "width") +
  geom_boxplot(width = 0.2, alpha = 0.5, outlier.size = 0.2) +
  facet_grid(metric ~ celltype, scales = "free_y") +
  labs(title = "UCell scores by condition and cell type",
       x = "Condition", y = "UCell score") +
  theme_minimal(base_size = 9) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

ggsave(file.path(fig_dir, "UCell_by_condition_celltype.pdf"), p3,
       width = 16, height = 10)

# ============================================================
# ============================================================
log_info("========== V4: NDUFB8 表达分层下 16 基因签名 ==========")

if ("NDUFB8" %in% colnames(counts_qc) && length(genes16_in) >= 5) {
  ndufb8_expr <- as.numeric(counts_qc[, "NDUFB8"])
  sig16_score <- Matrix::rowMeans(counts_qc[, genes16_in, drop = FALSE])

  has_expr <- ndufb8_expr > 0
  if (sum(has_expr) >= 50) {
    ndufb8_has <- ndufb8_expr[has_expr]
    sig16_has <- sig16_score[has_expr]
    high <- ndufb8_has > median(ndufb8_has)

    wt <- wilcox.test(sig16_has[high], sig16_has[!high])
    log_info("16 基因签名: 高 vs 低 NDUFB8 细胞")
    log_info("  高均值: {round(mean(sig16_has[high]),4)}")
    log_info("  低均值: {round(mean(sig16_has[!high]),4)}")
    log_info("  p = {format(wt$p.value, scientific=TRUE)}")
    log_info("  log2_ratio = {round(log2(mean(sig16_has[high])/mean(sig16_has[!high])),3)}")

    other16 <- setdiff(genes16_in, "NDUFB8")
    cor_results <- data.frame()
    for (g in other16) {
      g_vec <- as.numeric(counts_qc[has_expr, g])
      if (sd(g_vec) == 0 || sd(ndufb8_has) == 0) next
      r <- cor(ndufb8_has, g_vec, method = "spearman")
      cor_results <- rbind(cor_results, data.frame(
        gene = g, spearman_r = r, stringsAsFactors = FALSE))
    }
    if (nrow(cor_results) > 0) {
      cor_results <- cor_results[order(-abs(cor_results$spearman_r)), ]
      write.csv(cor_results,
                file.path(tab_dir, "NDUFB8_16gene_correlation.csv"),
                row.names = FALSE)
      log_info("NDUFB8 与 16 基因签名其他成员的 Spearman 相关 (top 5):")
      print(head(cor_results, 5), row.names = FALSE)
    }
  } else {
    log_warn("NDUFB8 表达细胞数不足（{sum(has_expr)}），无法做高/低对比")
  }
}

# ============================================================
# ============================================================
log_info("========== V4b: RGC 各 Complex 单独 UCell 评分 ==========")

if (length(complex_genes) >= 2) {
  complex_genes_all <- unique(unlist(complex_genes))
  log_info("Complex UCell 输入基因数: {length(complex_genes_all)}")

  counts_complex <- as.matrix(counts_qc[, complex_genes_all, drop = FALSE])

  complex_ucell <- as.data.frame(ScoreSignatures_UCell(
    t(counts_complex),
    features = complex_genes,
    ncores   = 1
  ))

  stopifnot(nrow(complex_ucell) == nrow(anno_qc))
  log_info("Complex UCell 返回类型: {class(complex_ucell)[1]}")
  log_info("Complex UCell 返回维度: {nrow(complex_ucell)} cells × {ncol(complex_ucell)} signatures")
  log_info("Complex UCell 返回列名: {paste(colnames(complex_ucell), collapse=', ')}")

  for (cp in names(complex_genes)) {
    col_name <- pick_col(complex_ucell, cp)
    anno_qc[[paste0("UCell_", cp)]] <- complex_ucell[[col_name]]
  }

  expected_cols_anno <- paste0("UCell_", names(complex_genes))
  missing_anno <- setdiff(expected_cols_anno, colnames(anno_qc))
  if (length(missing_anno) > 0) {
    stop(sprintf("anno_qc 缺少列: %s", paste(missing_anno, collapse = ", ")))
  }
  log_info("anno_qc Complex 列已加入: {paste(expected_cols_anno, collapse=', ')}")

  rgc_all <- anno_qc %>%
    filter(condition %in% c("Healthy", "dAMD"), celltype == "RGC")

  expected_cols <- expected_cols_anno
  missing_cols <- setdiff(expected_cols, colnames(rgc_all))
  if (length(missing_cols) > 0) {
    stop(sprintf("rgc_all 缺少列: %s", paste(missing_cols, collapse = ", ")))
  }
  log_info("rgc_all Complex 列: {paste(expected_cols, collapse=', ')}")

  if (min_n >= 10 && nrow(rgc_all) > 0) {
    set.seed(2026)
    rgc_sub_complex <- rgc_all %>%
      group_by(donor, condition) %>%
      sample_n(size = min(min_n, n()), replace = FALSE) %>%
      ungroup() %>%
      group_by(donor, condition) %>%
      summarise(
        across(all_of(expected_cols), mean),
        n_used = n(),
        .groups = "drop")

    log_info("rgc_sub_complex 列: {paste(colnames(rgc_sub_complex), collapse=', ')}")

    write.csv(rgc_sub_complex,
              file.path(tab_dir, "RGC_subsampled_complex_scores.csv"),
              row.names = FALSE)

    complex_res <- data.frame()
    for (cp in names(complex_genes)) {
      cn <- paste0("UCell_", cp)
      if (!cn %in% colnames(rgc_sub_complex)) next
      h <- rgc_sub_complex[[cn]][rgc_sub_complex$condition == "Healthy"]
      d <- rgc_sub_complex[[cn]][rgc_sub_complex$condition == "dAMD"]
      if (length(h) >= 3 && length(d) >= 3) {
        tt <- t.test(d, h)
        complex_res <- rbind(complex_res, data.frame(
          complex = cp,
          n_genes = length(complex_genes[[cp]]),
          mean_healthy = mean(h),
          mean_damd = mean(d),
          log2FC = log2(mean(d) / mean(h)),
          p = tt$p.value,
          stringsAsFactors = FALSE))
      }
    }

    if (nrow(complex_res) > 0) {
      complex_res <- complex_res[order(complex_res$p), ]
      write.csv(complex_res,
                file.path(tab_dir, "RGC_complex_breakdown.csv"),
                row.names = FALSE)

      log_info("RGC 各 Complex 拆解（单次亚采样，每 donor {min_n} 细胞）:")
      print(complex_res, row.names = FALSE, digits = 4)

      log_info("--- V4b: RGC 各 Complex 100 次亚采样稳定性 ---")

      complex_iter <- data.frame()
      for (i in seq_len(n_iter)) {
        set.seed(2026 + i)
        sub_i <- rgc_all %>%
          group_by(donor, condition) %>%
          sample_n(size = min(min_n, n()), replace = FALSE) %>%
          ungroup() %>%
          group_by(donor, condition) %>%
          summarise(across(all_of(expected_cols), mean),
                    .groups = "drop")

        row_i <- data.frame(iter = i)
        for (cp in names(complex_genes)) {
          cn <- paste0("UCell_", cp)
          if (!cn %in% colnames(sub_i)) next
          h_i <- sub_i[[cn]][sub_i$condition == "Healthy"]
          d_i <- sub_i[[cn]][sub_i$condition == "dAMD"]
          if (length(h_i) >= 3 && length(d_i) >= 3) {
            row_i[[paste0("log2FC_", cp)]] <- log2(mean(d_i) / mean(h_i))
            row_i[[paste0("p_", cp)]] <- t.test(d_i, h_i)$p.value
          }
        }
        complex_iter <- bind_rows(complex_iter, row_i)
      }

      write.csv(complex_iter,
                file.path(tab_dir, "RGC_complex_100iter.csv"),
                row.names = FALSE)

      complex_iter_summary <- data.frame()
      for (cp in names(complex_genes)) {
        lfc_col <- paste0("log2FC_", cp)
        p_col   <- paste0("p_", cp)
        if (lfc_col %in% colnames(complex_iter)) {
          complex_iter_summary <- rbind(complex_iter_summary, data.frame(
            complex = cp,
            n_genes = length(complex_genes[[cp]]),
            log2FC_median = median(complex_iter[[lfc_col]], na.rm = TRUE),
            log2FC_IQR_lo = quantile(complex_iter[[lfc_col]], 0.25, na.rm = TRUE),
            log2FC_IQR_hi = quantile(complex_iter[[lfc_col]], 0.75, na.rm = TRUE),
            pct_p_lt_005  = mean(complex_iter[[p_col]] < 0.05, na.rm = TRUE) * 100,
            pct_log2FC_pos = mean(complex_iter[[lfc_col]] > 0, na.rm = TRUE) * 100,
            stringsAsFactors = FALSE))
        }
      }

      write.csv(complex_iter_summary,
                file.path(tab_dir, "RGC_complex_100iter_summary.csv"),
                row.names = FALSE)

      log_info("RGC 各 Complex 100 次亚采样稳定性汇总:")
      print(complex_iter_summary, row.names = FALSE, digits = 4)

      # ---- 8c. visualization ----
      complex_iter_long <- complex_iter %>%
        select(iter, starts_with("log2FC_")) %>%
        pivot_longer(cols = starts_with("log2FC_"),
                     names_to = "complex", values_to = "log2FC") %>%
        mutate(complex = sub("^log2FC_", "", complex))

      med_order <- complex_iter_long %>%
        group_by(complex) %>%
        summarise(m = median(log2FC, na.rm = TRUE)) %>%
        arrange(m) %>%
        pull(complex)
      complex_iter_long$complex <- factor(complex_iter_long$complex,
                                           levels = med_order)

      p_complex <- ggplot(complex_iter_long,
                           aes(x = complex, y = log2FC, fill = complex)) +
        geom_violin(alpha = 0.5, scale = "width") +
        geom_boxplot(width = 0.15, alpha = 0.8, outlier.size = 0.3) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
        labs(title = sprintf("RGC Complex scores: %d subsampling iterations", n_iter),
             subtitle = sprintf("Each iteration: %d cells/donor", min_n),
             x = "Complex", y = "log2FC (dAMD / Healthy)") +
        theme_minimal(base_size = 12) +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              plot.subtitle = element_text(hjust = 0.5, size = 10),
              legend.position = "none")
      ggsave(file.path(fig_dir, "RGC_complex_100iter.pdf"), p_complex,
             width = 8, height = 5)

      bulk_path <- file.path(
        "E:/retina_tools/aging_rgc_axon/AMD",
        "results", "07_5c_OXPHOS_analysis", "tables",
        "OXPHOS13_complex_score.csv")
      if (file.exists(bulk_path)) {
        bulk_complex <- read.csv(bulk_path)
        log_info("Bulk 07_5c Complex 对照:")
        print(bulk_complex[, c("complex","n_subunits","beta","p","FDR")],
              row.names = FALSE, digits = 4)

        merged <- merge(
          complex_iter_summary[, c("complex","log2FC_median",
                                     "pct_p_lt_005","pct_log2FC_pos")],
          bulk_complex[, c("complex","beta","p","FDR")],
          by = "complex", all = TRUE)
        write.csv(merged,
                  file.path(tab_dir, "bulk_vs_RGC_complex_comparison.csv"),
                  row.names = FALSE)
        log_info("Bulk vs RGC Complex 对照:")
        print(merged, row.names = FALSE, digits = 4)
      } else {
        log_warn("07_5c 结果不存在，跳过 bulk-single-cell 对照")
      }
    }
  } else {
    log_warn("RGC 最小细胞数 < 10 或 RGC 数据为空，跳过 V4b")
  }
} else {
  log_warn("Complex 基因集不足（< 2），跳过 V4b")
}

# ============================================================
# 9. summary
# ============================================================
log_info("==================== 07_3 v10 完成 ====================")
log_info("输出目录: {out_dir}")

cat("\n", strrep("=", 70), "\n", sep = "")
cat("=== 07_3 NDUFB8 co-expression 分析汇总 (v10) ===\n")
cat(strrep("=", 70), "\n", sep = "")

cat("\n【V1】NDUFB8 各细胞类型表达\n")
cat("  → 见 tables/NDUFB8_by_celltype.csv, NDUFB8_sparsity_note.csv\n")

cat("\n【V2】NDUFB8 表达分层下 Complex I / OXPHOS 差异\n")
cat("  → 见 tables/NDUFB8_highlow_by_celltype.csv\n")

cat("\n【V3】UCell 评分 dAMD vs Healthy\n")
cat("  → 见 tables/donor_UCell_scores.csv, UCell_disease_by_celltype.csv\n")

cat("\n【V3b】pseudobulk 敏感性分析\n")
cat("  → 见 tables/pseudobulk_donor_scores.csv, pseudobulk_sensitivity.csv\n")

cat("\n【V3c】RGC n_cells 检查 + 单次亚采样\n")
cat("  → 见 tables/RGC_cells_per_donor.csv, RGC_subsampled_scores.csv\n")

cat("\n【V3d】RGC OXPHOS 100 次亚采样稳定性\n")
cat("  → 见 tables/RGC_subsampling_100iter_summary.csv\n")
cat("  → 见 figures/RGC_100iter_log2FC.pdf\n")

cat("\n【V4】NDUFB8 表达分层下 16 基因签名差异\n")
cat("  → 见 tables/NDUFB8_16gene_correlation.csv\n")

cat("\n【V4b】RGC 各 Complex（I-V）单独 UCell 评分\n")
cat("  → 见 tables/RGC_complex_breakdown.csv\n")
cat("  → 见 tables/RGC_complex_100iter_summary.csv\n")
cat("  → 见 tables/bulk_vs_RGC_complex_comparison.csv\n")
cat("  → 见 figures/RGC_complex_100iter.pdf\n")

log_info("=== 07_3 v10 done ===")

