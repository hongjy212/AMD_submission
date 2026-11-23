#!/usr/bin/env Rscript
# ============================================================
# 05_1_consensus_wgcna.r
#
# Output：results/05_1_consensus_wgcna/
#   data/consensusNet.rds
#   data/moduleMEs_list.rds
#   data/moduleColors_list.rds
#   data/moduleSizes_list.rds
#   tables/module_summary.csv
#   tables/consensus_module_genes_annotated.csv
#   tables/module_genes_<color>.csv
#   tables/moduleMEs_<grade>.csv
#   figures/soft_threshold_<grade>.pdf
#   figures/consensus_dendrogram*.pdf
#   figures/module_ME_correlation_<grade>.pdf
# ============================================================

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(WGCNA); library(dplyr); library(ggplot2)
})
allowWGCNAThreads()

DIR_THIS <- setup_output_dir("05_1_consensus_wgcna")
log_info("=== 05_1_consensus_wgcna.r start ===")

# ============================================================
# 0. Parameters
# ============================================================
MIN_MODULE_SIZE  <- 30    # 从 50 降到 30，适配 G4 只 58 样本
DEEP_SPLIT       <- 2
MERGE_CUT_HEIGHT <- 0.25
MAX_BLOCK_SIZE   <- 20000 # 15862 基因一 block 装下
COR_TYPE         <- "bicor"

# ============================================================
# ============================================================
log_info("--- 1. 载入 01_1b 输出 ---")

expr <- read_upstream("01_1b_batch_enhance", "data", "expr_wgcna.rds")
clin <- read_upstream("01_1b_batch_enhance", "data", "clinical_wgcna.rds")

if (nrow(expr) > 15000) {
  log_warn("expr 有 {nrow(expr)} 基因，未做 protein-coding 过滤")
  log_warn("建议先做 protein-coding 过滤，或用全部")
}

log_info("expr: {nrow(expr)} × {ncol(expr)}")
log_info("clinical: {nrow(clin)} 样本")

# align
common <- intersect(colnames(expr), rownames(clin))
expr <- expr[, common, drop = FALSE]
clin <- clin[common, , drop = FALSE]

log_info("对齐后: {ncol(expr)} 样本")

# ============================================================
# ============================================================
log_info("--- 2. 构建 multiExpr ---")

clin$grade_int <- as.integer(clin$grade)
grade_levels <- sort(unique(clin$grade_int))
log_info("条件: {paste('G' + grade_levels, collapse=', ')}")

multiExpr <- list()
sampleTable <- data.frame()

for (g in grade_levels) {
  samples_g <- rownames(clin)[clin$grade_int == g]
  if (length(samples_g) < 5) {
    log_warn("G{g}: {length(samples_g)} 样本 < 5，跳过")
    next
  }
  expr_g <- t(expr[, samples_g, drop = FALSE])
  cond_name <- paste0("G", g)
  multiExpr[[cond_name]] <- list(data = expr_g)
  log_info("{cond_name}: {nrow(expr_g)} 样本 × {ncol(expr_g)} 基因")

  sampleTable <- rbind(sampleTable, data.frame(
    Sample    = samples_g,
    Condition = cond_name,
    Grade     = g,
    stringsAsFactors = FALSE))
}
rownames(sampleTable) <- sampleTable$Sample

if (length(multiExpr) < 2) stop("multiExpr 少于 2 个条件")

# ============================================================
# ============================================================
log_info("--- 3. goodSamplesGenes 质控 ---")

for (name in names(multiExpr)) {
  log_info("检查 {name}")
  expr_t <- multiExpr[[name]]$data
  gsg <- goodSamplesGenes(expr_t, verbose = 0)
  if (!gsg$allOK) {
    log_warn("  {name}: 剔除 {sum(!gsg$goodSamples)} 样本, {sum(!gsg$goodGenes)} 基因")
    expr_t <- expr_t[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
    multiExpr[[name]]$data <- expr_t
  }
}

common_genes <- Reduce(intersect, lapply(multiExpr, function(x) colnames(x$data)))
log_info("共有基因: {length(common_genes)}")
for (name in names(multiExpr)) {
  multiExpr[[name]]$data <- multiExpr[[name]]$data[, common_genes, drop = FALSE]
  log_info("  {name}: {nrow(multiExpr[[name]]$data)} 样本 × {ncol(multiExpr[[name]]$data)} 基因")
}

# ============================================================
# ============================================================
log_info("--- 4. 软阈值评估 ---")

powers <- c(seq(1, 20, by = 2), seq(22, 30, by = 2))
soft_powers <- list()

for (name in names(multiExpr)) {
  log_info("软阈值 {name}")
  datExpr <- multiExpr[[name]]$data

  sft <- pickSoftThreshold(datExpr,
                            powerVector = powers,
                            networkType = "signed",
                            verbose = 0)
  power_est <- sft$powerEstimate
  if (is.na(power_est) || power_est < 1) power_est <- 6
  soft_powers[[name]] <- power_est
  log_info("  {name}: softPower = {power_est}")

  pdf(file.path(DIR_THIS, "figures", paste0("soft_threshold_", name, ".pdf")),
      width = 8, height = 6)
  par(mfrow = c(1, 2))
  plot(sft$fitIndices[, 1],
       -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       xlab = "Soft Threshold (power)",
       ylab = "Scale Free Topology Model Fit (signed R²)",
       type = "n", main = paste(name, "- Scale independence"))
  text(sft$fitIndices[, 1],
       -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
       labels = powers, col = "red")
  abline(h = 0.85, col = "blue")
  plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
       xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
       type = "n", main = paste(name, "- Mean connectivity"))
  text(sft$fitIndices[, 1], sft$fitIndices[, 5],
       labels = powers, col = "red")
  dev.off()
}

softPower <- max(unlist(soft_powers))
if (softPower < 6) softPower <- 6
log_info("consensus softPower: {softPower}")

# ============================================================
# 5. blockwiseConsensusModules
# ============================================================
log_info("--- 5. blockwiseConsensusModules ---")
log_info("参数: minModuleSize={MIN_MODULE_SIZE}, deepSplit={DEEP_SPLIT}, mergeCutHeight={MERGE_CUT_HEIGHT}")

t0 <- Sys.time()
consensusNet <- blockwiseConsensusModules(
  multiExpr       = multiExpr,
  power           = softPower,
  networkType     = "signed",
  TOMType         = "signed",
  corType         = COR_TYPE,
  minModuleSize   = MIN_MODULE_SIZE,
  deepSplit       = DEEP_SPLIT,
  mergeCutHeight  = MERGE_CUT_HEIGHT,
  maxBlockSize    = MAX_BLOCK_SIZE,
  saveDendrograms = TRUE,
  verbose         = 2
)
log_info("WGCNA 耗时: {round(as.numeric(difftime(Sys.time(), t0, units='mins')), 1)} 分钟")

n_modules_total <- length(unique(consensusNet$colors))
n_grey <- sum(consensusNet$colors == "grey")
log_info("模块总数（含 grey）: {n_modules_total}")
log_info("grey 基因数: {n_grey}")
save_rds(consensusNet, DIR_THIS, "consensusNet.rds")

# ============================================================
# ============================================================
log_info("--- 6. 绘制聚类树图 ---")

if (!is.null(consensusNet$dendrograms) && length(consensusNet$dendrograms) > 0) {
  nBlocks <- length(consensusNet$dendrograms)
  color_labels_all <- labels2colors(consensusNet$colors)
  names(color_labels_all) <- colnames(multiExpr[[1]]$data)

  for (b in seq_len(nBlocks)) {
    pdf(file.path(DIR_THIS, "figures",
                   paste0("consensus_dendrogram_block", b, ".pdf")),
        width = 14, height = 8)
    blockGenes <- consensusNet$blockGenes[[b]]
    plotDendroAndColors(consensusNet$dendrograms[[b]],
                         color_labels_all[blockGenes],
                         "Modules",
                         dendroLabels = FALSE, hang = 0.03,
                         addGuide = TRUE, guideHang = 0.05,
                         main = paste("Consensus Module Dendrogram - Block", b))
    dev.off()
  }
  log_info("树图已保存（{nBlocks} 个 block）")
}

# ============================================================
# ============================================================
log_info("--- 7. 提取各条件模块 MEs ---")

color_labels_all <- consensusNet$colors
names(color_labels_all) <- colnames(multiExpr[[1]]$data)

moduleColors_list <- list()
moduleMEs_list    <- list()
moduleSizes_list  <- list()

for (name in names(multiExpr)) {
  MEs <- moduleEigengenes(multiExpr[[name]]$data,
                          colors = color_labels_all)$eigengenes
  keep <- !grepl("^MEgrey", colnames(MEs))
  MEs <- MEs[, keep, drop = FALSE]

  moduleMEs_list[[name]]    <- MEs
  moduleColors_list[[name]] <- color_labels_all
  moduleSizes_list[[name]]  <- table(color_labels_all)

  log_info("{name}: {ncol(MEs)} 模块")
}

save_rds(moduleColors_list, DIR_THIS, "moduleColors_list.rds")
save_rds(moduleMEs_list,    DIR_THIS, "moduleMEs_list.rds")
save_rds(moduleSizes_list,  DIR_THIS, "moduleSizes_list.rds")

# ============================================================
# ============================================================
log_info("--- 8. 模块统计 ---")

summary_df <- data.frame()
for (name in names(moduleSizes_list)) {
  sizes <- moduleSizes_list[[name]]
  grey_total <- sum(sizes[names(sizes) == "grey"])
  n_modules  <- sum(names(sizes) != "grey")
  n_samples  <- nrow(multiExpr[[name]]$data)
  real_sizes <- sizes[names(sizes) != "grey"]
  max_size   <- if (length(real_sizes) > 0) max(real_sizes) else 0

  summary_df <- rbind(summary_df, data.frame(
    Condition       = name,
    N_Samples       = n_samples,
    N_Modules       = n_modules,
    Grey_Genes      = as.integer(grey_total),
    Grey_Pct        = round(grey_total / sum(sizes) * 100, 1),
    Max_Module_Size = max_size,
    stringsAsFactors = FALSE))
}
save_table(summary_df, DIR_THIS, "module_summary.csv")
print(summary_df, row.names = FALSE)

# ============================================================
# ============================================================
log_info("--- 9. ME 相关性热图 ---")

for (name in names(moduleMEs_list)) {
  MEs <- moduleMEs_list[[name]]
  if (is.null(MEs) || ncol(MEs) <= 1) next

  cor_MEs <- cor(MEs, use = "pairwise.complete.obs")
  pdf(file.path(DIR_THIS, "figures",
                 paste0("module_ME_correlation_", name, ".pdf")),
      width = 10, height = 8)
  par(cex = 0.9)
  WGCNA::labeledHeatmap(Matrix = cor_MEs,
                         xLabels = colnames(MEs),
                         yLabels = colnames(MEs),
                         colorLabels = TRUE,
                         colors = blueWhiteRed(50),
                         textMatrix = round(cor_MEs, 2),
                         setStdMargins = FALSE,
                         cex.text = 0.5,
                         main = paste("ME Correlation -", name))
  dev.off()
}

# ============================================================
# ============================================================
log_info("--- 10. 导出模块基因表 ---")

gene_names <- colnames(multiExpr[[1]]$data)
stopifnot(length(gene_names) == length(consensusNet$colors))


anno_path <- file.path(DIR_RESULT, "01_1_prepare_bulk_data", "data",
                        "gene_annotation.rds")
if (file.exists(anno_path)) {
  anno <- readRDS(anno_path)
} else {
  anno <- NULL
  log_warn("未找到 gene_annotation.rds")
}

module_df <- data.frame(
  ENSEMBL = gene_names,
  Module  = consensusNet$colors,
  stringsAsFactors = FALSE)

if (!is.null(anno) && all(c("ensembl", "symbol", "description") %in% colnames(anno))) {
  idx <- match(module_df$ENSEMBL, anno$ensembl)
  module_df$SYMBOL <- ifelse(is.na(idx), module_df$ENSEMBL,
                              ifelse(is.na(anno$symbol[idx]), module_df$ENSEMBL,
                                     anno$symbol[idx]))
  module_df$Description <- ifelse(is.na(idx), "No annotation",
                                   ifelse(is.na(anno$description[idx]), "No annotation",
                                          anno$description[idx]))
} else {
  module_df$SYMBOL <- module_df$ENSEMBL
  module_df$Description <- "No annotation"
}

module_df <- module_df[, c("ENSEMBL", "SYMBOL", "Module", "Description")]
save_table(module_df, DIR_THIS, "consensus_module_genes_annotated.csv")

for (mod in unique(module_df$Module)) {
  sub <- module_df[module_df$Module == mod, ]
  safe_mod <- gsub(" ", "_", mod)
  write.csv(sub,
            file.path(DIR_THIS, "tables", paste0("module_genes_", safe_mod, ".csv")),
            row.names = FALSE, fileEncoding = "UTF-8")
}

for (name in names(moduleMEs_list)) {
  MEs <- moduleMEs_list[[name]]
  if (!is.null(MEs)) {
    write.csv(as.data.frame(MEs),
              file.path(DIR_THIS, "tables", paste0("moduleMEs_", name, ".csv")),
              row.names = TRUE, fileEncoding = "UTF-8")
  }
}

# ============================================================
# 11. save sampleTable
# ============================================================
save_rds(sampleTable, DIR_THIS, "sample_table.rds")

# ============================================================
# 12. summary
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("=== 05_1 汇总 ===\n")
cat(strrep("=", 70), "\n", sep = "")

cat("\n【Q1】样本\n")
cat("  total:", ncol(expr), "\n")
for (name in names(multiExpr)) {
  cat(sprintf("  %s: %d\n", name, nrow(multiExpr[[name]]$data)))
}

cat("\n【Q2】consensus softPower:", softPower, "\n")

cat("\n【Q3】模块\n")
cat("  总数（含 grey）:", n_modules_total, "\n")
cat("  grey 基因数:     ", n_grey, "\n")
cat("  真实模块数:      ", n_modules_total - 1, "\n")

cat("\n【Q4】各条件模块统计\n")
print(summary_df, row.names = FALSE)

cat("\n【Q5】模块大小分布\n")
size_tbl <- table(consensusNet$colors)
print(sort(size_tbl, decreasing = TRUE))

cat("\n--- 输出位置 ---\n")
cat("目录:", DIR_THIS, "\n")

log_info("=== 05_1_consensus_wgcna.r done ===")


