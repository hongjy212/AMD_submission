# ============================================================
# 07_3_scTenifoldKnk_NDUFB3_RPE.r
# Purpose：scTenifoldKnk virtualknockout NDUFB3（RPE，GSE262151）
#
#
# Data source：E:/retina_tools/raw/GSE262151/snRNA_RPE.h5ad
# mappingfile：results/01_1_prepare_bulk_data/data/gene_annotation.rds
#
# ============================================================

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(zellkonverter)
  library(scTenifoldKnk)
  library(Matrix); library(dplyr); library(tidyr)
  library(ggplot2); library(fgsea)
})

# ============================================================
# 0. Parameters
# ============================================================
H5AD_FILE      <- "E:/retina_tools/raw/GSE262151/snRNA_RPE.h5ad"
GENE_ANNO_FILE <- "E:/retina_tools/aging_rgc_axon/AMD/results/01_1_prepare_bulk_data/data/gene_annotation.rds"

GENE_KO      <- "NDUFB3"
CELL_TYPE_KO <- "RPE"
N_CORES      <- 8
SEED         <- 42

DIR_THIS <- setup_output_dir(sprintf("07_3_scTenifoldKnk_%s_%s",
                                      GENE_KO, CELL_TYPE_KO))
log_info(sprintf("=== 07_3 scTenifoldKnk %s (%s) start ===",
                 GENE_KO, CELL_TYPE_KO))

# QC Parameters
QC_MIN_LIBSIZE <- 500
QC_MAX_MTRATIO <- 0.10
QC_MIN_PCT     <- 0.02

# scTenifoldKnk Parameters
NC_NNET     <- 5
NC_NCELLS   <- 300
NC_NCOMP    <- 3
TD_K        <- 3
TD_MAXITER  <- 500
TD_MAXERROR <- 1e-4

MIN_GENE_CELLS <- 15

# ============================================================
# ============================================================
log_info("--- 0. 读取外部基因集 ---")

GENE_SETS <- tryCatch(
  read_upstream("01_2_prepare_gene_sets", "data", "gene_sets.rds"),
  error = function(e) {
    log_warn(sprintf("gene_sets.rds 读取失败: %s", conditionMessage(e)))
    NULL
  }
)

CI_GENES <- character(0)
GENES_16 <- character(0)
GY_GENES <- character(0)
BL_GENES <- character(0)

if (!is.null(GENE_SETS) && "Mito_CI_Subunits" %in% names(GENE_SETS)) {
  CI_GENES <- toupper(GENE_SETS[["Mito_CI_Subunits"]]$genes)
  log_info(sprintf("Mito_CI_Subunits: %d 基因", length(CI_GENES)))
}

sig16_file <- file.path(DIR_RESULT, "02_1_differential_expression_main",
                         "tables", "sig16_genes.csv")
if (file.exists(sig16_file)) {
  GENES_16 <- toupper(read.csv(sig16_file)$gene)
  log_info(sprintf("16 基因签名: %d", length(GENES_16)))
}

gy_file <- file.path(DIR_RESULT, "05_1_consensus_wgcna", "tables",
                      "module_genes_greenyellow_symbolfix.csv")
if (file.exists(gy_file)) {
  g <- toupper(read.csv(gy_file)$SYMBOL)
  GY_GENES <- g[!is.na(g) & g != ""]
  log_info(sprintf("greenyellow: %d", length(GY_GENES)))
}

bl_file <- file.path(DIR_RESULT, "05_1_consensus_wgcna", "tables",
                      "module_genes_black_symbolfix.csv")
if (file.exists(bl_file)) {
  g <- toupper(read.csv(bl_file)$SYMBOL)
  BL_GENES <- g[!is.na(g) & g != ""]
  log_info(sprintf("black: %d", length(BL_GENES)))
}

# ============================================================
# 2. read h5ad
# ============================================================
log_info("--- 1. 读取 h5ad ---")

if (!file.exists(H5AD_FILE)) stop(sprintf("h5ad 不存在: %s", H5AD_FILE))

sce <- readH5AD(H5AD_FILE, reader = "R")
log_info(sprintf("h5ad: %d genes × %d cells", nrow(sce), ncol(sce)))
log_info(sprintf("assayNames: %s", paste(assayNames(sce), collapse = ", ")))

# ============================================================
# 3. Ensembl → Symbol mapping
# ============================================================
log_info("--- 2. Ensembl → Symbol 映射 ---")

if (!file.exists(GENE_ANNO_FILE)) {
  stop(sprintf("gene_annotation.rds 不存在: %s", GENE_ANNO_FILE))
}

gene_anno <- readRDS(GENE_ANNO_FILE)
log_info(sprintf("gene_annotation: %d 行 × %d 列",
                 nrow(gene_anno), ncol(gene_anno)))

ens_col <- grep("ensembl|^ens", colnames(gene_anno),
                ignore.case = TRUE, value = TRUE)[1]
sym_col <- grep("symbol|gene_name|^sym", colnames(gene_anno),
                ignore.case = TRUE, value = TRUE)[1]

if (is.na(ens_col)) stop("gene_annotation 里找不到 ensembl 列")
if (is.na(sym_col)) stop("gene_annotation 里找不到 symbol 列")

log_info(sprintf("Ensembl 列: %s, Symbol 列: %s", ens_col, sym_col))

ens_ids <- sub("\\..*$", "", as.character(gene_anno[[ens_col]]))
sym_ids <- as.character(gene_anno[[sym_col]])
valid <- !is.na(ens_ids) & !is.na(sym_ids) & ens_ids != "" & sym_ids != ""
symbol_map <- setNames(sym_ids[valid], ens_ids[valid])
log_info(sprintf("映射表: %d 条", length(symbol_map)))

rn <- sub("\\..*$", "", rownames(sce))
keep_gene <- rn %in% names(symbol_map)
log_info(sprintf("可映射基因: %d / %d (%.1f%%)",
                 sum(keep_gene), length(rn),
                 100 * sum(keep_gene) / length(rn)))

if (sum(keep_gene) < 1000) stop("可映射基因数太少")

sce_mapped <- sce[keep_gene, ]
new_rownames <- symbol_map[rn[keep_gene]]

expr_sum <- Matrix::rowSums(assay(sce_mapped, "X"))
dup_any <- duplicated(new_rownames) | duplicated(new_rownames, fromLast = TRUE)

if (sum(dup_any) > 0) {
  log_info(sprintf("重复 symbol: %d", sum(dup_any)))
  df_tmp <- data.frame(
    symbol = new_rownames, expr = expr_sum,
    idx = seq_along(new_rownames), stringsAsFactors = FALSE
  )
  df_tmp <- df_tmp[order(df_tmp$symbol, -df_tmp$expr), ]
  keep_idx <- df_tmp$idx[!duplicated(df_tmp$symbol)]
  sce_mapped <- sce_mapped[keep_idx, ]
  new_rownames <- new_rownames[keep_idx]
}

rownames(sce_mapped) <- make.unique(toupper(new_rownames))
log_info(sprintf("映射后基因数: %d", nrow(sce_mapped)))

# extractexpressionmatrix（genes × cells）
mtx_target <- assay(sce_mapped, "X")
if (!inherits(mtx_target, "dgCMatrix")) {
  mtx_target <- as(mtx_target, "CsparseMatrix")
}
rownames(mtx_target) <- rownames(sce_mapped)
colnames(mtx_target) <- make.unique(colnames(mtx_target))

log_info(sprintf("RPE 矩阵: %d genes × %d cells",
                 nrow(mtx_target), ncol(mtx_target)))

# ============================================================
# 4. genesfilter
# ============================================================
log_info("--- 3. 基因过滤 ---")

n_cells_per_gene <- Matrix::rowSums(mtx_target > 0)
keep <- n_cells_per_gene >= MIN_GENE_CELLS

force_keep <- unique(c(GENE_KO, CI_GENES, GENES_16, GY_GENES, BL_GENES))
force_keep <- force_keep[force_keep %in% rownames(mtx_target)]
keep[rownames(mtx_target) %in% force_keep] <- TRUE

mtx_target_f <- mtx_target[keep, , drop = FALSE]
log_info(sprintf("基因过滤: %d → %d", nrow(mtx_target), nrow(mtx_target_f)))

if (!GENE_KO %in% rownames(mtx_target_f)) {
  stop(sprintf("种子基因 %s 不在过滤后矩阵中", GENE_KO))
}

n_expr_target <- sum(mtx_target_f[GENE_KO, ] > 0)
pct_expr      <- 100 * n_expr_target / ncol(mtx_target_f)
log_info(sprintf("%s 中 %s: %d / %d (%.2f%%)",
                 CELL_TYPE_KO, GENE_KO,
                 n_expr_target, ncol(mtx_target_f), pct_expr))

if (n_expr_target < 30) {
  stop(sprintf("%s 中 %s 表达细胞数 (%d) 太少",
               CELL_TYPE_KO, GENE_KO, n_expr_target))
}

# ============================================================
# 5. scTenifoldKnk
# ============================================================
log_info(sprintf("--- 4. scTenifoldKnk (%s, %s) ---",
                 GENE_KO, CELL_TYPE_KO))
log_info(sprintf("参数: nc_nNet=%d, nc_nCells=%d, td_K=%d",
                 NC_NNET, NC_NCELLS, TD_K))

set.seed(SEED)
t0 <- Sys.time()

result <- tryCatch({
  scTenifoldKnk(
    countMatrix     = mtx_target_f,
    gKO             = GENE_KO,
    transcriptomeWide = FALSE,
    qc              = TRUE,
    qc_minLibSize   = QC_MIN_LIBSIZE,
    qc_removeOutlierCells = TRUE,
    qc_minPCT       = QC_MIN_PCT,
    qc_maxMTratio   = QC_MAX_MTRATIO,
    nc_lambda       = 0,
    nc_nNet         = NC_NNET,
    nc_nCells       = min(NC_NCELLS, ncol(mtx_target_f) - 1),
    nc_nComp        = NC_NCOMP,
    nc_scaleScores  = TRUE,
    nc_symmetric    = FALSE,
    nc_q            = 0.9,
    td_K            = TD_K,
    td_maxIter      = TD_MAXITER,
    td_maxError     = TD_MAXERROR,
    td_nDecimal     = 3,
    ma_nDim         = 2,
    dr_empiricalNull = FALSE,
    nCores          = N_CORES
  )
}, error = function(e) {
  cat("\n★★★ scTenifoldKnk 错误 ★★★\n")
  cat(conditionMessage(e), "\n")
  NULL
})

elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
log_info(sprintf("耗时: %.1f 分钟", elapsed))

if (is.null(result)) stop("scTenifoldKnk 失败")

saveRDS(result, file.path(DIR_THIS,
        sprintf("scTenifoldKnk_%s_%s.rds", GENE_KO, CELL_TYPE_KO)))

# ============================================================
# 6. extract DR genes
# ============================================================
log_info("--- 5. 提取 DR 基因 ---")

dr <- result$diffRegulation
dr <- dr[order(-abs(dr$Z)), ]

save_table(dr, DIR_THIS,
           sprintf("%s_diffRegulation_%s.csv", GENE_KO, CELL_TYPE_KO))
save_table(head(dr, 100), DIR_THIS,
           sprintf("%s_top100_%s.csv", GENE_KO, CELL_TYPE_KO))

n_sig <- sum(dr$p.adj < 0.05, na.rm = TRUE)
log_info(sprintf("显著扰动基因 (FDR<0.05): %d / %d", n_sig, nrow(dr)))

cat(sprintf("\n=== %s KO top 30 (%s) ===\n", GENE_KO, CELL_TYPE_KO))
print(head(dr[, c("gene","Z","p.value","p.adj")], 30),
      row.names = FALSE, digits = 4)

# ============================================================
# 7. crossvalidation
# ============================================================
log_info("--- 6. 交叉验证 ---")

ci_res  <- data.frame()
gy_res  <- data.frame()
bl_res  <- data.frame()
sig16_res <- data.frame()

if (length(CI_GENES) > 0) {
  ci_res <- dr[dr$gene %in% CI_GENES, ]
  save_table(ci_res, DIR_THIS,
             sprintf("%s_vs_ComplexI_%s.csv", GENE_KO, CELL_TYPE_KO))
  log_info(sprintf("Complex I 命中 %d, 显著 %d",
                   nrow(ci_res), sum(ci_res$p.adj < 0.05, na.rm = TRUE)))
}

if (length(GY_GENES) > 0) {
  gy_res <- dr[dr$gene %in% GY_GENES, ]
  save_table(gy_res, DIR_THIS,
             sprintf("%s_vs_greenyellow_%s.csv", GENE_KO, CELL_TYPE_KO))
  log_info(sprintf("greenyellow 命中 %d, 显著 %d (%.1f%%)",
                   nrow(gy_res), sum(gy_res$p.adj < 0.05, na.rm = TRUE),
                   100 * mean(gy_res$p.adj < 0.05, na.rm = TRUE)))
}

if (length(BL_GENES) > 0) {
  bl_res <- dr[dr$gene %in% BL_GENES, ]
  save_table(bl_res, DIR_THIS,
             sprintf("%s_vs_black_%s.csv", GENE_KO, CELL_TYPE_KO))
  log_info(sprintf("black 命中 %d, 显著 %d (%.1f%%)",
                   nrow(bl_res), sum(bl_res$p.adj < 0.05, na.rm = TRUE),
                   100 * mean(bl_res$p.adj < 0.05, na.rm = TRUE)))
}

if (length(GENES_16) > 0) {
  sig16_res <- dr[dr$gene %in% GENES_16, ]
  save_table(sig16_res, DIR_THIS,
             sprintf("%s_vs_16genes_%s.csv", GENE_KO, CELL_TYPE_KO))
  log_info(sprintf("16 基因命中 %d, 显著 %d",
                   nrow(sig16_res), sum(sig16_res$p.adj < 0.05, na.rm = TRUE)))
}

# ============================================================
# 8. GSEA
# ============================================================
log_info("--- 7. GSEA ---")

gsea_res <- NULL
if (!is.null(GENE_SETS)) {
  rank_vec <- setNames(dr$Z, dr$gene)
  rank_vec <- sort(rank_vec[is.finite(rank_vec)], decreasing = TRUE)

  gs_list <- lapply(GENE_SETS, function(gs) {
    g <- toupper(gs$genes); g[g %in% names(rank_vec)]
  })
  gs_list <- gs_list[sapply(gs_list, length) >= 5]

  if (length(gs_list) > 0) {
    set.seed(SEED)
    gsea_res <- tryCatch(
      fgsea(pathways = gs_list, stats = rank_vec,
            minSize = 5, maxSize = 500, nperm = 10000),
      error = function(e) NULL)

    if (!is.null(gsea_res) && nrow(gsea_res) > 0) {
      gsea_res <- gsea_res[order(gsea_res$pval), ]
      gsea_res$leadingEdge <- sapply(gsea_res$leadingEdge, paste, collapse = ";")
      save_table(gsea_res, DIR_THIS,
                 sprintf("%s_GSEA_%s.csv", GENE_KO, CELL_TYPE_KO))

      cat("\n=== GSEA top 20 ===\n")
      print(head(gsea_res[, c("pathway","NES","pval","padj","size")], 20),
            row.names = FALSE, digits = 4)
    }
  }
}

# ============================================================
# 9. ORA
# ============================================================
log_info("--- 8. ORA ---")

if (!is.null(GENE_SETS)) {
  sig_genes <- dr$gene[dr$p.adj < 0.05]
  if (length(sig_genes) >= 10) {
    bg <- dr$gene
    ora_res <- data.frame()
    for (gs_name in names(GENE_SETS)) {
      gs <- toupper(GENE_SETS[[gs_name]]$genes)
      gs <- intersect(gs, bg)
      if (length(gs) < 5) next
      hit <- intersect(sig_genes, gs)
      if (length(hit) < 3) next
      p <- phyper(length(hit) - 1, length(gs), length(bg) - length(gs),
                  length(sig_genes), lower.tail = FALSE)
      ora_res <- rbind(ora_res, data.frame(
        pathway = gs_name, n_geneset = length(gs),
        n_hit = length(hit), p = p, stringsAsFactors = FALSE))
    }
    if (nrow(ora_res) > 0) {
      ora_res$FDR <- p.adjust(ora_res$p, "BH")
      ora_res <- ora_res[order(ora_res$p), ]
      save_table(ora_res, DIR_THIS,
                 sprintf("%s_ORA_%s.csv", GENE_KO, CELL_TYPE_KO))
      cat("\n=== ORA top 20 ===\n")
      print(head(ora_res, 20), row.names = FALSE, digits = 4)
    }
  } else {
    log_info("显著扰动基因 < 10，跳过 ORA")
  }
}

# ============================================================
# 10. visualization
# ============================================================
log_info("--- 9. 可视化 ---")

dr$neglog10p <- -log10(pmax(dr$p.value, 1e-300))
dr$highlight <- "other"
dr$highlight[dr$gene %in% CI_GENES] <- "ComplexI"
dr$highlight[dr$gene %in% GY_GENES] <- "greenyellow"
dr$highlight[dr$gene %in% BL_GENES] <- "black"
dr$highlight[dr$gene %in% GENES_16] <- "16genes"

p <- ggplot(dr, aes(x = Z, y = neglog10p, color = highlight)) +
  geom_point(alpha = 0.4, size = 0.7) +
  scale_color_manual(values = c(other = "grey75", ComplexI = "#C44E52",
                                 greenyellow = "#4C72B0",
                                 black = "#55A868",
                                 "16genes" = "#F0AD4E")) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme_bw(base_size = 12) +
  labs(title = sprintf("%s virtual KO (%s, %.1f min)",
                       GENE_KO, CELL_TYPE_KO, elapsed),
       x = "Z score", y = "-log10(p)", color = "")
save_fig(p, DIR_THIS,
         sprintf("%s_volcano_%s", GENE_KO, CELL_TYPE_KO),
         width = 8, height = 5)

if (!is.null(gsea_res) && nrow(gsea_res) > 0) {
  top <- head(gsea_res[gsea_res$padj < 0.25, ], 20)
  if (nrow(top) == 0) top <- head(gsea_res, 20)
  top$pathway <- factor(top$pathway, levels = top$pathway[order(top$NES)])
  p_gsea <- ggplot(top, aes(x = pathway, y = NES, fill = NES > 0)) +
    geom_col() + coord_flip() + theme_bw(base_size = 10) +
    scale_fill_manual(values = c("TRUE" = "#C44E52", "FALSE" = "#4C72B0"),
                       guide = "none") +
    labs(x = "", y = "NES",
         title = sprintf("%s KO: GSEA (%s)", GENE_KO, CELL_TYPE_KO))
  save_fig(p_gsea, DIR_THIS,
           sprintf("%s_GSEA_bar_%s", GENE_KO, CELL_TYPE_KO),
           width = 8, height = 6)
}

pdf(file.path(DIR_THIS, "figures",
              sprintf("%s_KO_network_%s.pdf", GENE_KO, CELL_TYPE_KO)),
    width = 10, height = 10)
tryCatch(plotKO(result, annotate = TRUE),
         error = function(e) plot(1, type = "n", main = "plotKO unavailable"))
dev.off()

# ============================================================
# 11. summary
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat(sprintf("=== 07_3 scTenifoldKnk %s (%s) 汇总 ===\n",
            GENE_KO, CELL_TYPE_KO))
cat(strrep("=", 70), "\n", sep = "")

cat(sprintf("\n【1】数据来源: GSE262151 snRNA_RPE.h5ad\n"))
cat(sprintf("【2】种子基因: %s\n", GENE_KO))
cat(sprintf("【3】细胞类型: %s\n", CELL_TYPE_KO))
cat(sprintf("【4】原始细胞数: %d\n", ncol(sce)))
cat(sprintf("【5】映射后基因数: %d\n", nrow(sce_mapped)))
cat(sprintf("【6】过滤后矩阵: %d genes × %d cells\n",
            nrow(mtx_target_f), ncol(mtx_target_f)))
cat(sprintf("【7】%s 表达: %d / %d (%.2f%%)\n",
            GENE_KO, n_expr_target, ncol(mtx_target_f), pct_expr))
cat(sprintf("【8】运行时间: %.1f 分钟\n", elapsed))
cat(sprintf("【9】显著扰动基因 (FDR<0.05): %d\n", n_sig))
cat(sprintf("【10】top gene: %s (Z = %.2f, p.adj = %.3g)\n",
            dr$gene[1], dr$Z[1], dr$p.adj[1]))

if (length(CI_GENES) > 0)
  cat(sprintf("【11】Complex I 命中: %d, 显著: %d\n",
              nrow(ci_res), sum(ci_res$p.adj < 0.05, na.rm = TRUE)))
if (length(GY_GENES) > 0)
  cat(sprintf("【12】greenyellow 显著: %d (%.1f%%)\n",
              sum(gy_res$p.adj < 0.05, na.rm = TRUE),
              100 * mean(gy_res$p.adj < 0.05, na.rm = TRUE)))
if (length(GENES_16) > 0)
  cat(sprintf("【13】16 基因命中: %d, 显著: %d\n",
              nrow(sig16_res), sum(sig16_res$p.adj < 0.05, na.rm = TRUE)))

if (!is.null(gsea_res) && nrow(gsea_res) > 0) {
  cat("\n【14】GSEA top 5:\n")
  print(head(gsea_res[, c("pathway","NES","padj")], 5),
        row.names = FALSE, digits = 3)
}

cat("\n--- 输出位置 ---\n"); cat("目录:", DIR_THIS, "\n")
log_info(sprintf("=== 07_3 scTenifoldKnk %s (%s) done ===",
                 GENE_KO, CELL_TYPE_KO))
