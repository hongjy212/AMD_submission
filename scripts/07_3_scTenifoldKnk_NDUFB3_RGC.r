# ============================================================
# 07_3_scTenifoldKnk_NDUFB3_RGC.r
#
#   Rscript 07_3_scTenifoldKnk_NDUFB3_RGC.r                    # default NDUFB3
# ============================================================

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(scTenifoldKnk)
  library(Matrix); library(dplyr); library(tidyr)
  library(ggplot2); library(fgsea)
})

# ============================================================
# ============================================================
DIR_SC  <- "E:/retina_tools/raw/GSE221042"
N_CORES <- 8
SEED    <- 42

args <- commandArgs(trailingOnly = TRUE)
GENE_KO <- if (length(args) >= 1) {
  args[1]
} else if (nzchar(Sys.getenv("KO_GENE"))) {
  Sys.getenv("KO_GENE")
} else {
  "NDUFB3"
}

DIR_THIS <- setup_output_dir(sprintf("07_3_scTenifoldKnk_%s_RGC", GENE_KO))
log_info(sprintf("=== 07_3 scTenifoldKnk %s (RGC) start ===", GENE_KO))

# QC Parameters
QC_MIN_LIBSIZE <- 500
QC_MAX_MTRATIO <- 0.10
QC_MIN_PCT     <- 0.02

NC_NNET     <- 5         # 网络数
NC_NCELLS   <- 300       # 每网络子采样细胞数
NC_NCOMP    <- 3
TD_K        <- 3         # 张量秩
TD_MAXITER  <- 500       # 最大迭代
TD_MAXERROR <- 1e-4      # 收敛阈值（比标准放宽 10 倍，比极限严 10 倍）

# genesfilter
MIN_GENE_CELLS <- 15

# ============================================================
# ============================================================
log_info("--- 0. 读取外部基因集 ---")

GENE_SETS <- read_upstream("01_2_prepare_gene_sets", "data", "gene_sets.rds")

CI_GENES <- toupper(GENE_SETS[["Mito_CI_Subunits"]]$genes)
log_info(sprintf("Mito_CI_Subunits: %d 基因", length(CI_GENES)))

sig16_file <- file.path(DIR_RESULT, "02_1_differential_expression_main",
                         "tables", "sig16_genes.csv")
GENES_16 <- if (file.exists(sig16_file)) {
  toupper(read.csv(sig16_file)$gene)
} else {
  log_warn("sig16_genes.csv 不存在")
  character(0)
}
log_info(sprintf("16 基因签名: %d", length(GENES_16)))

# greenyellow module
gy_file <- file.path(DIR_RESULT, "05_1_consensus_wgcna", "tables",
                      "module_genes_greenyellow_symbolfix.csv")
GY_GENES <- if (file.exists(gy_file)) {
  g <- toupper(read.csv(gy_file)$SYMBOL)
  g[!is.na(g) & g != ""]
} else {
  log_warn("greenyellow 模块文件不存在")
  character(0)
}
log_info(sprintf("greenyellow: %d", length(GY_GENES)))

# black module
bl_file <- file.path(DIR_RESULT, "05_1_consensus_wgcna", "tables",
                      "module_genes_black_symbolfix.csv")
BL_GENES <- if (file.exists(bl_file)) {
  g <- toupper(read.csv(bl_file)$SYMBOL)
  g[!is.na(g) & g != ""]
} else {
  character(0)
}
log_info(sprintf("black: %d", length(BL_GENES)))

# ============================================================
# 2. read GSE221042
# ============================================================
log_info("--- 1. 读入 GSE221042 ---")

mtx_raw  <- readMM(file.path(DIR_SC, "GSE221042_data.mtx.gz"))
genes    <- read.csv(file.path(DIR_SC, "GSE221042_genes.csv.gz"),
                     header = FALSE, stringsAsFactors = FALSE)
meta_raw <- read.csv(file.path(DIR_SC, "GSE221042_sample_annotations.csv.gz"),
                     check.names = FALSE, stringsAsFactors = FALSE)

gene_symbol <- sub(" \\(.*", "", genes[[2]])
gene_symbol <- make.unique(toupper(trimws(gene_symbol)))

mtx <- t(mtx_raw)
rownames(mtx) <- gene_symbol
colnames(mtx) <- make.unique(as.character(meta_raw[[1]]))
meta <- meta_raw
rownames(meta) <- colnames(mtx)

log_info(sprintf("矩阵: %d genes x %d cells", nrow(mtx), ncol(mtx)))

# ============================================================
# 3. RGC annotation
# ============================================================
log_info("--- 2. RGC 注释 ---")

MARKERS <- list(
  Rod  = c("RHO","GNAT1","PDE6B","NRL"),
  Cone = c("ARR3","GNAT2","PDE6H"),
  HC   = c("ONECUT1","ONECUT2","LHX1","TFAP2B"),
  BC   = c("VSX1","VSX2","GRM6"),
  AC   = c("GAD1","GAD2"),
  MG   = c("RLBP1","GLUL"),
  RGC  = c("RBPMS","NEFL","NEFM"),
  RPE  = c("RPE65","BEST1","TTR"))

libsize_all <- Matrix::colSums(mtx)
norm_all    <- log1p(t(t(mtx) / libsize_all * 1e4))
gene_ids    <- rownames(norm_all)

scores <- do.call(cbind, lapply(MARKERS, function(g) {
  r <- which(gene_ids %in% g)
  if (length(r) == 0) return(rep(0, ncol(norm_all)))
  Matrix::colMeans(norm_all[r, , drop = FALSE])
}))
colnames(scores) <- names(MARKERS)
meta$cell_type <- colnames(scores)[apply(scores, 1, which.max)]

log_info("cell_type 分布:")
print(table(meta$cell_type))

rgc_cells <- rownames(meta)[meta$cell_type == "RGC"]
log_info(sprintf("RGC 细胞数: %d", length(rgc_cells)))

mtx_rgc <- mtx[, rgc_cells, drop = FALSE]

# ============================================================
# ============================================================
log_info("--- 3. 基因过滤 ---")

n_cells_per_gene <- Matrix::rowSums(mtx_rgc > 0)
keep <- n_cells_per_gene >= MIN_GENE_CELLS

force_keep <- unique(c(GENE_KO, CI_GENES, GENES_16, GY_GENES, BL_GENES))
force_keep <- force_keep[force_keep %in% rownames(mtx_rgc)]
keep[rownames(mtx_rgc) %in% force_keep] <- TRUE

mtx_rgc_f <- mtx_rgc[keep, , drop = FALSE]
log_info(sprintf("基因过滤: %d → %d", nrow(mtx_rgc), nrow(mtx_rgc_f)))

if (!GENE_KO %in% rownames(mtx_rgc_f)) {
  stop(sprintf("种子基因 %s 不在过滤后矩阵中，请检查 MIN_GENE_CELLS = %d",
               GENE_KO, MIN_GENE_CELLS))
}

n_expr_rgc <- sum(mtx_rgc_f[GENE_KO, ] > 0)
pct_expr   <- 100 * n_expr_rgc / ncol(mtx_rgc_f)
log_info(sprintf("RGC 中 %s: %d / %d (%.2f%%)",
                 GENE_KO, n_expr_rgc, ncol(mtx_rgc_f), pct_expr))

if (n_expr_rgc < 30) {
  stop(sprintf("RGC 中 %s 表达细胞数 (%d) 太少，无法做 KO", GENE_KO, n_expr_rgc))
}

# ============================================================
# 5. scTenifoldKnk
# ============================================================
log_info(sprintf("--- 4. scTenifoldKnk (%s, RGC) ---", GENE_KO))
log_info(sprintf("参数: nc_nNet=%d, nc_nCells=%d, td_K=%d, td_maxError=%g",
                 NC_NNET, NC_NCELLS, TD_K, TD_MAXERROR))

set.seed(SEED)
t0 <- Sys.time()

result <- tryCatch({
  scTenifoldKnk(
    countMatrix     = mtx_rgc_f,
    gKO             = GENE_KO,
    transcriptomeWide = FALSE,
    qc              = TRUE,
    qc_minLibSize   = QC_MIN_LIBSIZE,
    qc_removeOutlierCells = TRUE,
    qc_minPCT       = QC_MIN_PCT,
    qc_maxMTratio   = QC_MAX_MTRATIO,
    nc_lambda       = 0,
    nc_nNet         = NC_NNET,
    nc_nCells       = NC_NCELLS,
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

saveRDS(result, file.path(DIR_THIS, sprintf("scTenifoldKnk_%s_RGC.rds", GENE_KO)))
log_info("saved: scTenifoldKnk_RGC.rds")

# ============================================================
# 6. extract DR genes
# ============================================================
log_info("--- 5. 提取 DR 基因 ---")

dr <- result$diffRegulation
dr <- dr[order(-abs(dr$Z)), ]

save_table(dr, DIR_THIS, sprintf("%s_diffRegulation_RGC.csv", GENE_KO))
save_table(head(dr, 100), DIR_THIS, sprintf("%s_top100_RGC.csv", GENE_KO))

n_sig <- sum(dr$p.adj < 0.05, na.rm = TRUE)
log_info(sprintf("显著扰动基因 (FDR<0.05): %d / %d", n_sig, nrow(dr)))

cat(sprintf("\n=== %s KO top 30 (RGC) ===\n", GENE_KO))
print(head(dr[, c("gene","Z","p.value","p.adj")], 30),
      row.names = FALSE, digits = 4)

# ============================================================
# ============================================================
log_info("--- 6. 交叉验证 ---")

if (length(CI_GENES) > 0) {
  ci_res <- dr[dr$gene %in% CI_GENES, ]
  save_table(ci_res, DIR_THIS, sprintf("%s_vs_ComplexI_RGC.csv", GENE_KO))
  log_info(sprintf("Complex I 命中 %d, 显著 %d",
                   nrow(ci_res), sum(ci_res$p.adj < 0.05, na.rm = TRUE)))
}

if (length(GY_GENES) > 0) {
  gy_res <- dr[dr$gene %in% GY_GENES, ]
  save_table(gy_res, DIR_THIS, sprintf("%s_vs_greenyellow_RGC.csv", GENE_KO))
  log_info(sprintf("greenyellow 命中 %d, 显著 %d (%.1f%%)",
                   nrow(gy_res), sum(gy_res$p.adj < 0.05, na.rm = TRUE),
                   100 * mean(gy_res$p.adj < 0.05, na.rm = TRUE)))
}

if (length(BL_GENES) > 0) {
  bl_res <- dr[dr$gene %in% BL_GENES, ]
  save_table(bl_res, DIR_THIS, sprintf("%s_vs_black_RGC.csv", GENE_KO))
  log_info(sprintf("black 命中 %d, 显著 %d (%.1f%%)",
                   nrow(bl_res), sum(bl_res$p.adj < 0.05, na.rm = TRUE),
                   100 * mean(bl_res$p.adj < 0.05, na.rm = TRUE)))
}

if (length(GENES_16) > 0) {
  sig16_res <- dr[dr$gene %in% GENES_16, ]
  save_table(sig16_res, DIR_THIS, sprintf("%s_vs_16genes_RGC.csv", GENE_KO))
  log_info(sprintf("16 基因命中 %d", nrow(sig16_res)))
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
      save_table(gsea_res, DIR_THIS, sprintf("%s_GSEA_RGC.csv", GENE_KO))

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
      save_table(ora_res, DIR_THIS, sprintf("%s_ORA_RGC.csv", GENE_KO))
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
  labs(title = sprintf("%s virtual KO (RGC, %.1f min)", GENE_KO, elapsed),
       x = "Z score", y = "-log10(p)", color = "")
save_fig(p, DIR_THIS, sprintf("%s_volcano_RGC", GENE_KO), width = 8, height = 5)

if (!is.null(gsea_res) && nrow(gsea_res) > 0) {
  top <- head(gsea_res[gsea_res$padj < 0.25, ], 20)
  if (nrow(top) == 0) top <- head(gsea_res, 20)
  top$pathway <- factor(top$pathway, levels = top$pathway[order(top$NES)])
  p_gsea <- ggplot(top, aes(x = pathway, y = NES, fill = NES > 0)) +
    geom_col() + coord_flip() + theme_bw(base_size = 10) +
    scale_fill_manual(values = c("TRUE" = "#C44E52", "FALSE" = "#4C72B0"),
                       guide = "none") +
    labs(x = "", y = "NES",
         title = sprintf("%s KO: GSEA (RGC)", GENE_KO))
  save_fig(p_gsea, DIR_THIS, sprintf("%s_GSEA_bar_RGC", GENE_KO),
           width = 8, height = 6)
}

# KO networkfigure
pdf(file.path(DIR_THIS, "figures", sprintf("%s_KO_network_RGC.pdf", GENE_KO)),
    width = 10, height = 10)
tryCatch(plotKO(result, annotate = TRUE),
         error = function(e) plot(1, type = "n", main = "plotKO unavailable"))
dev.off()

# ============================================================
# 11. summary
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat(sprintf("=== 07_3 scTenifoldKnk %s (RGC) 汇总 ===\n", GENE_KO))
cat(strrep("=", 70), "\n", sep = "")

cat(sprintf("\n【1】种子基因: %s\n", GENE_KO))
cat(sprintf("【2】RGC 细胞数: %d\n", length(rgc_cells)))
cat(sprintf("【3】过滤后基因数: %d\n", nrow(mtx_rgc_f)))
cat(sprintf("【4】%s 表达: %d / %d (%.2f%%)\n",
            GENE_KO, n_expr_rgc, ncol(mtx_rgc_f), pct_expr))
cat(sprintf("【5】运行时间: %.1f 分钟\n", elapsed))
cat(sprintf("【6】显著扰动基因: %d\n", n_sig))
cat(sprintf("【7】top gene: %s (Z = %.2f, p.adj = %.3g)\n",
            dr$gene[1], dr$Z[1], dr$p.adj[1]))

if (length(CI_GENES) > 0)
  cat(sprintf("【8】Complex I 显著: %d\n",
              sum(ci_res$p.adj < 0.05, na.rm = TRUE)))
if (length(GY_GENES) > 0)
  cat(sprintf("【9】greenyellow 显著: %d (%.1f%%)\n",
              sum(gy_res$p.adj < 0.05, na.rm = TRUE),
              100 * mean(gy_res$p.adj < 0.05, na.rm = TRUE)))
if (length(GENES_16) > 0)
  cat(sprintf("【10】16 基因命中: %d\n", nrow(sig16_res)))

if (!is.null(gsea_res) && nrow(gsea_res) > 0) {
  cat("\n【11】GSEA top 5:\n")
  print(head(gsea_res[, c("pathway","NES","padj")], 5),
        row.names = FALSE, digits = 3)
}

cat("\n--- 输出位置 ---\n"); cat("目录:", DIR_THIS, "\n")
log_info(sprintf("=== 07_3 scTenifoldKnk %s (RGC) done ===", GENE_KO))
