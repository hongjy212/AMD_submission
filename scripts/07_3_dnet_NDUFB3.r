# ============================================================
# 07_3_dnet_NDUFB3.r
# ============================================================

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(Matrix); library(dplyr); library(tidyr)
  library(ggplot2); library(fgsea)
  library(dnet); library(igraph)
})

DIR_THIS <- setup_output_dir("07_3_dnet_NDUFB3")
log_info("=== 07_3 dnet RWR NDUFB3 (full) start ===")

# ============================================================
# 0. Parameters
# ============================================================
DIR_SC    <- "E:/retina_tools/raw/GSE221042"
MIN_CELLS <- 10
TOP_K     <- 100
RESTART   <- 0.75

args <- commandArgs(trailingOnly = TRUE)
GENE_SEED <- if (length(args) >= 1) args[1] else "NDUFB3"
log_info(sprintf("种子基因: %s", GENE_SEED))

GENE_SETS <- read_upstream("01_2_prepare_gene_sets", "data", "gene_sets.rds")
CI_GENES  <- toupper(GENE_SETS[["Mito_CI_Subunits"]]$genes)

sig16_file <- file.path(DIR_RESULT, "02_1_differential_expression_main",
                         "tables", "sig16_genes.csv")
GENES_16 <- if (file.exists(sig16_file)) toupper(read.csv(sig16_file)$gene) else character(0)

gy_file <- file.path(DIR_RESULT, "05_1_consensus_wgcna", "tables",
                      "module_genes_greenyellow_symbolfix.csv")
GY_GENES <- if (file.exists(gy_file)) {
  g <- toupper(read.csv(gy_file)$SYMBOL); g[!is.na(g) & g != ""]
} else character(0)

log_info(sprintf("Mito_CI_Subunits: %d", length(CI_GENES)))
log_info(sprintf("16 基因: %d", length(GENES_16)))
log_info(sprintf("greenyellow: %d", length(GY_GENES)))

# ============================================================
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

log_info(sprintf("矩阵: %d genes x %d cells", nrow(mtx), ncol(mtx)))

# ============================================================
# 2. genesfilter
# ============================================================
log_info("--- 2. 基因过滤 ---")

n_cells_per_gene <- Matrix::rowSums(mtx > 0)
keep <- n_cells_per_gene >= MIN_CELLS

force_keep <- unique(c(GENE_SEED, CI_GENES, GENES_16, GY_GENES))
force_keep <- force_keep[force_keep %in% rownames(mtx)]
keep[rownames(mtx) %in% force_keep] <- TRUE

mtx_f <- mtx[keep, , drop = FALSE]
log_info(sprintf("基因过滤: %d → %d", nrow(mtx), nrow(mtx_f)))

stopifnot(GENE_SEED %in% rownames(mtx_f))
n_seed_expr <- sum(mtx_f[GENE_SEED, ] > 0)
log_info(sprintf("%s 表达: %d / %d (%.2f%%)",
                 GENE_SEED, n_seed_expr, ncol(mtx_f),
                 100 * n_seed_expr / ncol(mtx_f)))

# ============================================================
# 3. normalization
# ============================================================
log_info("--- 3. 归一化 ---")

libsize <- Matrix::colSums(mtx_f)
norm <- log1p(t(t(mtx_f) / libsize * 1e4))
norm <- as.matrix(norm)
log_info(sprintf("归一化矩阵: %d x %d", nrow(norm), ncol(norm)))

# ============================================================
# ============================================================
log_info("--- 4. 计算共表达网络 ---")

calc_cor_topK <- function(mat, top_k, chunk_size = 500) {
  n_genes <- nrow(mat)
  gene_names <- rownames(mat)
  neighbors <- vector("list", n_genes)
  weights   <- vector("list", n_genes)
  mat_c <- mat - rowMeans(mat)

  for (i in seq(1, n_genes, by = chunk_size)) {
    idx_end <- min(i + chunk_size - 1, n_genes)
    idx <- i:idx_end
    chunk <- mat_c[idx, , drop = FALSE]
    cor_mat <- (chunk %*% t(mat_c)) / ncol(mat)
    for (j in seq_along(idx)) {
      r <- cor_mat[j, ]
      r[idx[j]] <- -Inf
      top_idx <- order(-abs(r))[1:top_k]
      neighbors[[idx[j]]] <- gene_names[top_idx]
      weights[[idx[j]]]   <- r[top_idx]
    }
    if (i %% 2000 == 0) log_info(sprintf("  完成 %d / %d", idx_end, n_genes))
  }
  list(neighbors = neighbors, weights = weights)
}

t0 <- Sys.time()
net <- calc_cor_topK(norm, TOP_K)
elapsed_cor <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
log_info(sprintf("共表达网络完成，耗时 %.1f 分钟", elapsed_cor))

# ============================================================
# 5. build igraph
# ============================================================
log_info("--- 5. 构建 igraph ---")

n_genes    <- nrow(norm)
gene_names <- rownames(norm)

rows <- rep(seq_len(n_genes), each = TOP_K)
cols <- unlist(lapply(net$neighbors, function(x) match(x, gene_names)))
vals <- unlist(net$weights)

ok <- !is.na(cols) & is.finite(vals) & vals > 0
A <- sparseMatrix(i = rows[ok], j = cols[ok], x = vals[ok],
                   dims = c(n_genes, n_genes))
A <- pmax(A, t(A))

g <- graph_from_adjacency_matrix(A, mode = "undirected",
                                  weighted = TRUE, diag = FALSE)
log_info(sprintf("igraph: %d nodes, %d edges", vcount(g), ecount(g)))

# ============================================================
# 6. dRWR
# ============================================================
log_info("--- 6. dRWR ---")

stopifnot(GENE_SEED %in% V(g)$name)

setSeeds <- matrix(0, nrow = vcount(g), ncol = 1)
rownames(setSeeds) <- V(g)$name
setSeeds[GENE_SEED, 1] <- 1

t0 <- Sys.time()
rwr_result <- dRWR(g = g,
                    normalise = "laplacian",
                    setSeeds = setSeeds,
                    restart = RESTART,
                    parallel = TRUE,
                    multicores = 8)
elapsed_rwr <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
log_info(sprintf("dRWR 完成，耗时 %.1f 秒", elapsed_rwr))

rwr_vec <- as.numeric(rwr_result)
names(rwr_vec) <- rownames(rwr_result)

prop_df <- data.frame(gene = names(rwr_vec),
                       propagation_score = rwr_vec,
                       stringsAsFactors = FALSE)
prop_df <- prop_df[order(-prop_df$propagation_score), ]
prop_df$rank <- seq_len(nrow(prop_df))

save_table(prop_df, DIR_THIS,
            sprintf("%s_dnet_propagation_scores.csv", GENE_SEED))
save_table(head(prop_df, 100), DIR_THIS,
            sprintf("%s_dnet_top100.csv", GENE_SEED))

cat(sprintf("\n=== %s dnet 传播 top 30 ===\n", GENE_SEED))
print(head(prop_df, 30), row.names = FALSE, digits = 6)

# ============================================================
# 7. crossvalidation
# ============================================================
log_info("--- 7. 交叉验证 ---")

if (length(CI_GENES) > 0) {
  ci_res <- prop_df[prop_df$gene %in% CI_GENES, ]
  save_table(ci_res, DIR_THIS, sprintf("%s_vs_ComplexI.csv", GENE_SEED))
  log_info(sprintf("Complex I 命中 %d, 中位排名 %d",
                   nrow(ci_res), median(ci_res$rank)))
}

if (length(GY_GENES) > 0) {
  gy_res <- prop_df[prop_df$gene %in% GY_GENES, ]
  save_table(gy_res, DIR_THIS, sprintf("%s_vs_greenyellow.csv", GENE_SEED))
  log_info(sprintf("greenyellow 命中 %d, 中位排名 %d",
                   nrow(gy_res), median(gy_res$rank)))
}

if (length(GENES_16) > 0) {
  sig16_res <- prop_df[prop_df$gene %in% GENES_16, ]
  save_table(sig16_res, DIR_THIS, sprintf("%s_vs_16genes.csv", GENE_SEED))
  log_info(sprintf("16 基因命中 %d", nrow(sig16_res)))
}

# ============================================================
# 8. GSEA
# ============================================================
log_info("--- 8. GSEA ---")

gsea_res <- NULL
if (!is.null(GENE_SETS)) {
  rank_vec <- setNames(log10(prop_df$propagation_score + 1e-10), prop_df$gene)
  rank_vec <- sort(rank_vec, decreasing = TRUE)

  gs_list <- lapply(GENE_SETS, function(gs) {
    g <- toupper(gs$genes); g[g %in% names(rank_vec)]
  })
  gs_list <- gs_list[sapply(gs_list, length) >= 5]

  if (length(gs_list) > 0) {
    set.seed(42)
    gsea_res <- tryCatch(
      fgsea(pathways = gs_list, stats = rank_vec,
            minSize = 5, maxSize = 500, nperm = 10000),
      error = function(e) NULL)

    if (!is.null(gsea_res) && nrow(gsea_res) > 0) {
      gsea_res <- gsea_res[order(gsea_res$pval), ]
      gsea_res$leadingEdge <- sapply(gsea_res$leadingEdge, paste, collapse = ";")
      save_table(gsea_res, DIR_THIS, sprintf("%s_GSEA.csv", GENE_SEED))

      cat("\n=== GSEA top 20 ===\n")
      print(head(gsea_res[, c("pathway","NES","pval","padj","size")], 20),
            row.names = FALSE, digits = 4)
    }
  }
}

# ============================================================
# 9. ORA
# ============================================================
log_info("--- 9. ORA ---")

if (!is.null(GENE_SETS)) {
  n_top <- ceiling(n_genes * 0.05)
  sig_genes <- head(prop_df$gene, n_top)
  bg <- prop_df$gene

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
    save_table(ora_res, DIR_THIS, sprintf("%s_ORA.csv", GENE_SEED))
    cat("\n=== ORA top 20 ===\n")
    print(head(ora_res, 20), row.names = FALSE, digits = 4)
  }
}

# ============================================================
# 10. visualization
# ============================================================
log_info("--- 10. 可视化 ---")

prop_df$neglog10 <- -log10(prop_df$propagation_score + 1e-10)
prop_df$highlight <- "other"
prop_df$highlight[prop_df$gene %in% CI_GENES] <- "ComplexI"
prop_df$highlight[prop_df$gene %in% GY_GENES] <- "greenyellow"
prop_df$highlight[prop_df$gene %in% GENES_16] <- "16genes"

p <- ggplot(prop_df[prop_df$rank <= 5000, ],
            aes(x = rank, y = neglog10, color = highlight)) +
  geom_point(alpha = 0.4, size = 0.7) +
  scale_color_manual(values = c(other = "grey75", ComplexI = "#C44E52",
                                 greenyellow = "#4C72B0",
                                 "16genes" = "#55A868")) +
  theme_bw(base_size = 12) +
  labs(title = sprintf("%s dnet propagation (top 5000)", GENE_SEED),
       x = "Rank", y = "-log10(propagation score)", color = "")
save_fig(p, DIR_THIS, sprintf("%s_dnet_propagation_rank", GENE_SEED),
         width = 10, height = 5)

if (!is.null(gsea_res) && nrow(gsea_res) > 0) {
  top <- head(gsea_res[gsea_res$padj < 0.25, ], 20)
  if (nrow(top) == 0) top <- head(gsea_res, 20)
  top$pathway <- factor(top$pathway, levels = top$pathway[order(top$NES)])
  p_gsea <- ggplot(top, aes(x = pathway, y = NES, fill = NES > 0)) +
    geom_col() + coord_flip() + theme_bw(base_size = 10) +
    scale_fill_manual(values = c("TRUE" = "#C44E52", "FALSE" = "#4C72B0"),
                       guide = "none") +
    labs(x = "", y = "NES",
         title = sprintf("%s dnet propagation: GSEA", GENE_SEED))
  save_fig(p_gsea, DIR_THIS, sprintf("%s_dnet_GSEA_bar", GENE_SEED),
           width = 8, height = 6)
}

# ============================================================
# 11. summary
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat(sprintf("=== %s dnet propagation 汇总 ===\n", GENE_SEED))
cat(strrep("=", 70), "\n", sep = "")

cat(sprintf("\n【1】种子基因: %s\n", GENE_SEED))
cat(sprintf("【2】细胞数: %d\n", ncol(norm)))
cat(sprintf("【3】过滤后基因数: %d\n", n_genes))
cat(sprintf("【4】%s 表达: %d (%.2f%%)\n",
            GENE_SEED, n_seed_expr, 100 * n_seed_expr / ncol(norm)))
cat(sprintf("【5】网络参数: top_k = %d, restart = %.2f\n", TOP_K, RESTART))
cat(sprintf("【6】共表达耗时: %.1f 分钟, dRWR 耗时: %.1f 秒\n",
            elapsed_cor, elapsed_rwr))

cat("\n【7】传播 top 10:\n")
print(head(prop_df[, c("gene","propagation_score","rank")], 10),
      row.names = FALSE, digits = 6)

if (length(CI_GENES) > 0)
  cat(sprintf("\n【8】Complex I 中位排名: %d\n", median(ci_res$rank)))
if (length(GY_GENES) > 0)
  cat(sprintf("【9】greenyellow 中位排名: %d\n", median(gy_res$rank)))
if (length(GENES_16) > 0)
  cat(sprintf("【10】16 基因中位排名: %d\n", median(sig16_res$rank)))

if (!is.null(gsea_res) && nrow(gsea_res) > 0) {
  cat("\n【11】GSEA top 5:\n")
  print(head(gsea_res[, c("pathway","NES","padj")], 5),
        row.names = FALSE, digits = 3)
}

cat("\n--- 输出位置 ---\n"); cat("目录:", DIR_THIS, "\n")
log_info(sprintf("=== %s dnet propagation done ===", GENE_SEED))
