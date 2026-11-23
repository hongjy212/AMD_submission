# ============================================================
# 02_7_functional_enrichment.r
# Input：
#   results/01_1_prepare_bulk_data/data/{expr_corrected,clinical,gene_annotation}.rds
#   results/01_2_prepare_gene_sets/data/gene_sets.rds
#   results/02_6_coexpression_network/tables/coexp_<gene>_full.csv × 16
#   tables/ora_<gene>_top20.csv       —— eachgenes top20 ORA
#   tables/ora_<gene>_top50.csv       —— eachgenes top50 ORA
#   tables/gsea_<gene>.csv            —— eachgenes GSEA
#   tables/all16_enrichment_summary.csv
#   tables/btbd19_pos.csv / btbd19_neg.csv
#   tables/ora_btbd19_pos.csv / ora_btbd19_neg.csv
#   figures/S100B_top20_bar.pdf/.png
#   figures/BTBD19_top30_bar.pdf/.png
# ============================================================

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(limma); library(splines); library(ggplot2)
  library(tidyr)
})

DIR_THIS <- setup_output_dir("02_7_functional_enrichment")
log_info("=== 02_7_functional_enrichment.r start ===")

# ============================================================
# ============================================================
log_info("--- 1. 载入数据 ---")

expr <- read_upstream("01_1_prepare_bulk_data", "data", "expr_corrected.rds")
clin <- read_upstream("01_1_prepare_bulk_data", "data", "clinical.rds")
anno <- read_upstream("01_1_prepare_bulk_data", "data", "gene_annotation.rds")
GENE_SETS <- read_upstream("01_2_prepare_gene_sets", "data", "gene_sets.rds")

log_info("gene sets (01_2): {length(GENE_SETS)} 个")

gs_inventory <- data.frame(
  geneset = names(GENE_SETS),
  group = sapply(GENE_SETS, function(x) x$group %||% NA),
  source = sapply(GENE_SETS, function(x) x$source %||% NA),
  n_genes = sapply(GENE_SETS, function(x) length(x$genes)),
  stringsAsFactors = FALSE)
save_table(gs_inventory, DIR_THIS, "gene_sets_inventory.csv")
print(gs_inventory, row.names = FALSE)

common <- intersect(colnames(expr), rownames(clin))
expr <- expr[, common, drop = FALSE]
clin <- clin[common, , drop = FALSE]

anno_sub <- anno[anno$ensembl %in% rownames(expr), ]
anno_sub <- anno_sub[!is.na(anno_sub$symbol), ]
expr <- expr[anno_sub$ensembl, , drop = FALSE]
expr <- expr[!duplicated(anno_sub$symbol), , drop = FALSE]
rownames(expr) <- anno_sub$symbol[!duplicated(anno_sub$symbol)]

clin$grade <- as.integer(clin$grade)
rownames(clin) <- colnames(expr)
UNIVERSE <- rownames(expr)

sig16_file <- file.path(DIR_RESULT,
                         "02_1_differential_expression_main",
                         "tables", "sig16_genes.csv")
sig16 <- read.csv(sig16_file, stringsAsFactors = FALSE)$gene
log_info("16 基因: {length(sig16)} 个")

log_info("表达矩阵: {nrow(expr)} × {ncol(expr)}")

# ============================================================
# ============================================================
`%||%` <- function(a, b) if (is.null(a)) b else a

run_ora <- function(genes, universe, gene_sets,
                    min_hit = 1, min_panel = 5) {
  genes <- intersect(genes, universe)
  n_sig <- length(genes)
  if (n_sig < 3) return(NULL)
  out <- data.frame()
  for (gs in names(gene_sets)) {
    panel <- intersect(gene_sets[[gs]]$genes, universe)
    if (length(panel) < min_panel) next
    hits <- intersect(genes, panel)
    if (length(hits) < min_hit) next
    p <- phyper(length(hits) - 1, length(panel),
                length(universe) - length(panel), n_sig,
                lower.tail = FALSE)
    expected <- n_sig * length(panel) / length(universe)
    out <- rbind(out, data.frame(
      geneset = gs, group = gene_sets[[gs]]$group,
      n_hit = length(hits), n_sig = n_sig, n_panel = length(panel),
      fold = round(length(hits) / max(expected, 1e-6), 2),
      p = p, hits = paste(hits, collapse = ","),
      stringsAsFactors = FALSE))
  }
  if (nrow(out) == 0) return(NULL)
  out$FDR <- p.adjust(out$p, "BH")
  out[order(out$p), ]
}

run_gsea <- function(ranks, gene_sets, minSize = 10, maxSize = 500,
                      nperm = 5000) {
  if (!requireNamespace("fgsea", quietly = TRUE)) return(NULL)
  library(fgsea)
  pathways <- lapply(gene_sets, function(gs) intersect(gs$genes, names(ranks)))
  pathways <- pathways[sapply(pathways, length) >= minSize &
                       sapply(pathways, length) <= maxSize]
  if (length(pathways) == 0) return(NULL)
  set.seed(SEED)
  res <- tryCatch(
    fgsea(pathways = pathways, stats = ranks,
          minSize = minSize, maxSize = maxSize, nperm = nperm),
    error = function(e) NULL)
  if (is.null(res)) return(NULL)
  res <- as.data.frame(res)
  res$leadingEdge <- NULL
  res[order(res$padj), ]
}

# ============================================================
# ============================================================
log_info("--- 2. 读共表达结果 ---")

DIR_02_6 <- file.path(DIR_RESULT, "02_6_coexpression_network", "tables")
coexp_list <- list()
for (g in sig16) {
  f <- file.path(DIR_02_6, paste0("coexp_", g, "_full.csv"))
  if (file.exists(f)) {
    coexp_list[[g]] <- read.csv(f, stringsAsFactors = FALSE)
  }
}
log_info("读入 {length(coexp_list)} 个共表达结果")

# ============================================================
# ============================================================
log_info("--- 3. 逐一 ORA + GSEA ---")

ora_summary <- data.frame()

for (g in names(coexp_list)) {
  log_info("富集 {g} ...")
  ce <- coexp_list[[g]]
  ce <- ce[order(-abs(ce$r)), ]

  top20 <- head(ce, 20)$gene
  top50 <- head(ce, 50)$gene

  ora20 <- run_ora(top20, UNIVERSE, GENE_SETS)
  ora50 <- run_ora(top50, UNIVERSE, GENE_SETS)

  if (!is.null(ora20)) {
    ora20$query <- g; ora20$set <- "top20"
    save_table(ora20, DIR_THIS, paste0("ora_", g, "_top20.csv"))
  }
  if (!is.null(ora50)) {
    ora50$query <- g; ora50$set <- "top50"
    save_table(ora50, DIR_THIS, paste0("ora_", g, "_top50.csv"))
  }

  ranks <- ce$r
  names(ranks) <- ce$gene
  ranks <- sort(ranks, decreasing = TRUE)
  gs_res <- run_gsea(ranks, GENE_SETS)
  if (!is.null(gs_res)) {
    gs_res$query <- g
    save_table(gs_res, DIR_THIS, paste0("gsea_", g, ".csv"))
  }

  n_ora_sig <- if (!is.null(ora50)) sum(ora50$FDR < 0.05) else 0
  n_gsea_sig <- if (!is.null(gs_res)) sum(gs_res$padj < 0.05) else 0
  top_path <- if (!is.null(ora50) && n_ora_sig > 0) {
    paste(head(ora50$geneset[ora50$FDR < 0.05], 3), collapse = "; ")
  } else ""

  ora_summary <- rbind(ora_summary, data.frame(
    gene = g,
    n_pos = sum(ce$FDR < 0.01 & ce$r > 0.5),
    n_neg = sum(ce$FDR < 0.01 & ce$r < -0.5),
    n_ora_fdr05 = n_ora_sig,
    n_gsea_fdr05 = n_gsea_sig,
    top_pathways = top_path,
    stringsAsFactors = FALSE))
}

ora_summary <- ora_summary[order(-ora_summary$n_gsea_fdr05,
                                  -ora_summary$n_ora_fdr05), ]
save_table(ora_summary, DIR_THIS, "all16_enrichment_summary.csv")

cat("\n=== 16 基因富集汇总 ===\n")
print(ora_summary, row.names = FALSE)

# ============================================================
# ============================================================
log_info("--- 4. BTBD19 正/负 ORA ---")

if ("BTBD19" %in% names(coexp_list)) {
  btbd_ce <- coexp_list[["BTBD19"]]
  btbd_pos <- btbd_ce$gene[btbd_ce$r > 0.5 & btbd_ce$FDR < 0.01]
  btbd_neg <- btbd_ce$gene[btbd_ce$r < -0.5 & btbd_ce$FDR < 0.01]

  save_table(data.frame(gene = btbd_pos), DIR_THIS, "btbd19_pos.csv")
  save_table(data.frame(gene = btbd_neg), DIR_THIS, "btbd19_neg.csv")

  ora_btbd_pos <- run_ora(btbd_pos, UNIVERSE, GENE_SETS)
  ora_btbd_neg <- run_ora(btbd_neg, UNIVERSE, GENE_SETS)

  if (!is.null(ora_btbd_pos)) {
    save_table(ora_btbd_pos, DIR_THIS, "ora_btbd19_pos.csv")
    cat("\n=== BTBD19 正相关 top 5 ===\n")
    print(head(ora_btbd_pos[, c("geneset","group","n_hit","fold","FDR")], 5),
          row.names = FALSE)
  }
  if (!is.null(ora_btbd_neg)) {
    save_table(ora_btbd_neg, DIR_THIS, "ora_btbd19_neg.csv")
    cat("\n=== BTBD19 负相关 top 5 ===\n")
    print(head(ora_btbd_neg[, c("geneset","group","n_hit","fold","FDR")], 5),
          row.names = FALSE)
  }
}

# ============================================================
# 6. figure
# ============================================================
log_info("--- 5. 图 ---")

if ("S100B" %in% names(coexp_list)) {
  ce <- coexp_list[["S100B"]]
  top20 <- head(ce[order(-abs(ce$r)), ], 20)
  top20$gene <- factor(top20$gene, levels = rev(top20$gene))
  p <- ggplot(top20, aes(x = gene, y = r, fill = r > 0)) +
    geom_col() + coord_flip() + theme_bw(base_size = 11) +
    scale_fill_manual(values = c("TRUE" = "#C44E52", "FALSE" = "#4C72B0"),
                       labels = c("Neg", "Pos"), guide = "none") +
    labs(x = "", y = "Pearson r with S100B",
         title = "S100B top 20 coexpressed genes")
  save_fig(p, DIR_THIS, "S100B_top20_bar", width = 7, height = 5)
}

if ("BTBD19" %in% names(coexp_list)) {
  ce <- coexp_list[["BTBD19"]]
  top30 <- head(ce[order(-abs(ce$r)), ], 30)
  top30$gene <- factor(top30$gene, levels = rev(top30$gene))
  p <- ggplot(top30, aes(x = gene, y = r, fill = r > 0)) +
    geom_col() + coord_flip() + theme_bw(base_size = 10) +
    scale_fill_manual(values = c("TRUE" = "#C44E52", "FALSE" = "#4C72B0"),
                       labels = c("Neg", "Pos"), name = "") +
    labs(x = "", y = "Pearson r with BTBD19",
         title = "BTBD19 top 30 coexpressed genes")
  save_fig(p, DIR_THIS, "BTBD19_top30_bar", width = 7, height = 6)
}

# ============================================================
# 7. summary
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("=== 02_7 功能富集 汇总 ===\n")
cat(strrep("=", 70), "\n", sep = "")

cat("\nGene sets（来自 01_2）:", length(GENE_SETS), "\n\n")

cat("=== 16 基因富集汇总 ===\n")
print(ora_summary, row.names = FALSE)

cat("\n--- 输出位置 ---\n")
cat("目录:", DIR_THIS, "\n")
cat("  tables/gene_sets_inventory.csv\n")
cat("  tables/ora_<gene>_top20.csv × 16\n")
cat("  tables/ora_<gene>_top50.csv × 16\n")
cat("  tables/gsea_<gene>.csv × 16\n")
cat("  tables/all16_enrichment_summary.csv\n")
cat("  tables/btbd19_pos.csv / btbd19_neg.csv\n")
cat("  tables/ora_btbd19_pos.csv / ora_btbd19_neg.csv\n")
cat("  figures/S100B_top20_bar.pdf/.png\n")
cat("  figures/BTBD19_top30_bar.pdf/.png\n")

log_info("=== 02_7_functional_enrichment.r done ===")

