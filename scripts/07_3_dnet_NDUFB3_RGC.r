# ============================================================
# 07_3_dnet_NDUFB3_RGC.r
# ============================================================

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(Matrix); library(dplyr); library(tidyr)
  library(ggplot2); library(fgsea); library(dnet); library(igraph)
})

DIR_THIS <- setup_output_dir("07_3_dnet_NDUFB3_RGC")
log_info("=== 07_3 dnet RWR NDUFB3 (RGC) start ===")

# ============================================================
# ============================================================
DIR_SC      <- "E:/retina_tools/raw/GSE221042"
MIN_CELLS   <- 10
TOP_K       <- 100
RESTART     <- 0.75   # dnet 推荐的默认重启概率

args <- commandArgs(trailingOnly = TRUE)
GENE_SEED <- if (length(args) >= 1) args[1] else "NDUFB3"

GENE_SETS <- read_upstream("01_2_prepare_gene_sets", "data", "gene_sets.rds")
CI_GENES  <- toupper(GENE_SETS[["Mito_CI_Subunits"]]$genes)

sig16_file <- file.path(DIR_RESULT, "02_1_differential_expression_main", "tables", "sig16_genes.csv")
GENES_16 <- if (file.exists(sig16_file)) toupper(read.csv(sig16_file)$gene) else character(0)

gy_file <- file.path(DIR_RESULT, "05_1_consensus_wgcna", "tables", "module_genes_greenyellow_symbolfix.csv")
GY_GENES <- if (file.exists(gy_file)) {
  g <- toupper(read.csv(gy_file)$SYMBOL); g[!is.na(g) & g != ""]
} else character(0)

# ============================================================
# ============================================================


g <- graph_from_adjacency_matrix(W, mode = "undirected", weighted = TRUE, diag = FALSE)

# ============================================================
# ============================================================
log_info("--- 使用 dRWR 进行网络传播 ---")

setSeeds <- matrix(0, nrow = vcount(g), ncol = 1)
rownames(setSeeds) <- V(g)$name
if (!GENE_SEED %in% V(g)$name) stop("种子基因不在网络中")
setSeeds[GENE_SEED, 1] <- 1

rwr_scores <- dRWR(g = g, 
                   normalise = "laplacian", 
                   setSeeds = setSeeds, 
                   restart = RESTART,
                   parallel = TRUE,   # 启用并行
                   multicores = 8)    # 根据你的CPU调整

rwr_vec <- as.numeric(rwr_scores)
names(rwr_vec) <- rownames(rwr_scores)

# ============================================================
# ============================================================
prop_df <- data.frame(
  gene = names(rwr_vec),
  propagation_score = rwr_vec,
  stringsAsFactors = FALSE
)
prop_df <- prop_df[order(-prop_df$propagation_score), ]
prop_df$rank <- seq_len(nrow(prop_df))

save_table(prop_df, DIR_THIS, sprintf("%s_dnet_propagation_scores_RGC.csv", GENE_SEED))
save_table(head(prop_df, 100), DIR_THIS, sprintf("%s_dnet_top100_RGC.csv", GENE_SEED))

cat(sprintf("\n=== %s dnet 传播 top 30 (RGC) ===\n", GENE_SEED))
print(head(prop_df, 30), row.names = FALSE, digits = 6)

# ============================================================
# ============================================================
