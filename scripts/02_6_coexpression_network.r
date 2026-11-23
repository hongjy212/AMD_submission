# ============================================================
# 02_6_coexpression_network.r
# Input：
#   results/01_1_prepare_bulk_data/data/{expr_corrected,clinical,gene_annotation}.rds
#   results/02_1_differential_expression_main/tables/sig16_genes.csv
# ============================================================

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(limma); library(splines); library(ggplot2)
  library(pheatmap)
})

DIR_THIS <- setup_output_dir("02_6_coexpression_network")
log_info("=== 02_6_coexpression_network.r start ===")

R_THRESHOLD   <- 0.5
FDR_THRESHOLD <- 0.01

# ============================================================
# ============================================================
log_info("--- 1. 载入数据 ---")

expr <- read_upstream("01_1_prepare_bulk_data", "data", "expr_corrected.rds")
clin <- read_upstream("01_1_prepare_bulk_data", "data", "clinical.rds")
anno <- read_upstream("01_1_prepare_bulk_data", "data", "gene_annotation.rds")

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
if (!file.exists(sig16_file)) {
  stop("请先运行 02_1_differential_expression_main.r", call. = FALSE)
}
sig16 <- read.csv(sig16_file, stringsAsFactors = FALSE)$gene
log_info("读入 16 基因: {paste(sig16, collapse=', ')}")
log_info("表达矩阵: {nrow(expr)} × {ncol(expr)}")

# ============================================================
# ============================================================
log_info("--- 2. 残差化 ---")

sp <- splines::ns(clin$age_n, df = 3)
covar_df <- data.frame(
  age_s1 = sp[,1], age_s2 = sp[,2], age_s3 = sp[,3],
  sex_M = as.integer(clin$sex == "Male"),
  rin = clin$rin_n, pmi = clin$pmi_n
)
rownames(covar_df) <- rownames(clin)
design_cov <- model.matrix(
  ~ grade + age_s1 + age_s2 + age_s3 + sex_M + rin + pmi,
  data = data.frame(grade = clin$grade, covar_df))
expr_sub <- expr[, rownames(covar_df), drop = FALSE]
fit_cov <- lmFit(expr_sub, design_cov)
resid_expr <- expr_sub - fit_cov$coefficients %*% t(design_cov)
log_info("残差矩阵: {nrow(resid_expr)} × {ncol(resid_expr)}")

# ============================================================
# ============================================================
log_info("--- 3. 共表达网络 ---")
summary_tbl <- data.frame()
all_pos_genes <- list()
all_neg_genes <- list()

for (g in sig16) {
  log_info("共表达 {g} ...")

  g_v <- resid_expr[g, ]
  n <- length(g_v)
  r_all <- cor(t(resid_expr), g_v, method = "pearson")[, 1]
  t_stat <- r_all * sqrt((n - 2) / (1 - r_all^2))
  p_all <- 2 * pt(-abs(t_stat), df = n - 2)
  FDR_all <- p.adjust(p_all, "BH")

  coexp <- data.frame(gene = rownames(resid_expr),
                      r = r_all, p = p_all, FDR = FDR_all,
                      stringsAsFactors = FALSE)
  coexp <- coexp[coexp$gene != g, ]
  coexp <- coexp[order(-abs(coexp$r)), ]
  save_table(coexp, DIR_THIS, paste0("coexp_", g, "_full.csv"))

  pos <- coexp$gene[coexp$FDR < FDR_THRESHOLD & coexp$r > R_THRESHOLD]
  neg <- coexp$gene[coexp$FDR < FDR_THRESHOLD & coexp$r < -R_THRESHOLD]
  all_pos_genes[[g]] <- pos
  all_neg_genes[[g]] <- neg

  summary_tbl <- rbind(summary_tbl, data.frame(
    gene = g,
    SD = round(sd(resid_expr[g, ]), 4),
    n_pos = length(pos), n_neg = length(neg),
    stringsAsFactors = FALSE))
}

summary_tbl <- summary_tbl[order(-summary_tbl$n_pos), ]
save_table(summary_tbl, DIR_THIS, "all16_summary.csv")

# ============================================================
# ============================================================
log_info("--- 4. 网络相似度 ---")

jmat <- matrix(NA, length(sig16), length(sig16),
                dimnames = list(sig16, sig16))
for (i in seq_along(sig16)) {
  for (j in seq_along(sig16)) {
    a <- all_pos_genes[[sig16[i]]]; b <- all_pos_genes[[sig16[j]]]
    u <- length(union(a, b))
    jmat[i, j] <- if (u == 0) 0 else length(intersect(a, b)) / u
  }
}
save_table(data.frame(
  gene1 = rownames(jmat)[row(jmat)],
  gene2 = colnames(jmat)[col(jmat)],
  jaccard = round(as.vector(jmat), 3)),
  DIR_THIS, "pos_network_jaccard.csv")

# ============================================================
# ============================================================
log_info("--- 5. 共识正相关基因 ---")

all_pos_flat <- unlist(all_pos_genes)
pos_freq <- sort(table(all_pos_flat), decreasing = TRUE)
pos_freq_tbl <- data.frame(
  gene = names(pos_freq),
  n_queries = as.integer(pos_freq),
  stringsAsFactors = FALSE)
save_table(pos_freq_tbl, DIR_THIS, "pos_gene_frequency.csv")

core_pos <- pos_freq_tbl$gene[pos_freq_tbl$n_queries >= 8]
save_table(data.frame(gene = core_pos), DIR_THIS, "core_pos_genes.csv")
log_info("≥8 query 核心正相关: {length(core_pos)}")

# ============================================================
# 6. SD vs n_pos diagnosis
# ============================================================
log_info("--- 6. SD vs n_pos ---")

coexp_check <- summary_tbl
coexp_check$mean_expr <- rowMeans(expr[sig16, ])
coexp_check$mean_resid <- rowMeans(resid_expr[sig16, ])

cat("\n=== SD vs n_pos ===\n")
print(coexp_check, row.names = FALSE)
cat("\ncor(SD, n_pos) =",
    round(cor(coexp_check$SD, coexp_check$n_pos), 3), "\n")
cat("cor(mean_expr, n_pos) =",
    round(cor(coexp_check$mean_expr, coexp_check$n_pos), 3), "\n")

# ============================================================
# 7. figure
# ============================================================
log_info("--- 7. 图 ---")

sum_long <- rbind(
  data.frame(gene = summary_tbl$gene, n = summary_tbl$n_pos,
             direction = "pos", stringsAsFactors = FALSE),
  data.frame(gene = summary_tbl$gene, n = summary_tbl$n_neg,
             direction = "neg", stringsAsFactors = FALSE))
p1 <- ggplot(sum_long, aes(x = reorder(gene, n), y = n, fill = direction)) +
  geom_col(position = "dodge") + coord_flip() + theme_bw(base_size = 11) +
  scale_fill_manual(values = c(pos = "#C44E52", neg = "#4C72B0")) +
  labs(x = "", y = "Coexpressed genes (|r|>0.5, FDR<0.01)", fill = "",
       title = "16 gene coexpression network size")
save_fig(p1, DIR_THIS, "all16_n_coexp", width = 7, height = 6)

pheatmap(jmat, cluster_rows = TRUE, cluster_cols = TRUE,
         display_numbers = TRUE, number_format = "%.2f",
         main = "16 gene network Jaccard",
         filename = file.path(DIR_THIS, "figures", "all16_jaccard.png"),
         width = 9, height = 8)

# ============================================================
# 8. summary
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("=== 02_6 共表达网络 汇总 ===\n")
cat(strrep("=", 70), "\n", sep = "")
print(summary_tbl, row.names = FALSE)
cat("\ncor(SD, n_pos) =", round(cor(summary_tbl$SD, summary_tbl$n_pos), 3), "\n")
cat("≥8 query 核心:", length(core_pos), "\n")

cat("\n--- 输出位置 ---\n")
cat("目录:", DIR_THIS, "\n")
cat("  tables/coexp_<gene>_full.csv × 16\n")
cat("  tables/all16_summary.csv\n")
cat("  tables/pos_network_jaccard.csv\n")
cat("  tables/pos_gene_frequency.csv\n")
cat("  tables/core_pos_genes.csv\n")
cat("  figures/all16_n_coexp.pdf/.png\n")
cat("  figures/all16_jaccard.pdf/.png\n")

log_info("=== 02_6_coexpression_network.r done ===")
