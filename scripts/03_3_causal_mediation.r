# ============================================================
# 03_3_causal_mediation.r (v5 final)
#
# v5 Scope：
# ============================================================
rm(list = ls(all.names = TRUE)); gc()
source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(limma); library(splines); library(ggplot2)
})

DIR_THIS <- setup_output_dir("03_3_causal_mediation")
log_info("=== 03_3_causal_mediation.r (v5 final) start ===")

# ============================================================
# ============================================================
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
log_info("表达矩阵: {nrow(expr)} × {ncol(expr)}")

# ============================================================
# ============================================================
med_raw <- NULL
med_source <- NA_character_
f1 <- file.path(DIR_RESULT, "03_1_celltype_composition", "data",
                "celltype_scores.rds")
if (is.null(med_raw) && file.exists(f1)) {
  s <- readRDS(f1); med_raw <- t(s)
  med_source <- "03_1/celltype_scores.rds"
}
f2 <- file.path(DIR_RESULT, "03_2_deconvolution", "data",
                "cell_proportions.rds")
if (is.null(med_raw) && file.exists(f2)) {
  p <- readRDS(f2); med_raw <- as.matrix(p)
  med_source <- "03_2/cell_proportions.rds"
}
if (is.null(med_raw)) stop("找不到中介变量文件")
log_info("中介变量来源: {med_source}")
log_info("中介变量: {nrow(med_raw)} × {ncol(med_raw)}")

# ============================================================
# 3. samplesalign
# ============================================================
common2 <- Reduce(intersect, list(rownames(med_raw), colnames(expr)))
med_raw <- med_raw[common2, , drop = FALSE]
expr    <- expr[, common2, drop = FALSE]
clin    <- clin[common2, , drop = FALSE]
log_info("对齐后样本: {length(common2)}")

# ============================================================
# 4. mediationvariable PCA
# ============================================================
med_pca <- prcomp(med_raw, scale. = TRUE)
n_pc <- min(5, ncol(med_pca$x))
med_use <- med_pca$x[, 1:n_pc, drop = FALSE]
colnames(med_use) <- paste0("cellPC", 1:n_pc)
var_exp <- round(100 * med_pca$sdev[1:n_pc]^2 / sum(med_pca$sdev^2), 1)
log_info("cellPC 方差解释: {paste(var_exp, '%', collapse=', ')}")

# ============================================================
# ============================================================
sp <- splines::ns(clin$age_n, df = 3)
covar_df <- data.frame(
  age_s1 = sp[,1], age_s2 = sp[,2], age_s3 = sp[,3],
  sex_M  = as.integer(clin$sex == "Male"),
  rin    = clin$rin_n,
  pmi    = clin$pmi_n
)
rownames(covar_df) <- rownames(clin)

design_df <- cbind(data.frame(grade = clin$grade), med_use, covar_df)
design_df <- design_df[complete.cases(design_df), , drop = FALSE]
expr_sub  <- expr[, rownames(design_df), drop = FALSE]

poly_mat <- poly(design_df$grade, 3)
p3_component <- poly_mat[, 3]

M <- as.matrix(design_df[, colnames(med_use), drop = FALSE])
C <- as.matrix(design_df[, colnames(covar_df), drop = FALSE])
Y <- t(expr_sub)
n_med <- ncol(M); n_genes <- ncol(Y)
log_info("样本: {nrow(design_df)}, 基因: {n_genes}, 中介: {n_med} cellPC")

# ============================================================
# ============================================================
fit_stats <- function(X, Y_mat, coef_idx) {
  n <- nrow(X); p <- ncol(X)
  qr_X <- qr(X)
  beta <- qr.coef(qr_X, Y_mat)
  resid <- Y_mat - X %*% beta
  sigma2 <- colSums(resid^2) / (n - p)
  XtX_inv <- chol2inv(chol(crossprod(X)))
  se_diag <- diag(XtX_inv)[coef_idx]
  se_mat  <- sqrt(outer(se_diag, sigma2, "*"))
  t_mat   <- beta[coef_idx, , drop = FALSE] / se_mat
  p_mat   <- 2 * pt(-abs(t_mat), df = n - p)
  list(beta = beta, se = se_mat, t = t_mat, p = p_mat)
}

# ============================================================
# ============================================================
log_info("--- A→M 诊断 ---")
X_A <- cbind(1, p3_component, C)
fit_M <- fit_stats(X_A, M, 2)
a_k <- fit_M$beta[2, ]
a_p <- fit_M$p[1, ]

am_tbl <- data.frame(
  mediator = colnames(M),
  a_coef = round(a_k, 4),
  a_p = round(a_p, 4),
  a_FDR = round(p.adjust(a_p, "BH"), 4),
  stringsAsFactors = FALSE
)
save_table(am_tbl, DIR_THIS, "A_to_M_diagnostic.csv")

cat("\n=== A → M 诊断 ===\n")
print(am_tbl, row.names = FALSE)

am_cor <- data.frame(
  mediator = colnames(M),
  pearson_r = round(cor(p3_component, M), 3),
  stringsAsFactors = FALSE
)
save_table(am_cor, DIR_THIS, "A_to_M_correlation.csv")
cat("\nA 与 M 的 Pearson 相关:\n")
print(am_cor, row.names = FALSE)

n_sig_am <- sum(a_p < 0.05)
A_to_M_ok <- n_sig_am >= 1
log_info("A→M 显著数: {n_sig_am} / {n_med}")
if (!A_to_M_ok) {
  log_warn("★ A→M 不显著，中介分析前提不成立")
} else {
  log_info("A→M 部分显著，中介分析可以继续")
}

# ============================================================
# ============================================================
log_info("--- TE 模型 ---")
fit_TE <- fit_stats(X_A, Y, 2)
TE   <- fit_TE$beta[2, ]
TE_p <- fit_TE$p[1, ]
TE_FDR <- p.adjust(TE_p, "BH")

te_tbl <- data.frame(
  gene = rownames(expr),
  TE = TE, TE_p = TE_p, TE_FDR = TE_FDR,
  stringsAsFactors = FALSE
)
save_table(te_tbl, DIR_THIS, "TE_full.csv")

n_TE_sig <- sum(TE_FDR < 0.05)
log_info("TE FDR<0.05: {n_TE_sig} / {n_genes}")

# ============================================================
# ============================================================
log_info("--- 中介分解（前提不成立，仅供参考）---")

X_Y <- cbind(1, p3_component, M, C)
m_idx <- 3:(2 + n_med)
fit_Y <- fit_stats(X_Y, Y, c(2, m_idx))
NDE   <- fit_Y$beta[2, ]
NDE_p <- fit_Y$p[1, ]
b_k   <- fit_Y$beta[m_idx, , drop = FALSE]
se_b  <- fit_Y$se[2:(1 + n_med), , drop = FALSE]

se_a <- fit_M$se[1, ]
NIE_k <- sweep(b_k, 1, a_k, "*")
NIE   <- colSums(NIE_k)
SE_NIE_k <- sqrt(sweep(se_b^2, 1, a_k^2, "*") +
                 sweep(b_k^2, 1, se_a^2, "*"))
SE_NIE <- sqrt(colSums(SE_NIE_k))
z_NIE <- NIE / SE_NIE
p_NIE <- 2 * pnorm(-abs(z_NIE))

FDR_NDE <- p.adjust(NDE_p, "BH")
FDR_NIE <- p.adjust(p_NIE, "BH")
prop_med <- NIE / TE
prop_med[!is.finite(prop_med)] <- NA

mediation_tbl <- data.frame(
  gene = rownames(expr),
  TE = TE, TE_p = TE_p, TE_FDR = TE_FDR,
  NDE = NDE, NDE_p = NDE_p, NDE_FDR = FDR_NDE,
  NIE = NIE, NIE_SE = SE_NIE, NIE_p = p_NIE, NIE_FDR = FDR_NIE,
  prop_mediated = prop_med,
  stringsAsFactors = FALSE
)
mediation_tbl$category <- "none"
mediation_tbl$category[mediation_tbl$TE_FDR < 0.05 &
                       mediation_tbl$NDE_FDR >= 0.05 &
                       mediation_tbl$NIE_FDR < 0.05] <- "fully_mediated"
mediation_tbl$category[mediation_tbl$TE_FDR < 0.05 &
                       mediation_tbl$NDE_FDR < 0.05 &
                       mediation_tbl$NIE_FDR < 0.05] <- "partially_mediated"
mediation_tbl$category[mediation_tbl$TE_FDR < 0.05 &
                       mediation_tbl$NDE_FDR < 0.05 &
                       mediation_tbl$NIE_FDR >= 0.05] <- "direct_only"
mediation_tbl$category[mediation_tbl$TE_FDR >= 0.05 &
                       mediation_tbl$NDE_FDR < 0.05] <- "suppressed"

save_table(mediation_tbl, DIR_THIS, "mediation_full.csv")

cat_tbl <- as.data.frame(table(mediation_tbl$category))
colnames(cat_tbl) <- c("category", "n")
cat_tbl <- cat_tbl[order(-cat_tbl$n), ]
save_table(cat_tbl, DIR_THIS, "mediation_summary.csv")

log_info("TE FDR<0.05    : {sum(mediation_tbl$TE_FDR < 0.05)}")
log_info("NDE FDR<0.05   : {sum(mediation_tbl$NDE_FDR < 0.05)}")
log_info("NIE FDR<0.05   : {sum(mediation_tbl$NIE_FDR < 0.05)}")

# ============================================================
# ============================================================
KEY_GENES <- c("HTRA1","ARMS2","CFH","C3","CFB","APOE","VEGFA","BCL11A",
                "NOTCH1","NOTCH2","HIF1A","EPAS1","TGFB1","IL6","TNF",
                "S100B","BTBD19","CCDC3","NDUFB8","UGGT1")
kg <- mediation_tbl[mediation_tbl$gene %in% KEY_GENES, ]
if (nrow(kg) > 0) {
  save_table(kg, DIR_THIS, "keygene_effects.csv")
  cat("\n=== 关键基因效应分解 ===\n")
  print(kg[, c("gene","TE","TE_FDR","NDE","NDE_p",
                "NIE","NIE_p","prop_mediated","category")],
        row.names = FALSE, digits = 4)
}

# ============================================================
# ============================================================
log_info("--- 图 ---")

plot_df <- mediation_tbl
plot_df$sig_type <- "ns"
plot_df$sig_type[plot_df$NDE_FDR < 0.05] <- "direct"
plot_df$sig_type[plot_df$category == "fully_mediated"] <- "fully_mediated"
plot_df$sig_type[plot_df$category == "suppressed"] <- "suppressed"

p1 <- ggplot(plot_df, aes(x = TE, y = NDE, color = sig_type)) +
  geom_point(alpha = 0.3, size = 0.5) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = 0, color = "gray70") +
  geom_vline(xintercept = 0, color = "gray70") +
  scale_color_manual(values = c(ns = "gray80", direct = "#C44E52",
                                 fully_mediated = "#4C72B0",
                                 suppressed = "#55A868")) +
  theme_bw(base_size = 12) +
  labs(x = "Total effect TE (p3)", y = "Direct effect NDE (p3)",
       color = "", title = "TE vs NDE (p3 component)")
save_fig(p1, DIR_THIS, "TE_vs_NDE", width = 7, height = 6)

# A→M diagnosisfigure
p2 <- ggplot(am_tbl, aes(x = reorder(mediator, -abs(a_coef)),
                          y = a_coef, fill = a_p < 0.05)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(aes(label = sprintf("p=%.2g", a_p)),
            vjust = ifelse(am_tbl$a_coef > 0, -0.5, 1.2), size = 3) +
  scale_fill_manual(values = c("TRUE" = "#C44E52", "FALSE" = "#B0B0B0"),
                    guide = "none") +
  theme_bw(base_size = 11) +
  labs(x = "", y = "A→M coefficient (p3 → cellPC)",
       title = "A→M path: p3 component to cell composition PCs")
save_fig(p2, DIR_THIS, "A_to_M_diagnostic", width = 7, height = 5)

top20_te <- head(te_tbl[order(te_tbl$TE_p), ], 20)
top20_te$gene <- factor(top20_te$gene,
                         levels = top20_te$gene[order(top20_te$TE)])
p3 <- ggplot(top20_te, aes(x = gene, y = TE, fill = TE > 0)) +
  geom_col() + coord_flip() + theme_bw(base_size = 10) +
  scale_fill_manual(values = c("TRUE" = "#C44E52", "FALSE" = "#4C72B0"),
                    guide = "none") +
  labs(x = "", y = "Total effect TE (p3)",
       title = "Top 20 genes by p3 total effect")
save_fig(p3, DIR_THIS, "TE_top20", width = 7, height = 6)

# ============================================================
# 12. summary
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("=== 03_3 v5 final 汇总 ===\n")
cat(strrep("=", 70), "\n", sep = "")
cat("中介变量来源:", med_source, "\n")
cat("中介变量:", n_med, "cellPC\n")
cat("样本数:", nrow(design_df), "  基因数:", n_genes, "\n\n")

cat("--- A→M 诊断（前提检查）---\n")
print(am_tbl, row.names = FALSE)
cat("\nA→M 显著数:", n_sig_am, "/", n_med, "\n")
if (!A_to_M_ok) {
  cat("★ A→M 不显著 → 中介分析前提不成立\n")
  cat("   阴性结果含义：p3 相变不改变细胞组成 → 支持 cell-intrinsic\n")
}

cat("\n--- TE（主结果）---\n")
cat("TE FDR<0.05:", n_TE_sig, "/", n_genes, "\n")

cat("\n--- 中介分解（前提不成立，仅供记录）---\n")
cat("NDE FDR<0.05:", sum(mediation_tbl$NDE_FDR < 0.05), "\n")
cat("NIE FDR<0.05:", sum(mediation_tbl$NIE_FDR < 0.05), "\n")
print(cat_tbl, row.names = FALSE)

if (nrow(kg) > 0) {
  cat("\n--- 关键基因效应分解 ---\n")
  print(kg[, c("gene","TE","TE_FDR","NDE","NDE_p",
                "NIE","NIE_p","prop_mediated","category")],
        row.names = FALSE, digits = 4)
}

cat("\n输出目录:", DIR_THIS, "\n")
log_info("=== 03_3 v5 final done ===")

