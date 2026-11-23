# ============================================================
# 02_4_single_gene_method_sweep.r
# ============================================================

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(limma); library(edgeR); library(MASS); library(splines)
  library(ggplot2); library(tidyr)
})

DIR_THIS <- setup_output_dir("02_4_single_gene_method_sweep")
log_info("=== 02_4_single_gene_method_sweep.r start ===")

# ============================================================
# ============================================================
log_info("--- 1. 载入数据 ---")

expr <- read_upstream("01_1_prepare_bulk_data", "data", "expr_corrected.rds")
counts <- read_upstream("01_1_prepare_bulk_data", "data", "expr_count.rds")
clin <- read_upstream("01_1_prepare_bulk_data", "data", "clinical.rds")
anno <- read_upstream("01_1_prepare_bulk_data", "data", "gene_annotation.rds")

# samplesalign
common <- intersect(colnames(expr), rownames(clin))
expr <- expr[, common, drop = FALSE]
counts <- counts[, common, drop = FALSE]
clin <- clin[common, , drop = FALSE]

anno_sub <- anno[anno$ensembl %in% rownames(expr), ]
anno_sub <- anno_sub[!is.na(anno_sub$symbol), ]
anno_sub <- anno_sub[!duplicated(anno_sub$symbol), ]

ensembl_keep <- anno_sub$ensembl      # ← 保留 ensembl 索引
symbol_keep  <- anno_sub$symbol       # ← 保留 symbol

counts <- counts[ensembl_keep, , drop = FALSE]

expr <- expr[ensembl_keep, , drop = FALSE]
rownames(expr) <- symbol_keep
rownames(counts) <- symbol_keep       # ★ counts 也 symbol 化

counts[is.na(counts)] <- 0
counts <- round(counts)

log_info("表达矩阵: {nrow(expr)} × {ncol(expr)}")
log_info("counts 对齐: {nrow(counts)} × {ncol(counts)}")
log_info("counts 范围: {min(counts)} ~ {max(counts)}")

clin$grade <- as.integer(clin$grade)
rownames(clin) <- colnames(expr)

keep <- clin$grade %in% c(1, 4)
expr_41 <- expr[, keep, drop = FALSE]
counts_41 <- counts[, keep, drop = FALSE]
clin_41 <- clin[keep, , drop = FALSE]
clin_41$group <- factor(ifelse(clin_41$grade == 4, "g4", "g1"),
                         levels = c("g1", "g4"))

log_info("g4_vs_g1: g1={sum(clin_41$grade==1)}, g4={sum(clin_41$grade==4)}")

# covariatematrix
sp <- splines::ns(clin_41$age_n, df = 3)
covar_41 <- data.frame(
  age_s1 = sp[,1], age_s2 = sp[,2], age_s3 = sp[,3],
  sex_M = as.integer(clin_41$sex == "Male"),
  rin = clin_41$rin_n, pmi = clin_41$pmi_n
)
design <- model.matrix(~ group + age_s1 + age_s2 + age_s3 + sex_M + rin + pmi,
                       data = data.frame(group = clin_41$group, covar_41))

# ============================================================
# 2. evaluatefunction
# ============================================================
eval_result <- function(res, method, effect_col = "logFC") {
  res <- res[!is.na(res$p), ]
  res$FDR <- p.adjust(res$p, "BH")
  sig <- res[res$FDR < 0.05, ]
  data.frame(
    method = method, n_total = nrow(res),
    n_p05 = sum(res$p < 0.05, na.rm = TRUE),
    n_fdr05 = nrow(sig),
    n_fdr10 = sum(res$FDR < 0.10, na.rm = TRUE),
    median_effect = round(median(abs(sig[[effect_col]]), na.rm = TRUE), 4),
    min_p = signif(min(res$p, na.rm = TRUE), 3),
    stringsAsFactors = FALSE)
}

all_results <- list()
summary_tbl <- data.frame()

# ============================================================
# 3. limma eBayes
# ============================================================
log_info("--- 方法 1: limma eBayes ---")
fit <- eBayes(lmFit(expr_41, design))
tt <- topTable(fit, coef = "groupg4", number = Inf, sort.by = "none")
res1 <- data.frame(gene = rownames(tt), logFC = tt$logFC, p = tt$P.Value,
                   stringsAsFactors = FALSE)
all_results$limma_eBayes <- res1
summary_tbl <- rbind(summary_tbl, eval_result(res1, "limma_eBayes"))

# ============================================================
# 4. limma robust
# ============================================================
log_info("--- 方法 2: limma robust ---")
fit_r <- tryCatch(eBayes(lmFit(expr_41, design, method = "robust")),
                   error = function(e) NULL)
if (!is.null(fit_r)) {
  tt_r <- topTable(fit_r, coef = "groupg4", number = Inf, sort.by = "none")
  res2 <- data.frame(gene = rownames(tt_r), logFC = tt_r$logFC, p = tt_r$P.Value,
                     stringsAsFactors = FALSE)
  all_results$limma_robust <- res2
  summary_tbl <- rbind(summary_tbl, eval_result(res2, "limma_robust"))
}

# ============================================================
# 5. limma-trend
# ============================================================
log_info("--- 方法 3: limma-trend ---")
log_cpm <- edgeR::cpm(counts_41, log = TRUE, prior.count = 2)
keep_expr <- rowSums(log_cpm > 0) >= 0.5 * ncol(log_cpm)
log_cpm <- log_cpm[keep_expr, , drop = FALSE]
fit_t <- eBayes(lmFit(log_cpm, design, trend = TRUE))
tt_t <- topTable(fit_t, coef = "groupg4", number = Inf, sort.by = "none")
res3 <- data.frame(gene = rownames(tt_t), logFC = tt_t$logFC, p = tt_t$P.Value,
                   stringsAsFactors = FALSE)
all_results$limma_trend <- res3
summary_tbl <- rbind(summary_tbl, eval_result(res3, "limma_trend"))

# ============================================================
# 6. voomWithQualityWeights
# ============================================================
log_info("--- 方法 4: voom quality ---")
dge <- DGEList(counts = counts_41)
dge <- edgeR::calcNormFactors(dge, method = "TMM")
vw <- voomWithQualityWeights(dge, design = design, plot = FALSE)
fit_q <- eBayes(lmFit(vw, design))
tt_q <- topTable(fit_q, coef = "groupg4", number = Inf, sort.by = "none")
res4 <- data.frame(gene = rownames(tt_q), logFC = tt_q$logFC, p = tt_q$P.Value,
                   stringsAsFactors = FALSE)
all_results$voom_quality <- res4
summary_tbl <- rbind(summary_tbl, eval_result(res4, "voom_quality"))

# ============================================================
# 7. edgeR QLF
# ============================================================
log_info("--- 方法 5: edgeR QLF ---")
dge2 <- DGEList(counts = counts_41)
dge2 <- edgeR::calcNormFactors(dge2, method = "TMM")
dge2 <- estimateDisp(dge2, design, robust = TRUE)
fit_e <- glmQLFit(dge2, design, robust = TRUE)
qlf <- glmQLFTest(fit_e, coef = "groupg4")
tt_e <- topTags(qlf, n = Inf, sort.by = "none")$table
res5 <- data.frame(gene = rownames(tt_e), logFC = tt_e$logFC, p = tt_e$PValue,
                   stringsAsFactors = FALSE)
all_results$edgeR_QLF <- res5
summary_tbl <- rbind(summary_tbl, eval_result(res5, "edgeR_QLF"))

# ============================================================
# 8. Wilcoxon
# ============================================================
log_info("--- 方法 6: Wilcoxon ---")
design_cov <- model.matrix(~ age_s1 + age_s2 + age_s3 + sex_M + rin + pmi,
                            data = covar_41)
fit_cov <- lmFit(expr_41, design_cov)
resid_expr <- expr_41 - fit_cov$coefficients %*% t(design_cov)

p_wilcox <- numeric(nrow(resid_expr))
for (i in seq_len(nrow(resid_expr))) {
  p_wilcox[i] <- tryCatch(
    wilcox.test(resid_expr[i, ] ~ clin_41$group)$p.value,
    error = function(e) NA_real_)
}
res7 <- data.frame(gene = rownames(resid_expr),
                    logFC = rowMeans(resid_expr[, clin_41$group == "g4"]) -
                            rowMeans(resid_expr[, clin_41$group == "g1"]),
                    p = p_wilcox, stringsAsFactors = FALSE)
all_results$Wilcoxon <- res7
summary_tbl <- rbind(summary_tbl, eval_result(res7, "Wilcoxon"))

# ============================================================
# 9. rlm
# ============================================================
log_info("--- 方法 7: rlm ---")
design_rlm <- data.frame(group_g4 = as.integer(clin_41$group == "g4"), covar_41)
beta_rlm <- numeric(nrow(expr_41))
p_rlm <- numeric(nrow(expr_41))
t0 <- Sys.time()
for (i in seq_len(nrow(expr_41))) {
  df <- design_rlm
  df$y <- expr_41[i, ]
  fit <- tryCatch(MASS::rlm(y ~ ., data = df, maxit = 50),
                   error = function(e) NULL)
  if (is.null(fit)) { beta_rlm[i] <- NA; p_rlm[i] <- NA; next }
  co <- summary(fit)$coefficients
  if (!"group_g4" %in% rownames(co)) { beta_rlm[i] <- NA; p_rlm[i] <- NA; next }
  beta_rlm[i] <- co["group_g4", "Value"]
  p_rlm[i] <- 2 * pt(-abs(co["group_g4", "t value"]), df = nrow(df) - ncol(df))
  if (i %% 3000 == 0) {
    dt <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
    log_info("  rlm {i}/{nrow(expr_41)} ({dt}s)")
  }
}
res9 <- data.frame(gene = rownames(expr_41), logFC = beta_rlm, p = p_rlm,
                   stringsAsFactors = FALSE)
res9 <- res9[!is.na(res9$p), ]
all_results$rlm <- res9
summary_tbl <- rbind(summary_tbl, eval_result(res9, "rlm"))

# ============================================================
# ============================================================
if (requireNamespace("ROTS", quietly = TRUE)) {
  log_info("--- 方法 8: ROTS ---")
  rots_res <- tryCatch({
    ROTS::ROTS(data = as.matrix(expr_41),
               groups = clin_41$group,
               B = 500, K = floor(ncol(expr_41) / 4),
               seed = SEED, verbose = FALSE)
  }, error = function(e) {
    log_warn("ROTS 失败: {conditionMessage(e)}"); NULL
  })
  if (!is.null(rots_res)) {
    res8 <- data.frame(gene = rownames(expr_41),
                        logFC = rots_res$logfc,
                        p = rots_res$pvalue,
                        stringsAsFactors = FALSE)
    all_results$ROTS <- res8
    summary_tbl <- rbind(summary_tbl, eval_result(res8, "ROTS"))
  }
} else {
  log_warn("ROTS 未安装，跳过方法 8")
}

# ============================================================
# 11. save
# ============================================================
log_info("--- 保存所有方法结果 ---")
for (nm in names(all_results)) {
  save_table(all_results[[nm]], DIR_THIS, paste0("method_", nm, ".csv"))
}
save_table(summary_tbl, DIR_THIS, "method_sweep_summary.csv")

cat("\n=== 各方法汇总 ===\n")
summary_tbl <- summary_tbl[order(-summary_tbl$n_fdr05), ]
print(summary_tbl[, c("method","n_fdr05","n_p05",
                       "median_effect","min_p")],
      row.names = FALSE)

# ============================================================
# ============================================================
log_info("--- 方法间一致性 ---")

sig_sets <- lapply(all_results, function(r) {
  r$FDR <- p.adjust(r$p, "BH")
  r$gene[!is.na(r$FDR) & r$FDR < 0.05]
})
sig_sets_nz <- sig_sets[sapply(sig_sets, length) >= 3]

if (length(sig_sets_nz) >= 2) {
  n <- length(sig_sets_nz)
  jmat <- matrix(NA, n, n,
                  dimnames = list(names(sig_sets_nz), names(sig_sets_nz)))
  for (i in 1:n) for (j in 1:n) {
    a <- sig_sets_nz[[i]]; b <- sig_sets_nz[[j]]
    u <- length(union(a, b))
    jmat[i, j] <- if (u == 0) 0 else length(intersect(a, b)) / u
  }
  save_table(data.frame(
    method1 = rownames(jmat)[row(jmat)],
    method2 = colnames(jmat)[col(jmat)],
    jaccard = round(as.vector(jmat), 3)),
    DIR_THIS, "method_overlap.csv")

  if (requireNamespace("pheatmap", quietly = TRUE)) {
    pheatmap::pheatmap(jmat, cluster_rows = TRUE, cluster_cols = TRUE,
                       display_numbers = TRUE, number_format = "%.2f",
                       main = "Method overlap (FDR<0.05)",
                       filename = file.path(DIR_THIS, "figures",
                                            "method_overlap.png"),
                       width = 8, height = 7)
  }
}

if (length(sig_sets_nz) >= 3) {
  gene_counts <- table(unlist(sig_sets_nz))
  consensus <- names(gene_counts[gene_counts >= 3])
  log_info("≥3 方法共识基因: {length(consensus)}")
  save_table(data.frame(gene = consensus,
                         n_methods = as.integer(gene_counts[consensus])),
              DIR_THIS, "consensus_genes.csv")

  cat("\n=== ≥3 方法共识基因 top 20 ===\n")
  cons_tbl <- data.frame(gene = consensus,
                          n_methods = as.integer(gene_counts[consensus]))
  cons_tbl <- cons_tbl[order(-cons_tbl$n_methods), ]
  print(head(cons_tbl, 20), row.names = FALSE)
}

# ============================================================
# 13. figure
# ============================================================
log_info("--- 图 ---")

if (nrow(summary_tbl) > 0) {
  summary_tbl$method <- factor(summary_tbl$method,
                                levels = summary_tbl$method[order(summary_tbl$n_fdr05)])
  p1 <- ggplot(summary_tbl, aes(x = method, y = n_fdr05)) +
    geom_col(fill = "#4C72B0") +
    geom_text(aes(label = n_fdr05), hjust = -0.2) +
    coord_flip() + theme_bw(base_size = 12) +
    labs(x = "", y = "FDR<0.05 基因数",
         title = "各方法显著基因数 (g4_vs_g1)")
  save_fig(p1, DIR_THIS, "methods_n_sig", width = 7, height = 5)
}

# ============================================================
# 14. summary
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("=== 02_4 汇总 ===\n")
cat(strrep("=", 70), "\n", sep = "")
cat("对比: g4_vs_g1  样本: ", ncol(expr_41),
    " (g1=", sum(clin_41$group=="g1"),
    ", g4=", sum(clin_41$group=="g4"), ")\n\n", sep = "")
cat("--- 各方法汇总 ---\n")
print(summary_tbl[, c("method","n_fdr05","n_fdr10","n_p05",
                       "median_effect","min_p")],
      row.names = FALSE)

cat("\n--- 输出位置 ---\n")
cat("目录:", DIR_THIS, "\n")

log_info("=== 02_4_single_gene_method_sweep.r done ===")


