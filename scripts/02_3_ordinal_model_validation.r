# ============================================================
# 02_3_ordinal_model_validation.r
#   A. MASS::polr (Wald + LRT)
# Input：
#   results/01_1_prepare_bulk_data/data/{expr_corrected,clinical,gene_annotation}.rds
#   tables/polr_results.csv
# ============================================================

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(limma); library(MASS); library(splines)
})

DIR_THIS <- setup_output_dir("02_3_ordinal_model_validation")
log_info("=== 02_3_ordinal_model_validation.r start ===")

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

log_info("表达矩阵: {nrow(expr)} × {ncol(expr)}")

# ============================================================
# 2. covariate
# ============================================================
sp <- splines::ns(clin$age_n, df = 3)
covar_df <- data.frame(
  age_s1 = sp[,1], age_s2 = sp[,2], age_s3 = sp[,3],
  sex_M = as.integer(clin$sex == "Male"),
  rin = clin$rin_n, pmi = clin$pmi_n
)
rownames(covar_df) <- rownames(clin)

grade_ord <- factor(clin$grade, ordered = TRUE)

# ============================================================
# ============================================================
log_info("--- 3. 候选基因 ---")

design_df <- data.frame(grade = factor(clin$grade), covar_df)
rownames(design_df) <- rownames(clin)
expr_sub <- expr[, rownames(design_df), drop = FALSE]

design_multi <- model.matrix(
  ~ 0 + grade + age_s1 + age_s2 + age_s3 + sex_M + rin + pmi,
  data = design_df)
colnames(design_multi)[grep("^grade", colnames(design_multi))] <- paste0("g", 1:4)

cm <- limma::makeContrasts(g3_vs_g2 = g3 - g2, g4_vs_g1 = g4 - g1,
                            levels = design_multi)
fit <- eBayes(contrasts.fit(lmFit(expr_sub, design_multi), cm))
tt_g3v2 <- topTable(fit, coef = "g3_vs_g2", number = Inf, sort.by = "none")
tt_g4v1 <- topTable(fit, coef = "g4_vs_g1", number = Inf, sort.by = "none")

cand_genes <- unique(c(
  rownames(tt_g3v2)[tt_g3v2$P.Value < 0.01],
  rownames(tt_g4v1)[tt_g4v1$P.Value < 0.01]
))
log_info("候选基因 (limma p<0.01): {length(cand_genes)}")

# ============================================================
# 4. polr Wald + LRT
# ============================================================
log_info("--- 4. polr (Wald + LRT) ---")

fit_polr_one <- function(y, grade_ord, covar_df) {
  df <- data.frame(grade = grade_ord, y = as.numeric(y),
                   covar_df, stringsAsFactors = FALSE)
  df <- df[complete.cases(df), , drop = FALSE]
  if (nrow(df) < 30) return(NULL)
  tryCatch({
    m0 <- MASS::polr(grade ~ age_s1 + age_s2 + age_s3 + sex_M + rin + pmi,
                     data = df, Hess = TRUE)
    m1 <- MASS::polr(grade ~ y + age_s1 + age_s2 + age_s3 + sex_M + rin + pmi,
                     data = df, Hess = TRUE)
    a <- anova(m0, m1, test = "Chisq")
    p_lrt <- a[["Pr(Chi)"]][2]
    cf <- coef(summary(m1))
    if (!"y" %in% rownames(cf)) return(NULL)
    z <- cf["y", "t value"]
    list(beta = cf["y", "Value"],
         se = cf["y", "Std. Error"],
         z = z,
         p_wald = 2 * pnorm(-abs(z)),
         p_lrt = p_lrt,
         n = nrow(df))
  }, error = function(e) NULL)
}

log_info("polr 拟合 {length(cand_genes)} 个基因...")
t0 <- Sys.time()
polr_res <- vector("list", length(cand_genes))
for (i in seq_along(cand_genes)) {
  polr_res[[i]] <- fit_polr_one(expr[cand_genes[i], ], grade_ord, covar_df)
  if (i %% 500 == 0) {
    dt <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
    log_info("  polr {i}/{length(cand_genes)} ({dt}s)")
  }
}
log_info("polr 总用时: {round(as.numeric(difftime(Sys.time(), t0, units='secs')), 1)}s")

res_polr <- data.frame(
  gene = cand_genes,
  beta = sapply(polr_res, function(r) if (is.null(r)) NA_real_ else r$beta),
  se = sapply(polr_res, function(r) if (is.null(r)) NA_real_ else r$se),
  z = sapply(polr_res, function(r) if (is.null(r)) NA_real_ else r$z),
  p_wald = sapply(polr_res, function(r) if (is.null(r)) NA_real_ else r$p_wald),
  p_lrt = sapply(polr_res, function(r) if (is.null(r)) NA_real_ else r$p_lrt),
  n = sapply(polr_res, function(r) if (is.null(r)) NA_integer_ else r$n),
  stringsAsFactors = FALSE
)
res_polr$FDR_wald <- p.adjust(res_polr$p_wald, "BH")
res_polr$FDR_lrt <- p.adjust(res_polr$p_lrt, "BH")
save_table(res_polr, DIR_THIS, "polr_results.csv")

log_info("polr Wald FDR<0.05: {sum(res_polr$FDR_wald < 0.05, na.rm=TRUE)}")
log_info("polr LRT FDR<0.05: {sum(res_polr$FDR_lrt < 0.05, na.rm=TRUE)}")
log_info("polr LRT p<0.001: {sum(res_polr$p_lrt < 0.001, na.rm=TRUE)}")

cat("\n=== polr top 20 (LRT) ===\n")
print(head(res_polr[order(res_polr$p_lrt),
                    c("gene","beta","p_wald","p_lrt","FDR_lrt")], 20),
      row.names = FALSE)

# ============================================================
# ============================================================
if (requireNamespace("ordinal", quietly = TRUE)) {
  log_info("--- 5. clm (ordinal 包) ---")
  library(ordinal)

  fit_clm_lrt <- function(y, grade_ord, covar_df) {
    df <- data.frame(grade = grade_ord, y = as.numeric(y),
                     covar_df, stringsAsFactors = FALSE)
    df <- df[complete.cases(df), , drop = FALSE]
    if (nrow(df) < 30) return(NULL)
    tryCatch({
      m0 <- clm(grade ~ age_s1 + age_s2 + age_s3 + sex_M + rin + pmi, data = df)
      m1 <- clm(grade ~ y + age_s1 + age_s2 + age_s3 + sex_M + rin + pmi, data = df)
      a <- anova(m0, m1, test = "LRT")
      p_lrt <- a[["Pr(>Chi)"]][2]
      cf <- coef(summary(m1))
      list(beta = cf["y", "Estimate"],
           z = cf["y", "z value"],
           p_lrt = p_lrt)
    }, error = function(e) NULL)
  }

  clm_res <- vector("list", length(cand_genes))
  for (i in seq_along(cand_genes)) {
    clm_res[[i]] <- fit_clm_lrt(expr[cand_genes[i], ], grade_ord, covar_df)
  }
  res_clm <- data.frame(
    gene = cand_genes,
    beta = sapply(clm_res, function(r) if (is.null(r)) NA_real_ else r$beta),
    z = sapply(clm_res, function(r) if (is.null(r)) NA_real_ else r$z),
    p_lrt = sapply(clm_res, function(r) if (is.null(r)) NA_real_ else r$p_lrt),
    stringsAsFactors = FALSE)
  res_clm$FDR_lrt <- p.adjust(res_clm$p_lrt, "BH")
  save_table(res_clm, DIR_THIS, "clm_LRT_results.csv")
  log_info("clm LRT FDR<0.05: {sum(res_clm$FDR_lrt < 0.05, na.rm=TRUE)}")
  log_info("clm LRT p<0.001: {sum(res_clm$p_lrt < 0.001, na.rm=TRUE)}")
} else {
  log_warn("ordinal 未安装，跳过 clm")
}

# ============================================================
# 6. summary
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("=== 02_3 汇总 ===\n")
cat(strrep("=", 70), "\n", sep = "")
cat("候选基因:", length(cand_genes), "\n")
cat("polr 有效结果:", sum(!is.na(res_polr$p_lrt)), "/", nrow(res_polr), "\n")
cat("polr Wald FDR<0.05:", sum(res_polr$FDR_wald < 0.05, na.rm=TRUE), "\n")
cat("polr LRT  FDR<0.05:", sum(res_polr$FDR_lrt  < 0.05, na.rm=TRUE), "\n")

cat("\n--- 输出位置 ---\n")
cat("目录:", DIR_THIS, "\n")
cat("  tables/polr_results.csv\n")
if (requireNamespace("ordinal", quietly = TRUE)) {
  cat("  tables/clm_LRT_results.csv\n")
}

log_info("=== 02_3_ordinal_model_validation.r done ===")

