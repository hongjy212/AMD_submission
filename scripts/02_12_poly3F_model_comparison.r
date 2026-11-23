# ============================================================
# 02_12_poly3F_model_comparison.r
#       - crossvalidation R²
# Output：results/02_12_poly3F_model_comparison/tables/
# ============================================================

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(splines)
})

ROOT <- DIR_ROOT
OUT  <- file.path(DIR_RESULT, "02_12_poly3F_model_comparison")
TAB  <- file.path(OUT, "tables")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
set.seed(SEED)
log_info("=== 02_12 poly3_F model comparison ===")

# ------------------------------------------------------------
# ------------------------------------------------------------
expr <- readRDS(file.path(ROOT, "results/01_1b_batch_enhance/data/expr_wgcna.rds"))
clin <- readRDS(file.path(ROOT, "results/01_1b_batch_enhance/data/clinical_wgcna.rds"))
anno <- readRDS(file.path(ROOT, "results/01_1_prepare_bulk_data/data/gene_annotation.rds"))

if ("ensembl" %in% colnames(anno)) {
  anno$ensembl_use <- as.character(anno$ensembl)
} else if ("ensembl_no_ver" %in% colnames(anno)) {
  anno$ensembl_use <- as.character(anno$ensembl_no_ver)
} else {
  anno$ensembl_use <- rownames(anno)
}
anno$ensembl_use <- sub("\\..*$", "", anno$ensembl_use)
anno$symbol <- as.character(anno$symbol)
keep <- anno$ensembl_use %in% rownames(expr) & !is.na(anno$symbol) & anno$symbol != ""
keep[is.na(keep)] <- FALSE
map_df <- unique(anno[keep, c("symbol", "ensembl_use")])
colnames(map_df) <- c("symbol", "ensembl")

common <- intersect(colnames(expr), clin$sample_id)
expr <- expr[, common, drop = FALSE]
clin <- clin[match(common, clin$sample_id), ]

base_df <- data.frame(
  grade = suppressWarnings(as.numeric(as.character(clin$grade))),
  age   = suppressWarnings(as.numeric(clin$age)),
  sex   = as.factor(clin$sex),
  rin   = suppressWarnings(as.numeric(clin$rin)),
  pmi   = suppressWarnings(as.numeric(clin$pmi)),
  stringsAsFactors = FALSE
)
complete_idx <- complete.cases(base_df)
base_df <- base_df[complete_idx, , drop = FALSE]
expr    <- expr[, complete_idx, drop = FALSE]
log_info(sprintf("complete samples: %d", nrow(base_df)))

# 16 genes
sig16_sym <- read_csv(file.path(DIR_RESULT,
  "02_1_differential_expression_main/tables/sig16_genes.csv"),
  show_col_types = FALSE)$gene
sig16_ens <- map_df$ensembl[match(sig16_sym, map_df$symbol)]
sig16_ens <- sig16_ens[!is.na(sig16_ens) & sig16_ens %in% rownames(expr)]
log_info(sprintf("16-gene mapped: %d", length(sig16_ens)))

# ------------------------------------------------------------
# ------------------------------------------------------------
fit_models <- function(y, dat) {
  f_null <- lm(y ~ ns(age, df = 3) + sex + rin + pmi, data = dat)
  f_lm   <- lm(y ~ grade + ns(age, df = 3) + sex + rin + pmi, data = dat)
  f_q    <- lm(y ~ poly(grade, 2) + ns(age, df = 3) + sex + rin + pmi, data = dat)
  f_c    <- lm(y ~ poly(grade, 3) + ns(age, df = 3) + sex + rin + pmi, data = dat)

  # AIC/BIC
  aic <- sapply(list(null = f_null, linear = f_lm, quadratic = f_q, cubic = f_c), AIC)
  bic <- sapply(list(null = f_null, linear = f_lm, quadratic = f_q, cubic = f_c), BIC)

  lrt_lin <- anova(f_null, f_lm)$`Pr(>F)`[2]
  lrt_q   <- anova(f_lm, f_q)$`Pr(>F)`[2]
  lrt_c   <- anova(f_q, f_c)$`Pr(>F)`[2]
  lrt_c_joint <- anova(f_null, f_c)$`Pr(>F)`[2]

  set.seed(SEED)
  folds <- sample(rep(1:5, length.out = nrow(dat)))
  cv_r2 <- function(formula) {
    preds <- numeric(nrow(dat))
    for (k in 1:5) {
      tr <- dat[folds != k, , drop = FALSE]
      te <- dat[folds == k, , drop = FALSE]
      tr$y <- y[folds != k]
      te$y <- y[folds == k]
      fit <- tryCatch(lm(formula, data = tr), error = function(e) NULL)
      if (is.null(fit)) { preds[folds == k] <- NA; next }
      preds[folds == k] <- tryCatch(predict(fit, newdata = te), error = function(e) NA)
    }
    ss_res <- sum((y - preds)^2, na.rm = TRUE)
    ss_tot <- sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)
    max(0, 1 - ss_res / ss_tot)
  }
  r2_null <- cv_r2(y ~ ns(age, df = 3) + sex + rin + pmi)
  r2_lm   <- cv_r2(y ~ grade + ns(age, df = 3) + sex + rin + pmi)
  r2_q    <- cv_r2(y ~ poly(grade, 2) + ns(age, df = 3) + sex + rin + pmi)
  r2_c    <- cv_r2(y ~ poly(grade, 3) + ns(age, df = 3) + sex + rin + pmi)

  data.frame(
    AIC_null = aic["null"], AIC_linear = aic["linear"],
    AIC_quadratic = aic["quadratic"], AIC_cubic = aic["cubic"],
    BIC_null = bic["null"], BIC_linear = bic["linear"],
    BIC_quadratic = bic["quadratic"], BIC_cubic = bic["cubic"],
    dAIC_cubic_vs_quad = aic["quadratic"] - aic["cubic"],
    dBIC_cubic_vs_quad = bic["quadratic"] - bic["cubic"],
    p_lrt_linear_vs_null = lrt_lin,
    p_lrt_quad_vs_linear = lrt_q,
    p_lrt_cubic_vs_quad = lrt_c,
    p_joint_cubic_vs_null = lrt_c_joint,
    cvR2_null = r2_null, cvR2_linear = r2_lm,
    cvR2_quad = r2_q, cvR2_cubic = r2_c,
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------
# ------------------------------------------------------------
log_info("--- fitting 4 nested models per gene ---")
res_list <- list()
for (ens in sig16_ens) {
  sym <- map_df$symbol[match(ens, map_df$ensembl)]
  y <- as.numeric(expr[ens, ])
  dat <- base_df
  res <- tryCatch(fit_models(y, dat), error = function(e) NULL)
  if (is.null(res)) next
  res$gene <- sym
  res$ensembl <- ens
  res_list[[ens]] <- res
}
model_cmp <- bind_rows(res_list) %>%
  select(gene, ensembl, everything())
print(model_cmp, row.names = FALSE)
write_csv(model_cmp, file.path(TAB, "poly3F_model_comparison.csv"))

# ------------------------------------------------------------
# ------------------------------------------------------------
log_info("--- model selection summary ---")
summary_df <- model_cmp %>%
  summarise(
    n_genes = n(),
    best_AIC_null = sum(AIC_null == pmin(AIC_null, AIC_linear, AIC_quadratic, AIC_cubic)),
    best_AIC_linear = sum(AIC_linear == pmin(AIC_null, AIC_linear, AIC_quadratic, AIC_cubic)),
    best_AIC_quad = sum(AIC_quadratic == pmin(AIC_null, AIC_linear, AIC_quadratic, AIC_cubic)),
    best_AIC_cubic = sum(AIC_cubic == pmin(AIC_null, AIC_linear, AIC_quadratic, AIC_cubic)),
    best_BIC_null = sum(BIC_null == pmin(BIC_null, BIC_linear, BIC_quadratic, BIC_cubic)),
    best_BIC_linear = sum(BIC_linear == pmin(BIC_null, BIC_linear, BIC_quadratic, BIC_cubic)),
    best_BIC_quad = sum(BIC_quadratic == pmin(BIC_null, BIC_linear, BIC_quadratic, BIC_cubic)),
    best_BIC_cubic = sum(BIC_cubic == pmin(BIC_null, BIC_linear, BIC_quadratic, BIC_cubic)),
    median_dAIC_cubic_vs_quad = median(dAIC_cubic_vs_quad, na.rm = TRUE),
    median_dBIC_cubic_vs_quad = median(dBIC_cubic_vs_quad, na.rm = TRUE),
    n_p_lrt_cubic_vs_quad_lt_005 = sum(p_lrt_cubic_vs_quad < 0.05, na.rm = TRUE),
    n_p_joint_cubic_vs_null_lt_005 = sum(p_joint_cubic_vs_null < 0.05, na.rm = TRUE),
    median_cvR2_null = median(cvR2_null, na.rm = TRUE),
    median_cvR2_linear = median(cvR2_linear, na.rm = TRUE),
    median_cvR2_quad = median(cvR2_quad, na.rm = TRUE),
    median_cvR2_cubic = median(cvR2_cubic, na.rm = TRUE)
  )
print(summary_df, row.names = FALSE)
write_csv(summary_df, file.path(TAB, "poly3F_model_selection_summary.csv"))

# ------------------------------------------------------------
# ------------------------------------------------------------
if (requireNamespace("qvalue", quietly = TRUE)) {
  library(qvalue)
  all_ens <- rownames(expr)
  all_p <- numeric(length(all_ens))
  for (i in seq_along(all_ens)) {
    y <- as.numeric(expr[all_ens[i], ])
    dat <- base_df; dat$y <- y
    fit_c <- tryCatch(lm(y ~ poly(grade, 3) + ns(age, df = 3) + sex + rin + pmi,
                         data = dat), error = function(e) NULL)
    fit_n <- tryCatch(lm(y ~ ns(age, df = 3) + sex + rin + pmi,
                         data = dat), error = function(e) NULL)
    if (is.null(fit_c) || is.null(fit_n)) { all_p[i] <- NA; next }
    all_p[i] <- tryCatch(anova(fit_n, fit_c)$`Pr(>F)`[2], error = function(e) NA)
  }
  valid <- !is.na(all_p)
  q <- qvalue(all_p[valid])
  q_df <- data.frame(
    ensembl = all_ens[valid],
    symbol = map_df$symbol[match(all_ens[valid], map_df$ensembl)],
    p_poly3F = all_p[valid],
    q_value = q$qvalues,
    pi0 = q$pi0
  )
  q_df_sig16 <- q_df %>% filter(symbol %in% sig16_sym)
  cat("\n=== 16-gene q-value ===\n")
  print(q_df_sig16, row.names = FALSE)
  write_csv(q_df, file.path(TAB, "poly3F_all_genes_qvalue.csv"))
  write_csv(q_df_sig16, file.path(TAB, "poly3F_sig16_qvalue.csv"))
} else {
  log_warn("qvalue package not installed; skipping local FDR")
}

# ------------------------------------------------------------
# ------------------------------------------------------------
log_info("--- batch correction sensitivity ---")
rin_models <- read_or_null <- function(path) {
  if (file.exists(path)) read_csv(path, show_col_types = FALSE) else NULL
}
rin_file <- file.path(DIR_RESULT, "diag/diag_rin_sensitivity/tables/rin_robust_summary.csv")
if (file.exists(rin_file)) {
  rin_df <- read_csv(rin_file, show_col_types = FALSE)
  cat("\n=== RIN sensitivity ===\n")
  print(rin_df, row.names = FALSE)
}

log_info("=== 02_12 done ===")
