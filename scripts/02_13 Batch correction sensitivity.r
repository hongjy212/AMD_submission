# ============================================================
# 02_13_batch_correction_sensitivity.r
#       - C：ComBat with grade as covariate
# Output：results/02_13_batch_correction_sensitivity/tables/
# ============================================================

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr); library(splines)
  library(limma); library(sva)
})

ROOT <- DIR_ROOT
OUT  <- file.path(DIR_RESULT, "02_13_batch_correction_sensitivity")
TAB  <- file.path(OUT, "tables")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
set.seed(SEED)
log_info("=== 02_13 batch correction sensitivity ===")

# ------------------------------------------------------------
# ------------------------------------------------------------
voom_full  <- readRDS(file.path(DIR_RESULT,
  "01_1_prepare_bulk_data/data/voom_full.rds"))
clin_wgcna <- readRDS(file.path(DIR_RESULT,
  "01_1b_batch_enhance/data/clinical_wgcna.rds"))
expr_wgcna <- readRDS(file.path(DIR_RESULT,
  "01_1b_batch_enhance/data/expr_wgcna.rds"))
anno       <- readRDS(file.path(DIR_RESULT,
  "01_1_prepare_bulk_data/data/gene_annotation.rds"))

common <- intersect(colnames(voom_full), clin_wgcna$sample_id)
common <- intersect(common, colnames(expr_wgcna))
voom_full  <- voom_full[, common, drop = FALSE]
expr_wgcna <- expr_wgcna[, common, drop = FALSE]
clin_wgcna <- clin_wgcna[match(common, clin_wgcna$sample_id), ]
log_info(sprintf("samples: %d", length(common)))

# ------------------------------------------------------------
# ------------------------------------------------------------
batch_cols <- intersect(c("batch", "lib_prep_batch", "rna_isolation_batch"),
                        colnames(clin_wgcna))
log_info(sprintf("batch columns: %s", paste(batch_cols, collapse = ", ")))

base_df <- data.frame(
  sample_id = clin_wgcna$sample_id,
  grade = suppressWarnings(as.numeric(as.character(clin_wgcna$grade))),
  age   = suppressWarnings(as.numeric(clin_wgcna$age)),
  sex   = as.factor(clin_wgcna$sex),
  rin   = suppressWarnings(as.numeric(clin_wgcna$rin)),
  pmi   = suppressWarnings(as.numeric(clin_wgcna$pmi)),
  stringsAsFactors = FALSE
)
for (bc in batch_cols) {
  base_df[[bc]] <- as.factor(clin_wgcna[[bc]])
}
complete_idx <- complete.cases(base_df)
base_df <- base_df[complete_idx, , drop = FALSE]
voom_full  <- voom_full[, complete_idx, drop = FALSE]
expr_wgcna <- expr_wgcna[, complete_idx, drop = FALSE]
log_info(sprintf("complete samples: %d", nrow(base_df)))

# ------------------------------------------------------------
# 3. genesannotationmapping
# ------------------------------------------------------------
if ("ensembl" %in% colnames(anno)) {
  anno$ensembl_use <- as.character(anno$ensembl)
} else if ("ensembl_no_ver" %in% colnames(anno)) {
  anno$ensembl_use <- as.character(anno$ensembl_no_ver)
} else {
  anno$ensembl_use <- rownames(anno)
}
anno$ensembl_use <- sub("\\..*$", "", anno$ensembl_use)
anno$symbol      <- as.character(anno$symbol)
keep <- anno$ensembl_use %in% rownames(voom_full) &
        !is.na(anno$symbol) & anno$symbol != ""
keep[is.na(keep)] <- FALSE
map_df <- unique(anno[keep, c("symbol", "ensembl_use")])
colnames(map_df) <- c("symbol", "ensembl")
log_info(sprintf("gene annotation mapped: %d", nrow(map_df)))

# ------------------------------------------------------------
# ------------------------------------------------------------

complexIII_sym <- c("UQCRQ", "UQCRFS1", "UQCRH")
complexIII_ens <- map_df$ensembl[match(complexIII_sym, map_df$symbol)]
complexIII_ens <- complexIII_ens[!is.na(complexIII_ens) &
                                  complexIII_ens %in% rownames(voom_full)]
log_info(sprintf("Complex III: %d genes", length(complexIII_ens)))

oxphos_file_rds <- file.path(DIR_RESULT,
  "01_2_prepare_gene_sets/data/gene_sets.rds")
oxphos_sym <- NULL
oxphos_ens <- NULL

if (file.exists(oxphos_file_rds)) {
  gs <- readRDS(oxphos_file_rds)
  log_info(sprintf("gene_sets.rds: %d gene sets", length(gs)))

  candidate_names <- c("OXPHOS", "KEGG_OXPHOS", "GO_OXPHOS",
                       "Mito_OXPHOS", "Reactome_OXPHOS")
  hit_name <- intersect(candidate_names, names(gs))[1]

  if (!is.na(hit_name)) {
    entry <- gs[[hit_name]]
    if (is.list(entry) && "genes" %in% names(entry)) {
      oxphos_sym <- toupper(entry$genes)
    } else if (is.character(entry)) {
      oxphos_sym <- toupper(entry)
    }
    log_info(sprintf("OXPHOS gene set '%s': %d genes",
                     hit_name, length(oxphos_sym)))

    oxphos_ens <- map_df$ensembl[match(oxphos_sym, map_df$symbol)]
    oxphos_ens <- oxphos_ens[!is.na(oxphos_ens) &
                             oxphos_ens %in% rownames(voom_full)]
    log_info(sprintf("OXPHOS mapped to Ensembl: %d genes", length(oxphos_ens)))
  } else {
    log_warn(sprintf("no OXPHOS in gene_sets.rds; available names (first 20): %s",
                     paste(head(names(gs), 20), collapse = ", ")))
  }
} else {
  log_warn(sprintf("gene_sets.rds not found at %s", oxphos_file_rds))
}

if (is.null(oxphos_ens) || length(oxphos_ens) < 20) {
  log_warn("using prefix-based OXPHOS fallback")
  oxphos_idx <- grep("^(NDUF|UQCR|CYC1|COX[0-9]|COX4I|COX5|COX6|COX7|COX8|ATP5|ATP6V)",
                     map_df$symbol)
  oxphos_ens <- map_df$ensembl[oxphos_idx]
  oxphos_ens <- oxphos_ens[!is.na(oxphos_ens) &
                           oxphos_ens %in% rownames(voom_full)]
  log_info(sprintf("prefix-based OXPHOS: %d genes", length(oxphos_ens)))
}

# --- 16-gene signature ---
sig16_sym <- read_csv(file.path(DIR_RESULT,
  "02_1_differential_expression_main/tables/sig16_genes.csv"),
  show_col_types = FALSE)$gene
sig16_ens <- map_df$ensembl[match(sig16_sym, map_df$symbol)]
sig16_ens <- sig16_ens[!is.na(sig16_ens) & sig16_ens %in% rownames(voom_full)]
log_info(sprintf("Complex III: %d | OXPHOS: %d | 16-gene: %d",
                 length(complexIII_ens),
                 length(oxphos_ens),
                 length(sig16_ens)))

# ------------------------------------------------------------
# 5. three batch correction
# ------------------------------------------------------------
correct_limma <- function(expr_mat, clin_df, batch_cols, protect_grade = TRUE) {
  if (protect_grade) {
    rhs <- c("grade", "age", "sex", "rin", "pmi", batch_cols)
  } else {
    rhs <- c("age", "sex", "rin", "pmi", batch_cols)
  }
  formula_str <- paste("~", paste(rhs, collapse = " + "))
  design <- model.matrix(as.formula(formula_str), data = clin_df)

  batch_coefs <- grep(paste0("^(", paste(batch_cols, collapse = "|"), ")"),
                      colnames(design), value = TRUE)
  log_info(sprintf("protect_grade=%s, design cols=%d, batch coefs=%d",
                   protect_grade, ncol(design), length(batch_coefs)))

  fit <- lmFit(expr_mat, design)
  corrected <- expr_mat
  for (bc in batch_coefs) {
    corrected <- corrected - fit$coefficients[, bc, drop = FALSE] %*%
      t(design[, bc, drop = FALSE])
  }
  corrected
}

correct_combat <- function(expr_mat, clin_df, batch_cols) {
  batch <- clin_df[[batch_cols[1]]]
  mod <- model.matrix(~ grade + age + sex + rin + pmi, data = clin_df)
  ComBat(expr_mat, batch = batch, mod = mod, par.prior = TRUE,
         prior.plots = FALSE)
}

log_info("--- Method A: limma, protect grade ---")
expr_A <- correct_limma(voom_full, base_df, batch_cols, protect_grade = TRUE)

log_info("--- Method B: limma, no grade protection ---")
expr_B <- correct_limma(voom_full, base_df, batch_cols, protect_grade = FALSE)

log_info("--- Method C: ComBat with grade covariate ---")
expr_C <- tryCatch(
  correct_combat(voom_full, base_df, batch_cols),
  error = function(e) {
    log_warn(sprintf("ComBat failed: %s", conditionMessage(e)))
    NULL
  })

# ------------------------------------------------------------
# ------------------------------------------------------------
compute_signal <- function(expr_mat, label) {
  if (is.null(expr_mat)) return(NULL)
  expr_mat <- as.matrix(expr_mat)
  z_mat <- t(scale(t(expr_mat)))

  if (length(complexIII_ens) >= 2) {
    cIII_score <- colMeans(z_mat[complexIII_ens, , drop = FALSE], na.rm = TRUE)
    dat <- base_df; dat$score <- cIII_score
    fit_cIII <- lm(score ~ grade + age + sex + rin + pmi, data = dat)
    beta_cIII <- coef(fit_cIII)["grade"]
    p_cIII <- summary(fit_cIII)$coefficients["grade", "Pr(>|t|)"]
  } else {
    beta_cIII <- NA; p_cIII <- NA
  }

  if (length(oxphos_ens) >= 3) {
    ox_score <- colMeans(z_mat[oxphos_ens, , drop = FALSE], na.rm = TRUE)
    dat <- base_df; dat$score <- ox_score
    fit_ox <- lm(score ~ grade + age + sex + rin + pmi, data = dat)
    beta_ox <- coef(fit_ox)["grade"]
    p_ox <- summary(fit_ox)$coefficients["grade", "Pr(>|t|)"]
  } else {
    beta_ox <- NA; p_ox <- NA
  }

  poly3_p <- sapply(sig16_ens, function(ens) {
    y <- as.numeric(expr_mat[ens, ])
    dat <- base_df; dat$y <- y
    fit_c <- tryCatch(lm(y ~ poly(grade, 3) + ns(age, df = 3) + sex + rin + pmi,
                         data = dat), error = function(e) NULL)
    fit_n <- tryCatch(lm(y ~ ns(age, df = 3) + sex + rin + pmi,
                         data = dat), error = function(e) NULL)
    if (is.null(fit_c) || is.null(fit_n)) return(NA_real_)
    tryCatch(anova(fit_n, fit_c)$`Pr(>F)`[2], error = function(e) NA_real_)
  })
  poly3_fdr <- p.adjust(poly3_p, method = "BH")
  n_poly3_sig <- sum(poly3_fdr < 0.05, na.rm = TRUE)

  data.frame(
    method            = label,
    beta_ComplexIII   = round(beta_cIII, 4),
    p_ComplexIII      = signif(p_cIII, 3),
    beta_OXPHOS       = round(beta_ox, 4),
    p_OXPHOS          = signif(p_ox, 3),
    n_poly3_sig       = n_poly3_sig,
    n_poly3_valid     = sum(!is.na(poly3_p)),
    median_FDR_16gene = signif(median(poly3_fdr, na.rm = TRUE), 3),
    stringsAsFactors = FALSE
  )
}

res <- bind_rows(
  compute_signal(expr_A,     "A_limma_protect_grade"),
  compute_signal(expr_B,     "B_limma_no_protection"),
  compute_signal(expr_C,     "C_ComBat_grade_covariate"),
  compute_signal(expr_wgcna, "D_current_published")
)
cat("\n=== Batch correction sensitivity summary ===\n")
print(res, row.names = FALSE)
write_csv(res, file.path(TAB, "batch_correction_sensitivity.csv"))

# ------------------------------------------------------------
# ------------------------------------------------------------
fit_16gene_fdr <- function(expr_mat) {
  if (is.null(expr_mat)) return(rep(NA_real_, length(sig16_ens)))
  sapply(sig16_ens, function(ens) {
    y <- as.numeric(expr_mat[ens, ])
    dat <- base_df; dat$y <- y
    fit_c <- tryCatch(lm(y ~ poly(grade, 3) + ns(age, df = 3) + sex + rin + pmi,
                         data = dat), error = function(e) NULL)
    fit_n <- tryCatch(lm(y ~ ns(age, df = 3) + sex + rin + pmi,
                         data = dat), error = function(e) NULL)
    if (is.null(fit_c) || is.null(fit_n)) return(NA_real_)
    tryCatch(anova(fit_n, fit_c)$`Pr(>F)`[2], error = function(e) NA_real_)
  })
}

cmp16 <- data.frame(
  gene  = map_df$symbol[match(sig16_ens, map_df$ensembl)],
  FDR_A = p.adjust(fit_16gene_fdr(expr_A), "BH"),
  FDR_B = p.adjust(fit_16gene_fdr(expr_B), "BH"),
  FDR_C = p.adjust(fit_16gene_fdr(expr_C), "BH"),
  FDR_D = p.adjust(fit_16gene_fdr(expr_wgcna), "BH"),
  stringsAsFactors = FALSE
)
cmp16$sig_A <- cmp16$FDR_A < 0.05
cmp16$sig_B <- cmp16$FDR_B < 0.05
cmp16$sig_C <- cmp16$FDR_C < 0.05
cmp16$sig_D <- cmp16$FDR_D < 0.05
cmp16$consistent_AB <- cmp16$sig_A == cmp16$sig_B

cat("\n=== 16-gene poly3_F FDR across methods ===\n")
print(cmp16, row.names = FALSE)
write_csv(cmp16, file.path(TAB, "batch_correction_16gene_fdr.csv"))

# ------------------------------------------------------------
# 8. summary
# ------------------------------------------------------------
cat("\n", strrep("=", 60), "\n", sep = "")
cat("=== 02_13 汇总 ===\n")
cat(strrep("=", 60), "\n", sep = "")

cat("\n[1] 核心 OXPHOS 信号 across 4 种方法：\n")
print(res, row.names = FALSE)

cat("\n[2] 16-gene poly3_F: A vs B 显著性一致率：\n")
cat(sprintf("   一致：%d / %d\n",
            sum(cmp16$consistent_AB, na.rm = TRUE), nrow(cmp16)))
cat(sprintf("   n sig in A (protect):    %d\n", sum(cmp16$sig_A)))
cat(sprintf("   n sig in B (no protect): %d\n", sum(cmp16$sig_B)))
cat(sprintf("   n sig in C (ComBat):     %d\n", sum(cmp16$sig_C)))
cat(sprintf("   n sig in D (published):  %d\n", sum(cmp16$sig_D)))
cat(sprintf("   Median FDR_A = %.4f, FDR_B = %.4f, FDR_C = %.4f, FDR_D = %.4f\n",
            median(cmp16$FDR_A, na.rm = TRUE),
            median(cmp16$FDR_B, na.rm = TRUE),
            median(cmp16$FDR_C, na.rm = TRUE),
            median(cmp16$FDR_D, na.rm = TRUE)))

log_info("=== 02_13 done ===")

