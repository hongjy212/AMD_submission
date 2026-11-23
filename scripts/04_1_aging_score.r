#!/usr/bin/env Rscript
# ============================================================
# 04_1_aging_score.r (v2: 419)
#
#
#   - SenMayo (Saul 2022)
#   - hUSI (Human Universal Senescence Index)
#   - ML Aging
#   - CoreScence (DeepScence)
#
# Input：
#   results/01_1b_batch_enhance/data/expr_wgcna.rds        ★ 419
#   results/01_1b_batch_enhance/data/clinical_wgcna.rds    ★ 419
#   results/01_1_prepare_bulk_data/data/gene_annotation.rds
#   data/aging_scores.rds               —— 4 × 419 scorematrix
#   figures/aging_score_poly3F_bar.pdf/.png
#   figures/aging_score_grade_curves.pdf/.png
# ============================================================

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(limma); library(splines); library(ggplot2)
  library(dplyr); library(tidyr); library(data.table)
  library(readxl); library(stringr)
})

DIR_THIS <- setup_output_dir("04_1_aging_score")
log_info("=== 04_1_aging_score.r v2 (419) start ===")

# ============================================================
# ============================================================
AGING_DIR <- "E:/retina_tools/pubdata/aging"
SENMAYO_FILE    <- file.path(AGING_DIR, "SenMayo/41467_2022_32552_MOESM4_ESM.xlsx")
HUSI_SCRIPT     <- file.path(AGING_DIR, "hUSI/HUSI-main/sc/hUSI.R")
ML_FILE         <- file.path(AGING_DIR, "ML_Aging/mmc12.txt")
CORESCENCE_FILE <- file.path(AGING_DIR, "DeepScence/mmc2.xlsx")

# ============================================================
# ============================================================
log_info("--- 1. 载入数据 ---")

expr <- read_upstream("01_1b_batch_enhance", "data", "expr_wgcna.rds")
gene_anno <- read_upstream("01_1_prepare_bulk_data", "data", "gene_annotation.rds")
clinical <- read_upstream("01_1b_batch_enhance", "data", "clinical_wgcna.rds")

common <- intersect(colnames(expr), rownames(clinical))
expr <- expr[, common, drop = FALSE]
clinical <- clinical[common, , drop = FALSE]

# Ensembl → Symbol
anno_sub <- gene_anno[gene_anno$ensembl %in% rownames(expr), ]
anno_sub <- anno_sub[!is.na(anno_sub$symbol), ]
expr <- expr[anno_sub$ensembl, , drop = FALSE]
expr <- expr[!duplicated(anno_sub$symbol), , drop = FALSE]
rownames(expr) <- toupper(anno_sub$symbol[!duplicated(anno_sub$symbol)])

clinical$grade <- as.integer(clinical$grade)
clinical$age <- as.numeric(clinical$age_n)
rownames(clinical) <- colnames(expr)

log_info("expr: {nrow(expr)} × {ncol(expr)}")
log_info("clinical: {nrow(clinical)} 样本")
log_info("grade: {paste(names(table(clinical$grade)), table(clinical$grade), sep='=', collapse=' ')}")

# ============================================================
# ============================================================
log_info("--- 2. 加载衰老基因集 ---")

load_gene_set <- function(path, sheet = NULL, gene_col = NULL,
                          filter_col = NULL, filter_threshold = NULL) {
  if (!file.exists(path)) return(NULL)

  if (grepl("\\.xlsx$", path, ignore.case = TRUE)) {
    df <- read_excel(path, sheet = sheet)
  } else if (grepl("\\.txt$|\\.tsv$|\\.csv$", path, ignore.case = TRUE)) {
    df <- fread(path, header = TRUE, data.table = FALSE)
  } else {
    return(NULL)
  }

  if (is.null(gene_col)) {
    possible <- grep("gene|symbol|id|Name", colnames(df),
                       ignore.case = TRUE, value = TRUE)
    if (length(possible) == 0) return(NULL)
    gene_col <- possible[1]
  }
  if (!gene_col %in% colnames(df)) return(NULL)

  genes <- as.character(df[[gene_col]]) %>% na.omit() %>% str_trim()
  genes <- genes[genes != ""]

  if (!is.null(filter_col) && !is.null(filter_threshold)) {
    if (filter_col %in% colnames(df)) {
      vals <- as.numeric(df[[filter_col]])
      genes <- genes[!is.na(vals) & vals >= filter_threshold]
    }
  }
  return(unique(genes))
}

compute_score <- function(expr, genes, standardize = TRUE) {
  genes_present <- intersect(genes, rownames(expr))
  if (length(genes_present) < 5) return(rep(NA, ncol(expr)))
  scores <- colMeans(expr[genes_present, , drop = FALSE], na.rm = TRUE)
  if (standardize) scores <- scale(scores)[, 1]
  return(scores)
}

# ---- SenMayo ----
senmayo_genes <- NULL
if (file.exists(SENMAYO_FILE)) {
  senmayo_genes <- load_gene_set(SENMAYO_FILE, sheet = "human",
                                   gene_col = "Gene(human)")
  senmayo_genes <- toupper(senmayo_genes)
  log_info("SenMayo: {length(senmayo_genes)} 基因")
} else {
  log_warn("SenMayo 不存在")
}

# ---- hUSI ----
husi_available <- FALSE
if (file.exists(HUSI_SCRIPT)) {
  tryCatch({
    source(HUSI_SCRIPT, local = TRUE)
    husi_available <- TRUE
    log_info("hUSI 函数已加载")
  }, error = function(e) log_warn("hUSI source 失败"))
}

# ---- ML Aging ----
ml_up <- NULL; ml_down <- NULL
if (file.exists(ML_FILE)) {
  lines <- readLines(ML_FILE, warn = FALSE)
  lines <- trimws(lines[!grepl("^\\s*$", lines)])
  gene_idx <- which(lines == "Gene")[1]
  if (!is.na(gene_idx) && (length(lines) - gene_idx) >= 4) {
    data_lines <- lines[(gene_idx + 4):length(lines)]
    n_records <- length(data_lines) %/% 4
    records <- split(data_lines[1:(n_records * 4)], rep(1:n_records, each = 4))
    ml_list <- lapply(records, function(block) {
      if (length(block) < 4) return(NULL)
      weight_num <- as.numeric(str_extract(block[2],
                            "[-+]?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?"))
      if (is.na(weight_num)) return(NULL)
      data.frame(gene = block[1], weight = weight_num,
                 class = block[3], stringsAsFactors = FALSE)
    })
    ml_weights <- bind_rows(ml_list) %>%
      filter(!is.na(weight), class %in% c("young", "old"))
    ml_up <- toupper(ml_weights %>% filter(class == "old", weight > 0) %>%
                       pull(gene) %>% unique())
    ml_down <- toupper(ml_weights %>% filter(class == "young", weight > 0) %>%
                         pull(gene) %>% unique())
    log_info("ML Aging: up={length(ml_up)}, down={length(ml_down)}")
  }
}

# ---- CoreScence ----
corescence_genes <- NULL
if (file.exists(CORESCENCE_FILE)) {
  mmc2 <- read_excel(CORESCENCE_FILE)
  if ("Occurrence" %in% colnames(mmc2)) {
    occ_threshold <- max(2, median(as.numeric(mmc2$Occurrence), na.rm = TRUE))
    corescence_genes <- load_gene_set(CORESCENCE_FILE,
                                        gene_col = "Gene symbols",
                                        filter_col = "Occurrence",
                                        filter_threshold = occ_threshold)
  } else {
    corescence_genes <- load_gene_set(CORESCENCE_FILE, gene_col = "Gene symbols")
  }
  corescence_genes <- toupper(corescence_genes)
  log_info("CoreScence: {length(corescence_genes)} 基因")
}

# ============================================================
# ============================================================
log_info("--- 3. 计算衰老评分 ---")

score_list <- list()

if (!is.null(senmayo_genes) && length(senmayo_genes) >= 5) {
  score_list$SenMayo <- compute_score(expr, senmayo_genes, standardize = TRUE)
  log_info("SenMayo: 非 NA = {sum(!is.na(score_list$SenMayo))}")
}

if (husi_available) {
  score_list$hUSI <- tryCatch({
    scores <- cal_hUSI(as.matrix(expr))
    if (length(scores) == ncol(expr)) scale(scores)[, 1] else
      rep(NA, ncol(expr))
  }, error = function(e) rep(NA, ncol(expr)))
  log_info("hUSI: 非 NA = {sum(!is.na(score_list$hUSI))}")
}

if (!is.null(ml_up) && !is.null(ml_down) &&
    length(ml_up) >= 5 && length(ml_down) >= 5) {
  up_score <- compute_score(expr, ml_up, standardize = FALSE)
  down_score <- compute_score(expr, ml_down, standardize = FALSE)
  score_list$ML_Aging <- scale(up_score - down_score)[, 1]
  log_info("ML_Aging: 非 NA = {sum(!is.na(score_list$ML_Aging))}")
}

if (!is.null(corescence_genes) && length(corescence_genes) >= 5) {
  score_list$CoreScence <- compute_score(expr, corescence_genes, standardize = TRUE)
  log_info("CoreScence: 非 NA = {sum(!is.na(score_list$CoreScence))}")
}

if (length(score_list) == 0) stop("无可用衰老评分", call. = FALSE)

score_mat <- do.call(rbind, score_list)
colnames(score_mat) <- colnames(expr)
log_info("aging score 矩阵: {nrow(score_mat)} × {ncol(score_mat)}")
log_info("评分: {paste(rownames(score_mat), collapse=', ')}")

save_rds(score_mat, DIR_THIS, "aging_scores.rds")

# save CSV
score_df <- as.data.frame(t(score_mat))
score_df$sample_id <- rownames(score_df)
save_table(score_df, DIR_THIS, "aging_scores.csv")

# ============================================================
# 4. covariate
# ============================================================
sp <- splines::ns(clinical$age_n, df = 3)
covar_df <- data.frame(
  age_s1 = sp[,1], age_s2 = sp[,2], age_s3 = sp[,3],
  sex_M = as.integer(clinical$sex == "Male"),
  rin = clinical$rin_n, pmi = clinical$pmi_n
)
rownames(covar_df) <- rownames(clinical)

# ============================================================
# ============================================================
log_info("--- 4. limma multi (7 对比) ---")

design_df <- data.frame(grade = factor(clinical$grade), covar_df)
rownames(design_df) <- rownames(clinical)
design_df <- design_df[complete.cases(design_df), , drop = FALSE]

score_sub <- score_mat[, rownames(design_df), drop = FALSE]

design_multi <- model.matrix(
  ~ 0 + grade + age_s1 + age_s2 + age_s3 + sex_M + rin + pmi,
  data = design_df)
colnames(design_multi)[grep("^grade", colnames(design_multi))] <- paste0("g", 1:4)

cm <- limma::makeContrasts(
  g2_vs_g1      = g2 - g1,
  g3_vs_g2      = g3 - g2,
  g4_vs_g3      = g4 - g3,
  g4_vs_g1      = g4 - g1,
  anyAMD_vs_g1  = (g2 + g3 + g4)/3 - g1,
  late_vs_early = (g3 + g4)/2 - (g1 + g2)/2,
  g4_vs_rest    = g4 - (g1 + g2 + g3)/3,
  levels = design_multi)

fit_multi <- eBayes(contrasts.fit(lmFit(score_sub, design_multi), cm))

multi_res <- data.frame()
for (ct in colnames(cm)) {
  tt <- topTable(fit_multi, coef = ct, number = Inf, sort.by = "none")
  multi_res <- rbind(multi_res, data.frame(
    score = rownames(tt), contrast = ct,
    logFC = tt$logFC, t = tt$t, p = tt$P.Value, FDR = tt$adj.P.Val,
    stringsAsFactors = FALSE))
}
save_table(multi_res, DIR_THIS, "aging_score_limma_multi.csv")

cat("\n=== limma multi 各对比 FDR<0.05 ===\n")
for (ct in colnames(cm)) {
  n_sig <- sum(multi_res$FDR[multi_res$contrast == ct] < 0.05, na.rm = TRUE)
  cat(sprintf("  %-15s  %d / %d\n", ct, n_sig, nrow(score_mat)))
}

cat("\n=== limma multi 完整结果 ===\n")
print(multi_res[order(multi_res$p), ], row.names = FALSE, digits = 4)

# ============================================================
# ============================================================
log_info("--- 5. poly3_F ---")

poly_mat <- poly(clinical$grade, 3)
colnames(poly_mat) <- c("grade_p1", "grade_p2", "grade_p3")
rownames(poly_mat) <- rownames(clinical)
design_poly_df <- data.frame(poly_mat, covar_df)
rownames(design_poly_df) <- rownames(clinical)
design_poly_df <- design_poly_df[complete.cases(design_poly_df), , drop = FALSE]
score_poly <- score_mat[, rownames(design_poly_df), drop = FALSE]

design_poly <- model.matrix(
  ~ grade_p1 + grade_p2 + grade_p3 + age_s1 + age_s2 + age_s3 + sex_M + rin + pmi,
  data = design_poly_df)
fit_poly <- eBayes(lmFit(score_poly, design_poly))

coef_poly <- grep("^grade_p", colnames(design_poly), value = TRUE)
tt_F <- topTable(fit_poly, coef = coef_poly, number = Inf, sort.by = "none")
polyF <- data.frame(
  score = rownames(tt_F),
  F_stat = tt_F$F, p = tt_F$P.Value, FDR = tt_F$adj.P.Val,
  stringsAsFactors = FALSE)
polyF <- polyF[order(polyF$p), ]
save_table(polyF, DIR_THIS, "aging_score_poly3F.csv")

cat("\n=== poly3_F ===\n")
print(polyF, row.names = FALSE, digits = 4)

# ============================================================
# 7. agediagnosis（age=none / linear / spline）
# ============================================================
log_info("--- 6. 年龄诊断 ---")

run_aging_age <- function(score_mat, clin, age_mode) {
  covar <- data.frame(sex_M = as.integer(clin$sex == "Male"),
                       rin = clin$rin_n, pmi = clin$pmi_n)
  if (age_mode == "linear") {
    covar$age <- clin$age_n
  } else if (age_mode == "spline") {
    sp <- splines::ns(clin$age_n, df = 3)
    colnames(sp) <- c("age_s1", "age_s2", "age_s3")
    covar <- cbind(covar, sp)
  }
  rownames(covar) <- rownames(clin)

  design_df <- data.frame(grade = factor(clin$grade), covar)
  rownames(design_df) <- rownames(clin)
  design_df <- design_df[complete.cases(design_df), , drop = FALSE]
  score_sub <- score_mat[, rownames(design_df), drop = FALSE]

  fml <- as.formula(paste("~ 0 + grade +",
                          paste(colnames(covar), collapse = " + ")))
  design <- model.matrix(fml, data = design_df)
  colnames(design)[grep("^grade", colnames(design))] <- paste0("g", 1:4)

  cm <- limma::makeContrasts(
    g3_vs_g2 = g3 - g2,
    g4_vs_g1 = g4 - g1,
    levels = design)
  fit <- eBayes(contrasts.fit(lmFit(score_sub, design), cm))

  out <- data.frame()
  for (ct in colnames(cm)) {
    tt <- topTable(fit, coef = ct, number = Inf, sort.by = "none")
    out <- rbind(out, data.frame(
      score = rownames(tt), age_mode = age_mode, contrast = ct,
      logFC = tt$logFC, p = tt$P.Value, FDR = tt$adj.P.Val,
      stringsAsFactors = FALSE))
  }
  out
}

res_none   <- run_aging_age(score_mat, clinical, "none")
res_linear <- run_aging_age(score_mat, clinical, "linear")
res_spline <- run_aging_age(score_mat, clinical, "spline")

age_cmp <- rbind(res_none, res_linear, res_spline)
age_cmp <- age_cmp[order(age_cmp$score, age_cmp$contrast, age_cmp$age_mode), ]
save_table(age_cmp, DIR_THIS, "aging_score_age_comparison.csv")

cat("\n=== 年龄处理对 logFC 的影响 ===\n")
wide <- tidyr::pivot_wider(age_cmp,
                             id_cols = c("score", "contrast"),
                             names_from = age_mode,
                             values_from = c("logFC", "p"))
print(as.data.frame(wide), row.names = FALSE, digits = 4)

# ============================================================
# ============================================================
log_info("--- 7. 年龄双模型分类 ---")

age_class_tbl <- data.frame()
for (i in seq_len(nrow(score_mat))) {
  sc_name <- rownames(score_mat)[i]
  y <- score_mat[i, ]

  df <- data.frame(y = y, grade = clinical$grade, age = clinical$age_n,
                    sex = clinical$sex, rin = clinical$rin_n,
                    pmi = clinical$pmi_n)
  df <- df[complete.cases(df), ]

  sp_g <- splines::ns(df$age, df = 3)

  fit1 <- lm(y ~ poly(grade, 3) + sp_g + sex + rin + pmi, data = df)
  p_g <- anova(fit1)["poly(grade, 3)", "Pr(>F)"]

  fit2 <- lm(y ~ sp_g + poly(grade, 3) + sex + rin + pmi, data = df)
  p_a <- anova(fit2)["sp_g", "Pr(>F)"]

  r_age <- cor(df$y, df$age, method = "pearson")

  age_class_tbl <- rbind(age_class_tbl, data.frame(
    score = sc_name,
    p_grade = p_g,
    p_age = p_a,
    r_age = round(r_age, 3),
    stringsAsFactors = FALSE))
}

age_class_tbl$type <- "无关"
age_class_tbl$type[age_class_tbl$p_grade < 0.05 &
                   age_class_tbl$p_age > 0.10] <- "纯疾病"
age_class_tbl$type[age_class_tbl$p_grade < 0.05 &
                   age_class_tbl$p_age < 0.05] <- "年龄+疾病"
age_class_tbl$type[age_class_tbl$p_grade > 0.10 &
                   age_class_tbl$p_age < 0.05] <- "纯年龄"

save_table(age_class_tbl, DIR_THIS, "aging_score_age_classification.csv")

cat("\n=== 年龄双模型分类 ===\n")
print(age_class_tbl, row.names = FALSE, digits = 4)

# ============================================================
# 9. figure
# ============================================================
log_info("--- 8. 图 ---")

plot_df <- polyF
plot_df$score <- factor(plot_df$score,
                         levels = plot_df$score[order(plot_df$F_stat)])
plot_df$sig <- ifelse(plot_df$FDR < 0.05, "FDR<0.05",
                ifelse(plot_df$p < 0.05, "p<0.05", "ns"))

p1 <- ggplot(plot_df, aes(x = score, y = F_stat, fill = sig)) +
  geom_col() + coord_flip() + theme_bw(base_size = 12) +
  scale_fill_manual(values = c("FDR<0.05" = "#C44E52",
                                "p<0.05" = "#FDB462",
                                "ns" = "#B0B0B0")) +
  labs(x = "", y = "poly3_F F statistic", fill = "",
       title = "Aging score poly3_F (non-monotonic, 419)")
save_fig(p1, DIR_THIS, "aging_score_poly3F_bar", width = 6, height = 4)

curve_data <- data.frame()
for (sc_name in rownames(score_mat)) {
  for (g in 1:4) {
    idx <- clinical$grade == g
    curve_data <- rbind(curve_data, data.frame(
      score = sc_name,
      grade = g,
      mean_score = mean(score_mat[sc_name, idx], na.rm = TRUE),
      se = sd(score_mat[sc_name, idx], na.rm = TRUE) / sqrt(sum(idx)),
      stringsAsFactors = FALSE))
  }
}

p2 <- ggplot(curve_data, aes(x = grade, y = mean_score,
                               color = score, group = score)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean_score - se, ymax = mean_score + se),
                 width = 0.15, alpha = 0.5) +
  theme_bw(base_size = 12) +
  scale_x_continuous(breaks = 1:4, labels = paste0("G", 1:4)) +
  labs(x = "AMD grade", y = "Aging score (z, mean ± SE)",
       title = "Aging scores by grade (419)", color = "Score")
save_fig(p2, DIR_THIS, "aging_score_grade_curves", width = 7, height = 5)

# ============================================================
# 10. summary
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("=== 04_1 aging score 汇总 (v2, 419) ===\n")
cat(strrep("=", 70), "\n", sep = "")
cat("评分数量:", nrow(score_mat), "\n")
cat("样本数:", ncol(score_mat), "\n")
cat("评分:", paste(rownames(score_mat), collapse = ", "), "\n\n")

cat("--- limma multi 各对比 FDR<0.05 ---\n")
for (ct in colnames(cm)) {
  n_sig <- sum(multi_res$FDR[multi_res$contrast == ct] < 0.05, na.rm = TRUE)
  cat(sprintf("  %-15s  %d\n", ct, n_sig))
}

cat("\n--- poly3_F ---\n")
print(polyF, row.names = FALSE, digits = 4)

cat("\n--- 年龄双模型分类 ---\n")
print(age_class_tbl, row.names = FALSE, digits = 4)

cat("\n--- 输出位置 ---\n")
cat("目录:", DIR_THIS, "\n")

log_info("=== 04_1_aging_score.r v2 (419) done ===")

