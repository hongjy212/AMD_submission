#!/usr/bin/env Rscript
# ============================================================
# diag_tissue.r
#
#
#
# Input：
#   results/01_1_prepare_bulk_data/data/{clinical,gene_annotation}.rds
#   results/01_2_prepare_gene_sets/data/gene_sets.rds
#   results/04_2_pathway_score/data/pathway_scores_ssgsea.rds
#
#   tables/column_probe.csv
#   tables/cramers_v.csv
#   tables/contingency_*.csv
#   tables/tissue_covariate_cross.csv
#   tables/tissue_pathway_kw.csv
#   tables/tissue_adjusted_effects.csv
#   tables/confounding_classification.csv
#   figures/tissue_grade_bar.pdf/.png
#   figures/cramers_v_bar.pdf/.png
# ============================================================

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(limma); library(splines); library(ggplot2)
})

DIR_THIS <- setup_output_dir("diag_tissue", subdir = "diag")
log_info("=== diag_tissue.r start ===")

# ============================================================
# ============================================================
log_info("--- 1. 载入 ---")

clin <- read_upstream("01_1_prepare_bulk_data", "data", "clinical.rds")
GENE_SETS <- read_upstream("01_2_prepare_gene_sets", "data", "gene_sets.rds")
score_mat <- read_upstream("04_2_pathway_score", "data",
                            "pathway_scores_ssgsea.rds")

clin$grade <- as.integer(clin$grade)
clin$age_v <- as.numeric(clin$age_n)
clin$pmi_v <- as.numeric(clin$pmi_n)
clin$rin_v <- as.numeric(clin$rin_n)

log_info("n samples: {nrow(clin)}")
log_info("cols: {paste(colnames(clin), collapse=', ')}")

# ============================================================
# ============================================================
log_info("--- 2. 探测候选分类列 ---")

CANDIDATE_COLS <- c(
  "tissue", "region", "macula", "peripheral", "sample_type",
  "study", "dataset", "source", "GSE", "geo", "batch", "cohort",
  "mgs_level", "platform", "site"
)
CANDIDATE_COLS <- intersect(CANDIDATE_COLS, colnames(clin))

probe <- data.frame()
for (cc in CANDIDATE_COLS) {
  v <- clin[[cc]]
  n_uniq <- length(unique(v[!is.na(v)]))
  top_val <- names(sort(table(v), decreasing = TRUE))[1]
  top_prop <- round(max(table(v)) / sum(!is.na(v)), 3)
  probe <- rbind(probe, data.frame(
    column = cc,
    n_unique = n_uniq,
    top_value = top_val,
    top_proportion = top_prop,
    n_missing = sum(is.na(v)),
    stringsAsFactors = FALSE))
}
save_table(probe, DIR_THIS, "column_probe.csv")

cat("\n=== 候选列探测 ===\n")
print(probe, row.names = FALSE)

usable_cols <- probe$column[probe$n_unique >= 2 & probe$n_unique <= 50]

single_cols <- probe$column[probe$n_unique == 1]
if (length(single_cols) > 0) {
  cat("\n=== 单值列（无变异，跳过）===\n")
  for (cc in single_cols) {
    cat(sprintf("  %s: 全部 = '%s'\n", cc, clin[[cc]][1]))
  }
}

grade_dup <- sapply(usable_cols, function(cc) {
  tab <- table(clin$grade, clin[[cc]])
  all(rowSums(tab > 0) == 1) && all(colSums(tab > 0) == 1)
})
if (any(grade_dup)) {
  cat("\n=== 与 grade 完全共线的列（跳过）===\n")
  for (cc in usable_cols[grade_dup]) cat("  ", cc, "\n")
  usable_cols <- usable_cols[!grade_dup]
}

log_info("可用分类列: {paste(usable_cols, collapse=', ')}")
if (length(usable_cols) == 0) {
  cat("\n无可用分类列，diag_tissue 结束。\n")
  quit(save = "no")
}

# ============================================================
# 3. Cramér's V
# ============================================================
log_info("--- 3. Cramér's V vs grade ---")

cramers_v <- function(tab) {
  chi <- suppressWarnings(chisq.test(tab, correct = FALSE))
  n <- sum(tab)
  k <- min(nrow(tab), ncol(tab))
  if (k <= 1) return(NA_real_)
  v <- sqrt(as.numeric(chi$statistic) / (n * (k - 1)))
  p <- chi$p.value
  list(V = round(v, 3), p = signif(p, 3))
}

cv_tbl <- data.frame()
for (cc in usable_cols) {
  tab <- table(clin$grade, clin[[cc]])
  cv <- cramers_v(tab)
  n_empty <- sum(tab == 0)
  n_cells <- length(tab)
  empty_prop <- round(n_empty / n_cells, 3)
  max_prop_per_grade <- max(apply(tab, 1, function(x) max(x) / sum(x)))
  cv_tbl <- rbind(cv_tbl, data.frame(
    column = cc,
    n_categories = ncol(tab),
    cramers_V = cv$V,
    chi_p = cv$p,
    empty_cells = n_empty,
    empty_proportion = empty_prop,
    max_prop_per_grade = round(max_prop_per_grade, 3),
    stringsAsFactors = FALSE))
}
cv_tbl <- cv_tbl[order(-cv_tbl$cramers_V), ]
save_table(cv_tbl, DIR_THIS, "cramers_v.csv")

cat("\n=== Cramér's V (vs grade) ===\n")
print(cv_tbl, row.names = FALSE)

classify <- function(v, empty_prop) {
  if (is.na(v)) return("NA")
  if (v >= 0.6 || empty_prop > 0.5) return("完全混杂")
  if (v >= 0.4) return("强混杂")
  if (v >= 0.2) return("中度混杂")
  return("弱混杂")
}
cv_tbl$confounding <- sapply(seq_len(nrow(cv_tbl)), function(i) {
  classify(cv_tbl$cramers_V[i], cv_tbl$empty_proportion[i])
})
save_table(cv_tbl, DIR_THIS, "confounding_classification.csv")

cat("\n=== 混杂等级 ===\n")
print(table(cv_tbl$confounding))

# ============================================================
# ============================================================
log_info("--- 4. 列联表 ---")

for (cc in usable_cols) {
  tab <- table(clin$grade, clin[[cc]])
  tab_df <- as.data.frame.matrix(tab)
  tab_df$grade <- rownames(tab_df)
  tab_df <- tab_df[, c("grade", setdiff(colnames(tab_df), "grade"))]
  save_table(tab_df, DIR_THIS, paste0("contingency_", cc, ".csv"))
}

# ============================================================
# ============================================================
log_info("--- 5. 分类列 vs 协变量 ---")

cross_out <- data.frame()
for (cc in usable_cols) {
  v <- clin[[cc]]
  # vs age
  r_age <- suppressWarnings(cor(as.integer(factor(v)),
                                 clin$age_v, method = "spearman",
                                 use = "complete.obs"))
  # vs PMI
  r_pmi <- suppressWarnings(cor(as.integer(factor(v)),
                                 clin$pmi_v, method = "spearman",
                                 use = "complete.obs"))
  # vs RIN
  r_rin <- suppressWarnings(cor(as.integer(factor(v)),
                                 clin$rin_v, method = "spearman",
                                 use = "complete.obs"))
  tab_sex <- table(v, clin$sex)
  p_sex <- suppressWarnings(chisq.test(tab_sex)$p.value)
  cross_out <- rbind(cross_out, data.frame(
    column = cc,
    rho_age = round(r_age, 3),
    rho_pmi = round(r_pmi, 3),
    rho_rin = round(r_rin, 3),
    sex_chi_p = signif(p_sex, 3),
    stringsAsFactors = FALSE))
}
save_table(cross_out, DIR_THIS, "tissue_covariate_cross.csv")

cat("\n=== 分类列 vs 协变量 ===\n")
print(cross_out, row.names = FALSE)

# ============================================================
# ============================================================
log_info("--- 6. 通路评分 vs 分类列 ---")

common_samp <- intersect(colnames(score_mat), rownames(clin))
kw_out <- data.frame()
for (cc in usable_cols) {
  v <- clin[common_samp, cc]
  if (length(unique(v[!is.na(v)])) < 2) next
  for (gs in rownames(score_mat)) {
    y <- score_mat[gs, common_samp]
    df <- data.frame(y = y, g = v)
    df <- df[complete.cases(df), ]
    if (length(unique(df$g)) < 2) next
    kw <- suppressWarnings(kruskal.test(y ~ g, data = df))
    # eta²: (H - k + 1) / (n - k)
    k <- length(unique(df$g))
    n <- nrow(df)
    eta2 <- max(0, (as.numeric(kw$statistic) - k + 1) / (n - k))
    kw_out <- rbind(kw_out, data.frame(
      column = cc, pathway = gs,
      kw_p = signif(kw$p.value, 3),
      eta2 = round(eta2, 4),
      stringsAsFactors = FALSE))
  }
}
save_table(kw_out, DIR_THIS, "tissue_pathway_kw.csv")

cat("\n=== 每个分类列 top 10 通路 (by eta²) ===\n")
for (cc in usable_cols) {
  sub <- kw_out[kw_out$column == cc, ]
  if (nrow(sub) == 0) next
  sub <- sub[order(-sub$eta2), ]
  cat(sprintf("\n-- %s --\n", cc))
  print(head(sub[, c("pathway","kw_p","eta2")], 10),
        row.names = FALSE)
}

# ============================================================
# ============================================================
log_info("--- 7. 控分类列前后效应变化 ---")

score_sub <- score_mat[, common_samp, drop = FALSE]

run_model <- function(score_mat, clin, model_name, use_cols = character(0),
                       use_age = TRUE) {
  covar <- data.frame(
    sex_M = as.integer(clin$sex == "Male"),
    rin = clin$rin_v, pmi = clin$pmi_v,
    stringsAsFactors = FALSE)
  if (use_age) {
    sp <- splines::ns(clin$age_v, df = 3)
    colnames(sp) <- c("age_s1", "age_s2", "age_s3")
    covar <- cbind(covar, sp)
  }
  for (cc in use_cols) {
    covar[[cc]] <- factor(clin[[cc]])
  }
  rownames(covar) <- rownames(clin)
  design_df <- data.frame(grade = factor(clin$grade), covar)
  rownames(design_df) <- rownames(clin)
  design_df <- design_df[complete.cases(design_df), , drop = FALSE]
  score_use <- score_mat[, rownames(design_df), drop = FALSE]

  fml <- as.formula(paste("~ 0 + grade +",
                          paste(setdiff(colnames(covar), ""),
                                collapse = " + ")))
  design <- model.matrix(fml, data = design_df)
  colnames(design)[grep("^grade", colnames(design))] <- paste0("g", 1:4)

  cm <- limma::makeContrasts(
    g3_vs_g2 = g3 - g2,
    g4_vs_g3 = g4 - g3,
    levels = design)
  fit <- eBayes(contrasts.fit(lmFit(score_use, design), cm))

  out <- data.frame()
  for (ct in colnames(cm)) {
    tt <- topTable(fit, coef = ct, number = Inf, sort.by = "none")
    out <- rbind(out, data.frame(
      pathway = rownames(tt), model = model_name, contrast = ct,
      logFC = tt$logFC, p = tt$P.Value, FDR = tt$adj.P.Val,
      stringsAsFactors = FALSE))
  }
  out
}

models <- list()
models[[length(models) + 1]] <- list(name = "M1_age_PMI_RIN",
                                       cols = character(0))
for (cc in usable_cols) {
  models[[length(models) + 1]] <- list(
    name = paste0("M1_plus_", cc), cols = cc)
}

all_models <- data.frame()
for (m in models) {
  r <- run_model(score_sub, clin, m$name, use_cols = m$cols)
  all_models <- rbind(all_models, r)
}

gs_used_path <- file.path(DIR_RESULT, "04_2_pathway_score", "tables",
                          "gene_sets_used.csv")
if (file.exists(gs_used_path)) {
  gs_used <- read.csv(gs_used_path, stringsAsFactors = FALSE)
  all_models$group <- gs_used$group[match(all_models$pathway, gs_used$geneset)]
}
save_table(all_models, DIR_THIS, "tissue_adjusted_effects.csv")

cat("\n=== 各模型下 g3_vs_g2 FDR<0.05 数量 ===\n")
for (m in unique(all_models$model)) {
  n <- sum(all_models$FDR[all_models$model == m &
                          all_models$contrast == "g3_vs_g2"] < 0.05,
           na.rm = TRUE)
  cat(sprintf("  %-30s  %d\n", m, n))
}

cat("\n=== 控某列后消失的 g3_vs_g2 通路 ===\n")
m1 <- all_models[all_models$model == "M1_age_PMI_RIN" &
                  all_models$contrast == "g3_vs_g2",
                  c("pathway","logFC","FDR")]
for (cc in usable_cols) {
  mn <- paste0("M1_plus_", cc)
  m2 <- all_models[all_models$model == mn &
                    all_models$contrast == "g3_vs_g2",
                    c("pathway","logFC","FDR")]
  cmp <- merge(m1, m2, by = "pathway", suffixes = c("_M1","_Mx"))
  driven <- cmp$FDR_M1 < 0.05 & cmp$FDR_Mx >= 0.05
  robust <- cmp$FDR_M1 < 0.05 & cmp$FDR_Mx < 0.05
  cat(sprintf("  %s: M1 显著 %d，控后消失 %d，控后稳健 %d\n",
              cc, sum(cmp$FDR_M1 < 0.05), sum(driven), sum(robust)))
  if (sum(driven) > 0) {
    cat("    消失的通路:", paste(cmp$pathway[driven], collapse = ", "), "\n")
  }
}

# ============================================================
# 8. figure
# ============================================================
log_info("--- 8. 图 ---")

for (cc in usable_cols) {
  tab_df <- as.data.frame(table(grade = clin$grade, cat = clin[[cc]]))
  p <- ggplot(tab_df, aes(x = factor(grade), y = Freq, fill = cat)) +
    geom_col(position = "stack") + theme_bw(base_size = 11) +
    labs(x = "AMD grade", y = "Count", fill = cc,
         title = sprintf("%s vs grade", cc))
  save_fig(p, DIR_THIS, paste0("tissue_grade_bar_", cc),
           width = 7, height = 5)
}

if (nrow(cv_tbl) > 0) {
  cv_tbl$column_f <- factor(cv_tbl$column,
                             levels = cv_tbl$column[order(cv_tbl$cramers_V)])
  p <- ggplot(cv_tbl, aes(x = column_f, y = cramers_V,
                            fill = confounding)) +
    geom_col(alpha = 0.8) +
    geom_hline(yintercept = c(0.2, 0.4, 0.6), linetype = "dashed",
                color = "gray40") +
    coord_flip() + theme_bw(base_size = 11) +
    labs(x = "", y = "Cramér's V (vs grade)", fill = "混杂等级",
         title = "Confounding of categorical columns with grade",
         subtitle = "0.2 weak, 0.4 strong, 0.6 complete")
  save_fig(p, DIR_THIS, "cramers_v_bar", width = 7, height = 5)
}

# ============================================================
# 9. summary
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("=== diag_tissue.r 汇总 ===\n")
cat(strrep("=", 70), "\n", sep = "")

cat("\n可用分类列:", paste(usable_cols, collapse = ", "), "\n")
if (length(single_cols) > 0) {
  cat("单值列（跳过）:", paste(single_cols, collapse = ", "), "\n")
}
if (exists("grade_dup") && any(grade_dup)) {
  cat("共线列（跳过）:",
      paste(names(grade_dup)[grade_dup], collapse = ", "), "\n")
}

cat("\n混杂等级:\n")
print(table(cv_tbl$confounding))

cat("\n最需要关注的列:\n")
print(head(cv_tbl[, c("column","cramers_V","chi_p","confounding")], 5),
      row.names = FALSE)

cat("\n--- 输出位置 ---\n")
cat("目录:", DIR_THIS, "\n")
cat("  tables/column_probe.csv\n")
cat("  tables/cramers_v.csv\n")
cat("  tables/confounding_classification.csv\n")
cat("  tables/tissue_covariate_cross.csv\n")
cat("  tables/tissue_pathway_kw.csv\n")
cat("  tables/tissue_adjusted_effects.csv\n")
cat("  figures/tissue_grade_bar_*.pdf/.png\n")
cat("  figures/cramers_v_bar.pdf/.png\n")

log_info("=== diag_tissue.r done ===")

