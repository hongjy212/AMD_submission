#!/usr/bin/env Rscript
# ============================================================
# diag_batch_sensitivity.r
#
#   C. 02_1 16 genes poly3_F
#
#
# Input：
#   results/01_1_prepare_bulk_data/data/{expr_corrected,clinical,gene_annotation}.rds
#   results/01_2_prepare_gene_sets/data/gene_sets.rds
#   results/02_1_differential_expression_main/data/...（16 genes）
#   results/04_1_aging_score/data/aging_scores.rds
#   results/04_2_pathway_score/data/pathway_scores_ssgsea.rds
#
#   figures/batch_effect_summary.pdf/.png
# ============================================================

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(limma); library(splines); library(ggplot2); library(tidyr)
})

DIR_THIS <- setup_output_dir("diag_batch_sensitivity", subdir = "diag")
log_info("=== diag_batch_sensitivity.r start ===")

# ============================================================
# ============================================================
log_info("--- 1. 载入 ---")

expr <- read_upstream("01_1_prepare_bulk_data", "data", "expr_corrected.rds")
clin <- read_upstream("01_1_prepare_bulk_data", "data", "clinical.rds")
anno <- read_upstream("01_1_prepare_bulk_data", "data", "gene_annotation.rds")
GENE_SETS <- read_upstream("01_2_prepare_gene_sets", "data", "gene_sets.rds")
score_mat <- read_upstream("04_2_pathway_score", "data",
                            "pathway_scores_ssgsea.rds")

common <- intersect(colnames(expr), rownames(clin))
expr <- expr[, common, drop = FALSE]
clin <- clin[common, , drop = FALSE]

anno_sub <- anno[anno$ensembl %in% rownames(expr), ]
anno_sub <- anno_sub[!is.na(anno_sub$symbol), ]
expr <- expr[anno_sub$ensembl, , drop = FALSE]
expr <- expr[!duplicated(anno_sub$symbol), , drop = FALSE]
rownames(expr) <- toupper(anno_sub$symbol[!duplicated(anno_sub$symbol)])

clin$grade <- as.integer(clin$grade)
clin$age_v <- as.numeric(clin$age_n)
clin$pmi_v <- as.numeric(clin$pmi_n)
clin$rin_v <- as.numeric(clin$rin_n)
clin$batch_f <- factor(clin$batch)
rownames(clin) <- colnames(expr)

log_info("expr: {nrow(expr)} × {ncol(expr)}")
log_info("score_mat: {nrow(score_mat)} × {ncol(score_mat)}")

# ============================================================
# ============================================================
run_limma <- function(score_mat, clin, use_batch = FALSE,
                       contrasts = c("g3_vs_g2", "g4_vs_g3", "g4_vs_g1")) {
  covar <- data.frame(
    age_s1 = splines::ns(clin$age_v, df = 3)[, 1],
    age_s2 = splines::ns(clin$age_v, df = 3)[, 2],
    age_s3 = splines::ns(clin$age_v, df = 3)[, 3],
    sex_M = as.integer(clin$sex == "Male"),
    rin = clin$rin_v, pmi = clin$pmi_v,
    stringsAsFactors = FALSE)
  if (use_batch) {
    covar$batch_f <- clin$batch_f
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
    g2_vs_g1 = g2 - g1,
    g3_vs_g2 = g3 - g2,
    g4_vs_g3 = g4 - g3,
    g4_vs_g1 = g4 - g1,
    levels = design)
  fit <- eBayes(contrasts.fit(lmFit(score_sub, design), cm))

  out <- data.frame()
  for (ct in contrasts) {
    if (!ct %in% colnames(cm)) next
    tt <- topTable(fit, coef = ct, number = Inf, sort.by = "none")
    out <- rbind(out, data.frame(
      feature = rownames(tt), contrast = ct,
      logFC = tt$logFC, p = tt$P.Value, FDR = tt$adj.P.Val,
      stringsAsFactors = FALSE))
  }
  out
}

run_poly3F <- function(score_mat, clin, use_batch = FALSE) {
  poly_mat <- poly(clin$grade, 3)
  colnames(poly_mat) <- c("grade_p1", "grade_p2", "grade_p3")
  rownames(poly_mat) <- rownames(clin)

  covar <- data.frame(
    poly_mat,
    age_s1 = splines::ns(clin$age_v, df = 3)[, 1],
    age_s2 = splines::ns(clin$age_v, df = 3)[, 2],
    age_s3 = splines::ns(clin$age_v, df = 3)[, 3],
    sex_M = as.integer(clin$sex == "Male"),
    rin = clin$rin_v, pmi = clin$pmi_v,
    stringsAsFactors = FALSE)
  if (use_batch) {
    covar$batch_f <- clin$batch_f
  }
  rownames(covar) <- rownames(clin)
  covar <- covar[complete.cases(covar), , drop = FALSE]
  score_sub <- score_mat[, rownames(covar), drop = FALSE]

  fml <- as.formula(paste("~ grade_p1 + grade_p2 + grade_p3 +",
                          paste(setdiff(colnames(covar),
                                        c("grade_p1","grade_p2","grade_p3")),
                                collapse = " + ")))
  design <- model.matrix(fml, data = covar)
  fit <- eBayes(lmFit(score_sub, design))
  coef_poly <- grep("^grade_p", colnames(design), value = TRUE)
  tt <- topTable(fit, coef = coef_poly, number = Inf, sort.by = "none")
  data.frame(feature = rownames(tt),
             F_stat = tt$F, p = tt$P.Value, FDR = tt$adj.P.Val,
             stringsAsFactors = FALSE)
}

# ============================================================
# ============================================================
log_info("--- 3. 模块 A: 通路 ---")

pathway_nobatch <- run_limma(score_mat, clin, use_batch = FALSE)
pathway_batch   <- run_limma(score_mat, clin, use_batch = TRUE)

pathway_cmp <- merge(
  pathway_nobatch[pathway_nobatch$contrast == "g3_vs_g2",
                   c("feature","logFC","FDR")],
  pathway_batch[pathway_batch$contrast == "g3_vs_g2",
                 c("feature","logFC","FDR")],
  by = "feature", suffixes = c("_nobatch","_batch"))
pathway_cmp$batch_driven <- pathway_cmp$FDR_nobatch < 0.05 &
                             pathway_cmp$FDR_batch >= 0.05
pathway_cmp$robust <- pathway_cmp$FDR_nobatch < 0.05 &
                       pathway_cmp$FDR_batch < 0.05
pathway_cmp <- pathway_cmp[order(pathway_cmp$FDR_nobatch), ]
save_table(pathway_cmp, DIR_THIS, "pathway_batch_effects.csv")

pathway_summary <- data.frame(
  module = "04_2_pathway",
  nobatch_sig = sum(pathway_cmp$FDR_nobatch < 0.05, na.rm = TRUE),
  batch_sig = sum(pathway_cmp$FDR_batch < 0.05, na.rm = TRUE),
  batch_driven = sum(pathway_cmp$batch_driven, na.rm = TRUE),
  robust = sum(pathway_cmp$robust, na.rm = TRUE))

cat("\n=== 04_2 通路 ===\n")
print(pathway_summary, row.names = FALSE)
cat("\n消失的通路:\n")
print(pathway_cmp$feature[pathway_cmp$batch_driven], row.names = FALSE)

# ============================================================
# ============================================================
log_info("--- 4. 模块 B: hUSI ---")

husi_file <- file.path(DIR_RESULT, "04_1_aging_score", "data",
                       "aging_scores.rds")
if (file.exists(husi_file)) {
  aging_scores <- readRDS(husi_file)
  common_samp <- intersect(colnames(aging_scores), rownames(clin))
  aging_sub <- aging_scores[, common_samp, drop = FALSE]

  husi_nobatch <- run_limma(aging_sub, clin[common_samp, ],
                             use_batch = FALSE,
                             contrasts = c("g3_vs_g2","g4_vs_g3"))
  husi_batch   <- run_limma(aging_sub, clin[common_samp, ],
                             use_batch = TRUE,
                             contrasts = c("g3_vs_g2","g4_vs_g3"))

  husi_cmp <- merge(husi_nobatch, husi_batch,
                     by = c("feature","contrast"),
                     suffixes = c("_nobatch","_batch"))
  husi_cmp <- husi_cmp[order(husi_cmp$contrast, husi_cmp$p_nobatch), ]
  save_table(husi_cmp, DIR_THIS, "hUSI_batch_effects.csv")

  cat("\n=== 04_1 hUSI ===\n")
  print(husi_cmp[, c("feature","contrast",
                      "logFC_nobatch","FDR_nobatch",
                      "logFC_batch","FDR_batch")],
        row.names = FALSE, digits = 4)

  husi_g3 <- husi_cmp[husi_cmp$contrast == "g3_vs_g2", ]
  husi_summary <- data.frame(
    module = "04_1_hUSI_g3v2",
    nobatch_sig = sum(husi_g3$FDR_nobatch < 0.05, na.rm = TRUE),
    batch_sig = sum(husi_g3$FDR_batch < 0.05, na.rm = TRUE),
    batch_driven = sum(husi_g3$FDR_nobatch < 0.05 &
                        husi_g3$FDR_batch >= 0.05, na.rm = TRUE),
    robust = sum(husi_g3$FDR_nobatch < 0.05 &
                  husi_g3$FDR_batch < 0.05, na.rm = TRUE))
} else {
  log_warn("04_1 aging_scores.rds 不存在，跳过")
  husi_summary <- NULL
}

# ============================================================
# 5. module C: 02_1 16 genes poly3_F
# ============================================================
log_info("--- 5. 模块 C: 16 基因 ---")

GENES_16 <- c("FCF1","UGGT1","CCDC3","S100B","BPNT1","BRI3BP","BTBD19",
              "DPH7","AKR1B1","ZNF623","ZNF587","RSPH3","TAF3","NDUFB8",
              "ZNF417","C6orf120")
genes_in <- intersect(GENES_16, rownames(expr))
log_info("16 基因命中: {length(genes_in)}")

if (length(genes_in) >= 5) {
  expr_16 <- expr[genes_in, , drop = FALSE]
  common_samp16 <- intersect(colnames(expr_16), rownames(clin))
  expr_16 <- expr_16[, common_samp16, drop = FALSE]
  clin_16 <- clin[common_samp16, ]

  # poly3_F
  poly_nobatch <- run_poly3F(expr_16, clin_16, use_batch = FALSE)
  poly_batch   <- run_poly3F(expr_16, clin_16, use_batch = TRUE)
  poly_cmp <- merge(poly_nobatch, poly_batch,
                     by = "feature",
                     suffixes = c("_nobatch","_batch"))
  poly_cmp <- poly_cmp[order(poly_cmp$p_nobatch), ]
  save_table(poly_cmp, DIR_THIS, "poly3F_16genes_batch.csv")

  cat("\n=== 02_1 16 基因 poly3_F ===\n")
  print(poly_cmp[, c("feature","F_stat_nobatch","FDR_nobatch",
                      "F_stat_batch","FDR_batch")],
        row.names = FALSE, digits = 4)

  limma_nobatch <- run_limma(expr_16, clin_16, use_batch = FALSE,
                              contrasts = c("g3_vs_g2","g4_vs_g3"))
  limma_batch   <- run_limma(expr_16, clin_16, use_batch = TRUE,
                              contrasts = c("g3_vs_g2","g4_vs_g3"))
  limma_cmp <- merge(limma_nobatch, limma_batch,
                      by = c("feature","contrast"),
                      suffixes = c("_nobatch","_batch"))
  save_table(limma_cmp, DIR_THIS, "limma_16genes_batch.csv")

  g3 <- limma_cmp[limma_cmp$contrast == "g3_vs_g2", ]
  poly_summary <- data.frame(
    module = "02_1_16genes",
    nobatch_sig = sum(poly_cmp$FDR_nobatch < 0.05, na.rm = TRUE),
    batch_sig = sum(poly_cmp$FDR_batch < 0.05, na.rm = TRUE),
    batch_driven = sum(poly_cmp$FDR_nobatch < 0.05 &
                        poly_cmp$FDR_batch >= 0.05, na.rm = TRUE),
    robust = sum(poly_cmp$FDR_nobatch < 0.05 &
                  poly_cmp$FDR_batch < 0.05, na.rm = TRUE))
} else {
  log_warn("16 基因命中不足")
  poly_summary <- NULL
}

# ============================================================
# ============================================================
log_info("--- 6. 模块 D: celltype ---")

celltype_file <- file.path(DIR_RESULT, "03_1_celltype_marker_score", "data",
                            "celltype_scores.rds")
if (file.exists(celltype_file)) {
  celltype_scores <- readRDS(celltype_file)
  common_samp_ct <- intersect(colnames(celltype_scores), rownames(clin))
  celltype_sub <- celltype_scores[, common_samp_ct, drop = FALSE]

  ct_nobatch <- run_limma(celltype_sub, clin[common_samp_ct, ],
                           use_batch = FALSE,
                           contrasts = c("g3_vs_g2","g4_vs_g3"))
  ct_batch   <- run_limma(celltype_sub, clin[common_samp_ct, ],
                           use_batch = TRUE,
                           contrasts = c("g3_vs_g2","g4_vs_g3"))
  ct_cmp <- merge(ct_nobatch, ct_batch,
                   by = c("feature","contrast"),
                   suffixes = c("_nobatch","_batch"))
  save_table(ct_cmp, DIR_THIS, "celltype_batch_effects.csv")

  ct_g3 <- ct_cmp[ct_cmp$contrast == "g3_vs_g2", ]
  ct_summary <- data.frame(
    module = "03_1_celltype_g3v2",
    nobatch_sig = sum(ct_g3$FDR_nobatch < 0.05, na.rm = TRUE),
    batch_sig = sum(ct_g3$FDR_batch < 0.05, na.rm = TRUE),
    batch_driven = sum(ct_g3$FDR_nobatch < 0.05 &
                        ct_g3$FDR_batch >= 0.05, na.rm = TRUE),
    robust = sum(ct_g3$FDR_nobatch < 0.05 &
                  ct_g3$FDR_batch < 0.05, na.rm = TRUE))
} else {
  log_warn("03_1 celltype_scores.rds 不存在，跳过")
  ct_summary <- NULL
}

# ============================================================
# ============================================================
log_info("--- 7. 总表 ---")

batch_summary <- do.call(rbind, Filter(Negate(is.null),
  list(pathway_summary, husi_summary, poly_summary, ct_summary)))
save_table(batch_summary, DIR_THIS, "batch_summary.csv")

cat("\n=== batch 敏感性总表 ===\n")
print(batch_summary, row.names = FALSE)

robust_conclusions <- data.frame()
if (!is.null(husi_summary) && husi_summary$robust > 0) {
  robust_conclusions <- rbind(robust_conclusions, data.frame(
    module = "hUSI",
    conclusion = "hUSI 在 g3_vs_g2 显著且控 batch 后稳健",
    evidence = sprintf("%d/%d 稳健", husi_summary$robust,
                        husi_summary$nobatch_sig)))
}
if (!is.null(poly_summary) && poly_summary$robust > 0) {
  robust_conclusions <- rbind(robust_conclusions, data.frame(
    module = "16_genes",
    conclusion = "16 基因 poly3_F 控 batch 后稳健",
    evidence = sprintf("%d/%d 稳健", poly_summary$robust,
                        poly_summary$nobatch_sig)))
}
if (pathway_summary$batch_driven > 0) {
  robust_conclusions <- rbind(robust_conclusions, data.frame(
    module = "pathway",
    conclusion = sprintf("04_2 中 %d 个通路是 batch artifact",
                          pathway_summary$batch_driven),
    evidence = "控 batch 后消失"))
}
if (!is.null(ct_summary) && ct_summary$robust == 0 && ct_summary$nobatch_sig == 0) {
  robust_conclusions <- rbind(robust_conclusions, data.frame(
    module = "celltype",
    conclusion = "celltype 层面本无显著信号（控 batch 前后一致）",
    evidence = "0/0"))
}
save_table(robust_conclusions, DIR_THIS, "robust_conclusions.csv")

cat("\n=== 稳健结论清单 ===\n")
print(robust_conclusions, row.names = FALSE)

# ============================================================
# 8. figure
# ============================================================
log_info("--- 8. 图 ---")

bs_long <- batch_summary[, c("module","nobatch_sig","batch_sig")]
bs_long <- tidyr::pivot_longer(bs_long,
                                cols = c("nobatch_sig","batch_sig"),
                                names_to = "model", values_to = "n_sig")
bs_long$model <- factor(bs_long$model,
                         levels = c("nobatch_sig","batch_sig"),
                         labels = c("Without batch","With batch"))

p <- ggplot(bs_long, aes(x = module, y = n_sig, fill = model)) +
  geom_col(position = "dodge", alpha = 0.8) +
  theme_bw(base_size = 11) +
  labs(x = "", y = "Number of FDR<0.05 features",
       fill = "", title = "Batch sensitivity: without vs with batch covariate")
save_fig(p, DIR_THIS, "batch_effect_summary", width = 8, height = 5)

pathway_cmp$color <- ifelse(pathway_cmp$batch_driven, "batch-driven",
                     ifelse(pathway_cmp$robust, "robust", "ns"))
p2 <- ggplot(pathway_cmp,
              aes(x = -log10(FDR_nobatch + 1e-10),
                  y = -log10(FDR_batch + 1e-10),
                  color = color)) +
  geom_point(size = 1.5, alpha = 0.7) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "gray40") +
  scale_color_manual(values = c("batch-driven" = "#C44E52",
                                 "robust" = "#4C72B0",
                                 "ns" = "#B0B0B0")) +
  theme_bw(base_size = 11) +
  labs(x = "-log10(FDR) without batch",
       y = "-log10(FDR) with batch",
       color = "",
       title = "Pathway batch sensitivity (g3_vs_g2)")
save_fig(p2, DIR_THIS, "pathway_batch_scatter", width = 7, height = 6)

# ============================================================
# 9. summary
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("=== diag_batch_sensitivity.r 汇总 ===\n")
cat(strrep("=", 70), "\n", sep = "")

cat("\n各模块 batch 敏感性:\n")
print(batch_summary, row.names = FALSE)

cat("\n--- 输出位置 ---\n")
cat("目录:", DIR_THIS, "\n")

log_info("=== diag_batch_sensitivity.r done ===")
