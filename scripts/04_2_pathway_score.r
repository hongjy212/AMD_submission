#!/usr/bin/env Rscript
# ============================================================
# 04_2_pathway_score.r (v3: 419, layer + batch adjusted)
#
#
#
# Methods：
#
# Input：
#   results/01_1b_batch_enhance/data/{expr_wgcna,clinical_wgcna}.rds   ★ 419
#   results/01_1_prepare_bulk_data/data/gene_annotation.rds
#   results/01_2_prepare_gene_sets/data/gene_sets.rds
#   data/pathway_scores_ssgsea.rds
#   data/pathway_scores_zscore.rds
#   tables/gene_sets_used.csv
#   tables/ssgsea_vs_zscore_correlation.csv
#   tables/pathway_limma_multi_adjusted.csv
#   tables/pathway_age_comparison_adjusted.csv
#   tables/pathway_age_classification_adjusted.csv
#   tables/pathway_hUSI_correlation.csv
#   tables/batch_within_g3v2.csv
#   tables/mito_subfunc_summary_adjusted.csv
#   tables/cgast_amp_summary_adjusted.csv
#   tables/key_pathway_focus_adjusted.csv
#   figures/pathway_poly3F_bar_adjusted.pdf/.png
#   figures/pathway_grade_curves_adjusted.pdf/.png
#   figures/pathway_heatmap_adjusted.pdf/.png
#   figures/mito_subfunc_curves_adjusted.pdf/.png
#   figures/cgast_amp_curves_adjusted.pdf/.png
#   figures/batch_within_g3v2.pdf/.png
# ============================================================

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(limma); library(splines); library(ggplot2)
  library(GSVA); library(pheatmap); library(tidyr)
})

DIR_THIS <- setup_output_dir("04_2_pathway_score")
log_info("=== 04_2_pathway_score.r (v3, 419, layer+batch adjusted) start ===")

# Parameters
MIN_GENES <- 10
MAX_GENES <- 500

# ============================================================
# ============================================================
log_info("--- 1. 载入数据 ---")

expr <- read_upstream("01_1b_batch_enhance", "data", "expr_wgcna.rds")
clin <- read_upstream("01_1b_batch_enhance", "data", "clinical_wgcna.rds")
anno <- read_upstream("01_1_prepare_bulk_data", "data", "gene_annotation.rds")
GENE_SETS <- read_upstream("01_2_prepare_gene_sets", "data", "gene_sets.rds")

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

sid <- as.character(clin$sample_id)
clin$layer <- sub(".*_(L[0-9]+)$", "\\1", sid)
clin$layer <- factor(clin$layer, levels = paste0("L", 1:8))

rownames(clin) <- colnames(expr)

log_info(paste0("expr: ", nrow(expr), " × ", ncol(expr)))
log_info(paste0("gene sets: ", length(GENE_SETS)))
log_info(paste0("grade: ", paste(names(table(clin$grade)), table(clin$grade), sep="=", collapse=" ")))
log_info(paste0("batch: ", paste(names(table(clin$batch_f)), table(clin$batch_f), sep="=", collapse=" ")))
log_info(paste0("layer: ", paste(names(table(clin$layer)), table(clin$layer), sep="=", collapse=" ")))

# ============================================================
# 2. filter gene sets
# ============================================================
log_info(paste0("--- 2. 过滤 gene sets (", MIN_GENES, " ≤ n ≤ ", MAX_GENES, ") ---"))

gs_keep <- list()
gs_meta <- data.frame()

for (gs_name in names(GENE_SETS)) {
  gs <- GENE_SETS[[gs_name]]
  genes_in_expr <- intersect(gs$genes, rownames(expr))
  n <- length(genes_in_expr)

  if (n >= MIN_GENES && n <= MAX_GENES) {
    gs_keep[[gs_name]] <- genes_in_expr
    gs_meta <- rbind(gs_meta, data.frame(
      geneset = gs_name,
      group = if (!is.null(gs$group)) gs$group else NA,
      source = if (!is.null(gs$source)) gs$source else NA,
      source_id = if (!is.null(gs$source_id)) gs$source_id else NA,
      n_genes = n,
      stringsAsFactors = FALSE))
  }
}

log_info(paste0("保留: ", length(gs_keep), " / ", length(GENE_SETS)))
log_info("按 group:")
print(table(gs_meta$group))
save_table(gs_meta, DIR_THIS, "gene_sets_used.csv")

# ============================================================
# ============================================================
log_info("--- 3. ssGSEA 打分 ---")

set.seed(SEED)
ssgsea_param <- gsvaParam(
  exprData = as.matrix(expr),
  geneSets = gs_keep,
  kcdf = "Gaussian",
  minSize = MIN_GENES,
  maxSize = MAX_GENES
)
ssgsea_scores <- gsva(ssgsea_param)
log_info(paste0("ssGSEA 矩阵: ", nrow(ssgsea_scores), " × ", ncol(ssgsea_scores)))
save_rds(ssgsea_scores, DIR_THIS, "pathway_scores_ssgsea.rds")

# ============================================================
# ============================================================
log_info("--- 4. z-score 平均打分 ---")

expr_z <- t(scale(t(expr)))
expr_z[is.na(expr_z)] <- 0

zscore_mat <- do.call(rbind, lapply(gs_keep, function(genes) {
  colMeans(expr_z[genes, , drop = FALSE])
}))
save_rds(zscore_mat, DIR_THIS, "pathway_scores_zscore.rds")

common_gs <- intersect(rownames(ssgsea_scores), rownames(zscore_mat))
if (length(common_gs) >= 5) {
  r_cmp <- sapply(common_gs, function(gs) {
    cor(ssgsea_scores[gs, ], zscore_mat[gs, ])
  })
  save_table(data.frame(geneset = common_gs, r_ssgsea_zscore = round(r_cmp, 3)),
              DIR_THIS, "ssgsea_vs_zscore_correlation.csv")
  log_info(paste0("ssGSEA vs z-score 中位相关: ", round(median(r_cmp, na.rm=TRUE), 3)))
}

score_mat <- ssgsea_scores

# ============================================================
# ============================================================
make_covar <- function(clin) {
  sp <- splines::ns(clin$age_v, df = 3)
  covar <- data.frame(
    age_s1 = sp[, 1], age_s2 = sp[, 2], age_s3 = sp[, 3],
    sex_M = as.integer(clin$sex == "Male"),
    rin = clin$rin_v, pmi = clin$pmi_v,
    layer = clin$layer,
    batch = clin$batch_f,
    stringsAsFactors = FALSE)
  rownames(covar) <- rownames(clin)
  covar
}

run_limma_multi <- function(score_mat, clin, contrasts = NULL) {
  covar <- make_covar(clin)
  design_df <- data.frame(grade = factor(clin$grade), covar)
  rownames(design_df) <- rownames(clin)
  design_df <- design_df[complete.cases(design_df), , drop = FALSE]
  score_sub <- score_mat[, rownames(design_df), drop = FALSE]

  design_multi <- model.matrix(
    ~ 0 + grade + age_s1 + age_s2 + age_s3 + sex_M + rin + pmi + layer + batch,
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

  out <- data.frame()
  for (ct in colnames(cm)) {
    tt <- topTable(fit_multi, coef = ct, number = Inf, sort.by = "none")
    out <- rbind(out, data.frame(
      pathway = rownames(tt), contrast = ct,
      logFC = tt$logFC, t = tt$t, p = tt$P.Value, FDR = tt$adj.P.Val,
      stringsAsFactors = FALSE))
  }
  out
}

# ============================================================
# ============================================================
log_info("--- 6. limma multi (7 对比, layer+batch adjusted) ---")

gs_info <- gs_meta[match(rownames(score_mat), gs_meta$geneset), ]
multi_res <- run_limma_multi(score_mat, clin)
multi_res$group <- gs_info$group[match(multi_res$pathway, gs_info$geneset)]
multi_res$source_id <- gs_info$source_id[match(multi_res$pathway, gs_info$geneset)]
save_table(multi_res, DIR_THIS, "pathway_limma_multi_adjusted.csv")

cat("\n=== limma multi 各对比 FDR<0.05 (layer+batch adjusted) ===\n")
for (ct in unique(multi_res$contrast)) {
  n <- sum(multi_res$FDR[multi_res$contrast == ct] < 0.05, na.rm = TRUE)
  cat(sprintf("  %-15s  %d / %d\n", ct, n, nrow(score_mat)))
}

sub <- multi_res[multi_res$contrast == "g3_vs_g2" & multi_res$FDR < 0.05, ]
if (nrow(sub) > 0) {
  cat("\n=== g3_vs_g2 FDR<0.05 明细 ===\n")
  print(sub[order(sub$FDR), c("pathway","group","logFC","p","FDR")],
        row.names = FALSE, digits = 4)
}

# ============================================================
# ============================================================
log_info("--- 7. poly3_F (layer+batch adjusted) ---")

covar <- make_covar(clin)
log_info(paste0("covar: ", nrow(covar), " × ", ncol(covar)))

poly_mat <- poly(clin$grade, 3)
colnames(poly_mat) <- c("grade_p1", "grade_p2", "grade_p3")
rownames(poly_mat) <- rownames(clin)

design_poly_df <- data.frame(poly_mat, covar)
rownames(design_poly_df) <- rownames(clin)
design_poly_df <- design_poly_df[complete.cases(design_poly_df), , drop = FALSE]

score_poly <- score_mat[, rownames(design_poly_df), drop = FALSE]

design_poly <- model.matrix(
  ~ 0 + grade_p1 + grade_p2 + grade_p3 +
      age_s1 + age_s2 + age_s3 + sex_M + rin + pmi + layer + batch,
  data = design_poly_df)
rownames(design_poly) <- rownames(design_poly_df)

log_info(paste0("poly3 design: ", nrow(design_poly), " × ", ncol(design_poly)))

# 7.2 limma fit
fit_poly <- eBayes(lmFit(score_poly, design_poly))

coef_poly <- grep("^grade_p", colnames(design_poly), value = TRUE)
tt_F <- topTable(fit_poly, coef = coef_poly, number = Inf, sort.by = "none")

tt_p1 <- topTable(fit_poly, coef = "grade_p1", number = Inf, sort.by = "none")
tt_p2 <- topTable(fit_poly, coef = "grade_p2", number = Inf, sort.by = "none")
tt_p3 <- topTable(fit_poly, coef = "grade_p3", number = Inf, sort.by = "none")

polyF <- data.frame(
  pathway = rownames(tt_F),
  F_stat  = tt_F$F,
  p       = tt_F$P.Value,
  FDR     = tt_F$adj.P.Val,
  p1      = tt_p1$P.Value,
  p2      = tt_p2$P.Value,
  p3      = tt_p3$P.Value,
  p1_FDR  = p.adjust(tt_p1$P.Value, method = "BH"),
  p2_FDR  = p.adjust(tt_p2$P.Value, method = "BH"),
  p3_FDR  = p.adjust(tt_p3$P.Value, method = "BH"),
  p1_coef = tt_p1$logFC,
  p2_coef = tt_p2$logFC,
  p3_coef = tt_p3$logFC,
  stringsAsFactors = FALSE
)

polyF$pattern <- with(polyF,
  ifelse(!is.na(p3) & p3 < 0.05 & p1 > 0.05 & p2 > 0.05, "p3_only",
  ifelse(!is.na(p3) & p3 < 0.05, "p3+other",
  ifelse(!is.na(p1) & p1 < 0.05, "p1_sig",
  ifelse(!is.na(p2) & p2 < 0.05, "p2_sig", "none")))))

# 7.7 sort + save
polyF <- polyF[order(polyF$p3), ]

save_table(polyF, DIR_THIS, "pathway_poly3F_adjusted.csv")
save_table(polyF[, c("pathway", "F_stat", "p", "FDR")],
           DIR_THIS, "pathway_poly3F_adjusted_original.csv")

# 7.8 summary
log_info(paste0("poly3_F 整体 FDR<0.05: ",
                sum(polyF$FDR < 0.05, na.rm = TRUE), " / ", nrow(polyF)))

n_p3_only <- sum(polyF$pattern == "p3_only", na.rm = TRUE)
n_p1_sig  <- sum(polyF$pattern == "p1_sig",  na.rm = TRUE)
log_info(paste0("p3_only: ", n_p3_only, " / ", nrow(polyF)))
log_info(paste0("p1_sig : ", n_p1_sig,  " / ", nrow(polyF)))

axon_paths <- c("Axonal_Transport", "Axonal_Growth_Cone",
                "Growth_Cone", "Optic_Nerve_Development")
mito_paths <- grep("OXPHOS|Mito_CI|Mito_CII|Mito_CIII|Mito_CIV|Mito_CV|Mito_TCA",
                    polyF$pathway, value = TRUE)

cat("\n=== 轴突通路 poly3 分项 ===\n")
print(polyF[polyF$pathway %in% axon_paths,
             c("pathway", "p1", "p2", "p3", "pattern")], row.names = FALSE)

cat("\n=== 线粒体/OXPHOS 通路 poly3 分项 ===\n")
print(polyF[polyF$pathway %in% mito_paths,
             c("pathway", "p1", "p2", "p3", "pattern")], row.names = FALSE)

cat("\n=== 全部 p3_only 通路 ===\n")
print(polyF[polyF$pattern == "p3_only",
             c("pathway", "p1", "p2", "p3")], row.names = FALSE)

cat("\n=== 全部 p1_sig 通路 ===\n")
print(polyF[polyF$pattern == "p1_sig",
             c("pathway", "p1", "p2", "p3")], row.names = FALSE)

# ============================================================
# 8. agediagnosis
# ============================================================
log_info("--- 8. 年龄诊断 ---")

run_pathway_age <- function(score_mat, clin, age_mode) {
  covar <- data.frame(
    sex_M = as.integer(clin$sex == "Male"),
    rin = clin$rin_v, pmi = clin$pmi_v,
    layer = clin$layer, batch = clin$batch_f,
    stringsAsFactors = FALSE)
  if (age_mode == "linear") {
    covar$age <- clin$age_v
  } else if (age_mode == "spline") {
    sp <- splines::ns(clin$age_v, df = 3)
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
      pathway = rownames(tt), age_mode = age_mode, contrast = ct,
      logFC = tt$logFC, p = tt$P.Value, FDR = tt$adj.P.Val,
      stringsAsFactors = FALSE))
  }
  out
}

age_cmp <- rbind(
  run_pathway_age(score_mat, clin, "none"),
  run_pathway_age(score_mat, clin, "linear"),
  run_pathway_age(score_mat, clin, "spline"))
age_cmp$group <- gs_info$group[match(age_cmp$pathway, gs_info$geneset)]
save_table(age_cmp, DIR_THIS, "pathway_age_comparison_adjusted.csv")

age_g3v2 <- age_cmp[age_cmp$contrast == "g3_vs_g2", ]
stable_pathways <- age_g3v2$pathway[
  age_g3v2$age_mode == "spline" & age_g3v2$p < 0.05]
log_info(paste0("g3_vs_g2 控年龄后 p<0.05 的 pathway: ", length(stable_pathways)))

if (length(stable_pathways) > 0) {
  cat("\n=== g3_vs_g2 控年龄后显著 pathway ===\n")
  sig_tbl <- age_g3v2[age_g3v2$pathway %in% stable_pathways &
                        age_g3v2$age_mode == "spline", ]
  print(sig_tbl[order(sig_tbl$p), c("pathway","group","logFC","p","FDR")],
        row.names = FALSE, digits = 4)
}

# ============================================================
# ============================================================
log_info("--- 9. 年龄双模型分类 ---")

age_class_tbl <- data.frame()
for (i in seq_len(nrow(score_mat))) {
  sc_name <- rownames(score_mat)[i]
  y <- score_mat[i, ]

  df <- data.frame(y = y, grade = clin$grade, age = clin$age_v,
                    sex = clin$sex, rin = clin$rin_v, pmi = clin$pmi_v,
                    layer = clin$layer, batch = clin$batch_f)
  df <- df[complete.cases(df), ]

  sp_g <- splines::ns(df$age, df = 3)
  fit1 <- lm(y ~ poly(grade, 3) + sp_g + sex + rin + pmi + layer + batch,
              data = df)
  p_g <- anova(fit1)["poly(grade, 3)", "Pr(>F)"]

  fit2 <- lm(y ~ sp_g + poly(grade, 3) + sex + rin + pmi + layer + batch,
              data = df)
  p_a <- anova(fit2)["sp_g", "Pr(>F)"]

  age_class_tbl <- rbind(age_class_tbl, data.frame(
    pathway = sc_name, p_grade = p_g, p_age = p_a,
    stringsAsFactors = FALSE))
}

age_class_tbl$type <- "无关"
age_class_tbl$type[age_class_tbl$p_grade < 0.05 &
                   age_class_tbl$p_age > 0.10] <- "纯疾病"
age_class_tbl$type[age_class_tbl$p_grade < 0.05 &
                   age_class_tbl$p_age < 0.05] <- "年龄+疾病"
age_class_tbl$type[age_class_tbl$p_grade > 0.10 &
                   age_class_tbl$p_age < 0.05] <- "纯年龄"
age_class_tbl$group <- gs_info$group[match(age_class_tbl$pathway, gs_info$geneset)]
save_table(age_class_tbl, DIR_THIS, "pathway_age_classification_adjusted.csv")

cat("\n=== 年龄双模型分类 (layer+batch adjusted) ===\n")
print(table(age_class_tbl$type))

pure_disease <- age_class_tbl[age_class_tbl$type == "纯疾病", ]
if (nrow(pure_disease) > 0) {
  cat("\n=== 纯疾病 pathway ===\n")
  print(pure_disease, row.names = FALSE, digits = 4)
}

# ============================================================
# ============================================================
log_info("--- 10. 与 hUSI 相关 ---")

husi_file <- file.path(DIR_RESULT, "04_1_aging_score", "data", "aging_scores.rds")
if (file.exists(husi_file)) {
  aging_scores <- readRDS(husi_file)
  common_samples <- intersect(colnames(score_mat), colnames(aging_scores))

  if (length(common_samples) > 30) {
    hUSI_v <- aging_scores["hUSI", common_samples]
    cor_tbl <- data.frame()
    for (gs in rownames(score_mat)) {
      r <- cor(score_mat[gs, common_samples], hUSI_v, use = "complete.obs")
      cor_tbl <- rbind(cor_tbl, data.frame(
        pathway = gs,
        group = gs_info$group[match(gs, gs_info$geneset)],
        r_hUSI = round(r, 3),
        stringsAsFactors = FALSE))
    }
    cor_tbl <- cor_tbl[order(-abs(cor_tbl$r_hUSI)), ]
    save_table(cor_tbl, DIR_THIS, "pathway_hUSI_correlation.csv")

    cat("\n=== 与 hUSI 相关 top 20 ===\n")
    print(head(cor_tbl, 20), row.names = FALSE)
  }
} else {
  log_warn("04_1 aging_scores.rds 不存在，跳过")
}

# ============================================================
# ============================================================
log_info("--- 11. batch 内验证 ---")

run_within_batch <- function(score_mat, clin, batch_val) {
  idx <- clin$batch_f == batch_val
  sub_clin <- clin[idx, ]
  sub_score <- score_mat[, rownames(sub_clin), drop = FALSE]

  sp <- splines::ns(sub_clin$age_v, df = 3)
  covar <- data.frame(
    age_s1 = sp[, 1], age_s2 = sp[, 2], age_s3 = sp[, 3],
    sex_M = as.integer(sub_clin$sex == "Male"),
    rin = sub_clin$rin_v, pmi = sub_clin$pmi_v,
    layer = sub_clin$layer,
    stringsAsFactors = FALSE)
  rownames(covar) <- rownames(sub_clin)
  design_df <- data.frame(grade = factor(sub_clin$grade), covar)
  rownames(design_df) <- rownames(sub_clin)
  design_df <- design_df[complete.cases(design_df), , drop = FALSE]
  score_use <- sub_score[, rownames(design_df), drop = FALSE]

  design <- model.matrix(
    ~ 0 + grade + age_s1 + age_s2 + age_s3 + sex_M + rin + pmi + layer,
    data = design_df)
  colnames(design)[grep("^grade", colnames(design))] <- paste0("g", 1:4)
  cm <- limma::makeContrasts(g3_vs_g2 = g3 - g2, levels = design)
  fit <- eBayes(contrasts.fit(lmFit(score_use, design), cm))
  tt <- topTable(fit, coef = "g3_vs_g2", number = Inf, sort.by = "none")
  data.frame(feature = rownames(tt), batch = batch_val,
             logFC = tt$logFC, p = tt$P.Value, FDR = tt$adj.P.Val,
             stringsAsFactors = FALSE)
}

batches <- levels(clin$batch_f)
batch_within <- do.call(rbind, lapply(batches, function(b) {
  tryCatch(run_within_batch(score_mat, clin, b),
           error = function(e) { log_warn(paste0(b, ": ", e$message)); NULL })
}))
save_table(batch_within, DIR_THIS, "batch_within_g3v2.csv")

wide <- reshape(batch_within[, c("feature","batch","logFC","p")],
                 idvar = "feature", timevar = "batch", direction = "wide")
wide$direction_consistent <- sign(wide[[grep("^logFC", colnames(wide))[1]]]) ==
                              sign(wide[[grep("^logFC", colnames(wide))[2]]])
save_table(wide, DIR_THIS, "batch_within_g3v2_wide.csv")

cat("\n=== batch 内 g3_vs_g2 方向一致性 ===\n")
cat(sprintf("  两个 batch 内方向一致: %d / %d\n",
            sum(wide$direction_consistent, na.rm = TRUE), nrow(wide)))

wide$group <- gs_info$group[match(wide$feature, gs_info$geneset)]
cat("\n按 group 的方向一致性:\n")
print(table(wide$group, wide$direction_consistent))

# ============================================================
# ============================================================
log_info("--- 12. 专项汇总 ---")

mito_gs <- gs_meta$geneset[gs_meta$group == "mito_subfunc"]
if (length(mito_gs) > 0) {
  mito_summary <- data.frame()
  for (gs in mito_gs) {
    p3 <- polyF[polyF$pathway == gs, ]
    g3 <- age_g3v2[age_g3v2$pathway == gs & age_g3v2$age_mode == "spline", ]
    g4g3 <- multi_res[multi_res$pathway == gs & multi_res$contrast == "g4_vs_g3", ]
    w <- wide[wide$feature == gs, ]
    mito_summary <- rbind(mito_summary, data.frame(
      pathway = gs,
      n_genes = gs_meta$n_genes[gs_meta$geneset == gs],
      poly3F_p = if (nrow(p3) > 0) p3$p else NA,
      poly3F_FDR = if (nrow(p3) > 0) p3$FDR else NA,
      poly3F_p3 = if (nrow(p3) > 0) p3$p3 else NA,
      poly3F_pattern = if (nrow(p3) > 0) p3$pattern else NA,
      g3v2_spline_logFC = if (nrow(g3) > 0) g3$logFC else NA,
      g3v2_spline_p = if (nrow(g3) > 0) g3$p else NA,
      g4v3_logFC = if (nrow(g4g3) > 0) g4g3$logFC else NA,
      g4v3_FDR = if (nrow(g4g3) > 0) g4g3$FDR else NA,
      direction_consistent = if (nrow(w) > 0) w$direction_consistent else NA,
      stringsAsFactors = FALSE))
  }
  mito_summary <- mito_summary[order(mito_summary$poly3F_p), ]
  save_table(mito_summary, DIR_THIS, "mito_subfunc_summary_adjusted.csv")
  cat("\n=== MitoCarta 子功能 (layer+batch adjusted) ===\n")
  print(mito_summary, row.names = FALSE, digits = 4)
}

# ---- cGAS-STING ----
cgast_gs <- gs_meta$geneset[gs_meta$group == "cgast_amp"]
if (length(cgast_gs) > 0) {
  cgast_summary <- data.frame()
  for (gs in cgast_gs) {
    p3 <- polyF[polyF$pathway == gs, ]
    g3 <- age_g3v2[age_g3v2$pathway == gs & age_g3v2$age_mode == "spline", ]
    g4g3 <- multi_res[multi_res$pathway == gs & multi_res$contrast == "g4_vs_g3", ]
    cgast_summary <- rbind(cgast_summary, data.frame(
      pathway = gs,
      n_genes = gs_meta$n_genes[gs_meta$geneset == gs],
      poly3F_p = if (nrow(p3) > 0) p3$p else NA,
      poly3F_p3 = if (nrow(p3) > 0) p3$p3 else NA,
      poly3F_pattern = if (nrow(p3) > 0) p3$pattern else NA,
      g3v2_spline_logFC = if (nrow(g3) > 0) g3$logFC else NA,
      g3v2_spline_p = if (nrow(g3) > 0) g3$p else NA,
      g4v3_logFC = if (nrow(g4g3) > 0) g4g3$logFC else NA,
      g4v3_FDR = if (nrow(g4g3) > 0) g4g3$FDR else NA,
      stringsAsFactors = FALSE))
  }
  save_table(cgast_summary, DIR_THIS, "cgast_amp_summary_adjusted.csv")
  cat("\n=== cGAS-STING (layer+batch adjusted) ===\n")
  print(cgast_summary, row.names = FALSE, digits = 4)
}

key_groups <- c("amd_pathology", "cgast_amp", "visual", "axon",
                "artifact", "metabolism", "mito_subfunc")
key_genesets <- gs_meta$geneset[gs_meta$group %in% key_groups]
key_focus <- data.frame()
for (gs in key_genesets) {
  p3 <- polyF[polyF$pathway == gs, ]
  g3 <- age_g3v2[age_g3v2$pathway == gs & age_g3v2$age_mode == "spline", ]
  g4g3 <- multi_res[multi_res$pathway == gs & multi_res$contrast == "g4_vs_g3", ]
  w <- wide[wide$feature == gs, ]
  key_focus <- rbind(key_focus, data.frame(
    pathway = gs,
    group = gs_meta$group[gs_meta$geneset == gs],
    n_genes = gs_meta$n_genes[gs_meta$geneset == gs],
    poly3F_p = if (nrow(p3) > 0) p3$p else NA,
    poly3F_FDR = if (nrow(p3) > 0) p3$FDR else NA,
    p1 = if (nrow(p3) > 0) p3$p1 else NA,
    p2 = if (nrow(p3) > 0) p3$p2 else NA,
    p3 = if (nrow(p3) > 0) p3$p3 else NA,
    pattern = if (nrow(p3) > 0) p3$pattern else NA,
    g3v2_spline_logFC = if (nrow(g3) > 0) g3$logFC else NA,
    g3v2_spline_p = if (nrow(g3) > 0) g3$p else NA,
    g4v3_logFC = if (nrow(g4g3) > 0) g4g3$logFC else NA,
    g4v3_FDR = if (nrow(g4g3) > 0) g4g3$FDR else NA,
    direction_consistent = if (nrow(w) > 0) w$direction_consistent else NA,
    stringsAsFactors = FALSE))
}
key_focus <- key_focus[order(key_focus$poly3F_p), ]
save_table(key_focus, DIR_THIS, "key_pathway_focus_adjusted.csv")
cat("\n=== 关键通路聚焦 top 30 (layer+batch adjusted) ===\n")
print(head(key_focus, 30), row.names = FALSE, digits = 4)

# ============================================================
# 13. figure
# ============================================================
log_info("--- 13. 图 ---")

top20 <- head(polyF, 20)
top20$pathway <- factor(top20$pathway,
                         levels = top20$pathway[order(top20$F_stat)])
top20$sig <- ifelse(top20$FDR < 0.05, "FDR<0.05",
             ifelse(top20$p < 0.05, "p<0.05", "ns"))
p1 <- ggplot(top20, aes(x = pathway, y = F_stat, fill = sig)) +
  geom_col() + coord_flip() + theme_bw(base_size = 10) +
  scale_fill_manual(values = c("FDR<0.05" = "#C44E52",
                                "p<0.05" = "#FDB462",
                                "ns" = "#B0B0B0")) +
  labs(x = "", y = "poly3_F F statistic", fill = "",
       title = "Pathway poly3_F top 20 (layer+batch adjusted)")
save_fig(p1, DIR_THIS, "pathway_poly3F_bar_adjusted", width = 8, height = 7)

top30 <- head(polyF$pathway, 30)
curve_data <- data.frame()
for (gs in top30) {
  if (!gs %in% rownames(score_mat)) next
  for (g in 1:4) {
    idx <- clin$grade == g
    curve_data <- rbind(curve_data, data.frame(
      pathway = gs, grade = g,
      mean_score = mean(score_mat[gs, idx], na.rm = TRUE),
      se = sd(score_mat[gs, idx], na.rm = TRUE) / sqrt(sum(idx)),
      stringsAsFactors = FALSE))
  }
}
p2 <- ggplot(curve_data, aes(x = grade, y = mean_score,
                               color = pathway, group = pathway)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_score - se, ymax = mean_score + se),
                 width = 0.15, alpha = 0.5) +
  theme_bw(base_size = 11) +
  scale_x_continuous(breaks = 1:4, labels = paste0("G", 1:4)) +
  labs(x = "AMD grade", y = "ssGSEA score (mean ± SE)",
       title = "Top 30 pathways by poly3_F (layer+batch adjusted)",
       color = "Pathway")
save_fig(p2, DIR_THIS, "pathway_grade_curves_adjusted", width = 12, height = 8)

top20_names <- head(polyF$pathway, 20)
top20_names <- top20_names[top20_names %in% rownames(score_mat)]
heat_mat <- matrix(NA_real_, nrow = length(top20_names), ncol = 4,
                    dimnames = list(top20_names, paste0("G", 1:4)))
for (gs in top20_names) {
  for (g in 1:4) {
    heat_mat[gs, g] <- mean(score_mat[gs, clin$grade == g])
  }
}
pheatmap(heat_mat, cluster_rows = TRUE, cluster_cols = FALSE,
         display_numbers = TRUE, number_format = "%.2f",
         color = colorRampPalette(c("#4C72B0","#FFFFFF","#C44E52"))(50),
         main = "Top 20 pathways × grade (layer+batch adjusted)",
         filename = file.path(DIR_THIS, "figures",
                              "pathway_heatmap_adjusted.png"),
         width = 7, height = 8)

if (length(mito_gs) > 0) {
  mito_curve <- data.frame()
  for (gs in mito_gs) {
    if (!gs %in% rownames(score_mat)) next
    for (g in 1:4) {
      idx <- clin$grade == g
      mito_curve <- rbind(mito_curve, data.frame(
        pathway = gs, grade = g,
        mean_score = mean(score_mat[gs, idx], na.rm = TRUE),
        se = sd(score_mat[gs, idx], na.rm = TRUE) / sqrt(sum(idx)),
        stringsAsFactors = FALSE))
    }
  }
  p_mito <- ggplot(mito_curve, aes(x = grade, y = mean_score,
                                     color = pathway, group = pathway)) +
    geom_line(linewidth = 0.8) + geom_point(size = 1.8) +
    geom_errorbar(aes(ymin = mean_score - se, ymax = mean_score + se),
                   width = 0.15, alpha = 0.5) +
    theme_bw(base_size = 10) +
    scale_x_continuous(breaks = 1:4, labels = paste0("G", 1:4)) +
    labs(x = "AMD grade", y = "ssGSEA score (mean ± SE)",
         title = "MitoCarta subfunction pathways (layer+batch adjusted)",
         color = "Pathway")
  save_fig(p_mito, DIR_THIS, "mito_subfunc_curves_adjusted",
           width = 12, height = 8)
}

if (length(cgast_gs) > 0) {
  cgast_curve <- data.frame()
  for (gs in cgast_gs) {
    if (!gs %in% rownames(score_mat)) next
    for (g in 1:4) {
      idx <- clin$grade == g
      cgast_curve <- rbind(cgast_curve, data.frame(
        pathway = gs, grade = g,
        mean_score = mean(score_mat[gs, idx], na.rm = TRUE),
        se = sd(score_mat[gs, idx], na.rm = TRUE) / sqrt(sum(idx)),
        stringsAsFactors = FALSE))
    }
  }
  p_cgast <- ggplot(cgast_curve, aes(x = grade, y = mean_score,
                                       color = pathway, group = pathway)) +
    geom_line(linewidth = 1) + geom_point(size = 2.5) +
    geom_errorbar(aes(ymin = mean_score - se, ymax = mean_score + se),
                   width = 0.15, alpha = 0.5) +
    theme_bw(base_size = 11) +
    scale_x_continuous(breaks = 1:4, labels = paste0("G", 1:4)) +
    labs(x = "AMD grade", y = "ssGSEA score (mean ± SE)",
         title = "cGAS-STING pathways (layer+batch adjusted)",
         color = "Pathway")
  save_fig(p_cgast, DIR_THIS, "cgast_amp_curves_adjusted", width = 8, height = 5)
}

wide_plot <- wide[!is.na(wide$direction_consistent), ]
if (nrow(wide_plot) > 0) {
  colnames_logFC <- grep("^logFC", colnames(wide_plot), value = TRUE)
  if (length(colnames_logFC) >= 2) {
    wide_plot$x <- wide_plot[[colnames_logFC[1]]]
    wide_plot$y <- wide_plot[[colnames_logFC[2]]]
    wide_plot$color <- ifelse(wide_plot$direction_consistent,
                                "consistent", "opposite")
    p_batch <- ggplot(wide_plot, aes(x = x, y = y, color = color)) +
      geom_point(size = 1.5, alpha = 0.7) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray60") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
      scale_color_manual(values = c("consistent" = "#4C72B0",
                                     "opposite" = "#C44E52")) +
      theme_bw(base_size = 11) +
      labs(x = paste0("logFC in ", batches[1]),
           y = paste0("logFC in ", batches[2]),
           color = "",
           title = "Batch within-batch g3_vs_g2 (direction consistency)")
    save_fig(p_batch, DIR_THIS, "batch_within_g3v2", width = 7, height = 6)
  }
}

# ============================================================
# 14. summary
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("=== 04_2 pathway score 汇总 (v3, 419, layer+batch adjusted) ===\n")
cat(strrep("=", 70), "\n", sep = "")
cat("pathway 数:", nrow(score_mat), "\n")
cat("样本数:", ncol(score_mat), "\n\n")

cat("--- limma multi 各对比 FDR<0.05 ---\n")
for (ct in unique(multi_res$contrast)) {
  n <- sum(multi_res$FDR[multi_res$contrast == ct] < 0.05, na.rm = TRUE)
  cat(sprintf("  %-15s  %d\n", ct, n))
}

cat("\n--- poly3_F 分项汇总 ---\n")
cat("整体 FDR<0.05:", sum(polyF$FDR < 0.05, na.rm = TRUE), " / ", nrow(polyF), "\n")
cat("p3_only      :", sum(polyF$pattern == "p3_only", na.rm = TRUE), " / ", nrow(polyF), "\n")
cat("p3+other     :", sum(polyF$pattern == "p3+other", na.rm = TRUE), " / ", nrow(polyF), "\n")
cat("p1_sig       :", sum(polyF$pattern == "p1_sig", na.rm = TRUE), " / ", nrow(polyF), "\n")
cat("p2_sig       :", sum(polyF$pattern == "p2_sig", na.rm = TRUE), " / ", nrow(polyF), "\n")
cat("none         :", sum(polyF$pattern == "none", na.rm = TRUE), " / ", nrow(polyF), "\n")

cat("\n--- 年龄双模型分类 ---\n")
print(table(age_class_tbl$type))

cat("\n--- batch 内方向一致性 ---\n")
cat(sprintf("  一致: %d / %d\n",
            sum(wide$direction_consistent, na.rm = TRUE), nrow(wide)))

cat("\n--- 输出位置 ---\n")
cat("目录:", DIR_THIS, "\n")
log_info("=== 04_2_pathway_score.r (v3, 419) done ===")

