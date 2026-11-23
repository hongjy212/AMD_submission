#!/usr/bin/env Rscript
# ============================================================
#
#
# Input：
#   results/01_1_prepare_bulk_data/data/{expr_corrected,clinical,gene_annotation}.rds
# Output：
#   results/02_1_453_discovery/tables/
#     sig16_poly3_detail.csv
#     sig16_genes.csv
#     poly3_F_full.csv
#     genes16_poly3F_adjusted.csv
# ============================================================

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(limma); library(splines); library(dplyr)
})

ROOT <- "E:/retina_tools/aging_rgc_axon/AMD"
OUT  <- file.path(ROOT, "results", "02_1_453_discovery")
TAB  <- file.path(OUT, "tables")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)

log_info("=== 02_1_453_discovery.r v2 start ===")

# ============================================================
# ============================================================
expr <- readRDS(file.path(ROOT, "results", "01_1_prepare_bulk_data",
                           "data", "expr_corrected.rds"))
clin <- readRDS(file.path(ROOT, "results", "01_1_prepare_bulk_data",
                           "data", "clinical.rds"))
anno <- readRDS(file.path(ROOT, "results", "01_1_prepare_bulk_data",
                           "data", "gene_annotation.rds"))

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
log_info("grade: {paste(names(table(clin$grade)), table(clin$grade), sep='=', collapse=' ')}")

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

# ============================================================
# 3. poly3_F
# ============================================================
poly_mat <- poly(clin$grade, 3)
colnames(poly_mat) <- c("grade_p1", "grade_p2", "grade_p3")
rownames(poly_mat) <- rownames(clin)
design_poly_df <- data.frame(poly_mat, covar_df)
rownames(design_poly_df) <- rownames(clin)
design_poly_df <- design_poly_df[complete.cases(design_poly_df), , drop = FALSE]
expr_poly <- expr[, rownames(design_poly_df), drop = FALSE]

design_poly <- model.matrix(
  ~ grade_p1 + grade_p2 + grade_p3 + age_s1 + age_s2 + age_s3 + sex_M + rin + pmi,
  data = design_poly_df)
fit_poly <- eBayes(lmFit(expr_poly, design_poly))

coef_poly <- grep("^grade_p", colnames(design_poly), value = TRUE)
tt_F <- topTable(fit_poly, coef = coef_poly, number = Inf, sort.by = "none")
polyF <- data.frame(
  gene = rownames(tt_F),
  t = tt_F$F, p = tt_F$P.Value, FDR = tt_F$adj.P.Val,
  stringsAsFactors = FALSE)

write.csv(polyF, file.path(TAB, "poly3_F_full.csv"), row.names = FALSE)
log_info("saved: poly3_F_full.csv")

log_info("poly3_F (453): FDR<0.05 = {sum(polyF$FDR < 0.05, na.rm=TRUE)}")

# ============================================================
# ============================================================
sig16 <- polyF$gene[polyF$FDR < 0.05]
sig16_detail <- polyF[polyF$gene %in% sig16, ]
sig16_detail <- sig16_detail[order(sig16_detail$FDR), ]

write.csv(data.frame(gene = sig16),
          file.path(TAB, "sig16_genes.csv"), row.names = FALSE)
write.csv(sig16_detail,
          file.path(TAB, "sig16_poly3_detail.csv"), row.names = FALSE)

log_info("poly3_F FDR<0.05: {length(sig16)} 基因")
cat("\n=== 453 上的 16 基因 ===\n")
print(sig16_detail, row.names = FALSE)

# ============================================================
# ============================================================
run_genes16_poly3F <- function(expr, clin, model_name, use_layer_batch) {
  covar <- data.frame(
    age_s1 = splines::ns(clin$age_n, df = 3)[, 1],
    age_s2 = splines::ns(clin$age_n, df = 3)[, 2],
    age_s3 = splines::ns(clin$age_n, df = 3)[, 3],
    sex_M = as.integer(clin$sex == "Male"),
    rin = clin$rin_n, pmi = clin$pmi_n,
    stringsAsFactors = FALSE)
  if (use_layer_batch) {
    sid <- as.character(clin$sample_id)
    layer <- sub(".*_(L[0-9]+)$", "\\1", sid)
    covar$layer <- factor(layer, levels = paste0("L", 1:8))
    covar$batch <- factor(clin$batch)
  }
  rownames(covar) <- rownames(clin)

  poly_mat <- poly(clin$grade, 3)
  colnames(poly_mat) <- c("grade_p1", "grade_p2", "grade_p3")
  rownames(poly_mat) <- rownames(clin)

  design_df <- data.frame(poly_mat, covar)
  rownames(design_df) <- rownames(clin)
  design_df <- design_df[complete.cases(design_df), , drop = FALSE]

  expr_16 <- expr[sig16, rownames(design_df), drop = FALSE]

  fml <- as.formula(paste("~", paste(colnames(design_df), collapse = " + ")))
  design <- model.matrix(fml, data = design_df)
  fit <- eBayes(lmFit(expr_16, design))

  coef_poly <- c("grade_p1", "grade_p2", "grade_p3")
  tt_F <- topTable(fit, coef = coef_poly, number = Inf, sort.by = "none")

  coef_p1 <- topTable(fit, coef = "grade_p1", number = Inf, sort.by = "none")
  coef_p2 <- topTable(fit, coef = "grade_p2", number = Inf, sort.by = "none")
  coef_p3 <- topTable(fit, coef = "grade_p3", number = Inf, sort.by = "none")

  data.frame(
    gene = rownames(tt_F),
    model = model_name,
    F_stat = tt_F$F,
    p = tt_F$P.Value,
    FDR = tt_F$adj.P.Val,
    p1_logFC = coef_p1$logFC, p1_p = coef_p1$P.Value, p1_FDR = coef_p1$adj.P.Val,
    p2_logFC = coef_p2$logFC, p2_p = coef_p2$P.Value, p2_FDR = coef_p2$adj.P.Val,
    p3_logFC = coef_p3$logFC, p3_p = coef_p3$P.Value, p3_FDR = coef_p3$adj.P.Val,
    stringsAsFactors = FALSE)
}

genes16_poly3F <- rbind(
  run_genes16_poly3F(expr, clin, "M1_age_PMI_RIN", FALSE),
  run_genes16_poly3F(expr, clin, "M1_plus_layer_batch", TRUE)
)
write.csv(genes16_poly3F,
          file.path(TAB, "genes16_poly3F_adjusted.csv"), row.names = FALSE)
log_info("saved: genes16_poly3F_adjusted.csv")

log_info("=== 02_1_453_discovery.r v2 done ===")


