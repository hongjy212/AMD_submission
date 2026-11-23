#!/usr/bin/env Rscript
# =============================================================================
# 04_2c_pathway_poly3_cellprop.r (v4 final)
#
# =============================================================================

rm(list = ls()); gc()
setwd("E:/retina_tools/aging_rgc_axon/AMD")
set.seed(2026)

suppressPackageStartupMessages({
  library(limma); library(splines); library(dplyr)
  library(tidyr); library(ggplot2); library(logger)
})

out_dir <- "results/04_2c_cellprop_adjusted"
tab_dir <- file.path(out_dir, "tables")
fig_dir <- file.path(out_dir, "figures")
for (d in c(out_dir, tab_dir, fig_dir))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

log_appender(appender_tee(file.path(out_dir, "04_2c_v4.log")))
log_layout(layout_glue_generator(format = '{time} [{level}] {msg}'))
log_info("=== 04_2c v4 final ===")

# ============================================================
# ============================================================
score_mat <- readRDS("results/04_2_pathway_score/data/pathway_scores_ssgsea.rds")
cell_prop <- read.csv("results/03_2_deconvolution/tables/cell_proportions.csv")
clin <- readRDS("results/01_1_prepare_bulk_data/data/clinical.rds")

common <- Reduce(intersect, list(colnames(score_mat), rownames(clin),
                                  cell_prop$sample_id))
score_mat <- score_mat[, common, drop = FALSE]
clin_sub <- clin[common, , drop = FALSE]
cell_prop_sub <- cell_prop[match(common, cell_prop$sample_id), ]
rownames(cell_prop_sub) <- common

clin_sub$grade   <- as.integer(clin_sub$grade)
clin_sub$age_v   <- as.numeric(clin_sub$age_n)
clin_sub$rin_v   <- as.numeric(clin_sub$rin_n)
clin_sub$pmi_v   <- as.numeric(clin_sub$pmi_n)
clin_sub$batch_f <- factor(clin_sub$batch)
clin_sub$sex_f   <- factor(clin_sub$sex)
clin_sub$layer_f <- factor(sub("^.*_(L[0-9]+)$", "\\1",
                                as.character(clin_sub$sample_id)),
                            levels = paste0("L", 1:8))

# ---- CLR ----
cell_types <- setdiff(colnames(cell_prop_sub), "sample_id")
prop_mat <- as.matrix(cell_prop_sub[, cell_types, drop = FALSE])
prop_pc <- prop_mat + 0.5
gm <- exp(rowMeans(log(prop_pc)))
clr_mat <- log(prop_pc) - log(gm)
cell_types_use <- setdiff(cell_types, "Rod")
clr_df <- as.data.frame(clr_mat[, cell_types_use, drop = FALSE])
colnames(clr_df) <- paste0("clr_", cell_types_use)

sp_age <- splines::ns(clin_sub$age_v, df = 3)

design_base <- data.frame(
  grade_f = factor(clin_sub$grade),
  age_s1 = sp_age[, 1], age_s2 = sp_age[, 2], age_s3 = sp_age[, 3],
  sex_f = clin_sub$sex_f,
  rin_v = clin_sub$rin_v,
  pmi_v = clin_sub$pmi_v,
  layer_f = clin_sub$layer_f,
  batch_f = clin_sub$batch_f,
  stringsAsFactors = FALSE
)
design_adj <- cbind(design_base, clr_df)

# ============================================================
# ============================================================
build_poly_design <- function(design_df) {
  pm <- poly(as.integer(design_df$grade_f), 3)
  colnames(pm) <- c("p1", "p2", "p3")
  covar <- design_df[, setdiff(colnames(design_df), "grade_f"), drop = FALSE]
  model.matrix(~ 0 + p1 + p2 + p3 + .,
                data = data.frame(pm, covar, check.names = FALSE))
}

diag_design <- function(X, tag) {
  sv <- svd(X)$d
  tol <- max(dim(X)) * max(sv) * .Machine$double.eps
  num_rank <- sum(sv > tol)
  cond_num <- max(sv) / min(sv)
  log_info("[{tag}] ncol={ncol(X)} 数值秩={num_rank} cond={round(cond_num, 1)}")
  if (num_rank < ncol(X)) log_warn("[{tag}] 数值秩不足")
  if (cond_num > 1000) log_warn("[{tag}] cond > 1000，解释需谨慎")
  invisible(list(rank = num_rank, cond = cond_num))
}

X_orig <- build_poly_design(design_base)
X_adj  <- build_poly_design(design_adj)
diag_design(X_orig, "orig")
diag_design(X_adj,  "adj+CLR")

# ============================================================
# ============================================================
fit_and_extract <- function(X, score_mat) {
  fit <- eBayes(lmFit(score_mat, X))
  tt1 <- topTable(fit, coef = "p1", number = Inf, sort.by = "none")
  tt2 <- topTable(fit, coef = "p2", number = Inf, sort.by = "none")
  tt3 <- topTable(fit, coef = "p3", number = Inf, sort.by = "none")
  data.frame(
    pathway = rownames(tt1),
    p1 = tt1$P.Value, p1_FDR = tt1$adj.P.Val, p1_coef = tt1$logFC,
    p2 = tt2$P.Value, p2_FDR = tt2$adj.P.Val, p2_coef = tt2$logFC,
    p3 = tt3$P.Value, p3_FDR = tt3$adj.P.Val, p3_coef = tt3$logFC,
    stringsAsFactors = FALSE)
}

res_orig <- fit_and_extract(X_orig, score_mat)
res_adj  <- fit_and_extract(X_adj,  score_mat)

# ============================================================
# ============================================================
add_pattern <- function(res) {
  res$pattern_p <- with(res,
    ifelse(!is.na(p3) & p3 < 0.05 & p1 > 0.05 & p2 > 0.05, "p3_only",
    ifelse(!is.na(p3) & p3 < 0.05, "p3+other",
    ifelse(!is.na(p1) & p1 < 0.05, "p1_sig",
    ifelse(!is.na(p2) & p2 < 0.05, "p2_sig", "none")))))
  res$pattern_fdr <- with(res,
    ifelse(!is.na(p3_FDR) & p3_FDR < 0.05 & p1_FDR > 0.10 & p2_FDR > 0.10, "p3_only",
    ifelse(!is.na(p3_FDR) & p3_FDR < 0.05, "p3+other",
    ifelse(!is.na(p1_FDR) & p1_FDR < 0.05, "p1_sig",
    ifelse(!is.na(p2_FDR) & p2_FDR < 0.05, "p2_sig", "none")))))
  res
}
res_orig <- add_pattern(res_orig)
res_adj  <- add_pattern(res_adj)

# ============================================================
# ============================================================
merged <- merge(
  res_orig[, c("pathway","p3","p3_FDR","p3_coef","pattern_p","pattern_fdr")],
  res_adj[,  c("pathway","p3","p3_FDR","p3_coef","pattern_p","pattern_fdr")],
  by = "pathway", suffixes = c("_orig","_adj"))

for (pt in c("p","fdr")) {
  pc <- paste0("pattern_", pt)
  st <- paste0("still_", pt); lo <- paste0("lost_", pt); nw <- paste0("new_", pt)
  merged[[st]] <- merged[[paste0(pc,"_orig")]] == "p3_only" &
                  merged[[paste0(pc,"_adj")]]  == "p3_only"
  merged[[lo]] <- merged[[paste0(pc,"_orig")]] == "p3_only" &
                  merged[[paste0(pc,"_adj")]]  != "p3_only"
  merged[[nw]] <- merged[[paste0(pc,"_orig")]] != "p3_only" &
                  merged[[paste0(pc,"_adj")]]  == "p3_only"
}

merged$coef_sign_consistent <- sign(merged$p3_coef_orig) == sign(merged$p3_coef_adj)
merged$coef_sign_reversed <- !merged$coef_sign_consistent &
                              merged$pattern_p_orig == "p3_only" &
                              merged$pattern_p_adj  == "p3_only"

write.csv(merged, file.path(tab_dir, "pathway_poly3_cellprop_adj.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# ============================================================
# 6. summary
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("=== 04_2c v4 final 汇总 ===\n")
cat(strrep("=", 70), "\n", sep = "")

cat("\n【设计矩阵诊断】\n")
cat("orig: cond =", round(svd(X_orig)$d[1]/min(svd(X_orig)$d), 1), "\n")
cat("adj : cond =", round(svd(X_adj)$d[1]/min(svd(X_adj)$d), 1), "\n")

cat("\n【p3_only 对比：p 版（与 04_2 一致）】\n")
cat("  原始:", sum(merged$pattern_p_orig == "p3_only"), "\n")
cat("  校正后:", sum(merged$pattern_p_adj == "p3_only"), "\n")
cat("  保留:", sum(merged$still_p), "\n")
cat("  丢失:", sum(merged$lost_p), "\n")
cat("  新增:", sum(merged$new_p), "\n")
cat("  系数方向反转:", sum(merged$coef_sign_reversed), "\n")

cat("\n【p3_only 对比：FDR 版（附加）】\n")
cat("  原始:", sum(merged$pattern_fdr_orig == "p3_only"), "\n")
cat("  校正后:", sum(merged$pattern_fdr_adj == "p3_only"), "\n")
cat("  保留:", sum(merged$still_fdr), "\n")

cat("\n【判读】\n")
n_still_p <- sum(merged$still_p)
if (n_still_p >= 15) {
  cat("→ 相变不是细胞组成驱动的（cell-intrinsic）\n")
} else if (n_still_p >= 8) {
  cat("→ 部分由细胞组成驱动，需在 Discussion 说明\n")
} else {
  cat("→ 主要是细胞组成假象，需重写叙事\n")
}

cat("\n丢失的 p3_only 通路 (p 版):\n")
if (sum(merged$lost_p) > 0) {
  print(merged[merged$lost_p, c("pathway","p3_orig","p3_adj",
                                 "pattern_p_orig","pattern_p_adj")],
        row.names = FALSE, digits = 4)
} else cat("  无\n")

cat("\n新增的 p3_only 通路 (p 版):\n")
if (sum(merged$new_p) > 0) {
  print(merged[merged$new_p, c("pathway","p3_orig","p3_adj",
                                "pattern_p_orig","pattern_p_adj")],
        row.names = FALSE, digits = 4)
} else cat("  无\n")

log_info("=== 04_2c v4 final done ===")

