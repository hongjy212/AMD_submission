#!/usr/bin/env Rscript
# ============================================================
# 01_1b_batch_enhance.r
#
# Input：
#   results/01_1_prepare_bulk_data/data/{expr_corrected,clinical,gene_annotation}.rds
#   E:/retina_tools/raw/GSE115828/NIHMS1518564-supplement-2.xlsx
#   E:/retina_tools/raw/GSE115828/GSE115828_series_matrix.txt.gz
#   data/clinical_wgcna.rds
#   tables/batch_audit.csv
#   tables/r2_diagnostic.csv
#   figures/PCA_before_after.pdf/.png
# ============================================================

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(data.table); library(readxl); library(GEOquery)
  library(limma); library(splines); library(ggplot2); library(dplyr)
})

DIR_THIS <- setup_output_dir("01_1b_batch_enhance")
log_info("=== 01_1b_batch_enhance.r start ===")

# ============================================================
# ============================================================
log_info("--- 1. 载入 01_1 输出 ---")

expr <- read_upstream("01_1_prepare_bulk_data", "data", "expr_corrected.rds")
clin <- read_upstream("01_1_prepare_bulk_data", "data", "clinical.rds")
anno <- read_upstream("01_1_prepare_bulk_data", "data", "gene_annotation.rds")

common <- intersect(colnames(expr), rownames(clin))
expr <- expr[, common, drop = FALSE]
clin <- clin[common, , drop = FALSE]

log_info("expr: {nrow(expr)} × {ncol(expr)}")
log_info("clinical: {nrow(clin)} 样本")

# ============================================================
# ============================================================
log_info("--- 2. 解析 sample_id ---")

sid <- as.character(clin$sample_id)
clin$batch     <- sub("^(R\\d+)-.*", "\\1", sid)
clin$pf_id     <- sub("^R\\d+-([0-9]+pf)_.*", "\\1", sid)
clin$donor_num <- sub("^R\\d+-[0-9]+pf_([0-9]+)-.*", "\\1", sid)
clin$layer     <- sub("^.*_(L[0-9]+)$", "\\1", sid)

log_info("batch 唯一值: {paste(unique(clin$batch), collapse=', ')}")
log_info("donor_num 样例: {paste(head(clin$donor_num, 5), collapse=', ')}")
log_info("layer 唯一值: {paste(sort(unique(clin$layer)), collapse=', ')}")

# ============================================================
# ============================================================
log_info("--- 3. 读补充材料 ---")

supp_file <- "E:/retina_tools/raw/GSE115828/NIHMS1518564-supplement-2.xlsx"
if (!file.exists(supp_file)) {
  stop("补充材料不存在: ", supp_file)
}

supp_extract_donor <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("[[:space:]]", "", x)
  sub("_[0-9]+$", "", x)
}

supp_raw <- readxl::read_excel(supp_file, sheet = 1, col_names = FALSE)
header_row <- which(apply(supp_raw, 1, function(r) {
  any(grepl("^r_id$", trimws(as.character(r)), ignore.case = TRUE))
}))[1]
if (is.na(header_row)) header_row <- 5
log_info("补充材料表头行: {header_row}")

supp_meta <- readxl::read_excel(supp_file, sheet = 1, skip = header_row - 1)
supp_meta <- as.data.frame(supp_meta, stringsAsFactors = FALSE)
log_info("补充材料维度: {nrow(supp_meta)} × {ncol(supp_meta)}")
log_info("列名: {paste(colnames(supp_meta), collapse=', ')}")

supp_meta$donor_num_key <- supp_extract_donor(supp_meta$r_id)
supp_meta <- supp_meta[!duplicated(supp_meta$donor_num_key), , drop = FALSE]

# ============================================================
# ============================================================
log_info("--- 4. 合并补充材料 ---")

idx <- match(clin$donor_num, supp_meta$donor_num_key)
log_info("匹配样本数: {sum(!is.na(idx))} / {length(idx)}")

assign_from_supp <- function(col_out, col_in) {
  if (col_in %in% colnames(supp_meta)) {
    clin[[col_out]] <<- supp_meta[[col_in]][idx]
    log_info("  {col_out} <- {col_in}")
  } else {
    log_warn("  补充材料缺列: {col_in}，{col_out} 保持 NA")
    clin[[col_out]] <<- NA
  }
}

assign_from_supp("rna_isolation_batch", "rna_isolation_batch")
assign_from_supp("lib_prep_batch",      "lib_prep_batch")
assign_from_supp("library_prepper",     "library_prepper")
assign_from_supp("seq_date",            "library_sequenced_date")
assign_from_supp("rna_iso_date",        "rna_isolation_date")
assign_from_supp("rnaseq_reads",        "rnaseq_reads")
assign_from_supp("smoking",             "smoking_history")
assign_from_supp("hypertension",        "hypertension")
assign_from_supp("heart_disease",       "heart_disease")
assign_from_supp("cholesterol",         "cholesterol")
assign_from_supp("death_category",      "death_category")
assign_from_supp("donor",               "donor")

n_matched <- sum(!is.na(clin$rna_isolation_batch))
log_info("批次信息匹配: {n_matched} / {nrow(clin)}")

save_rds(clin, DIR_THIS, "clinical_453_full.rds")

# ============================================================
# ============================================================
log_info("--- 5. 剔除混杂样本 ---")

clin$batch_f <- factor(clin$batch)

batch_check_vars <- c("rna_isolation_batch", "lib_prep_batch", "library_prepper")
has_batch <- rep(TRUE, nrow(clin))
for (vn in batch_check_vars) {
  if (vn %in% colnames(clin)) {
    has_batch <- has_batch & !is.na(clin[[vn]])
  }
}
log_info("批次信息缺失: {sum(!has_batch)} 个")

# 5b. library_prepper == 3
prep_val <- suppressWarnings(as.numeric(as.character(clin$library_prepper)))
is_drop_prepper <- !is.na(prep_val) & prep_val == 3
log_info("library_prepper==3: {sum(is_drop_prepper)} 个")

keep_mask <- has_batch & !is_drop_prepper
log_info("合计剔除: {sum(!keep_mask)} 个 → 保留 {sum(keep_mask)} 个")

drop_tbl <- data.frame(
  sample_id = clin$sample_id[!keep_mask],
  grade = clin$grade[!keep_mask],
  reason = ifelse(!has_batch[!keep_mask], "missing_batch_info",
                   "library_prepper==3"),
  stringsAsFactors = FALSE)
save_table(drop_tbl, DIR_THIS, "dropped_samples.csv")

expr_sub <- expr[, keep_mask, drop = FALSE]
clin_sub <- clin[keep_mask, , drop = FALSE]

log_info("剔除后样本数: {nrow(clin_sub)}")
print(table(clin_sub$grade))
print(table(clin_sub$batch_f, clin_sub$grade))

# ============================================================
# ============================================================
log_info("--- 6. 批次审计 ---")

batch_audit <- data.frame()
for (vn in c("batch", "rna_isolation_batch", "lib_prep_batch",
              "library_prepper", "seq_date")) {
  if (vn %in% colnames(clin_sub)) {
    tab <- table(clin_sub[[vn]], useNA = "ifany")
    batch_audit <- rbind(batch_audit, data.frame(
      variable = vn, level = names(tab), n = as.integer(tab),
      stringsAsFactors = FALSE))
  }
}
save_table(batch_audit, DIR_THIS, "batch_audit.csv")

# rna_isolation_batch × grade
cat("\n=== rna_isolation_batch × grade ===\n")
print(table(clin_sub$rna_isolation_batch, clin_sub$grade, useNA = "ifany"))

cat("\n=== lib_prep_batch × grade ===\n")
print(table(clin_sub$lib_prep_batch, clin_sub$grade, useNA = "ifany"))

# ============================================================
# 7. recompute PC1-PC5
# ============================================================
log_info("--- 7. 重算 PC1-PC5 ---")

pca_new <- prcomp(t(expr_sub), center = TRUE, scale. = FALSE)
for (i in 1:5) {
  clin_sub[[paste0("PC", i)]] <- pca_new$x[, i]
}
log_info("PC1 方差解释: {round(summary(pca_new)$importance[2,1]*100, 1)}%")

# ============================================================
# ============================================================
log_info("--- 8. 构造 design_full ---")

clin_sub$grade_int <- as.integer(clin_sub$grade)
clin_sub$age_v     <- as.numeric(clin_sub$age_n)
clin_sub$rin_v     <- as.numeric(clin_sub$rin_n)
clin_sub$pmi_v     <- as.numeric(clin_sub$pmi_n)
clin_sub$sex_f     <- factor(clin_sub$sex)
clin_sub$layer_f   <- factor(clin_sub$layer, levels = paste0("L", 1:8))

design_vars <- c(
  "grade_factor",                         # 保护
  "age_v", "sex_f", "rin_v", "pmi_v",     # 保护
  "batch_f",                              # 扣
  "lib_prep_batch",                       # 扣
  "rna_isolation_batch",                  # 扣
  paste0("PC", 1:5)                       # 扣
)

if ("rnaseq_reads" %in% colnames(clin_sub)) {
  reads_num <- suppressWarnings(as.numeric(clin_sub$rnaseq_reads))
  if (sum(is.finite(reads_num)) > nrow(clin_sub) * 0.9) {
    clin_sub$log_reads <- log10(reads_num + 1)
    design_vars <- c(design_vars, "log_reads")
    log_info("加入 log_reads")
  }
}

clin_sub$lib_prep_batch      <- factor(clin_sub$lib_prep_batch)
clin_sub$rna_isolation_batch <- factor(clin_sub$rna_isolation_batch)
grade_factor <- factor(clin_sub$grade_int, levels = 1:4)

# design matrix
design_df <- clin_sub
design_df$grade_factor <- grade_factor
design_full <- model.matrix(
  as.formula(paste("~ 0 +", paste(design_vars, collapse = " + "))),
  data = design_df
)
colnames(design_full)[grep("^grade_factor", colnames(design_full))] <- paste0("G", 1:4)

log_info("design_full: {nrow(design_full)} × {ncol(design_full)}")
log_info("列名: {paste(colnames(design_full), collapse=', ')}")

rank_val <- qr(design_full)$rank
if (rank_val < ncol(design_full)) {
  log_warn("design_full 秩亏: rank={rank_val}, ncol={ncol(design_full)}")
  svd_d <- svd(design_full)
  log_warn("最小奇异值: {round(min(svd_d$d), 4)}")
  log_warn("尝试去掉 rna_isolation_batch...")
  design_vars <- setdiff(design_vars, "rna_isolation_batch")
  design_full <- model.matrix(
    as.formula(paste("~ 0 +", paste(design_vars, collapse = " + "))),
    data = design_df
  )
  colnames(design_full)[grep("^grade_factor", colnames(design_full))] <- paste0("G", 1:4)
  rank_val <- qr(design_full)$rank
  log_info("修正后 design_full rank = {rank_val}")
}
log_info("design_full rank = {rank_val}")

# ============================================================
# ============================================================
log_info("--- 9. lmFit 扣批次 ---")

fit <- lmFit(expr_sub, design_full)

coef_remove <- grep("^batch_f|^lib_prep_batch|^rna_isolation_batch|^PC[0-9]+$",
                     colnames(design_full), value = TRUE)
log_info("扣掉的系数: {paste(coef_remove, collapse=', ')}")

if (length(coef_remove) == 0) {
  log_warn("没有可扣的系数，expr_wgcna = expr_sub")
  expr_wgcna <- expr_sub
} else {
  effect <- fit$coefficients[, coef_remove, drop = FALSE] %*%
            t(design_full[, coef_remove, drop = FALSE])
  expr_wgcna <- expr_sub - effect
}
log_info("expr_wgcna: {nrow(expr_wgcna)} × {ncol(expr_wgcna)}")

# ============================================================
# 10. R² diagnosis
# ============================================================
log_info("--- 10. R² 诊断 ---")

pca_before <- prcomp(t(expr_sub), center = TRUE, scale. = FALSE)
pca_after  <- prcomp(t(expr_wgcna), center = TRUE, scale. = FALSE)

r2_diag <- data.frame()

for (vn in c("batch_f", "lib_prep_batch", "rna_isolation_batch",
              "library_prepper", "layer_f", "grade_factor")) {
  if (vn %in% colnames(clin_sub)) {
    x <- factor(clin_sub[[vn]])
    if (nlevels(x) > 1 && sum(!is.na(x)) > 3) {
      r2b <- summary(aov(pca_before$x[, 1] ~ x))[[1]]
      r2a <- summary(aov(pca_after$x[, 1]  ~ x))[[1]]
      r2_diag <- rbind(r2_diag, data.frame(
        variable = vn,
        type = "categorical",
        R2_before = round(r2b[1, "Sum Sq"] / sum(r2b[, "Sum Sq"]), 3),
        R2_after  = round(r2a[1, "Sum Sq"] / sum(r2a[, "Sum Sq"]), 3),
        stringsAsFactors = FALSE))
    }
  }
}

for (vn in c("age_v", "rin_v", "pmi_v", "log_reads")) {
  if (vn %in% colnames(clin_sub)) {
    x <- as.numeric(clin_sub[[vn]])
    if (sum(is.finite(x)) > 3) {
      r_b <- cor(pca_before$x[, 1], x, use = "complete.obs")
      r_a <- cor(pca_after$x[, 1],  x, use = "complete.obs")
      r2_diag <- rbind(r2_diag, data.frame(
        variable = vn,
        type = "continuous",
        R2_before = round(r_b^2, 3),
        R2_after  = round(r_a^2, 3),
        stringsAsFactors = FALSE))
    }
  }
}

r2_diag$delta <- r2_diag$R2_after - r2_diag$R2_before
save_table(r2_diag, DIR_THIS, "r2_diagnostic.csv")

cat("\n=== PC1 校正前后 R² ===\n")
print(r2_diag, row.names = FALSE)

# ============================================================
# ============================================================
log_info("--- 11. PCA 对比图 ---")

pca_df <- rbind(
  data.frame(PC1 = pca_before$x[, 1], PC2 = pca_before$x[, 2],
              batch = clin_sub$batch_f, stage = "Before"),
  data.frame(PC1 = pca_after$x[, 1],  PC2 = pca_after$x[, 2],
              batch = clin_sub$batch_f, stage = "After")
)

p1 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = batch)) +
  geom_point(size = 1.2, alpha = 0.7) +
  facet_wrap(~ stage, scales = "free") +
  theme_bw(base_size = 11) +
  labs(title = "PCA: before vs after batch correction",
       color = "batch")
save_fig(p1, DIR_THIS, "PCA_before_after_batch", width = 10, height = 5)

pca_df$rna_iso_batch <- rep(clin_sub$rna_isolation_batch, 2)
pca_df$lib_prep_batch <- rep(clin_sub$lib_prep_batch, 2)

p2 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = lib_prep_batch)) +
  geom_point(size = 1.2, alpha = 0.7) +
  facet_wrap(~ stage, scales = "free") +
  theme_bw(base_size = 11) +
  labs(title = "PCA: lib_prep_batch",
       color = "lib_prep_batch")
save_fig(p2, DIR_THIS, "PCA_before_after_lib_batch", width = 10, height = 5)

p3 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = rna_iso_batch)) +
  geom_point(size = 1.2, alpha = 0.7) +
  facet_wrap(~ stage, scales = "free") +
  theme_bw(base_size = 11) +
  labs(title = "PCA: rna_isolation_batch",
       color = "rna_iso_batch")
save_fig(p3, DIR_THIS, "PCA_before_after_iso_batch", width = 10, height = 5)

# grade
pca_df$grade <- rep(factor(clin_sub$grade_int), 2)
p4 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = grade)) +
  geom_point(size = 1.2, alpha = 0.7) +
  facet_wrap(~ stage, scales = "free") +
  theme_bw(base_size = 11) +
  labs(title = "PCA: grade",
       color = "grade")
save_fig(p4, DIR_THIS, "PCA_before_after_grade", width = 10, height = 5)

# ============================================================
# 12. save
# ============================================================
log_info("--- 12. 保存 ---")

save_rds(expr_wgcna, DIR_THIS, "expr_wgcna.rds")
save_rds(clin_sub,   DIR_THIS, "clinical_wgcna.rds")

# ============================================================
# 13. summary
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("=== 01_1b 汇总 ===\n")
cat(strrep("=", 70), "\n", sep = "")

cat("\n【Q1】样本过滤\n")
cat("  before:", 453, "\n")
cat("  after: ", nrow(clin_sub), "\n")

cat("\n【Q2】批次列匹配\n")
cat("  rna_isolation_batch 非 NA:", sum(!is.na(clin_sub$rna_isolation_batch)), "\n")
cat("  lib_prep_batch 非 NA:      ", sum(!is.na(clin_sub$lib_prep_batch)), "\n")

cat("\n【Q3】PC1 R² 变化\n")
print(r2_diag, row.names = FALSE)

cat("\n【Q4】关键判据\n")
batch_r2_after <- r2_diag$R2_after[r2_diag$variable == "batch_f"]
if (length(batch_r2_after) > 0) {
  cat(sprintf("  batch_f R² after: %.3f  →  %s\n",
              batch_r2_after,
              ifelse(batch_r2_after < 0.02, "✓ 校正成功", "⚠ 仍残留")))
}
iso_r2_after <- r2_diag$R2_after[r2_diag$variable == "rna_isolation_batch"]
if (length(iso_r2_after) > 0) {
  cat(sprintf("  rna_isolation_batch R² after: %.3f\n", iso_r2_after))
}

cat("\n--- 输出位置 ---\n")
cat("目录:", DIR_THIS, "\n")

log_info("=== 01_1b_batch_enhance.r done ===")
