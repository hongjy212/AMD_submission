#!/usr/bin/env Rscript
# ============================================================
# 03_2_deconvolution.r (v3, self-contained)
#
#
# Output：results/03_2_deconvolution/
# ============================================================

rm(list = ls(all.names = TRUE)); gc()
setwd("E:/retina_tools/aging_rgc_axon/AMD")
set.seed(20260925)

suppressPackageStartupMessages({
  library(cellGeometry); library(SingleCellExperiment)
  library(zellkonverter); library(MASS); library(dplyr)
  library(tidyr); library(ggplot2); library(Matrix)
})

# ---- directory ----
DIR_ROOT   <- "E:/retina_tools/aging_rgc_axon/AMD"
DIR_RESULT <- file.path(DIR_ROOT, "results")
DIR_THIS   <- file.path(DIR_RESULT, "03_2_deconvolution")
for (d in c(DIR_THIS,
            file.path(DIR_THIS, "tables"),
            file.path(DIR_THIS, "figures"),
            file.path(DIR_THIS, "data"))) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

cat("=== 03_2_deconvolution.r start ===\n")

# ---- Parameters ----
CELLTYPE_LEVEL <- "majorclass"
SEED <- 20260925
COL_GRADE <- c("1" = "#4C72B0", "2" = "#55A868",
               "3" = "#C44E52", "4" = "#8172B3")

# ============================================================
# ============================================================
cat("--- 1. 载入 bulk ---\n")

expr_bulk <- readRDS(file.path(DIR_RESULT, "01_1_prepare_bulk_data",
                                "data", "expr_corrected.rds"))
gene_anno <- readRDS(file.path(DIR_RESULT, "01_1_prepare_bulk_data",
                                "data", "gene_annotation.rds"))
clinical  <- readRDS(file.path(DIR_RESULT, "01_1_prepare_bulk_data",
                                "data", "clinical.rds"))

cat("bulk 矩阵:", nrow(expr_bulk), "×", ncol(expr_bulk), "\n")
cat("clinical:", nrow(clinical), "样本\n")

# Ensembl → Symbol
if (mean(grepl("^ENSG", rownames(expr_bulk))) > 0.5) {
  cat("Ensembl 检测，映射到 Symbol...\n")
  map <- setNames(toupper(gene_anno$symbol), gene_anno$ensembl)
  new_ids <- map[rownames(expr_bulk)]
  new_ids[is.na(new_ids)] <- rownames(expr_bulk)[is.na(new_ids)]

  if (any(duplicated(new_ids))) {
    keep <- !duplicated(new_ids)
    expr_bulk <- expr_bulk[keep, ]
    new_ids <- new_ids[keep]
  }
  rownames(expr_bulk) <- new_ids
  cat("映射后基因数:", nrow(expr_bulk), "\n")
}

# voom logCPM → CPM
cat("将 voom logCPM 转回 CPM 空间...\n")
expr_cpm <- 2^expr_bulk - 0.5
expr_cpm[expr_cpm < 0] <- 0
cat("CPM 范围: [", round(min(expr_cpm), 2), ", ",
    round(max(expr_cpm), 2), "]\n")

keep <- apply(expr_cpm, 1, max) > 1
expr_cpm <- expr_cpm[keep, ]
cat("过滤后 CPM:", nrow(expr_cpm), "×", ncol(expr_cpm), "\n")

# samplesalign
common_samples <- intersect(colnames(expr_cpm), clinical$sample_id)
expr_cpm <- expr_cpm[, common_samples, drop = FALSE]
clinical <- clinical[match(common_samples, clinical$sample_id), ]
cat("对齐后:", ncol(expr_cpm), "样本\n")

clinical$amd_grade <- as.integer(clinical$grade)
clinical$age <- as.numeric(clinical$age_n)
clinical$sex <- as.factor(clinical$sex)

cat("AMD 等级分布:\n"); print(table(clinical$amd_grade))

clinical$batch_f <- factor(clinical$batch)
sid <- as.character(clinical$sample_id)
clinical$layer_f <- factor(sub("^.*_(L[0-9]+)$", "\\1", sid),
                            levels = paste0("L", 1:8))
cat("batch 分布:\n"); print(table(clinical$batch_f))
cat("layer 分布:\n"); print(table(clinical$layer_f))
cat("layer NA 数:", sum(is.na(clinical$layer_f)), "\n")

# ============================================================
# ============================================================
cat("--- 2. 载入单细胞 h5ad ---\n")

ref_file <- "E:/retina_tools/raw/GSE265801/HRCA_snRNA_subsampled_300.h5ad"
if (!file.exists(ref_file)) {
  alt <- list.files("E:/retina_tools/raw/GSE265801",
                    pattern = "\\.h5ad$", full.names = TRUE)
  if (length(alt) > 0) {
    ref_file <- alt[1]
    cat("自动找到:", ref_file, "\n")
  } else stop("缺 h5ad")
}

sce <- zellkonverter::readH5AD(ref_file, use_hdf5 = TRUE)
cat("单细胞:", ncol(sce), "cells ×", nrow(sce), "genes\n")
cat("meta 列:", paste(colnames(colData(sce)), collapse=', '), "\n")

celltype_col <- grep("celltype|cell_type|CellType",
                      colnames(colData(sce)), value = TRUE, ignore.case = TRUE)[1]
majorclass_col <- grep("majorclass|major_class|major",
                        colnames(colData(sce)), value = TRUE, ignore.case = TRUE)[1]

if (is.na(celltype_col) && is.na(majorclass_col)) {
  stop("未找到 celltype 或 majorclass 列")
}

cat("celltype 列:", celltype_col, "\n")
cat("majorclass 列:", majorclass_col, "\n")

if (!is.na(majorclass_col)) {
  cat("=== majorclass 分布 ===\n")
  print(table(colData(sce)[[majorclass_col]]))
}

# ============================================================
# ============================================================
cat("--- 3. 采样 ---\n")

coldata <- as.data.frame(colData(sce))

level_col <- if (!is.na(majorclass_col)) majorclass_col else celltype_col
cat("使用层级:", level_col, "\n")

labels_all <- coldata[[level_col]]
subtype_counts <- table(labels_all)
subtype_counts <- subtype_counts[!is.na(names(subtype_counts)) & names(subtype_counts) != ""]

max_per_subtype <- 2000
min_per_subtype <- 10
max_total <- 150000

set.seed(SEED)
selected <- c()
for (st in names(subtype_counts)) {
  cells <- which(labels_all == st)
  if (length(cells) > max_per_subtype) cells <- sample(cells, max_per_subtype)
  if (length(cells) >= min_per_subtype) selected <- c(selected, cells)
}
selected <- unique(selected)
if (length(selected) > max_total) selected <- sample(selected, max_total)
cat("采样细胞数:", format(length(selected), big.mark=','), "\n")

assay_names <- assayNames(sce)
orig_assay <- if ("counts" %in% assay_names) "counts" else assay_names[1]
counts_mat <- assay(sce, orig_assay)[, selected, drop = FALSE]

cell_labels <- coldata[[level_col]][selected]
valid <- !is.na(cell_labels) & cell_labels != ""
counts_mat <- counts_mat[, valid, drop = FALSE]
cell_labels <- cell_labels[valid]

cat("采样后:", nrow(counts_mat), "genes ×", ncol(counts_mat), "cells\n")

if (mean(grepl("^ENSG", rownames(counts_mat))) > 0.5) {
  map <- setNames(toupper(gene_anno$symbol), gene_anno$ensembl)
  new_ids <- map[rownames(counts_mat)]
  new_ids[is.na(new_ids)] <- rownames(counts_mat)[is.na(new_ids)]

  if (any(duplicated(new_ids))) {
    keep <- !duplicated(new_ids)
    counts_mat <- counts_mat[keep, ]
    new_ids <- new_ids[keep]
  }
  rownames(counts_mat) <- new_ids
}

rm(sce); gc()

common_genes <- intersect(rownames(counts_mat), rownames(expr_cpm))
counts_mat <- counts_mat[common_genes, , drop = FALSE]
expr_cpm <- expr_cpm[common_genes, , drop = FALSE]
cat("共同基因:", length(common_genes), "\n")

# ============================================================
# ============================================================
cat("--- 4. cellMarkers ---\n")

mk <- cellGeometry::cellMarkers(
  scdata = counts_mat,
  subclass = cell_labels,
  cellgroup = NULL,
  dual_mean = TRUE,
  min_cells = 10,
  cores = 1,
  bulkdata = expr_cpm,
  verbose = TRUE
)
cat("签名基因数:", nrow(mk$signature), "\n")
saveRDS(mk, file.path(DIR_THIS, "data", "cellMarkers.rds"))

cat("--- 5. 反卷积 ---\n")
t0 <- Sys.time()
fit <- tryCatch({
  cellGeometry::deconvolute(
    mk = mk,
    test = expr_cpm,
    logged_bulk = FALSE,
    count_space = TRUE,
    arith_mean = FALSE,
    use_filter = TRUE,
    weight_method = "equal",
    check_comp = TRUE,
    group_comp_amount = 0
  )
}, error = function(e) {
  cat("反卷积失败:", conditionMessage(e), "\n")
  NULL
})

if (is.null(fit)) stop("反卷积失败")
dt <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
cat("反卷积用时:", dt, "min\n")

prop_percent <- fit$subclass$percent
rownames(prop_percent) <- colnames(expr_cpm)
cat("比例矩阵:", paste(dim(prop_percent), collapse=' × '), "\n")

saveRDS(prop_percent, file.path(DIR_THIS, "data", "cell_proportions.rds"))
prop_df <- as.data.frame(prop_percent)
prop_df$sample_id <- rownames(prop_df)
write.csv(prop_df, file.path(DIR_THIS, "tables", "cell_proportions.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# ============================================================
# ============================================================
cat("--- 6. 关联 grade ---\n")

common_samples <- intersect(rownames(prop_percent), clinical$sample_id)
prop_percent <- prop_percent[common_samples, , drop = FALSE]
clinical_sub <- clinical[match(common_samples, clinical$sample_id), ]
clinical_sub$grade_f <- factor(clinical_sub$amd_grade, ordered = TRUE)

if (!"batch_f" %in% colnames(clinical_sub)) {
  clinical_sub$batch_f <- factor(clinical_sub$batch)
}
if (!"layer_f" %in% colnames(clinical_sub)) {
  sid2 <- as.character(clinical_sub$sample_id)
  clinical_sub$layer_f <- factor(sub("^.*_(L[0-9]+)$", "\\1", sid2),
                                  levels = paste0("L", 1:8))
}

ct_cols <- colnames(prop_percent)
cat("celltype 数:", length(ct_cols), "\n")
cat("celltype:", paste(ct_cols, collapse=', '), "\n")

polr_res <- data.frame()
for (ct in ct_cols) {
  y <- prop_percent[, ct]
  df <- data.frame(y = y, grade_f = clinical_sub$grade_f,
                   age = clinical_sub$age, sex = clinical_sub$sex,
                   rin = clinical_sub$rin_n, pmi = clinical_sub$pmi_n,
                   batch = clinical_sub$batch_f,
                   layer = clinical_sub$layer_f)
  df <- df[complete.cases(df), ]

  fit_p <- tryCatch(
    MASS::polr(grade_f ~ y + age + sex + rin + pmi + batch + layer,
               data = df, Hess = TRUE),
    error = function(e) NULL)
  if (is.null(fit_p)) next

  co <- summary(fit_p)$coefficients
  if (!"y" %in% rownames(co)) next

  est <- co["y", "Value"]
  tval <- co["y", "t value"]
  polr_res <- rbind(polr_res, data.frame(
    celltype = ct,
    Estimate = round(est, 4),
    OR = round(exp(est), 4),
    t = round(tval, 3),
    p = 2 * pnorm(-abs(tval)),
    stringsAsFactors = FALSE))
}

if (nrow(polr_res) > 0) {
  polr_res$FDR <- p.adjust(polr_res$p, "BH")
  polr_res <- polr_res[order(polr_res$p), ]
  write.csv(polr_res,
            file.path(DIR_THIS, "tables", "celltype_vs_grade_polr.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
  cat("\n=== polr (含 batch+layer) ===\n")
  print(polr_res, row.names = FALSE, digits = 4)
}

lm_res <- data.frame()
for (ct in ct_cols) {
  y <- prop_percent[, ct]
  df <- data.frame(y = y, grade = clinical_sub$amd_grade,
                   age = clinical_sub$age, sex = clinical_sub$sex,
                   rin = clinical_sub$rin_n, pmi = clinical_sub$pmi_n,
                   batch = clinical_sub$batch_f,
                   layer = clinical_sub$layer_f)
  df <- df[complete.cases(df), ]

  fit_l <- tryCatch(lm(y ~ grade + age + sex + rin + pmi + batch + layer,
                        data = df),
                     error = function(e) NULL)
  if (is.null(fit_l)) next

  co <- summary(fit_l)$coefficients
  lm_res <- rbind(lm_res, data.frame(
    celltype = ct,
    beta_grade = round(co["grade", "Estimate"], 5),
    t = round(co["grade", "t value"], 3),
    p = co["grade", "Pr(>|t|)"],
    stringsAsFactors = FALSE))
}

if (nrow(lm_res) > 0) {
  lm_res$FDR <- p.adjust(lm_res$p, "BH")
  lm_res <- lm_res[order(lm_res$p), ]
  write.csv(lm_res,
            file.path(DIR_THIS, "tables", "celltype_vs_grade_lm.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
  cat("\n=== lm (含 batch+layer) ===\n")
  print(lm_res, row.names = FALSE, digits = 4)
}

# ============================================================
# 7. figure
# ============================================================
cat("--- 7. 图 ---\n")

prop_long <- as.data.frame(prop_percent)
prop_long$sample_id <- rownames(prop_long)
prop_long <- pivot_longer(prop_long, cols = all_of(ct_cols),
                            names_to = "celltype", values_to = "proportion")
prop_long$grade <- clinical_sub$amd_grade[match(prop_long$sample_id,
                                                  clinical_sub$sample_id)]
prop_long$grade <- factor(prop_long$grade)

p1 <- ggplot(prop_long, aes(x = grade, y = proportion, fill = grade)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  facet_wrap(~ celltype, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = COL_GRADE) +
  theme_bw(base_size = 11) +
  labs(x = "AMD grade", y = "Estimated proportion (%)",
       fill = "Grade", title = "Cell proportions by AMD grade")
ggsave(file.path(DIR_THIS, "figures", "proportions_by_grade.pdf"),
       p1, width = 11, height = 8)
ggsave(file.path(DIR_THIS, "figures", "proportions_by_grade.png"),
       p1, width = 11, height = 8, dpi = 300)

if (nrow(polr_res) > 0) {
  plot_df <- polr_res
  plot_df$celltype <- factor(plot_df$celltype,
                              levels = plot_df$celltype[order(plot_df$OR)])
  plot_df$log_OR <- log(plot_df$OR)

  p2 <- ggplot(plot_df, aes(x = celltype, y = log_OR, fill = p < 0.05)) +
    geom_col() + coord_flip() + theme_bw(base_size = 11) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = c("TRUE" = "#C44E52", "FALSE" = "#B0B0B0"),
                      labels = c("TRUE" = "p<0.05", "FALSE" = "NS")) +
    labs(x = "", y = "log(OR)", fill = "",
         title = "polr: AMD grade ~ cell proportion (+ batch + layer)")
  ggsave(file.path(DIR_THIS, "figures", "polr_logOR.pdf"),
         p2, width = 7, height = 5)
  ggsave(file.path(DIR_THIS, "figures", "polr_logOR.png"),
         p2, width = 7, height = 5, dpi = 300)
}

# ============================================================
# 8. summary
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("=== 03_2 反卷积 汇总 (含 batch + layer) ===\n")
cat(strrep("=", 70), "\n", sep = "")
cat("样本数:", nrow(prop_percent), "\n")
cat("celltype 数:", ncol(prop_percent), "\n")
cat("celltype:", paste(colnames(prop_percent), collapse=', '), "\n\n")

if (nrow(polr_res) > 0) {
  cat("--- polr ---\n")
  print(polr_res, row.names = FALSE, digits = 4)
}
if (nrow(lm_res) > 0) {
  cat("\n--- lm ---\n")
  print(lm_res, row.names = FALSE, digits = 4)
}

cat("\n=== 03_2_deconvolution.r done ===\n")

