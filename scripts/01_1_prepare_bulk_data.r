# ============================================================
# 01_1_prepare_bulk_data.r
#   data/expr_count.rds
#   data/voom_full.rds
#   data/expr_corrected.rds
#   data/clinical.rds
#   data/gene_annotation.rds
#   tables/*.csv
#   figures/*.png / *.pdf
# ============================================================

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

DIR_THIS <- setup_output_dir("01_1_prepare_bulk_data")
log_info("=== 01_1_prepare_bulk_data.r start ===")

DIR_GSE <- file.path(DIR_RAW, "GSE115828")
stopifnot(dir.exists(DIR_GSE))

log_info("--- 1. Normalized_CPM 白名单 ---")
cpm_file <- file.path(DIR_GSE, "GSE115828_Normalized_CPM.txt.gz")
stopifnot(file.exists(cpm_file))

con <- gzfile(cpm_file, "r")
header_cpm <- readLines(con, n = 1)
close(con)
cols_cpm <- strsplit(header_cpm, "\t")[[1]]
samples_keep <- cols_cpm[-1]
log_info("Normalized_CPM 样本数: {length(samples_keep)}")

con <- gzfile(cpm_file, "r")
lines_cpm <- readLines(con, warn = FALSE)
close(con)
genes_cpm <- vapply(strsplit(lines_cpm[-1], "\t"),
                    function(x) x[1], character(1))
log_info("Normalized_CPM 基因数: {length(genes_cpm)}")

# ---------- 2. raw count ----------
log_info("--- 2. raw count ---")
count_file <- file.path(DIR_GSE, "GSE115828_RSEM_gene_counts.tsv.gz")
stopifnot(file.exists(count_file))

counts_all <- read.delim(gzfile(count_file), row.names = 1,
                          check.names = FALSE)
log_info("raw count: {nrow(counts_all)} x {ncol(counts_all)}")

samples_use <- intersect(colnames(counts_all), samples_keep)
log_info("样本交集: {length(samples_use)} (raw={ncol(counts_all)}, CPM={length(samples_keep)})")
if (length(samples_use) == 0) stop("样本交集为空")
counts_all <- counts_all[, samples_use, drop = FALSE]

# ---------- 3. series matrix → clinical ----------
log_info("--- 3. series matrix ---")
sm_file <- file.path(DIR_GSE, "GSE115828_series_matrix.txt.gz")
sm_lines <- readLines(gzfile(sm_file), warn = FALSE)
sm_meta  <- sm_lines[grepl("^!Sample", sm_lines)]

extract_field <- function(pattern) {
  l <- sm_meta[grepl(pattern, sm_meta)]
  if (length(l) == 0) return(character(0))
  v <- strsplit(l[1], "\t")[[1]][-1]
  gsub('^"|"$', '', v)
}

title <- extract_field("!Sample_title")
geo   <- extract_field("!Sample_geo_accession")
char_lines <- sm_meta[grepl("^!Sample_characteristics_ch1", sm_meta)]

clinical <- data.frame(sample_id = title, geo = geo, stringsAsFactors = FALSE)
for (line in char_lines) {
  v <- strsplit(line, "\t")[[1]][-1]
  v <- gsub('^"|"$', '', v)
  key <- sub(":.*$", "", v[1])
  val <- sub("^[^:]+:\\s*", "", v)
  clinical[[key]] <- val
}

clinical$grade <- suppressWarnings(as.integer(clinical$mgs_level))
clinical$age_n <- suppressWarnings(as.numeric(clinical$age))
clinical$rin_n <- suppressWarnings(as.numeric(clinical$rin))
clinical$pmi_n <- suppressWarnings(as.numeric(clinical$pmi))
clinical$sex   <- ifelse(clinical$Sex == "M", "Male",
                  ifelse(clinical$Sex == "F", "Female", NA))
clinical$batch <- sub("^([A-Za-z0-9]+)-.*$", "\\1", clinical$sample_id)

log_info("clinical 维度: {nrow(clinical)} x {ncol(clinical)}")
log_info("grade 分布:")
print(table(clinical$grade, useNA = "ifany"))

# ---------- 4. align ----------
log_info("--- 4. 对齐 ---")
common <- intersect(colnames(counts_all), clinical$sample_id)
log_info("count ∩ clinical: {length(common)}")
counts_all <- counts_all[, common, drop = FALSE]
clinical   <- clinical[match(common, clinical$sample_id), , drop = FALSE]
rownames(clinical) <- clinical$sample_id

keep <- !is.na(clinical$grade) & !is.na(clinical$age_n)
log_info("过滤前 {nrow(clinical)}, 过滤后 {sum(keep)}")
clinical   <- clinical[keep, , drop = FALSE]
counts_all <- counts_all[, keep, drop = FALSE]

log_info("--- 5. 基因白名单 ---")
rn  <- rownames(counts_all)
rn0 <- sub("\\..*$", "", rn)
gk0 <- sub("\\..*$", "", genes_cpm)
keep_gene <- rn0 %in% gk0
log_info("基因交集: {sum(keep_gene)} (raw={length(rn)}, CPM={length(genes_cpm)})")
counts_all <- counts_all[keep_gene, , drop = FALSE]

# ---------- 6. X/Y exclude ----------
log_info("--- 6. X/Y 排除 ---")
gene_anno <- data.frame(ensembl = rownames(counts_all), stringsAsFactors = FALSE)
gene_anno$ensembl_no_ver <- sub("\\..*$", "", gene_anno$ensembl)
gene_anno$chr <- suppressMessages(
  mapIds(org.Hs.eg.db,
         keys    = gene_anno$ensembl_no_ver,
         column  = "CHR",
         keytype = "ENSEMBL",
         multiVals = "first")
)
log_info("X: {sum(gene_anno$chr == 'X', na.rm=TRUE)}, Y: {sum(gene_anno$chr == 'Y', na.rm=TRUE)}")
keep_xy <- !(gene_anno$chr %in% c("X", "Y"))
counts_all <- counts_all[keep_xy, , drop = FALSE]
gene_anno  <- gene_anno[keep_xy, , drop = FALSE]
log_info("X/Y 排除后: {nrow(counts_all)} 基因")

# ---------- 7. voom ----------
log_info("--- 7. voom ---")
dge <- DGEList(counts = as.matrix(counts_all))
dge <- calcNormFactors(dge, method = "TMM")

clinical$grade_factor <- factor(clinical$grade)
design_simple <- model.matrix(~ grade_factor, data = clinical)

voom_obj <- voom(dge, design = design_simple, plot = FALSE)
voom_full <- voom_obj$E
log_info("voom_full: {nrow(voom_full)} x {ncol(voom_full)}")

keep_expr <- rowSums(voom_full > 0) >= 0.5 * ncol(voom_full)
voom_full <- voom_full[keep_expr, , drop = FALSE]
gene_anno <- gene_anno[keep_expr, , drop = FALSE]
counts_all <- counts_all[keep_expr, , drop = FALSE]
log_info("表达过滤后: {nrow(voom_full)} 基因")

# ---------- 8. PCA ----------
log_info("--- 8. PCA ---")
gene_var  <- apply(voom_full, 1, var)
top_genes <- names(sort(gene_var, decreasing = TRUE))[1:min(5000, nrow(voom_full))]
pca <- prcomp(t(voom_full[top_genes, , drop = FALSE]), scale. = TRUE)

var_pct <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 2)
log_info("PC1-10 方差: {paste(var_pct[1:10], collapse = ', ')}")

pc_df <- as.data.frame(pca$x[, 1:5])
colnames(pc_df) <- paste0("PC", 1:5)
pc_df$sample_id <- rownames(pc_df)
clinical <- merge(clinical, pc_df, by = "sample_id", all.x = TRUE)
rownames(clinical) <- clinical$sample_id

pc_diag <- data.frame()
for (i in 1:5) {
  r_b <- cor(as.numeric(factor(clinical$batch)),
             clinical[[paste0("PC", i)]], use = "complete.obs")
  r_g <- cor(clinical$grade,
             clinical[[paste0("PC", i)]], use = "complete.obs")
  pc_diag <- rbind(pc_diag, data.frame(
    PC = paste0("PC", i),
    r_batch = round(r_b, 3),
    r_grade = round(r_g, 3),
    var_pct = var_pct[i],
    stringsAsFactors = FALSE))
  log_info("PC{i}: r(batch)={round(r_b,3)}, r(grade)={round(r_g,3)}")
}
save_table(pc_diag, DIR_THIS, "pc_covariate.csv")

log_info("--- 9. 批次校正 ---")
clinical$sex_factor <- factor(clinical$sex)

protect_vars <- c("grade_factor", "age_n", "sex_factor", "rin_n", "pmi_n")
correct_vars <- c("batch", paste0("PC", 1:5))

fml <- as.formula(paste("~ 0 +",
                        paste(c(protect_vars, correct_vars), collapse = " + ")))
design_full <- model.matrix(fml, data = clinical)
log_info("design: {nrow(design_full)} x {ncol(design_full)}")

fit <- lmFit(voom_full, design_full)
coef_remove <- grep("^batch|^PC[0-9]", colnames(design_full), value = TRUE)
log_info("扣掉系数: {paste(coef_remove, collapse = ', ')}")

expr_corrected <- voom_full
if (length(coef_remove) > 0) {
  beta_remove <- fit$coefficients[, coef_remove, drop = FALSE]
  X_remove    <- design_full[, coef_remove, drop = FALSE]
  expr_corrected <- voom_full - beta_remove %*% t(X_remove)
}
log_info("校正后矩阵: {nrow(expr_corrected)} x {ncol(expr_corrected)}")

# ---------- 10. Symbol annotation ----------
log_info("--- 10. Symbol 注释 ---")
gene_anno$symbol <- suppressMessages(
  mapIds(org.Hs.eg.db,
         keys    = gene_anno$ensembl_no_ver,
         column  = "SYMBOL",
         keytype = "ENSEMBL",
         multiVals = "first")
)
log_info("Symbol 注释: {sum(!is.na(gene_anno$symbol))} / {nrow(gene_anno)}")

log_info("--- 11. 保存 ---")
save_rds(counts_all,     DIR_THIS, "expr_count.rds")
save_rds(voom_full,      DIR_THIS, "voom_full.rds")
save_rds(expr_corrected, DIR_THIS, "expr_corrected.rds")
save_rds(clinical,       DIR_THIS, "clinical.rds")
save_rds(gene_anno,      DIR_THIS, "gene_annotation.rds")

# ---------- 12. QC summary ----------
qc_summary <- data.frame(
  item = c("n_samples", "n_genes_raw", "n_genes_after_XY",
            "n_genes_after_expr", "grade_1", "grade_2",
            "grade_3", "grade_4",
            "batch_R42015", "batch_R42016",
            "age_median", "rin_median", "pmi_median"),
  value = c(nrow(clinical), nrow(counts_all), nrow(gene_anno),
            nrow(voom_full),
            sum(clinical$grade == 1), sum(clinical$grade == 2),
            sum(clinical$grade == 3), sum(clinical$grade == 4),
            sum(clinical$batch == "R42015"), sum(clinical$batch == "R42016"),
            median(clinical$age_n, na.rm = TRUE),
            median(clinical$rin_n, na.rm = TRUE),
            median(clinical$pmi_n, na.rm = TRUE)),
  stringsAsFactors = FALSE)
save_table(qc_summary, DIR_THIS, "qc_summary.csv")

# ---------- 13. figure ----------
log_info("--- 13. 图 ---")

pca_before <- pca$x[, 1:2]
pca_after_obj <- prcomp(t(expr_corrected[top_genes, , drop = FALSE]), scale. = TRUE)
pca_after <- pca_after_obj$x[, 1:2]

df_before <- data.frame(PC1 = pca_before[,1], PC2 = pca_before[,2],
                         batch = factor(clinical$batch),
                         stage = "Before")
df_after <- data.frame(PC1 = pca_after[,1], PC2 = pca_after[,2],
                        batch = factor(clinical$batch),
                        stage = "After")
df_pca <- rbind(df_before, df_after)

p_pca <- ggplot(df_pca, aes(x = PC1, y = PC2, color = batch)) +
  geom_point(alpha = 0.6, size = 1.2) +
  facet_wrap(~ stage, scales = "free") +
  theme_bw(base_size = 12) +
  labs(title = "PCA before/after batch correction", color = "Batch")
save_fig(p_pca, DIR_THIS, "Fig1_PCA_before_after", width = 8, height = 4)

p_age <- ggplot(clinical, aes(x = factor(grade), y = age_n)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, fill = "#B0B0B0") +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
  theme_bw(base_size = 12) +
  labs(x = "AMD grade", y = "Age (years)",
       title = "Age distribution by AMD grade")
save_fig(p_age, DIR_THIS, "Fig2_age_by_grade", width = 5, height = 4)

# ---------- 14. summary ----------
cat("\n========================================\n")
cat("=== 01_1_prepare_bulk_data.r 汇总 ===\n")
cat("========================================\n")
cat("样本数:", nrow(clinical), "\n")
cat("基因数:", nrow(voom_full), "\n\n")
cat("grade 分布:\n"); print(table(clinical$grade))
cat("\nbatch 分布:\n"); print(table(clinical$batch))
cat("\nRIN 中位数:", median(clinical$rin_n, na.rm = TRUE), "\n")
cat("age 中位数:", median(clinical$age_n, na.rm = TRUE), "\n\n")
cat("PC1-5 方差占比:", paste(var_pct[1:5], collapse = ", "), "\n\n")
cat("PC 诊断:\n"); print(pc_diag, row.names = FALSE)

cat("\n--- 输出位置 ---\n")
cat("目录:", DIR_THIS, "\n")
cat("  data/expr_count.rds\n")
cat("  data/voom_full.rds\n")
cat("  data/expr_corrected.rds\n")
cat("  data/clinical.rds\n")
cat("  data/gene_annotation.rds\n")

log_info("=== 01_1_prepare_bulk_data.r done ===")
