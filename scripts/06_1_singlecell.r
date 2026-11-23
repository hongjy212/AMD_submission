# ============================================================
# 06_1_singlecell.r (v6 final)
#
# ============================================================

rm(list = ls()); gc()
setwd("E:/retina_tools/aging_rgc_axon/AMD")
set.seed(2026)

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(Matrix); library(data.table)
  library(ggplot2); library(dplyr); library(tidyr)
})

if (!exists("setup_output_dir")) {
  setup_output_dir <- function(name) {
    d <- file.path("results", name)
    dir.create(file.path(d, "tables"),  recursive = TRUE, showWarnings = FALSE)
    dir.create(file.path(d, "figures"), recursive = TRUE, showWarnings = FALSE)
    if (exists("log_info")) log_info("输出目录: {d}")
    d
  }
}
if (!exists("save_table")) {
  save_table <- function(df, dir, name) {
    write.csv(df, file.path(dir, "tables", name),
              row.names = FALSE, fileEncoding = "UTF-8")
  }
}
if (!exists("save_fig")) {
  save_fig <- function(p, dir, name, w = 7, h = 5) {
    ggsave(file.path(dir, "figures", paste0(name, ".pdf")), p,
           width = w, height = h)
  }
}

DIR_THIS <- setup_output_dir("06_1_singlecell")
log_info("=== 06_1 singlecell v6 final start ===")

DIR_SC  <- "E:/retina_tools/raw/GSE221042"
mtx_f   <- file.path(DIR_SC, "GSE221042_data.mtx.gz")
genes_f <- file.path(DIR_SC, "GSE221042_genes.csv.gz")
meta_f  <- file.path(DIR_SC, "GSE221042_sample_annotations.csv.gz")
stopifnot(file.exists(mtx_f), file.exists(genes_f), file.exists(meta_f))

# ============================================================
# 1. read
# ============================================================
log_info("--- 1. 读入 ---")

genes <- fread(genes_f, header = FALSE, data.table = FALSE)
colnames(genes)[1:2] <- c("idx", "gene_full")
genes$symbol <- toupper(trimws(sub(" \\(.*$", "", genes$gene_full)))
log_info("genes: {nrow(genes)}")

meta <- fread(meta_f, header = TRUE, data.table = FALSE)
log_info("meta 原始列名: {paste(colnames(meta), collapse=', ')}")
if (ncol(meta) >= 5) {
  colnames(meta)[1:5] <- c("cell_id","donor","condition","n_genes","total_counts")
}
log_info("meta: {nrow(meta)} × {ncol(meta)}")
log_info("condition: {paste(names(table(meta$condition)), table(meta$condition), sep='=', collapse=' ')}")

mtx <- readMM(mtx_f)
log_info("mtx: {nrow(mtx)} × {ncol(mtx)}")

if (nrow(mtx) == nrow(meta) && ncol(mtx) == nrow(genes)) {
  log_info("★ mtx = cells × genes")
} else if (ncol(mtx) == nrow(meta) && nrow(mtx) == nrow(genes)) {
  log_info("★ mtx = genes × cells → 转置")
  mtx <- t(mtx)
} else {
  stop(sprintf("mtx 与 meta/genes 不匹配: mtx %d×%d, meta %d, genes %d",
               nrow(mtx), ncol(mtx), nrow(meta), nrow(genes)))
}
rownames(mtx) <- meta$cell_id
colnames(mtx) <- genes$symbol
mtx <- as(mtx, "CsparseMatrix")

# ============================================================
# 2. QC
# ============================================================
log_info("--- 2. QC ---")
total_counts <- Matrix::rowSums(mtx)
n_genes_cell <- Matrix::rowSums(mtx > 0)
mt_genes <- grep("^MT-", colnames(mtx), value = TRUE)
mt_pct <- if (length(mt_genes) > 0) {
  Matrix::rowSums(mtx[, mt_genes, drop = FALSE]) / total_counts * 100
} else rep(0, nrow(mtx))

keep <- n_genes_cell >= 200 & n_genes_cell <= 7500 & mt_pct <= 10
mtx  <- mtx[keep, ]
meta <- meta[keep, ]
log_info("QC 后: {nrow(mtx)} cells")

# ============================================================
# 3. normalization
# ============================================================
libsize <- Matrix::rowSums(mtx)
norm <- mtx / libsize * 1e4
norm <- log1p(norm)

# ============================================================
# ============================================================
log_info("--- 4. Marker 打分 ---")
MARKERS <- list(
  Rod  = c("RHO","GNAT1","PDE6B","NRL"),
  Cone = c("ARR3","GNAT2","PDE6H"),
  HC   = c("ONECUT1","ONECUT2","LHX1","TFAP2B"),
  BC   = c("VSX1","VSX2","GRM6"),
  AC   = c("GAD1","GAD2"),
  MG   = c("RLBP1","GLUL"),
  RGC  = c("RBPMS","NEFL","NEFM"),
  RPE  = c("RPE65","BEST1","TTR")
)
MARKERS <- lapply(MARKERS, function(g) intersect(g, colnames(norm)))
MARKERS <- MARKERS[sapply(MARKERS, length) >= 1]
log_info("有效 cell type: {length(MARKERS)}")

scores <- sapply(MARKERS, function(g) {
  cols <- which(colnames(norm) %in% g)
  if (length(cols) == 0) return(rep(0, nrow(norm)))
  Matrix::rowMeans(norm[, cols, drop = FALSE])
})
if (is.null(dim(scores))) {
  scores <- matrix(scores, ncol = length(MARKERS),
                   dimnames = list(NULL, names(MARKERS)))
}
stopifnot(nrow(scores) == nrow(norm))

max_score <- apply(scores, 1, max)
max_type  <- colnames(scores)[apply(scores, 1, which.max)]
sorted2   <- t(apply(scores, 1, function(x) sort(x, decreasing = TRUE)[1:2]))
gap       <- sorted2[, 1] - sorted2[, 2]
is_ambig  <- gap < 0.05
is_low    <- max_score < 0.1

meta$cell_type <- max_type
meta$cell_type[is_ambig] <- "Ambiguous"
meta$cell_type[is_low]   <- "Low_quality"
meta$marker_max <- max_score
meta$marker_gap <- gap
log_info("cell type 分布:"); print(table(meta$cell_type))

# ============================================================
# ============================================================
log_info("--- 5. donor × celltype 比例 ---")
meta_clean <- meta[!meta$cell_type %in% c("Ambiguous","Low_quality"), ]
log_info("清理后 cells: {nrow(meta_clean)}")

donor_pct <- meta_clean %>%
  group_by(donor, condition, celltype = cell_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(donor) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()
save_table(donor_pct, DIR_THIS, "donor_celltype_pct.csv")

# ============================================================
# 6. pathwayscore
# ============================================================
log_info("--- 6. 通路评分 ---")
GENE_SETS <- readRDS("results/01_2_prepare_gene_sets/data/gene_sets.rds")
oxphos_genes     <- toupper(GENE_SETS[["OXPHOS"]]$genes)
complement_genes <- toupper(GENE_SETS[["Complement"]]$genes)
genes16 <- toupper(read.csv(
  "results/02_1_differential_expression_main/tables/sig16_genes.csv")$gene)

pathway_list <- list(
  OXPHOS      = intersect(oxphos_genes, colnames(norm)),
  Complement  = intersect(complement_genes, colnames(norm)),
  Signature16 = intersect(genes16, colnames(norm))
)
log_info("OXPHOS: {length(pathway_list$OXPHOS)}  Complement: {length(pathway_list$Complement)}  Sig16: {length(pathway_list$Signature16)}")

pathway_scores <- sapply(pathway_list, function(g) {
  cols <- which(colnames(norm) %in% g)
  if (length(cols) == 0) return(rep(0, nrow(norm)))
  Matrix::rowMeans(norm[, cols, drop = FALSE])
})
if (is.null(dim(pathway_scores))) {
  pathway_scores <- matrix(pathway_scores, ncol = length(pathway_list),
                           dimnames = list(NULL, names(pathway_list)))
}
stopifnot(nrow(pathway_scores) == nrow(norm))

idx <- match(meta_clean$cell_id, rownames(norm))
stopifnot(sum(is.na(idx)) == 0)
meta_clean$OXPHOS      <- pathway_scores[idx, "OXPHOS"]
meta_clean$Complement  <- pathway_scores[idx, "Complement"]
meta_clean$Signature16 <- pathway_scores[idx, "Signature16"]
log_info("meta_clean OXPHOS NA 数: {sum(is.na(meta_clean$OXPHOS))}")

donor_ct_pathway <- meta_clean %>%
  group_by(donor, condition, celltype = cell_type) %>%
  summarise(
    n_cells          = n(),
    mean_OXPHOS      = mean(OXPHOS, na.rm = TRUE),
    mean_Complement  = mean(Complement, na.rm = TRUE),
    mean_Signature16 = mean(Signature16, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(n_cells >= 20)
save_table(donor_ct_pathway, DIR_THIS, "celltype_pathway_scores.csv")
log_info("donor_ct_pathway: {nrow(donor_ct_pathway)} 行")

# ============================================================
# 7. condition compare
# ============================================================
log_info("--- 7. condition 比较 ---")

ct_test <- data.frame()
for (ct in unique(donor_pct$celltype)) {
  sub <- donor_pct[donor_pct$celltype == ct, ]
  h <- na.omit(sub$pct[sub$condition %in% c("Healthy","Control")])
  d <- na.omit(sub$pct[sub$condition == "dAMD"])
  if (length(h) >= 3 && length(d) >= 3) {
    tt <- tryCatch(t.test(d, h), error = function(e) NULL)
    if (is.null(tt)) next
    ct_test <- rbind(ct_test, data.frame(
      celltype = ct, n_healthy = length(h), n_damd = length(d),
      mean_healthy = mean(h), mean_damd = mean(d),
      delta = mean(d) - mean(h), p = tt$p.value,
      stringsAsFactors = FALSE))
  }
}
if (nrow(ct_test) > 0) {
  ct_test$FDR <- p.adjust(ct_test$p, method = "BH")
  ct_test <- ct_test[order(ct_test$p), ]
}
save_table(ct_test, DIR_THIS, "celltype_by_condition_ttest.csv")
log_info("celltype 比例差异:"); print(ct_test, row.names = FALSE, digits = 4)

pt_test <- data.frame()
for (ct in unique(donor_ct_pathway$celltype)) {
  for (metric in c("mean_OXPHOS","mean_Complement","mean_Signature16")) {
    sub <- donor_ct_pathway[donor_ct_pathway$celltype == ct, ]
    h <- na.omit(sub[[metric]][sub$condition %in% c("Healthy","Control")])
    d <- na.omit(sub[[metric]][sub$condition == "dAMD"])
    if (length(h) >= 3 && length(d) >= 3) {
      tt <- tryCatch(t.test(d, h), error = function(e) NULL)
      if (is.null(tt)) next
      pt_test <- rbind(pt_test, data.frame(
        celltype = ct, metric = metric,
        n_healthy = length(h), n_damd = length(d),
        mean_healthy = mean(h), mean_damd = mean(d),
        delta = mean(d) - mean(h), p = tt$p.value,
        stringsAsFactors = FALSE))
    }
  }
}
if (nrow(pt_test) > 0) {
  pt_test$FDR <- p.adjust(pt_test$p, method = "BH")
  pt_test <- pt_test[order(pt_test$p), ]
}
save_table(pt_test, DIR_THIS, "celltype_pathway_condition_test.csv")
log_info("细胞类型内通路差异:"); print(pt_test, row.names = FALSE, digits = 4)

# ============================================================
# 8. figure
# ============================================================
log_info("--- 8. 图 ---")

if (nrow(donor_pct) > 0) {
  p1 <- ggplot(donor_pct, aes(x = celltype, y = pct, fill = condition)) +
    geom_boxplot(outlier.size = 0.5, linewidth = 0.3) +
    labs(title = "Cell type proportions by condition",
         x = "", y = "% of cells") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_fig(p1, DIR_THIS, "celltype_by_condition", 8, 5)
}

# 8.2 OXPHOS delta
if (nrow(pt_test) > 0) {
  oxphos_test <- pt_test[pt_test$metric == "mean_OXPHOS", ]
  if (nrow(oxphos_test) > 0) {
    oxphos_test$celltype <- factor(oxphos_test$celltype,
      levels = oxphos_test$celltype[order(oxphos_test$delta)])
    p2 <- ggplot(oxphos_test, aes(x = celltype, y = delta, fill = delta < 0)) +
      geom_col(width = 0.6) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_text(aes(label = sprintf("p=%.2g", p)),
                vjust = ifelse(oxphos_test$delta > 0, -0.5, 1.2), size = 2.5) +
      scale_fill_manual(values = c("TRUE" = "#3B6DB5", "FALSE" = "#C44E52"),
                         guide = "none") +
      labs(title = "OXPHOS score change (dAMD - Healthy) by cell type",
           x = "", y = "Delta OXPHOS") +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    save_fig(p2, DIR_THIS, "celltype_OXPHOS_change", 7, 4)
  }

  # 8.3 Complement delta
  compl_test <- pt_test[pt_test$metric == "mean_Complement", ]
  if (nrow(compl_test) > 0) {
    compl_test$celltype <- factor(compl_test$celltype,
      levels = compl_test$celltype[order(compl_test$delta)])
    p3 <- ggplot(compl_test, aes(x = celltype, y = delta, fill = delta > 0)) +
      geom_col(width = 0.6) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      geom_text(aes(label = sprintf("p=%.2g", p)),
                vjust = ifelse(compl_test$delta > 0, -0.5, 1.2), size = 2.5) +
      scale_fill_manual(values = c("TRUE" = "#C44E52", "FALSE" = "#3B6DB5"),
                         guide = "none") +
      labs(title = "Complement score change (dAMD - Healthy) by cell type",
           x = "", y = "Delta Complement") +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    save_fig(p3, DIR_THIS, "celltype_Complement_change", 7, 4)
  }
}

# ============================================================
# 9. summary
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("=== 06_1 v6 final 汇总 ===\n")
cat(strrep("=", 70), "\n", sep = "")
cat("\n细胞数:", nrow(mtx), "\n")
cat("清理后:", nrow(meta_clean), "\n")
cat("细胞类型:\n"); print(table(meta$cell_type))
cat("\ncondition 细胞比例差异:\n")
print(ct_test, row.names = FALSE, digits = 4)
cat("\n细胞类型内通路评分差异:\n")
print(pt_test, row.names = FALSE, digits = 4)
cat("\n输出目录:", DIR_THIS, "\n")
log_info("=== 06_1 v6 final done ===")

