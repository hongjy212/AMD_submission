#!/usr/bin/env Rscript
# ============================================================
# 03_1_celltype_marker_score.r
#
#
# Input：
#   results/01_1_prepare_bulk_data/data/{expr_corrected,clinical,gene_annotation}.rds
#   results/02_1_differential_expression_main/tables/poly3_F_full.csv
#   results/02_1_differential_expression_main/tables/limma_multi_full.csv
#   E:/retina_tools/pubdata/retina_fingerprint/41467_2024_50853_MOESM6_ESM.xlsx
#
#   data/celltype_marker_scores.rds
#   tables/celltype_marker_scores.csv
#   tables/marker_score_correlation.csv
#   tables/polr_by_celltype.csv
#   tables/poly3F_by_celltype.csv
#   tables/poly3_components_by_celltype.csv
#   tables/celltype_zigzag_check.csv
#   tables/lm_Rod_Cone_age_grade.csv
#   tables/celltype_gsea_poly3F.csv
#   tables/celltype_gsea_g3v2.csv
#   tables/celltype_gsea_g4v1.csv
#   tables/celltype_gsea_combined.csv
#   figures/polr_logOR.pdf/.png
#   figures/poly3F_bar.pdf/.png
#   figures/celltype_poly3_curves.pdf/.png
#   figures/celltype_scores_by_grade.pdf/.png
#   figures/celltype_gsea_heatmap.pdf/.png
#   figures/celltype_gsea_poly3F_bar.pdf/.png
#   figures/celltype_gsea_g3v2_bar.pdf/.png
# ============================================================
rm(list = ls(all.names = TRUE)); gc()
source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(ggplot2)
  library(MASS); library(splines); library(stringr)
  library(fgsea); library(pheatmap)
})

DIR_THIS <- setup_output_dir("03_1_celltype_marker_score")
log_info("=== 03_1_celltype_marker_score.r start ===")

# Parameters
TOP_N_SCORE  <- 20    # marker score 用 top 20（加权平均）
TOP_N_MARKERS <- 100  # GSEA 用 top 100（集合更大）

# ============================================================
# ============================================================
log_info("--- 1. 加载 bulk 数据 ---")

expr_corrected <- read_upstream("01_1_prepare_bulk_data", "data",
                                 "expr_corrected.rds")
gene_anno <- read_upstream("01_1_prepare_bulk_data", "data",
                            "gene_annotation.rds")
clinical <- read_upstream("01_1_prepare_bulk_data", "data", "clinical.rds")

clinical$amd_grade <- as.integer(clinical$grade)
clinical$age       <- clinical$age_n
clinical$RIN       <- clinical$rin_n
clinical$PMI       <- clinical$pmi_n

log_info("expr: {nrow(expr_corrected)} × {ncol(expr_corrected)}")
log_info("clinical: {nrow(clinical)} 样本")

# ============================================================
# ============================================================
log_info("--- 2. 载入 02_1 的排序结果 ---")

poly3F_file <- file.path(DIR_RESULT, "02_1_differential_expression_main",
                          "tables", "poly3_F_full.csv")
limma_file <- file.path(DIR_RESULT, "02_1_differential_expression_main",
                         "tables", "limma_multi_full.csv")

if (!file.exists(poly3F_file)) stop("缺 poly3_F_full.csv，请先跑 02_1")
if (!file.exists(limma_file)) stop("缺 limma_multi_full.csv，请先跑 02_1")

poly3F <- read.csv(poly3F_file, stringsAsFactors = FALSE)
limma_multi <- read.csv(limma_file, stringsAsFactors = FALSE)
log_info("poly3_F: {nrow(poly3F)} 基因")
log_info("limma_multi: {nrow(limma_multi)} 行")

# ============================================================
# ============================================================
log_info("--- 3. 读 MOESM6 ---")

celltype_file <- "E:/retina_tools/pubdata/retina_fingerprint/41467_2024_50853_MOESM6_ESM.xlsx"
if (!file.exists(celltype_file)) stop("缺 MOESM6")

sheet_names <- readxl::excel_sheets(celltype_file)
log_info("表单: {paste(sheet_names, collapse=', ')}")

all_sheets <- list()
for (sheet in sheet_names) {
  df <- readxl::read_excel(celltype_file, sheet = sheet)
  if (!all(c("names", "scores") %in% colnames(df))) next
  df <- df[!duplicated(df$names), ]
  df <- df[order(-df$scores), ]
  all_sheets[[sheet]] <- df
}

# ============================================================
# ============================================================
log_info("--- 4. marker score (TOP_N_SCORE = {TOP_N_SCORE}) ---")

compute_weighted_score <- function(df, expr, gene_anno, top_n) {
  df <- head(df, top_n)

  ensembl <- gene_anno$ensembl[match(toupper(df$names),
                                       toupper(gene_anno$symbol))]
  valid <- !is.na(ensembl) & ensembl %in% rownames(expr)
  ensembl <- ensembl[valid]
  weights <- df$scores[valid]

  if (length(ensembl) < 5) return(NULL)

  sub <- expr[ensembl, , drop = FALSE]
  sub_z <- t(scale(t(sub)))
  score <- colSums(sub_z * weights) / sum(weights)
  attr(score, "n_matched") <- length(ensembl)
  return(score)
}

celltype_scores_list <- list()
celltype_n_genes <- list()
for (ct in names(all_sheets)) {
  score <- compute_weighted_score(all_sheets[[ct]], expr_corrected,
                                    gene_anno, TOP_N_SCORE)
  if (!is.null(score)) {
    celltype_scores_list[[ct]] <- score
    celltype_n_genes[[ct]] <- attr(score, "n_matched")
    log_info("  {ct}: 匹配 {attr(score, 'n_matched')} / {TOP_N_SCORE}")
  }
}

if (length(celltype_scores_list) == 0) stop("无 celltype 数据")

# Rod vs Cone validation
if (all(c("Rod", "Cone") %in% names(celltype_scores_list))) {
  r_rc <- cor(celltype_scores_list$Rod, celltype_scores_list$Cone)
  log_info("Rod vs Cone: {round(r_rc, 4)}")
}

# build marker_score data.frame
marker_score_df <- data.frame(
  sample_id = names(celltype_scores_list[[1]]),
  stringsAsFactors = FALSE)
for (ct in names(celltype_scores_list)) {
  marker_score_df[[paste0(ct, "_score")]] <-
    as.numeric(celltype_scores_list[[ct]])
}

merged <- merge(marker_score_df, clinical, by = "sample_id")
merged$grade_factor <- factor(merged$amd_grade)
merged$grade <- merged$amd_grade

cell_types <- names(celltype_scores_list)

# save score
save_rds(marker_score_df, DIR_THIS, "celltype_marker_scores.rds")
save_table(marker_score_df, DIR_THIS, "celltype_marker_scores.csv")

score_mat <- as.matrix(marker_score_df[, -1, drop = FALSE])
cor_m <- cor(score_mat, use = "pairwise.complete.obs")
save_table(as.data.frame(round(cor_m, 4)), DIR_THIS,
           "marker_score_correlation.csv")

# ============================================================
# ============================================================
log_info("--- 5. polr（单调） ---")

polr_results <- data.frame()
for (ct in cell_types) {
  score_col <- paste0(ct, "_score")
  fml <- as.formula(paste(
    "grade_factor ~", score_col,
    "+ age + sex + batch + RIN + PMI + PC1 + PC2 + PC3 + PC4 + PC5"))
  fit <- tryCatch(polr(fml, data = merged, Hess = TRUE),
                   error = function(e) NULL)
  if (is.null(fit)) next

  co <- summary(fit)$coefficients
  if (!score_col %in% rownames(co)) next

  est <- co[score_col, "Value"]
  tval <- co[score_col, "t value"]
  pval <- 2 * pnorm(-abs(tval))

  polr_results <- rbind(polr_results, data.frame(
    celltype = ct,
    n_genes = celltype_n_genes[[ct]],
    Estimate = round(est, 4),
    OR = round(exp(est), 4),
    t = round(tval, 3),
    p = pval,
    stringsAsFactors = FALSE))
}
polr_results$FDR <- p.adjust(polr_results$p, "BH")
polr_results <- polr_results[order(polr_results$p), ]
save_table(polr_results, DIR_THIS, "polr_by_celltype.csv")

cat("\n=== polr（单调） ===\n")
print(polr_results, row.names = FALSE, digits = 4)

# ============================================================
# ============================================================
log_info("--- 6. poly3_F（非单调） ---")

poly_mat <- poly(merged$grade, 3)
colnames(poly_mat) <- c("grade_p1", "grade_p2", "grade_p3")
merged <- cbind(merged, poly_mat)

sp <- splines::ns(merged$age_n, df = 3)
merged$age_s1 <- sp[,1]; merged$age_s2 <- sp[,2]; merged$age_s3 <- sp[,3]

poly3_results <- data.frame()
poly3_components <- data.frame()

for (ct in cell_types) {
  score_col <- paste0(ct, "_score")
  y <- merged[[score_col]]

  df <- data.frame(
    y = y,
    grade_p1 = merged$grade_p1,
    grade_p2 = merged$grade_p2,
    grade_p3 = merged$grade_p3,
    age_s1 = merged$age_s1,
    age_s2 = merged$age_s2,
    age_s3 = merged$age_s3,
    sex_M = as.integer(merged$sex == "Male"),
    rin = merged$rin_n,
    pmi = merged$pmi_n)
  df <- df[complete.cases(df), ]

  fit <- lm(y ~ grade_p1 + grade_p2 + grade_p3 +
              age_s1 + age_s2 + age_s3 + sex_M + rin + pmi, data = df)
  fit_null <- lm(y ~ age_s1 + age_s2 + age_s3 + sex_M + rin + pmi, data = df)
  aov_cmp <- anova(fit_null, fit)

  poly3_results <- rbind(poly3_results, data.frame(
    celltype = ct,
    F_poly3 = round(aov_cmp$F[2], 3),
    p_F = aov_cmp$`Pr(>F)`[2],
    stringsAsFactors = FALSE))

  co <- summary(fit)$coefficients
  poly3_components <- rbind(poly3_components, data.frame(
    celltype = ct,
    t_p1 = round(co["grade_p1", "t value"], 3),
    p_p1 = co["grade_p1", "Pr(>|t|)"],
    t_p2 = round(co["grade_p2", "t value"], 3),
    p_p2 = co["grade_p2", "Pr(>|t|)"],
    t_p3 = round(co["grade_p3", "t value"], 3),
    p_p3 = co["grade_p3", "Pr(>|t|)"],
    stringsAsFactors = FALSE))
}

poly3_results$FDR <- p.adjust(poly3_results$p_F, "BH")
poly3_results <- poly3_results[order(poly3_results$p_F), ]
save_table(poly3_results, DIR_THIS, "poly3F_by_celltype.csv")

poly3_components$FDR_p1 <- p.adjust(poly3_components$p_p1, "BH")
poly3_components$FDR_p2 <- p.adjust(poly3_components$p_p2, "BH")
poly3_components$FDR_p3 <- p.adjust(poly3_components$p_p3, "BH")
save_table(poly3_components, DIR_THIS, "poly3_components_by_celltype.csv")

cat("\n=== poly3_F（非单调） ===\n")
print(poly3_results, row.names = FALSE, digits = 4)

# ============================================================
# ============================================================
log_info("--- 7. 锯齿检查 ---")

curve_data <- data.frame()
for (ct in cell_types) {
  score_col <- paste0(ct, "_score")
  for (g in 1:4) {
    sub <- merged[merged$grade == g, score_col]
    curve_data <- rbind(curve_data, data.frame(
      celltype = ct, grade = g,
      mean_score = mean(sub, na.rm = TRUE),
      se = sd(sub, na.rm = TRUE) / sqrt(length(sub)),
      stringsAsFactors = FALSE))
  }
}

zigzag_tbl <- data.frame()
for (ct in cell_types) {
  cd <- curve_data[curve_data$celltype == ct, ]
  cd <- cd[order(cd$grade), ]
  zigzag_tbl <- rbind(zigzag_tbl, data.frame(
    celltype = ct,
    g1 = round(cd$mean_score[cd$grade == 1], 4),
    g2 = round(cd$mean_score[cd$grade == 2], 4),
    g3 = round(cd$mean_score[cd$grade == 3], 4),
    g4 = round(cd$mean_score[cd$grade == 4], 4),
    is_zigzag = cd$mean_score[cd$grade == 2] > cd$mean_score[cd$grade == 3],
    stringsAsFactors = FALSE))
}
save_table(zigzag_tbl, DIR_THIS, "celltype_zigzag_check.csv")

cat("\n=== 锯齿检查 ===\n")
print(zigzag_tbl, row.names = FALSE)

# ============================================================
# 8. lm（Rod/Cone ~ age + grade）
# ============================================================
log_info("--- 8. lm Rod/Cone ---")

lm_results <- data.frame()
for (ct in c("Rod", "Cone")) {
  score_col <- paste0(ct, "_score")
  if (!score_col %in% colnames(merged)) next

  fit <- lm(as.formula(paste(
    score_col, "~ age + grade_factor + sex + batch + RIN + PMI",
    "+ PC1 + PC2 + PC3 + PC4 + PC5")), data = merged)
  co <- summary(fit)$coefficients

  lm_results <- rbind(lm_results, data.frame(
    celltype = ct,
    age_beta = round(co["age", "Estimate"], 6),
    age_p = co["age", "Pr(>|t|)"],
    grade3_beta = round(co["grade_factor3", "Estimate"], 6),
    grade3_p = co["grade_factor3", "Pr(>|t|)"],
    grade4_beta = round(co["grade_factor4", "Estimate"], 6),
    grade4_p = co["grade_factor4", "Pr(>|t|)"],
    stringsAsFactors = FALSE))
}
save_table(lm_results, DIR_THIS, "lm_Rod_Cone_age_grade.csv")

# ============================================================
# ============================================================
log_info("--- 9. GSEA (top {TOP_N_MARKERS} markers) ---")

# build marker sets（top 100）
marker_sets <- list()
for (ct in names(all_sheets)) {
  df <- head(all_sheets[[ct]], TOP_N_MARKERS)
  genes <- toupper(trimws(df$names))
  genes <- genes[!is.na(genes) & genes != ""]
  marker_sets[[ct]] <- unique(genes)
  log_info("  {ct}: {length(unique(genes))} markers")
}

ranks_poly3F <- poly3F$t
names(ranks_poly3F) <- poly3F$gene
ranks_poly3F <- sort(ranks_poly3F[!is.na(ranks_poly3F)], decreasing = TRUE)

lm_g3v2 <- limma_multi[limma_multi$contrast == "g3_vs_g2", ]
ranks_g3v2 <- lm_g3v2$t
names(ranks_g3v2) <- lm_g3v2$gene
ranks_g3v2 <- sort(ranks_g3v2[!is.na(ranks_g3v2)], decreasing = TRUE)

lm_g4v1 <- limma_multi[limma_multi$contrast == "g4_vs_g1", ]
ranks_g4v1 <- lm_g4v1$t
names(ranks_g4v1) <- lm_g4v1$gene
ranks_g4v1 <- sort(ranks_g4v1[!is.na(ranks_g4v1)], decreasing = TRUE)

run_gsea_one <- function(ranks, marker_sets, ranking_name) {
  set.seed(SEED)
  res <- tryCatch({
    fgsea::fgsea(pathways = marker_sets, stats = ranks,
                  minSize = 10, maxSize = 500, nperm = 10000)
  }, error = function(e) NULL)
  if (is.null(res)) return(NULL)
  res <- as.data.frame(res)
  res$leadingEdge <- NULL
  res$ranking <- ranking_name
  res[order(res$padj), ]
}

gsea_poly3F <- run_gsea_one(ranks_poly3F, marker_sets, "poly3_F")
gsea_g3v2   <- run_gsea_one(ranks_g3v2,   marker_sets, "g3_vs_g2")
gsea_g4v1   <- run_gsea_one(ranks_g4v1,   marker_sets, "g4_vs_g1")

# save
if (!is.null(gsea_poly3F)) {
  save_table(gsea_poly3F, DIR_THIS, "celltype_gsea_poly3F.csv")
  cat("\n=== GSEA poly3_F ===\n")
  print(gsea_poly3F[, c("pathway","NES","pval","padj","size")],
        row.names = FALSE, digits = 4)
}
if (!is.null(gsea_g3v2)) {
  save_table(gsea_g3v2, DIR_THIS, "celltype_gsea_g3v2.csv")
  cat("\n=== GSEA g3_vs_g2 ===\n")
  print(gsea_g3v2[, c("pathway","NES","pval","padj","size")],
        row.names = FALSE, digits = 4)
}
if (!is.null(gsea_g4v1)) {
  save_table(gsea_g4v1, DIR_THIS, "celltype_gsea_g4v1.csv")
  cat("\n=== GSEA g4_vs_g1 ===\n")
  print(gsea_g4v1[, c("pathway","NES","pval","padj","size")],
        row.names = FALSE, digits = 4)
}

all_gsea <- rbind(
  if (!is.null(gsea_poly3F)) gsea_poly3F else NULL,
  if (!is.null(gsea_g3v2))   gsea_g3v2   else NULL,
  if (!is.null(gsea_g4v1))   gsea_g4v1   else NULL)

if (!is.null(all_gsea) && nrow(all_gsea) > 0) {
  combined <- tidyr::pivot_wider(
    all_gsea[, c("pathway", "ranking", "NES", "padj")],
    names_from = ranking,
    values_from = c(NES, padj))
  save_table(as.data.frame(combined), DIR_THIS, "celltype_gsea_combined.csv")

  cat("\n=== 三排序合并 ===\n")
  print(as.data.frame(combined), row.names = FALSE, digits = 3)
}

# ============================================================
# 10. figure
# ============================================================
log_info("--- 10. 图 ---")

plot_df <- polr_results
plot_df$celltype <- factor(plot_df$celltype,
                            levels = plot_df$celltype[order(plot_df$OR)])
plot_df$log_OR <- log(plot_df$OR)

p1 <- ggplot(plot_df, aes(x = celltype, y = log_OR, fill = p < 0.05)) +
  geom_col() + coord_flip() + theme_bw(base_size = 11) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("TRUE" = "#C44E52", "FALSE" = "#B0B0B0"),
                    labels = c("TRUE" = "p<0.05", "FALSE" = "NS")) +
  labs(x = "", y = "log(OR)", fill = "",
       title = "polr: AMD grade ~ celltype marker score (monotonic)")
save_fig(p1, DIR_THIS, "polr_logOR", width = 7, height = 5)

plot_df2 <- poly3_results
plot_df2$celltype <- factor(plot_df2$celltype,
                              levels = plot_df2$celltype[order(plot_df2$F_poly3)])
plot_df2$sig <- ifelse(plot_df2$FDR < 0.05, "FDR<0.05",
                 ifelse(plot_df2$p_F < 0.05, "p<0.05", "ns"))

p2 <- ggplot(plot_df2, aes(x = celltype, y = F_poly3, fill = sig)) +
  geom_col() + coord_flip() + theme_bw(base_size = 11) +
  scale_fill_manual(values = c("FDR<0.05" = "#C44E52",
                                "p<0.05" = "#FDB462",
                                "ns" = "#B0B0B0")) +
  labs(x = "", y = "poly3_F F statistic", fill = "",
       title = "poly3_F: celltype marker score ~ grade (non-monotonic)")
save_fig(p2, DIR_THIS, "poly3F_bar", width = 7, height = 5)

p3 <- ggplot(curve_data, aes(x = grade, y = mean_score,
                               color = celltype, group = celltype)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean_score - se, ymax = mean_score + se),
                 width = 0.15, alpha = 0.5) +
  theme_bw(base_size = 12) +
  scale_x_continuous(breaks = 1:4, labels = paste0("G", 1:4)) +
  labs(x = "AMD grade", y = "Marker score (mean ± SE)",
       title = "Celltype marker score by grade", color = "Cell type")
save_fig(p3, DIR_THIS, "celltype_poly3_curves", width = 9, height = 6)

score_cols <- grep("_score$", colnames(merged), value = TRUE)
score_long <- merged %>%
  dplyr::select(sample_id, grade_factor, all_of(score_cols)) %>%
  pivot_longer(cols = all_of(score_cols),
                names_to = "celltype", values_to = "score") %>%
  mutate(celltype = sub("_score$", "", celltype)) %>%
  filter(!is.na(score))

p4 <- ggplot(score_long, aes(x = grade_factor, y = score, fill = grade_factor)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.3) +
  facet_wrap(~ celltype, scales = "free_y", ncol = 3) +
  scale_fill_brewer(palette = "RdYlBu") +
  labs(title = "Cell type marker scores by AMD grade (top 20)",
       x = "AMD grade", y = "Marker score (z)") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none")
save_fig(p4, DIR_THIS, "celltype_scores_by_grade", width = 11, height = 8)

if (exists("combined")) {
  nes_cols <- grep("^NES_", colnames(combined), value = TRUE)
  if (length(nes_cols) > 0) {
    nes_mat <- as.matrix(combined[, nes_cols, drop = FALSE])
    rownames(nes_mat) <- combined$pathway
    colnames(nes_mat) <- sub("^NES_", "", nes_cols)

    if ("poly3_F" %in% colnames(nes_mat)) {
      nes_mat <- nes_mat[order(nes_mat[, "poly3_F"]), ]
    }

    pheatmap(nes_mat,
             cluster_rows = FALSE, cluster_cols = FALSE,
             display_numbers = TRUE, number_format = "%.2f",
             color = colorRampPalette(c("#4C72B0", "#FFFFFF", "#C44E52"))(50),
             main = "Celltype marker GSEA NES (3 rankings)",
             filename = file.path(DIR_THIS, "figures",
                                   "celltype_gsea_heatmap.png"),
             width = 6, height = 5)
  }
}

if (!is.null(gsea_poly3F)) {
  plot_df <- gsea_poly3F
  plot_df$pathway <- factor(plot_df$pathway,
                              levels = plot_df$pathway[order(plot_df$NES)])
  plot_df$sig <- ifelse(plot_df$padj < 0.05, "padj<0.05",
                  ifelse(plot_df$pval < 0.05, "p<0.05", "ns"))

  p <- ggplot(plot_df, aes(x = pathway, y = NES, fill = sig)) +
    geom_col() + coord_flip() + theme_bw(base_size = 11) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = c("padj<0.05" = "#C44E52",
                                  "p<0.05" = "#FDB462",
                                  "ns" = "#B0B0B0")) +
    labs(x = "", y = "NES (poly3_F ranking)", fill = "",
         title = "Celltype marker GSEA (poly3_F)")
  save_fig(p, DIR_THIS, "celltype_gsea_poly3F_bar", width = 7, height = 5)
}

if (!is.null(gsea_g3v2)) {
  plot_df <- gsea_g3v2
  plot_df$pathway <- factor(plot_df$pathway,
                              levels = plot_df$pathway[order(plot_df$NES)])
  plot_df$sig <- ifelse(plot_df$padj < 0.05, "padj<0.05",
                  ifelse(plot_df$pval < 0.05, "p<0.05", "ns"))

  p <- ggplot(plot_df, aes(x = pathway, y = NES, fill = sig)) +
    geom_col() + coord_flip() + theme_bw(base_size = 11) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = c("padj<0.05" = "#C44E52",
                                  "p<0.05" = "#FDB462",
                                  "ns" = "#B0B0B0")) +
    labs(x = "", y = "NES (g3_vs_g2 ranking)", fill = "",
         title = "Celltype marker GSEA (g3_vs_g2)")
  save_fig(p, DIR_THIS, "celltype_gsea_g3v2_bar", width = 7, height = 5)
}

# ============================================================
# 11. summary
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("=== 03_1 celltype marker score 汇总 ===\n")
cat(strrep("=", 70), "\n", sep = "")
cat("TOP_N_SCORE =", TOP_N_SCORE, "（marker score 用）\n")
cat("TOP_N_MARKERS =", TOP_N_MARKERS, "（GSEA 用）\n")
cat("celltype 数:", length(cell_types), "\n")
cat("样本数:", nrow(merged), "\n\n")

cat("--- polr（单调）---\n")
print(polr_results, row.names = FALSE, digits = 4)

cat("\n--- poly3_F（非单调）---\n")
print(poly3_results, row.names = FALSE, digits = 4)

cat("\n--- 锯齿检查 ---\n")
print(zigzag_tbl, row.names = FALSE)

cat("\n--- GSEA poly3_F ---\n")
if (!is.null(gsea_poly3F)) {
  print(gsea_poly3F[, c("pathway","NES","pval","padj","size")],
        row.names = FALSE, digits = 4)
}

cat("\n--- GSEA g3_vs_g2 ---\n")
if (!is.null(gsea_g3v2)) {
  print(gsea_g3v2[, c("pathway","NES","pval","padj","size")],
        row.names = FALSE, digits = 4)
}

cat("\n--- 输出位置 ---\n")
cat("目录:", DIR_THIS, "\n")

log_info("=== 03_1_celltype_marker_score.r done ===")

