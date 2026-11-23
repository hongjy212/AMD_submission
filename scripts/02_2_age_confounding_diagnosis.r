# ============================================================
# 02_2_age_confounding_diagnosis.r
# Input：
#   results/01_1_prepare_bulk_data/data/{expr_corrected,clinical,gene_annotation,gene_sets}.rds
#   results/02_1_differential_expression_main/tables/sig16_genes.csv
#   tables/age_by_grade.csv
#   tables/age_pairwise_diff.csv
#   tables/multi_none.csv
#   tables/multi_spline.csv
#   tables/sig_overlap.csv
#   tables/enrichment_g4v1_ageNone.csv
#   tables/enrichment_g3v2_ageSpline.csv
#   tables/sig16_layered.csv
#   tables/unique3_grade_means.csv
#   tables/overlap13_zigzag.csv
#   figures/age_by_grade.pdf/.png
#   figures/unique3_grade_curves.pdf/.png
#   figures/overlap13_curves.pdf/.png
# ============================================================

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(limma); library(splines); library(ggplot2)
  library(tidyr); library(dplyr)
})

DIR_THIS <- setup_output_dir("02_2_age_confounding_diagnosis")
log_info("=== 02_2_age_confounding_diagnosis.r start ===")

# ============================================================
# ============================================================
log_info("--- 1. 载入数据 ---")

expr <- read_upstream("01_1_prepare_bulk_data", "data", "expr_corrected.rds")
clin <- read_upstream("01_1_prepare_bulk_data", "data", "clinical.rds")
anno <- read_upstream("01_1_prepare_bulk_data", "data", "gene_annotation.rds")
GENE_SETS <- read_upstream("01_2_prepare_gene_sets", "data", "gene_sets.rds")

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

sig16_file <- file.path(DIR_RESULT,
                         "02_1_differential_expression_main",
                         "tables", "sig16_genes.csv")
if (!file.exists(sig16_file)) {
  stop("请先运行 02_1_differential_expression_main.r", call. = FALSE)
}
sig16 <- read.csv(sig16_file, stringsAsFactors = FALSE)$gene
log_info("读入 16 基因: {paste(sig16, collapse=', ')}")

# ============================================================
# 2. age vs grade diagnosis
# ============================================================
log_info("--- 2. age vs grade 诊断 ---")

age_tbl <- aggregate(age_n ~ grade, data = clin,
                     FUN = function(x) c(
                       n = length(x),
                       median = median(x, na.rm = TRUE),
                       mean = mean(x, na.rm = TRUE),
                       sd = sd(x, na.rm = TRUE)))
age_df <- data.frame(
  grade = age_tbl$grade,
  n     = age_tbl$age_n[, "n"],
  median = round(age_tbl$age_n[, "median"], 1),
  mean   = round(age_tbl$age_n[, "mean"], 1),
  sd     = round(age_tbl$age_n[, "sd"], 1)
)
save_table(age_df, DIR_THIS, "age_by_grade.csv")
print(age_df, row.names = FALSE)

aov_fit <- aov(age_n ~ factor(grade), data = clin)
aov_p <- summary(aov_fit)[[1]][["Pr(>F)"]][1]
sp_test <- cor.test(clin$grade, clin$age_n, method = "spearman")
log_info("ANOVA p = {signif(aov_p, 3)}")
log_info("Spearman rho = {round(sp_test$estimate, 3)}, p = {signif(sp_test$p.value, 3)}")

if (abs(sp_test$estimate) > 0.3 && sp_test$p.value < 0.01) {
  log_warn("|rho|>0.3 且 p<0.01 → age 是强混杂")
} else if (abs(sp_test$estimate) > 0.15 && sp_test$p.value < 0.05) {
  log_info("|rho| 中等 → age 是弱混杂")
} else {
  log_info("|rho|<0.15 → age 与 grade 弱关联")
}

grade_pairs <- list(c(1,2), c(2,3), c(3,4), c(1,4))
pair_tbl <- data.frame()
for (p in grade_pairs) {
  a <- clin$age_n[clin$grade == p[1]]
  b <- clin$age_n[clin$grade == p[2]]
  tt <- t.test(b, a)
  pair_tbl <- rbind(pair_tbl, data.frame(
    pair = paste0("g", p[2], "_vs_g", p[1]),
    n1 = length(a), n2 = length(b),
    median_diff = round(median(b, na.rm = TRUE) - median(a, na.rm = TRUE), 1),
    t_p = signif(tt$p.value, 3),
    stringsAsFactors = FALSE))
}
save_table(pair_tbl, DIR_THIS, "age_pairwise_diff.csv")
print(pair_tbl, row.names = FALSE)

# ============================================================
# ============================================================
log_info("--- 3. age=none vs age=spline ---")

run_limma_multi <- function(expr, clin, age_mode) {
  covar <- data.frame(sex_M = as.integer(clin$sex == "Male"),
                       rin = clin$rin_n, pmi = clin$pmi_n)
  if (age_mode == "linear") {
    covar$age <- clin$age_n
  } else if (age_mode == "spline") {
    sp <- splines::ns(clin$age_n, df = 3)
    colnames(sp) <- c("age_s1", "age_s2", "age_s3")
    covar <- cbind(covar, sp)
  }
  rownames(covar) <- rownames(clin)

  design_df <- data.frame(grade = factor(clin$grade), covar)
  rownames(design_df) <- rownames(clin)
  design_df <- design_df[complete.cases(design_df), , drop = FALSE]
  expr_sub <- expr[, rownames(design_df), drop = FALSE]

  fml <- as.formula(paste("~ 0 + grade +",
                          paste(colnames(covar), collapse = " + ")))
  design <- model.matrix(fml, data = design_df)
  colnames(design)[grep("^grade", colnames(design))] <- paste0("g", 1:4)

  cm <- limma::makeContrasts(
    g2_vs_g1 = g2 - g1, g3_vs_g2 = g3 - g2, g4_vs_g3 = g4 - g3,
    g4_vs_g1 = g4 - g1,
    anyAMD_vs_g1 = (g2 + g3 + g4)/3 - g1,
    late_vs_early = (g3 + g4)/2 - (g1 + g2)/2,
    g4_vs_rest = g4 - (g1 + g2 + g3)/3,
    levels = design)

  fit <- eBayes(contrasts.fit(lmFit(expr_sub, design), cm))
  out <- data.frame()
  for (ct in colnames(cm)) {
    tt <- topTable(fit, coef = ct, number = Inf, sort.by = "none")
    out <- rbind(out, data.frame(
      gene = rownames(tt), contrast = ct,
      logFC = tt$logFC, p = tt$P.Value, FDR = tt$adj.P.Val,
      stringsAsFactors = FALSE))
  }
  out
}

res_none <- run_limma_multi(expr, clin, "none")
res_spline <- run_limma_multi(expr, clin, "spline")
save_table(res_none, DIR_THIS, "multi_none.csv")
save_table(res_spline, DIR_THIS, "multi_spline.csv")

# ============================================================
# ============================================================
log_info("--- 4. 显著基因对比 ---")

sig_none_g4v1 <- res_none$gene[res_none$contrast == "g4_vs_g1" &
                                !is.na(res_none$FDR) & res_none$FDR < 0.05]
sig_spline_g3v2 <- res_spline$gene[res_spline$contrast == "g3_vs_g2" &
                                    !is.na(res_spline$FDR) & res_spline$FDR < 0.05]

log_info("age=none   g4_vs_g1: {length(sig_none_g4v1)}")
log_info("age=spline g3_vs_g2: {length(sig_spline_g3v2)}")

overlap <- intersect(sig_none_g4v1, sig_spline_g3v2)
overlap_tbl <- data.frame(
  set1 = "g4_vs_g1_ageNone", set2 = "g3_vs_g2_ageSpline",
  n_set1 = length(sig_none_g4v1), n_set2 = length(sig_spline_g3v2),
  n_overlap = length(overlap),
  jaccard = round(length(overlap) / length(union(sig_none_g4v1, sig_spline_g3v2)), 4))
save_table(overlap_tbl, DIR_THIS, "sig_overlap.csv")
print(overlap_tbl, row.names = FALSE)

# ============================================================
# ============================================================
log_info("--- 5. 富集对比 ---")
UNIVERSE <- rownames(expr)

run_ora <- function(genes, universe, gene_sets, min_hit = 2, min_panel = 5) {
  genes <- intersect(genes, universe)
  n_sig <- length(genes)
  if (n_sig < 3) return(NULL)
  out <- data.frame()
  for (gs in names(gene_sets)) {
    panel <- intersect(gene_sets[[gs]]$genes, universe)
    if (length(panel) < min_panel) next
    hits <- intersect(genes, panel)
    if (length(hits) < min_hit) next
    p <- phyper(length(hits)-1, length(panel),
                length(universe)-length(panel), n_sig, lower.tail = FALSE)
    expected <- n_sig * length(panel) / length(universe)
    out <- rbind(out, data.frame(
      geneset = gs, group = gene_sets[[gs]]$group,
      n_hit = length(hits), n_sig = n_sig,
      fold = round(length(hits)/max(expected, 1e-6), 2),
      p = p, hits = paste(hits, collapse = ","),
      stringsAsFactors = FALSE))
  }
  if (nrow(out) == 0) return(NULL)
  out$FDR <- p.adjust(out$p, "BH")
  out[order(out$p), ]
}

ora_none <- run_ora(sig_none_g4v1, UNIVERSE, GENE_SETS)
ora_spline <- run_ora(sig_spline_g3v2, UNIVERSE, GENE_SETS)

if (!is.null(ora_none)) {
  save_table(ora_none, DIR_THIS, "enrichment_g4v1_ageNone.csv")
  cat("\n=== age=None g4_vs_g1 富集 ===\n")
  print(head(ora_none[, c("geneset","group","n_hit","fold","FDR")], 10),
        row.names = FALSE)
}
if (!is.null(ora_spline)) {
  save_table(ora_spline, DIR_THIS, "enrichment_g3v2_ageSpline.csv")
  cat("\n=== age=Spline g3_vs_g2 富集 ===\n")
  print(head(ora_spline[, c("geneset","group","n_hit","fold","FDR")], 10),
        row.names = FALSE)
}

# ============================================================
# ============================================================
log_info("--- 6. 16 基因内部分层 ---")

overlap13 <- intersect(sig16, sig_spline_g3v2)
unique3   <- setdiff(sig16, sig_spline_g3v2)

log_info("16 ∩ 86 = {length(overlap13)}")
log_info("poly3_F 独有 = {length(unique3)}")

cat("\n=== 16 基因分层 ===\n")
cat("与 limma g3_vs_g2 重叠 (", length(overlap13), "):\n", sep = "")
cat(paste(overlap13, collapse = ", "), "\n\n")
cat("poly3_F 独有 (", length(unique3), "):\n", sep = "")
cat(paste(unique3, collapse = ", "), "\n\n")

save_table(data.frame(
  gene = sig16,
  category = ifelse(sig16 %in% overlap13,
                    "limma_overlap", "poly3F_unique"),
  stringsAsFactors = FALSE),
  DIR_THIS, "sig16_layered.csv")

if (length(unique3) >= 2) {
  expr_grade_u3 <- sapply(1:4, function(g) {
    rowMeans(expr[unique3, clin$grade == g, drop = FALSE])
  })
  rownames(expr_grade_u3) <- unique3
  colnames(expr_grade_u3) <- paste0("g", 1:4)
  save_table(as.data.frame(expr_grade_u3),
             DIR_THIS, "unique3_grade_means.csv")

  cat("=== poly3_F 独有基因的 grade 均值 ===\n")
  print(round(expr_grade_u3, 3))

  expr_long_u3 <- data.frame(
    gene = rep(unique3, each = 4),
    grade = rep(1:4, times = length(unique3)),
    expr = as.numeric(expr_grade_u3))

  colors_u3 <- c("#4C72B0", "#C44E52", "#55A868",
                 "#8172B2", "#FDB462")[1:length(unique3)]
  p_u3 <- ggplot(expr_long_u3,
                  aes(x = grade, y = expr, color = gene)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 3.5) +
    theme_bw(base_size = 13) +
    scale_x_continuous(breaks = 1:4, labels = paste0("G", 1:4)) +
    scale_color_manual(values = colors_u3) +
    labs(x = "AMD grade", y = "Mean expression (log2 CPM)",
         title = "Non-monotonic pattern in poly3F-unique genes",
         subtitle = "1 < 2 > 3 < 4",
         color = "")
  save_fig(p_u3, DIR_THIS, "unique3_grade_curves", width = 7, height = 4.5)
}

if (length(overlap13) >= 3) {
  expr_grade_13 <- sapply(1:4, function(g) {
    rowMeans(expr[overlap13, clin$grade == g, drop = FALSE])
  })
  rownames(expr_grade_13) <- overlap13
  colnames(expr_grade_13) <- paste0("g", 1:4)

  zigzag_flag <- expr_grade_13[, "g2"] > expr_grade_13[, "g3"]

  zigzag_tbl <- data.frame(
    gene = overlap13,
    g1 = round(expr_grade_13[, "g1"], 3),
    g2 = round(expr_grade_13[, "g2"], 3),
    g3 = round(expr_grade_13[, "g3"], 3),
    g4 = round(expr_grade_13[, "g4"], 3),
    is_zigzag = zigzag_flag,
    stringsAsFactors = FALSE)
  save_table(zigzag_tbl, DIR_THIS, "overlap13_zigzag.csv")

  cat("\n=== 13 个重叠基因的锯齿检查 ===\n")
  cat("g2 > g3 的基因数:", sum(zigzag_flag), "/",
      length(overlap13), "\n")
  print(zigzag_tbl, row.names = FALSE)

  expr_long_13 <- data.frame(
    gene = rep(overlap13, each = 4),
    grade = rep(1:4, times = length(overlap13)),
    expr = as.numeric(expr_grade_13))
  expr_long_13$zigzag <- rep(zigzag_flag, each = 4)

  p_13 <- ggplot(expr_long_13,
                  aes(x = grade, y = expr,
                      color = zigzag, group = gene)) +
    geom_line(linewidth = 0.8, alpha = 0.7) +
    geom_point(size = 1.5) +
    theme_bw(base_size = 12) +
    scale_x_continuous(breaks = 1:4, labels = paste0("G", 1:4)) +
    scale_color_manual(
      values = c("TRUE" = "#C44E52", "FALSE" = "#B0B0B0"),
      labels = c("TRUE" = "g2>g3 (zigzag)",
                 "FALSE" = "g2<g3")) +
    labs(x = "AMD grade", y = "Mean expression",
         title = "13 limma-overlap genes: zigzag pattern check",
         color = "")
  save_fig(p_13, DIR_THIS, "overlap13_curves", width = 7, height = 4.5)
}

# ============================================================
# 7. figure（age_by_grade）
# ============================================================
log_info("--- 7. 图 ---")

p1 <- ggplot(clin, aes(x = factor(grade), y = age_n)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, fill = "#B0B0B0") +
  geom_jitter(width = 0.15, alpha = 0.3, size = 0.8) +
  stat_summary(fun = median, geom = "line", aes(group = 1),
               color = "red", size = 0.8) +
  stat_summary(fun = median, geom = "point", color = "red", size = 2) +
  theme_bw(base_size = 12) +
  labs(x = "AMD grade", y = "Age (years)",
       title = sprintf("Age by grade (Spearman rho=%.3f, p=%.2g)",
                       sp_test$estimate, sp_test$p.value),
       subtitle = sprintf("ANOVA p=%.3g", aov_p))
save_fig(p1, DIR_THIS, "age_by_grade", width = 5, height = 4)

# ============================================================
# 8. summary
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("=== 02_2 汇总 ===\n")
cat(strrep("=", 70), "\n", sep = "")

cat("\n--- 1. age vs grade ---\n")
cat(sprintf("Spearman rho = %.3f, p = %.3g\n",
            sp_test$estimate, sp_test$p.value))
cat(sprintf("ANOVA p = %.3g\n", aov_p))
cat("\n各 grade 年龄中位数:\n")
print(age_df, row.names = FALSE)

cat("\n--- 2. age=none vs age=spline ---\n")
cat(sprintf("age=None   g4_vs_g1: %d 基因\n", length(sig_none_g4v1)))
cat(sprintf("age=Spline g3_vs_g2: %d 基因\n", length(sig_spline_g3v2)))
cat(sprintf("Jaccard = %.4f\n", overlap_tbl$jaccard))

cat("\n--- 3. 富集对比 ---\n")
if (!is.null(ora_none)) {
  cat("age=None g4_vs_g1 top5:\n")
  print(head(ora_none[, c("geneset","n_hit","fold","FDR")], 5),
        row.names = FALSE)
}
if (!is.null(ora_spline)) {
  cat("\nage=Spline g3_vs_g2 top5:\n")
  print(head(ora_spline[, c("geneset","n_hit","fold","FDR")], 5),
        row.names = FALSE)
}

cat("\n--- 4. 16 基因内部分层 ---\n")
cat(sprintf("limma 重叠: %d\n", length(overlap13)))
cat(sprintf("poly3F 独有: %d\n", length(unique3)))
cat("poly3F 独有基因:", paste(unique3, collapse = ", "), "\n")

if (exists("zigzag_tbl")) {
  cat(sprintf("\n13 重叠基因中 g2>g3 (锯齿) 的: %d / %d\n",
              sum(zigzag_tbl$is_zigzag), nrow(zigzag_tbl)))
}

cat("\n--- 输出位置 ---\n")
cat("目录:", DIR_THIS, "\n")
cat("  tables/age_by_grade.csv\n")
cat("  tables/age_pairwise_diff.csv\n")
cat("  tables/multi_none.csv\n")
cat("  tables/multi_spline.csv\n")
cat("  tables/sig_overlap.csv\n")
cat("  tables/enrichment_g4v1_ageNone.csv\n")
cat("  tables/enrichment_g3v2_ageSpline.csv\n")
cat("  tables/sig16_layered.csv\n")
cat("  tables/unique3_grade_means.csv\n")
cat("  tables/overlap13_zigzag.csv\n")
cat("  figures/age_by_grade.pdf/.png\n")
cat("  figures/unique3_grade_curves.pdf/.png\n")
cat("  figures/overlap13_curves.pdf/.png\n")

log_info("=== 02_2_age_confounding_diagnosis.r done ===")

