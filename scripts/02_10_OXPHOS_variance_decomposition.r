# ============================================================
# Output：results/02_10_OXPHOS_variance_decomposition/
# ============================================================

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr)
  library(ggplot2); library(patchwork)
})

has_vp <- requireNamespace("variancePartition", quietly = TRUE)

ROOT <- "E:/retina_tools/aging_rgc_axon/AMD"
OUT  <- file.path(ROOT, "results", "02_10_OXPHOS_variance_decomposition")
TAB  <- file.path(OUT, "tables"); FIG <- file.path(OUT, "figures")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

log_info("=== 02_10 (v4) ===")
set.seed(42)

# ------------------------------------------------------------
# ------------------------------------------------------------
expr <- readRDS(file.path(ROOT, "results/01_1b_batch_enhance/data/expr_wgcna.rds"))
clin <- readRDS(file.path(ROOT, "results/01_1b_batch_enhance/data/clinical_wgcna.rds"))
anno <- readRDS(file.path(ROOT, "results/01_1_prepare_bulk_data/data/gene_annotation.rds"))

if ("ensembl" %in% colnames(anno)) {
  anno$ensembl_use <- as.character(anno$ensembl)
} else if ("ensembl_no_ver" %in% colnames(anno)) {
  anno$ensembl_use <- as.character(anno$ensembl_no_ver)
} else {
  anno$ensembl_use <- rownames(anno)
}
anno$ensembl_use <- sub("\\..*$", "", anno$ensembl_use)
anno$symbol <- as.character(anno$symbol)
keep <- anno$ensembl_use %in% rownames(expr) & !is.na(anno$symbol) & anno$symbol != ""
keep[is.na(keep)] <- FALSE
map_df <- unique(anno[keep, c("symbol", "ensembl_use")])
colnames(map_df) <- c("symbol", "ensembl")

common <- intersect(colnames(expr), clin$sample_id)
expr <- expr[, common, drop = FALSE]
clin <- clin[match(common, clin$sample_id), ]

base_df <- data.frame(
  grade = suppressWarnings(as.numeric(as.character(clin$grade))),
  age   = suppressWarnings(as.numeric(clin$age)),
  sex   = as.factor(clin$sex),
  rin   = suppressWarnings(as.numeric(clin$rin)),
  pmi   = suppressWarnings(as.numeric(clin$pmi)),
  stringsAsFactors = FALSE
)
complete_idx <- complete.cases(base_df)
base_df <- base_df[complete_idx, , drop = FALSE]
expr    <- expr[, complete_idx, drop = FALSE]
log_info(paste0("完整样本：", nrow(base_df)))

# ------------------------------------------------------------
# ------------------------------------------------------------
oxphos_patterns <- list(
  Complex_I   = "^NDUF",
  Complex_II  = "^(SDH[ABCD]|SDHAF)",
  Complex_III = "^(UQCR|CYC1)",
  Complex_IV  = "^(COX[0-9]|COX4I|COX5[AB]|COX6[ABC]|COX7[ABC]|COX8A)",
  Complex_V   = "^(ATP5|ATP6|ATP6V)"
)
oxphos_by_complex <- lapply(oxphos_patterns, function(p) {
  grep(p, map_df$symbol, value = TRUE)
})
oxphos_all <- unique(unlist(oxphos_by_complex))

sig16 <- read_csv(file.path(ROOT,
  "results/02_1_differential_expression_main/tables/sig16_genes.csv"),
  show_col_types = FALSE)$gene
gy <- read_csv(file.path(ROOT,
  "results/07_5b_GSE231536_attachments/tables/greenyellow_eqtm_genes.csv"),
  show_col_types = FALSE)$gene

to_ens <- function(syms) {
  e <- map_df$ensembl[match(syms, map_df$symbol)]
  e[!is.na(e) & e %in% rownames(expr)]
}
oxphos_ens <- to_ens(oxphos_all)
sig16_ens  <- to_ens(sig16)
gy_ens     <- to_ens(gy)

# ------------------------------------------------------------
# ------------------------------------------------------------
calc_set_score <- function(ens, label) {
  if (length(ens) < 3) return(NULL)
  mat   <- expr[ens, , drop = FALSE]
  score <- rowMeans(scale(t(mat)), na.rm = TRUE)
  dat <- base_df; dat$score <- score
  dat <- dat[complete.cases(dat), , drop = FALSE]
  fit_full <- lm(score ~ grade + age + sex + rin + pmi, data = dat)
  fit_null <- lm(score ~ age + sex + rin + pmi, data = dat)
  r2_full <- summary(fit_full)$r.squared
  r2_null <- summary(fit_null)$r.squared
  data.frame(
    set = label,
    n_genes = length(ens),
    n_samples = nrow(dat),
    grade_partial_R2 = (r2_full - r2_null) / (1 - r2_null),
    grade_p = anova(fit_null, fit_full)$`Pr(>F)`[2],
    grade_coef = coef(fit_full)["grade"],
    row.names = NULL
  )
}

calc_single_gene_dir <- function(ens, label) {
  if (length(ens) < 1) return(NULL)
  res <- do.call(rbind, lapply(ens, function(e) {
    y <- as.numeric(expr[e, ])
    fit <- lm(y ~ grade + age + sex + rin + pmi, data = base_df)
    coefs <- summary(fit)$coefficients
    data.frame(
      set = label,
      symbol = map_df$symbol[match(e, map_df$ensembl)],
      ensembl = e,
      grade_coef = coefs["grade", "Estimate"],
      grade_p    = coefs["grade", "Pr(>|t|)"],
      stringsAsFactors = FALSE)
  }))
  res$direction <- ifelse(res$grade_coef > 0, "UP", "DOWN")
  res[order(res$grade_p), ]
}

dir_binom <- function(dir_vec) {
  n_up <- sum(dir_vec == "UP"); n_down <- sum(dir_vec == "DOWN")
  n <- n_up + n_down
  if (n == 0) return(c(n_up = NA, n_down = NA, n = 0, binom_p = NA))
  p <- binom.test(n_up, n, p = 0.5)$p.value
  c(n_up = n_up, n_down = n_down, n = n, binom_p = p)
}

# ------------------------------------------------------------
# ------------------------------------------------------------
log_info("--- 6. Gene-set scores ---")

sig16_dir_prelim <- calc_single_gene_dir(sig16_ens, "sig16")
strong_up <- sig16_dir_prelim %>%
  filter(grade_p < 0.05 & direction == "UP") %>%
  pull(symbol)
sig16_other <- setdiff(sig16, strong_up)
log_info(paste0("sig16 强驱动 UP：", paste(strong_up, collapse = ", ")))
log_info(paste0("sig16 其余 ", length(sig16_other), " 个：",
                paste(sig16_other, collapse = ", ")))

score_results <- rbind(
  calc_set_score(oxphos_ens,         "OXPHOS_129"),
  calc_set_score(sig16_ens,          "sig16_all16"),
  calc_set_score(sig16_ens[sig16 %in% strong_up], "sig16_strong_up"),
  calc_set_score(sig16_ens[sig16 %in% sig16_other], "sig16_other13"),
  calc_set_score(gy_ens,             "greenyellow_13")
)
cat("\n=== Gene-set 评分 ===\n"); print(score_results, row.names = FALSE)
write_csv(score_results, file.path(TAB, "gene_set_scores.csv"))

# ------------------------------------------------------------
# ------------------------------------------------------------
log_info("--- 7. Single-gene direction ---")
sig16_dir  <- calc_single_gene_dir(sig16_ens,  "sig16")
gy_dir     <- calc_single_gene_dir(gy_ens,     "greenyellow_13")
oxphos_dir <- calc_single_gene_dir(oxphos_ens, "OXPHOS_129")

sig16_binom  <- dir_binom(sig16_dir$direction)
gy_binom     <- dir_binom(gy_dir$direction)
oxphos_binom <- dir_binom(oxphos_dir$direction)

write_csv(sig16_dir,  file.path(TAB, "sig16_single_gene_direction.csv"))
write_csv(gy_dir,     file.path(TAB, "greenyellow_single_gene_direction.csv"))
write_csv(oxphos_dir, file.path(TAB, "oxphos_single_gene_direction.csv"))

# ------------------------------------------------------------
# ------------------------------------------------------------
log_info("--- 7b. OXPHOS by complex ---")

complex_results <- list()
for (cn in names(oxphos_by_complex)) {
  cn_sym <- oxphos_by_complex[[cn]]
  cn_ens <- to_ens(cn_sym)
  if (length(cn_ens) < 3) next
  dir_df <- calc_single_gene_dir(cn_ens, cn)
  binom  <- dir_binom(dir_df$direction)
  complex_results[[cn]] <- data.frame(
    complex  = cn,
    n_genes  = length(cn_ens),
    n_up     = binom["n_up"],
    n_down   = binom["n_down"],
    pct_down = binom["n_down"] / (binom["n_up"] + binom["n_down"]),
    binom_p  = binom["binom_p"],
    median_coef = median(dir_df$grade_coef, na.rm = TRUE),
    n_sig_dn = sum(dir_df$grade_p < 0.05 & dir_df$direction == "DOWN"),
    n_sig_up = sum(dir_df$grade_p < 0.05 & dir_df$direction == "UP"),
    row.names = NULL
  )
}
complex_df <- bind_rows(complex_results)
cat("\n=== OXPHOS 按复合物拆分 ===\n")
print(complex_df, row.names = FALSE)
write_csv(complex_df, file.path(TAB, "oxphos_by_complex.csv"))

complex_scores <- list()
for (cn in names(oxphos_by_complex)) {
  cn_ens <- to_ens(oxphos_by_complex[[cn]])
  if (length(cn_ens) < 3) next
  s <- calc_set_score(cn_ens, cn)
  if (!is.null(s)) complex_scores[[cn]] <- s
}
if (length(complex_scores) > 0) {
  complex_score_df <- bind_rows(complex_scores)
  cat("\n=== 复合物 gene-set 评分 ===\n")
  print(complex_score_df, row.names = FALSE)
  write_csv(complex_score_df, file.path(TAB, "oxphos_by_complex_scores.csv"))
}

# ------------------------------------------------------------
# ------------------------------------------------------------
log_info("--- 7c. Methylation-expression-grade three-way consistency ---")

cand_files <- c(
  file.path(ROOT, "results/07_5_GSE231536_methylation/tables/pathway_MGS1_vs_MGS3_eQTM.csv"),
  file.path(ROOT, "results/07_5_GSE231536_methylation/tables/MGS1_vs_MGS3_deltaM.csv"),
  file.path(ROOT, "results/07_5_GSE231536_methylation/tables/cpg_poly2_results.csv")
)
cand_files <- cand_files[file.exists(cand_files)]

threeway_df <- NULL
for (f in cand_files) {
  log_info(paste0("  尝试：", basename(f)))
  d <- tryCatch(read_csv(f, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(d)) next
  cn <- colnames(d)
  cat("    列名：", paste(cn, collapse = ", "), "\n")

  gene_col <- intersect(c("gene_sym", "gene", "symbol"), cn)[1]
  if (is.na(gene_col)) { log_info("    跳过：无基因列"); next }

  d_gy <- d %>% filter(!!sym(gene_col) %in% gy)
  cat("    greenyellow 匹配：", nrow(d_gy), " 行\n")
  if (nrow(d_gy) < 3) next

  delta_col <- intersect(c("delta_M", "deltaM", "delta_beta"), cn)[1]
  beta_col  <- intersect(c("beta", "eQTM_beta", "eqtm_beta"), cn)[1]
  cons_col  <- intersect(c("consistent", "concordant"), cn)[1]

  if (is.na(delta_col)) {
    log_info("    无 delta_M 列，跳过"); next
  }

  expr_dir_all <- rbind(sig16_dir, gy_dir, oxphos_dir) %>%
    filter(set == "greenyellow_13") %>%
    dplyr::select(symbol, grade_coef, grade_p, direction)

  merged <- d_gy %>%
    rename(symbol = !!sym(gene_col),
           delta_M = !!sym(delta_col)) %>%
    inner_join(expr_dir_all, by = "symbol")

  if (is.na(cons_col) && !is.na(beta_col)) {
    merged <- merged %>%
      mutate(consistent = sign(delta_M) == -sign(!!sym(beta_col)))
  } else if (!is.na(cons_col)) {
    merged <- merged %>% rename(consistent = !!sym(cons_col))
  } else {
    log_info("    无 eQTM 一致性列，跳过"); next
  }

  merged <- merged %>%
    mutate(
      meth_expr_concordant = consistent,
      meth_grade_concordant = sign(delta_M) == sign(grade_coef),
      threeway = meth_expr_concordant & meth_grade_concordant
    )

  threeway_df <- merged
  cat("\n=== 三级一致性明细（greenyellow 13）===\n")
  print(threeway_df[, c("symbol", "delta_M", "grade_coef",
                        "grade_p", "consistent",
                        "meth_grade_concordant", "threeway")],
        row.names = FALSE)
  write_csv(threeway_df, file.path(TAB, "threeway_consistency_greenyellow.csv"))

  # summary
  n_total <- nrow(threeway_df)
  n_me    <- sum(threeway_df$meth_expr_concordant, na.rm = TRUE)
  n_mg    <- sum(threeway_df$meth_grade_concordant, na.rm = TRUE)
  n_3way  <- sum(threeway_df$threeway, na.rm = TRUE)
  binom_me <- binom.test(n_me, n_total, 0.5)$p.value
  binom_3way <- binom.test(n_3way, n_total, 0.5)$p.value
  cat("\n三级一致性汇总：\n")
  cat("  甲基化-表达一致：", n_me, "/", n_total,
      " (binomial p =", signif(binom_me, 3), ")\n")
  cat("  甲基化-grade 同向：", n_mg, "/", n_total, "\n")
  cat("  三级一致：", n_3way, "/", n_total,
      " (binomial p =", signif(binom_3way, 3), ")\n")

  summary_df <- data.frame(
    set = "greenyellow_13",
    n_total = n_total,
    n_meth_expr = n_me,
    pct_meth_expr = n_me / n_total,
    binom_p_meth_expr = binom_me,
    n_meth_grade = n_mg,
    n_threeway = n_3way,
    pct_threeway = n_3way / n_total,
    binom_p_threeway = binom_3way,
    row.names = NULL
  )
  write_csv(summary_df, file.path(TAB, "threeway_consistency_summary.csv"))
  break  # 用第一个匹配的文件
}

if (is.null(threeway_df)) {
  log_info("  没有可用的 CpG 级甲基化文件，跳过三级一致性")
}

# ------------------------------------------------------------
# 8. variancePartition
# ------------------------------------------------------------
if (has_vp) {
  log_info("--- 8. variancePartition ---")
  suppressPackageStartupMessages(library(variancePartition))
  form <- ~ grade + age + sex + rin + pmi
  vp_summary_list <- list()
  for (item in list(
    list(ens = oxphos_ens, label = "OXPHOS_129"),
    list(ens = sig16_ens,  label = "sig16"),
    list(ens = gy_ens,     label = "greenyellow_13")
  )) {
    expr_sub <- expr[item$ens, , drop = FALSE]
    keep_g <- apply(expr_sub, 1, function(x) var(x, na.rm = TRUE) > 0)
    expr_sub <- expr_sub[keep_g, , drop = FALSE]
    if (nrow(expr_sub) < 3) next
    vp <- tryCatch(
      fitExtractVarPartModel(expr_sub, form, base_df, quiet = TRUE),
      error = function(e) NULL)
    if (is.null(vp)) next
    vp_df <- as.data.frame(vp)
    vp_df$ensembl <- rownames(vp_df)
    vp_df$symbol  <- map_df$symbol[match(vp_df$ensembl, map_df$ensembl)]
    vp_df$set     <- item$label
    vp_summary_list[[item$label]] <- vp_df
  }
  if (length(vp_summary_list) > 0) {
    vp_all <- bind_rows(vp_summary_list)
    write_csv(vp_all, file.path(TAB, "variance_partition_all_sets.csv"))
    median_by_set <- vp_all %>%
      group_by(set) %>%
      summarise(across(c(grade, age, sex, rin, pmi, Residuals),
                       ~ round(median(.x, na.rm = TRUE), 4)),
                .groups = "drop")
    cat("\n=== variancePartition ===\n")
    print(as.data.frame(median_by_set), row.names = FALSE)
  }
}

# ------------------------------------------------------------
# ------------------------------------------------------------
log_info("--- 9. Figure S14 ---")

plot_dir_bar <- function(dat, title) {
  ggplot(dat, aes(x = reorder(symbol, grade_coef),
                  y = grade_coef, fill = direction)) +
    geom_col() + geom_hline(yintercept = 0, color = "grey30") +
    scale_fill_manual(values = c("UP" = "firebrick", "DOWN" = "steelblue")) +
    labs(x = "", y = "Grade coefficient", title = title, fill = "Direction") +
    theme_classic(base_size = 9) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          plot.title = element_text(face = "bold"),
          legend.position = "top")
}

sig16_dir$group <- ifelse(sig16_dir$symbol %in% strong_up,
                          "Strong UP", "Other")
pA <- ggplot(sig16_dir, aes(x = reorder(symbol, grade_coef),
                             y = grade_coef, fill = direction)) +
  geom_col() +
  geom_hline(yintercept = 0, color = "grey30") +
  geom_text(aes(label = ifelse(grade_p < 0.05, "*", "")),
            vjust = ifelse(sig16_dir$grade_coef > 0, -0.5, 1.2), size = 3) +
  scale_fill_manual(values = c("UP" = "firebrick", "DOWN" = "steelblue")) +
  labs(x = "", y = "Grade coefficient",
       title = "A  16-gene signature (3 strong UP *)",
       fill = "Direction") +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        plot.title = element_text(face = "bold"),
        legend.position = "top")

# B: greenyellow 13
pB <- plot_dir_bar(gy_dir, "B  greenyellow 13 (12 DOWN)")

pC <- ggplot(complex_df, aes(x = factor(complex, levels = complex_df$complex),
                              y = pct_down, fill = binom_p < 0.05)) +
  geom_col() +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  geom_text(aes(label = paste0(n_down, "/", n_genes)),
            vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("FALSE" = "grey60", "TRUE" = "firebrick"),
                    labels = c("P ≥ 0.05", "P < 0.05"),
                    name = "Binomial") +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "", y = "Proportion DOWN",
       title = "C  OXPHOS by complex (94/129 DOWN)") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        plot.title = element_text(face = "bold"))

score_plot_df <- score_results %>%
  mutate(sig_label = ifelse(grade_p < 0.05, "*", ""))
pD <- ggplot(score_plot_df,
             aes(x = reorder(set, grade_coef), y = grade_coef,
                 fill = grade_p < 0.05)) +
  geom_col() +
  geom_hline(yintercept = 0, color = "grey30") +
  geom_text(aes(label = sig_label),
            vjust = ifelse(score_plot_df$grade_coef > 0, -0.5, 1.2), size = 4) +
  scale_fill_manual(values = c("FALSE" = "grey60", "TRUE" = "firebrick"),
                    guide = "none") +
  coord_flip() +
  labs(x = "", y = "Grade coefficient (gene-set score)",
       title = "D  Gene-set scores (* p<0.05)") +
  theme_classic(base_size = 10) +
  theme(plot.title = element_text(face = "bold"))

fig_s14 <- (pA | pB) / (pC | pD) +
  plot_annotation(
    title = "G2–G3 transcriptional programs",
    subtitle = "16-gene (heterogeneous) | OXPHOS by complex (coordinated down) | greenyellow 13",
    theme = theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
                  plot.subtitle = element_text(size = 10, hjust = 0.5,
                                                color = "grey30")))
ggsave(file.path(FIG, "FigureS14_variance_direction.pdf"), fig_s14,
       width = 12, height = 9, device = cairo_pdf)
ggsave(file.path(FIG, "FigureS14_variance_direction.png"), fig_s14,
       width = 12, height = 9, dpi = 300)

if (!is.null(threeway_df) && nrow(threeway_df) > 0) {
  tw_long <- threeway_df %>%
    dplyr::select(symbol, delta_M, grade_coef) %>%
    pivot_longer(-symbol, names_to = "measure", values_to = "value") %>%
    mutate(sign_val = ifelse(value > 0, "UP", "DOWN"))
  pE <- ggplot(tw_long, aes(x = measure, y = symbol, fill = sign_val)) +
    geom_tile(color = "white") +
    scale_fill_manual(values = c("UP" = "firebrick", "DOWN" = "steelblue")) +
    labs(x = "", y = "", fill = "Direction",
         title = "E  Three-way consistency (greenyellow)") +
    theme_minimal(base_size = 10) +
    theme(plot.title = element_text(face = "bold"))
  ggsave(file.path(FIG, "FigureS14E_threeway_heatmap.pdf"), pE,
         width = 6, height = 4, device = cairo_pdf)
  ggsave(file.path(FIG, "FigureS14E_threeway_heatmap.png"), pE,
         width = 6, height = 4, dpi = 300)
}

# ------------------------------------------------------------
# 10. summary
# ------------------------------------------------------------
cat("\n", strrep("=", 72), "\n", sep = "")
cat("=== 02_10 summary (v4) ===\n")
cat(strrep("=", 72), "\n", sep = "")

cat("\n[1] Gene-set scores（含 sig16 拆分）\n")
print(score_results, row.names = FALSE)

cat("\n[1b] sig16 拆分：强 UP vs 其余\n")
cat("  强 UP（", length(strong_up), "）：",
    paste(strong_up, collapse = ", "), "\n", sep = "")
cat("  其余（", length(sig16_other), "）：",
    paste(sig16_other, collapse = ", "), "\n", sep = "")

cat("\n[2] 单基因方向分布\n")
dir_summary <- data.frame(
  set = c("sig16", "greenyellow_13", "OXPHOS_129"),
  n = c(nrow(sig16_dir), nrow(gy_dir), nrow(oxphos_dir)),
  n_up = c(sig16_binom["n_up"], gy_binom["n_up"], oxphos_binom["n_up"]),
  n_down = c(sig16_binom["n_down"], gy_binom["n_down"], oxphos_binom["n_down"]),
  binom_p = c(sig16_binom["binom_p"], gy_binom["binom_p"], oxphos_binom["binom_p"]),
  row.names = NULL
)
print(dir_summary, row.names = FALSE)
write_csv(dir_summary, file.path(TAB, "direction_summary.csv"))

cat("\n[2b] OXPHOS 按复合物\n")
print(complex_df, row.names = FALSE)

if (!is.null(threeway_df) && nrow(threeway_df) > 0) {
  cat("\n[3] 三级一致性（greenyellow 13）\n")
  cat("  甲基化-表达一致：", n_me, "/", n_total, "\n", sep = "")
  cat("  甲基化-grade 同向：", n_mg, "/", n_total, "\n", sep = "")
  cat("  三级一致：", n_3way, "/", n_total,
      " (p = ", signif(binom_3way, 3), ")\n", sep = "")
}

cat("\n[4] 结论\n")
cat("  - sig16: 11 UP / 5 DOWN, 但 3 个强驱动 UP (BTBD19, UGGT1, S100B)\n")
cat("  - OXPHOS 129: 94/129 DOWN, 通路协调下调, 最强复合物 = ",
    complex_df$complex[which.min(complex_df$binom_p)], "\n", sep = "")
cat("  - greenyellow 13: 12/13 DOWN, 是 OXPHOS 下调的极端子集\n")
cat("  - 两个信号不重叠（16 ∩ greenyellow = 0）\n")

cat("\n[Files saved]\n"); print(list.files(TAB)); print(list.files(FIG))
log_info("=== 02_10 done ===")

