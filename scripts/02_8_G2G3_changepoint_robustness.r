# ============================================================
# Output：results/02_8_G2G3_changepoint_robustness/
# ============================================================

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr)
  library(splines); library(ggplot2); library(patchwork)
})

has_segmented <- requireNamespace("segmented", quietly = TRUE)
has_mcp       <- requireNamespace("mcp",       quietly = TRUE)
has_bcp       <- requireNamespace("bcp",       quietly = TRUE)
if (has_segmented) suppressPackageStartupMessages(library(segmented))
if (has_mcp)       suppressPackageStartupMessages(library(mcp))
if (has_bcp)       suppressPackageStartupMessages(library(bcp))

ROOT <- "E:/retina_tools/aging_rgc_axon/AMD"
OUT  <- file.path(ROOT, "results", "02_8_G2G3_changepoint_robustness")
TAB  <- file.path(OUT, "tables"); FIG <- file.path(OUT, "figures")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

log_info("=== 02_8 G2-G3 changepoint robustness (v3) ===")
log_info(paste0("Packages: segmented=", has_segmented,
                ", mcp=", has_mcp, ", bcp=", has_bcp))
set.seed(42)

# ------------------------------------------------------------
# ------------------------------------------------------------
expr_file    <- file.path(ROOT, "results", "01_1b_batch_enhance", "data", "expr_wgcna.rds")
clin_file    <- file.path(ROOT, "results", "01_1b_batch_enhance", "data", "clinical_wgcna.rds")
anno_file    <- file.path(ROOT, "results", "01_1_prepare_bulk_data", "data", "gene_annotation.rds")
sig16_file   <- file.path(ROOT, "results", "02_1_differential_expression_main", "tables", "sig16_genes.csv")
gene_p3_file <- file.path(ROOT, "results", "06_1_TRADE_vs_adjacent", "tables", "poly3F_gene_level.csv")

expr <- readRDS(expr_file); clin <- readRDS(clin_file); anno <- readRDS(anno_file)

# ------------------------------------------------------------
# 2. samplesalign
# ------------------------------------------------------------
common <- intersect(colnames(expr), clin$sample_id)
expr <- expr[, common, drop = FALSE]
clin <- clin[match(common, clin$sample_id), ]
rownames(clin) <- clin$sample_id
log_info(paste0("对齐后样本数：", length(common)))

# ------------------------------------------------------------
# 3. Ensembl → Symbol mapping
# ------------------------------------------------------------
if ("ensembl" %in% colnames(anno)) {
  anno$ensembl_use <- as.character(anno$ensembl)
} else if ("ensembl_no_ver" %in% colnames(anno)) {
  anno$ensembl_use <- as.character(anno$ensembl_no_ver)
} else {
  anno$ensembl_use <- rownames(anno)
}
anno$ensembl_use <- sub("\\..*$", "", anno$ensembl_use)
anno$symbol      <- as.character(anno$symbol)

keep <- anno$ensembl_use %in% rownames(expr) &
        !is.na(anno$symbol) & anno$symbol != ""
keep[is.na(keep)] <- FALSE
map_df <- unique(anno[keep, c("symbol", "ensembl_use")])
colnames(map_df) <- c("symbol", "ensembl")

# ------------------------------------------------------------
# ------------------------------------------------------------
sig16_sym <- read_csv(sig16_file, show_col_types = FALSE)$gene
sig16_map <- merge(data.frame(symbol = sig16_sym, stringsAsFactors = FALSE),
                   map_df, by = "symbol", all.x = FALSE)
sig16_map <- sig16_map[!duplicated(sig16_map$symbol), ]
sig16_ens <- sig16_map$ensembl
log_info(paste0("16 基因中可评估：", length(sig16_ens), " / ", length(sig16_sym)))
missing_sym <- setdiff(sig16_sym, sig16_map$symbol)
if (length(missing_sym) > 0) {
  log_info(paste0("未映射基因：", paste(missing_sym, collapse = ", ")))
}

gene2ens <- setNames(sig16_map$ensembl, sig16_map$symbol)
ens2gene <- setNames(sig16_map$symbol, sig16_map$ensembl)

# ------------------------------------------------------------
# ------------------------------------------------------------
log_info(paste0("clin$grade 类型：", class(clin$grade)[1]))
log_info(paste0("clin$grade 唯一值：", paste(sort(unique(clin$grade)), collapse = ", ")))

grade_num <- suppressWarnings(as.numeric(as.character(clin$grade)))
if (any(is.na(grade_num) & !is.na(clin$grade))) {
  stop("grade 无法转换为数值，请检查原始编码（可能是 G1/G2 字符串）")
}

base_df <- data.frame(
  grade = grade_num,
  age   = suppressWarnings(as.numeric(clin$age)),
  sex   = as.factor(clin$sex),
  rin   = suppressWarnings(as.numeric(clin$rin)),
  pmi   = suppressWarnings(as.numeric(clin$pmi)),
  stringsAsFactors = FALSE
)

na_summary <- sapply(base_df, function(x) sum(is.na(x)))
log_info(paste0("协变量 NA 数：", paste(names(na_summary), na_summary,
                                       sep = "=", collapse = ", ")))

complete_idx <- complete.cases(base_df)
log_info(paste0("完整案例样本数：", sum(complete_idx), " / ", nrow(base_df)))
if (sum(complete_idx) < nrow(base_df)) {
  base_df <- base_df[complete_idx, , drop = FALSE]
  expr    <- expr[, complete_idx, drop = FALSE]
}

grade_tab <- table(base_df$grade)
log_info(paste0("grade 分布：",
                paste(names(grade_tab), as.integer(grade_tab),
                      sep = "=", collapse = ", ")))

get_expr <- function(ens) as.numeric(expr[ens, ])

fit_p3_safe <- function(y, dat) {
  tryCatch({
    fit_full <- lm(y ~ poly(grade, 3) + ns(age, df = 3) + sex + rin + pmi, data = dat)
    fit_null <- lm(y ~ ns(age, df = 3) + sex + rin + pmi, data = dat)
    aov_res  <- anova(fit_null, fit_full)
    p_poly   <- if (nrow(aov_res) >= 2) aov_res$`Pr(>F)`[2] else NA_real_
    list(p3_p = p_poly, fit = fit_full, df_resid = fit_full$df.residual)
  }, error = function(e) {
    list(p3_p = NA_real_, fit = NULL, df_resid = NA_integer_)
  })
}

# ============================================================
# ============================================================
log_info("--- 6. Natural spline model ---")
spline_results <- list()
for (ens in sig16_ens) {
  g <- ens2gene[ens]
  dat <- base_df; dat$y <- get_expr(ens)
  res <- tryCatch({
    fit_full <- lm(y ~ ns(grade, df = 3) + ns(age, df = 3) + sex + rin + pmi, data = dat)
    fit_null <- lm(y ~ ns(age, df = 3) + sex + rin + pmi, data = dat)
    aov_res <- anova(fit_null, fit_full)
    spline_F <- aov_res$F[2]; spline_p <- aov_res$`Pr(>F)`[2]
    r2_full <- summary(fit_full)$r.squared
    r2_null <- summary(fit_null)$r.squared
    partial_r2 <- (r2_full - r2_null) / (1 - r2_null)
    data.frame(gene = g, ensembl = ens,
               spline_F = spline_F, spline_p = spline_p,
               partial_R2 = partial_r2, stringsAsFactors = FALSE)
  }, error = function(e) {
    data.frame(gene = g, ensembl = ens, spline_F = NA_real_,
               spline_p = NA_real_, partial_R2 = NA_real_,
               stringsAsFactors = FALSE)
  })
  spline_results[[ens]] <- res
}
spline_df <- bind_rows(spline_results) %>%
  mutate(spline_FDR = p.adjust(spline_p, method = "BH")) %>%
  arrange(spline_p)
write_csv(spline_df, file.path(TAB, "spline_ns_results.csv"))
log_info(paste0("Spline: ", sum(spline_df$spline_FDR < 0.05, na.rm = TRUE),
                " / ", nrow(spline_df), " genes FDR < 0.05"))

# ============================================================
# ============================================================
log_info("--- 7. Segmented regression ---")
seg_df <- NULL
if (has_segmented) {
  seg_results <- list()
  for (ens in sig16_ens) {
    g <- ens2gene[ens]; dat <- base_df; dat$y <- get_expr(ens)
    res <- tryCatch({
      fit_lm <- lm(y ~ grade + age + sex + rin + pmi, data = dat)
      seg <- segmented(fit_lm, seg.Z = ~grade, psi = list(grade = 2.5),
                       control = seg.control(it.max = 50, fix.npsi = TRUE))
      bp <- seg$psi[1, "Est."]; bp_se <- seg$psi[1, "St.Err"]
      davies_p <- tryCatch(davies.test(seg, seg.Z = ~grade,
                                       alternative = "two.sided")$p.value,
                           error = function(e) NA_real_)
      data.frame(gene = g, ensembl = ens, breakpoint = bp,
                 bp_se = bp_se, davies_p = davies_p, stringsAsFactors = FALSE)
    }, error = function(e) NULL)
    if (!is.null(res)) seg_results[[ens]] <- res
  }
  if (length(seg_results) > 0) {
    seg_df <- bind_rows(seg_results) %>%
      mutate(davies_FDR = p.adjust(davies_p, method = "BH")) %>%
      arrange(davies_p)
    write_csv(seg_df, file.path(TAB, "segmented_breakpoint_results.csv"))
    log_info(paste0("Segmented: median breakpoint = ",
                    round(median(seg_df$breakpoint, na.rm = TRUE), 2),
                    " (n = ", nrow(seg_df), ")"))
    log_info(paste0("Segmented: Davies FDR < 0.05 = ",
                    sum(seg_df$davies_FDR < 0.05, na.rm = TRUE),
                    " / ", nrow(seg_df),
                    " → 若为 0 则结果不可用（grade 只有 4 个离散值）"))
  } else {
    log_info("Segmented: 无基因成功拟合")
  }
}

# ============================================================
# ============================================================
log_info("--- 8. Bayesian changepoint ---")

mcp_df <- NULL
if (has_mcp) {
  log_info("  尝试 mcp ...")
  mcp_results <- list()
  mcp_errors <- 0
  for (ens in sig16_ens) {
    g <- ens2gene[ens]
    dat <- data.frame(y = get_expr(ens), grade = base_df$grade)
    fit_mcp <- tryCatch({
      mcp(list(y ~ 1 + grade, ~ 1 + grade), data = dat, par_x = "grade",
          iter = 4000, chains = 2, adapt = 1000, cores = 1)
    }, error = function(e) {
      if (mcp_errors == 0) {
        message("  [mcp first error] ", conditionMessage(e))
      }
      mcp_errors <<- mcp_errors + 1
      NULL
    })
    if (is.null(fit_mcp)) next
    s <- tryCatch(summary(fit_mcp), error = function(e) NULL)
    if (is.null(s)) next
    cp_row <- s[s$name == "cp_1", ]
    if (nrow(cp_row) == 0) next
    mcp_results[[ens]] <- data.frame(
      gene = g, ensembl = ens,
      cp_mean = cp_row$mean, cp_lower = cp_row$lower,
      cp_upper = cp_row$upper,
      rhat = if ("Rhat" %in% colnames(cp_row)) cp_row$Rhat else NA_real_,
      stringsAsFactors = FALSE)
  }
  if (length(mcp_results) > 0) {
    mcp_df <- bind_rows(mcp_results)
    write_csv(mcp_df, file.path(TAB, "bayesian_changepoint_mcp.csv"))
    log_info(paste0("  mcp 成功：", nrow(mcp_df), " / ", length(sig16_ens),
                    " genes, median cp = ",
                    round(median(mcp_df$cp_mean, na.rm = TRUE), 2),
                    ", max Rhat = ",
                    round(max(mcp_df$rhat, na.rm = TRUE), 3)))
  } else {
    log_info(paste0("  mcp 全部失败（", mcp_errors, " errors），改用 bcp"))
  }
} else {
  log_info("  mcp 未安装")
}

bcp_df <- NULL
if (is.null(mcp_df) && has_bcp) {
  log_info("  运行 bcp ...")
  bcp_results <- list()
  for (ens in sig16_ens) {
    g <- ens2gene[ens]
    y <- get_expr(ens)
    ord <- order(base_df$grade)
    y_sorted <- y[ord]
    res <- tryCatch(
      bcp(y_sorted, return.samples = FALSE, burnin = 500, mcmc = 2000),
      error = function(e) NULL)
    if (is.null(res)) next
    pp <- res$posterior.prob
    cp_idx <- which.max(pp)
    bcp_results[[ens]] <- data.frame(
      gene = g, ensembl = ens,
      cp_index = cp_idx,
      cp_grade = base_df$grade[ord][cp_idx],
      cp_prob = pp[cp_idx],
      stringsAsFactors = FALSE)
  }
  if (length(bcp_results) > 0) {
    bcp_df <- bind_rows(bcp_results)
    write_csv(bcp_df, file.path(TAB, "bayesian_changepoint_bcp.csv"))
    log_info(paste0("  bcp 成功：", nrow(bcp_df), " / ", length(sig16_ens),
                    " genes, median cp_grade = ",
                    round(median(bcp_df$cp_grade, na.rm = TRUE), 2),
                    ", median cp_prob = ",
                    round(median(bcp_df$cp_prob, na.rm = TRUE), 3)))
  } else {
    log_info("  bcp 也失败")
  }
}

# ============================================================
# ============================================================
log_info("--- 9. Effect size ---")
effect_results <- list()
for (ens in sig16_ens) {
  g <- ens2gene[ens]; y <- get_expr(ens)
  res <- tryCatch({
    dat <- data.frame(y = y, base_df)
    fit_full <- lm(y ~ poly(grade, 3) + ns(age, df = 3) + sex + rin + pmi, data = dat)
    fit_null <- lm(y ~ ns(age, df = 3) + sex + rin + pmi, data = dat)
    r2_full <- summary(fit_full)$r.squared
    r2_null <- summary(fit_null)$r.squared
    partial_r2 <- (r2_full - r2_null) / (1 - r2_null)
    g2 <- y[base_df$grade == 2]; g3 <- y[base_df$grade == 3]
    cohens_d <- (mean(g3) - mean(g2)) / sqrt((var(g2) + var(g3)) / 2)
    data.frame(gene = g, ensembl = ens, partial_R2 = partial_r2,
               cohens_d_G2vsG3 = cohens_d,
               mean_G2 = mean(g2), mean_G3 = mean(g3),
               stringsAsFactors = FALSE)
  }, error = function(e) {
    data.frame(gene = g, ensembl = ens, partial_R2 = NA_real_,
               cohens_d_G2vsG3 = NA_real_, mean_G2 = NA_real_,
               mean_G3 = NA_real_, stringsAsFactors = FALSE)
  })
  effect_results[[ens]] <- res
}
effect_df <- bind_rows(effect_results) %>% arrange(desc(partial_R2))
write_csv(effect_df, file.path(TAB, "effect_size_results.csv"))
log_info(paste0("Effect size: median partial R2 = ",
                round(median(effect_df$partial_R2, na.rm = TRUE), 4)))

# ============================================================
# 10. Bootstrap
# ============================================================
log_info("--- 10. Bootstrap p3 stability ---")
n_boot <- 1000
boot_results <- list()
for (i_ens in seq_along(sig16_ens)) {
  ens <- sig16_ens[i_ens]
  g <- ens2gene[ens]; y <- get_expr(ens)
  set.seed(42 + i_ens)
  p3_sig_count <- 0; n_valid <- 0
  for (b in seq_len(n_boot)) {
    idx <- sample(length(y), replace = TRUE)
    dat_b <- data.frame(
      y = y[idx], grade = base_df$grade[idx], age = base_df$age[idx],
      sex = base_df$sex[idx], rin = base_df$rin[idx], pmi = base_df$pmi[idx])
    res_b <- fit_p3_safe(dat_b$y, dat_b)
    if (!is.na(res_b$p3_p)) {
      n_valid <- n_valid + 1
      if (res_b$p3_p < 0.05) p3_sig_count <- p3_sig_count + 1
    }
  }
  boot_results[[ens]] <- data.frame(
    gene = g, ensembl = ens,
    pct_sig_p3 = ifelse(n_valid > 0, p3_sig_count / n_valid, NA_real_),
    n_valid = n_valid, stringsAsFactors = FALSE)
}
boot_df <- bind_rows(boot_results) %>% arrange(desc(pct_sig_p3))
write_csv(boot_df, file.path(TAB, "bootstrap_p3_16genes.csv"))
log_info(paste0("Bootstrap: median pct_sig_p3 = ",
                round(median(boot_df$pct_sig_p3, na.rm = TRUE), 3)))

# ============================================================
# 11. Leave-one-out
# ============================================================
log_info("--- 11. Leave-one-out ---")
target_sym <- intersect(c("FCF1", "NDUFB8", "UGGT1", "S100B"), sig16_map$symbol)
target_ens <- gene2ens[target_sym]
log_info(paste0("LOO 目标基因：", paste(target_sym, collapse = ", ")))

loo_sample_results <- list()
for (ens in target_ens) {
  g <- ens2gene[ens]; y <- get_expr(ens); n <- length(y)
  p3_pvals <- numeric(n)
  for (i in seq_len(n)) {
    dat_i <- data.frame(
      y = y[-i], grade = base_df$grade[-i], age = base_df$age[-i],
      sex = base_df$sex[-i], rin = base_df$rin[-i], pmi = base_df$pmi[-i])
    res_i <- fit_p3_safe(dat_i$y, dat_i)
    p3_pvals[i] <- res_i$p3_p
  }
  n_valid_loo <- sum(!is.na(p3_pvals))
  log_info(paste0("  ", g, ": valid LOO p = ", n_valid_loo, " / ", n))
  loo_sample_results[[ens]] <- data.frame(
    gene = g, ensembl = ens,
    n_valid_loo = n_valid_loo,
    p3_p_loo_median = median(p3_pvals, na.rm = TRUE),
    p3_p_loo_max = if (all(is.na(p3_pvals))) NA_real_ else max(p3_pvals, na.rm = TRUE),
    pct_loo_sig = mean(p3_pvals < 0.05, na.rm = TRUE),
    stringsAsFactors = FALSE)
}
loo_sample_df <- bind_rows(loo_sample_results)
write_csv(loo_sample_df, file.path(TAB, "FCF1_leave_one_out.csv"))
log_info(paste0("Sample-level LOO: ",
                paste(loo_sample_df$gene, round(loo_sample_df$pct_loo_sig, 2),
                      sep = "=", collapse = ", ")))

sig16_noFCF1_sym <- setdiff(sig16_map$symbol, "FCF1")
sig16_noFCF1_ens <- gene2ens[sig16_noFCF1_sym]
poly3_noFCF1 <- list()
for (ens in sig16_noFCF1_ens) {
  g <- ens2gene[ens]
  dat <- data.frame(y = get_expr(ens), grade = base_df$grade,
                    age = base_df$age, sex = base_df$sex,
                    rin = base_df$rin, pmi = base_df$pmi)
  res <- fit_p3_safe(dat$y, dat)
  poly3_noFCF1[[ens]] <- data.frame(gene = g, ensembl = ens,
                                    p3_p = res$p3_p, stringsAsFactors = FALSE)
}
poly3_noFCF1_df <- bind_rows(poly3_noFCF1) %>%
  mutate(p3_FDR = p.adjust(p3_p, method = "BH")) %>% arrange(p3_p)
write_csv(poly3_noFCF1_df, file.path(TAB, "poly3F_without_FCF1.csv"))
log_info(paste0("Without FCF1: ",
                sum(poly3_noFCF1_df$p3_FDR < 0.05, na.rm = TRUE),
                " / ", nrow(poly3_noFCF1_df), " genes FDR < 0.05"))

# ============================================================
# ============================================================
log_info("--- 12. Panel direction balance ---")
panel_dir_df <- NULL; panel_sign_df <- NULL
if (file.exists(gene_p3_file)) {
  p3_genes_raw <- read_csv(gene_p3_file, show_col_types = FALSE)

  log_info(paste0("poly3F_gene_level.csv 列名：",
                  paste(colnames(p3_genes_raw), collapse = ", ")))

  sym_col <- intersect(c("SYMBOL","symbol","gene_name"), colnames(p3_genes_raw))[1]
  lfc_col <- intersect(c("p3_logFC","p3_coef","p3_beta"), colnames(p3_genes_raw))[1]

  if (!is.na(sym_col) && !is.na(lfc_col)) {
    p3_genes <- p3_genes_raw %>%
      rename(SYMBOL = !!sym(sym_col), p3_logFC = !!lfc_col) %>%
      filter(!is.na(SYMBOL), !is.na(p3_logFC))

    log_info(paste0("p3_logFC 总体：n = ", nrow(p3_genes),
                    ", mean = ", round(mean(p3_genes$p3_logFC), 4),
                    ", median = ", round(median(p3_genes$p3_logFC), 4),
                    ", pct_pos = ", round(mean(p3_genes$p3_logFC > 0), 3)))
    set.seed(1)
    sample_rows <- p3_genes[sample(nrow(p3_genes), min(5, nrow(p3_genes))), ]
    log_info("  随机 5 行（用于核对方向定义）：")
    print(sample_rows)

    panels <- list(
      OXPHOS = grep("^(NDUF|SDH|UQCR|COX|ATP5)", p3_genes$SYMBOL, value = TRUE),
      Complement = grep("^(C[0-9]|CF[BDHI]|SERPING)", p3_genes$SYMBOL, value = TRUE),
      IEG = c("FOS","FOSB","JUN","JUNB","EGR1","EGR2","EGR3",
              "ATF3","ZFP36","DUSP1","IER2","SOCS3"),
      Photoreceptor = grep("^(RHO|NR2E3|NRL|GNAT1|PDE6)", p3_genes$SYMBOL, value = TRUE),
      Apoptosis = grep("^(BAX|BAK1|CASP|BCL2|BID|PMAIP1)", p3_genes$SYMBOL, value = TRUE),
      UPR = grep("^(ATF4|ATF6|ERN1|EIF2AK3|HSPA5|DDIT3|XBP1)", p3_genes$SYMBOL, value = TRUE),
      Spliceosome = grep("^(SNRP|SRSF|HNRNP|RBM)", p3_genes$SYMBOL, value = TRUE))

    dir_rows <- list(); sign_rows <- list()
    for (pn in names(panels)) {
      sub <- p3_genes %>% filter(SYMBOL %in% panels[[pn]])
      log_info(paste0("  Panel ", pn, ": n = ", nrow(sub),
                      ", mean_logFC = ", round(mean(sub$p3_logFC, na.rm = TRUE), 4)))
      if (nrow(sub) < 5) next
      n_pos <- sum(sub$p3_logFC > 0, na.rm = TRUE)
      n_neg <- sum(sub$p3_logFC < 0, na.rm = TRUE)
      dir_rows[[pn]] <- data.frame(
        panel = pn, n_genes = nrow(sub),
        n_p3_positive = n_pos, n_p3_negative = n_neg,
        pct_positive = n_pos / nrow(sub),
        mean_p3_logFC = mean(sub$p3_logFC, na.rm = TRUE),
        sd_p3_logFC = sd(sub$p3_logFC, na.rm = TRUE),
        cancel_ratio = min(n_pos, n_neg) / max(n_pos, n_neg),
        stringsAsFactors = FALSE)
      binom_res <- binom.test(n_pos, nrow(sub), p = 0.5,
                              alternative = "two.sided")
      sign_rows[[pn]] <- data.frame(panel = pn, n_pos = n_pos,
        n_total = nrow(sub), binom_p = binom_res$p.value,
        stringsAsFactors = FALSE)
    }
    if (length(dir_rows) > 0) {
      panel_dir_df <- bind_rows(dir_rows)
      write_csv(panel_dir_df, file.path(TAB, "panel_p3_direction.csv"))
    }
    if (length(sign_rows) > 0) {
      panel_sign_df <- bind_rows(sign_rows) %>%
        mutate(binom_FDR = p.adjust(binom_p, method = "BH"))
      write_csv(panel_sign_df, file.path(TAB, "panel_p3_sign_test.csv"))
    }
    log_info(paste0("Panel direction: ", nrow(panel_dir_df), " panels"))
  }
}

# ============================================================
# ============================================================
log_info("--- 13. Figure S13 ---")
pA <- ggplot(spline_df, aes(x = reorder(gene, spline_p), y = -log10(spline_p))) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(x = "", y = "-log10(P)", title = "A  Natural spline model") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        plot.title = element_text(face = "bold"))

pB <- ggplot(boot_df, aes(x = reorder(gene, pct_sig_p3), y = pct_sig_p3)) +
  geom_col(fill = "darkgreen") +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
  labs(x = "", y = "Proportion bootstrap F-test P<0.05",
       title = "B  Bootstrap stability (1000x)") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        plot.title = element_text(face = "bold"))

pC <- if (!is.null(panel_dir_df) && nrow(panel_dir_df) > 0) {
  ggplot(panel_dir_df, aes(x = reorder(panel, -pct_positive), y = pct_positive)) +
    geom_col(fill = "orange") +
    geom_hline(yintercept = 0.5, linetype = "dashed") + ylim(0, 1) +
    labs(x = "", y = "Proportion p3-positive",
         title = "C  Panel direction balance") +
    theme_classic(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          plot.title = element_text(face = "bold"))
} else ggplot() + theme_void() + ggtitle("C  Panel (not available)")

pD <- ggplot(effect_df, aes(x = reorder(gene, cohens_d_G2vsG3),
                            y = cohens_d_G2vsG3,
                            fill = cohens_d_G2vsG3 > 0)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey30") +
  scale_fill_manual(values = c("TRUE" = "firebrick", "FALSE" = "steelblue"),
                    guide = "none") +
  labs(x = "", y = "Cohen's d (G3 vs G2)",
       title = "D  Effect size G2 vs G3") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        plot.title = element_text(face = "bold"))

fig_s13 <- (pA | pB) / (pC | pD) +
  plot_annotation(
    title = "G2-G3 signal robustness assessment",
    subtitle = "Natural spline / bootstrap / panel balance / effect size",
    theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
                  plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey30")))
ggsave(file.path(FIG, "FigureS13_G2G3_robustness.pdf"), fig_s13,
       width = 12, height = 9, device = cairo_pdf)
ggsave(file.path(FIG, "FigureS13_G2G3_robustness.png"), fig_s13,
       width = 12, height = 9, dpi = 300)

# ============================================================
# 14. summary
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("=== 02_8 G2-G3 robustness summary (v3) ===\n")
cat(strrep("=", 70), "\n", sep = "")
cat("\n[0] 数据概况:\n")
cat("    Samples: ", nrow(base_df), "\n", sep = "")
cat("    16-gene mapped: ", length(sig16_ens), " / 16\n", sep = "")
cat("    grade 分布: ",
    paste(names(table(base_df$grade)), as.integer(table(base_df$grade)),
          sep = "=", collapse = ", "), "\n", sep = "")

cat("\n[1] Natural spline:\n")
cat("    FDR < 0.05: ", sum(spline_df$spline_FDR < 0.05, na.rm = TRUE),
    " / ", nrow(spline_df), "\n", sep = "")
cat("    Median partial R2: ",
    round(median(spline_df$partial_R2, na.rm = TRUE), 4), "\n", sep = "")

if (!is.null(seg_df)) {
  cat("\n[2] Segmented regression (仅记录，不建议用于论文):\n")
  cat("    Fitted: ", nrow(seg_df), " genes\n", sep = "")
  cat("    Median breakpoint: ",
      round(median(seg_df$breakpoint, na.rm = TRUE), 2), "\n", sep = "")
  cat("    Davies FDR < 0.05: ",
      sum(seg_df$davies_FDR < 0.05, na.rm = TRUE),
      " / ", nrow(seg_df), "  ← 若为 0 则无统计证据\n", sep = "")
}

cat("\n[2b] Bayesian changepoint:\n")
if (!is.null(mcp_df)) {
  cat("    mcp: ", nrow(mcp_df), " genes, median cp = ",
      round(median(mcp_df$cp_mean, na.rm = TRUE), 2),
      ", max Rhat = ", round(max(mcp_df$rhat, na.rm = TRUE), 3), "\n", sep = "")
} else if (!is.null(bcp_df)) {
  cat("    bcp: ", nrow(bcp_df), " genes, median cp_grade = ",
      round(median(bcp_df$cp_grade, na.rm = TRUE), 2),
      ", median cp_prob = ",
      round(median(bcp_df$cp_prob, na.rm = TRUE), 3), "\n", sep = "")
} else {
  cat("    mcp 和 bcp 均不可用，跳过\n")
}

cat("\n[3] Effect size:\n")
cat("    Median partial R2: ",
    round(median(effect_df$partial_R2, na.rm = TRUE), 4), "\n", sep = "")
cat("    Median |Cohen's d| (G2 vs G3): ",
    round(median(abs(effect_df$cohens_d_G2vsG3), na.rm = TRUE), 3), "\n", sep = "")

cat("\n[4] Bootstrap stability:\n")
cat("    Genes with > 80% sig: ",
    sum(boot_df$pct_sig_p3 > 0.8, na.rm = TRUE),
    " / ", nrow(boot_df), "\n", sep = "")
cat("    Median pct_sig: ",
    round(median(boot_df$pct_sig_p3, na.rm = TRUE), 3), "\n", sep = "")

cat("\n[5] Sample-level LOO:\n")
for (i in seq_len(nrow(loo_sample_df))) {
  cat("    ", loo_sample_df$gene[i],
      ": valid = ", loo_sample_df$n_valid_loo[i],
      ", pct_loo_sig = ", round(loo_sample_df$pct_loo_sig[i], 3),
      ", median P = ", round(loo_sample_df$p3_p_loo_median[i], 4), "\n", sep = "")
}

cat("\n[6] Without FCF1:\n")
cat("    FDR < 0.05: ",
    sum(poly3_noFCF1_df$p3_FDR < 0.05, na.rm = TRUE),
    " / ", nrow(poly3_noFCF1_df), "\n", sep = "")

if (!is.null(panel_dir_df)) {
  cat("\n[7] Panel direction balance:\n")
  for (i in seq_len(nrow(panel_dir_df))) {
    cat("    ", panel_dir_df$panel[i],
        ": n = ", panel_dir_df$n_genes[i],
        ", mean_logFC = ", round(panel_dir_df$mean_p3_logFC[i], 4),
        ", pct_pos = ", round(panel_dir_df$pct_positive[i], 3),
        "\n", sep = "")
  }
  if (!is.null(panel_sign_df)) {
    cat("\n[8] Panel sign test (two-sided):\n")
    print(as.data.frame(panel_sign_df), row.names = FALSE)
  }
  cat("\n[9] ⚠️ 方向核对提醒：\n")
  cat("    请确认 poly3F_gene_level.csv 的 p3_logFC 是 G3-G2 的 logFC，\n")
  cat("    还是 poly(grade,3) 的最高次项系数。后者正负无生物学意义。\n")
}

cat("\n[Files saved]\n"); print(list.files(TAB)); print(list.files(FIG))
cat("\nOUT:", OUT, "\n")
log_info("=== 02_8 done ===")


