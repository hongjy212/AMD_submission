# ============================================================
# 07_5c_OXPHOS_analysis.r (v8: 418, outlier removed + per-SD beta
#                            + OXPHOS pathway per-SD)
# Purpose: comprehensive analysis of the OXPHOS genes from 07_5b
#   Part A: Fisher enrichment (greenyellow / p3 / modules vs eQTM)
#   Part B: single-gene regression + FDR
#   Part C: complex-level summary
#   Part D: composite score + Stouffer + binom + sample-level CSV
#   Part D2: per-1-SD beta (adjusted and unadjusted) for composite
#   Part E: figures (all English labels)
#
#   - Output OXPHOS_pathway_perSD.csv
#
#   - Output OXPHOS13_composite_perSD.csv
#
#     samples ID: R42015-506pf_294-IR_L3
#
#
# Input:
#   results/07_5b_GSE231536_attachments/tables/greenyellow_eqtm_genes.csv
#   results/01_1b_batch_enhance/data/{expr_wgcna,clinical_wgcna}.rds
#   results/01_1_prepare_bulk_data/data/gene_annotation.rds
#   results/05_1_consensus_wgcna/data/consensusNet.rds
#   results/01_2_prepare_gene_sets/data/gene_sets.rds
#   results/04_2_pathway_score/tables/pathway_poly3F_adjusted.csv
#   GSE231536 MOESM12/14 (eQTM targets)
# Output:
#   results/07_5c_OXPHOS_analysis/
# ============================================================
source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(readxl)
})

ROOT <- "E:/retina_tools/aging_rgc_axon/AMD"
OUT  <- file.path(ROOT, "results", "07_5c_OXPHOS_analysis")
TAB  <- file.path(OUT, "tables")
FIG  <- file.path(OUT, "figures")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

write_tab <- function(df, name) {
  p <- file.path(TAB, name)
  write.csv(df, p, row.names = FALSE)
  log_info("saved: {name} ({nrow(df)} rows)")
  invisible(p)
}

log_info("=== 07_5c OXPHOS analysis v8 (418 + per-SD + pathway per-SD) ===")

# ============================================================
# 0. Load data
# ============================================================
expr <- readRDS(file.path(ROOT, "results", "01_1b_batch_enhance",
                           "data", "expr_wgcna.rds"))
clin <- readRDS(file.path(ROOT, "results", "01_1b_batch_enhance",
                           "data", "clinical_wgcna.rds"))
anno <- readRDS(file.path(ROOT, "results", "01_1_prepare_bulk_data",
                           "data", "gene_annotation.rds"))
net  <- readRDS(file.path(ROOT, "results", "05_1_consensus_wgcna",
                           "data", "consensusNet.rds"))
GENE_SETS <- readRDS(file.path(ROOT, "results", "01_2_prepare_gene_sets",
                                "data", "gene_sets.rds"))

log_info("expr: {nrow(expr)} × {ncol(expr)}")
log_info("clin: {nrow(clin)} samples")

ens2sym <- setNames(toupper(anno$symbol), sub("\\..*$", "", anno$ensembl))
expr_sym <- expr
rownames(expr_sym) <- ens2sym[rownames(expr_sym)]
common <- intersect(colnames(expr_sym), clin$sample_id)
expr_sym <- expr_sym[, common, drop = FALSE]
clin2 <- clin[match(common, clin$sample_id), ]
log_info("common samples: {length(common)}")

gy_file <- file.path(ROOT, "results", "07_5b_GSE231536_attachments",
                      "tables", "greenyellow_eqtm_genes.csv")
if (file.exists(gy_file)) {
  gy_genes <- toupper(read.csv(gy_file)$gene)
  log_info("Loaded from 07_5b: {length(gy_genes)} genes")
} else {
  stop("07_5b not run; please run 07_5b first")
}

gy_expr <- expr_sym[rownames(expr_sym) %in% gy_genes, , drop = FALSE]
log_info("Matched genes: {nrow(gy_expr)} / {length(gy_genes)}")

# ============================================================
# ============================================================
gy_mean <- colMeans(gy_expr, na.rm = TRUE)
gy_z <- (gy_mean - mean(gy_mean, na.rm = TRUE)) / sd(gy_mean, na.rm = TRUE)
outlier_samples <- names(gy_z)[gy_z < -5]

if (length(outlier_samples) > 0) {
  log_warn(sprintf("Excluding %d outlier sample(s) with 13-gene mean z < -5: %s",
                   length(outlier_samples),
                   paste(outlier_samples, collapse = ", ")))
  gy_expr <- gy_expr[, !colnames(gy_expr) %in% outlier_samples, drop = FALSE]

  clin2 <- clin2[clin2$sample_id %in% colnames(gy_expr), ]
  clin2 <- clin2[match(colnames(gy_expr), clin2$sample_id), ]

  log_info(sprintf("After outlier exclusion: %d samples, %d genes",
                   ncol(gy_expr), nrow(gy_expr)))
} else {
  log_info("No outlier samples detected")
}

cat("\n=== Sample distribution after exclusion ===\n")
print(table(clin2$grade))

# ============================================================
# Part A: Fisher enrichment
# ============================================================
log_info("--- Part A: Fisher ---")

DIR_MET <- "E:/retina_tools/raw/GSE231536"

extract_targets <- function(df) {
  col <- grep("Significant.*target genes|eQTM target genes",
              colnames(df), ignore.case = TRUE)[1]
  if (is.na(col)) return(character(0))
  s <- paste(df[[col]], collapse = ",")
  g <- trimws(unlist(strsplit(s, "[,;]")))
  g <- g[!is.na(g) & g != "" & g != "-"]
  toupper(unique(g))
}

m12 <- read_excel(file.path(DIR_MET, "41467_2024_46063_MOESM12_ESM.xlsx"),
                  sheet = "eQTM_target_genes_GO")
m14 <- read_excel(file.path(DIR_MET, "41467_2024_46063_MOESM14_ESM.xlsx"))
eqtm_targets <- unique(c(extract_targets(m12), extract_targets(m14)))
log_info("eQTM targets: {length(eqtm_targets)}")

sig16 <- toupper(read.csv(file.path(
  ROOT, "results", "02_1_differential_expression_main",
  "tables", "sig16_genes.csv"))$gene)

polyF <- read.csv(file.path(
  ROOT, "results", "04_2_pathway_score",
  "tables", "pathway_poly3F_adjusted.csv"))
p3_paths <- polyF$pathway[polyF$pattern == "p3_only"]
p3_genes <- unique(toupper(unlist(lapply(
  intersect(p3_paths, names(GENE_SETS)),
  function(p) GENE_SETS[[p]]$genes))))

mod_sym <- list(
  greenyellow = toupper(na.omit(ens2sym[sub("\\..*$", "",
                   names(net$colors)[net$colors == "greenyellow"])])),
  black       = toupper(na.omit(ens2sym[sub("\\..*$", "",
                   names(net$colors)[net$colors == "black"])])),
  green       = toupper(na.omit(ens2sym[sub("\\..*$", "",
                   names(net$colors)[net$colors == "green"])]))
)

targets <- list(
  sig16           = sig16,
  p3_pathway      = p3_genes,
  mod_greenyellow = mod_sym$greenyellow,
  mod_black       = mod_sym$black,
  mod_green       = mod_sym$green
)

N_BG <- 20000
fisher_enrich <- function(genes, ref_set, N = N_BG) {
  a <- length(intersect(genes, ref_set))
  b <- length(genes) - a
  c <- length(ref_set) - a
  d <- N - a - b - c
  if (d < 0) d <- 0
  m <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  ft <- tryCatch(fisher.test(m, alternative = "greater",
                              simulate.p.value = TRUE, B = 1e5),
                  error = function(e) NULL)
  if (is.null(ft)) return(list(OR = NA, p = NA))
  list(OR = as.numeric(ft$estimate), p = ft$p.value)
}

fisher_res <- data.frame()
for (tgt in names(targets)) {
  g <- targets[[tgt]]
  n_g <- length(g)
  n_hit <- length(intersect(g, eqtm_targets))
  fr <- fisher_enrich(g, eqtm_targets)
  fisher_res <- rbind(fisher_res, data.frame(
    target    = tgt,
    n_total   = n_g,
    n_hit     = n_hit,
    pct_hit   = round(100 * n_hit / n_g, 2),
    n_ref     = length(eqtm_targets),
    fold      = round((n_hit / n_g) / (length(eqtm_targets) / N_BG), 2),
    OR        = round(fr$OR, 2),
    p_fisher  = fr$p))
}
fisher_res$FDR <- p.adjust(fisher_res$p_fisher, method = "BH")
fisher_res <- fisher_res[order(fisher_res$p_fisher), ]
rownames(fisher_res) <- NULL
write_tab(fisher_res, "fisher_eqtm_enrichment.csv")

# ============================================================
# Part B: single-gene regression
# ============================================================
log_info("--- Part B: single-gene regression ---")

res <- data.frame()
for (g in rownames(gy_expr)) {
  y <- gy_expr[g, ]
  df <- data.frame(y = y, grade = clin2$grade, age = clin2$age_n,
                   sex = clin2$sex, rin = clin2$rin_n, pmi = clin2$pmi_n)
  df <- df[complete.cases(df), ]
  fit <- lm(y ~ grade + age + sex + rin + pmi, data = df)
  co <- summary(fit)$coefficients
  res <- rbind(res, data.frame(
    gene = g,
    beta = co["grade","Estimate"],
    se   = co["grade","Std. Error"],
    t    = co["grade","t value"],
    p    = co["grade","Pr(>|t|)"]))
}
res$FDR <- p.adjust(res$p, method = "BH")
res <- res[order(res$p), ]
rownames(res) <- NULL
write_tab(res, "OXPHOS13_gene_vs_grade.csv")

cat("\n=== Part B: genes vs grade ===\n")
print(res, row.names = FALSE, digits = 3)
cat("\nFDR<0.05:", sum(res$FDR < 0.05), "/", nrow(res), "\n")
cat("beta<0:  ", sum(res$beta < 0),    "/", nrow(res), "\n")

# ============================================================
# Part C: complex-level summary
# ============================================================
log_info("--- Part C: complex ---")

complex_map <- c(
  NDUFAB1 = "I", NDUFB1 = "I", NDUFB2 = "I", NDUFB3 = "I",
  NDUFC1  = "I", NDUFA12 = "I",
  SDHB    = "II",
  UQCRQ   = "III", UQCRFS1 = "III", UQCRH = "III",
  COX7A2  = "IV", COX7C = "IV",
  GABARAPL2 = "Autophagy"
)
res$complex <- complex_map[res$gene]

comp_summary <- data.frame()
for (cp in unique(res$complex)) {
  sub <- res[res$complex == cp, ]
  comp_summary <- rbind(comp_summary, data.frame(
    complex = cp, n = nrow(sub),
    mean_beta = mean(sub$beta), min_p = min(sub$p)))
}
write_tab(comp_summary, "OXPHOS13_complex_summary.csv")

complexes <- setdiff(unique(res$complex), "Autophagy")
comp_res <- data.frame()
for (cp in complexes) {
  genes <- names(complex_map)[complex_map == cp]
  genes <- intersect(genes, rownames(gy_expr))
  if (length(genes) < 2) next
  sub <- gy_expr[genes, , drop = FALSE]
  z <- t(scale(t(sub)))
  sc <- colMeans(z, na.rm = TRUE)
  df <- data.frame(y = as.numeric(sc), grade = clin2$grade,
                   age = clin2$age_n, sex = clin2$sex,
                   rin = clin2$rin_n, pmi = clin2$pmi_n)
  df <- df[complete.cases(df), ]
  fit <- lm(y ~ grade + age + sex + rin + pmi, data = df)
  co <- summary(fit)$coefficients
  rho <- cor.test(as.numeric(sc), clin2$grade,
                   method = "spearman", exact = FALSE)
  comp_res <- rbind(comp_res, data.frame(
    complex = cp, n_subunits = length(genes),
    beta = co["grade","Estimate"], p = co["grade","Pr(>|t|)"],
    rho = rho$estimate, rho_p = rho$p.value))
}
comp_res$FDR <- p.adjust(comp_res$p, "BH")
comp_res <- comp_res[order(comp_res$p), ]
rownames(comp_res) <- NULL
write_tab(comp_res, "OXPHOS13_complex_score.csv")

cat("\n=== Part C: complex score ===\n")
print(comp_res, row.names = FALSE, digits = 3)

# ============================================================
# Part D: composite score + Stouffer + binom + sample-level CSV
# ============================================================
log_info("--- Part D: composite score ---")

z_all <- t(scale(t(gy_expr)))
oxphos_score <- colMeans(z_all, na.rm = TRUE)

df_score <- data.frame(
  score = oxphos_score, grade = clin2$grade,
  age = clin2$age_n, sex = clin2$sex,
  rin = clin2$rin_n, pmi = clin2$pmi_n)
df_score <- df_score[complete.cases(df_score), ]

fit_score <- lm(score ~ grade + age + sex + rin + pmi, data = df_score)
co_s <- summary(fit_score)$coefficients
rho_s <- cor.test(df_score$score, df_score$grade,
                   method = "spearman", exact = FALSE)

z_st <- abs(qnorm(res$p/2)) * sign(res$beta)
Z <- sum(z_st) / sqrt(length(z_st))
p_st <- pnorm(Z)

log_info("Stouffer diagnostic: n_pos_z={sum(z_st > 0)}, n_neg_z={sum(z_st < 0)}, Z={round(Z,3)}, p={format(p_st, scientific=TRUE)}")

binom_p <- binom.test(sum(res$beta < 0), nrow(res), 0.5)$p.value

cat("\n=== Part D: composite score ===\n")
cat(sprintf("  beta = %.4f, p = %.4f\n",
            co_s["grade","Estimate"], co_s["grade","Pr(>|t|)"]))
cat(sprintf("  rho = %.3f, p = %.4f\n",
            rho_s$estimate, rho_s$p.value))
cat(sprintf("  Stouffer Z = %.3f, p = %.4g\n", Z, p_st))
cat(sprintf("  binom p = %.4g\n", binom_p))

summary_row <- data.frame(
  n_genes = nrow(res),
  n_beta_neg = sum(res$beta < 0),
  n_p05 = sum(res$p < 0.05),
  n_fdr05 = sum(res$FDR < 0.05),
  score_beta = co_s["grade","Estimate"],
  score_p = co_s["grade","Pr(>|t|)"],
  score_rho = rho_s$estimate,
  score_rho_p = rho_s$p.value,
  stouffer_Z = Z,
  stouffer_p = p_st,
  binom_p = binom_p)
write_tab(summary_row, "OXPHOS13_summary.csv")

# ============================================================
# ============================================================
log_info("--- Part D2: per-SD beta (composite score) ---")

df_score_z <- df_score %>%
  mutate(score_z = as.numeric(scale(score)),
         grade_z = as.numeric(scale(grade)))

fit_z_adj <- lm(score_z ~ grade_z + age + sex + rin + pmi, data = df_score_z)
co_z_adj <- summary(fit_z_adj)$coefficients

beta_perSD_adj <- co_z_adj["grade_z", "Estimate"]
se_perSD_adj   <- co_z_adj["grade_z", "Std. Error"]
p_perSD_adj    <- co_z_adj["grade_z", "Pr(>|t|)"]
ci_lo_adj      <- beta_perSD_adj - 1.96 * se_perSD_adj
ci_hi_adj      <- beta_perSD_adj + 1.96 * se_perSD_adj

fit_z_unadj <- lm(score_z ~ grade_z, data = df_score_z)
co_z_unadj <- summary(fit_z_unadj)$coefficients

beta_perSD_unadj <- co_z_unadj["grade_z", "Estimate"]
se_perSD_unadj   <- co_z_unadj["grade_z", "Std. Error"]
p_perSD_unadj    <- co_z_unadj["grade_z", "Pr(>|t|)"]
ci_lo_unadj      <- beta_perSD_unadj - 1.96 * se_perSD_unadj
ci_hi_unadj      <- beta_perSD_unadj + 1.96 * se_perSD_unadj

cat("\n=== Part D2: per-1-SD beta (composite) ===\n")
cat(sprintf("  Adjusted   (age, sex, RIN, PMI): beta = %.4f (95%% CI %.4f to %.4f), p = %.4f\n",
            beta_perSD_adj, ci_lo_adj, ci_hi_adj, p_perSD_adj))
cat(sprintf("  Unadjusted                     : beta = %.4f (95%% CI %.4f to %.4f), p = %.4f\n",
            beta_perSD_unadj, ci_lo_unadj, ci_hi_unadj, p_perSD_unadj))

perSD_row <- data.frame(
  model       = c("adjusted", "unadjusted"),
  beta_per_SD = c(beta_perSD_adj, beta_perSD_unadj),
  se          = c(se_perSD_adj, se_perSD_unadj),
  ci_lower    = c(ci_lo_adj, ci_lo_unadj),
  ci_upper    = c(ci_hi_adj, ci_hi_unadj),
  p_value     = c(p_perSD_adj, p_perSD_unadj),
  n           = c(nrow(df_score_z), nrow(df_score_z)),
  covariates  = c("age + sex + rin + pmi", "none"),
  stringsAsFactors = FALSE)

write_tab(perSD_row, "OXPHOS13_composite_perSD.csv")

# ============================================================
# ============================================================
log_info("--- Part D3: OXPHOS pathway score per-SD beta ---")

pathway_file <- file.path(ROOT, "results", "04_2_pathway_score",
                           "data", "pathway_scores_ssgsea.rds")
if (file.exists(pathway_file)) {
  pathway_scores <- readRDS(pathway_file)
  cat("pathway_scores 类型:", class(pathway_scores), "\n")
  cat("pathway_scores 维度:", paste(dim(pathway_scores), collapse = " x "), "\n")

  if (is.matrix(pathway_scores)) {
    ps <- as.data.frame(t(pathway_scores))
    ps$sample_id <- rownames(ps)
  } else {
    ps <- as.data.frame(pathway_scores)
    if (!"sample_id" %in% colnames(ps)) {
      ps$sample_id <- rownames(ps)
    }
  }

  if ("OXPHOS" %in% colnames(ps)) {
    ps_use <- ps %>%
      select(sample_id, OXPHOS) %>%
      left_join(clin2 %>% select(sample_id, grade,
                                  age = age_n, sex,
                                  rin = rin_n, pmi = pmi_n),
                by = "sample_id") %>%
      filter(!is.na(OXPHOS), !is.na(grade))

    cat("OXPHOS pathway score 样本数:", nrow(ps_use), "\n")

    ps_use <- ps_use %>%
      mutate(score_z = as.numeric(scale(OXPHOS)),
             grade_z = as.numeric(scale(grade)))

    fit_ox <- lm(score_z ~ grade_z + age + sex + rin + pmi, data = ps_use)
    co_ox <- summary(fit_ox)$coefficients

    beta_ox <- co_ox["grade_z", "Estimate"]
    se_ox   <- co_ox["grade_z", "Std. Error"]
    p_ox    <- co_ox["grade_z", "Pr(>|t|)"]
    ci_lo_ox <- beta_ox - 1.96 * se_ox
    ci_hi_ox <- beta_ox + 1.96 * se_ox

    fit_ox_unadj <- lm(score_z ~ grade_z, data = ps_use)
    co_ox_unadj <- summary(fit_ox_unadj)$coefficients

    beta_ox_unadj <- co_ox_unadj["grade_z", "Estimate"]
    se_ox_unadj   <- co_ox_unadj["grade_z", "Std. Error"]
    p_ox_unadj    <- co_ox_unadj["grade_z", "Pr(>|t|)"]
    ci_lo_ox_unadj <- beta_ox_unadj - 1.96 * se_ox_unadj
    ci_hi_ox_unadj <- beta_ox_unadj + 1.96 * se_ox_unadj

    cat(sprintf("\nOXPHOS pathway score per 1 SD (adjusted):\n"))
    cat(sprintf("  beta = %.4f (95%% CI %.4f to %.4f), p = %.4f\n",
                beta_ox, ci_lo_ox, ci_hi_ox, p_ox))
    cat(sprintf("OXPHOS pathway score per 1 SD (unadjusted):\n"))
    cat(sprintf("  beta = %.4f (95%% CI %.4f to %.4f), p = %.4f\n",
                beta_ox_unadj, ci_lo_ox_unadj, ci_hi_ox_unadj, p_ox_unadj))

    oxphos_perSD <- data.frame(
      model       = c("adjusted", "unadjusted"),
      beta_per_SD = c(beta_ox, beta_ox_unadj),
      se          = c(se_ox, se_ox_unadj),
      ci_lower    = c(ci_lo_ox, ci_lo_ox_unadj),
      ci_upper    = c(ci_hi_ox, ci_hi_ox_unadj),
      p_value     = c(p_ox, p_ox_unadj),
      n           = c(nrow(ps_use), nrow(ps_use)),
      covariates  = c("age + sex + rin + pmi", "none"),
      stringsAsFactors = FALSE)

    write_tab(oxphos_perSD, "OXPHOS_pathway_perSD.csv")
  } else {
    log_warn("Part D3: OXPHOS column not found in pathway_scores_ssgsea.rds")
    log_warn("Available columns: {paste(colnames(ps), collapse=', ')}")
  }
} else {
  log_warn("Part D3: pathway_scores_ssgsea.rds not found; skipped")
}

# ============================================================
# Part D4: sample-level composite score CSV
# ============================================================
composite_sample <- data.frame(
  sample_id = colnames(gy_expr),
  score     = as.numeric(oxphos_score),
  n_genes   = colSums(!is.na(z_all)),
  grade     = as.integer(clin2$grade),
  stringsAsFactors = FALSE
) %>%
  arrange(sample_id)

write_tab(composite_sample, "OXPHOS13_composite_score.csv")
log_info("composite sample-level: {nrow(composite_sample)} samples, n_genes range = {min(composite_sample$n_genes)}-{max(composite_sample$n_genes)}")
cat("\nsamples:", nrow(composite_sample), "\n")
print(table(composite_sample$grade))

# ============================================================
# Part E: figures
# ============================================================
log_info("--- Part E: figures ---")

res_plot <- res
res_plot$gene <- factor(res_plot$gene,
                         levels = res_plot$gene[order(res_plot$beta)])
p1 <- ggplot(res_plot, aes(x = gene, y = beta,
                             fill = FDR < 0.05)) +
  geom_col() +
  geom_errorbar(aes(ymin = beta - se, ymax = beta + se), width = 0.2) +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_fill_manual(values = c("TRUE" = "#C44E52", "FALSE" = "#B0B0B0"),
                     name = "FDR<0.05") +
  coord_flip() + theme_bw(12) +
  labs(x = "", y = "beta (per grade)",
       title = "OXPHOS candidate genes vs AMD grade")
ggsave(file.path(FIG, "Fig1_gene_beta.pdf"), p1, width = 6, height = 5)

p2 <- ggplot(df_score, aes(x = factor(grade), y = score)) +
  geom_boxplot(outlier.shape = NA, width = 0.6,
                fill = "#4C72B0", alpha = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE,
               color = "#C44E52", alpha = 0.2) +
  theme_bw(12) +
  labs(x = "AMD grade", y = "OXPHOS score",
       title = sprintf("OXPHOS score vs grade (rho=%.3f, p=%.3f)",
                       rho_s$estimate, rho_s$p.value))
ggsave(file.path(FIG, "Fig2_oxphos_score.pdf"), p2, width = 6, height = 4)

if (nrow(comp_res) > 0) {
  comp_res$complex <- factor(comp_res$complex,
                              levels = comp_res$complex[order(comp_res$beta)])
  p3 <- ggplot(comp_res, aes(x = complex, y = beta,
                               fill = FDR < 0.05)) +
    geom_col() +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_fill_manual(values = c("TRUE" = "#C44E52", "FALSE" = "#B0B0B0"),
                       name = "FDR<0.05") +
    theme_bw(12) +
    labs(x = "Complex", y = "beta (per grade)",
         title = "OXPHOS complex score vs AMD grade")
  ggsave(file.path(FIG, "Fig3_complex_beta.pdf"), p3, width = 5, height = 4)
}

if (nrow(fisher_res) > 0) {
  fisher_res$target <- factor(fisher_res$target,
                               levels = fisher_res$target[order(fisher_res$OR)])
  p4 <- ggplot(fisher_res, aes(x = target, y = OR,
                                 fill = FDR < 0.05)) +
    geom_col() +
    geom_hline(yintercept = 1, linetype = 2) +
    scale_fill_manual(values = c("TRUE" = "#C44E52", "FALSE" = "#B0B0B0"),
                       name = "FDR<0.05") +
    coord_flip() + theme_bw(12) +
    labs(x = "", y = "Odds Ratio",
         title = "eQTM target enrichment by target")
  ggsave(file.path(FIG, "Fig4_fisher_OR.pdf"), p4, width = 6, height = 4)
}

# ============================================================
# Summary
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("=== 07_5c summary (v8, 418 + per-SD + pathway per-SD) ===\n")
cat(strrep("=", 70), "\n", sep = "")

cat("\n[Part A] Fisher enrichment:\n")
print(fisher_res[, c("target","n_total","n_hit","fold","OR",
                       "p_fisher","FDR")], row.names = FALSE)

cat("\n[Part B] genes:\n")
cat("  beta<0:  ", sum(res$beta < 0), "/", nrow(res), "\n")
cat("  p<0.05:  ", sum(res$p < 0.05), "/", nrow(res), "\n")
cat("  FDR<0.05:", sum(res$FDR < 0.05), "/", nrow(res), "\n")

cat("\n[Part C] Complexes:\n")
print(comp_res[, c("complex","n_subunits","beta","p","FDR")],
      row.names = FALSE, digits = 3)

cat("\n[Part D] Composite score:\n")
cat(sprintf("  beta = %.4f, p = %.4f\n",
            co_s["grade","Estimate"], co_s["grade","Pr(>|t|)"]))
cat(sprintf("  rho = %.3f, p = %.4f\n",
            rho_s$estimate, rho_s$p.value))
cat(sprintf("  Stouffer Z = %.3f, p = %.4g\n", Z, p_st))
cat(sprintf("  binom p = %.4g\n", binom_p))

cat("\n[Part D2] Per-1-SD beta (composite):\n")
print(perSD_row, row.names = FALSE, digits = 4)

if (exists("oxphos_perSD")) {
  cat("\n[Part D3] Per-1-SD beta (OXPHOS pathway):\n")
  print(oxphos_perSD, row.names = FALSE, digits = 4)
}

log_info("=== 07_5c v8 done ===")

