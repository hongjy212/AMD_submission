# ============================================================
# 07_5d_metabolism_pathways.r (v2: 419)
# Purpose: compare metabolic pathway scores vs AMD grade
#   - Key pathways: OXPHOS / Glycolysis / Hypoxia / mTORC1 /
#     FAO / Mito_biogenesis / Mito_mitophagy / Mito_antioxidant /
#     Mito_dynamics / Mito_mtDNA_rep
#   - Validate OXPHOS downregulation at pathway level
#   - Test Warburg hypothesis (classic: Glycolysis up + OXPHOS down)
#   - Compare with greenyellow module (OXPHOS module)
#
#
# Input:
#   results/05_1_consensus_wgcna/data/consensusNet.rds
#   results/01_1b_batch_enhance/data/{expr_wgcna,clinical_wgcna}.rds   ★ 419
#   results/01_1_prepare_bulk_data/data/gene_annotation.rds
#   results/01_2_prepare_gene_sets/data/gene_sets.rds
# Output:
#   results/07_5d_metabolism_pathways/
# ============================================================
source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
})

ROOT <- "E:/retina_tools/aging_rgc_axon/AMD"
OUT  <- file.path(ROOT, "results", "07_5d_metabolism_pathways")
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

log_info("=== 07_5d metabolism pathways v2 (419) ===")

# ============================================================
# 0. Load  ★ v2：419
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

ens2sym <- setNames(toupper(anno$symbol), sub("\\..*$", "", anno$ensembl))
expr_sym <- expr
rownames(expr_sym) <- ens2sym[rownames(expr_sym)]
common <- intersect(colnames(expr_sym), clin$sample_id)
expr_sym <- expr_sym[, common, drop = FALSE]
clin2 <- clin[match(common, clin$sample_id), ]
log_info("Samples: {length(common)}, genes: {nrow(expr_sym)}")

# ============================================================
# 0b. Module availability check
# ============================================================
mod_colors <- net$colors
available_mods <- sort(unique(mod_colors))
log_info("WGCNA module colors: {paste(available_mods, collapse=', ')}")

requested_mods <- c("greenyellow", "darkorange")
available_requested <- intersect(requested_mods, available_mods)
missing_requested <- setdiff(requested_mods, available_mods)
if (length(missing_requested) > 0) {
  log_warn("Modules requested but NOT in WGCNA: {paste(missing_requested, collapse=', ')}")
  log_warn("  -> will be skipped; check consensusNet.rds if this is unexpected")
}

# ============================================================
# 1. Key pathways
# ============================================================
key_sets <- c(
  "OXPHOS",
  "Glycolysis",
  "Hypoxia",
  "mTORC1",
  "Fatty_acid_oxidation",
  "Mito_biogenesis",
  "Mito_mitophagy",
  "Mito_antioxidant",
  "Mito_dynamics",
  "Mito_mtDNA_rep"
)

pathway_res <- data.frame()
for (sn in key_sets) {
  if (!sn %in% names(GENE_SETS)) {
    log_warn("Missing gene set: {sn}")
    next
  }
  sg <- toupper(GENE_SETS[[sn]]$genes)
  sg_in <- intersect(sg, rownames(expr_sym))
  if (length(sg_in) < 5) {
    log_warn("{sn}: only {length(sg_in)} genes in expr, skip")
    next
  }

  sub <- expr_sym[sg_in, , drop = FALSE]
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

  pathway_res <- rbind(pathway_res, data.frame(
    pathway = sn,
    n_genes = length(sg_in),
    beta = co["grade","Estimate"],
    se = co["grade","Std. Error"],
    p = co["grade","Pr(>|t|)"],
    rho = rho$estimate,
    rho_p = rho$p.value,
    dir = ifelse(co["grade","Estimate"] > 0, "up", "down")))
}
pathway_res$FDR <- p.adjust(pathway_res$p, "BH")
pathway_res <- pathway_res[order(pathway_res$p), ]
rownames(pathway_res) <- NULL

write_tab(pathway_res, "metabolic_pathway_vs_grade.csv")

cat("\n=== Metabolic pathways vs AMD grade (419) ===\n")
print(pathway_res, row.names = FALSE, digits = 3)

# ============================================================
# 2. Warburg hypothesis (classic definition)
# ============================================================
log_info("--- Warburg test ---")

get_path <- function(nm) pathway_res[pathway_res$pathway == nm, ]
gly  <- get_path("Glycolysis")
hyp  <- get_path("Hypoxia")
oxp  <- get_path("OXPHOS")
mtor <- get_path("mTORC1")

is_sig <- function(df, direction = c("up", "down")) {
  direction <- match.arg(direction)
  if (nrow(df) == 0) return(FALSE)
  if (direction == "up")   return(df$beta > 0 && df$p < 0.05)
  if (direction == "down") return(df$beta < 0 && df$p < 0.05)
}

gly_up   <- is_sig(gly,  "up")
hyp_up   <- is_sig(hyp,  "up")
oxp_down <- is_sig(oxp,  "down")
mtor_up  <- is_sig(mtor, "up")

cat("\n=== Warburg hypothesis (classic definition) ===\n")
cat("  Glycolysis UP (p<0.05):", gly_up, "\n")
cat("  Hypoxia UP    (p<0.05):", hyp_up, "\n")
cat("  OXPHOS DOWN   (p<0.05):", oxp_down, "\n")
cat("  mTORC1 UP     (p<0.05):", mtor_up, "\n")

if (gly_up && oxp_down) {
  cat("\n-> WARBURG-LIKE: glycolysis up + OXPHOS down\n")
  if (hyp_up)  cat("   + hypoxia up (hypoxic Warburg)\n")
  if (mtor_up) cat("   + mTORC1 up (anabolic reprogramming)\n")
} else if (oxp_down && !gly_up && !hyp_up) {
  cat("\n-> OXPHOS-SELECTIVE DOWNREGULATION (not Warburg)\n")
} else if (hyp_up && oxp_down && !gly_up) {
  cat("\n-> HYPOXIA-DRIVEN OXPHOS SUPPRESSION (not classic Warburg)\n")
} else if (oxp_down) {
  cat("\n-> PARTIAL OXPHOS DOWNREGULATION (glycolysis not up)\n")
} else {
  cat("\n-> NO CLEAR METABOLIC SIGNATURE\n")
}

evidence_score <- data.frame(
  gly_up = gly_up, hyp_up = hyp_up,
  oxp_down = oxp_down, mtor_up = mtor_up,
  n_evidence = sum(gly_up, hyp_up, oxp_down),
  classic_warburg = gly_up && oxp_down)
write_tab(evidence_score, "warburg_evidence.csv")

# ============================================================
# 3. Module score analysis (greenyellow + optional darkorange)
# ============================================================
log_info("--- Module score analysis ---")

module_analysis <- function(mod_name, net_colors, ens2sym, expr_sym, clin2) {
  ens <- names(net_colors)[net_colors == mod_name]
  if (length(ens) == 0) {
    log_warn("{mod_name}: 0 genes in WGCNA, skip")
    return(NULL)
  }
  sym <- toupper(na.omit(ens2sym[sub("\\..*$", "", ens)]))
  sym_in <- intersect(sym, rownames(expr_sym))
  if (length(sym_in) < 5) {
    log_warn("{mod_name}: only {length(sym_in)} genes in expr, skip")
    return(NULL)
  }
  log_info("{mod_name}: {length(ens)} WGCNA genes, {length(sym_in)} matched")

  sub <- expr_sym[sym_in, , drop = FALSE]
  z <- t(scale(t(sub)))
  score <- colMeans(z, na.rm = TRUE)

  df <- data.frame(score = as.numeric(score), grade = clin2$grade,
                   age = clin2$age_n, sex = clin2$sex,
                   rin = clin2$rin_n, pmi = clin2$pmi_n)
  df <- df[complete.cases(df), ]

  fit <- lm(score ~ grade + age + sex + rin + pmi, data = df)
  co <- summary(fit)$coefficients
  rho <- cor.test(df$score, df$grade, method = "spearman", exact = FALSE)

  list(
    module = mod_name,
    n_genes = length(sym_in),
    beta = co["grade","Estimate"],
    p = co["grade","Pr(>|t|)"],
    rho = rho$estimate,
    rho_p = rho$p.value,
    score = score,
    df = df
  )
}

gy_res <- module_analysis("greenyellow", mod_colors, ens2sym, expr_sym, clin2)
do_res <- module_analysis("darkorange",  mod_colors, ens2sym, expr_sym, clin2)

if (!is.null(gy_res)) {
  cat("\n=== greenyellow module score ===\n")
  cat(sprintf("  beta = %.4f, p = %.4f\n", gy_res$beta, gy_res$p))
  cat(sprintf("  rho = %.3f, p = %.4f\n",  gy_res$rho, gy_res$rho_p))
  cat("  Mean score by grade:\n")
  print(tapply(gy_res$df$score, gy_res$df$grade, mean))
}
if (!is.null(do_res)) {
  cat("\n=== darkorange module score ===\n")
  cat(sprintf("  beta = %.4f, p = %.4f\n", do_res$beta, do_res$p))
  cat(sprintf("  rho = %.3f, p = %.4f\n",  do_res$rho, do_res$rho_p))
  cat("  Mean score by grade:\n")
  print(tapply(do_res$df$score, do_res$df$grade, mean))
}

module_summary <- dplyr::bind_rows(
  if (!is.null(gy_res)) data.frame(
    module = "greenyellow",
    n_genes = gy_res$n_genes,
    beta = gy_res$beta,
    p = gy_res$p,
    rho = gy_res$rho,
    rho_p = gy_res$rho_p) else NULL,
  if (!is.null(do_res)) data.frame(
    module = "darkorange",
    n_genes = do_res$n_genes,
    beta = do_res$beta,
    p = do_res$p,
    rho = do_res$rho,
    rho_p = do_res$rho_p) else NULL
)
if (nrow(module_summary) > 0) {
  write_tab(module_summary, "module_gene_vs_grade.csv")
} else {
  log_warn("No module results to save (all requested modules missing)")
}

# ============================================================
# 4. Figures
# ============================================================
log_info("--- Figures ---")

if (nrow(pathway_res) > 0) {
  pathway_res$pathway <- factor(pathway_res$pathway,
                                 levels = pathway_res$pathway[order(pathway_res$beta)])
  p1 <- ggplot(pathway_res, aes(x = pathway, y = beta, fill = dir)) +
    geom_col() +
    geom_errorbar(aes(ymin = beta - se, ymax = beta + se), width = 0.2) +
    geom_hline(yintercept = 0, linetype = 2) +
    scale_fill_manual(values = c("up" = "#C44E52", "down" = "#4C72B0"),
                       name = "direction") +
    coord_flip() +
    theme_bw(12) +
    labs(x = "", y = "beta (per grade)",
         title = "Metabolic pathway score vs AMD grade (419)") +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(file.path(FIG, "Fig1_metabolic_pathways.pdf"), p1,
         width = 6, height = 5)
}

if (nrow(oxp) > 0) {
  oxp_genes <- intersect(toupper(GENE_SETS$OXPHOS$genes), rownames(expr_sym))
  oxp_expr <- expr_sym[oxp_genes, , drop = FALSE]
  oxp_score <- colMeans(t(scale(t(oxp_expr))), na.rm = TRUE)
  df_oxp <- data.frame(score = oxp_score, grade = clin2$grade)
  p2 <- ggplot(df_oxp, aes(x = factor(grade), y = score)) +
    geom_boxplot(outlier.shape = NA, width = 0.6,
                  fill = "#4C72B0", alpha = 0.5) +
    geom_jitter(width = 0.15, alpha = 0.4, size = 1) +
    geom_smooth(aes(group = 1), method = "lm", se = TRUE,
                 color = "#C44E52", alpha = 0.2) +
    theme_bw(12) +
    labs(x = "AMD grade", y = "OXPHOS score",
         title = "OXPHOS pathway score by AMD grade (419)",
         subtitle = sprintf("beta=%.4f, p=%.4f, FDR=%.4f",
                            oxp$beta, oxp$p, oxp$FDR)) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(file.path(FIG, "Fig2_OXPHOS_score.pdf"), p2,
         width = 6, height = 4)
}

mod_list <- list()
if (!is.null(gy_res)) mod_list[["greenyellow"]] <- gy_res$score
if (!is.null(do_res)) mod_list[["darkorange"]]  <- do_res$score

if (length(mod_list) > 0) {
  df_dual <- data.frame(
    sample_id = common,
    grade = clin2$grade,
    as.data.frame(mod_list, check.names = FALSE))
  df_long <- tidyr::pivot_longer(df_dual, cols = names(mod_list),
                                  names_to = "score_type",
                                  values_to = "score")
  p3 <- ggplot(df_long, aes(x = factor(grade), y = score, fill = score_type)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.8) +
    facet_wrap(~ score_type, scales = "free_y") +
    scale_fill_manual(values = c("greenyellow" = "#ADFF2F",
                                  "darkorange" = "#FF8C00")) +
    theme_bw(12) +
    labs(x = "AMD grade", y = "score",
         title = "Module scores vs AMD grade (419)") +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5))
  ggsave(file.path(FIG, "Fig3_module_scores.pdf"), p3,
         width = 4 * length(mod_list), height = 4)
}

# ============================================================
# Summary
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("=== 07_5d summary (v2, 419) ===\n")
cat(strrep("=", 70), "\n", sep = "")

cat("\n[Pathways vs AMD grade]:\n")
print(pathway_res[, c("pathway","n_genes","beta","p","FDR","dir")],
      row.names = FALSE, digits = 3)

if (!is.null(gy_res)) {
  cat("\n[greenyellow module]:\n")
  cat(sprintf("  beta = %.4f, p = %.4f, rho = %.3f\n",
              gy_res$beta, gy_res$p, gy_res$rho))
}
if (!is.null(do_res)) {
  cat("\n[darkorange module]:\n")
  cat(sprintf("  beta = %.4f, p = %.4f, rho = %.3f\n",
              do_res$beta, do_res$p, do_res$rho))
}

cat("\n[Warburg hypothesis]:\n")
cat("  Glycolysis UP:", gly_up, "\n")
cat("  Hypoxia UP:   ", hyp_up, "\n")
cat("  OXPHOS DOWN:  ", oxp_down, "\n")
if (gly_up && oxp_down) {
  cat("  -> WARBURG-LIKE\n")
} else if (oxp_down && !gly_up && !hyp_up) {
  cat("  -> OXPHOS-selective downregulation (NOT Warburg)\n")
} else {
  cat("  -> No clear Warburg signature\n")
}

cat("\nOUT =", OUT, "\n")
log_info("=== 07_5d v2 (419) done ===")

