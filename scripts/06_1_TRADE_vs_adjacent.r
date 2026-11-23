#!/usr/bin/env Rscript
# ============================================================
# 06_1_TRADE_vs_adjacent.r (v2)
#
# v2 Fix：
#
#   04_x → pathway/score
#   05_x → WGCNA module
#
#   tables/TRADE_impact_3modes.csv
#   tables/TRADE_impact_decomposition.csv
#   tables/DEG_pure_<comp>.csv
#   tables/panel_perturbation_3modes.csv
#   tables/poly3F_panel_scores.csv
#   tables/poly3F_gene_level.csv
#   tables/poly3F_panel_summary.csv
#   tables/poly3F_vs_PPF.csv
#   figures/TRADE_impact_comparison.pdf/.png
#   figures/panel_PPF_3v2.pdf/.png
#   figures/poly3F_panel_curves.pdf/.png
#   figures/poly3F_vs_PPF_scatter.pdf
# ============================================================

if (exists(".AMD_SETUP_LOADED", envir = .GlobalEnv)) {
  rm(".AMD_SETUP_LOADED", envir = .GlobalEnv)
}
source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")

suppressPackageStartupMessages({
  library(limma); library(TRADEtools); library(GSVA)
  library(dplyr); library(tidyr); library(ggplot2)
  library(pheatmap); library(splines)
})

DIR_THIS <- setup_output_dir("06_1_TRADE_vs_adjacent")
log_info("=== 06_1_TRADE_vs_adjacent.r (v2) start ===")

# ============================================================
# ============================================================
log_info("--- 1. 载入 ---")

expr <- read_upstream("01_1b_batch_enhance", "data", "expr_wgcna.rds")
clin <- read_upstream("01_1b_batch_enhance", "data", "clinical_wgcna.rds")
anno <- read_upstream("01_1_prepare_bulk_data", "data", "gene_annotation.rds")
GENE_SETS <- read_upstream("01_2_prepare_gene_sets", "data", "gene_sets.rds")

common <- intersect(colnames(expr), rownames(clin))
expr <- expr[, common, drop = FALSE]
clin <- clin[common, , drop = FALSE]

log_info("expr: {nrow(expr)} × {ncol(expr)}")
log_info("gene sets: {length(GENE_SETS)} 个")

# ============================================================
# ============================================================
log_info("--- 2. 整理变量 ---")

clin$grade <- as.integer(clin$grade)
clin$age_v <- as.numeric(clin$age_n)
clin$rin_v <- as.numeric(clin$rin_n)
clin$pmi_v <- as.numeric(clin$pmi_n)
clin$batch_f <- factor(clin$batch)
sid <- as.character(clin$sample_id)
clin$layer_f <- factor(sub("^.*_(L[0-9]+)$", "\\1", sid),
                        levels = paste0("L", 1:8))
clin$sex_f <- factor(clin$sex)
rownames(clin) <- colnames(expr)

log_info("grade: {paste(names(table(clin$grade)), table(clin$grade), sep='=', collapse=' ')}")
log_info("batch: {paste(names(table(clin$batch_f)), table(clin$batch_f), sep='=', collapse=' ')}")
log_info("layer NA: {sum(is.na(clin$layer_f))}")

# ============================================================
# ============================================================
comparisons <- list(
  list(label = "2_vs_1", high = 2, low = 1),
  list(label = "3_vs_2", high = 3, low = 2),
  list(label = "4_vs_3", high = 4, low = 3)
)

design_functions <- list(
  total = function(cs) {
    grade_f <- factor(cs$grade)
    model.matrix(~ 0 + grade_f, data = cs)
  },
  adjusted = function(cs) {
    grade_f <- factor(cs$grade)
    model.matrix(~ 0 + grade_f + age_v + sex_f, data = cs)
  },
  pure = function(cs) {
    cs <- as.data.frame(cs)
    cs$grade_f <- factor(cs$grade)
    cs$batch_f <- droplevels(factor(cs$batch_f))
    cs$layer_f <- droplevels(factor(cs$layer_f))

    pure_vars <- c("grade_f", "age_v", "sex_f")
    if (nlevels(cs$batch_f) > 1) pure_vars <- c(pure_vars, "batch_f")
    if (nlevels(cs$layer_f) > 1) pure_vars <- c(pure_vars, "layer_f")
    pure_vars <- c(pure_vars, "rin_v", "pmi_v")

    for (vn in pure_vars) {
      if (anyNA(cs[[vn]])) stop("pure design 变量 ", vn, " 含 NA")
    }

    model.matrix(
      as.formula(paste("~ 0 +", paste(pure_vars, collapse = " + "))),
      data = cs)
  }
)

# ============================================================
# ============================================================
log_info("--- 4. TRADE 分析 ---")

trade_summary_all <- data.frame()
deg_all <- list()
lfcSE_diag <- data.frame()

for (design_name in names(design_functions)) {
  design_fn <- design_functions[[design_name]]
  log_info("Design: {design_name}")

  for (comp in comparisons) {
    label <- comp$label
    high  <- comp$high
    low   <- comp$low

    keep_idx <- which(clin$grade %in% c(high, low))
    clinical_sub <- clin[keep_idx, , drop = FALSE]
    expr_sub <- expr[, keep_idx, drop = FALSE]

    design_sub <- design_fn(clinical_sub)
    colnames(design_sub)[1:2] <- c(paste0("G", low), paste0("G", high))

    stopifnot(nrow(design_sub) == ncol(expr_sub))
    rank_sub <- qr(design_sub)$rank
    if (rank_sub < ncol(design_sub)) {
      log_warn("  [{design_name}/{label}] design 秩亏，跳过")
      next
    }

    contrasts_mat <- limma::makeContrasts(
      contrast = get(paste0("G", high)) - get(paste0("G", low)),
      levels = design_sub)

    fit_init     <- lmFit(expr_sub, design_sub)
    fit_contrast <- contrasts.fit(fit_init, contrasts_mat)
    fit_eb       <- eBayes(fit_contrast, trend = TRUE)

    se_check  <- fit_init$stdev.unscaled * fit_init$sigma
    se_median <- median(sqrt(se_check[, paste0("G", low)]^2 +
                              se_check[, paste0("G", high)]^2))

    logFC_vec <- fit_eb$coefficients[, "contrast"]
    se_vec    <- fit_eb$stdev.unscaled[, "contrast"] * fit_eb$sigma
    pval_vec  <- fit_eb$p.value[, "contrast"]

    res <- data.frame(
      Gene      = rownames(fit_eb$coefficients),
      logFC     = logFC_vec,
      lfcSE     = se_vec,
      t         = logFC_vec / se_vec,
      P.Value   = pval_vec,
      adj.P.Val = p.adjust(pval_vec, method = "BH"),
      stringsAsFactors = FALSE)
    res$symbol <- anno$symbol[match(res$Gene, anno$ensembl)]
    res <- res[order(res$P.Value), ]

    if (design_name == "pure") {
      write.csv(res,
                file.path(DIR_THIS, "tables", paste0("DEG_pure_", label, ".csv")),
                row.names = FALSE)
      deg_all[[label]] <- res
    }

    n_sig_05 <- sum(res$adj.P.Val < 0.05, na.rm = TRUE)
    n_sig_10 <- sum(res$adj.P.Val < 0.10, na.rm = TRUE)

    lfcSE_diag <- rbind(lfcSE_diag, data.frame(
      Design = design_name, Comparison = label,
      n_genes = length(se_vec),
      se_median = round(median(se_vec, na.rm = TRUE), 6),
      se_max    = round(max(se_vec, na.rm = TRUE), 6),
      lfc_median = round(median(abs(logFC_vec), na.rm = TRUE), 6),
      design_rank = rank_sub, design_ncol = ncol(design_sub),
      stringsAsFactors = FALSE))

    trade_input <- data.frame(
      log2FoldChange = logFC_vec,
      lfcSE          = se_vec,
      pvalue         = pval_vec,
      row.names      = rownames(fit_eb$coefficients))
    trade_input <- trade_input[complete.cases(trade_input), ]
    trade_input <- trade_input[trade_input$lfcSE > 0, ]

    trade_res <- tryCatch({
      TRADE(
        mode = "univariate", results1 = trade_input,
        annot_table = NULL,
        log2FoldChange = "log2FoldChange",
        lfcSE = "lfcSE", pvalue = "pvalue",
        model_significant = TRUE, n_sample = 2000, verbose = FALSE)
    }, error = function(e) {
      log_warn("  [{design_name}/{label}] TRADE 失败: {e$message}")
      NULL
    })

    if (!is.null(trade_res)) {
      impact <- trade_res$distribution_summary$transcriptome_wide_impact
      Me     <- trade_res$distribution_summary$Me

      trade_summary_all <- rbind(trade_summary_all, data.frame(
        Design = design_name, Comparison = label, Stage = high,
        TranscriptomeImpact = impact, EffectiveDEGs = Me,
        n_DEG_FDR05 = n_sig_05, n_DEG_FDR10 = n_sig_10,
        n_samples = nrow(clinical_sub), se_median = se_median,
        stringsAsFactors = FALSE))
    }
  }
}

write.csv(trade_summary_all,
          file.path(DIR_THIS, "tables", "TRADE_impact_3modes.csv"),
          row.names = FALSE)
write.csv(lfcSE_diag,
          file.path(DIR_THIS, "tables", "lfcSE_diagnostic.csv"),
          row.names = FALSE)

log_info("TRADE 结果:")
print(trade_summary_all, row.names = FALSE, digits = 4)

# ============================================================
# ============================================================
impact_wide <- trade_summary_all %>%
  dplyr::select(Design, Comparison, TranscriptomeImpact) %>%
  pivot_wider(names_from = Design, values_from = TranscriptomeImpact)
impact_wide <- impact_wide %>%
  mutate(delta_age_sex = total - adjusted,
         delta_full    = total - pure,
         pure_pct      = round(100 * pure / total, 2))
write.csv(impact_wide,
          file.path(DIR_THIS, "tables", "TRADE_impact_decomposition.csv"),
          row.names = FALSE)

log_info("Impact 分解:")
print(impact_wide, row.names = FALSE, digits = 4)

# ============================================================
# ============================================================
log_info("--- 6. Panel 扰动分解 ---")

key_panels <- list(
  IEG               = GENE_SETS$IEG$genes,
  OXPHOS            = GENE_SETS$OXPHOS$genes,
  Phototransduction = GENE_SETS$Phototransduction$genes,
  Photoreceptor_OS  = GENE_SETS$Photoreceptor_OuterSegment$genes,
  Photoreceptor_IS  = GENE_SETS$Photoreceptor_InnerSegment$genes,
  Photoreceptor_CC  = GENE_SETS$Photoreceptor_ConnectingCilium$genes,
  Notch             = GENE_SETS$Notch$genes,
  UPR               = GENE_SETS$ER_stress$genes,
  Complement        = GENE_SETS$Complement$genes,
  Spliceosome       = GENE_SETS$Spliceosome$genes,
  Mito_biogenesis   = GENE_SETS$Mito_biogenesis$genes
)

log_info("key_panels: {length(key_panels)} 个 panel")

panel_gene_counts <- data.frame()
sym_of_expr <- toupper(anno$symbol[match(rownames(expr), anno$ensembl)])
sym_of_expr <- sym_of_expr[!is.na(sym_of_expr)]

for (pn in names(key_panels)) {
  genes <- toupper(key_panels[[pn]])
  n_total   <- length(genes)
  n_in_data <- sum(genes %in% sym_of_expr)
  panel_gene_counts <- rbind(panel_gene_counts, data.frame(
    Panel             = pn,
    n_genes_in_set    = n_total,
    n_genes_in_data   = n_in_data,
    pct_in_data       = round(100 * n_in_data / n_total, 1),
    stringsAsFactors  = FALSE))
}
write.csv(panel_gene_counts,
          file.path(DIR_THIS, "tables", "panel_gene_counts.csv"),
          row.names = FALSE)

log_info("Panel 基因数:")
print(panel_gene_counts, row.names = FALSE)

panel_perturbation <- data.frame()
for (label in names(deg_all)) {
  deg_df <- deg_all[[label]]
  deg_df$symbol_upper <- toupper(deg_df$symbol)
  if (nrow(deg_df) == 0) next
  global_mean_abs_lfc <- mean(abs(deg_df$logFC), na.rm = TRUE)

  for (panel_name in names(key_panels)) {
    genes <- toupper(key_panels[[panel_name]])
    sub <- deg_df[deg_df$symbol_upper %in% genes, ]
    if (nrow(sub) < 3) next

    panel_mean_abs_lfc <- mean(abs(sub$logFC), na.rm = TRUE)
    ppf <- panel_mean_abs_lfc / global_mean_abs_lfc
    panel_mean_signed <- mean(sub$logFC, na.rm = TRUE)
    direction <- ifelse(panel_mean_signed > 0, "Up", "Down")

    panel_perturbation <- rbind(panel_perturbation, data.frame(
      Comparison = label, Panel = panel_name,
      n_genes_in_panel = nrow(sub),
      mean_abs_logFC = round(panel_mean_abs_lfc, 5),
      mean_signed_logFC = round(panel_mean_signed, 5),
      PPF = round(ppf, 3), direction = direction,
      stringsAsFactors = FALSE))
  }
}
panel_perturbation <- panel_perturbation[
  order(panel_perturbation$Comparison, -panel_perturbation$PPF), ]

write.csv(panel_perturbation,
          file.path(DIR_THIS, "tables", "panel_perturbation_3modes.csv"),
          row.names = FALSE)

ppf_3v2 <- panel_perturbation[panel_perturbation$Comparison == "3_vs_2", ]
log_info("Panel 扰动（3_vs_2）:")
print(ppf_3v2, row.names = FALSE)

# ============================================================
# ============================================================
log_info("--- 7. Panel 层面 poly3_F ---")

panel_ensembl <- list()
for (pn in names(key_panels)) {
  pn_genes <- toupper(key_panels[[pn]])
  hit <- rownames(expr)[sym_of_expr %in% pn_genes]
  if (length(hit) >= 5) {
    panel_ensembl[[pn]] <- hit
  }
}
log_info("可用 panel: {length(panel_ensembl)} / {length(key_panels)}")

set.seed(SEED)
ssgsea_param <- gsvaParam(
  exprData = as.matrix(expr),
  geneSets = panel_ensembl,
  kcdf = "Gaussian",
  minSize = 5,
  maxSize = 5000
)
panel_scores <- gsva(ssgsea_param)
log_info("panel scores: {nrow(panel_scores)} × {ncol(panel_scores)}")

poly_mat <- poly(clin$grade, 3)
colnames(poly_mat) <- c("grade_p1", "grade_p2", "grade_p3")
rownames(poly_mat) <- rownames(clin)

sp_age <- splines::ns(clin$age_v, df = 3)
poly_df <- data.frame(
  poly_mat,
  age_s1  = sp_age[, 1],
  age_s2  = sp_age[, 2],
  age_s3  = sp_age[, 3],
  sex_f   = clin$sex_f,
  batch_f = clin$batch_f,
  layer_f = clin$layer_f,
  rin_v   = clin$rin_v,
  pmi_v   = clin$pmi_v,
  stringsAsFactors = FALSE)
rownames(poly_df) <- rownames(clin)

design_poly <- model.matrix(
  ~ grade_p1 + grade_p2 + grade_p3 +
    age_s1 + age_s2 + age_s3 + sex_f + batch_f + layer_f + rin_v + pmi_v,
  data = poly_df)

log_info("poly design: {nrow(design_poly)} × {ncol(design_poly)}")
log_info("design rank: {qr(design_poly)$rank}")

fit_poly <- eBayes(lmFit(panel_scores, design_poly))

tt_F  <- topTable(fit_poly, coef = c("grade_p1","grade_p2","grade_p3"),
                   number = Inf, sort.by = "none")
tt_p1 <- topTable(fit_poly, coef = "grade_p1", number = Inf, sort.by = "none")
tt_p2 <- topTable(fit_poly, coef = "grade_p2", number = Inf, sort.by = "none")
tt_p3 <- topTable(fit_poly, coef = "grade_p3", number = Inf, sort.by = "none")

poly3F_panel <- data.frame(
  Panel       = rownames(tt_F),
  F_stat      = tt_F$F,
  F_p         = tt_F$P.Value,
  F_FDR       = tt_F$adj.P.Val,
  p1_logFC    = tt_p1$logFC,
  p1_p        = tt_p1$P.Value,
  p1_FDR      = tt_p1$adj.P.Val,
  p2_logFC    = tt_p2$logFC,
  p2_p        = tt_p2$P.Value,
  p2_FDR      = tt_p2$adj.P.Val,
  p3_logFC    = tt_p3$logFC,
  p3_p        = tt_p3$P.Value,
  p3_FDR      = tt_p3$adj.P.Val,
  stringsAsFactors = FALSE)
poly3F_panel <- poly3F_panel[order(poly3F_panel$F_p), ]

write.csv(poly3F_panel,
          file.path(DIR_THIS, "tables", "poly3F_panel_scores.csv"),
          row.names = FALSE)

log_info("Panel 层面 poly3_F:")
print(poly3F_panel[, c("Panel","F_stat","p1_FDR","p2_FDR","p3_FDR",
                        "p3_logFC")],
      row.names = FALSE, digits = 4)

# ============================================================
# ============================================================
log_info("--- 8. 基因层面 poly3_F ---")

fit_gene <- eBayes(lmFit(expr, design_poly))

coef_poly <- c("grade_p1","grade_p2","grade_p3")
tt_gene_F  <- topTable(fit_gene, coef = coef_poly,
                        number = Inf, sort.by = "none")
tt_gene_p1 <- topTable(fit_gene, coef = "grade_p1",
                        number = Inf, sort.by = "none")
tt_gene_p2 <- topTable(fit_gene, coef = "grade_p2",
                        number = Inf, sort.by = "none")
tt_gene_p3 <- topTable(fit_gene, coef = "grade_p3",
                        number = Inf, sort.by = "none")

gene_poly3F <- data.frame(
  ENSEMBL     = rownames(tt_gene_F),
  SYMBOL      = anno$symbol[match(rownames(tt_gene_F), anno$ensembl)],
  F_stat      = tt_gene_F$F,
  F_p         = tt_gene_F$P.Value,
  F_FDR       = tt_gene_F$adj.P.Val,
  p1_logFC    = tt_gene_p1$logFC,
  p1_FDR      = tt_gene_p1$adj.P.Val,
  p2_logFC    = tt_gene_p2$logFC,
  p2_FDR      = tt_gene_p2$adj.P.Val,
  p3_logFC    = tt_gene_p3$logFC,
  p3_FDR      = tt_gene_p3$adj.P.Val,
  stringsAsFactors = FALSE)

write.csv(gene_poly3F,
          file.path(DIR_THIS, "tables", "poly3F_gene_level.csv"),
          row.names = FALSE)

log_info("基因层面 poly3_F:")
log_info("  F_FDR<0.05:  {sum(gene_poly3F$F_FDR < 0.05, na.rm=TRUE)} / {nrow(gene_poly3F)}")
log_info("  p1_FDR<0.05: {sum(gene_poly3F$p1_FDR < 0.05, na.rm=TRUE)}")
log_info("  p2_FDR<0.05: {sum(gene_poly3F$p2_FDR < 0.05, na.rm=TRUE)}")
log_info("  p3_FDR<0.05: {sum(gene_poly3F$p3_FDR < 0.05, na.rm=TRUE)}")

gene_poly3F$SYMBOL_upper <- toupper(gene_poly3F$SYMBOL)
bg_p3 <- mean(gene_poly3F$p3_FDR < 0.05, na.rm = TRUE)
bg_F  <- mean(gene_poly3F$F_FDR < 0.05, na.rm = TRUE)

panel_poly3F_summary <- data.frame()
for (pn in names(key_panels)) {
  pn_genes <- toupper(key_panels[[pn]])
  sub <- gene_poly3F[gene_poly3F$SYMBOL_upper %in% pn_genes, ]
  n_total <- nrow(sub)
  if (n_total < 5) next

  n_p1 <- sum(sub$p1_FDR < 0.05, na.rm = TRUE)
  n_p2 <- sum(sub$p2_FDR < 0.05, na.rm = TRUE)
  n_p3 <- sum(sub$p3_FDR < 0.05, na.rm = TRUE)
  n_F  <- sum(sub$F_FDR < 0.05, na.rm = TRUE)

  panel_poly3F_summary <- rbind(panel_poly3F_summary, data.frame(
    Panel = pn,
    N_genes = n_total,
    N_p1_sig = n_p1, N_p2_sig = n_p2, N_p3_sig = n_p3, N_F_sig = n_F,
    Pct_p1 = round(100 * n_p1 / n_total, 2),
    Pct_p2 = round(100 * n_p2 / n_total, 2),
    Pct_p3 = round(100 * n_p3 / n_total, 2),
    Pct_F  = round(100 * n_F / n_total, 2),
    Enrich_p3 = round((n_p3 / n_total) / bg_p3, 2),
    Enrich_F  = round((n_F / n_total) / bg_F, 2),
    stringsAsFactors = FALSE))
}
panel_poly3F_summary <- panel_poly3F_summary[
  order(-panel_poly3F_summary$Enrich_p3), ]

write.csv(panel_poly3F_summary,
          file.path(DIR_THIS, "tables", "poly3F_panel_summary.csv"),
          row.names = FALSE)

log_info("Panel 聚合（基因层面）:")
print(panel_poly3F_summary, row.names = FALSE)

# ============================================================
# ============================================================
log_info("--- 9. poly3F vs PPF 对比 ---")

cmp <- merge(
  poly3F_panel[, c("Panel","p3_FDR","p3_logFC","F_FDR")],
  panel_poly3F_summary[, c("Panel","Pct_p3","Enrich_p3")],
  by = "Panel", all = TRUE)
cmp <- merge(cmp, ppf_3v2[, c("Panel","PPF","direction")],
              by = "Panel", all = TRUE)
cmp <- cmp[order(cmp$p3_FDR), ]

cmp$pattern <- "other"
cmp$pattern[cmp$p3_FDR < 0.05 & cmp$PPF > 1.2] <- "both"
cmp$pattern[cmp$p3_FDR < 0.05 & cmp$PPF <= 1.2] <- "poly3F_only"
cmp$pattern[cmp$p3_FDR >= 0.05 & cmp$PPF > 1.2] <- "PPF_only"

write.csv(cmp,
          file.path(DIR_THIS, "tables", "poly3F_vs_PPF.csv"),
          row.names = FALSE)

log_info("poly3F vs PPF:")
print(cmp, row.names = FALSE, digits = 4)

cat("\n=== poly3F vs PPF 模式 ===\n")
print(table(cmp$pattern))

# ============================================================
# 10. figure
# ============================================================
log_info("--- 10. 图 ---")

p_impact <- ggplot(trade_summary_all,
                   aes(x = Stage, y = TranscriptomeImpact,
                       color = Design, group = Design)) +
  geom_line(linewidth = 1.2) + geom_point(size = 4) +
  geom_text(aes(label = sprintf("%.5f", TranscriptomeImpact)),
            vjust = -1, size = 3, show.legend = FALSE) +
  scale_color_manual(values = c("total" = "#E41A1C",
                                "adjusted" = "#FF7F00",
                                "pure" = "#377EB8")) +
  scale_x_continuous(breaks = c(2, 3, 4),
                     labels = c("G2 vs G1", "G3 vs G2", "G4 vs G3")) +
  labs(title = "Transcriptome perturbation across adjacent stages",
       subtitle = "total=grade only, adjusted=+age+sex, pure=+batch+layer+RIN+PMI",
       x = "Adjacent-stage comparison",
       y = "TRADE Impact Score", color = "Design") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 9),
        legend.position = "bottom")
save_fig(p_impact, DIR_THIS, "TRADE_impact_comparison", width = 9, height = 6)

# 10b. PPF barplot (3_vs_2)
if (nrow(ppf_3v2) > 0) {
  ppf_3v2$Panel <- factor(ppf_3v2$Panel,
                           levels = rev(ppf_3v2$Panel[order(ppf_3v2$PPF)]))
  p_ppf <- ggplot(ppf_3v2, aes(x = Panel, y = PPF, fill = direction)) +
    geom_col(width = 0.7) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
    geom_text(aes(label = sprintf("%.2f", PPF)), hjust = -0.1, size = 3) +
    scale_fill_manual(values = c("Up" = "#c0392b", "Down" = "#2980b9")) +
    coord_flip() +
    labs(title = "Panel Perturbation Factor (G3_vs_G2)",
         x = NULL, y = "PPF") +
    theme_minimal(base_size = 11)
  save_fig(p_ppf, DIR_THIS, "panel_PPF_3v2", width = 8, height = 6)
}

p_fdr <- poly3F_panel[, c("Panel","p1_FDR","p2_FDR","p3_FDR")]
rownames(p_fdr) <- p_fdr$Panel
p_fdr$Panel <- NULL
p_fdr$p1_FDR[is.na(p_fdr$p1_FDR)] <- 1
p_fdr$p2_FDR[is.na(p_fdr$p2_FDR)] <- 1
p_fdr$p3_FDR[is.na(p_fdr$p3_FDR)] <- 1

pdf(file.path(DIR_THIS, "figures", "poly3F_panel_FDR_heatmap.pdf"),
    width = 7, height = 6)
pheatmap(as.matrix(p_fdr),
         cluster_rows = TRUE, cluster_cols = FALSE,
         display_numbers = TRUE, number_format = "%.3f",
         color = colorRampPalette(c("#c0392b","#FFFFFF"))(50),
         main = "Panel poly3_F: p1/p2/p3 FDR",
         fontsize = 10)
dev.off()

top6 <- head(poly3F_panel$Panel, 6)
curve_df <- data.frame()
for (pn in top6) {
  for (g in 1:4) {
    idx <- clin$grade == g
    curve_df <- rbind(curve_df, data.frame(
      Panel = pn, Grade = g,
      mean_score = mean(panel_scores[pn, idx], na.rm = TRUE),
      se = sd(panel_scores[pn, idx], na.rm = TRUE) / sqrt(sum(idx)),
      stringsAsFactors = FALSE))
  }
}

p_curve <- ggplot(curve_df, aes(x = Grade, y = mean_score,
                                  color = Panel, group = Panel)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean_score - se, ymax = mean_score + se),
                 width = 0.15, alpha = 0.5) +
  theme_bw(base_size = 11) +
  scale_x_continuous(breaks = 1:4, labels = paste0("G", 1:4)) +
  labs(x = "AMD grade", y = "ssGSEA score (mean ± SE)",
       title = "Panel poly3_F: grade curves")
save_fig(p_curve, DIR_THIS, "poly3F_panel_curves", width = 9, height = 5)

p_sc <- ggplot(cmp, aes(x = PPF, y = -log10(p3_FDR))) +
  geom_point(aes(color = pattern), size = 3, alpha = 0.8) +
  geom_text(aes(label = Panel), size = 3, vjust = -0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = 1.2, linetype = "dashed") +
  scale_color_manual(values = c("both" = "#c0392b",
                                 "poly3F_only" = "#3498db",
                                 "PPF_only" = "#f39c12",
                                 "other" = "grey70")) +
  labs(x = "PPF (G3_vs_G2)", y = "-log10(p3_FDR)",
       title = "Panel: poly3_F vs PPF", color = "Pattern") +
  theme_bw(base_size = 11)
save_fig(p_sc, DIR_THIS, "poly3F_vs_PPF_scatter", width = 8, height = 6)

# ============================================================
# 11. summary
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("=== 06_1 汇总 (v2) ===\n")
cat(strrep("=", 70), "\n", sep = "")

cat("\n【Q1】TRADE impact 分解\n")
print(impact_wide, row.names = FALSE, digits = 4)

cat("\n【Q2】Panel 扰动（3_vs_2）\n")
print(ppf_3v2, row.names = FALSE, digits = 4)

cat("\n【Q3】Panel 层面 poly3_F\n")
print(poly3F_panel[, c("Panel","p1_FDR","p2_FDR","p3_FDR","p3_logFC")],
      row.names = FALSE, digits = 4)

cat("\n【Q4】基因层面 poly3_F\n")
cat(sprintf("  F_FDR<0.05:  %d / %d\n",
            sum(gene_poly3F$F_FDR < 0.05, na.rm=TRUE), nrow(gene_poly3F)))
cat(sprintf("  p1_FDR<0.05: %d\n", sum(gene_poly3F$p1_FDR < 0.05, na.rm=TRUE)))
cat(sprintf("  p2_FDR<0.05: %d\n", sum(gene_poly3F$p2_FDR < 0.05, na.rm=TRUE)))
cat(sprintf("  p3_FDR<0.05: %d\n", sum(gene_poly3F$p3_FDR < 0.05, na.rm=TRUE)))

cat("\n【Q5】Panel 基因层面 p3 富集 top 5\n")
print(head(panel_poly3F_summary[, c("Panel","N_genes","Pct_p3","Enrich_p3")], 5),
      row.names = FALSE)

cat("\n【Q6】poly3F vs PPF 模式\n")
print(table(cmp$pattern))

cat("\n--- 输出位置 ---\n")
cat("目录:", DIR_THIS, "\n")

log_info("=== 06_1_TRADE_vs_adjacent.r (v2) done ===")

