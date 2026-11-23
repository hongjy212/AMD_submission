# 10_1C_supp_figures.R (v10.7)
# Supplementary figures S1-S23 (S6 retired; S7-S24 renumbered to S6-S23)

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr); library(tidyr); library(readr)
  library(patchwork); library(pheatmap); library(tibble)
  library(grid); library(gridExtra)
})

results_dir <- DIR_RESULT
fig_dir     <- file.path(results_dir, "10_1C_supp_figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
log_info("=== 10_1C Supplementary figures (v10.7) ===")

theme_pub <- theme_classic(base_size = 11) +
  theme(axis.text = element_text(color = "black"),
        axis.title = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 11))

save_fig <- function(p, name, w = 6, h = 4) {
  suppressMessages({
    ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = w, height = h)
    ggsave(file.path(fig_dir, paste0(name, ".tiff")), p, width = w, height = h, dpi = 300)
  })
  log_info(sprintf("saved figure: %s (%dx%d)", name, w, h))
}

read_required <- function(path) {
  if (!file.exists(path)) stop(sprintf("REQUIRED FILE MISSING: %s", path))
  read_csv(path, show_col_types = FALSE, progress = FALSE)
}

fdr_legend_guides <- guides(
  fill = guide_legend(
    override.aes = list(fill = c("firebrick", "grey70"))
  )
)

# Figure S1
gene_vs_grade <- read_required(file.path(results_dir,
  "07_5c_OXPHOS_analysis/tables/OXPHOS13_gene_vs_grade.csv"))

gene_vs_grade <- gene_vs_grade %>%
  mutate(sig = factor(ifelse(FDR < 0.05, "FDR < 0.05", "n.s."),
                      levels = c("FDR < 0.05", "n.s.")))

gene_vs_grade$gene <- factor(gene_vs_grade$gene,
                              levels = gene_vs_grade$gene[order(gene_vs_grade$beta)])

p_s1 <- ggplot(gene_vs_grade, aes(x = beta, y = gene, fill = sig)) +
  geom_col(show.legend = TRUE) +
  geom_errorbarh(aes(xmin = beta - se, xmax = beta + se), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(
    values = c("FDR < 0.05" = "firebrick", "n.s." = "grey70"),
    drop = FALSE,
    name = "") +
  fdr_legend_guides +
  labs(x = "beta (per grade)", y = "",
       title = "S1  OXPHOS13 single-gene breakdown") +
  theme_pub
save_fig(p_s1, "FigureS1", w = 7, h = 5)

# Figure S2
bn_noage <- read_required(file.path(results_dir,
  "07_6a_bayesian_network/tables/bn_strength_no_age.csv"))
bn_noage_sig <- bn_noage %>% filter(strength >= 0.85) %>% arrange(desc(strength))
log_info(sprintf("Figure S2A: %d strong edges after excluding age", nrow(bn_noage_sig)))

p_s2a <- ggplot(bn_noage_sig,
                aes(x = strength, y = reorder(paste(from, "->", to), strength))) +
  geom_col(fill = "steelblue") +
  geom_vline(xintercept = 0.85, linetype = "dashed", color = "red") +
  labs(x = "Bootstrap strength", y = "",
       title = "A  Stable edges after excluding age") +
  theme_pub +
  scale_y_discrete(expand = expansion(add = c(0.5, 0.5))) +
  theme(axis.text.y = element_text(size = 9))

bn_grade_noage <- read_required(file.path(results_dir,
  "07_6a_bayesian_network/tables/bn_edges_into_grade_no_age.csv")) %>%
  arrange(desc(strength))
log_info(sprintf("Figure S2B: %d edges into grade after excluding age", nrow(bn_grade_noage)))

p_s2b <- ggplot(bn_grade_noage, aes(x = strength, y = reorder(from, strength))) +
  geom_col(fill = "grey60") +
  geom_vline(xintercept = 0.85, linetype = "dashed", color = "red") +
  labs(x = "Bootstrap strength", y = "",
       title = "B  Edges into AMD grade (all < 0.05)") +
  theme_pub +
  scale_y_discrete(expand = expansion(add = c(0.5, 0.5))) +
  theme(axis.text.y = element_text(size = 9))

fig_s2 <- p_s2a | p_s2b
save_fig(fig_s2, "FigureS2", w = 14, h = 6)

# Figure S3
sem <- read_required(file.path(results_dir,
  "07_6b_lavaan_mediation/tables/model_comparison.csv"))

sem$sig_ci <- factor(
  ifelse(sem$ci_low > 0 | sem$ci_high < 0, "Significant", "Non-significant"),
  levels = c("Significant", "Non-significant"))

p_s3 <- ggplot(sem, aes(x = estimate, y = model, color = sig_ci)) +
  geom_point(size = 2.5, show.legend = TRUE) +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(
    values = c("Significant" = "firebrick",
               "Non-significant" = "grey50"),
    drop = FALSE,
    name = "") +
  guides(color = guide_legend(
    override.aes = list(color = c("firebrick", "grey50"))
  )) +
  labs(x = "Standardized estimate", y = "",
       title = "S3  SEM mediation models") +
  theme_pub + theme(legend.position = "bottom")
save_fig(p_s3, "FigureS3", w = 8, h = 6)

# Figure S4
cmp <- read_required(file.path(results_dir,
  "07_7b_MR_plink/MR_cross_cohort_comparison.csv"))
forest_df <- cmp %>%
  select(symbol, OR_all, p_all, OR_early, p_early) %>%
  pivot_longer(cols = c(OR_all, OR_early),
               names_to = "outcome", values_to = "OR") %>%
  mutate(outcome = recode(outcome,
                           OR_all = "All-stage AMD (FinnGen)",
                           OR_early = "Early AMD (IAMDGC)"))

p_s4 <- ggplot(forest_df, aes(x = OR, y = reorder(symbol, OR), color = outcome)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10() +
  scale_color_manual(
    values = c("All-stage AMD (FinnGen)" = "#2166AC",
               "Early AMD (IAMDGC)" = "#B2182B"),
    name = "") +
  labs(x = "Odds ratio (log scale)", y = "",
       title = "S4  MR forest: OXPHOS genes vs AMD risk") +
  theme_pub + theme(legend.position = "bottom")
save_fig(p_s4, "FigureS4", w = 8, h = 6)

# Figure S5
ml <- read_required(file.path(results_dir,
  "07_8_ML_prediction_model/Table5c_external_AUC_all_cohorts.csv"))
p_s5 <- ggplot(ml, aes(x = cohort, y = auc, fill = model)) +
  geom_col(position = "dodge") +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                position = position_dodge(width = 0.9), width = 0.2) +
  geom_hline(yintercept = 0.5, linetype = "dashed") +
  labs(x = "", y = "AUC",
       title = "S5  ML performance across cohorts") +
  theme_pub + theme(axis.text.x = element_text(angle = 30, hjust = 1))
save_fig(p_s5, "FigureS5", w = 9, h = 5)

# Figure S6 -- retired (redundant with main Figure 6D); kept for record
# cell_prop <- read_required(file.path(results_dir,
#   "06_1_singlecell/tables/celltype_by_condition_ttest.csv"))
# cell_prop_long <- cell_prop %>%
#   select(celltype, mean_healthy, mean_damd) %>%
#   pivot_longer(-celltype, names_to = "condition", values_to = "mean") %>%
#   mutate(condition = recode(condition,
#                              mean_healthy = "Healthy",
#                              mean_damd = "dAMD"))
# p_s6 <- ggplot(cell_prop_long, aes(x = celltype, y = mean, fill = condition)) +
#   geom_col(position = "dodge") +
#   scale_fill_manual(values = c("Healthy" = "grey60", "dAMD" = "firebrick"),
#                     name = "") +
#   labs(x = "", y = "Mean proportion (%)",
#        title = "S6  Cell-type proportions (dAMD vs Healthy)") +
#   theme_pub + theme(axis.text.x = element_text(angle = 45, hjust = 1))
# save_fig(p_s6, "FigureS6", w = 9, h = 5)

# Figure S6 (was S7)
panel_ppf <- read_required(file.path(results_dir,
  "06_1_TRADE_vs_adjacent/tables/panel_perturbation_3modes.csv"))

panel_g3v2 <- panel_ppf %>%
  filter(Comparison == "3_vs_2") %>%
  arrange(desc(PPF))

panel_g3v2$Panel <- factor(panel_g3v2$Panel,
                            levels = panel_g3v2$Panel[order(panel_g3v2$PPF)])

p_s7 <- ggplot(panel_g3v2, aes(x = PPF, y = Panel, fill = direction)) +
  geom_col() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(
    values = c("Up" = "firebrick", "Down" = "steelblue"),
    name = "Direction") +
  labs(x = "Panel perturbation factor (PPF)",
       y = "",
       title = "S6  TRADE panel-level perturbation at G3 vs G2") +
  theme_pub
save_fig(p_s7, "FigureS6", w = 8, h = 6)

# Figure S7 (was S8)
mci <- read_required(file.path(results_dir,
  "06_2_BioTIP/tables/MCI_summary_across.csv"))
p_s8 <- ggplot(mci, aes(x = State, y = MCI, fill = Comparison)) +
  geom_col(position = "dodge") +
  labs(x = "State", y = "MCI",
       title = "S7  BioTIP MCI scores") +
  theme_pub
save_fig(p_s8, "FigureS7", w = 8, h = 5)

# Figure S8 (was S9)
perm_dist <- read_required(file.path(results_dir,
  "06_2_permutation_test/tables/permutation_dist_16genes.csv"))
perm_sum <- read_required(file.path(results_dir,
  "06_2_permutation_test/tables/permutation_summary.csv"))

obs_16 <- perm_sum %>% filter(object == "16_genes") %>% pull(real_p3_sig)

p_s9a <- ggplot(perm_dist, aes(x = n_p3_sig)) +
  geom_histogram(bins = 30, fill = "grey70", color = "black") +
  geom_vline(xintercept = obs_16, color = "red", linetype = "dashed", linewidth = 1) +
  coord_cartesian(xlim = c(0, 16)) +
  labs(x = "p3-significant genes (permutation)", y = "Frequency",
       title = "A  16-gene signature") +
  theme_pub
save_fig(p_s9a, "FigureS8", w = 7, h = 5)

# Figure S9 (was S10)
top_axis <- read_required(file.path(results_dir,
  "06_3_cellchat/tables/top_axis.csv"))
p_s10 <- ggplot(top_axis, aes(x = source, y = target, fill = n_interactions)) +
  geom_tile(color = "white") +
  geom_text(aes(label = n_interactions), size = 3) +
  scale_fill_gradient(low = "white", high = "firebrick") +
  labs(x = "", y = "",
       title = "S9  CellChat communication count (70,378 cells)") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_fig(p_s10, "FigureS9", w = 8, h = 6)

# Figure S10 (was S11)
grid <- read_required(file.path(results_dir,
  "07_7b_MR_plink/NDUFB3_sensitivity_grid.csv"))
grid_plot <- grid %>%
  filter(!is.na(OR_MR), outcome == "earlyAMD") %>%
  mutate(logOR = log(OR_MR))
p_s11 <- ggplot(grid_plot, aes(x = factor(r2_thr), y = factor(p_thr), fill = logOR)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", OR_MR)), size = 2.6) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, name = "log(OR)") +
  labs(x = "LD r2 threshold", y = "eQTL P-value threshold",
       title = "S10  NDUFB3 sensitivity grid (early AMD)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank())
save_fig(p_s11, "FigureS10", w = 8, h = 5)

# Figure S11 (was S12)
eqtl <- read_required(file.path(results_dir,
  "07_7a_NDUFB3/NDUFB3_multi_tissue_eQTL.csv")) %>%
  mutate(tissue_label = gsub("_", " ", tissue),
         neg_log10p = -log10(pval),
         gws = ifelse(pval < 5e-8, "Genome-wide significant", "Nominal")) %>%
  arrange(nes)
p_s12 <- ggplot(eqtl, aes(x = nes, y = reorder(tissue_label, nes))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(aes(size = neg_log10p, fill = gws),
             shape = 21, color = "black", stroke = 0.3) +
  scale_fill_manual(values = c("Genome-wide significant" = "firebrick",
                                "Nominal" = "grey70"), name = "") +
  scale_size_continuous(range = c(2, 6), name = expression(-log[10](P))) +
  labs(x = "Normalised effect size (NES)", y = "",
       title = "S11  NDUFB3 multi-tissue eQTL (GTEx V10)") +
  theme_pub + theme(legend.position = "bottom")
save_fig(p_s12, "FigureS11", w = 9, h = 7)

# Figure S12 (was S13)
spline <- read_required(file.path(results_dir,
  "02_8_G2G3_changepoint_robustness/tables/spline_ns_results.csv"))
bootstrap <- read_required(file.path(results_dir,
  "02_8_G2G3_changepoint_robustness/tables/bootstrap_p3_16genes.csv"))

p_s13a <- ggplot(spline, aes(x = reorder(gene, partial_R2), y = partial_R2)) +
  geom_col(fill = "steelblue") +
  labs(x = "", y = "Partial R²",
       title = "A  Natural spline (partial R²)") +
  theme_pub + theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_s13b <- ggplot(bootstrap, aes(x = reorder(gene, pct_sig_p3), y = pct_sig_p3)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
  labs(x = "", y = "p3 significance (%)",
       title = "B  Bootstrap stability (p3 significance %)") +
  theme_pub + theme(axis.text.x = element_text(angle = 45, hjust = 1))

fig_s13 <- p_s13a / p_s13b
save_fig(fig_s13, "FigureS12_G2G3_robustness", w = 9, h = 8)

# Figure S13 (was S14)
oxphos_cplx <- read_required(file.path(results_dir,
  "02_10_OXPHOS_variance_decomposition/tables/oxphos_by_complex.csv"))
oxphos_cplx$complex <- factor(oxphos_cplx$complex,
                               levels = oxphos_cplx$complex[order(oxphos_cplx$median_coef)])
p_s14 <- ggplot(oxphos_cplx, aes(x = complex, y = median_coef,
                                  fill = binom_p < 0.05)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(
    values = c("FALSE" = "grey70", "TRUE" = "firebrick"),
    breaks = c("FALSE", "TRUE"),
    labels = c("Not significant", "Significant"),
    name = "Binomial test") +
  labs(x = "OXPHOS complex", y = "Median coefficient",
       title = "S13  OXPHOS by complex") +
  theme_pub
save_fig(p_s14, "FigureS13_OXPHOS_variance", w = 7, h = 5)

# Figure S14 (was S15)
eqtm_bg <- read_required(file.path(results_dir,
  "02_9_eQTM_background_sensitivity/tables/eqtm_background_sensitivity.csv"))
eqtm_bg_long <- eqtm_bg %>%
  select(target, OR_20000, OR_14233) %>%
  pivot_longer(cols = c(OR_20000, OR_14233),
               names_to = "background", values_to = "OR") %>%
  mutate(background = recode(background,
                              OR_20000 = "20,000-gene background",
                              OR_14233 = "14,233-expressed-gene background")) %>%
  filter(OR > 0)

n_dropped <- sum(eqtm_bg$OR_20000 == 0 | eqtm_bg$OR_14233 == 0)
log_info(sprintf("Figure S14: %d rows dropped (OR = 0); %d rows plotted",
                 n_dropped, nrow(eqtm_bg_long)))

p_s15 <- ggplot(eqtm_bg_long, aes(x = target, y = OR, fill = background)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_y_log10() +
  labs(x = "", y = "Odds ratio (log scale)",
       title = "S14  eQTM background sensitivity",
       caption = "sig16 (OR = 0, no eQTM hit) is not shown on the log scale.") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        plot.caption = element_text(size = 7, color = "grey30", hjust = 0))
save_fig(p_s15, "FigureS14_eQTM_background", w = 8, h = 5)

# Figure S15 (was S16)
ko_gsea_rgc <- read_required(file.path(results_dir,
  "07_3_scTenifoldKnk_NDUFB3_RGC/tables/NDUFB3_GSEA_RGC.csv"))
ko_16_rgc <- read_required(file.path(results_dir,
  "07_3_scTenifoldKnk_NDUFB3_RGC/tables/NDUFB3_vs_16genes_RGC.csv"))

top_rgc <- ko_gsea_rgc %>% arrange(padj) %>% head(20)
top_rgc$pathway <- factor(top_rgc$pathway,
                           levels = top_rgc$pathway[order(top_rgc$NES)])
p_s16a <- ggplot(top_rgc, aes(x = NES, y = pathway, fill = NES > 0)) +
  geom_col() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(
    values = c("TRUE" = "firebrick", "FALSE" = "steelblue"),
    labels = c("TRUE" = "Up", "FALSE" = "Down"),
    name = "") +
  labs(x = "NES", y = "",
       title = "A  NDUFB3 KO GSEA (top 20) — RGC") +
  theme_pub

ko_16_rgc <- ko_16_rgc %>%
  mutate(sig = factor(ifelse(p.adj < 0.05, "FDR < 0.05", "n.s."),
                      levels = c("FDR < 0.05", "n.s.")))

p_s16b <- ggplot(ko_16_rgc, aes(x = Z, y = reorder(gene, Z), fill = sig)) +
  geom_col(show.legend = TRUE) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(
    values = c("FDR < 0.05" = "firebrick", "n.s." = "grey70"),
    drop = FALSE,
    name = "") +
  fdr_legend_guides +
  labs(x = "Z-score", y = "",
       title = "B  16-gene candidate set perturbation — RGC") +
  theme_pub

ko_gsea_rod <- read_required(file.path(results_dir,
  "07_3_scTenifoldKnk_NDUFB3_ROD/tables/NDUFB3_GSEA_Rod.csv"))
ko_16_rod <- read_required(file.path(results_dir,
  "07_3_scTenifoldKnk_NDUFB3_ROD/tables/NDUFB3_vs_16genes_Rod.csv"))

top_rod <- ko_gsea_rod %>% arrange(padj) %>% head(20)
top_rod$pathway <- factor(top_rod$pathway,
                           levels = top_rod$pathway[order(top_rod$NES)])
p_s16c <- ggplot(top_rod, aes(x = NES, y = pathway, fill = NES > 0)) +
  geom_col() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(
    values = c("TRUE" = "firebrick", "FALSE" = "steelblue"),
    labels = c("TRUE" = "Up", "FALSE" = "Down"),
    name = "") +
  labs(x = "NES", y = "",
       title = "C  NDUFB3 KO GSEA (top 20) — ROD") +
  theme_pub

ko_16_rod <- ko_16_rod %>%
  mutate(sig = factor(ifelse(p.adj < 0.05, "FDR < 0.05", "n.s."),
                      levels = c("FDR < 0.05", "n.s.")))

p_s16d <- ggplot(ko_16_rod, aes(x = Z, y = reorder(gene, Z), fill = sig)) +
  geom_col(show.legend = TRUE) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(
    values = c("FDR < 0.05" = "firebrick", "n.s." = "grey70"),
    drop = FALSE,
    name = "") +
  fdr_legend_guides +
  labs(x = "Z-score", y = "",
       title = "D  16-gene candidate set perturbation — ROD") +
  theme_pub

fig_s16 <- (p_s16a | p_s16b) / (p_s16c | p_s16d)
save_fig(fig_s16, "FigureS15_NDUFB3_KO_RGC_ROD", w = 14, h = 12)

# Figure S16 (was S17)
atac <- read_required(file.path(results_dir,
  "07_1_atac_validation/tables/gene_set_ATAC_summary.csv"))

atac <- atac %>%
  mutate(sig = factor(ifelse(wilcox_p_zero < 0.05, "wilcox p<0.05", "n.s."),
                      levels = c("wilcox p<0.05", "n.s.")))

p_s17 <- ggplot(atac, aes(x = reorder(gene_set, cohen_d), y = cohen_d,
                           fill = sig)) +
  geom_col(show.legend = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(
    values = c("wilcox p<0.05" = "firebrick", "n.s." = "grey70"),
    drop = FALSE,
    name = "") +
  guides(fill = guide_legend(
    override.aes = list(fill = c("firebrick", "grey70"))
  )) +
  facet_wrap(~ tissue, ncol = 2, scales = "free_y") +
  coord_flip() +
  labs(x = "", y = "Cohen's d",
       title = "S16  ATAC-seq promoter accessibility") +
  theme_pub
save_fig(p_s17, "FigureS16a_ATAC_cohen_d", w = 12, h = 8)

atac_rna <- read_required(file.path(results_dir,
  "07_1_atac_validation/tables/ATAC_vs_RNA_p3_consistency.csv"))

p_s17b <- ggplot(atac_rna, aes(x = tissue, y = spearman_r, fill = tissue)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(aes(label = sprintf("r=%.2f\npct_concord=%.1f%%",
                                 spearman_r, pct_same_dir)),
            vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("RPE" = "#2166AC", "Retina" = "#B2182B"),
                    guide = "none") +
  labs(x = "Tissue", y = "Spearman r (ATAC vs RNA)",
       title = "S16b  ATAC vs RNA consistency (tissue-level)") +
  coord_cartesian(ylim = c(0, max(atac_rna$spearman_r) * 1.3)) +
  theme_pub
save_fig(p_s17b, "FigureS16b_ATAC_vs_RNA", w = 6, h = 5)

# Figure S17 (was S18)
composite_perSD <- read_required(file.path(results_dir,
  "07_5c_OXPHOS_analysis/tables/OXPHOS13_composite_perSD.csv"))
pathway_perSD <- read_required(file.path(results_dir,
  "07_5c_OXPHOS_analysis/tables/OXPHOS_pathway_perSD.csv"))

perSD_df <- bind_rows(
  composite_perSD %>% mutate(score = "13-gene composite"),
  pathway_perSD   %>% mutate(score = "OXPHOS pathway (187 genes)")
) %>%
  filter(model == "adjusted") %>%
  mutate(
    score = factor(score,
                   levels = c("13-gene composite",
                              "OXPHOS pathway (187 genes)")),
    sig_group = factor(
      ifelse(p_value < 0.05, "Significant (p < 0.05)", "Not significant"),
      levels = c("Significant (p < 0.05)", "Not significant")),
    label = sprintf("beta = %.3f (95%% CI %.3f to %.3f), p = %.3f",
                    beta_per_SD, ci_lower, ci_upper, p_value)
  )

p_s18 <- ggplot(perSD_df, aes(x = beta_per_SD, y = score, color = sig_group)) +
  geom_point(size = 4, show.legend = TRUE) +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey30") +
  geom_text(aes(label = label),
            vjust = -1.2, size = 3, color = "black") +
  scale_color_manual(
    values = c("Significant (p < 0.05)" = "firebrick",
               "Not significant" = "grey50"),
    drop = FALSE,
    name = "") +
  guides(color = guide_legend(
    override.aes = list(color = c("firebrick", "grey50"))
  )) +
  labs(x = "Beta per 1 SD (adjusted for age, sex, RIN, PMI)", y = "",
       title = "S17  Per-1-SD observational effect sizes",
       caption = "Both estimates are directionally negative but not statistically significant, consistent with the MR null results.") +
  theme_pub +
  theme(legend.position = "bottom",
        plot.caption = element_text(size = 8, color = "grey30", hjust = 0))
save_fig(p_s18, "FigureS17_per_SD_effect_sizes", w = 8, h = 4)

# Figure S18 (was S19)
ssgsea_file <- file.path(results_dir,
  "05_3_module_de/data/module_scores_ssgsea.rds")
clin_file <- file.path(results_dir,
  "01_1b_batch_enhance/data/clinical_wgcna.rds")
if (!file.exists(ssgsea_file) || !file.exists(clin_file)) {
  stop("Figure S18: ssGSEA or clinical file missing")
}
ssgsea <- readRDS(ssgsea_file)
clin <- readRDS(clin_file)

ssgsea_df <- as.data.frame(t(ssgsea)) %>%
  tibble::rownames_to_column("sample_id") %>%
  left_join(clin %>% select(sample_id, grade), by = "sample_id")

gy_by_grade <- ssgsea_df %>%
  group_by(grade) %>%
  summarise(mean_ssgsea = mean(greenyellow, na.rm = TRUE),
            sd_ssgsea = sd(greenyellow, na.rm = TRUE),
            n = n(),
            .groups = "drop")

log_info("greenyellow ssGSEA means by grade:")
print(gy_by_grade)

p_s19 <- ggplot(gy_by_grade, aes(x = factor(grade), y = mean_ssgsea,
                                  fill = factor(grade))) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = mean_ssgsea - sd_ssgsea / sqrt(n),
                    ymax = mean_ssgsea + sd_ssgsea / sqrt(n)),
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("#2166AC", "#67A9CF", "#EF8A62", "#B2182B"),
                    guide = "none") +
  labs(x = "AMD grade", y = "greenyellow ssGSEA (mean ± SE)",
       title = "S18  greenyellow module ssGSEA by AMD grade") +
  theme_pub
save_fig(p_s19, "FigureS18_greenyellow_ssGSEA", w = 6, h = 4)

# Figure S19 (was S20)
perm_df <- read_required(file.path(results_dir,
  "06_3_cellchat/tables/RPE_vs_RodCone_permutation.csv"))
final_df <- read_required(file.path(results_dir,
  "06_3_cellchat/tables/RPE_vs_RodCone_final.csv"))

plot_df <- perm_df %>%
  mutate(
    sig = case_when(
      FDR < 0.05 ~ "FDR < 0.05",
      FDR < 0.10 ~ "FDR < 0.10",
      TRUE ~ "NS"
    ),
    neglog10fdr = -log10(pmax(FDR, 1e-10))
  )

p_s20a <- ggplot(plot_df, aes(x = log2_ratio, y = neglog10fdr, color = sig)) +
  geom_point(aes(size = pmax(n_RPE + n_RodCone, 1)), alpha = 0.7) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("FDR < 0.05" = "#C44E52",
                                 "FDR < 0.10" = "#F0AD4E",
                                 "NS" = "grey70")) +
  scale_size_continuous(name = "Total pairs", range = c(2, 8)) +
  coord_cartesian(xlim = c(-3, 3)) +
  theme_pub +
  labs(x = "log2(RPE / RodCone)", y = "-log10(FDR)",
       title = "A  RPE vs Rod/Cone pathway usage (permutation test)",
       color = "") +
  theme(legend.position = "bottom")

forest_df <- final_df %>%
  filter(!is.na(boot_ci_low) & !is.na(boot_ci_high)) %>%
  arrange(desc(log2_ratio)) %>%
  head(20)
forest_df$pathway <- factor(forest_df$pathway,
                             levels = forest_df$pathway[order(forest_df$log2_ratio)])

p_s20b <- ggplot(forest_df, aes(x = log2_ratio, y = pathway)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = boot_ci_low, xmax = boot_ci_high),
                 height = 0.2, color = "grey40") +
  geom_point(aes(color = FDR < 0.05), size = 3) +
  scale_color_manual(values = c("TRUE" = "#C44E52", "FALSE" = "grey60"),
                     name = "FDR < 0.05") +
  scale_x_continuous(breaks = c(-15, -10, -5, 0, 5, 10, 15)) +
  theme_pub +
  labs(x = "log2(RPE / RodCone)", y = "",
       title = "B  Top pathways (bootstrap 95% CI)") +
  theme(legend.position = "bottom")

fig_s20 <- p_s20a / p_s20b +
  plot_layout(heights = c(1, 1.2))
save_fig(fig_s20, "FigureS19_CellChat_RPE", w = 9, h = 10)

# Figure S20 (was S21)
ko_dir_rpe <- file.path(results_dir, "07_3_scTenifoldKnk_NDUFB3_RPE")

ko_gsea_rpe  <- read_required(file.path(ko_dir_rpe,
  "tables/NDUFB3_GSEA_RPE.csv"))
ko_16_rpe    <- read_required(file.path(ko_dir_rpe,
  "tables/NDUFB3_vs_16genes_RPE.csv"))
ko_ci_rpe    <- read_required(file.path(ko_dir_rpe,
  "tables/NDUFB3_vs_ComplexI_RPE.csv"))

top_rpe <- ko_gsea_rpe %>% arrange(padj) %>% head(20)
top_rpe$pathway <- factor(top_rpe$pathway,
                           levels = top_rpe$pathway[order(top_rpe$NES)])
p_s21a <- ggplot(top_rpe, aes(x = NES, y = pathway, fill = NES > 0)) +
  geom_col() +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(
    values = c("TRUE" = "firebrick", "FALSE" = "steelblue"),
    labels = c("TRUE" = "Up", "FALSE" = "Down"),
    name = "") +
  labs(x = "NES", y = "",
       title = "A  NDUFB3 KO GSEA (top 20) — RPE") +
  theme_pub

ko_16_rpe <- ko_16_rpe %>%
  mutate(sig = factor(ifelse(p.adj < 0.05, "FDR < 0.05", "n.s."),
                      levels = c("FDR < 0.05", "n.s.")))

p_s21b <- ggplot(ko_16_rpe, aes(x = Z, y = reorder(gene, Z), fill = sig)) +
  geom_col(show.legend = TRUE) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(
    values = c("FDR < 0.05" = "firebrick", "n.s." = "grey70"),
    drop = FALSE,
    name = "") +
  fdr_legend_guides +
  labs(x = "Z-score", y = "",
       title = "B  16-gene candidate set perturbation — RPE") +
  theme_pub

ko_ci_rpe <- ko_ci_rpe %>%
  mutate(sig = factor(ifelse(p.adj < 0.05, "FDR < 0.05", "n.s."),
                      levels = c("FDR < 0.05", "n.s."))) %>%
  arrange(desc(abs(Z))) %>%
  head(20)

p_s21c <- ggplot(ko_ci_rpe, aes(x = Z, y = reorder(gene, Z), fill = sig)) +
  geom_col(show.legend = TRUE) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(
    values = c("FDR < 0.05" = "firebrick", "n.s." = "grey70"),
    drop = FALSE,
    name = "") +
  fdr_legend_guides +
  labs(x = "Z-score", y = "",
       title = "C  Complex I perturbation (top 20) — RPE") +
  theme_pub

fig_s21 <- p_s21a | (p_s21b / p_s21c)
save_fig(fig_s21, "FigureS20_NDUFB3_KO_RPE", w = 14, h = 9)

# Figure S21 (was S22)
island <- read_required(file.path(results_dir,
  "07_10c_cpg_island_distribution/tables/cpg_island_distribution.csv"))

island <- island %>%
  mutate(
    relation_island = factor(relation_island,
      levels = c("Island", "N_Shore", "S_Shore",
                 "N_Shelf", "S_Shelf", "OpenSea"))
  )

island_long <- island %>%
  select(relation_island, pct_target, pct_bg) %>%
  pivot_longer(cols = c(pct_target, pct_bg),
               names_to = "set", values_to = "pct") %>%
  mutate(set = recode(set,
    pct_target = "greenyellow eQTM CpGs (n = 101)",
    pct_bg     = "All eQTM CpGs (background)"))

p_s22a <- ggplot(island_long,
                 aes(x = relation_island, y = pct, fill = set)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("greenyellow eQTM CpGs (n = 101)" = "#7FBC41",
                               "All eQTM CpGs (background)" = "grey70"),
                    name = "") +
  labs(x = "Relation to CpG island", y = "Percentage of CpGs (%)",
       title = "A  CpG island distribution") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "bottom")

hit <- read_required(file.path(results_dir,
  "07_10c_cpg_island_distribution/tables/motif_hit_rate_descriptive.csv"))

hit_top <- hit %>%
  arrange(desc(fold_neg)) %>%
  head(12) %>%
  select(motif, pct_neg, pct_pos, pct_bg) %>%
  pivot_longer(cols = c(pct_neg, pct_pos, pct_bg),
               names_to = "set", values_to = "pct") %>%
  mutate(set = recode(set,
    pct_neg = "Negative-beta CpGs (n = 55)",
    pct_pos = "Positive-beta CpGs (n = 46)",
    pct_bg  = "Background (n = 2000)"))

hit_top$motif <- factor(hit_top$motif,
                         levels = hit$motif[order(-hit$fold_neg)][1:12])

p_s22b <- ggplot(hit_top, aes(x = motif, y = pct, fill = set)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("Negative-beta CpGs (n = 55)" = "#B2182B",
                               "Positive-beta CpGs (n = 46)" = "#2166AC",
                               "Background (n = 2000)" = "grey70"),
                    name = "") +
  labs(x = "OXPHOS-related motif (JASPAR2024 CORE)",
       y = "CpGs with motif hit (%)",
       title = "B  Descriptive motif hit rates (not enrichment)") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "bottom")

fig_s22 <- p_s22a / p_s22b +
  plot_layout(heights = c(1, 1.2)) &
  theme(plot.margin = margin(6, 8, 6, 10))

save_fig(fig_s22, "FigureS21_CpG_island_and_motif_hit", w = 10, h = 10)

# Figure S22 (was S23)
motif_primary <- read_required(file.path(results_dir,
  "07_1b_atac_motif_enrichment/tables/motif_enrichment_primary.csv"))
motif_robust <- read_required(file.path(results_dir,
  "07_1b_atac_motif_enrichment/tables/motif_robust_summary.csv"))

top_m <- head(motif_primary[order(motif_primary$FDR), ], 20)
top_m$motif_name <- factor(top_m$motif_name, levels = rev(top_m$motif_name))
top_m$sig <- factor(ifelse(top_m$FDR < 0.05, "FDR<0.05", "n.s."),
                    levels = c("FDR<0.05", "n.s."))

p_s23a <- ggplot(top_m, aes(x = motif_name, y = -log10(FDR), fill = sig)) +
  geom_col(show.legend = TRUE) + coord_flip() + theme_pub +
  scale_fill_manual(values = c("FDR<0.05" = "firebrick", "n.s." = "grey70"),
                    drop = FALSE, name = "") +
  guides(fill = guide_legend(
    override.aes = list(fill = c("firebrick", "grey70"))
  )) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(x = "", y = "-log10(FDR)",
       title = "A  Motif enrichment, fixed background (primary)")

sens_top <- head(motif_robust[order(-motif_robust$freq_sig), ], 20)
sens_top$motif_name <- factor(sens_top$motif_name,
                              levels = rev(sens_top$motif_name))

p_s23b <- ggplot(sens_top, aes(x = motif_name, y = freq_sig)) +
  geom_col(fill = "steelblue") + coord_flip() + theme_pub +
  geom_hline(yintercept = 0.80, linetype = "dashed", color = "firebrick") +
  labs(x = "", y = "Fraction of 100 random backgrounds with FDR < 0.05",
       title = "B  Sensitivity to background choice")

fig_s23 <- p_s23a | p_s23b
save_fig(fig_s23, "FigureS22_greenyellow_motif_enrichment", w = 12, h = 6)

# Figure S23 (was S24): Conventional two-group contrasts vs poly3_F
log_info("--- Figure S23: conventional vs poly3_F ---")

de_multi <- read_required(file.path(results_dir,
  "02_1_differential_expression_main/tables/limma_multi_full.csv"))
poly3f <- read_required(file.path(results_dir,
  "02_1_453_discovery/tables/genes16_poly3F_adjusted.csv"))
sig16  <- read_required(file.path(results_dir,
  "02_1_453_discovery/tables/sig16_genes.csv"))

de_plot <- de_multi %>%
  filter(contrast %in% c("g4_vs_g1", "g3_vs_g2", "g2_vs_g1")) %>%
  mutate(
    contrast = recode(contrast,
      g4_vs_g1 = "G4 vs G1",
      g3_vs_g2 = "G3 vs G2",
      g2_vs_g1 = "G2 vs G1"),
    neglog10fdr = -log10(pmax(FDR, 1e-10))
  )

cat("Panel A rows:", nrow(de_plot), "\n")

p_s24a <- ggplot(de_plot, aes(x = logFC, y = neglog10fdr)) +
  geom_point(color = "grey70", alpha = 0.5, size = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  facet_wrap(~ contrast, ncol = 3) +
  labs(x = "logFC", y = "-log10(FDR)",
       title = "A  Conventional two-group contrasts") +
  theme_pub +
  theme(strip.text = element_text(face = "bold"))

poly3f_p3 <- poly3f %>%
  filter(model == "M1_age_PMI_RIN") %>%
  mutate(
    neglog10fdr = -log10(pmax(p3_FDR, 1e-10)),
    is_16gene = gene %in% sig16$gene
  )

cat("Panel B rows:", nrow(poly3f_p3), "\n")

p_s24b <- ggplot(poly3f_p3, aes(x = p3_logFC, y = neglog10fdr)) +
  geom_point(data = poly3f_p3 %>% filter(!is_16gene),
             color = "grey70", alpha = 0.5, size = 0.8) +
  geom_point(data = poly3f_p3 %>% filter(is_16gene),
             color = "firebrick", size = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
  labs(x = "p3 coefficient (logFC)", y = "-log10(FDR)",
       title = "B  poly3_F non-linear detection (p3 term)") +
  theme_pub

fig_s24 <- p_s24a / p_s24b +
  plot_layout(heights = c(1, 1.2)) &
  theme(plot.margin = margin(6, 8, 6, 10))

save_fig(fig_s24, "FigureS23_conventional_vs_poly3F", w = 12, h = 8)

log_info("=== 10_1C v10.7 done ===")

