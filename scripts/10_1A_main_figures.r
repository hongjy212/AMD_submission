# 10_1A_main_figures.R (v4.11)
# Main figures 1-7
# v4.11: Figure 2 restructured (2A left/right; 2C removed);
#        greenyellow permutation moved to Figure 4B (4B/4C/4D -> 4C/4D/4E)

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr); library(tidyr); library(readr)
  library(patchwork); library(pheatmap); library(tibble)
  library(grid); library(gridExtra)
})

results_dir <- DIR_RESULT
fig_dir     <- file.path(results_dir, "10_1A_main_figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
log_info("=== 10_1A Main figures (v4.11) ===")

read_required <- function(path) {
  if (!file.exists(path)) stop(sprintf("REQUIRED FILE MISSING: %s", path))
  read_csv(path, show_col_types = FALSE)
}

check_cols <- function(df, required, name) {
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    stop(sprintf("%s: missing columns: %s", name, paste(missing, collapse = ", ")))
  }
  invisible(TRUE)
}

theme_pub <- theme_classic(base_size = 11) +
  theme(
    axis.text       = element_text(color = "black"),
    axis.title      = element_text(face = "bold"),
    plot.title      = element_text(face = "bold", size = 11),
    plot.margin     = margin(6, 8, 6, 10),
    axis.text.y     = element_text(margin = margin(r = 4)),
    axis.text.x     = element_text(margin = margin(t = 4))
  )

save_fig <- function(p, name, w = 7, h = 5, dpi = 300) {
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p,
         width = w, height = h, device = "pdf")
  ggsave(file.path(fig_dir, paste0(name, ".tiff")), p,
         width = w, height = h, dpi = dpi, compression = "lzw")
  log_info(sprintf("saved figure: %s (%gx%g)", name, w, h))
}

normalize_complex <- function(x) {
  x <- toupper(trimws(x))
  x <- gsub("^COMPLEX\\s*", "", x)
  x <- gsub("^C\\s*", "", x)
  x <- trimws(x)
  num_to_roman <- c("1" = "I", "2" = "II", "3" = "III", "4" = "IV", "5" = "V")
  x <- ifelse(x %in% names(num_to_roman), num_to_roman[x], x)
  paste0("Complex ", x)
}

pathway_labels <- c(
  "Mito_CI_Subunits"   = "Complex I subunits",
  "Mito_CV_Subunits"   = "Complex V subunits",
  "Mito_CIV_Subunits"  = "Complex IV subunits",
  "Mito_CIII_Subunits" = "Complex III subunits",
  "Mito_Complex_I"     = "Complex I",
  "Mito_Complex_II"    = "Complex II",
  "Mito_Complex_III"   = "Complex III",
  "Mito_Complex_IV"    = "Complex IV",
  "Mito_Complex_V"     = "Complex V",
  "OXPHOS"             = "OXPHOS",
  "Mito_TCA"           = "TCA cycle",
  "Mito_TCA_associated"= "TCA-associated",
  "Mito_Calcium"       = "Calcium",
  "Mito_CalciumCycle"  = "Calcium cycle",
  "Fatty_acid_oxidation" = "Fatty acid oxidation",
  "Mito_FattyAcidOxidation" = "Mito fatty acid oxidation",
  "Mito_Chaperones"    = "Chaperones",
  "Mito_FeS_Cluster"   = "Fe-S cluster",
  "Mito_Ribosome"      = "Mito ribosome",
  "Mito_antioxidant"   = "Mito antioxidant",
  "Mito_ROS_Glutathione" = "ROS/glutathione",
  "Mito_Dynamics_Surveillance" = "Dynamics/surveillance",
  "Mito_mitophagy"     = "Mitophagy",
  "Mito_biogenesis"    = "Mito biogenesis",
  "Mito_dynamics"      = "Mito dynamics",
  "Mito_mtDNA_rep"     = "mtDNA replication",
  "Mito_mtDNA_replication" = "mtDNA replication"
)

bulk_cplx <- read_required(file.path(results_dir,
  "07_5c_OXPHOS_analysis/tables/OXPHOS13_complex_score.csv")) %>%
  mutate(complex_std = normalize_complex(complex))

rgc_cplx <- read_required(file.path(results_dir,
  "07_3_NDUFB8_coexpression/tables/RGC_complex_100iter_summary.csv")) %>%
  mutate(complex_std = normalize_complex(complex))

cross_complex <- full_join(
  bulk_cplx %>% select(complex_std, bulk_n_subunits = n_subunits,
                       bulk_beta = beta, bulk_p = p, bulk_FDR = FDR),
  rgc_cplx %>% select(complex_std, rgc_n_genes = n_genes,
                      rgc_log2FC = log2FC_median,
                      rgc_IQR_lo = log2FC_IQR_lo,
                      rgc_IQR_hi = log2FC_IQR_hi,
                      pct_p_lt_005, pct_log2FC_pos),
  by = "complex_std"
) %>%
  mutate(
    complex_std = factor(complex_std,
                         levels = c("Complex I", "Complex II",
                                    "Complex III", "Complex IV", "Complex V")),
    bulk_direction = ifelse(bulk_beta < 0, "down", "up"),
    rgc_direction  = ifelse(rgc_log2FC < 0, "down", "up"),
    concordant     = bulk_direction == rgc_direction
  ) %>%
  arrange(complex_std)

# ============================================================
# Figure 1
# ============================================================
p1a <- ggplot() +
  annotate("rect", xmin = 0.02, xmax = 0.32, ymin = 0.55, ymax = 0.95,
           fill = "#E8F1F8", color = "grey40") +
  annotate("text", x = 0.17, y = 0.85, label = "Bulk retina\n(n = 419)",
           size = 3.2, fontface = "bold") +
  annotate("text", x = 0.17, y = 0.65,
           label = "GSE115828\nG1-G4", size = 2.8) +
  annotate("rect", xmin = 0.35, xmax = 0.65, ymin = 0.55, ymax = 0.95,
           fill = "#F8E8E8", color = "grey40") +
  annotate("text", x = 0.50, y = 0.85, label = "Single-cell retina\n(n = 70,394)",
           size = 3.2, fontface = "bold") +
  annotate("text", x = 0.50, y = 0.65,
           label = "GSE221042\n6 Healthy / 4 dAMD", size = 2.8) +
  annotate("rect", xmin = 0.68, xmax = 0.98, ymin = 0.55, ymax = 0.95,
           fill = "#F0F8E8", color = "grey40") +
  annotate("text", x = 0.83, y = 0.85, label = "DNA methylation\n+ eQTL + GWAS",
           size = 3.2, fontface = "bold") +
  annotate("text", x = 0.83, y = 0.65,
           label = "GSE231536\neQTLGen / GTEx\nFinnGen / GCST010723",
           size = 2.8, lineheight = 0.9) +
  annotate("segment", x = 0.17, xend = 0.17, y = 0.55, yend = 0.40,
           arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8) +
  annotate("segment", x = 0.50, xend = 0.50, y = 0.55, yend = 0.40,
           arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8) +
  annotate("segment", x = 0.83, xend = 0.83, y = 0.55, yend = 0.40,
           arrow = arrow(length = unit(0.2, "cm")), linewidth = 0.8) +
  annotate("rect", xmin = 0.02, xmax = 0.98, ymin = 0.05, ymax = 0.38,
           fill = "#F7F7F7", color = "grey40") +
  annotate("text", x = 0.50, y = 0.30,
           label = "poly3_F  |  WGCNA  |  ssGSEA  |  UCell  |  cellGeometry",
           size = 3.0, fontface = "bold") +
  annotate("text", x = 0.50, y = 0.18,
           label = "Bayesian network  |  SEM mediation  |  Mendelian randomization  |  scTenifoldKnk",
           size = 2.8, fontface = "bold") +
  annotate("text", x = 0.50, y = 0.07,
           label = "Multi-scale causal inference framework",
           size = 3.2, fontface = "bold.italic", color = "firebrick") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  labs(title = "A  Study design") +
  theme_void(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5))

cmp_all <- read_required(file.path(results_dir,
  "07_7b_MR_plink/MR_cross_cohort_comparison.csv"))
check_cols(cmp_all, c("symbol", "OR_all", "p_all", "OR_early", "p_early"),
           "MR_cross_cohort_comparison")

comp_genes <- c("CFH", "C3", "CD59")
comp_plot_b <- cmp_all %>%
  filter(symbol %in% comp_genes) %>%
  select(symbol, OR_all, p_all, OR_early, p_early) %>%
  pivot_longer(cols = c(OR_all, OR_early),
               names_to = "outcome", values_to = "OR") %>%
  mutate(outcome = recode(outcome,
                          OR_all   = "All-stage AMD (FinnGen)",
                          OR_early = "Early AMD (IAMDGC)"))

comp_fdr <- read_required(file.path(results_dir,
  "07_7e_MR_unified_fdr/tables/complement_positive_control_unified_fdr.csv"))
check_cols(comp_fdr, c("symbol", "outcome", "p_unified_fdr"),
           "complement_positive_control_unified_fdr")
cfh_fdr <- comp_fdr %>% filter(symbol == "CFH", outcome == "earlyAMD") %>%
  pull(p_unified_fdr)
c3_p_all <- comp_plot_b %>% filter(symbol == "C3") %>% pull(p_all) %>% unique()
log_info(sprintf("Figure 1B: CFH earlyAMD unified FDR = %.4f", cfh_fdr))
log_info(sprintf("Figure 1B: C3 all-stage nominal p = %.4f", c3_p_all[1]))

p1b <- ggplot(comp_plot_b,
              aes(x = OR, y = symbol, color = outcome, shape = outcome)) +
  geom_point(size = 4, position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey30") +
  scale_x_log10() +
  scale_color_manual(values = c("All-stage AMD (FinnGen)" = "#2166AC",
                                "Early AMD (IAMDGC)" = "#B2182B"),
                     name = "") +
  scale_shape_manual(values = c("All-stage AMD (FinnGen)" = 16,
                                "Early AMD (IAMDGC)" = 17),
                     name = "") +
  labs(x = "Odds ratio (log scale)", y = "",
       title = "B  Complement positive controls",
       caption = sprintf("CFH early AMD unified FDR = %.4f; C3 all-stage nominal p = %.4f",
                         cfh_fdr, c3_p_all[1])) +
  theme_pub +
  scale_y_discrete(expand = expansion(add = c(0.5, 0.5))) +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 10, face = "bold"),
        plot.caption = element_text(size = 7, color = "grey30", hjust = 0))

fig1 <- p1a | p1b
save_fig(fig1, "Figure1", w = 12, h = 5.5)

# ============================================================
# Figure 2
# v4.11: 2A = left (permutation) + right (16-gene p3 signal);
#        2C moved to Figure 4B; layout = 2A / (2B | 2D)
# ============================================================
perm_dist <- read_required(file.path(results_dir,
  "06_2_permutation_test/tables/permutation_dist_16genes.csv"))
perm_sum  <- read_required(file.path(results_dir,
  "06_2_permutation_test/tables/permutation_summary.csv"))
trade     <- read_required(file.path(results_dir,
  "06_1_TRADE_vs_adjacent/tables/TRADE_impact_3modes.csv"))
perm_mod  <- read_required(file.path(results_dir,
  "06_2_permutation_test/tables/permutation_dist_modules.csv"))

obs_16 <- perm_sum %>% filter(object == "16_genes") %>% pull(real_p3_sig)
obs_gy <- perm_sum %>% filter(object == "greenyellow") %>% pull(real_p3_sig)
log_info(sprintf("Figure 2A: obs_16 = %d", obs_16))
log_info(sprintf("Figure 4B: obs_gy = %d", obs_gy))

# ---- Figure 2A left: permutation histogram ----
p2a_left <- ggplot(perm_dist, aes(x = n_p3_sig)) +
  geom_histogram(bins = 30, fill = "grey70", color = "black") +
  geom_vline(xintercept = obs_16, color = "red",
             linetype = "dashed", linewidth = 1) +
  coord_cartesian(xlim = c(0, max(perm_dist$n_p3_sig, obs_16) + 1)) +
  labs(x = "p3-significant genes (permutation)", y = "Frequency",
       title = "A  Permutation test: 16-gene signature") +
  theme_pub

# ---- Figure 2A right: per-gene p3 signal for the 16 candidates ----
gene16_p3 <- read_required(file.path(results_dir,
  "02_1_453_discovery/tables/genes16_poly3F_adjusted.csv")) %>%
  filter(model == "M1_age_PMI_RIN") %>%
  mutate(neglog10fdr = -log10(pmax(p3_FDR, 1e-10)),
         gene = factor(gene, levels = gene[order(p3_FDR, decreasing = TRUE)]))

log_info(sprintf("Figure 2A right: %d genes, max -log10(p3 FDR) = %.2f",
                 nrow(gene16_p3), max(gene16_p3$neglog10fdr)))

p2a_right <- ggplot(gene16_p3, aes(x = neglog10fdr, y = gene)) +
  geom_col(fill = "firebrick", width = 0.7) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "grey30") +
  labs(x = expression(-log[10]*"(p3 FDR)"), y = "",
       title = "Per-gene p3 signal (16-gene signature)") +
  theme_pub +
  theme(axis.text.y = element_text(size = 8))

p2a_combined <- p2a_left | p2a_right

# ---- Figure 2B: TRADE ----
p2b <- ggplot(trade, aes(x = Comparison, y = TranscriptomeImpact, fill = Design)) +
  geom_col(position = "dodge") +
  labs(x = "", y = "Transcriptome impact",
       title = "B  TRADE: G3 vs G2 shows the highest impact under adjusted and pure designs (total design peaks at G4 vs G3)") +
  theme_pub +
  theme(legend.position = "bottom")

# ---- Figure 2C: robustness panels (C1-C4) ----
spline <- read_required(file.path(results_dir,
  "02_8_G2G3_changepoint_robustness/tables/spline_ns_results.csv"))
bootstrap <- read_required(file.path(results_dir,
  "02_8_G2G3_changepoint_robustness/tables/bootstrap_p3_16genes.csv"))
loo <- read_required(file.path(results_dir,
  "02_8_G2G3_changepoint_robustness/tables/FCF1_leave_one_out.csv"))
batch <- read_required(file.path(results_dir,
  "02_13_batch_correction_sensitivity/tables/batch_correction_sensitivity.csv"))

p2c_spline <- ggplot(spline, aes(x = reorder(gene, partial_R2), y = partial_R2)) +
  geom_col(fill = "steelblue") +
  labs(x = "", y = "Partial R²", title = "C1  Natural spline") +
  theme_pub + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

p2c_boot <- ggplot(bootstrap, aes(x = reorder(gene, pct_sig_p3), y = pct_sig_p3)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
  labs(x = "", y = "p3 significance (%)", title = "C2  Bootstrap stability") +
  theme_pub + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

p2c_loo <- ggplot(loo, aes(x = reorder(gene, pct_loo_sig), y = pct_loo_sig)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
  labs(x = "", y = "LOO significance (%)", title = "C3  Leave-one-out") +
  theme_pub + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

batch <- batch %>%
  mutate(method = recode(method,
    "A_limma_protect_grade"    = "limma (grade protected)",
    "B_limma_no_protection"    = "limma (no protection)",
    "C_ComBat_grade_covariate" = "ComBat",
    "D_current_published"      = "Main pipeline"))

p2c_batch <- ggplot(batch, aes(x = reorder(method, beta_ComplexIII),
                               y = beta_ComplexIII,
                               fill = p_ComplexIII < 0.05)) +
  geom_col() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("grey70", "firebrick"), name = "p<0.05") +
  labs(x = "", y = "beta ComplexIII",
       title = "C4  Batch correction sensitivity") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
        legend.position = "none")

p2c <- (p2c_spline | p2c_boot) / (p2c_loo | p2c_batch)

# ---- Assemble Figure 2: 2A / 2B / 2C ----
fig2 <- p2a_combined /
        p2b /
        p2c +
  plot_layout(heights = c(1, 0.7, 1.6)) &
  theme(plot.margin = margin(6, 8, 6, 10))

save_fig(fig2, "Figure2", w = 14, h = 15)

# ============================================================
# Figure 3
# ============================================================
ssgsea_path <- file.path(results_dir, "05_3_module_de/data/module_scores_ssgsea.rds")
clin_path   <- file.path(results_dir, "01_1b_batch_enhance/data/clinical_wgcna.rds")
if (!file.exists(ssgsea_path) || !file.exists(clin_path)) {
  stop("Figure 3B: ssGSEA or clinical file missing")
}
ssgsea_mat <- readRDS(ssgsea_path)
clin_df    <- readRDS(clin_path)
gy_ssgsea <- as.data.frame(t(ssgsea_mat)) %>%
  tibble::rownames_to_column("sample_id") %>%
  left_join(clin_df %>% select(sample_id, grade), by = "sample_id") %>%
  group_by(grade) %>%
  summarise(mean_ssgsea = mean(greenyellow, na.rm = TRUE), .groups = "drop")

path_limma <- read_required(file.path(results_dir,
  "04_2_pathway_score/tables/pathway_limma_multi_adjusted.csv"))
check_cols(path_limma, c("contrast", "pathway", "logFC", "t", "FDR"),
           "pathway_limma_multi_adjusted")
g3v2_sig <- path_limma %>%
  filter(contrast == "g3_vs_g2") %>%
  filter(grepl("Mito_|OXPHOS|Fatty_acid|TCA", pathway, ignore.case = TRUE)) %>%
  filter(FDR < 0.05) %>%
  arrange(logFC)
log_info(sprintf("Figure 3A: %d pathways with FDR < 0.05", nrow(g3v2_sig)))

p3a <- ggplot(g3v2_sig, aes(x = logFC, y = reorder(pathway, logFC))) +
  geom_point(size = 3.5, color = "darkred") +
  geom_errorbar(aes(xmin = logFC - 1.96*abs(logFC/t),
                    xmax = logFC + 1.96*abs(logFC/t)),
                orientation = "y", height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "logFC (G3 vs G2)", y = "",
       title = "A  Mitochondrial pathways at G3 vs G2 (FDR < 0.05)") +
  theme_pub +
  scale_y_discrete(expand = expansion(add = c(0.4, 0.4)),
                   labels = pathway_labels) +
  theme(axis.text.y = element_text(size = 9),
        plot.margin = margin(6, 8, 6, 40))

metab <- read_required(file.path(results_dir,
  "07_5d_metabolism_pathways/tables/metabolic_pathway_vs_grade.csv"))
check_cols(metab, c("pathway", "beta", "FDR"), "metabolic_pathway_vs_grade")

pathway_perSD <- read_required(file.path(results_dir,
  "07_5c_OXPHOS_analysis/tables/OXPHOS_pathway_perSD.csv")) %>%
  filter(model == "adjusted")

oxphos_metab <- metab %>% filter(pathway == "OXPHOS")
p3b_caption <- sprintf(
  paste0(
    "Per-grade model: OXPHOS beta = %.3f (nominal p = %.3f; FDR = %.3f).\n",
    "Per 1 SD: beta = %.3f (95%% CI %.3f to %.3f), p = %.3f.\n",
    "greenyellow ssGSEA (not plotted): G1 = %.3f, G2 = %.3f, G3 = %.3f, G4 = %.3f."),
  oxphos_metab$beta, oxphos_metab$p, oxphos_metab$FDR,
  pathway_perSD$beta_per_SD, pathway_perSD$ci_lower,
  pathway_perSD$ci_upper, pathway_perSD$p_value,
  gy_ssgsea$mean_ssgsea[gy_ssgsea$grade == 1],
  gy_ssgsea$mean_ssgsea[gy_ssgsea$grade == 2],
  gy_ssgsea$mean_ssgsea[gy_ssgsea$grade == 3],
  gy_ssgsea$mean_ssgsea[gy_ssgsea$grade == 4])
log_info(sprintf("Figure 3B caption: %s", p3b_caption))

p3b <- ggplot(metab, aes(x = beta, y = reorder(pathway, beta), fill = FDR < 0.05)) +
  geom_col() +
  scale_fill_manual(values = c("grey70", "firebrick"), name = "FDR<0.05") +
  scale_y_discrete(expand = expansion(add = c(0.4, 0.4)),
                   labels = pathway_labels) +
  labs(x = "Beta per grade", y = "",
       title = "B  Metabolic pathways: OXPHOS-selective",
       caption = p3b_caption) +
  theme_pub +
  theme(legend.position = "none",
        axis.text.y = element_text(size = 9),
        plot.caption = element_text(size = 7, color = "grey30", hjust = 0))

cmp <- read_required(file.path(results_dir,
  "07_7b_MR_plink/MR_cross_cohort_comparison.csv"))
check_cols(cmp, c("symbol", "OR_all", "OR_early"), "MR_cross_cohort_comparison")

mr_double <- cmp %>%
  select(symbol, OR_all, OR_early, p_all, p_early) %>%
  pivot_longer(cols = c(OR_all, OR_early),
               names_to = "outcome", values_to = "OR") %>%
  mutate(outcome = recode(outcome,
                          OR_all   = "All-stage AMD",
                          OR_early = "Early AMD"),
         gene_type = case_when(
           grepl("^NDUF|^COX|^UQCR|^SDH|^GABARAPL2", symbol) ~ "OXPHOS",
           grepl("^CFH|^CFHR|^CD59|^C3$|^CFB$", symbol) ~ "Complement",
           TRUE ~ "Other"),
         highlight = case_when(
           symbol == "NDUFB3" ~ "NDUFB3",
           gene_type == "Complement" ~ "Complement",
           TRUE ~ "Other OXPHOS"))

x_lo <- max(0.001, min(mr_double$OR, na.rm = TRUE) * 0.7)
x_hi <- max(mr_double$OR, na.rm = TRUE) * 1.3

p3c <- ggplot(mr_double,
              aes(x = OR, y = reorder(symbol, OR), color = highlight,
                  shape = outcome)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey30") +
  scale_x_log10(limits = c(x_lo, x_hi)) +
  scale_color_manual(values = c("NDUFB3" = "firebrick",
                                "Complement" = "#2166AC",
                                "Other OXPHOS" = "grey50"),
                     name = "") +
  scale_shape_manual(values = c("All-stage AMD" = 16, "Early AMD" = 17),
                     name = "") +
  labs(title = "C  MR: OXPHOS genes and complement positive controls",
       x = "Odds ratio (log scale)", y = "") +
  theme_pub +
  scale_y_discrete(expand = expansion(add = c(0.4, 0.4))) +
  theme(legend.position = "bottom",
        legend.box = "horizontal",
        legend.spacing.x = unit(0.15, "cm"),
        legend.key.width = unit(0.4, "cm"),
        legend.text = element_text(size = 7.5))

bn <- read_required(file.path(results_dir,
  "07_6a_bayesian_network/tables/bn_edges_into_grade.csv"))
sem <- read_required(file.path(results_dir,
  "07_6b_lavaan_mediation/tables/model_comparison.csv"))

bn_data_only <- bn %>%
  mutate(edge = paste0(from, " -> grade"),
         edge_group = ifelse(grepl("^age", edge), "Age-grade", "Molecular")) %>%
  select(edge, strength_data, edge_group)

x_max   <- max(bn_data_only$strength_data, na.rm = TRUE)
x_upper <- if (x_max > 0.90) 1.05 else 0.90

p3d_bn <- ggplot(bn_data_only,
                 aes(x = strength_data,
                     y = reorder(edge, strength_data),
                     fill = edge_group)) +
  geom_col() +
  scale_fill_manual(values = c("Age-grade" = "grey60",
                               "Molecular" = "steelblue"), name = "") +
  geom_vline(xintercept = 0.85, linetype = "dashed", color = "red") +
  labs(x = "Bootstrap edge strength", y = "",
       title = "D1  BN edges into grade (data-driven only)",
       caption = "Edges ordered by bootstrap strength; strong edges in Supplementary Table S8.") +
  theme_pub +
  scale_y_discrete(expand = expansion(add = c(0.5, 0.5))) +
  scale_x_continuous(limits = c(0, x_upper), breaks = seq(0, x_upper, 0.2)) +
  theme(axis.text.y = element_text(size = 8),
        legend.position = "right",
        legend.justification = c(0, 0.5),
        plot.caption = element_text(size = 7, color = "grey30", hjust = 0))

p3d_sem <- ggplot(sem, aes(x = estimate, y = model)) +
  geom_point(size = 2.5, color = "darkred") +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high),
                orientation = "y", height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Standardized estimate", y = "",
       title = "D2  SEM mediation") +
  theme_pub +
  theme(axis.text.y = element_text(size = 8),
        plot.margin = margin(6, 8, 6, 60))

p3d <- p3d_bn / p3d_sem

fig3 <- (p3a | p3b) / (p3c | p3d) +
  plot_layout(heights = c(1, 1.4)) &
  theme(plot.margin = margin(6, 8, 6, 10))

save_fig(fig3, "Figure3", w = 14.5, h = 12)

# ============================================================
# Figure 4
# v4.11: 4B = greenyellow permutation (moved from old 2C);
#        old 4B/4C/4D -> 4C/4D/4E
# ============================================================
mod_trait <- read_required(file.path(results_dir,
  "05_2_module_trait/tables/module_trait_correlation.csv"))
mod_sum <- read_required(file.path(results_dir,
  "05_2_module_trait/tables/module_summary_annotated.csv"))

mat_df <- mod_trait %>%
  select(Module, Trait, r) %>%
  pivot_wider(names_from = Trait, values_from = r)

col_order <- c("grade", "hUSI", "ML_Aging", "genes16", "mito",
               "age", "RIN", "PMI", "batch", "layer", "sex")
col_order <- intersect(col_order, colnames(mat_df))
mat_df <- mat_df %>% select(Module, all_of(col_order))

mat <- as.matrix(mat_df[, -1])
rownames(mat) <- mat_df$Module
mat[is.na(mat)] <- 0

row_order <- mod_sum %>%
  filter(Module %in% rownames(mat)) %>%
  arrange(FDR_grade) %>%
  pull(Module)
row_order <- c(row_order, setdiff(rownames(mat), row_order))
mat <- mat[row_order, , drop = FALSE]

p4a_grob <- pheatmap(
  mat, cluster_rows = FALSE, cluster_cols = FALSE,
  color = colorRampPalette(c("navy", "white", "firebrick"))(50),
  fontsize = 9, silent = TRUE
)$gtable

pdf(file.path(fig_dir, "Figure4A.pdf"), width = 8, height = 6)
grid.newpage(); grid.draw(p4a_grob); dev.off()
tiff(file.path(fig_dir, "Figure4A.tiff"), width = 8, height = 6,
     units = "in", res = 300, compression = "lzw")
grid.newpage(); grid.draw(p4a_grob); dev.off()
log_info("saved figure: Figure4A (PDF + TIFF)")

p4a <- wrap_elements(p4a_grob) +
  labs(title = "A  Module-trait correlations (full table: Supplementary S4)")

# ---- Figure 4B: greenyellow permutation ----
p4b_perm <- ggplot(perm_mod, aes(x = greenyellow)) +
  geom_histogram(bins = 30, fill = "grey70", color = "black") +
  geom_vline(xintercept = obs_gy, color = "red", linetype = "dashed") +
  coord_cartesian(xlim = c(0, max(perm_mod$greenyellow, obs_gy) + 1)) +
  labs(x = "p3-significant genes", y = "Frequency",
       title = "B  Permutation: greenyellow module") +
  theme_pub

eqtm_bg <- read_required(file.path(results_dir,
  "02_9_eQTM_background_sensitivity/tables/eqtm_background_sensitivity.csv"))
check_cols(eqtm_bg,
  c("target", "OR_20000", "OR_14233", "FDR"),
  "eqtm_background_sensitivity")

gy_eqtm <- eqtm_bg %>% filter(target == "mod_greenyellow")
gy_or_20000 <- round(gy_eqtm$OR_20000, 2)
gy_or_14233 <- round(gy_eqtm$OR_14233, 2)
gy_fdr      <- gy_eqtm$FDR
log_info(sprintf("Figure 4C: greenyellow OR_20000=%.2f OR_14233=%.2f FDR=%.3e",
                 gy_or_20000, gy_or_14233, gy_fdr))

p4b_caption <- sprintf(
  "20,000-gene background; 14,233-gene background: OR = %.2f", gy_or_14233)

# ---- Figure 4C: eQTM enrichment ----
p4c <- ggplot(eqtm_bg, aes(x = reorder(target, -OR_20000), y = OR_20000)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(x = "", y = "Odds ratio",
       title = sprintf("C  eQTM enrichment: greenyellow dominant\n(OR = %.2f, FDR = %.2e)",
                       gy_or_20000, gy_fdr),
       caption = p4b_caption) +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        plot.caption = element_text(size = 7, color = "grey30", hjust = 0))

cross <- read_required(file.path(results_dir,
  "07_5b_GSE231536_attachments/tables/cross_A_eqtm_egene.csv"))
cross_b <- read_required(file.path(results_dir,
  "07_5b_GSE231536_attachments/tables/cross_B_detail.csv"))

tier_rdh5 <- cross_b %>%
  filter(target == "p3_pathway", reference == "tier1") %>%
  pull(gene)
tier_cfi <- cross_b %>%
  filter(target == "p3_pathway", reference == "tier2") %>%
  pull(gene)
log_info(sprintf("Figure 4D: Tier1=%s, Tier2=%s",
                 paste(tier_rdh5, collapse = ","),
                 paste(tier_cfi, collapse = ",")))

# ---- Figure 4D: eQTM vs eGene ----
p4d <- ggplot(cross, aes(x = target, y = fold, fill = reference)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(x = "", y = "Fold enrichment",
       title = "D  eQTM vs eGene: two independent axes",
       caption = sprintf(
         "Within p3_pathway: %s (Tier1) and %s (Tier2) overlap AMD GWAS causal genes.",
         paste(tier_rdh5, collapse = ","),
         paste(tier_cfi, collapse = ","))) +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
        legend.position = "bottom",
        legend.box = "horizontal",
        plot.caption = element_text(size = 7, color = "grey30", hjust = 0))

atac_perm <- read_required(file.path(results_dir,
  "07_1_atac_validation/tables/permutation_results.csv"))
check_cols(atac_perm,
  c("gene_set", "tissue", "level", "FDR_coef_neg", "obs_mean_coef",
    "FDR_dar"),
  "permutation_results")

atac_gene <- atac_perm %>%
  filter(level == "gene",
         gene_set %in% c("16genes", "greenyellow")) %>%
  mutate(sig = ifelse(FDR_coef_neg < 0.05, "FDR<0.05", "n.s."),
         tissue = factor(tissue, levels = c("RPE", "Retina")))

gy_rpe_fdr <- atac_gene %>%
  filter(gene_set == "greenyellow", tissue == "RPE") %>%
  pull(FDR_coef_neg)
t16_rpe_fdr <- atac_gene %>%
  filter(gene_set == "16genes", tissue == "RPE") %>%
  pull(FDR_coef_neg)
t16_ret_fdr <- atac_gene %>%
  filter(gene_set == "16genes", tissue == "Retina") %>%
  pull(FDR_coef_neg)
gy_ret_fdr <- atac_gene %>%
  filter(gene_set == "greenyellow", tissue == "Retina") %>%
  pull(FDR_coef_neg)
log_info(sprintf("Figure 4E: greenyellow RPE FDR=%.4g; 16genes RPE FDR=%.4g; 16genes Retina FDR=%.4g; greenyellow Retina FDR=%.4g",
                 gy_rpe_fdr, t16_rpe_fdr, t16_ret_fdr, gy_ret_fdr))

p4d_caption <- sprintf(
  "greenyellow RPE FDR = %.4g; 16genes RPE FDR = %.4g; 16genes Retina FDR = %.4g; greenyellow Retina FDR = %.4g.",
  gy_rpe_fdr, t16_rpe_fdr, t16_ret_fdr, gy_ret_fdr)

# ---- Figure 4E: ATAC-seq ----
p4e <- ggplot(atac_gene, aes(x = gene_set, y = obs_mean_coef, fill = sig)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ tissue, ncol = 2) +
  scale_fill_manual(
    values = c("FDR<0.05" = "firebrick",
               "n.s."     = "grey70"),
    breaks = c("FDR<0.05", "n.s."),
    labels = c("FDR<0.05", "n.s."),
    name = "") +
  scale_y_continuous(limits = c(-0.75, 0.05),
                     breaks = seq(-0.6, 0, 0.2)) +
  labs(x = "", y = "Mean coefficient",
       title = "E  ATAC-seq promoter accessibility (RPE + Retina)",
       caption = p4d_caption) +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "bottom",
        legend.box = "horizontal",
        plot.caption = element_text(size = 7, color = "grey30", hjust = 0))

# ---- Assemble Figure 4: 4A / (4B | 4C) / (4D | 4E) ----
fig4 <- p4a /
        (p4b_perm | p4c) /
        (p4d | p4e) +
  plot_layout(heights = c(1.4, 1, 1)) &
  theme(plot.margin = margin(6, 8, 6, 10))

save_fig(fig4, "Figure4", w = 14, h = 16)

# ============================================================
# Figure 5
# ============================================================
mr_early <- read_required(file.path(results_dir,
  "07_7b_MR_plink/MR_earlyAMD.csv"))
grid <- read_required(file.path(results_dir,
  "07_7b_MR_plink/NDUFB3_sensitivity_grid.csv"))
adipose <- read_required(file.path(results_dir,
  "07_7c_GTEx_MR_NDUFB3/tables/MR_results.csv"))
eqtl_multi <- read_required(file.path(results_dir,
  "07_7a_NDUFB3/NDUFB3_multi_tissue_eQTL.csv"))
coloc <- read_required(file.path(results_dir,
  "07_7f_coloc_NDUFB3/coloc_NDUFB3_both_tissues.csv"))

blood_single <- mr_early %>%
  filter(symbol == "NDUFB3", method == "Wald ratio") %>%
  transmute(label = "Blood (eQTLGen, single-SNP Wald)",
            OR = OR_MR, lower = OR_lower, upper = OR_upper, p = p_MR)

blood_ivw <- grid %>%
  filter(outcome == "earlyAMD", n_iv >= 3, method == "IVW") %>%
  slice(1) %>%
  transmute(label = "Blood (eQTLGen, 3-SNP IVW)",
            OR = OR_MR, lower = OR_lower, upper = OR_upper, p = p_MR)
log_info(sprintf("Figure 5A: blood 3-SNP IVW OR=%.3f CI=%.3f-%.3f p=%.2e",
                 blood_ivw$OR, blood_ivw$lower, blood_ivw$upper, blood_ivw$p))

adipose_row <- adipose %>%
  transmute(label = "Adipose (GTEx, single-SNP Wald)",
            OR = or, lower = or_lci95, upper = or_uci95, p = pval)

ndufb3_replication <- bind_rows(blood_single, blood_ivw, adipose_row)

p5a <- ggplot(ndufb3_replication, aes(x = OR, y = label)) +
  geom_point(size = 4, color = "firebrick") +
  geom_errorbar(aes(xmin = lower, xmax = upper),
                orientation = "y", height = 0.2, color = "firebrick") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey30") +
  scale_x_log10() +
  labs(x = "Odds ratio (log scale)", y = "",
       title = "A  NDUFB3 replication in blood and adipose") +
  theme_pub +
  scale_y_discrete(expand = expansion(add = c(0.4, 0.4)))

grid_plot <- grid %>%
  filter(!is.na(OR_MR), outcome == "earlyAMD") %>%
  mutate(logOR = log(OR_MR),
         p_label = ifelse(p_MR < 0.001, "***",
                   ifelse(p_MR < 0.01, "**",
                   ifelse(p_MR < 0.05, "*", "ns"))))

p5b <- ggplot(grid_plot, aes(x = factor(r2_thr),
                             y = factor(p_thr, levels = c("1e-06", "1e-05", "1e-04")),
                             fill = logOR)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f\n%s", OR_MR, p_label)), size = 2.6) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, name = "log(OR)") +
  labs(x = "LD r2 threshold", y = "eQTL P-value threshold",
       title = "B  NDUFB3 sensitivity grid (early AMD)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank())

eqtl_plot <- eqtl_multi %>%
  mutate(tissue_label = gsub("_", " ", tissue),
         neg_log10p = -log10(pval),
         gws = ifelse(pval < 5e-8, "Genome-wide significant", "Nominal")) %>%
  arrange(nes)

p5c <- ggplot(eqtl_plot, aes(x = nes, y = reorder(tissue_label, nes))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(aes(size = neg_log10p, fill = gws),
             shape = 21, color = "black", stroke = 0.3) +
  scale_fill_manual(values = c("Genome-wide significant" = "firebrick",
                               "Nominal" = "grey70"), name = "") +
  scale_size_continuous(range = c(2, 6), name = expression(-log[10](P))) +
  labs(x = "Normalised effect size (NES)", y = "",
       title = "C  NDUFB3 multi-tissue eQTL (GTEx V10)") +
  theme_pub +
  scale_y_discrete(expand = expansion(add = c(0.4, 0.4))) +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 8))

coloc_plot <- coloc %>%
  mutate(tissue_label = recode(tissue,
                               "Adipose_Subcutaneous" = "Subcutaneous adipose\n(n = 581)",
                               "Whole_Blood" = "Whole blood\n(n = 670)"),
         r2_label = sprintf("r² = %.2f", r2_top_variants)) %>%
  select(tissue_label, PP_H4, PP_H3, PP_H1, r2_label) %>%
  pivot_longer(cols = c(PP_H4, PP_H3, PP_H1),
               names_to = "hypothesis", values_to = "PP") %>%
  mutate(hypothesis = recode(hypothesis,
                             PP_H4 = "H4 (shared causal)",
                             PP_H3 = "H3 (independent)",
                             PP_H1 = "H1 (eQTL only)"))

p5d <- ggplot(coloc_plot, aes(x = tissue_label, y = PP, fill = hypothesis)) +
  geom_col(position = "dodge") +
  geom_text(data = coloc_plot %>% filter(hypothesis == "H4 (shared causal)"),
            aes(label = r2_label), vjust = -0.5, size = 3, fontface = "bold") +
  scale_fill_manual(values = c("H4 (shared causal)" = "firebrick",
                               "H3 (independent)" = "grey50",
                               "H1 (eQTL only)" = "steelblue"),
                    name = "") +
  labs(x = "", y = "Posterior probability",
       title = "D  Colocalization: NDUFB3 cis-eQTL vs early AMD GWAS",
       subtitle = "GTEx v6p eQTL via COLSTATS v2.3") +
  theme_pub +
  theme(legend.position = "bottom",
        axis.text.x = element_text(size = 9))

fig5 <- (p5a | p5b) / (p5c | p5d) +
  plot_layout(heights = c(1, 1.2)) &
  theme(plot.margin = margin(6, 8, 6, 10))
save_fig(fig5, "Figure5", w = 13, h = 11)

# ============================================================
# Figure 6
# ============================================================
rgc <- cross_complex %>% filter(!is.na(rgc_log2FC))

p6a <- ggplot(rgc, aes(x = complex_std, y = rgc_log2FC, fill = rgc_log2FC > 0)) +
  geom_col() +
  geom_errorbar(aes(ymin = rgc_IQR_lo, ymax = rgc_IQR_hi), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(aes(
    y = ifelse(rgc_log2FC > 0, rgc_IQR_hi, rgc_IQR_lo),
    label = sprintf("%.3f", rgc_log2FC),
    vjust = ifelse(rgc_log2FC > 0, -0.5, 1.5)
  ), size = 3, fontface = "bold") +
  scale_fill_manual(values = c("steelblue", "firebrick"), name = "Up") +
  labs(x = "OXPHOS complex", y = "log2FC (dAMD vs Healthy)",
       title = "A  RGC: complex redistribution") +
  theme_pub

per_cplx <- rgc %>%
  select(complex_std, pct_log2FC_pos, pct_p_lt_005) %>%
  mutate(pct_p_lt_005 = pct_p_lt_005 * 100,
         complex_std = as.character(complex_std))

agg_oxphos <- read_required(file.path(results_dir,
  "07_3_NDUFB8_coexpression/tables/RGC_subsampling_100iter_summary.csv")) %>%
  filter(metric == "OXPHOS") %>%
  transmute(complex_std = "Aggregate OXPHOS",
            pct_log2FC_pos = pct_log2FC_pos,
            pct_p_lt_005   = pct_p_lt_005)

iter_summary <- bind_rows(per_cplx, agg_oxphos) %>%
  pivot_longer(cols = c(pct_log2FC_pos, pct_p_lt_005),
               names_to = "metric", values_to = "pct") %>%
  mutate(complex_std = factor(complex_std,
    levels = c("Complex I", "Complex II", "Complex III", "Complex IV",
               "Complex V", "Aggregate OXPHOS")))

p6b <- ggplot(iter_summary, aes(x = complex_std, y = pct, fill = metric)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = 50, linetype = "dashed", color = "grey30") +
  scale_fill_manual(values = c("pct_log2FC_pos" = "firebrick",
                               "pct_p_lt_005"   = "steelblue"),
                    labels = c("% iterations positive log2FC",
                               "% iterations p < 0.05"),
                    name = "") +
  labs(x = "OXPHOS complex", y = "Proportion of iterations (%)",
       title = "B  Subsampling consistency (100 iterations)",
       caption = paste0(
         "Aggregate OXPHOS: 84% of iterations p < 0.05 (median p = 0.042), ",
         "up-regulated (median log2FC = +0.0135, 100% positive).\n",
         "No individual complex reached p < 0.05 in any iteration.")) +
  theme_pub +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "bottom",
        plot.caption = element_text(size = 7, color = "grey30", hjust = 0))

cross_layer <- cross_complex %>%
  filter(!is.na(bulk_beta), !is.na(rgc_log2FC)) %>%
  mutate(
    bulk_z = as.numeric(scale(bulk_beta)),
    rgc_z  = as.numeric(scale(rgc_log2FC))
  ) %>%
  select(complex_std, bulk_z, rgc_z) %>%
  pivot_longer(-complex_std, names_to = "layer", values_to = "effect") %>%
  mutate(layer = recode(layer,
                        bulk_z = "Bulk retina (z)",
                        rgc_z  = "RGC (z)"))

p6c_caption <- {
  raw <- cross_complex %>%
    filter(complex_std %in% c("Complex I", "Complex III", "Complex IV")) %>%
    select(complex_std, bulk_beta, rgc_log2FC)
  raw <- raw[match(c("Complex I", "Complex III", "Complex IV"), raw$complex_std), ]
  sprintf(
    paste0(
      "Values are z-scored within each layer across the three complexes with ",
      "data in both layers (I, III, IV). Because z-scoring is relative within ",
      "each layer, a same-direction effect with a smaller magnitude can appear ",
      "reversed; direction should be read from the raw values.\n",
      "Complex I (bulk beta = %.4f; RGC log2FC = %.4f) and Complex IV ",
      "(bulk beta = %.4f; RGC log2FC = %.4f) are concordantly downregulated ",
      "in both layers; Complex III is discordant ",
      "(bulk beta = %.4f; RGC log2FC = %.4f)."),
    raw$bulk_beta[1], raw$rgc_log2FC[1],
    raw$bulk_beta[3], raw$rgc_log2FC[3],
    raw$bulk_beta[2], raw$rgc_log2FC[2]
  )
}

p6c <- ggplot(cross_layer, aes(x = complex_std, y = effect, fill = layer)) +
  geom_col(position = "dodge") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("Bulk retina (z)" = "#2166AC",
                               "RGC (z)" = "#B2182B"),
                    name = "") +
  labs(x = "OXPHOS complex", y = "Standardized effect (z-score)",
       title = "C  Cross-layer direction reversal for Complex III",
       caption = p6c_caption) +
  theme_pub +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "bottom",
        plot.caption = element_text(size = 7, color = "grey30", hjust = 0))

cell_prop <- read_required(file.path(results_dir,
  "06_1_singlecell/tables/celltype_by_condition_ttest.csv"))
cell_prop_long <- cell_prop %>%
  select(celltype, mean_healthy, mean_damd) %>%
  pivot_longer(-celltype, names_to = "condition", values_to = "mean") %>%
  mutate(condition = recode(condition,
                            mean_healthy = "Healthy",
                            mean_damd = "dAMD"))

p6d <- ggplot(cell_prop_long, aes(x = celltype, y = mean, fill = condition)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c("Healthy" = "grey60", "dAMD" = "firebrick"),
                    name = "") +
  labs(x = "", y = "Mean proportion (%)",
       title = "D  Cell-type proportions") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))

fig6 <- (p6a | p6b) / (p6c | p6d) +
  plot_layout(heights = c(1, 1)) &
  theme(plot.margin = margin(6, 8, 6, 10))
save_fig(fig6, "Figure6", w = 13, h = 10)

# ============================================================
# Figure 7
# ============================================================
gsea_rgc <- read_required(file.path(results_dir,
  "07_3_scTenifoldKnk_NDUFB3_RGC/tables/NDUFB3_GSEA_RGC.csv"))
gsea_rpe <- read_required(file.path(results_dir,
  "07_3_scTenifoldKnk_NDUFB3_RPE/tables/NDUFB3_GSEA_RPE.csv"))
gsea_rod <- read_required(file.path(results_dir,
  "07_3_scTenifoldKnk_NDUFB3_ROD/tables/NDUFB3_GSEA_Rod.csv"))

gsea_all <- bind_rows(
  gsea_rgc %>% mutate(celltype = "RGC"),
  gsea_rpe %>% mutate(celltype = "RPE"),
  gsea_rod %>% mutate(celltype = "ROD")
)

oxphos_pathways <- c("Mito_Complex_I", "Mito_Complex_II", "Mito_Complex_III",
                     "Mito_Complex_IV", "Mito_Complex_V", "OXPHOS")

gsea_oxphos <- gsea_all %>%
  filter(pathway %in% oxphos_pathways) %>%
  mutate(pathway = factor(pathway, levels = oxphos_pathways),
         celltype = factor(celltype, levels = c("RGC", "RPE", "ROD")),
         sig = factor(ifelse(padj < 0.05, "FDR < 0.05", "n.s."),
                      levels = c("FDR < 0.05", "n.s.")))

p7a <- ggplot(gsea_oxphos, aes(x = pathway, y = NES, fill = sig)) +
  geom_col() +
  facet_wrap(~ celltype, ncol = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("FDR < 0.05" = "firebrick", "n.s." = "grey70"),
                    breaks = c("FDR < 0.05", "n.s."),
                    labels = c("FDR < 0.05", "n.s."),
                    name = "",
                    drop = FALSE) +
  labs(x = "", y = "NES",
       title = "A  NDUFB3 virtual KO: OXPHOS complex gene sets (GSEA)") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
        legend.position = "bottom")

gene16_rgc <- read_required(file.path(results_dir,
  "07_3_scTenifoldKnk_NDUFB3_RGC/tables/NDUFB3_vs_16genes_RGC.csv"))
gene16_rpe <- read_required(file.path(results_dir,
  "07_3_scTenifoldKnk_NDUFB3_RPE/tables/NDUFB3_vs_16genes_RPE.csv"))
gene16_rod <- read_required(file.path(results_dir,
  "07_3_scTenifoldKnk_NDUFB3_ROD/tables/NDUFB3_vs_16genes_Rod.csv"))

gene16_all <- bind_rows(
  gene16_rgc %>% mutate(celltype = "RGC"),
  gene16_rpe %>% mutate(celltype = "RPE"),
  gene16_rod %>% mutate(celltype = "ROD")
) %>%
  mutate(celltype = factor(celltype, levels = c("RGC", "RPE", "ROD")),
         sig = factor(ifelse(p.adj < 0.05, "FDR < 0.05", "n.s."),
                      levels = c("FDR < 0.05", "n.s.")))

dummy_row <- gene16_all[1, ]
dummy_row$gene <- "DUMMY"
dummy_row$sig <- factor("FDR < 0.05", levels = c("FDR < 0.05", "n.s."))
dummy_row$Z <- NA

gene16_all_fixed <- bind_rows(gene16_all, dummy_row)

p7b <- ggplot(gene16_all_fixed, aes(x = Z, y = reorder(gene, Z), fill = sig)) +
  geom_col(na.rm = TRUE) +
  facet_wrap(~ celltype, ncol = 3, scales = "free_y") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(
    values = c("FDR < 0.05" = "firebrick",
               "n.s."     = "grey70"),
    breaks = c("FDR < 0.05", "n.s."),
    labels = c("FDR < 0.05", "n.s."),
    name = "",
    drop = FALSE
  ) +
  labs(x = "Z-score", y = "",
       title = "B  NDUFB3 network perturbation: 16-gene candidate set per cell type",
       caption = "Genes shown are the evaluable members of the 16-gene candidate set in each cell type; members not evaluable in a given cell type are omitted.") +
  theme_pub +
  theme(legend.position = "bottom",
        axis.text.y = element_text(size = 7),
        plot.caption = element_text(size = 7, color = "grey30", hjust = 0))

gsea_tile <- gsea_all %>%
  filter(pathway %in% oxphos_pathways) %>%
  select(pathway, celltype, NES, padj) %>%
  mutate(pathway = factor(pathway, levels = rev(oxphos_pathways)),
         celltype = factor(celltype, levels = c("RGC", "RPE", "ROD")),
         sig_label = ifelse(padj < 0.05, "*", ""),
         tile_label = sprintf("%.2f%s", NES, sig_label))

p7c <- ggplot(gsea_tile, aes(x = celltype, y = pathway, fill = NES)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = tile_label), size = 3, color = "black") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "firebrick",
                       midpoint = 0, name = "NES") +
  labs(x = "", y = "",
       title = "C  Cross-cell-type NES (OXPHOS complex gene sets)",
       caption = paste0(
         "* FDR < 0.05; missing tile = gene set not present in that cell type; ",
         "RPE column shows no significant enrichment (all FDR > 0.9)")) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        plot.caption = element_text(size = 7, color = "grey30", hjust = 0))

fig7 <- p7a / p7b / p7c +
  plot_layout(heights = c(1.4, 1, 0.9)) &
  theme(plot.margin = margin(6, 8, 6, 10))
save_fig(fig7, "Figure7", w = 12, h = 14)

log_info("=== 10_1A v4.11 done (7 figures) ===")


