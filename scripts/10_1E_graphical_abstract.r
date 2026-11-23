# 10_1E_graphical_abstract.R (v11.6)
# Graphical abstract

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressWarnings(suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr); library(tidyr); library(readr)
  library(patchwork); library(grid); library(tibble)
}))

results_dir <- DIR_RESULT
fig_dir     <- file.path(results_dir, "10_1E_graphical_abstract")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
log_info("=== 10_1E Graphical Abstract (v11.6) ===")

col_bulk <- "#2166AC"; col_rgc <- "#B2182B"; col_gy <- "#7FBC41"; col_grey <- "grey70"

theme_ga <- theme_void(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
        plot.margin = margin(4, 4, 4, 4))

normalize_complex <- function(x) {
  x <- toupper(trimws(x))
  x <- gsub("^COMPLEX\\s*", "", x)
  x <- gsub("^C\\s*", "", x)
  x <- trimws(x)
  num_to_roman <- c("1" = "I", "2" = "II", "3" = "III", "4" = "IV", "5" = "V")
  x <- ifelse(x %in% names(num_to_roman), num_to_roman[x], x)
  paste0("Complex ", x)
}

read_required <- function(path) {
  if (!file.exists(path)) stop(sprintf("REQUIRED FILE MISSING: %s", path))
  read_csv(path, show_col_types = FALSE, progress = FALSE)
}

# ------------------------------------------------------------
# Read data
# ------------------------------------------------------------
log_info("--- reading shared data ---")
bulk_cplx <- read_required(file.path(results_dir,
  "07_5c_OXPHOS_analysis/tables/OXPHOS13_complex_score.csv")) %>%
  mutate(complex_std = normalize_complex(complex))

rgc_cplx <- read_required(file.path(results_dir,
  "07_3_NDUFB8_coexpression/tables/RGC_complex_100iter_summary.csv")) %>%
  mutate(complex_std = normalize_complex(complex))

cross <- read_required(file.path(results_dir,
  "07_5b_GSE231536_attachments/tables/cross_A_eqtm_egene.csv"))
bn <- read_required(file.path(results_dir,
  "07_6a_bayesian_network/tables/bn_edges_into_grade.csv"))
cmp <- read_required(file.path(results_dir,
  "07_7b_MR_plink/MR_cross_cohort_comparison.csv"))

# --- Panel A: disease-stage framework ---
panel_a <- ggplot() +
  annotate("segment", x = 1, xend = 4, y = 0.5, yend = 0.5,
           arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
           linewidth = 1.2, color = "grey30") +
  annotate("point", x = 1:4, y = 0.5, size = 6, color = "grey30") +
  annotate("text", x = 1:4, y = 0.5, label = c("G1", "G2", "G3", "G4"),
           color = "white", fontface = "bold", size = 3.5) +
  annotate("rect", xmin = 1.5, xmax = 3.5, ymin = 0.15, ymax = 0.85,
           fill = "red", alpha = 0.15, color = "red", linetype = "dashed") +
  annotate("text", x = 2.5, y = 1.20,
           label = "G2-G3 transition:\nhighest disease-specific\nperturbation (pure = 69.8%;\nG2-G1 = 64.2%)",
           color = "red", size = 3.2, fontface = "bold", lineheight = 0.9) +
  annotate("text", x = 1, y = -0.15, label = "Normal", color = "grey30", size = 3) +
  annotate("text", x = 4, y = -0.15, label = "Advanced AMD", color = "grey30", size = 3) +
  coord_cartesian(xlim = c(0.5, 4.5), ylim = c(-0.4, 1.5)) +
  labs(title = "A  Disease-stage framework") +
  theme_ga

# --- Panel B: bulk OXPHOS downregulation ---
bulk_plot <- bulk_cplx %>%
  filter(!is.na(beta)) %>%
  mutate(complex_std = factor(complex_std,
                              levels = c("Complex I", "Complex III", "Complex IV")),
         sig = factor(ifelse(FDR < 0.05, "FDR<0.05", "n.s."),
                      levels = c("FDR<0.05", "n.s.")))

panel_b <- ggplot(bulk_plot, aes(x = complex_std, y = beta, fill = sig)) +
  geom_col(width = 0.6, show.legend = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  geom_text(aes(label = sprintf("%.3f", beta)),
            vjust = 1.5, size = 3, color = "white", fontface = "bold") +
  scale_fill_manual(
    values = c("FDR<0.05" = "firebrick", "n.s." = "grey70"),
    drop = FALSE,
    name = "") +
  guides(fill = guide_legend(
    override.aes = list(fill = c("firebrick", "grey70"))
  )) +
  labs(title = "B  Bulk retina: OXPHOS down",
       y = "Beta per grade", x = "") +
  theme_ga +
  theme(axis.text.x = element_text(size = 9, color = "black"),
        axis.title.y = element_text(size = 9),
        legend.position = "bottom",
        legend.text = element_text(size = 8))

# --- Panel C: causal re-positioning ---
bn_key <- bn %>%
  filter(grepl("OXPHOS|ComplexIII|UQCRQ|age|MEgreenyellow", from, ignore.case = TRUE)) %>%
  select(from, strength_data) %>%
  mutate(from_label = recode(from,
                             "age" = "Age",
                             "OXPHOS_pathway" = "OXPHOS pathway",
                             "OXPHOS13" = "OXPHOS13 score",
                             "ComplexIII_mean" = "Complex III",
                             "UQCRQ_z" = "UQCRQ",
                             "MEgreenyellow" = "greenyellow module",
                             .default = from))

panel_c <- ggplot(bn_key, aes(x = strength_data, y = reorder(from_label, strength_data))) +
  geom_col(fill = col_grey, width = 0.6) +
  geom_vline(xintercept = 0.85, linetype = "dashed", color = "red") +
  geom_text(aes(label = sprintf("%.3f", strength_data)),
            hjust = -0.15, size = 3) +
  labs(title = "C  Causal re-positioning: weak edges",
       x = "Bootstrap edge strength", y = "") +
  coord_cartesian(xlim = c(0, 1.15)) +
  theme_ga +
  theme(axis.text.x = element_text(size = 8, color = "black"),
        axis.text.y = element_text(size = 9, color = "black"),
        axis.title.x = element_text(size = 9))

# --- Panel D: MR forest ---
mr_plot <- cmp %>%
  filter(!is.na(OR_early)) %>%
  mutate(highlight = case_when(
    symbol == "NDUFB3" ~ "NDUFB3",
    grepl("^CFH|^CD59|^C3$", symbol) ~ "Complement",
    TRUE ~ "Other OXPHOS"
  ))

x_lo <- min(mr_plot$OR_early, na.rm = TRUE) * 0.7
x_hi <- max(mr_plot$OR_early, na.rm = TRUE) * 1.3

panel_d <- ggplot(mr_plot,
                  aes(x = OR_early, y = reorder(symbol, OR_early),
                      color = highlight)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey30") +
  scale_x_log10(limits = c(x_lo, x_hi)) +
  scale_color_manual(values = c("NDUFB3" = "firebrick",
                                "Complement" = "#2166AC",
                                "Other OXPHOS" = "grey50"),
                     name = "") +
  labs(title = "D  MR: OXPHOS genes and complement\ncontrols (early AMD)",
       x = "Odds ratio (log scale)", y = "") +
  theme_ga +
  theme(axis.text.x = element_text(size = 8, color = "black"),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.title.x = element_text(size = 9),
        plot.title = element_text(size = 10, lineheight = 0.9),
        legend.position = "bottom",
        legend.text = element_text(size = 8))

# --- Panel E: two-axis regulation ---
cross_wide <- cross %>%
  select(target, reference, fold) %>%
  pivot_wider(names_from = reference, values_from = fold)

eqtm_col  <- intersect(c("eQTM_targets", "eQTM", "eqtm"), colnames(cross_wide))[1]
egene_col <- intersect(c("eGenes", "eGene", "egene"), colnames(cross_wide))[1]
stopifnot(!is.na(eqtm_col), !is.na(egene_col))

cross_long <- cross_wide %>%
  filter(target %in% c("mod_greenyellow", "sig16")) %>%
  mutate(target_label = recode(target,
                               "mod_greenyellow" = "greenyellow\n(OXPHOS module)",
                               "sig16" = "16-gene\nsignature",
                               .default = target)) %>%
  pivot_longer(cols = c(all_of(eqtm_col), all_of(egene_col)),
               names_to = "axis_raw", values_to = "fold") %>%
  mutate(axis = recode(axis_raw,
                       "eQTM_targets" = "eQTM (epigenetic)",
                       "eQTM" = "eQTM (epigenetic)",
                       "eGenes" = "eGene (genetic)",
                       "eGene" = "eGene (genetic)"))

panel_e <- ggplot(cross_long, aes(x = target_label, y = fold, fill = axis)) +
  geom_col(position = "dodge", width = 0.6) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey30") +
  geom_text(aes(label = sprintf("%.1f", fold)),
            position = position_dodge(0.6), vjust = -0.5,
            size = 3.2, fontface = "bold") +
  scale_fill_manual(values = c("eQTM (epigenetic)" = col_gy,
                               "eGene (genetic)" = "grey40"), name = "") +
  labs(title = "E  Two-axis hypothesis (speculative)",
       y = "Fold enrichment", x = "") +
  theme_ga +
  theme(axis.text.x = element_text(size = 9, color = "black"),
        axis.title.y = element_text(size = 9),
        legend.position = "bottom",
        legend.text = element_text(size = 8))

# --- Panel F: single-cell complex redistribution ---
rgc_plot <- rgc_cplx %>%
  filter(!is.na(log2FC_median)) %>%
  mutate(complex_std = factor(complex_std,
                              levels = c("Complex I", "Complex II", "Complex III",
                                         "Complex IV", "Complex V")),
         direction = ifelse(log2FC_median > 0, "Up in dAMD", "Down in dAMD"),
         label_vjust = ifelse(log2FC_median > 0, -0.5, 1.5))

panel_f <- ggplot(rgc_plot, aes(x = complex_std, y = log2FC_median,
                                fill = direction)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey30") +
  geom_text(aes(label = sprintf("%.3f", log2FC_median),
                vjust = label_vjust),
            size = 3, color = "black", fontface = "bold") +
  scale_fill_manual(values = c("Up in dAMD" = col_rgc, "Down in dAMD" = col_bulk),
                    name = "") +
  labs(title = "F  Single-cell RGC: complex redistribution\n(exploratory; donor-level p > 0.05)",
       y = "log2FC (dAMD vs Healthy)", x = "") +
  theme_ga +
  theme(axis.text.x = element_text(size = 8, color = "black", angle = 30, hjust = 1),
        axis.title.y = element_text(size = 9),
        plot.title = element_text(size = 10, lineheight = 0.9),
        legend.position = c(0.02, 0.98),
        legend.justification = c(0, 1),
        legend.direction = "vertical",
        legend.key.size = unit(0.35, "cm"),
        legend.text = element_text(size = 7),
        legend.background = element_rect(fill = "white", color = NA),
        legend.margin = margin(1, 1, 1, 1))

# --- Panel G: key-conclusion text box ---
panel_g <- ggplot() +
  annotate("rect", xmin = 0, xmax = 1, ymin = 0, ymax = 1,
           fill = "#F7F7F7", color = "grey30", linewidth = 0.5) +
  annotate("text", x = 0.5, y = 0.97,
           label = "Key conclusions",
           fontface = "bold", size = 4.6, hjust = 0.5) +

  annotate("text", x = 0.05, y = 0.88,
           label = "1  OXPHOS downregulation is a reproducible\n     correlate, not a strong causal driver",
           size = 3.1, hjust = 0, lineheight = 1.02, fontface = "bold") +
  annotate("text", x = 0.07, y = 0.74,
           label = "Observational \u03b2 = -0.025 (composite) / -0.080\n(full pathway) per 1 SD, both p > 0.05;\nBN/SEM edges \u2264 0.014; G2\u2013G3 is the highest,\nbut not unique, perturbation (pure 69.8%)",
           size = 2.8, hjust = 0, lineheight = 0.98) +

  annotate("text", x = 0.05, y = 0.54,
           label = "2  Two-axis regulation of the OXPHOS axis",
           size = 3.1, hjust = 0, lineheight = 1.02, fontface = "bold") +
  annotate("text", x = 0.07, y = 0.40,
           label = "Genetic (16-gene, eGene-enriched) and epigenetic\n(greenyellow, eQTM OR = 8.91) axes are partly\nindependent; Complex III is nominally discordant\nacross layers (bulk down, RGC up; exploratory)",
           size = 2.8, hjust = 0, lineheight = 0.98) +

  annotate("text", x = 0.05, y = 0.22,
           label = "3  NDUFB3 is the only gene-level candidate",
           size = 3.1, hjust = 0, lineheight = 1.02, fontface = "bold") +
  annotate("text", x = 0.07, y = 0.10,
           label = "Stage-specific for early AMD, but instrument-fragile\n(blood IVW OR = 4.61; adipose OR = 1.87); requires\nexperimental validation",
           size = 2.8, hjust = 0, lineheight = 0.98) +

  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_void()

# --- Assemble final graphical abstract ---
log_info("--- assembling final graphical abstract ---")
top_row    <- panel_a + plot_layout(widths = 1)
mid_row    <- (panel_b | panel_c | panel_d) + plot_layout(widths = c(1, 1.1, 1))
bottom_row <- (panel_e | panel_f | panel_g) + plot_layout(widths = c(1, 1, 0.85))

final_ga <- top_row / mid_row / bottom_row +
  plot_layout(heights = c(0.75, 1, 1.60)) +
  plot_annotation(
    title = "Multi-scale causal inference in AMD progression",
    subtitle = "OXPHOS downregulation is a correlate, not a strong causal driver; the axis is split across genetic and epigenetic layers, and NDUFB3 is a provisional stage-specific candidate",
    theme = theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey30")
    )
  )

# --- Save figures ---
log_info("--- saving figures ---")

suppressMessages({
  ggsave(file.path(fig_dir, "GraphicalAbstract.pdf"),
         final_ga, width = 14, height = 10, device = cairo_pdf)

  tryCatch({
    ggsave(file.path(fig_dir, "GraphicalAbstract.tiff"),
           final_ga, width = 14, height = 10, dpi = 300, device = "tiff")
  }, error = function(e) {
    log_info("TIFF via device='tiff' failed; retrying via ragg::agg_tiff")
    ggsave(file.path(fig_dir, "GraphicalAbstract.tiff"),
           final_ga, width = 14, height = 10, dpi = 300, device = ragg::agg_tiff)
  })

  ggsave(file.path(fig_dir, "GraphicalAbstract.png"),
         final_ga, width = 14, height = 10, dpi = 200, device = "png")
})
log_info("saved graphical abstract (PDF + TIFF + PNG)")

log_info("=== 10_1E v11.6 done ===")