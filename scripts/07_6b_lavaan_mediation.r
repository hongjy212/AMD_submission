# ============================================================
# 07_6b_lavaan_mediation.r
# Purpose: SEM mediation
#   - OXPHOS: Model A (OXPHOS13), B (latent), C, D
#   - Complex III: Model A, B (latent), C, D for ComplexIII_mean
#   - UQCRQ single gene: Model A, C, D
# Input:  results/07_6_prepare_sem_bnlearn_data/data/df_cont.rds
# Output: results/07_6b_lavaan_mediation/
# ============================================================
source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(lavaan)
})

ROOT <- "E:/retina_tools/aging_rgc_axon/AMD"
OUT  <- file.path(ROOT, "results", "07_6b_lavaan_mediation")
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

safe_estimates <- function(fit, labels) {
  if (is.null(fit)) {
    return(data.frame(label = labels, est = NA_real_, se = NA_real_,
                       ci.lower = NA_real_, ci.upper = NA_real_,
                       pvalue = NA_real_))
  }
  pe <- parameterEstimates(fit, boot.ci.type = "bca.simple")
  pe[pe$label %in% labels,
     c("label","est","se","ci.lower","ci.upper","pvalue")]
}

log_info("=== 07_6b lavaan mediation ===")

# ============================================================
# 0. Load
# ============================================================
df_sem <- readRDS(file.path(ROOT, "results", "07_6_prepare_sem_bnlearn_data",
                             "data", "df_cont.rds"))
log_info("df_cont: {nrow(df_sem)} x {ncol(df_sem)}")

# Standardize continuous variables
z_vars <- c("age","RIN","PMI",
            "OXPHOS13","OXPHOS_pathway",
            "MEgreenyellow","MEblack",
            "ComplexIII_mean","UQCRQ_z")
for (v in z_vars) {
  if (v %in% colnames(df_sem)) {
    df_sem[[paste0("z_", v)]] <- as.numeric(scale(df_sem[[v]]))
  }
}

# Already z-scored indicators: use as-is
if (all(c("UQCRFS1_gene_z","UQCRH_gene_z") %in% colnames(df_sem))) {
  df_sem$UQCRQ_ind   <- df_sem$UQCRQ_z
  df_sem$UQCRFS1_ind <- df_sem$UQCRFS1_gene_z
  df_sem$UQCRH_ind   <- df_sem$UQCRH_gene_z
} else {
  df_sem$UQCRQ_ind <- df_sem$UQCRFS1_ind <- df_sem$UQCRH_ind <- NA_real_
}

df_sem$sex     <- factor(df_sem$sex, levels = c("Female", "Male"))
df_sem$z_sex   <- as.integer(df_sem$sex == "Male")
df_sem$z_grade <- as.numeric(scale(df_sem$grade))

labels_main <- c("a","b","cprime","indirect","total","prop_med")

# ============================================================
# 1. OXPHOS: Model A
# ============================================================
log_info("--- 1. OXPHOS Model A ---")

modelA <- '
  z_OXPHOS13 ~ a * z_age + z_sex + z_RIN + z_PMI
  z_grade    ~ b * z_OXPHOS13 + cprime * z_age + z_sex + z_RIN + z_PMI
  indirect := a * b
  total    := cprime + a * b
  prop_med := (a * b) / total
'
fitA <- sem(modelA, data = df_sem, se = "bootstrap", bootstrap = 1000)
selA <- safe_estimates(fitA, labels_main)
write_tab(selA, "model_A_estimates.csv")
print(selA, row.names = FALSE, digits = 4)

# ============================================================
# 2. OXPHOS: Model B (latent)
# ============================================================
log_info("--- 2. OXPHOS Model B ---")

modelB <- '
  OXPHOS =~ z_OXPHOS13 + z_OXPHOS_pathway + z_MEgreenyellow
  OXPHOS ~ a * z_age + z_sex + z_RIN + z_PMI
  z_grade ~ b * OXPHOS + cprime * z_age + z_sex + z_RIN + z_PMI
  indirect := a * b
  total    := cprime + a * b
  prop_med := (a * b) / total
'
fitB <- tryCatch(
  sem(modelB, data = df_sem, se = "bootstrap", bootstrap = 1000),
  warning = function(w) {
    log_warn("Model B warning: {conditionMessage(w)}")
    suppressWarnings(sem(modelB, data = df_sem,
                          se = "bootstrap", bootstrap = 1000))
  },
  error = function(e) { log_err("Model B failed: {conditionMessage(e)}"); NULL }
)
if (!is.null(fitB)) {
  fmB <- fitMeasures(fitB, c("cfi","rmsea","srmr","chisq","df","pvalue"))
  write_tab(data.frame(metric = names(fmB), value = as.numeric(fmB)),
            "model_B_fit.csv")
}
selB <- safe_estimates(fitB, labels_main)
write_tab(selB, "model_B_estimates.csv")
print(selB, row.names = FALSE, digits = 4)

# ============================================================
# 3. OXPHOS: Model C
# ============================================================
log_info("--- 3. OXPHOS Model C ---")
modelC <- '
  z_grade ~ b * z_OXPHOS13 + c * z_age + z_sex + z_RIN + z_PMI
'
fitC <- sem(modelC, data = df_sem, se = "bootstrap", bootstrap = 1000)
selC <- safe_estimates(fitC, c("b","c"))
write_tab(selC, "model_C_estimates.csv")
print(selC, row.names = FALSE, digits = 4)

# ============================================================
# 4. OXPHOS: Model D (reverse)
# ============================================================
log_info("--- 4. OXPHOS Model D ---")
modelD <- '
  z_OXPHOS13 ~ b * z_grade + z_age + z_sex + z_RIN + z_PMI
'
fitD <- sem(modelD, data = df_sem, se = "bootstrap", bootstrap = 1000)
selD <- safe_estimates(fitD, "b")
write_tab(selD, "model_D_reverse.csv")
print(selD, row.names = FALSE, digits = 4)

# ============================================================
# 5. ComplexIII_mean: Model A
# ============================================================
log_info("--- 5. ComplexIII_mean Model A ---")
modelA_cIII <- '
  z_ComplexIII_mean ~ a * z_age + z_sex + z_RIN + z_PMI
  z_grade           ~ b * z_ComplexIII_mean + cprime * z_age + z_sex + z_RIN + z_PMI
  indirect := a * b
  total    := cprime + a * b
  prop_med := (a * b) / total
'
fitA_cIII <- sem(modelA_cIII, data = df_sem, se = "bootstrap", bootstrap = 1000)
selA_cIII <- safe_estimates(fitA_cIII, labels_main)
write_tab(selA_cIII, "modelA_ComplexIII_mean.csv")
print(selA_cIII, row.names = FALSE, digits = 4)

# ============================================================
# 6. ComplexIII_mean: Model B (latent)
# ============================================================
log_info("--- 6. ComplexIII_mean Model B (latent) ---")

if (all(!is.na(df_sem$UQCRQ_ind)) &&
    all(!is.na(df_sem$UQCRFS1_ind)) &&
    all(!is.na(df_sem$UQCRH_ind))) {

  modelB_cIII <- '
    CIII =~ UQCRQ_ind + UQCRFS1_ind + UQCRH_ind
    CIII ~ a * z_age + z_sex + z_RIN + z_PMI
    z_grade ~ b * CIII + cprime * z_age + z_sex + z_RIN + z_PMI
    indirect := a * b
    total    := cprime + a * b
    prop_med := (a * b) / total
  '
  fitB_cIII <- tryCatch(
    sem(modelB_cIII, data = df_sem, se = "bootstrap", bootstrap = 1000),
    error = function(e) { log_err("Model B cIII failed: {conditionMessage(e)}"); NULL }
  )
  if (!is.null(fitB_cIII)) {
    fmB_cIII <- fitMeasures(fitB_cIII,
                            c("cfi","rmsea","srmr","chisq","df","pvalue"))
    write_tab(data.frame(metric = names(fmB_cIII),
                          value = as.numeric(fmB_cIII)),
              "modelB_fit_ComplexIII_latent.csv")
  }
  selB_cIII <- safe_estimates(fitB_cIII, labels_main)
} else {
  log_warn("Latent Complex III indicators unavailable; skipping Model B.")
  selB_cIII <- safe_estimates(NULL, labels_main)
}
write_tab(selB_cIII, "modelB_ComplexIII_latent.csv")
print(selB_cIII, row.names = FALSE, digits = 4)

# ============================================================
# 7. ComplexIII_mean: Model C
# ============================================================
log_info("--- 7. ComplexIII_mean Model C ---")
modelC_cIII <- '
  z_grade ~ b * z_ComplexIII_mean + c * z_age + z_sex + z_RIN + z_PMI
'
fitC_cIII <- sem(modelC_cIII, data = df_sem, se = "bootstrap", bootstrap = 1000)
selC_cIII <- safe_estimates(fitC_cIII, c("b","c"))
write_tab(selC_cIII, "modelC_ComplexIII_mean_to_grade.csv")
print(selC_cIII, row.names = FALSE, digits = 4)

# ============================================================
# 8. ComplexIII_mean: Model D (reverse)
# ============================================================
log_info("--- 8. ComplexIII_mean Model D ---")
modelD_cIII <- '
  z_ComplexIII_mean ~ b * z_grade + z_age + z_sex + z_RIN + z_PMI
'
fitD_cIII <- sem(modelD_cIII, data = df_sem, se = "bootstrap", bootstrap = 1000)
selD_cIII <- safe_estimates(fitD_cIII, "b")
write_tab(selD_cIII, "modelD_reverse_ComplexIII_mean.csv")
print(selD_cIII, row.names = FALSE, digits = 4)

# ============================================================
# 9. UQCRQ single-gene: Model A / C / D
# ============================================================
log_info("--- 9. UQCRQ single-gene models ---")

modelA_uqcrq <- '
  z_UQCRQ_z ~ a * z_age + z_sex + z_RIN + z_PMI
  z_grade   ~ b * z_UQCRQ_z + cprime * z_age + z_sex + z_RIN + z_PMI
  indirect := a * b
  total    := cprime + a * b
  prop_med := (a * b) / total
'
fitA_uqcrq <- sem(modelA_uqcrq, data = df_sem, se = "bootstrap", bootstrap = 1000)
selA_uqcrq <- safe_estimates(fitA_uqcrq, labels_main)
write_tab(selA_uqcrq, "modelA_UQCRQ.csv")

modelC_uqcrq <- '
  z_grade ~ b * z_UQCRQ_z + c * z_age + z_sex + z_RIN + z_PMI
'
fitC_uqcrq <- sem(modelC_uqcrq, data = df_sem, se = "bootstrap", bootstrap = 1000)
selC_uqcrq <- safe_estimates(fitC_uqcrq, c("b","c"))
write_tab(selC_uqcrq, "modelC_UQCRQ_to_grade.csv")

modelD_uqcrq <- '
  z_UQCRQ_z ~ b * z_grade + z_age + z_sex + z_RIN + z_PMI
'
fitD_uqcrq <- sem(modelD_uqcrq, data = df_sem, se = "bootstrap", bootstrap = 1000)
selD_uqcrq <- safe_estimates(fitD_uqcrq, "b")
write_tab(selD_uqcrq, "modelD_reverse_UQCRQ.csv")

# ============================================================
# 10. Comparison table
# ============================================================
log_info("--- 10. Comparison ---")

get_est <- function(df, lbl, col) {
  v <- df[[col]][df$label == lbl]
  if (length(v) == 0) NA_real_ else v[1]
}

comp <- data.frame(
  model = c("OXPHOS_A_indirect",
             "OXPHOS_B_indirect_latent",
             "OXPHOS_C_direct",
             "OXPHOS_D_reverse",
             "CIII_A_indirect",
             "CIII_B_indirect_latent",
             "CIII_C_direct",
             "CIII_D_reverse",
             "UQCRQ_A_indirect",
             "UQCRQ_C_direct",
             "UQCRQ_D_reverse"),
  estimate = c(get_est(selA, "indirect", "est"),
                get_est(selB, "indirect", "est"),
                get_est(selC, "b", "est"),
                get_est(selD, "b", "est"),
                get_est(selA_cIII, "indirect", "est"),
                get_est(selB_cIII, "indirect", "est"),
                get_est(selC_cIII, "b", "est"),
                get_est(selD_cIII, "b", "est"),
                get_est(selA_uqcrq, "indirect", "est"),
                get_est(selC_uqcrq, "b", "est"),
                get_est(selD_uqcrq, "b", "est")),
  ci_low  = c(get_est(selA, "indirect", "ci.lower"),
               get_est(selB, "indirect", "ci.lower"),
               get_est(selC, "b", "ci.lower"),
               get_est(selD, "b", "ci.lower"),
               get_est(selA_cIII, "indirect", "ci.lower"),
               get_est(selB_cIII, "indirect", "ci.lower"),
               get_est(selC_cIII, "b", "ci.lower"),
               get_est(selD_cIII, "b", "ci.lower"),
               get_est(selA_uqcrq, "indirect", "ci.lower"),
               get_est(selC_uqcrq, "b", "ci.lower"),
               get_est(selD_uqcrq, "b", "ci.lower")),
  ci_high = c(get_est(selA, "indirect", "ci.upper"),
               get_est(selB, "indirect", "ci.upper"),
               get_est(selC, "b", "ci.upper"),
               get_est(selD, "b", "ci.upper"),
               get_est(selA_cIII, "indirect", "ci.upper"),
               get_est(selB_cIII, "indirect", "ci.upper"),
               get_est(selC_cIII, "b", "ci.upper"),
               get_est(selD_cIII, "b", "ci.upper"),
               get_est(selA_uqcrq, "indirect", "ci.upper"),
               get_est(selC_uqcrq, "b", "ci.upper"),
               get_est(selD_uqcrq, "b", "ci.upper")),
  p       = c(get_est(selA, "indirect", "pvalue"),
               get_est(selB, "indirect", "pvalue"),
               get_est(selC, "b", "pvalue"),
               get_est(selD, "b", "pvalue"),
               get_est(selA_cIII, "indirect", "pvalue"),
               get_est(selB_cIII, "indirect", "pvalue"),
               get_est(selC_cIII, "b", "pvalue"),
               get_est(selD_cIII, "b", "pvalue"),
               get_est(selA_uqcrq, "indirect", "pvalue"),
               get_est(selC_uqcrq, "b", "pvalue"),
               get_est(selD_uqcrq, "b", "pvalue"))
)
comp$sig_ci <- ifelse(!is.na(comp$ci_low) & !is.na(comp$ci_high) &
                       comp$ci_low * comp$ci_high > 0, "yes", "no")
write_tab(comp, "model_comparison.csv")

cat("\n=== Model comparison ===\n")
print(comp, row.names = FALSE, digits = 4)

# ============================================================
# 11. Forest plot
# ============================================================
log_info("--- 11. Figure ---")

plot_df <- comp[!is.na(comp$estimate), ]
plot_df$model <- factor(plot_df$model, levels = rev(plot_df$model))

p <- ggplot(plot_df, aes(x = estimate, y = model, color = sig_ci)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high),
                  height = 0.2, linewidth = 0.8) +
  geom_point(size = 3.5) +
  scale_color_manual(values = c("yes" = "#C44E52", "no" = "#B0B0B0"),
                      name = "CI excludes 0") +
  theme_bw(12) +
  labs(x = "Standardized effect (95% bootstrap CI)", y = "",
       title = "SEM mediation: OXPHOS and Complex III in AMD") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave(file.path(FIG, "Fig1_sem_forest.pdf"), p, width = 9, height = 6)
ggsave(file.path(FIG, "Fig1_sem_forest.png"), p, width = 9, height = 6, dpi = 300)

# ============================================================
# 12. Summary
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("=== 07_6b summary ===\n")
cat(strrep("=", 70), "\n", sep = "")

cat("\n[OXPHOS Model A]\n"); print(selA, row.names = FALSE, digits = 4)
cat("\n[ComplexIII_mean Model A]\n"); print(selA_cIII, row.names = FALSE, digits = 4)
cat("\n[ComplexIII_mean Model C]\n"); print(selC_cIII, row.names = FALSE, digits = 4)
cat("\n[ComplexIII_mean Model D]\n"); print(selD_cIII, row.names = FALSE, digits = 4)
cat("\n[UQCRQ Model A]\n"); print(selA_uqcrq, row.names = FALSE, digits = 4)
cat("\n[UQCRQ Model C]\n"); print(selC_uqcrq, row.names = FALSE, digits = 4)
cat("\n[UQCRQ Model D]\n"); print(selD_uqcrq, row.names = FALSE, digits = 4)
cat("\n[Comparison]\n"); print(comp, row.names = FALSE, digits = 4)

cat("\nOUT:", OUT, "\n")
log_info("=== 07_6b done ===")

