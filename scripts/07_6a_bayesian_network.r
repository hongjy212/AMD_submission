# ============================================================
# 07_6a_bayesian_network.r
# Purpose: BN for OXPHOS and Complex III
#   - Analysis A: data-driven (no whitelist)
#   - Analysis B: with OXPHOS biology priors
#   - Analysis B2: with Complex III biology priors
#   - Sensitivity: without age
# Input:  results/07_6_prepare_sem_bnlearn_data/data/df_disc.rds
# Output: results/07_6a_bayesian_network/
# ============================================================
source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(ggplot2); library(bnlearn)
})

ROOT <- "E:/retina_tools/aging_rgc_axon/AMD"
OUT  <- file.path(ROOT, "results", "07_6a_bayesian_network")
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

log_info("=== 07_6a Bayesian network ===")

# ============================================================
# 0. Load
# ============================================================
df_disc <- readRDS(file.path(ROOT, "results", "07_6_prepare_sem_bnlearn_data",
                              "data", "df_disc.rds"))
log_info("df_disc: {nrow(df_disc)} x {ncol(df_disc)}")
log_info("Variables: {paste(colnames(df_disc), collapse=', ')}")

N_BOOT <- 1000
set.seed(SEED)

bl_age_sink <- data.frame(
  from = setdiff(colnames(df_disc), "age"),
  to   = rep("age", length(setdiff(colnames(df_disc), "age")))
)

# ============================================================
# 1. Analysis A: data-driven
# ============================================================
log_info("--- 1. Analysis A: data-driven ---")

dag_A <- hc(df_disc, score = "bic", blacklist = bl_age_sink,
             restart = 10, perturb = 5)
log_info("A arcs: {narcs(dag_A)}")

strength_A <- boot.strength(
  df_disc, R = N_BOOT, algorithm = "hc",
  algorithm.args = list(score = "bic", blacklist = bl_age_sink,
                         restart = 5, perturb = 3)
)
strength_A <- strength_A[order(-strength_A$strength), ]
rownames(strength_A) <- NULL
write_tab(strength_A, "bn_strength_data_driven.csv")

strong_A <- strength_A[strength_A$strength >= 0.85, ]
cat("\n=== A: strong edges (>=0.85) ===\n")
print(head(strong_A, 25), row.names = FALSE, digits = 3)

# ============================================================
# 2. Analysis B: with OXPHOS priors
# ============================================================
log_info("--- 2. Analysis B: OXPHOS priors ---")

wl <- data.frame(
  from = c("age", "age", "age",
            "MEgreenyellow", "MEgreenyellow",
            "OXPHOS_pathway",
            "OXPHOS_pathway",
            "OXPHOS13",
            "MEblack"),
  to   = c("grade", "OXPHOS13", "OXPHOS_pathway",
            "OXPHOS13", "OXPHOS_pathway",
            "OXPHOS13",
            "grade",
            "grade",
            "MEgreen")
)
log_info("Whitelist (OXPHOS): {nrow(wl)} edges")

bl_B <- rbind(
  bl_age_sink,
  data.frame(
    from = rep("grade", ncol(df_disc) - 1),
    to   = setdiff(colnames(df_disc), "grade")
  )
)

dag_B <- hc(df_disc, score = "bic", whitelist = wl, blacklist = bl_B,
             restart = 10, perturb = 5)
log_info("B arcs: {narcs(dag_B)}")

strength_B <- boot.strength(
  df_disc, R = N_BOOT, algorithm = "hc",
  algorithm.args = list(score = "bic", whitelist = wl, blacklist = bl_B,
                         restart = 5, perturb = 3)
)
strength_B <- strength_B[order(-strength_B$strength), ]
rownames(strength_B) <- NULL
write_tab(strength_B, "bn_strength_with_prior.csv")

strong_B <- strength_B[strength_B$strength >= 0.85, ]
cat("\n=== B: strong edges (>=0.85) ===\n")
print(head(strong_B, 25), row.names = FALSE, digits = 3)

# ============================================================
# 3. Analysis B2: with Complex III priors
# ============================================================
log_info("--- 3. Analysis B2: Complex III priors ---")

wl_cIII <- data.frame(
  from = c("age", "age", "age",
            "ComplexIII_mean", "ComplexIII_mean",
            "UQCRQ_z"),
  to   = c("grade", "ComplexIII_mean", "UQCRQ_z",
            "grade", "UQCRQ_z",
            "grade")
)
log_info("Whitelist (Complex III): {nrow(wl_cIII)} edges")

dag_B2 <- hc(df_disc, score = "bic", whitelist = wl_cIII, blacklist = bl_B,
              restart = 10, perturb = 5)
log_info("B2 arcs: {narcs(dag_B2)}")

strength_B2 <- boot.strength(
  df_disc, R = N_BOOT, algorithm = "hc",
  algorithm.args = list(score = "bic", whitelist = wl_cIII, blacklist = bl_B,
                         restart = 5, perturb = 3)
)
strength_B2 <- strength_B2[order(-strength_B2$strength), ]
rownames(strength_B2) <- NULL
write_tab(strength_B2, "bn_strength_prior_cIII.csv")

# ============================================================
# 4. Stable edges (OXPHOS)
# ============================================================
log_info("--- 4. Stable edges (OXPHOS) ---")

comp <- merge(
  strength_A[, c("from","to","strength","direction")],
  strength_B[, c("from","to","strength","direction")],
  by = c("from", "to"), all = TRUE,
  suffixes = c("_data", "_prior")
)
colnames(comp) <- c("from", "to",
                     "strength_data", "direction_data",
                     "strength_prior", "direction_prior")

comp$stable <- !is.na(comp$strength_data) & !is.na(comp$strength_prior) &
               comp$strength_data >= 0.85 & comp$strength_prior >= 0.85

comp <- comp[order(-comp$stable,
                    -pmax(comp$strength_data, comp$strength_prior,
                          na.rm = TRUE)), ]
write_tab(comp, "bn_stable_edges.csv")

stable_only <- comp[comp$stable, ]
cat("\n=== Stable edges (OXPHOS) ===\n")
if (nrow(stable_only) > 0) {
  print(stable_only[, c("from","to","strength_data","strength_prior",
                          "direction_data","direction_prior")],
        row.names = FALSE, digits = 3)
} else {
  cat("  (none)\n")
}

# ============================================================
# 4b. Stable edges (Complex III)
# ============================================================
log_info("--- 4b. Stable edges (Complex III) ---")

comp_cIII <- merge(
  strength_A[, c("from","to","strength","direction")],
  strength_B2[, c("from","to","strength","direction")],
  by = c("from", "to"), all = TRUE,
  suffixes = c("_data", "_prior")
)
colnames(comp_cIII) <- c("from", "to",
                          "strength_data", "direction_data",
                          "strength_prior", "direction_prior")

comp_cIII$stable <- !is.na(comp_cIII$strength_data) &
                    !is.na(comp_cIII$strength_prior) &
                    comp_cIII$strength_data >= 0.85 &
                    comp_cIII$strength_prior >= 0.85

comp_cIII <- comp_cIII[order(-comp_cIII$stable,
                              -pmax(comp_cIII$strength_data,
                                    comp_cIII$strength_prior,
                                    na.rm = TRUE)), ]
write_tab(comp_cIII, "bn_stable_edges_cIII.csv")

stable_cIII <- comp_cIII[comp_cIII$stable, ]
cat("\n=== Stable edges (Complex III) ===\n")
if (nrow(stable_cIII) > 0) {
  print(stable_cIII[, c("from","to","strength_data","strength_prior",
                          "direction_data","direction_prior")],
        row.names = FALSE, digits = 3)
} else {
  cat("  (none)\n")
}

# ============================================================
# 5. Edges into grade
# ============================================================
log_info("--- 5. Edges into grade ---")

grade_edges <- comp[comp$to == "grade", ]
grade_edges <- grade_edges[order(-grade_edges$strength_prior,
                                   -grade_edges$strength_data), ]
write_tab(grade_edges, "bn_edges_into_grade.csv")

cat("\n=== Edges into grade (OXPHOS priors) ===\n")
print(grade_edges[, c("from","strength_data","strength_prior",
                        "direction_data","direction_prior","stable")],
      row.names = FALSE, digits = 3)

grade_edges_cIII <- comp_cIII[comp_cIII$to == "grade", ]
grade_edges_cIII <- grade_edges_cIII[order(-grade_edges_cIII$strength_prior,
                                             -grade_edges_cIII$strength_data), ]
write_tab(grade_edges_cIII, "bn_edges_into_grade_cIII.csv")

cat("\n=== Edges into grade (Complex III priors) ===\n")
print(grade_edges_cIII[, c("from","strength_data","strength_prior",
                             "direction_data","direction_prior","stable")],
      row.names = FALSE, digits = 3)

# ============================================================
# 6. Averaged network (OXPHOS priors)
# ============================================================
log_info("--- 6. Averaged network (OXPHOS priors) ---")
avg_net <- averaged.network(strength_B, threshold = 0.5)
log_info("Averaged network: {narcs(avg_net)} arcs")

arcs_avg <- data.frame(arcs(avg_net))
if (nrow(arcs_avg) > 0) {
  colnames(arcs_avg) <- c("from", "to")
  key <- paste(arcs_avg$from, arcs_avg$to)
  arcs_avg$strength <- strength_B$strength[match(key,
                                                  paste(strength_B$from,
                                                        strength_B$to))]
  arcs_avg$direction <- strength_B$direction[match(key,
                                                    paste(strength_B$from,
                                                          strength_B$to))]
  arcs_avg <- arcs_avg[order(-arcs_avg$strength), ]
}
write_tab(arcs_avg, "bn_averaged_edges.csv")

# ============================================================
# 7. Sensitivity: EXCLUDE AGE from data
# ============================================================
log_info("--- 7. Sensitivity: exclude age ---")

df_noage <- df_disc[, setdiff(colnames(df_disc), "age"), drop = FALSE]
log_info("df_noage: {nrow(df_noage)} x {ncol(df_noage)}")

bl_noage <- data.frame(
  from = rep("grade", ncol(df_noage) - 1),
  to   = setdiff(colnames(df_noage), "grade")
)

dag_noage <- hc(df_noage, score = "bic", blacklist = bl_noage,
                 restart = 10, perturb = 5)
log_info("DAG (no age) arcs: {narcs(dag_noage)}")

strength_noage <- boot.strength(
  df_noage, R = N_BOOT, algorithm = "hc",
  algorithm.args = list(score = "bic", blacklist = bl_noage,
                         restart = 5, perturb = 3)
)
strength_noage <- strength_noage[order(-strength_noage$strength), ]
rownames(strength_noage) <- NULL
write_tab(strength_noage, "bn_strength_no_age.csv")

grade_edges_noage <- strength_noage[strength_noage$to == "grade", ]
grade_edges_noage <- grade_edges_noage[order(-grade_edges_noage$strength), ]
write_tab(grade_edges_noage, "bn_edges_into_grade_no_age.csv")

cat("\n=== Sensitivity: edges into grade WITHOUT age ===\n")
print(grade_edges_noage[, c("from","to","strength","direction")],
      row.names = FALSE, digits = 3)

with_age <- strength_A[strength_A$to == "grade",
                        c("from","strength","direction")]
cat("\n=== Comparison: edges into grade (with age, data-driven) ===\n")
print(with_age, row.names = FALSE, digits = 3)

# ============================================================
# 8. Figures
# ============================================================
log_info("--- 8. Figures ---")

# Fig1: OXPHOS priors vs data-driven, top 20
top_edges <- comp[order(-pmax(comp$strength_data, comp$strength_prior,
                                na.rm = TRUE)), ]
top_edges <- head(top_edges, 20)
top_edges$edge <- paste(top_edges$from, "->", top_edges$to)
top_edges$edge <- factor(top_edges$edge, levels = rev(top_edges$edge))

plot_long <- top_edges %>%
  select(edge, strength_data, strength_prior) %>%
  pivot_longer(-edge, names_to = "analysis", values_to = "strength")
plot_long$analysis <- factor(plot_long$analysis,
                              levels = c("strength_data", "strength_prior"),
                              labels = c("Data-driven", "With OXPHOS priors"))

p1 <- ggplot(plot_long, aes(x = strength, y = edge, fill = analysis)) +
  geom_col(position = "dodge") +
  geom_vline(xintercept = 0.85, linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("Data-driven" = "#4C72B0",
                                "With OXPHOS priors" = "#C44E52")) +
  theme_bw(12) +
  labs(x = "Bootstrap strength", y = "",
       title = "Top 20 edges: data-driven vs OXPHOS priors",
       fill = "") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave(file.path(FIG, "Fig1_strength_compare.pdf"), p1, width = 9, height = 6)

# Fig2: Complex III priors vs data-driven
top_edges_cIII <- comp_cIII[order(-pmax(comp_cIII$strength_data,
                                         comp_cIII$strength_prior,
                                         na.rm = TRUE)), ]
top_edges_cIII <- head(top_edges_cIII, 20)
top_edges_cIII$edge <- paste(top_edges_cIII$from, "->", top_edges_cIII$to)
top_edges_cIII$edge <- factor(top_edges_cIII$edge,
                               levels = rev(top_edges_cIII$edge))

plot_long_cIII <- top_edges_cIII %>%
  select(edge, strength_data, strength_prior) %>%
  pivot_longer(-edge, names_to = "analysis", values_to = "strength")
plot_long_cIII$analysis <- factor(plot_long_cIII$analysis,
                                    levels = c("strength_data", "strength_prior"),
                                    labels = c("Data-driven", "With Complex III priors"))

p2 <- ggplot(plot_long_cIII, aes(x = strength, y = edge, fill = analysis)) +
  geom_col(position = "dodge") +
  geom_vline(xintercept = 0.85, linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("Data-driven" = "#4C72B0",
                                "With Complex III priors" = "#55A868")) +
  theme_bw(12) +
  labs(x = "Bootstrap strength", y = "",
       title = "Top 20 edges: data-driven vs Complex III priors",
       fill = "") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave(file.path(FIG, "Fig2_strength_compare_cIII.pdf"), p2, width = 9, height = 6)

# Fig3: age-exclusion comparison
cmp_df <- rbind(
  data.frame(setting = "With age",
             from = with_age$from,
             strength = with_age$strength),
  data.frame(setting = "Without age",
             from = grade_edges_noage$from,
             strength = grade_edges_noage$strength)
)
cmp_df$setting <- factor(cmp_df$setting,
                          levels = c("With age", "Without age"))

p3 <- ggplot(cmp_df, aes(x = strength, y = reorder(from, strength),
                           fill = setting)) +
  geom_col(position = "dodge") +
  geom_vline(xintercept = 0.85, linetype = "dashed", color = "red") +
  scale_fill_manual(values = c("With age" = "#4C72B0",
                                "Without age" = "#C44E52")) +
  theme_bw(12) +
  labs(x = "Bootstrap strength", y = "Predictor",
       title = "Edges into grade: with vs without age",
       fill = "") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave(file.path(FIG, "Fig3_age_exclusion.pdf"), p3, width = 8, height = 5)

# ============================================================
# 9. Summary
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("=== 07_6a summary ===\n")
cat(strrep("=", 70), "\n", sep = "")

cat("\n[A: data-driven]\n")
cat("Arcs:", narcs(dag_A), "  Strong:", nrow(strong_A), "\n")

cat("\n[B: OXPHOS priors]\n")
cat("Arcs:", narcs(dag_B), "  Strong:", nrow(strong_B), "\n")

cat("\n[B2: Complex III priors]\n")
cat("Arcs:", narcs(dag_B2), "\n")

cat("\n[Stable edges: OXPHOS]\n")
if (nrow(stable_only) > 0) {
  print(stable_only[, c("from","to","strength_data","strength_prior")],
        row.names = FALSE, digits = 3)
} else {
  cat("  (none)\n")
}

cat("\n[Stable edges: Complex III]\n")
if (nrow(stable_cIII) > 0) {
  print(stable_cIII[, c("from","to","strength_data","strength_prior")],
        row.names = FALSE, digits = 3)
} else {
  cat("  (none)\n")
}

cat("\n[Edges into grade: OXPHOS priors]\n")
print(grade_edges[, c("from","strength_data","strength_prior","stable")],
      row.names = FALSE, digits = 3)

cat("\n[Edges into grade: Complex III priors]\n")
print(grade_edges_cIII[, c("from","strength_data","strength_prior","stable")],
      row.names = FALSE, digits = 3)

cat("\n[Edges into grade: without age]\n")
print(grade_edges_noage[, c("from","strength","direction")],
      row.names = FALSE, digits = 3)

cat("\nOUT:", OUT, "\n")
log_info("=== 07_6a done ===")
