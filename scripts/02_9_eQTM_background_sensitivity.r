# ============================================================
# 02_9_eQTM_background_sensitivity.r（v5）
# ============================================================
# ============================================================

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr)
})

ROOT <- "E:/retina_tools/aging_rgc_axon/AMD"
OUT  <- file.path(ROOT, "results", "02_9_eQTM_background_sensitivity")
TAB  <- file.path(OUT, "tables")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)

log_info("=== 02_9 eQTM background sensitivity (v5) ===")

# ------------------------------------------------------------
# 1. Backgroundsensitivity
# ------------------------------------------------------------
eqtm <- read_csv(file.path(ROOT,
  "results/07_5c_OXPHOS_analysis/tables/fisher_eqtm_enrichment.csv"),
  show_col_types = FALSE)
cat("\neQTM enrichment 原始结果：\n"); print(eqtm)

n_eqtm_targets <- 241
n_genes_expr   <- 14233
N_genome       <- 20000

eqtm_sens <- eqtm %>%
  mutate(
    a = round(pct_hit / 100 * n_total),
    b = n_total - a,
    c = n_eqtm_targets - a,
    d_20000 = N_genome - n_total - c,
    d_14233 = n_genes_expr - n_total - c,
    OR_20000 = (a * d_20000) / (b * c),
    OR_14233 = (a * d_14233) / (b * c),
    OR_ratio    = OR_20000 / OR_14233,
    log_OR_diff = log(OR_20000 / OR_14233)
  )

cat("\n背景敏感性（OR_ratio = OR_20000 / OR_14233）：\n")
print(eqtm_sens[, c("target", "n_total", "n_hit",
                    "OR_20000", "OR_14233",
                    "OR_ratio", "log_OR_diff")])

write_csv(eqtm_sens, file.path(TAB, "eqtm_background_sensitivity.csv"))

log_info(paste0("Median OR ratio = ",
                round(median(eqtm_sens$OR_ratio, na.rm = TRUE), 3)))
log_info(paste0("Max |log OR diff| = ",
                round(max(abs(eqtm_sens$log_OR_diff), na.rm = TRUE), 3)))

contingency <- data.frame(
  target                = eqtm_sens$target,
  `In query, eQTM+`     = eqtm_sens$a,
  `In query, eQTM-`     = eqtm_sens$b,
  `Not query, eQTM+`    = eqtm_sens$c,
  `Not query, eQTM-`    = eqtm_sens$d_20000,
  check.names = FALSE
)
write_csv(contingency, file.path(TAB, "eqtm_contingency_tables.csv"))
cat("\n完整列联表（20,000 背景）：\n"); print(contingency)

# ------------------------------------------------------------
# ------------------------------------------------------------
mod_file <- file.path(ROOT,
  "results/07_5_GSE231536_methylation/tables/module_direction_validation.csv")

if (file.exists(mod_file)) {
  mod_df <- read_csv(mod_file, show_col_types = FALSE)
  cat("\nmodule_direction_validation.csv：\n"); print(mod_df)

  gy_row <- mod_df %>% filter(module == "greenyellow")
  if (nrow(gy_row) == 1) {
    log_info(paste0("greenyellow 方向一致性: n = ", gy_row$n_cpg,
                    ", concordant = ", gy_row$n_consistent,
                    " (", round(gy_row$pct, 1), "%), binom_p = ",
                    signif(gy_row$binom_p, 3)))
  }

  write_csv(mod_df, file.path(TAB, "direction_consistency_by_module.csv"))
  cat("\n方向一致性（全部模块）：\n"); print(mod_df)

} else {
  log_info("module_direction_validation.csv 不存在，跳过方向一致性")
}

# ------------------------------------------------------------
# ------------------------------------------------------------
cpg_file_16 <- file.path(ROOT,
  "results/07_5_GSE231536_methylation/tables/gene16_MGS1_vs_MGS3_eQTM.csv")
if (file.exists(cpg_file_16)) {
  cpg16  <- read_csv(cpg_file_16, show_col_types = FALSE)
  n16    <- nrow(cpg16)
  conc16 <- sum(cpg16$consistent, na.rm = TRUE)
  binom16 <- binom.test(conc16, n16, p = 0.5)
  log_info(paste0("16-gene signature: n = ", n16,
                  ", concordant = ", conc16,
                  " (", round(conc16 / n16 * 100, 1), "%), p = ",
                  signif(binom16$p.value, 3)))
}

# ------------------------------------------------------------
# ------------------------------------------------------------
detail_file <- file.path(ROOT,
  "results/07_5_GSE231536_methylation/tables/module_direction_with_bulk.csv")
if (file.exists(detail_file)) {
  cat("\n=== module_direction_with_bulk.csv ===\n")
  detail <- read_csv(detail_file, show_col_types = FALSE)
  cat("列名：", paste(colnames(detail), collapse = ", "), "\n")
  cat("行数：", nrow(detail), "\n")
  print(head(detail, 5))

  if (all(c("module", "gene_sym", "consistent") %in% colnames(detail))) {
    gy_detail <- detail %>% filter(module == "greenyellow")
    if (nrow(gy_detail) > 0) {
      log_info(paste0("greenyellow CpG 明细: n = ", nrow(gy_detail),
                      ", 唯一基因 = ", length(unique(gy_detail$gene_sym))))

      if ("p_FDR" %in% colnames(gy_detail)) {
        gy_rep <- gy_detail %>%
          group_by(gene_sym) %>%
          slice_min(p_FDR, n = 1, with_ties = FALSE) %>%
          ungroup()
      } else {
        gy_rep <- gy_detail %>%
          group_by(gene_sym) %>% slice(1) %>% ungroup()
      }
      n_rep   <- nrow(gy_rep)
      conc_rep <- sum(gy_rep$consistent, na.rm = TRUE)
      binom_rep <- binom.test(conc_rep, n_rep, p = 0.5)
      log_info(paste0("Independent CpG (1 per gene): n = ", n_rep,
                      ", concordant = ", conc_rep,
                      " (", round(conc_rep / n_rep * 100, 1),
                      "%), p = ", signif(binom_rep$p.value, 3)))

      indep_df <- data.frame(
        level      = "Independent CpG (1 per gene)",
        n          = n_rep,
        concordant = conc_rep,
        pct        = conc_rep / n_rep * 100,
        binom_p    = binom_rep$p.value
      )
      write_csv(indep_df,
                file.path(TAB, "direction_consistency_independent_cpg.csv"))
    }
  }
}

log_info("=== 02_9 v5 done ===")

