# ============================================================
# 08_4_fdr_positive_control.r
# Output:  results/08_4_fdr_positive_control/
# ============================================================

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(data.table); library(dplyr)
})

ROOT     <- "E:/retina_tools/aging_rgc_axon/AMD"
DIR_THIS <- file.path(ROOT, "results", "08_4_fdr_positive_control")
DIR_FDR  <- file.path(DIR_THIS, "01_fdr")
DIR_POS  <- file.path(DIR_THIS, "02_positive_control")
for (d in c(DIR_FDR, DIR_POS)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

IN_MR  <- file.path(ROOT, "results/08_3_mr/01_mr")
IN_ROB <- file.path(ROOT, "results/08_3_mr/02_robust")

log_info("=== 08_4 fdr + positive control ===")

all_pvals <- list()
if (file.exists(file.path(IN_MR, "MR_allstage.csv"))) {
  all_pvals$allstage <- fread(file.path(IN_MR, "MR_allstage.csv"))[, source := "main_allstage"]
}
if (file.exists(file.path(IN_MR, "MR_early.csv"))) {
  all_pvals$early <- fread(file.path(IN_MR, "MR_early.csv"))[, source := "main_early"]
}
if (file.exists(file.path(IN_ROB, "robust_methods.csv"))) {
  rob <- fread(file.path(IN_ROB, "robust_methods.csv"))
  rob <- rob[method != "Inverse variance weighted"]
  rob[, `:=`(source = "robust", p_MR = pval, OR_MR = exp(b),
             OR_lower = exp(b - 1.96 * se), OR_upper = exp(b + 1.96 * se),
             outcome = "early")]
  all_pvals$robust <- rob[, .(gene, symbol, gene_type, outcome, method,
                               OR_MR, p_MR, OR_lower, OR_upper, source)]
}

pval_df <- rbindlist(all_pvals, fill = TRUE)
pval_df <- pval_df[!is.na(p_MR)]
pval_df[, p_unified_fdr := p.adjust(p_MR, method = "BH")]

log_info("Total tests: {nrow(pval_df)}")
log_info("  OXPHOS: {sum(pval_df$gene_type == 'OXPHOS')}")
log_info("  complement: {sum(pval_df$gene_type == 'complement_positive_control')}")

fwrite(pval_df, file.path(DIR_FDR, "all_pvalues_unified_fdr.csv"))
fwrite(pval_df[gene_type == "OXPHOS"], file.path(DIR_FDR, "OXPHOS_unified_fdr.csv"))
fwrite(pval_df[gene_type == "complement_positive_control"],
       file.path(DIR_FDR, "complement_unified_fdr.csv"))

gene_summary <- pval_df[, .(
  min_p = min(p_MR, na.rm = TRUE),
  min_fdr = min(p_unified_fdr, na.rm = TRUE),
  n_tests = .N
), by = .(gene, symbol, gene_type, outcome)][order(min_p)]
fwrite(gene_summary, file.path(DIR_FDR, "gene_unified_fdr_summary.csv"))

log_info("Top 20 by min p:")
print(head(gene_summary, 20))

comp <- pval_df[gene_type == "complement_positive_control"][order(p_MR)]
if (nrow(comp) == 0) {
  log_warn("No complement genes found")
} else {
  comp_min <- comp[, .SD[which.min(p_MR)], by = symbol][order(p_MR)]
  fwrite(comp_min, file.path(DIR_POS, "positive_control_summary.csv"))
  for (i in seq_len(nrow(comp_min))) {
    r <- comp_min[i]
    log_info("  {r$symbol} ({r$outcome}): OR = {round(r$OR_MR,3)}, p = {format(r$p_MR, scientific=TRUE)}, unified FDR = {format(r$p_unified_fdr, scientific=TRUE)}")
  }
  n_sig <- sum(comp_min$p_unified_fdr < 0.05)
  log_info("Positive controls FDR<0.05: {n_sig} / {nrow(comp_min)}")
}

oxphos <- gene_summary[gene_type == "OXPHOS"]
report <- data.table(
  metric = c("n_OXPHOS_genes_tested","n_nominal_p_lt_0.05",
             "n_unified_FDR_lt_0.05","min_p","min_FDR",
             "top_gene","top_outcome"),
  value = c(nrow(oxphos), sum(oxphos$min_p < 0.05, na.rm = TRUE),
            sum(oxphos$min_fdr < 0.05, na.rm = TRUE),
            min(oxphos$min_p, na.rm = TRUE), min(oxphos$min_fdr, na.rm = TRUE),
            oxphos[which.min(min_p), symbol], oxphos[which.min(min_p), outcome]))
fwrite(report, file.path(DIR_THIS, "summary_report.csv"))
print(report)

log_info("=== 08_4 done ===")
