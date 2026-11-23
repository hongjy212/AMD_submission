# 10_1D_supp_tables.R (v3.10)
# Supplementary tables S1-S16

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressWarnings(suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(openxlsx)
}))

results_dir <- DIR_RESULT
table_dir   <- file.path(results_dir, "10_1D_supp_tables")
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)
log_info("=== 10_1D Supplementary tables (v3.10, 16 workbooks S1-S16) ===")

read_or_null <- function(path) {
  if (file.exists(path)) read_csv(path, show_col_types = FALSE, progress = FALSE)
  else { log_warn(sprintf("missing file: %s", path)); NULL }
}

add_sheet_safe <- function(wb, sheet_name, df) {
  if (is.null(df) || nrow(df) == 0) {
    log_warn(sprintf("  skip sheet: %s (NULL or empty)", sheet_name))
    return(invisible(FALSE))
  }
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, df)
  log_info(sprintf("  + sheet: %s (%d x %d)", sheet_name, nrow(df), ncol(df)))
  invisible(TRUE)
}

save_wb <- function(wb, name) {
  sheets <- names(wb)
  if (length(sheets) == 0) {
    log_warn(sprintf("workbook %s has no sheets; skipped", name))
    return(invisible(FALSE))
  }
  f <- file.path(table_dir, name)
  saveWorkbook(wb, f, overwrite = TRUE)
  log_info(sprintf("saved workbook: %s (%d sheets: %s)",
                   name, length(sheets), paste(sheets, collapse = ", ")))
  invisible(TRUE)
}

# S1
log_info("--- S1. Cohort_and_data_summary ---")
wb <- createWorkbook()
add_sheet_safe(wb, "Cohorts", data.frame(
  cohort = c("GSE115828", "GSE221042", "GSE231536", "GSE29801",
             "GSE135092", "GSE265801", "GSE99287"),
  platform = c("Illumina RNA-seq", "10x 3' scRNA-seq", "Illumina methylation",
               "Affymetrix microarray", "Illumina microarray",
               "snRNA-seq", "ATAC-seq"),
  n_samples = c(453, 70378, "various", 62, 128, "reference", "reference"),
  role = c("Bulk discovery", "Single-cell", "Methylation + eQTM",
           "External replication", "ML external validation",
           "Deconvolution reference", "Promoter accessibility"),
  stringsAsFactors = FALSE
))
add_sheet_safe(wb, "ATAC_summary",
  read_or_null(file.path(results_dir,
    "07_1_atac_validation/tables/gene_set_ATAC_summary.csv")))
save_wb(wb, "S1_Cohort_and_data_summary.xlsx")

# S2
log_info("--- S2. Poly3F_model_diagnostics ---")
wb <- createWorkbook()
add_sheet_safe(wb, "Nested_model_comparison",
  read_or_null(file.path(results_dir,
    "02_12_poly3F_model_comparison/tables/poly3F_model_comparison.csv")))
add_sheet_safe(wb, "Model_selection_summary",
  read_or_null(file.path(results_dir,
    "02_12_poly3F_model_comparison/tables/poly3F_model_selection_summary.csv")))
add_sheet_safe(wb, "Local_FDR_qvalue",
  read_or_null(file.path(results_dir,
    "02_12_poly3F_model_comparison/tables/poly3F_all_genes_qvalue.csv")))
add_sheet_safe(wb, "Sig16_qvalue",
  read_or_null(file.path(results_dir,
    "02_12_poly3F_model_comparison/tables/poly3F_sig16_qvalue.csv")))
add_sheet_safe(wb, "Spline_NS",
  read_or_null(file.path(results_dir,
    "02_8_G2G3_changepoint_robustness/tables/spline_ns_results.csv")))
add_sheet_safe(wb, "Bootstrap_p3",
  read_or_null(file.path(results_dir,
    "02_8_G2G3_changepoint_robustness/tables/bootstrap_p3_16genes.csv")))
save_wb(wb, "S2_Poly3F_model_diagnostics.xlsx")

# S3
log_info("--- S3. G2G3_robustness ---")
wb <- createWorkbook()
add_sheet_safe(wb, "Leave_one_out",
  read_or_null(file.path(results_dir,
    "02_8_G2G3_changepoint_robustness/tables/FCF1_leave_one_out.csv")))
add_sheet_safe(wb, "No_FCF1",
  read_or_null(file.path(results_dir,
    "02_8_G2G3_changepoint_robustness/tables/poly3F_without_FCF1.csv")))
add_sheet_safe(wb, "Panel_direction",
  read_or_null(file.path(results_dir,
    "02_8_G2G3_changepoint_robustness/tables/panel_p3_direction.csv")))
add_sheet_safe(wb, "Panel_sign_test",
  read_or_null(file.path(results_dir,
    "02_8_G2G3_changepoint_robustness/tables/panel_p3_sign_test.csv")))
add_sheet_safe(wb, "Effect_size",
  read_or_null(file.path(results_dir,
    "02_8_G2G3_changepoint_robustness/tables/effect_size_results.csv")))
add_sheet_safe(wb, "Segmented_breakpoint",
  read_or_null(file.path(results_dir,
    "02_8_G2G3_changepoint_robustness/tables/segmented_breakpoint_results.csv")))
add_sheet_safe(wb, "Bayesian_changepoint",
  read_or_null(file.path(results_dir,
    "02_8_G2G3_changepoint_robustness/tables/bayesian_changepoint_mcp.csv")))
save_wb(wb, "S3_G2G3_robustness.xlsx")

# S4
log_info("--- S4. WGCNA_modules ---")
wb <- createWorkbook()
add_sheet_safe(wb, "Module_trait_full",
  read_or_null(file.path(results_dir,
    "05_2_module_trait/tables/module_summary_annotated.csv")))
add_sheet_safe(wb, "Module_DE",
  read_or_null(file.path(results_dir,
    "05_3_module_de/tables/module_summary_all.csv")))
add_sheet_safe(wb, "Cross_layer_complex",
  read_or_null(file.path(results_dir,
    "07_5c_OXPHOS_analysis/tables/OXPHOS13_complex_score.csv")))

mod <- read_or_null(file.path(results_dir,
  "05_2_module_trait/tables/module_summary_annotated.csv"))
mtc <- read_or_null(file.path(results_dir,
  "05_2_module_trait/tables/module_trait_correlation.csv"))
if (!is.null(mod)) {
  if (!is.null(mtc) && all(c("Module", "Trait", "r", "FDR") %in% colnames(mtc))) {
    fdr_wide <- mtc %>%
      select(Module, Trait, r, FDR) %>%
      pivot_wider(
        names_from  = Trait,
        values_from = c(r, FDR),
        names_glue  = "{.value}_{Trait}"
      )
    fdr_cols_new <- setdiff(
      grep("^FDR_", names(fdr_wide), value = TRUE),
      names(mod)
    )
    if (length(fdr_cols_new) > 0) {
      mod <- mod %>%
        left_join(fdr_wide %>% select(Module, all_of(fdr_cols_new)),
                  by = "Module")
    }
  }
  desired_cols <- c(
    "Module", "N_genes",
    "r_grade", "FDR_grade",
    "r_hUSI", "FDR_hUSI",
    "r_ML_Aging", "FDR_ML_Aging",
    "r_mito", "FDR_mito",
    "r_genes16", "FDR_genes16",
    "r_batch", "r_layer", "r_rin",
    "top_enrich", "top_enrich_FDR"
  )
  available_cols <- intersect(desired_cols, colnames(mod))
  part_B <- mod %>% select(all_of(available_cols))
  if ("FDR_grade" %in% colnames(part_B)) {
    part_B <- part_B %>% arrange(FDR_grade)
  }
  add_sheet_safe(wb, "Table2_PartB_WGCNA_module", part_B)
  log_info(sprintf("  S4 added Table2_PartB_WGCNA_module (%d x %d)",
                   nrow(part_B), ncol(part_B)))
}
save_wb(wb, "S4_WGCNA_modules.xlsx")

# S5
log_info("--- S5. eQTM_eGene ---")
wb <- createWorkbook()
add_sheet_safe(wb, "Fisher_enrichment",
  read_or_null(file.path(results_dir,
    "07_5c_OXPHOS_analysis/tables/fisher_eqtm_enrichment.csv")))
add_sheet_safe(wb, "Background_sensitivity",
  read_or_null(file.path(results_dir,
    "02_9_eQTM_background_sensitivity/tables/eqtm_background_sensitivity.csv")))
add_sheet_safe(wb, "Contingency_tables",
  read_or_null(file.path(results_dir,
    "02_9_eQTM_background_sensitivity/tables/eqtm_contingency_tables.csv")))
add_sheet_safe(wb, "Direction_by_module",
  read_or_null(file.path(results_dir,
    "02_9_eQTM_background_sensitivity/tables/direction_consistency_by_module.csv")))
add_sheet_safe(wb, "Independent_CpG",
  read_or_null(file.path(results_dir,
    "02_9_eQTM_background_sensitivity/tables/direction_consistency_independent_cpg.csv")))
add_sheet_safe(wb, "Block_permutation",
  read_or_null(file.path(results_dir,
    "02_11_eQTM_direction_permutation/tables/block_permutation_summary.csv")))
add_sheet_safe(wb, "CpG_poly2",
  read_or_null(file.path(results_dir,
    "07_5_GSE231536_methylation/tables/cpg_poly2_results.csv")))
add_sheet_safe(wb, "Threeway_greenyellow",
  read_or_null(file.path(results_dir,
    "02_10_OXPHOS_variance_decomposition/tables/threeway_consistency_greenyellow.csv")))

eqtm  <- read_or_null(file.path(results_dir,
  "07_5c_OXPHOS_analysis/tables/fisher_eqtm_enrichment.csv"))
cross <- read_or_null(file.path(results_dir,
  "07_5b_GSE231536_attachments/tables/cross_A_eqtm_egene.csv"))
if (!is.null(eqtm) && !is.null(cross)) {
  if (all(c("target", "n_total", "n_hit", "pct_hit", "n_ref",
            "fold", "OR", "p_fisher", "FDR") %in% colnames(eqtm)) &&
      all(c("target", "reference", "fold") %in% colnames(cross))) {

    cross_wide <- cross %>%
      select(target, reference, fold) %>%
      pivot_wider(names_from = reference, values_from = fold)

    eqtm_col  <- intersect(c("eQTM_targets", "eQTM", "eqtm"), names(cross_wide))[1]
    egene_col <- intersect(c("eGenes", "eGene", "egene"), names(cross_wide))[1]

    part_C <- eqtm %>%
      select(target, n_total, n_hit, pct_hit, n_ref,
             fold, OR, p_fisher, FDR) %>%
      left_join(cross_wide, by = "target")

    part_C <- part_C %>%
      mutate(
        eqtm_fold  = if (!is.na(eqtm_col))  .data[[eqtm_col]]  else NA_real_,
        egene_fold = if (!is.na(egene_col)) .data[[egene_col]] else NA_real_,
        axis = case_when(
          !is.na(eqtm_fold)  & FDR < 0.05 & eqtm_fold  > 1   ~ "epigenetic (eQTM)",
          !is.na(egene_fold) & egene_fold >= 1.5 &
            (is.na(FDR) | FDR >= 0.05)                      ~ "genetic (eGene)",
          TRUE                                              ~ "other"
        )
      ) %>%
      arrange(FDR)

    add_sheet_safe(wb, "Table2_PartC_eQTM_eGene", part_C)
    log_info(sprintf("  S5 added Table2_PartC_eQTM_eGene (%d x %d)",
                     nrow(part_C), ncol(part_C)))
  }
}
save_wb(wb, "S5_eQTM_eGene.xlsx")

# S6
log_info("--- S6. Metabolic_and_mito_pathways ---")
wb <- createWorkbook()
add_sheet_safe(wb, "Metabolic_pathways",
  read_or_null(file.path(results_dir,
    "07_5d_metabolism_pathways/tables/metabolic_pathway_vs_grade.csv")))
add_sheet_safe(wb, "Pathway_bootstrap",
  read_or_null(file.path(results_dir,
    "04_2_pathway_score/tables/bootstrap_p3_results.csv")))
add_sheet_safe(wb, "OXPHOS_variance_partition",
  read_or_null(file.path(results_dir,
    "02_10_OXPHOS_variance_decomposition/tables/variance_partition_all_sets.csv")))
add_sheet_safe(wb, "OXPHOS_by_complex",
  read_or_null(file.path(results_dir,
    "02_10_OXPHOS_variance_decomposition/tables/oxphos_by_complex.csv")))
add_sheet_safe(wb, "Gene_set_scores",
  read_or_null(file.path(results_dir,
    "02_10_OXPHOS_variance_decomposition/tables/gene_set_scores.csv")))
add_sheet_safe(wb, "OXPHOS_composite_perSD",
  read_or_null(file.path(results_dir,
    "07_5c_OXPHOS_analysis/tables/OXPHOS13_composite_perSD.csv")))
add_sheet_safe(wb, "OXPHOS_pathway_perSD",
  read_or_null(file.path(results_dir,
    "07_5c_OXPHOS_analysis/tables/OXPHOS_pathway_perSD.csv")))
add_sheet_safe(wb, "TRADE_impact_3modes",
  read_or_null(file.path(results_dir,
    "06_1_TRADE_vs_adjacent/tables/TRADE_impact_3modes.csv")))
add_sheet_safe(wb, "TRADE_impact_decomposition",
  read_or_null(file.path(results_dir,
    "06_1_TRADE_vs_adjacent/tables/TRADE_impact_decomposition.csv")))
save_wb(wb, "S6_Metabolic_and_mito_pathways.xlsx")

# S7
log_info("--- S7. Batch_correction_sensitivity ---")
wb <- createWorkbook()
add_sheet_safe(wb, "Batch_methods_summary",
  read_or_null(file.path(results_dir,
    "02_13_batch_correction_sensitivity/tables/batch_correction_sensitivity.csv")))
add_sheet_safe(wb, "16gene_FDR_by_method",
  read_or_null(file.path(results_dir,
    "02_13_batch_correction_sensitivity/tables/batch_correction_16gene_fdr.csv")))
add_sheet_safe(wb, "RIN_sensitivity",
  read_or_null(file.path(results_dir,
    "diag/diag_rin_sensitivity/tables/rin_robust_summary.csv")))
save_wb(wb, "S7_Batch_correction_sensitivity.xlsx")

# S8
log_info("--- S8. Causal_structure ---")
wb <- createWorkbook()
add_sheet_safe(wb, "BN_stable_edges",
  read_or_null(file.path(results_dir,
    "07_6a_bayesian_network/tables/bn_stable_edges.csv")))
add_sheet_safe(wb, "BN_edges_into_grade",
  read_or_null(file.path(results_dir,
    "07_6a_bayesian_network/tables/bn_edges_into_grade.csv")))
add_sheet_safe(wb, "SEM_model_comparison",
  read_or_null(file.path(results_dir,
    "07_6b_lavaan_mediation/tables/model_comparison.csv")))
save_wb(wb, "S8_Causal_structure.xlsx")

# S9
log_info("--- S9. MR_NDUFB3_main ---")
wb <- createWorkbook()
add_sheet_safe(wb, "Main_MR",
  read_or_null(file.path(results_dir, "07_7b_MR_plink/MR_earlyAMD.csv")))
add_sheet_safe(wb, "Cross_cohort_MR",
  read_or_null(file.path(results_dir,
    "07_7b_MR_plink/MR_cross_cohort_comparison.csv")))
add_sheet_safe(wb, "Sensitivity_grid",
  read_or_null(file.path(results_dir,
    "07_7b_MR_plink/NDUFB3_sensitivity_grid.csv")))
add_sheet_safe(wb, "Adipose_replication",
  read_or_null(file.path(results_dir,
    "07_7c_GTEx_MR_NDUFB3/tables/MR_results.csv")))
add_sheet_safe(wb, "Multi_tissue_eQTL",
  read_or_null(file.path(results_dir,
    "07_7a_NDUFB3/NDUFB3_multi_tissue_eQTL.csv")))
save_wb(wb, "S9_MR_NDUFB3_main.xlsx")

# S10
log_info("--- S10. MR_NDUFB3_sensitivity ---")
wb <- createWorkbook()
add_sheet_safe(wb, "Robust_MR",
  read_or_null(file.path(results_dir,
    "07_7d_MR_robust_sensitivity/tables/mr_robust_earlyAMD.csv")))
add_sheet_safe(wb, "Heterogeneity",
  read_or_null(file.path(results_dir,
    "07_7d_MR_robust_sensitivity/tables/heterogeneity_earlyAMD.csv")))
add_sheet_safe(wb, "Pleiotropy",
  read_or_null(file.path(results_dir,
    "07_7d_MR_robust_sensitivity/tables/pleiotropy_earlyAMD.csv")))
add_sheet_safe(wb, "Leave_one_out",
  read_or_null(file.path(results_dir,
    "07_7d_MR_robust_sensitivity/tables/leaveoneout_earlyAMD.csv")))
add_sheet_safe(wb, "Single_SNP",
  read_or_null(file.path(results_dir,
    "07_7d_MR_robust_sensitivity/tables/single_snp_earlyAMD.csv")))
add_sheet_safe(wb, "Unified_FDR_summary",
  read_or_null(file.path(results_dir,
    "07_7e_MR_unified_fdr/tables/NDUFB3_unified_fdr_summary.csv")))
add_sheet_safe(wb, "Unified_FDR_gene",
  read_or_null(file.path(results_dir,
    "07_7e_MR_unified_fdr/tables/gene_unified_fdr_summary.csv")))
add_sheet_safe(wb, "Unified_FDR_complement",
  read_or_null(file.path(results_dir,
    "07_7e_MR_unified_fdr/tables/complement_positive_control_unified_fdr.csv")))
add_sheet_safe(wb, "IV_strength",
  read_or_null(file.path(results_dir,
    "07_7d_MR_robust_sensitivity/tables/NDUFB3_IV_F_stat.csv")))
add_sheet_safe(wb, "Coloc_adipose_blood",
  read_or_null(file.path(results_dir,
    "07_7f_coloc_NDUFB3/coloc_NDUFB3_both_tissues.csv")))
add_sheet_safe(wb, "Retinal_eQTL_screening",
  read_or_null(file.path(results_dir,
    "07_7g_retinal_eqtl_screening_NDUFB3/tables/NDUFB3_retinal_eQTL_screening_summary.csv")))
save_wb(wb, "S10_MR_NDUFB3_sensitivity.xlsx")

# S11
log_info("--- S11. Single_cell ---")
wb <- createWorkbook()
add_sheet_safe(wb, "Cell_proportions",
  read_or_null(file.path(results_dir,
    "06_1_singlecell/tables/celltype_by_condition_ttest.csv")))
add_sheet_safe(wb, "Pathway_by_celltype",
  read_or_null(file.path(results_dir,
    "06_1_singlecell/tables/celltype_pathway_condition_test.csv")))
add_sheet_safe(wb, "RGC_complex_summary",
  read_or_null(file.path(results_dir,
    "07_3_NDUFB8_coexpression/tables/RGC_complex_100iter_summary.csv")))
add_sheet_safe(wb, "CellChat_matrix",
  read_or_null(file.path(results_dir,
    "06_3_cellchat/tables/top_axis.csv")))
add_sheet_safe(wb, "CellChat_communications",
  read_or_null(file.path(results_dir,
    "06_3_cellchat/tables/communication_count.csv")))
add_sheet_safe(wb, "CellChat_RPE_related",
  read_or_null(file.path(results_dir,
    "06_3_cellchat/tables/RPE_related_communications.csv")))
add_sheet_safe(wb, "CellChat_RPE_vs_RodCone",
  read_or_null(file.path(results_dir,
    "06_3_cellchat/tables/RPE_vs_RodCone_final.csv")))
save_wb(wb, "S11_Single_cell.xlsx")

# S12
log_info("--- S12. Virtual_KO_and_external ---")
wb <- createWorkbook()

add_sheet_safe(wb, "NDUFB3_KO_GSEA_RGC",
  read_or_null(file.path(results_dir,
    "07_3_scTenifoldKnk_NDUFB3_RGC/tables/NDUFB3_GSEA_RGC.csv")))
add_sheet_safe(wb, "NDUFB3_KO_16gene_RGC",
  read_or_null(file.path(results_dir,
    "07_3_scTenifoldKnk_NDUFB3_RGC/tables/NDUFB3_vs_16genes_RGC.csv")))
add_sheet_safe(wb, "NDUFB3_KO_GSEA_RPE",
  read_or_null(file.path(results_dir,
    "07_3_scTenifoldKnk_NDUFB3_RPE/tables/NDUFB3_GSEA_RPE.csv")))
add_sheet_safe(wb, "NDUFB3_KO_16gene_RPE",
  read_or_null(file.path(results_dir,
    "07_3_scTenifoldKnk_NDUFB3_RPE/tables/NDUFB3_vs_16genes_RPE.csv")))
add_sheet_safe(wb, "NDUFB3_KO_ComplexI_RPE",
  read_or_null(file.path(results_dir,
    "07_3_scTenifoldKnk_NDUFB3_RPE/tables/NDUFB3_vs_ComplexI_RPE.csv")))
add_sheet_safe(wb, "NDUFB3_KO_top100_RPE",
  read_or_null(file.path(results_dir,
    "07_3_scTenifoldKnk_NDUFB3_RPE/tables/NDUFB3_top100_RPE.csv")))
add_sheet_safe(wb, "NDUFB3_KO_GSEA_ROD",
  read_or_null(file.path(results_dir,
    "07_3_scTenifoldKnk_NDUFB3_ROD/tables/NDUFB3_GSEA_Rod.csv")))
add_sheet_safe(wb, "NDUFB3_KO_16gene_ROD",
  read_or_null(file.path(results_dir,
    "07_3_scTenifoldKnk_NDUFB3_ROD/tables/NDUFB3_vs_16genes_Rod.csv")))
add_sheet_safe(wb, "GSE29801_permutation",
  read_or_null(file.path(results_dir,
    "07_2_gse29801_replication/tables/permutation_results.csv")))
add_sheet_safe(wb, "ML_AUC",
  read_or_null(file.path(results_dir,
    "07_8_ML_prediction_model/Table5c_external_AUC_all_cohorts.csv")))
add_sheet_safe(wb, "LASSO_genes",
  read_or_null(file.path(results_dir,
    "07_8_ML_prediction_model/Table8_LASSO_annotated.csv")))

required_s12 <- c("NDUFB3_KO_GSEA_RGC", "NDUFB3_KO_16gene_RGC",
                  "NDUFB3_KO_GSEA_RPE", "NDUFB3_KO_16gene_RPE",
                  "NDUFB3_KO_ComplexI_RPE", "NDUFB3_KO_top100_RPE",
                  "NDUFB3_KO_GSEA_ROD", "NDUFB3_KO_16gene_ROD")
missing_s12 <- setdiff(required_s12, names(wb))
if (length(missing_s12) > 0) {
  log_warn(sprintf("S12 missing required sheets: %s",
                   paste(missing_s12, collapse = ", ")))
} else {
  log_info("S12: all 8 virtual KO sheets present (RGC / RPE / ROD)")
}
save_wb(wb, "S12_Virtual_KO_and_external.xlsx")

# S13
log_info("--- S13. CellChat RPE-related communications ---")
wb <- createWorkbook()
add_sheet_safe(wb, "communication_count",
  read_or_null(file.path(results_dir,
    "06_3_cellchat/tables/communication_count.csv")))
add_sheet_safe(wb, "pathway_prob",
  read_or_null(file.path(results_dir,
    "06_3_cellchat/tables/pathway_prob.csv")))
add_sheet_safe(wb, "RPE_related_communications",
  read_or_null(file.path(results_dir,
    "06_3_cellchat/tables/RPE_related_communications.csv")))
add_sheet_safe(wb, "RPE_pathway_count",
  read_or_null(file.path(results_dir,
    "06_3_cellchat/tables/RPE_pathway_count.csv")))
add_sheet_safe(wb, "RPE_vs_RodCone_observed",
  read_or_null(file.path(results_dir,
    "06_3_cellchat/tables/RPE_vs_RodCone_observed.csv")))
add_sheet_safe(wb, "RPE_vs_RodCone_permutation",
  read_or_null(file.path(results_dir,
    "06_3_cellchat/tables/RPE_vs_RodCone_permutation.csv")))
add_sheet_safe(wb, "RPE_vs_RodCone_bootstrap_CI",
  read_or_null(file.path(results_dir,
    "06_3_cellchat/tables/RPE_vs_RodCone_bootstrap_CI.csv")))
add_sheet_safe(wb, "RPE_vs_RodCone_final",
  read_or_null(file.path(results_dir,
    "06_3_cellchat/tables/RPE_vs_RodCone_final.csv")))
add_sheet_safe(wb, "RPE_marker_QC",
  read_or_null(file.path(results_dir,
    "06_3_cellchat/tables/RPE_marker_QC.csv")))
save_wb(wb, "S13_CellChat_RPE_related_communications.xlsx")

# S14
log_info("--- S14. Full_OXPHOS_MR_coloc ---")
wb <- createWorkbook()

add_sheet_safe(wb, "OXPHOS_gene_list",
  read_or_null(file.path(results_dir,
    "08_1_oxphos_gene_list/tables/oxphos_target_genes.csv")))
add_sheet_safe(wb, "OXPHOS_gene_summary",
  read_or_null(file.path(results_dir,
    "08_1_oxphos_gene_list/tables/oxphos_gene_summary.csv")))
add_sheet_safe(wb, "OXPHOS_gene_by_complex",
  read_or_null(file.path(results_dir,
    "08_1_oxphos_gene_list/tables/oxphos_gene_by_complex.csv")))
add_sheet_safe(wb, "IV_summary",
  read_or_null(file.path(results_dir,
    "08_2_eqtl_iv/02_iv/iv_summary.csv")))
add_sheet_safe(wb, "Per_gene_SNP_count",
  read_or_null(file.path(results_dir,
    "08_2_eqtl_iv/01_eqtl/per_gene_snp_count.csv")))
add_sheet_safe(wb, "MR_allstage",
  read_or_null(file.path(results_dir,
    "08_3_mr/01_mr/MR_allstage.csv")))
add_sheet_safe(wb, "MR_early",
  read_or_null(file.path(results_dir,
    "08_3_mr/01_mr/MR_early.csv")))
add_sheet_safe(wb, "Robust_MR",
  read_or_null(file.path(results_dir,
    "08_3_mr/02_robust/robust_methods.csv")))
add_sheet_safe(wb, "Heterogeneity",
  read_or_null(file.path(results_dir,
    "08_3_mr/02_robust/heterogeneity.csv")))
add_sheet_safe(wb, "Pleiotropy",
  read_or_null(file.path(results_dir,
    "08_3_mr/02_robust/pleiotropy.csv")))
add_sheet_safe(wb, "Leave_one_out",
  read_or_null(file.path(results_dir,
    "08_3_mr/02_robust/leaveoneout.csv")))
add_sheet_safe(wb, "Single_SNP",
  read_or_null(file.path(results_dir,
    "08_3_mr/02_robust/single_snp.csv")))
add_sheet_safe(wb, "Unified_FDR_all",
  read_or_null(file.path(results_dir,
    "08_4_fdr_positive_control/01_fdr/all_pvalues_unified_fdr.csv")))
add_sheet_safe(wb, "Unified_FDR_OXPHOS",
  read_or_null(file.path(results_dir,
    "08_4_fdr_positive_control/01_fdr/OXPHOS_unified_fdr.csv")))
add_sheet_safe(wb, "Unified_FDR_complement",
  read_or_null(file.path(results_dir,
    "08_4_fdr_positive_control/01_fdr/complement_unified_fdr.csv")))
add_sheet_safe(wb, "Unified_FDR_gene_summary",
  read_or_null(file.path(results_dir,
    "08_4_fdr_positive_control/01_fdr/gene_unified_fdr_summary.csv")))
add_sheet_safe(wb, "Positive_control_summary",
  read_or_null(file.path(results_dir,
    "08_4_fdr_positive_control/02_positive_control/positive_control_summary.csv")))
add_sheet_safe(wb, "Coloc_full_OXPHOS",
  read_or_null(file.path(results_dir,
    "08_5_coloc/coloc_all_genes.csv")))
save_wb(wb, "S14_Full_OXPHOS_MR_coloc.xlsx")

# S15
log_info("--- S15. Greenyellow_methylation_motif_hit ---")
wb <- createWorkbook()

add_sheet_safe(wb, "CpG_island_distribution",
  read_or_null(file.path(results_dir,
    "07_10c_cpg_island_distribution/tables/cpg_island_distribution.csv")))
add_sheet_safe(wb, "Motif_hit_rate_descriptive",
  read_or_null(file.path(results_dir,
    "07_10c_cpg_island_distribution/tables/motif_hit_rate_descriptive.csv")))
add_sheet_safe(wb, "Motif_enrichment_significant",
  read_or_null(file.path(results_dir,
    "07_10b_motif_enrichment/tables/motif_enrichment.csv")))

save_wb(wb, "S15_Greenyellow_methylation_motif_hit.xlsx")

# S16
log_info("--- S16. Greenyellow_motif_enrichment_monaLisa ---")
wb <- createWorkbook()

add_sheet_safe(wb, "Primary_enrichment",
  read_or_null(file.path(results_dir,
    "07_1b_atac_motif_enrichment/tables/motif_enrichment_primary.csv")))
add_sheet_safe(wb, "Robust_summary",
  read_or_null(file.path(results_dir,
    "07_1b_atac_motif_enrichment/tables/motif_robust_summary.csv")))
add_sheet_safe(wb, "Sensitivity_100_backgrounds",
  read_or_null(file.path(results_dir,
    "07_1b_atac_motif_enrichment/tables/motif_sensitivity_100_random_backgrounds.csv")))

save_wb(wb, "S16_Greenyellow_motif_enrichment_monaLisa.xlsx")

log_info("=== 10_1D v3.10 done (16 workbooks S1-S16) ===")

