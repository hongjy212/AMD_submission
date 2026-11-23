# ============================================================
# 07_7f2_coloc_NDUFB3.r
# Purpose: Assemble COLSTATS v2.3 colocalization results into a tidy
#          table for Table 3 and Figure 5D.
#
# Data provenance:
#   - Adipose_Subcutaneous (cs_4612_c, n = 581):
#       transcribed from results/07_7f_coloc_NDUFB3/NDUFB3_colstats.txt
#       (COLSTATS v2.3 raw web output, session 2026-10-11)
#   - Whole_Blood (cs_1067364_c, n = 670):
#       transcribed from COLSTATS v2.3 web output, session 2026-10-12;
#       web-session screenshots archived in
#       results/07_7f_coloc_NDUFB3/figure/ (not included in code package)
#
# Note: PP.H0-PP.H4, variant counts, top variants and r2 are transcribed
#       from the COLSTATS output, not computed locally. The complete
#       workflow is documented in Supplementary Code 1 and DATA_SOURCES.md.
# ============================================================

setwd("E:/retina_tools/aging_rgc_axon/AMD/")

library(data.table)

coloc_results <- data.table(
  analysis = c("NDUFB3_adipose_vs_earlyAMD", "NDUFB3_blood_vs_earlyAMD"),
  tissue = c("Adipose_Subcutaneous", "Whole_Blood"),
  n_samples = c(581, 670),
  region = c("2:201120000-203133000", "2:201120000-203133000"),
  n_variants_t1 = c(3429, 3429),
  n_variants_t2 = c(2773, 2874),
  n_shared = c(2495, 2538),
  PP_H0 = c(0.0168, 0.0384),
  PP_H1 = c(0.3644, 0.8339),
  PP_H2 = c(0.0042, 0.0048),
  PP_H3 = c(0.0911, 0.1043),
  PP_H4 = c(0.5235, 0.0185),
  top_variant_T1 = c("chr2_201135672_C_G", "chr2_201135672_C_G"),
  top_variant_T2 = c("chr2_201131954_CT_C", "chr2_201690409_T_C"),
  r2_top_variants = c(0.88, 0.00),
  direction_concordant = c(TRUE, TRUE)
)

fwrite(coloc_results, "E:/retina_tools/aging_rgc_axon/AMD/results/07_7f_coloc_NDUFB3/coloc_NDUFB3_both_tissues.csv")