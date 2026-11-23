# ============================================================
# 07_7f_coloc_NDUFB3.r
# Purpose: coloc.abf between GTEx v10 adipose/blood NDUFB3 cis-eQTL
#          and early AMD GWAS (ebi-a-GCST010723)
# ============================================================
setwd("E:/retina_tools/aging_rgc_axon/AMD/")
library(data.table)
gwas <- readRDS("E:/retina_tools/pubdata/eQTL/gwas/GCST010723_NDUFB3_region.rds")
fwrite(as.data.table(gwas), "E:/retina_tools/aging_rgc_axon/AMD/results/07_7f_coloc_NDUFB3/ndufb3_gwas_region.csv")


