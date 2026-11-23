# ============================================================
# 08_2_eqtl_iv.r
# Output:  results/08_2_eqtl_iv/
# ============================================================

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(ieugwasr)
})

ROOT     <- "E:/retina_tools/aging_rgc_axon/AMD"
DIR_THIS <- file.path(ROOT, "results", "08_2_eqtl_iv")
DIR_EQTL <- file.path(DIR_THIS, "01_eqtl")
DIR_IV   <- file.path(DIR_THIS, "02_iv")
for (d in c(DIR_EQTL, DIR_IV)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

EQTLGEN_FILE <- "E:/retina_tools/pubdata/eQTL/eQTLGen/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"
PLINK_BIN    <- "D:/tools/plink_win64_20260923/plink.exe"
BFILE        <- "E:/retina_tools/pubdata/eQTL/LD_reference/1kg.v3/EUR"
P_THRESHOLD  <- 1e-5
CLUMP_R2     <- 0.001
CLUMP_KB     <- 10000

OXPHOS_FULL_RDS <- file.path(ROOT,
  "results/08_1_oxphos_gene_list/data/oxphos_target_genes.rds")

COMPLEMENT_GENES <- data.frame(
  symbol = c("CFH","CD59","CFHR5","C3","CFB"),
  ensembl = c("ENSG00000000971","ENSG00000085063","ENSG00000134365",
              "ENSG00000125730","ENSG00000243649"),
  gene_type = "complement_positive_control", stringsAsFactors = FALSE)

log_info("=== 08_2 eqtl + iv ===")

if (!file.exists(OXPHOS_FULL_RDS)) stop("Run 08_1 first: ", OXPHOS_FULL_RDS)
ox <- as.data.frame(readRDS(OXPHOS_FULL_RDS))
ox <- ox[!is.na(ox$ensembl), ]
oxphos <- data.frame(symbol = ox$gene, ensembl = ox$ensembl,
                     gene_type = "OXPHOS", stringsAsFactors = FALSE)
log_info("Loaded {nrow(oxphos)} OXPHOS genes")

target_genes <- rbind(oxphos, COMPLEMENT_GENES)
log_info("Target genes: {nrow(target_genes)}")

log_info("Reading eQTLGen (2-5 min)...")
eqtl_all <- fread(EQTLGEN_FILE, sep = "\t", header = TRUE,
                  stringsAsFactors = FALSE, showProgress = TRUE)
log_info("Total rows: {nrow(eqtl_all)}")

eqtl <- eqtl_all[Gene %in% target_genes$ensembl, ]
rm(eqtl_all); gc()
log_info("After gene filter: {nrow(eqtl)}")

eqtl <- merge(eqtl, target_genes,
              by.x = "Gene", by.y = "ensembl", all.x = TRUE)

eqtl[, SNPPos := as.numeric(SNPPos)]
eqtl[, GenePos := as.numeric(GenePos)]
eqtl[, cis_dist := abs(SNPPos - GenePos)]
eqtl <- eqtl[cis_dist <= 1e6, ]
log_info("After cis filter: {nrow(eqtl)}")

eqtl[, Pvalue := as.numeric(Pvalue)]
eqtl_sig <- eqtl[Pvalue < P_THRESHOLD, ]
log_info("After P<{P_THRESHOLD}: {nrow(eqtl_sig)}")

fwrite(eqtl_sig, file.path(DIR_EQTL, "eqtl_cis.csv"))
per_gene <- eqtl_sig[, .N, by = .(symbol, gene_type)][order(-N)]
fwrite(per_gene, file.path(DIR_EQTL, "per_gene_snp_count.csv"))

log_info("--- LD clumping ---")
iv_list <- list()
for (g in unique(eqtl_sig$Gene)) {
  sub <- eqtl_sig[Gene == g, ]
  g_sym <- unique(sub$symbol)[1]
  if (nrow(sub) < 2) {
    iv_list[[g]] <- sub
    log_info("  [{g_sym}] 1 SNP (no clump)")
    next
  }
  clumped <- tryCatch({
    ieugwasr::ld_clump(
      dat = dplyr::tibble(rsid = sub$SNP, pval = sub$Pvalue),
      plink_bin = PLINK_BIN, bfile = BFILE,
      clump_kb = CLUMP_KB, clump_r2 = CLUMP_R2, clump_p = P_THRESHOLD)
  }, error = function(e) {
    log_warn("  [{g_sym}] clump fail: {conditionMessage(e)}"); NULL
  })
  if (is.null(clumped) || nrow(clumped) == 0) {
    iv_list[[g]] <- sub[0, ]
    log_warn("  [{g_sym}] 0 SNP after clump")
  } else {
    iv_list[[g]] <- sub[SNP %in% clumped$rsid, ]
    log_info("  [{g_sym}] {nrow(sub)} -> {nrow(iv_list[[g]])} SNP")
  }
}

iv <- rbindlist(iv_list)
if (nrow(iv) == 0) stop("No IVs survived clumping")

iv[, Zscore := as.numeric(Zscore)]
iv[, NrSamples := as.numeric(NrSamples)]
iv[, beta := Zscore / sqrt(NrSamples)]
iv[, se   := 1 / sqrt(NrSamples)]
iv[, F_stat := Zscore^2]

fwrite(iv, file.path(DIR_IV, "iv_all.csv"))
iv_summary <- iv[, .(
  n_iv = .N, max_F = max(F_stat, na.rm = TRUE),
  min_p = min(Pvalue, na.rm = TRUE)
), by = .(symbol, gene_type, Gene)][order(gene_type, -n_iv)]
fwrite(iv_summary, file.path(DIR_IV, "iv_summary.csv"))

log_info("Total IVs: {nrow(iv)}")
log_info("Genes with IV: {nrow(iv_summary)}")
log_info("=== 08_2 done ===")
