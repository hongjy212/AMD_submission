# ============================================================
# 08_5_coloc.r
# Output:  results/08_5_coloc/
# ============================================================

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(coloc)
  library(rtracklayer); library(GenomicRanges)
})

ROOT     <- "E:/retina_tools/aging_rgc_axon/AMD"
DIR_THIS <- file.path(ROOT, "results", "08_5_coloc")
dir.create(DIR_THIS, recursive = TRUE, showWarnings = FALSE)

IN_IV    <- file.path(ROOT, "results/08_2_eqtl_iv/02_iv/iv_all.csv")
GWAS_FILE <- "E:/retina_tools/pubdata/eQTL/eQTLGen/GCST010723_buildGRCh37.tsv"
CHAIN_FILE <- "E:/retina_tools/pubdata/eQTL/eQTLGen/hg19ToHg38.over.chain/hg19ToHg38.over.chain"

GTEX_DIR <- "E:/retina_tools/pubdata/eQTL/eQTLGen/GTEx_Analysis_v8_eQTL"
GTEX_FILES <- c(
  Adipose_Subcutaneous = file.path(GTEX_DIR, "Adipose_Subcutaneous.v8.signif_variant_gene_pairs.txt.gz"),
  Whole_Blood          = file.path(GTEX_DIR, "Whole_Blood.v8.signif_variant_gene_pairs.txt.gz")
)
GTEX_N <- c(Adipose_Subcutaneous = 581L, Whole_Blood = 670L)

N_GWAS_CASE <- 14034L
N_GWAS_CTRL <- 91214L
N_GWAS      <- N_GWAS_CASE + N_GWAS_CTRL
S_GWAS      <- N_GWAS_CASE / N_GWAS

log_info("=== 08_5 coloc ===")

iv <- fread(IN_IV)
genes <- unique(iv[, .(ensembl = Gene, symbol)])
genes[, ensembl_nover := sub("\\..*$", "", ensembl)]
log_info("Genes to coloc: {nrow(genes)}")

log_info("Reading GWAS...")
gwas <- fread(GWAS_FILE,
              select = c("variant_id","chromosome","base_pair_location",
                         "effect_allele","other_allele",
                         "effect_allele_frequency","beta","standard_error","p_value"))
log_info("GWAS rows raw: {nrow(gwas)}")

gwas <- gwas[!is.na(variant_id)]
gwas <- gwas[!is.na(effect_allele) & !is.na(other_allele)]
log_info("GWAS after NA filter: {nrow(gwas)}")

# ---- 2. liftover ----
gr <- GRanges(
  seqnames = paste0("chr", gwas$chromosome),
  ranges   = IRanges(start = gwas$base_pair_location, width = 1)
)
mcols(gr)$idx <- seq_len(nrow(gwas))

ch <- import.chain(CHAIN_FILE)
gr38 <- liftOver(gr, ch)

n_map <- elementNROWS(gr38)
idx_keep <- which(n_map == 1)
gr38_1 <- unlist(gr38[idx_keep])

idx_orig <- mcols(gr38_1)$idx
log_info("Mapped SNPs: {length(idx_orig)} / {nrow(gwas)}")

gwas38 <- copy(gwas[idx_orig])
gwas38$chr38 <- as.integer(sub("^chr", "", as.character(seqnames(gr38_1))))
gwas38$pos38 <- start(gr38_1)
gwas38$MAF <- pmin(gwas38$effect_allele_frequency,
                   1 - gwas38$effect_allele_frequency)
gwas38 <- gwas38[!is.na(gwas38$MAF) & gwas38$MAF > 0 & gwas38$MAF < 1, ]
gwas38 <- gwas38[!duplicated(paste(gwas38$chr38, gwas38$pos38)), ]
log_info("GWAS after liftover + dedup: {nrow(gwas38)}")

rm(gr, gr38, gr38_1); gc()

coloc_results <- list()

for (tissue in names(GTEX_FILES)) {
  f <- GTEX_FILES[[tissue]]
  if (!file.exists(f)) { log_warn("GTEx file not found: {f}"); next }

  log_info("--- {tissue} ---")
  eqtl_all <- fread(f)
  log_info("  eQTL rows: {nrow(eqtl_all)}")
  eqtl_all[, ensembl_nover := sub("\\..*$", "", gene_id)]

  eqtl_all <- eqtl_all[!is.na(slope) & !is.na(slope_se) &
                       !is.na(maf) & !is.na(variant_id)]

  for (i in seq_len(nrow(genes))) {
    ens <- genes$ensembl_nover[i]
    sym <- genes$symbol[i]

    eqtl <- eqtl_all[ensembl_nover == ens]
    if (nrow(eqtl) < 2) next

    parts <- tstrsplit(eqtl$variant_id, "_", fixed = TRUE)
    eqtl[, chr := as.integer(sub("^chr", "", parts[[1]]))]
    eqtl[, pos := as.numeric(parts[[2]])]

    med_pos <- median(eqtl$pos, na.rm = TRUE)
    gwas_sub <- gwas38[gwas38$chr38 == eqtl$chr[1] &
                       gwas38$pos38 >= med_pos - 1e6 &
                       gwas38$pos38 <= med_pos + 1e6, ]
    if (nrow(gwas_sub) < 2) next

    common <- merge(
      gwas_sub[, .(chr38, pos38, snp = variant_id,
                   beta_out = beta, se_out = standard_error,
                   pval_out = p_value, MAF_out = MAF,
                   ea_out = effect_allele, oa_out = other_allele)],
      eqtl[, .(chr, pos, beta_eq = slope, se_eq = slope_se,
               pval_eq = pval_nominal, maf_eq = maf)],
      by.x = c("chr38", "pos38"),
      by.y = c("chr", "pos")
    )
    if (nrow(common) < 2) next

    common <- common[!is.na(snp) & !is.na(beta_out) & !is.na(se_out) &
                     !is.na(beta_eq) & !is.na(se_eq) &
                     !is.na(MAF_out) & !is.na(maf_eq), ]
    if (nrow(common) < 2) next

    D1 <- list(beta = common$beta_out,
               varbeta = common$se_out^2,
               snp = common$snp, position = common$pos38,
               N = N_GWAS, type = "cc", s = S_GWAS,
               MAF = common$MAF_out)
    D2 <- list(beta = common$beta_eq,
               varbeta = common$se_eq^2,
               snp = common$snp, position = common$pos38,
               N = GTEX_N[[tissue]], type = "quant",
               MAF = common$maf_eq)

    res <- tryCatch(coloc.abf(D1, D2)$summary,
                    error = function(e) {
                      log_warn(sprintf("[%s|%s] coloc failed: %s",
                                       sym, tissue, conditionMessage(e))); NULL
                    })
    if (is.null(res)) next

    coloc_results[[length(coloc_results) + 1]] <- data.table(
      gene = sym, tissue = tissue,
      nsnps = res["nsnps"],
      PP_H0 = res["PP.H0.abf"], PP_H1 = res["PP.H1.abf"],
      PP_H2 = res["PP.H2.abf"], PP_H3 = res["PP.H3.abf"],
      PP_H4 = res["PP.H4.abf"])
    log_info(sprintf("  [%s] PP.H4 = %.3f", sym, res["PP.H4.abf"]))
  }
}

if (length(coloc_results) > 0) {
  coloc_dt <- rbindlist(coloc_results)
  fwrite(coloc_dt, file.path(DIR_THIS, "coloc_all_genes.csv"))
  log_info("Saved: coloc_all_genes.csv ({nrow(coloc_dt)} rows)")
  print(coloc_dt)
} else {
  log_warn("No coloc results generated")
}
log_info("=== 08_5 done ===")

