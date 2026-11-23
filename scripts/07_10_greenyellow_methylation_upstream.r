#!/usr/bin/env Rscript
# ============================================================
# 07_10_greenyellow_methylation_upstream.r (v3, English)
# Purpose: map greenyellow eQTM CpGs to hg19 coordinates and
#          export BED files for downstream motif / CpG island analysis
#
# Input :
#   E:/retina_tools/raw/GSE231536/eQTM_significant_pairs_retina.txt
#   results/05_1_consensus_wgcna/tables/consensus_module_genes_annotated.csv
#   results/01_1_prepare_bulk_data/data/gene_annotation.rds
#   IlluminaHumanMethylationEPICanno.ilm10b4.hg19 annotation
#
# Output: results/07_10_greenyellow_methylation_upstream/
#   tables/greenyellow_eqtm_cpg_with_coords.csv
#   bed/target_neg_beta.bed
#   bed/target_pos_beta.bed
#   bed/target_all_greenyellow.bed
#   bed/target_promoter_only.bed
#   bed/background_all_eqtm.bed
# ============================================================

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(minfi)
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  library(dplyr); library(tidyr)
})

DIR_THIS <- setup_output_dir("07_10_greenyellow_methylation_upstream")
dir.create(file.path(DIR_THIS, "bed"), showWarnings = FALSE)
log_info("=== 07_10_greenyellow_methylation_upstream.r v3 start ===")

strip_v <- function(x) sub("\\..*$", "", x)

# ------------------------------------------------------------
# 1. Load greenyellow gene set
# ------------------------------------------------------------
log_info("--- 1. Load greenyellow genes ---")

module_df  <- read_upstream("05_1_consensus_wgcna", "tables",
                             "consensus_module_genes_annotated.csv")
anno_gene  <- read_upstream("01_1_prepare_bulk_data", "data",
                             "gene_annotation.rds")

gy_ens <- unique(strip_v(module_df$ENSEMBL[module_df$Module == "greenyellow"]))
gy_ens <- gy_ens[!is.na(gy_ens) & gy_ens != ""]
gy_sym <- unique(anno_gene$symbol[strip_v(anno_gene$ensembl) %in% gy_ens])
gy_sym <- gy_sym[!is.na(gy_sym) & gy_sym != ""]
log_info("greenyellow module: {length(gy_ens)} genes, {length(gy_sym)} symbols")

# ------------------------------------------------------------
# 2. Load eQTM table and subset to greenyellow
# ------------------------------------------------------------
log_info("--- 2. Load eQTM table ---")

eqtm <- read.delim("E:/retina_tools/raw/GSE231536/eQTM_significant_pairs_retina.txt",
                    stringsAsFactors = FALSE)
log_info("eQTM pairs in total: {nrow(eqtm)}")

gy_eqtm <- eqtm[eqtm$trait %in% gy_ens, ]
gy_eqtm <- gy_eqtm[!duplicated(gy_eqtm$geno), ]
log_info("greenyellow eQTM CpGs: {nrow(gy_eqtm)}")
log_info("p.adj < 0.05: {sum(gy_eqtm$p.adj < 0.05)}")

# ------------------------------------------------------------
# 3. Add hg19 coordinates from EPIC annotation
# ------------------------------------------------------------
log_info("--- 3. Add EPIC coordinates ---")

ann    <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
ann_df <- as.data.frame(ann)

idx <- match(gy_eqtm$geno, ann_df$Name)
gy_eqtm$chr             <- ann_df$chr[idx]
gy_eqtm$pos             <- ann_df$pos[idx]
gy_eqtm$strand          <- ann_df$strand[idx]
gy_eqtm$gene            <- ann_df$UCSC_RefGene_Name[idx]
gy_eqtm$region          <- ann_df$UCSC_RefGene_Group[idx]
gy_eqtm$relation_island <- ann_df$Relation_to_Island[idx]

n_matched <- sum(!is.na(gy_eqtm$chr))
log_info("coordinates matched: {n_matched} / {nrow(gy_eqtm)}")

gy_eqtm <- gy_eqtm[!is.na(gy_eqtm$chr) & !is.na(gy_eqtm$pos), ]
log_info("CpGs with valid coordinates: {nrow(gy_eqtm)}")

sym_map <- setNames(anno_gene$symbol, strip_v(anno_gene$ensembl))
gy_eqtm$gene_sym <- toupper(sym_map[gy_eqtm$trait])

save_table(gy_eqtm, DIR_THIS, "greenyellow_eqtm_cpg_with_coords.csv")

# ------------------------------------------------------------
# 4. Split by eQTM beta sign
# ------------------------------------------------------------
neg_cpg <- gy_eqtm[gy_eqtm$beta < 0, ]
pos_cpg <- gy_eqtm[gy_eqtm$beta > 0, ]
log_info("beta < 0 (negative): {nrow(neg_cpg)}")
log_info("beta > 0 (positive): {nrow(pos_cpg)}")

# ------------------------------------------------------------
# 5. Export BED files
# ------------------------------------------------------------
write_bed <- function(df, file, extra_cols = NULL) {
  if (nrow(df) == 0) {
    log_warn("empty table, skipping: {basename(file)}")
    return(invisible(NULL))
  }
  bed <- data.frame(
    chr    = paste0("chr", sub("^chr", "", df$chr)),
    start  = pmax(1, df$pos - 51),
    end    = df$pos + 50,
    name   = df$geno,
    score  = 0,
    strand = ifelse(df$strand == "+", "+", "-"),
    stringsAsFactors = FALSE)
  if (!is.null(extra_cols)) {
    for (col in extra_cols) bed[[col]] <- df[[col]]
  }
  write.table(bed, file, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  log_info("saved: {basename(file)} ({nrow(bed)} CpGs)")
}

write_bed(neg_cpg, file.path(DIR_THIS, "bed", "target_neg_beta.bed"),
          extra_cols = c("beta", "p.adj"))
write_bed(pos_cpg, file.path(DIR_THIS, "bed", "target_pos_beta.bed"),
          extra_cols = c("beta", "p.adj"))
write_bed(gy_eqtm, file.path(DIR_THIS, "bed", "target_all_greenyellow.bed"),
          extra_cols = c("beta", "p.adj"))

eqtm_all <- eqtm[!duplicated(eqtm$geno), ]
idx2 <- match(eqtm_all$geno, ann_df$Name)
eqtm_all$chr    <- ann_df$chr[idx2]
eqtm_all$pos    <- ann_df$pos[idx2]
eqtm_all$strand <- ann_df$strand[idx2]
eqtm_all <- eqtm_all[!is.na(eqtm_all$chr) & !is.na(eqtm_all$pos), ]
write_bed(eqtm_all, file.path(DIR_THIS, "bed", "background_all_eqtm.bed"))
log_info("background eQTM CpGs: {nrow(eqtm_all)}")

# ------------------------------------------------------------
# 6. Stratify by TSS distance
# ------------------------------------------------------------
log_info("--- 6. Stratify by UCSC_RefGene_Group ---")

region_tab <- as.data.frame(table(gy_eqtm$region))
print(region_tab)

promoter_cpg <- gy_eqtm[grepl("TSS200|TSS1500", gy_eqtm$region), ]
log_info("promoter CpGs (TSS200/TSS1500): {nrow(promoter_cpg)}")

if (nrow(promoter_cpg) >= 20) {
  write_bed(promoter_cpg,
            file.path(DIR_THIS, "bed", "target_promoter_only.bed"),
            extra_cols = c("beta", "p.adj"))
}

# ------------------------------------------------------------
# 7. Summary
# ------------------------------------------------------------
cat("\n", strrep("=", 70), "\n", sep = "")
cat("=== 07_10 v3 summary ===\n")
cat(strrep("=", 70), "\n", sep = "")

cat("\n[1] greenyellow genes: ", length(gy_sym), "\n", sep = "")
cat("[2] eQTM CpGs: ", nrow(gy_eqtm), " (", nrow(neg_cpg), " negative, ",
    nrow(pos_cpg), " positive)\n", sep = "")
cat("[3] promoter CpGs: ", nrow(promoter_cpg), "\n", sep = "")

cat("\n[4] CpG count by UCSC_RefGene_Group:\n")
print(region_tab, row.names = FALSE)

cat("\n[5] BED files:\n")
for (f in list.files(file.path(DIR_THIS, "bed"), full.names = TRUE)) {
  cat("  ", basename(f), " (", file.info(f)$size, " bytes)\n", sep = "")
}

cat("\n[6] Motif enrichment input files (for downstream analysis):\n")
cat("  target_neg_beta.bed      : beta < 0, methylation up -> expression down\n")
cat("  target_pos_beta.bed      : beta > 0, methylation up -> expression up\n")
cat("  target_promoter_only.bed : TSS200/TSS1500 CpGs\n")
cat("  target_all_greenyellow.bed : all greenyellow eQTM CpGs\n")
cat("  background_all_eqtm.bed  : background (all eQTM CpGs)\n")

log_info("=== 07_10 v3 done ===")

