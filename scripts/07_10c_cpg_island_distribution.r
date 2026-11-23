#!/usr/bin/env Rscript
# ============================================================
# 07_10c_cpg_island_distribution.r
# Purpose:
#   (A) CpG island distribution of greenyellow eQTM CpGs
#       vs all eQTM CpGs (background)
#   (B) Descriptive motif hit-rate comparison for OXPHOS-related
#       motifs (no p-value claim; motif enrichment itself was
#       negative in 07_10b and is reported as a negative result)
#
# Input:
#   results/07_10_greenyellow_methylation_upstream/tables/
#     greenyellow_eqtm_cpg_with_coords.csv
#   results/07_10_greenyellow_methylation_upstream/bed/
#     target_neg_beta.bed, target_pos_beta.bed,
#     background_all_eqtm.bed
#
# Output: results/07_10c_cpg_island_distribution/
#   tables/cpg_island_distribution.csv
#   tables/motif_hit_rate_descriptive.csv
# ============================================================

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(minfi)
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  library(dplyr); library(tidyr)
})

DIR_THIS <- setup_output_dir("07_10c_cpg_island_distribution")
BED_DIR  <- file.path(DIR_RESULT, "07_10_greenyellow_methylation_upstream", "bed")
TAB_DIR  <- file.path(DIR_RESULT, "07_10_greenyellow_methylation_upstream", "tables")
log_info("=== 07_10c_cpg_island_distribution.r start ===")

# ============================================================
# Part A: CpG island distribution
# ============================================================
log_info("--- Part A: CpG island distribution ---")

gy_coords <- read.csv(file.path(TAB_DIR, "greenyellow_eqtm_cpg_with_coords.csv"),
                       stringsAsFactors = FALSE)
log_info("greenyellow eQTM CpGs: {nrow(gy_coords)}")

ann    <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
ann_df <- as.data.frame(ann)

eqtm <- read.delim("E:/retina_tools/raw/GSE231536/eQTM_significant_pairs_retina.txt",
                    stringsAsFactors = FALSE)
eqtm_all <- eqtm[!duplicated(eqtm$geno), ]
idx_bg <- match(eqtm_all$geno, ann_df$Name)
eqtm_all$relation_island <- ann_df$Relation_to_Island[idx_bg]
eqtm_all <- eqtm_all[!is.na(eqtm_all$relation_island), ]

target_dist <- as.data.frame(table(gy_coords$relation_island))
colnames(target_dist) <- c("relation_island", "n_target")
target_dist$pct_target <- round(100 * target_dist$n_target / nrow(gy_coords), 1)

bg_dist <- as.data.frame(table(eqtm_all$relation_island))
colnames(bg_dist) <- c("relation_island", "n_bg")
bg_dist$pct_bg <- round(100 * bg_dist$n_bg / nrow(eqtm_all), 1)

island_cmp <- merge(target_dist, bg_dist, by = "relation_island", all = TRUE)
island_cmp[is.na(island_cmp)] <- 0
island_cmp$fold <- round(island_cmp$pct_target / pmax(island_cmp$pct_bg, 0.1), 2)
island_cmp <- island_cmp[order(-island_cmp$pct_target), ]

save_table(island_cmp, DIR_THIS, "cpg_island_distribution.csv")

cat("\n=== CpG island distribution ===\n")
print(island_cmp, row.names = FALSE)

n_island   <- sum(gy_coords$relation_island %in% c("Island", "N_Shore", "S_Shore"))
n_open_sea <- sum(gy_coords$relation_island == "OpenSea")
log_info("In CpG island/shore: {n_island} / {nrow(gy_coords)} ({round(100*n_island/nrow(gy_coords),1)}%)")
log_info("In OpenSea: {n_open_sea} / {nrow(gy_coords)} ({round(100*n_open_sea/nrow(gy_coords),1)}%)")

# ============================================================
# Part B: Descriptive motif hit-rate comparison
# ============================================================
log_info("--- Part B: Descriptive motif hit-rate comparison ---")

has_motif_pkgs <- requireNamespace("universalmotif", quietly = TRUE) &&
                  requireNamespace("JASPAR2024", quietly = TRUE) &&
                  requireNamespace("TFBSTools", quietly = TRUE) &&
                  requireNamespace("RSQLite", quietly = TRUE) &&
                  requireNamespace("BSgenome.Hsapiens.UCSC.hg19", quietly = TRUE) &&
                  requireNamespace("GenomicRanges", quietly = TRUE) &&
                  requireNamespace("Biostrings", quietly = TRUE)

if (!has_motif_pkgs) {
  log_warn("Motif packages not installed; skipping Part B")
} else {
  suppressPackageStartupMessages({
    library(universalmotif); library(JASPAR2024); library(TFBSTools)
    library(RSQLite); library(BSgenome.Hsapiens.UCSC.hg19)
    library(GenomicRanges); library(IRanges); library(Biostrings)
  })

  read_bed_gr <- function(f) {
    d <- read.delim(f, header = FALSE, stringsAsFactors = FALSE)
    GRanges(seqnames = d$V1,
            ranges   = IRanges(start = d$V2 + 1, end = d$V3),
            strand   = d$V6, name = d$V4)
  }

  bed_files <- c(
    neg_beta   = "target_neg_beta.bed",
    pos_beta   = "target_pos_beta.bed",
    background = "background_all_eqtm.bed"
  )
  gr_list <- list()
  for (nm in names(bed_files)) {
    f <- file.path(BED_DIR, bed_files[[nm]])
    if (file.exists(f)) gr_list[[nm]] <- read_bed_gr(f)
  }

  get_seq <- function(gr) getSeq(BSgenome.Hsapiens.UCSC.hg19, gr)
  clean_seq <- function(dss) dss[!grepl("N", as.character(dss), fixed = TRUE)]

  set.seed(SEED)
  seq_bg_full <- clean_seq(get_seq(gr_list$background))
  seq_bg <- seq_bg_full[sample(length(seq_bg_full),
                                min(2000, length(seq_bg_full)))]
  log_info("background sequences: {length(seq_bg)}")

  seq_neg <- clean_seq(get_seq(gr_list$neg_beta))
  seq_pos <- clean_seq(get_seq(gr_list$pos_beta))
  log_info("neg_beta: {length(seq_neg)}, pos_beta: {length(seq_pos)}")

  jaspar <- JASPAR2024::JASPAR2024()
  sq24   <- RSQLite::dbConnect(RSQLite::SQLite(), db(jaspar))
  motifs <- TFBSTools::getMatrixSet(sq24,
              list(species = "Homo sapiens", collection = "CORE"))
  RSQLite::dbDisconnect(sq24)
  motifs_um <- universalmotif::convert_motifs(motifs)
  motif_names <- sapply(motifs_um, function(m) m@name)

  oxphos_tfs <- c("NRF1","NRF2","NFE2L2","ESRRA","ESRRB","PPARGC1A",
                  "SP1","SP2","SP3","SP4","KLF4","KLF9","KLF15",
                  "E2F1","E2F2","E2F4","YY1","MYC","SREBF1","SREBF2",
                  "GABPA","TFAM","PPARG","NR1H3","NRF")
  motifs_ox <- motifs_um[grepl(paste(oxphos_tfs, collapse = "|"),
                                motif_names, ignore.case = TRUE)]
  log_info("OXPHOS-related motifs: {length(motifs_ox)}")

  count_hits <- function(motifs_sub, seq_target) {
    sapply(motifs_sub, function(m) {
      hits <- scan_sequences(m, seq_target, threshold = 0.001, RC = TRUE)
      length(unique(hits$sequence))
    })
  }

  hit_neg <- count_hits(motifs_ox, seq_neg)
  hit_pos <- count_hits(motifs_ox, seq_pos)
  hit_bg  <- count_hits(motifs_ox, seq_bg)

  hit_df <- data.frame(
    motif       = sapply(motifs_ox, function(m) m@name),
    n_neg       = as.integer(hit_neg),
    pct_neg     = round(100 * hit_neg / length(seq_neg), 1),
    n_pos       = as.integer(hit_pos),
    pct_pos     = round(100 * hit_pos / length(seq_pos), 1),
    n_bg        = as.integer(hit_bg),
    pct_bg      = round(100 * hit_bg / length(seq_bg), 1),
    stringsAsFactors = FALSE)
  hit_df$fold_neg <- round(hit_df$pct_neg / pmax(hit_df$pct_bg, 0.1), 2)
  hit_df <- hit_df[order(-hit_df$fold_neg), ]

  save_table(hit_df, DIR_THIS, "motif_hit_rate_descriptive.csv")

  cat("\n=== Descriptive motif hit rates (top 15 by fold_neg) ===\n")
  print(head(hit_df, 15), row.names = FALSE)
}

# ============================================================
# Summary
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("=== 07_10c summary ===\n")
cat(strrep("=", 70), "\n", sep = "")

cat("\n[A] CpG island distribution:\n")
print(island_cmp, row.names = FALSE)

cat("\n[B] Descriptive motif hit rates:\n")
if (exists("hit_df")) {
  cat("  Top motif by fold_neg:", hit_df$motif[1],
      "(fold =", hit_df$fold_neg[1], ")\n")
  cat("  Full results: tables/motif_hit_rate_descriptive.csv\n")
} else {
  cat("  Skipped (motif packages not available)\n")
}

cat("\n--- Output files ---\n")
cat("  tables/cpg_island_distribution.csv\n")
cat("  tables/motif_hit_rate_descriptive.csv\n")

log_info("=== 07_10c done ===")