#!/usr/bin/env Rscript
# ============================================================
# 07_1b_atac_motif_enrichment.r  (v5, monaLisa implementation)
#
# Purpose:
#   TF motif enrichment in greenyellow_closed promoters (RPE, GSE99287)
#   using monaLisa, with a reproducible design:
#     - Primary analysis: fixed foreground vs fixed background bin
#       (82 greenyellow_closed promoters vs 82 randomly sampled
#        promoters from the GSE99287 promoter pool, seed = 42).
#     - Sensitivity analysis: 100 random background draws, recording
#       each motif's FDR and counting how often FDR < 0.05.
#   Only motifs that are significant in the primary analysis AND in
#   >= 80% of random-background draws are called "robust".
#
# Changes in v5 (relative to v4):
#   - The entire output directory is deleted at start, so no legacy
#     v2/v3 files (genome/promoter background tables, two-background
#     figures) remain.
#   - Header output list matches the actual v5 outputs.
#   - Results draft mentions the most reproducible nominal signal
#     (e.g. ESRRA) even when n_robust == 0.
#   - sens_summary groups by motif_id only; motif_name is re-attached
#     from pwms, avoiding NA-group warnings.
#
# Inputs:
#   results/07_1_atac_validation/data/promo_long_rpe.rds
#   results/05_1_consensus_wgcna/tables/consensus_module_genes_annotated.csv
#   results/01_1_prepare_bulk_data/data/gene_annotation.rds
#
# Outputs (results/07_1b_atac_motif_enrichment/):
#   bed/greenyellow_closed_RPE.hg19.bed
#   bed/random_background_82.bed
#   tables/motif_enrichment_primary.csv
#   tables/motif_sensitivity_100_random_backgrounds.csv
#   tables/motif_robust_summary.csv
#   figures/motif_enrichment_primary.pdf
#   figures/motif_sensitivity_frequency.pdf
#   text/methods_motif_draft.txt
#   text/results_motif_draft.txt
# ============================================================

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(rtracklayer)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(JASPAR2024)
  library(TFBSTools)
  library(RSQLite)
  library(monaLisa)
  library(SummarizedExperiment)
})

log_info("=== 07_1b_atac_motif_enrichment.r (v5, monaLisa) start ===")

# ------------------------------------------------------------
# 0. Paths and parameters
# ------------------------------------------------------------
RESULT_ROOT  <- "E:/retina_tools/aging_rgc_axon/AMD/results"
OUT_DIR      <- file.path(RESULT_ROOT, "07_1b_atac_motif_enrichment")

# --- Clean slate: remove the entire previous output directory ---
if (dir.exists(OUT_DIR)) {
  unlink(OUT_DIR, recursive = TRUE, force = TRUE)
  log_info("removed previous output directory: {OUT_DIR}")
}

BED_DIR      <- file.path(OUT_DIR, "bed")
TAB_DIR      <- file.path(OUT_DIR, "tables")
FIG_DIR      <- file.path(OUT_DIR, "figures")
TXT_DIR      <- file.path(OUT_DIR, "text")
for (d in c(BED_DIR, TAB_DIR, FIG_DIR, TXT_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

PROMO_RDS <- file.path(RESULT_ROOT, "07_1_atac_validation", "data",
                        "promo_long_rpe.rds")
MOD_CSV   <- file.path(RESULT_ROOT, "05_1_consensus_wgcna", "tables",
                        "consensus_module_genes_annotated.csv")
ANNO_RDS  <- file.path(RESULT_ROOT, "01_1_prepare_bulk_data", "data",
                        "gene_annotation.rds")

FDR_CLOSED   <- 0.05
MODULE       <- "greenyellow"
N_RANDOM_BG  <- 100
ROBUST_FREQ  <- 0.80
SEED         <- 42
set.seed(SEED)

stopifnot(file.exists(PROMO_RDS), file.exists(MOD_CSV), file.exists(ANNO_RDS))

# ------------------------------------------------------------
# 1. Load inputs
# ------------------------------------------------------------
log_info("--- 1. Load inputs ---")

promo_rpe <- as.data.table(readRDS(PROMO_RDS))
mod_df    <- read.csv(MOD_CSV)
anno      <- readRDS(ANNO_RDS)

log_info("promo_long_rpe: {nrow(promo_rpe)} rows, {uniqueN(promo_rpe$gene)} genes")

# ------------------------------------------------------------
# 2. Map module genes to symbols
# ------------------------------------------------------------
log_info("--- 2. Map module genes to symbols ---")

ens <- sub("\\..*$", "", mod_df$ENSEMBL[mod_df$Module == MODULE])
idx <- match(ens, sub("\\..*$", "", anno$ensembl))
gy_genes <- toupper(anno$symbol[idx])
gy_genes <- unique(gy_genes[!is.na(gy_genes) & gy_genes != ""])

log_info("greenyellow genes after mapping: {length(gy_genes)}")
log_info("found in promo_long_rpe: {sum(gy_genes %in% promo_rpe$gene)}")

# ------------------------------------------------------------
# 3. Representative closed peak per gene
# ------------------------------------------------------------
log_info("--- 3. Build greenyellow_closed representative peaks ---")

gy_closed <- promo_rpe[
  gene %in% gy_genes & coef < 0 & fdr < FDR_CLOSED,
  .SD[which.min(coef)],
  by = gene
]
N_PEAKS <- nrow(gy_closed)
log_info("greenyellow_closed representative peaks: {N_PEAKS}")
if (N_PEAKS == 0) stop("No greenyellow_closed peaks")

gy_bed <- data.frame(
  chr    = gy_closed$Chr,
  start  = gy_closed$Hg19_start,
  end    = gy_closed$Hg19_end,
  name   = gy_closed$gene,
  score  = 0,
  strand = "."
)
bed_main <- file.path(BED_DIR, "greenyellow_closed_RPE.hg19.bed")
write.table(gy_bed, bed_main, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
log_info("written: {bed_main} ({nrow(gy_bed)} rows)")

# ------------------------------------------------------------
# 4. Extract foreground sequences
# ------------------------------------------------------------
log_info("--- 4. Extract foreground sequences ---")

gy_gr <- GRanges(seqnames = gy_closed$Chr,
                 ranges = IRanges(start = gy_closed$Hg19_start,
                                  end   = gy_closed$Hg19_end))
names(gy_gr) <- gy_closed$gene
seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, gy_gr)
names(seqs) <- gy_closed$gene
log_info("foreground sequences: {length(seqs)}")

# ------------------------------------------------------------
# 5. JASPAR2024 human CORE PWMs
# ------------------------------------------------------------
log_info("--- 5. Load JASPAR2024 human CORE PWMs ---")

jaspar <- JASPAR2024()
sq24 <- RSQLite::dbConnect(RSQLite::SQLite(), db(jaspar))
pwms_pf <- TFBSTools::getMatrixSet(
  sq24,
  list(species = "Homo sapiens", collection = "CORE", all_versions = FALSE)
)
RSQLite::dbDisconnect(sq24)
pwms <- TFBSTools::toPWM(pwms_pf)
log_info("PWMs loaded: {length(pwms)}")

# ------------------------------------------------------------
# 6. Helpers
# ------------------------------------------------------------
extract_motif_res <- function(se, pwms, fg_level = "foreground") {
  res <- as.data.frame(assay(se, "negLog10Padj"))
  if (is.null(colnames(res))) {
    fg_col <- 1L
  } else {
    fg_col <- which(colnames(res) == fg_level)
    if (length(fg_col) == 0) fg_col <- ncol(res)
  }
  out <- data.frame(
    motif_id     = rownames(res),
    negLog10Padj = res[[fg_col]],
    stringsAsFactors = FALSE
  )
  out$motif_name <- name(pwms)[match(out$motif_id, names(pwms))]
  out$FDR        <- 10^(-out$negLog10Padj)
  out <- out[!is.na(out$motif_name), ]
  out[order(-out$negLog10Padj),
      c("motif_id","motif_name","negLog10Padj","FDR")]
}

fmt_fdr <- function(x) {
  if (length(x) == 0 || is.na(x)) return("NA")
  if (x < 0.001) sprintf("%.2e", x) else sprintf("%.3f", x)
}

run_fg_vs_bg <- function(fg_seqs, bg_seqs, pwms) {
  all_seqs <- c(fg_seqs, bg_seqs)
  all_bins <- factor(c(rep("foreground", length(fg_seqs)),
                       rep("background", length(bg_seqs))),
                     levels = c("background", "foreground"))
  se <- calcBinnedMotifEnrR(
    seqs = all_seqs,
    bins = all_bins,
    pwmL = pwms,
    background = "otherBins",
    verbose = FALSE
  )
  extract_motif_res(se, pwms, fg_level = "foreground")
}

# ------------------------------------------------------------
# 7. Fixed background (primary analysis)
# ------------------------------------------------------------
log_info("--- 7. Primary: fixed background (seed = {SEED}) ---")

all_promoter <- promo_rpe[Type == "promoter"]

set.seed(SEED)
bg_idx <- sample(nrow(all_promoter), N_PEAKS)
bg_bed <- data.frame(
  chr    = all_promoter$Chr[bg_idx],
  start  = all_promoter$Hg19_start[bg_idx],
  end    = all_promoter$Hg19_end[bg_idx],
  name   = paste0("bg_", seq_along(bg_idx)),
  score  = 0,
  strand = "."
)
bed_bg <- file.path(BED_DIR, "random_background_82.bed")
write.table(bg_bed, bed_bg, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
log_info("fixed background written: {bed_bg}")

bg_gr <- GRanges(seqnames = bg_bed$chr,
                 ranges = IRanges(start = bg_bed$start, end = bg_bed$end))
bg_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg19, bg_gr)
names(bg_seqs) <- bg_bed$name

res_primary <- run_fg_vs_bg(seqs, bg_seqs, pwms)
write.csv(res_primary,
          file.path(TAB_DIR, "motif_enrichment_primary.csv"),
          row.names = FALSE)

n_primary_sig <- sum(res_primary$FDR < 0.05, na.rm = TRUE)
log_info("primary: motifs FDR<0.05 = {n_primary_sig}")

# ------------------------------------------------------------
# 8. Sensitivity: N_RANDOM_BG random background draws
# ------------------------------------------------------------
log_info("--- 8. Sensitivity: {N_RANDOM_BG} random backgrounds ---")

sens_list <- vector("list", N_RANDOM_BG)
for (i in seq_len(N_RANDOM_BG)) {
  set.seed(SEED + i)
  idx_i <- sample(nrow(all_promoter), N_PEAKS)
  bg_gr_i <- GRanges(seqnames = all_promoter$Chr[idx_i],
                     ranges = IRanges(start = all_promoter$Hg19_start[idx_i],
                                      end   = all_promoter$Hg19_end[idx_i]))
  bg_seqs_i <- getSeq(BSgenome.Hsapiens.UCSC.hg19, bg_gr_i)
  names(bg_seqs_i) <- paste0("bg_", i, "_", seq_along(bg_gr_i))

  r <- run_fg_vs_bg(seqs, bg_seqs_i, pwms)
  r$iter <- i
  sens_list[[i]] <- r[, c("iter","motif_id","motif_name","FDR")]
  if (i %% 10 == 0) log_info("sensitivity {i}/{N_RANDOM_BG} done")
}
sens <- rbindlist(sens_list)
write.csv(sens,
          file.path(TAB_DIR, "motif_sensitivity_100_random_backgrounds.csv"),
          row.names = FALSE)

# Group by motif_id only to avoid NA-group warnings; re-attach motif_name
sens_summary <- sens[, .(
  n_iter     = .N,
  n_sig      = sum(FDR < 0.05, na.rm = TRUE),
  freq_sig   = mean(FDR < 0.05, na.rm = TRUE),
  median_FDR = median(FDR, na.rm = TRUE),
  min_FDR    = min(FDR, na.rm = TRUE)
), by = motif_id]
sens_summary$motif_name <- name(pwms)[match(sens_summary$motif_id, names(pwms))]
sens_summary <- sens_summary[, c("motif_id","motif_name","n_iter","n_sig",
                                 "freq_sig","median_FDR","min_FDR")]

# ------------------------------------------------------------
# 9. Robust motif summary
# ------------------------------------------------------------
log_info("--- 9. Robust motif summary ---")

robust <- merge(
  res_primary[, c("motif_id","motif_name","FDR")],
  sens_summary[, c("motif_id","motif_name","freq_sig","n_sig",
                   "median_FDR","min_FDR")],
  by = c("motif_id","motif_name"),
  all.x = TRUE
)
robust$primary_sig <- robust$FDR < 0.05
robust$robust      <- robust$primary_sig & robust$freq_sig >= ROBUST_FREQ
robust <- robust[order(-robust$freq_sig, robust$FDR), ]
write.csv(robust,
          file.path(TAB_DIR, "motif_robust_summary.csv"),
          row.names = FALSE)

n_robust <- sum(robust$robust, na.rm = TRUE)
log_info("robust motifs (primary FDR<0.05 & freq_sig>={ROBUST_FREQ}): {n_robust}")
if (n_robust > 0) {
  print(robust[robust$robust, ], row.names = FALSE, digits = 4)
}

# ------------------------------------------------------------
# 10. Figures
# ------------------------------------------------------------
log_info("--- 10. Figures ---")

# 10a. Primary barplot (top 20 by FDR)
top <- head(res_primary[order(-res_primary$negLog10Padj), ], 20)
top$motif_name <- factor(top$motif_name, levels = rev(top$motif_name))
top$sig <- ifelse(top$FDR < 0.05, "FDR<0.05", "ns")
pA <- ggplot(top, aes(x = motif_name, y = negLog10Padj, fill = sig)) +
  geom_col() + coord_flip() + theme_bw(base_size = 10) +
  scale_fill_manual(values = c("FDR<0.05" = "#C44E52", "ns" = "#B0B0B0")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(x = "", y = "-log10(FDR)", fill = "",
       title = "Motif enrichment, fixed background (primary)")
ggsave(file.path(FIG_DIR, "motif_enrichment_primary.pdf"),
       pA, width = 8, height = 6)

# 10b. Sensitivity frequency plot (top 20 by freq_sig)
top_sens <- head(sens_summary[order(-sens_summary$freq_sig), ], 20)
top_sens$motif_name <- factor(top_sens$motif_name,
                               levels = rev(top_sens$motif_name))
pB <- ggplot(top_sens, aes(x = motif_name, y = freq_sig)) +
  geom_col(fill = "#4C72B0") +
  coord_flip() + theme_bw(base_size = 10) +
  geom_hline(yintercept = ROBUST_FREQ, linetype = "dashed", color = "#C44E52") +
  labs(x = "", y = "Fraction of random backgrounds with FDR < 0.05",
       title = "Sensitivity of motif enrichment to background choice")
ggsave(file.path(FIG_DIR, "motif_sensitivity_frequency.pdf"),
       pB, width = 8, height = 6)

# ------------------------------------------------------------
# 11. Draft Methods / Results text
# ------------------------------------------------------------
log_info("--- 11. Draft text ---")

top_primary <- res_primary[which.min(res_primary$FDR), ]
top_sens    <- sens_summary[which.max(sens_summary$freq_sig), ]

methods_txt <- c(
  "Methods (motif enrichment).",
  "TF motif enrichment was performed with monaLisa (v1.18.0) using the",
  "JASPAR2024 human CORE PWM collection (720 motifs). Foreground sequences",
  sprintf("corresponded to the %d greenyellow-closed promoter peaks", N_PEAKS),
  "(representative peak per gene; coef < 0 and FDR < 0.05). To avoid the",
  "instability of genome-wide background sampling with small peak sets, we",
  "used a fixed foreground-vs-background design: the background bin comprised",
  sprintf("%d randomly sampled promoter peaks from the same GSE99287 promoter", N_PEAKS),
  sprintf("pool (seed = %d), written to bed/random_background_82.bed.", SEED),
  "Enrichment was tested by one-sided Fisher's exact test and P values were",
  "BH-adjusted across motifs. To assess robustness to background choice, the",
  sprintf("analysis was repeated with %d independent random background draws,", N_RANDOM_BG),
  "and motifs were called robust if primary FDR < 0.05 and FDR < 0.05",
  sprintf("in >= %d%% of draws.", round(100 * ROBUST_FREQ))
)

if (n_robust == 0) {
  results_txt <- c(
    "Results (motif enrichment).",
    sprintf("In the fixed-background primary analysis, %d motif(s) reached FDR < 0.05;",
            n_primary_sig),
    sprintf("the strongest was %s (FDR = %s).",
            top_primary$motif_name, fmt_fdr(top_primary$FDR)),
    sprintf("Across %d random background draws, no motif reached FDR < 0.05 in >= %d%% of draws.",
            N_RANDOM_BG, round(100 * ROBUST_FREQ)),
    sprintf("The most reproducible nominal signal was %s, significant in %d/%d draws,",
            top_sens$motif_name, top_sens$n_sig, N_RANDOM_BG),
    sprintf("which is below the prespecified %d%% robustness threshold.",
            round(100 * ROBUST_FREQ)),
    "No TF motif showed a robust enrichment signal for greenyellow promoter",
    "closure. This is consistent with a distributed, epigenetically mediated",
    "silencing rather than loss of a single TF-binding program."
  )
} else {
  results_txt <- c(
    "Results (motif enrichment).",
    sprintf("In the fixed-background primary analysis, %d motif(s) reached FDR < 0.05;",
            n_primary_sig),
    sprintf("the strongest was %s (FDR = %s).",
            top_primary$motif_name, fmt_fdr(top_primary$FDR)),
    sprintf("Across %d random background draws, %d motif(s) reached FDR < 0.05",
            N_RANDOM_BG, n_robust),
    sprintf("in >= %d%% of draws and in the primary analysis.",
            round(100 * ROBUST_FREQ)),
    sprintf("The most reproducible signal was %s, significant in %d/%d draws.",
            top_sens$motif_name, top_sens$n_sig, N_RANDOM_BG),
    "This motif should be interpreted as a candidate regulatory feature",
    "requiring independent validation, not as a causal driver."
  )
}

writeLines(methods_txt, file.path(TXT_DIR, "methods_motif_draft.txt"))
writeLines(results_txt, file.path(TXT_DIR, "results_motif_draft.txt"))

# ------------------------------------------------------------
# 12. Summary
# ------------------------------------------------------------
cat("\n", strrep("=", 64), "\n", sep = "")
cat("Summary: 07_1b motif enrichment (monaLisa, v5)\n")
cat(strrep("=", 64), "\n", sep = "")
cat(sprintf("greenyellow_closed peaks:                 %d\n", N_PEAKS))
cat(sprintf("PWMs tested:                              %d\n", length(pwms)))
cat(sprintf("primary (fixed background), FDR < 0.05:   %d\n", n_primary_sig))
cat(sprintf("random background draws:                  %d\n", N_RANDOM_BG))
cat(sprintf("robust (primary + >=%.0f%% draws):         %d\n",
            100 * ROBUST_FREQ, n_robust))
cat(sprintf("primary top motif:                        %s (FDR = %s)\n",
            top_primary$motif_name, fmt_fdr(top_primary$FDR)))
cat(sprintf("most reproducible motif:                  %s (%d/%d draws)\n",
            top_sens$motif_name, top_sens$n_sig, N_RANDOM_BG))
cat(sprintf("output directory:                         %s\n", OUT_DIR))
cat(strrep("=", 64), "\n", sep = "")

log_info("=== 07_1b_atac_motif_enrichment.r (v5, monaLisa) done ===")

