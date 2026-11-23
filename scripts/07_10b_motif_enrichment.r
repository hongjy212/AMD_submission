#!/usr/bin/env Rscript
# ============================================================
# 07_10b_motif_enrichment.r (slim, archived, logger, v5)
#
# STATUS: MIXED RESULT
#   1 GC-rich motif (KLF15) significant in pos_beta,
#   attributable to CpG-island background, not TF-specific.
#
# History:
#   - v3 (flank = 50 bp):  0 significant motifs
#   - v4 (flank = 250 bp): 0 significant motifs (DFrame bug masked 1 row)
#   - v5 (flank = 250 bp): 1 significant motif in pos_beta (KLF15)
#
# Test: universalmotif::enrich_motifs() with JASPAR2024 human CORE
#       (720 motifs), background = 2000 eQTM CpGs,
#       max.p = 0.05, max.q = 0.20, threshold = 0.001, RC = TRUE
#
# Targets:
#   - 55 negative-beta CpGs  (methylation up -> expression down)
#   - 46 positive-beta CpGs  (methylation up -> expression up)
#   - 46 promoter CpGs       (TSS200/TSS1500)   [not tested in v5]
#   - 101 all greenyellow CpGs                  [not tested in v5]
#
# Result:
#   neg_beta : 0 significant motifs
#   pos_beta : 1 significant motif (KLF15, fold = 1.27, Q = 6.8e-05)
#              KLF15 consensus CCCCGCCC is purely GC-rich;
#              the same motif is present in 81.8% of negative-beta
#              CpGs (07_10c), indicating CpG-island background
#              rather than specific KLF15 binding.
#
# Conclusion:
#   greenyellow methylation is not driven by a single TF motif.
#   The KLF15 signal is attributable to CpG-island GC content.
#   Reported as a mixed (GC-rich) result in Supplementary Results S4.
#   The CpG island distribution (07_10c) provides the positive
#   functional annotation (43.6% in CpG islands vs 34.4% background).
#
# Input : results/07_10_greenyellow_methylation_upstream/bed/*.bed
# Output: results/07_10b_motif_enrichment/
#           tables/motif_enrichment.csv   (only if nrow > 0)
#           07_10b_summary.txt
#           07_10b_motif_enrichment.log
# ============================================================

# ------------------------------------------------------------
# 0. Setup: logger + paths (do NOT source 00_setup.r to avoid
#    log_info name clash with logger)
# ------------------------------------------------------------
suppressPackageStartupMessages({
  library(logger)
  library(universalmotif); library(JASPAR2024); library(TFBSTools)
  library(RSQLite); library(BSgenome.Hsapiens.UCSC.hg19)
  library(GenomicRanges); library(IRanges); library(Biostrings)
  library(dplyr)
})

DIR_ROOT   <- "E:/retina_tools/aging_rgc_axon/AMD"
DIR_RESULT <- file.path(DIR_ROOT, "results")
DIR_THIS   <- file.path(DIR_RESULT, "07_10b_motif_enrichment")
BED_DIR    <- file.path(DIR_RESULT,
                         "07_10_greenyellow_methylation_upstream", "bed")
SEED       <- 20260925

for (d in c(DIR_THIS,
            file.path(DIR_THIS, "tables"),
            file.path(DIR_THIS, "figures"))) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

log_file <- file.path(DIR_THIS, "07_10b_motif_enrichment.log")
log_appender(appender_tee(log_file))
log_layout(layout_glue_generator(format = '{time} [{level}] {msg}'))

log_info("=== 07_10b v5 (mixed result) start ===")
log_info("DIR_THIS = {DIR_THIS}")
log_info("BED_DIR  = {BED_DIR}")
log_info("SEED     = {SEED}")

# ------------------------------------------------------------
# 1. Parameters
# ------------------------------------------------------------
FLANK       <- 250
N_BG_SAMPLE <- 2000
THRESHOLD   <- 0.001
MAX_P       <- 0.05
MAX_Q       <- 0.20

log_info("FLANK = {FLANK} bp, N_BG_SAMPLE = {N_BG_SAMPLE}, ",
         "THRESHOLD = {THRESHOLD}, MAX_P = {MAX_P}, MAX_Q = {MAX_Q}")

# ------------------------------------------------------------
# 2. Read BED and extract sequences (hg19)
# ------------------------------------------------------------
read_bed_gr <- function(f) {
  d <- read.delim(f, header = FALSE, stringsAsFactors = FALSE)
  GRanges(seqnames = d$V1,
          ranges   = IRanges(start = d$V2 + 1, end = d$V3),
          strand   = d$V6, name = d$V4)
}

expand_gr <- function(gr, flank) {
  start(gr) <- pmax(1, start(gr) - flank)
  end(gr)   <- end(gr) + flank
  gr
}

bed_files <- c(
  neg_beta        = "target_neg_beta.bed",
  pos_beta        = "target_pos_beta.bed",
  promoter_only   = "target_promoter_only.bed",
  all_greenyellow = "target_all_greenyellow.bed",
  background      = "background_all_eqtm.bed"
)

gr_list <- list()
for (nm in names(bed_files)) {
  f <- file.path(BED_DIR, bed_files[[nm]])
  if (!file.exists(f)) {
    log_warn("missing BED: {f}")
    next
  }
  gr_list[[nm]] <- read_bed_gr(f)
  log_info("{nm}: {length(gr_list[[nm]])} regions")
}
stopifnot("background" %in% names(gr_list))

gr_list <- lapply(gr_list, expand_gr, flank = FLANK)

get_seq   <- function(gr) getSeq(BSgenome.Hsapiens.UCSC.hg19, gr)
clean_seq <- function(dss) dss[!grepl("N", as.character(dss), fixed = TRUE)]

set.seed(SEED)
seq_bg_full <- clean_seq(get_seq(gr_list$background))
seq_bg      <- seq_bg_full[sample(length(seq_bg_full),
                                    min(N_BG_SAMPLE, length(seq_bg_full)))]
log_info("background sequences: {length(seq_bg_full)} -> {length(seq_bg)}")

seq_targets <- lapply(
  gr_list[c("neg_beta","pos_beta","promoter_only","all_greenyellow")],
  function(gr) clean_seq(get_seq(gr)))
for (nm in names(seq_targets)) {
  log_info("{nm}: {length(seq_targets[[nm]])} sequences")
}

# ------------------------------------------------------------
# 3. JASPAR2024 human CORE motifs
# ------------------------------------------------------------
jaspar <- JASPAR2024::JASPAR2024()
sq24   <- RSQLite::dbConnect(RSQLite::SQLite(), db(jaspar))
motifs <- TFBSTools::getMatrixSet(sq24,
            list(species = "Homo sapiens", collection = "CORE"))
RSQLite::dbDisconnect(sq24)
motifs_um <- universalmotif::convert_motifs(motifs)
log_info("JASPAR2024 CORE motifs: {length(motifs_um)}")

oxphos_tfs <- c("NRF1","NRF2","NFE2L2","ESRRA","ESRRB","PPARGC1A",
                "SP1","SP2","SP3","SP4","KLF4","KLF9","KLF15",
                "E2F1","E2F2","E2F4","YY1","MYC","SREBF1","SREBF2",
                "GABPA","TFAM","PPARG","NR1H3","NRF")
motif_names <- sapply(motifs_um, function(m) m@name)
n_ox <- sum(grepl(paste(oxphos_tfs, collapse = "|"),
                    motif_names, ignore.case = TRUE))
log_info("OXPHOS-related motifs present in JASPAR2024 CORE: {n_ox}")

# ------------------------------------------------------------
# 4. Motif enrichment
# ------------------------------------------------------------
run_enrich <- function(motifs_sub, seq_target, seq_bg, label) {
  log_info(">>> {label}: {length(motifs_sub)} motifs x {length(seq_target)} targets")
  t0 <- Sys.time()
  res <- tryCatch(
    enrich_motifs(
      motifs        = motifs_sub,
      sequences     = seq_target,
      bkg.sequences = seq_bg,
      max.p         = MAX_P,
      max.q         = MAX_Q,
      threshold     = THRESHOLD,
      RC            = TRUE,
      verbose       = 0
    ),
    error = function(e) {
      log_error("failed: {label} - {conditionMessage(e)}")
      NULL
    }
  )
  dt <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
  if (is.null(res)) {
    log_info("<<< {label}: failed after {dt}s")
    return(NULL)
  }
  n_sig <- if (is.null(nrow(res))) 0L else as.integer(nrow(res))
  log_info("<<< {label}: {n_sig} significant motifs ({dt}s)")
  if (n_sig == 0L) return(NULL)
  res$target_set <- label
  res
}

enrich_results <- list()
for (nm in c("neg_beta", "pos_beta")) {
  r <- run_enrich(motifs_um, seq_targets[[nm]], seq_bg, paste0(nm, "_all"))
  if (!is.null(r)) enrich_results[[nm]] <- r
}

# ------------------------------------------------------------
# 5. Save (only if non-empty)
# ------------------------------------------------------------
if (length(enrich_results) > 0) {
  combined <- bind_rows(lapply(enrich_results, as.data.frame))
  utils::write.csv(combined,
                   file.path(DIR_THIS, "tables", "motif_enrichment.csv"),
                   row.names = FALSE, fileEncoding = "UTF-8")
  log_info("saved: tables/motif_enrichment.csv ({nrow(combined)} rows)")
} else {
  log_info("No significant motif in any target set")
  log_info("No CSV written (no rows to report)")
}

# ------------------------------------------------------------
# 6. Interpretation: flag GC-rich motifs explicitly
# ------------------------------------------------------------
interpretation <- character(0)

if (!is.null(enrich_results$pos_beta)) {
  pos_df <- as.data.frame(enrich_results$pos_beta)
  for (i in seq_len(nrow(pos_df))) {
    m_name <- pos_df$motif[i]
    m_cons <- pos_df$motif.consensus[i]
    fold   <- pos_df$target.enrichment[i]
    qval   <- pos_df$Qval[i]
    # GC-rich flag: consensus composed of >= 80% G/C
    gc_frac <- nchar(gsub("[^GCgc]", "", m_cons)) / nchar(m_cons)
    is_gc_rich <- gc_frac >= 0.8
    interpretation <- c(interpretation, sprintf(
      "pos_beta: %s (consensus %s, fold = %.2f, Q = %.2g) - %s",
      m_name, m_cons, fold, qval,
      if (is_gc_rich) "GC-rich motif; likely CpG-island background, not TF-specific"
      else "non-GC-rich motif; consider as candidate TF"
    ))
  }
}

if (is.null(enrich_results$neg_beta)) {
  interpretation <- c(interpretation, "neg_beta: no significant motif")
}
if (is.null(enrich_results$pos_beta)) {
  interpretation <- c(interpretation, "pos_beta: no significant motif")
}

for (line in interpretation) log_info(line)

# ------------------------------------------------------------
# 7. Write summary file for archival
# ------------------------------------------------------------
summary_txt <- c(
  "07_10b_motif_enrichment - mixed result",
  paste0("Date: ", Sys.Date()),
  "",
  "Test:",
  "  universalmotif::enrich_motifs",
  "  JASPAR2024 human CORE (720 motifs)",
  "  background = 2000 eQTM CpGs",
  paste0("  max.p = ", MAX_P, ", max.q = ", MAX_Q,
         ", threshold = ", THRESHOLD, ", RC = TRUE"),
  paste0("  flank = ", FLANK, " bp"),
  "",
  "Targets:",
  paste0("  neg_beta        : ", length(seq_targets$neg_beta), " CpGs"),
  paste0("  pos_beta        : ", length(seq_targets$pos_beta), " CpGs"),
  paste0("  promoter_only   : ", length(seq_targets$promoter_only),
         " CpGs (not tested in v5)"),
  paste0("  all_greenyellow : ", length(seq_targets$all_greenyellow),
         " CpGs (not tested in v5)"),
  "",
  "Result:",
  if (!is.null(enrich_results$neg_beta)) {
    paste0("  neg_beta        : ", nrow(as.data.frame(enrich_results$neg_beta)),
           " significant motif(s)")
  } else {
    "  neg_beta        : 0 significant motifs"
  },
  if (!is.null(enrich_results$pos_beta)) {
    pos_df <- as.data.frame(enrich_results$pos_beta)
    c(
      paste0("  pos_beta        : ", nrow(pos_df), " significant motif(s)"),
      paste0("                    ", pos_df$motif, " (consensus ",
             pos_df$motif.consensus, ", fold = ",
             round(pos_df$target.enrichment, 2), ", Q = ",
             signif(pos_df$Qval, 3), ")"),
      "                    KLF15 consensus CCCCGCCC is purely GC-rich;",
      "                    the same motif is present in 81.8% of negative-beta",
      "                    CpGs (07_10c), indicating CpG-island background",
      "                    rather than specific KLF15 binding."
    )
  } else {
    "  pos_beta        : 0 significant motifs"
  },
  "",
  "Conclusion:",
  "  greenyellow methylation is not driven by a single TF motif.",
  "  The KLF15 signal is attributable to CpG-island GC content.",
  "  Reported as a mixed (GC-rich) result in Supplementary Results S4.",
  "  The CpG island distribution (07_10c) provides the positive",
  "  functional annotation (43.6% in CpG islands vs 34.4% background)."
)
writeLines(summary_txt, file.path(DIR_THIS, "07_10b_summary.txt"))
log_info("saved: 07_10b_summary.txt")

log_info("=== 07_10b v5 done ===")

