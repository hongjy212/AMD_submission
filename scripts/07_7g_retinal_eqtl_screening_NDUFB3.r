# ============================================================
# 07_7g_retinal_eqtl_screening_NDUFB3.r  (v2.1)
# Purpose: Screen four independent retinal eQTL resources
#          (Orozco eye-eQTL, Strunz mega-analysis, EyeGEx,
#           Wang et al. 2023 sc-multiomics) for NDUFB3 cis-eQTL,
#          and produce a unified summary table.
#
# v2.1 Fix：
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(dplyr)
})

ROOT <- "E:/retina_tools/aging_rgc_axon/AMD"
OUT  <- file.path(ROOT, "results", "07_7g_retinal_eqtl_screening_NDUFB3")
TAB  <- file.path(OUT, "tables")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)

log_info <- function(...) cat(sprintf("[%s] [INFO] %s\n",
  format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)))
log_warn <- function(...) cat(sprintf("[%s] [WARN] %s\n",
  format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)))

GENE_ENSEMBL <- "ENSG00000119013"
GENE_SYMBOL  <- "NDUFB3"

NDUFB3_CIS_CHR   <- 2L
NDUFB3_CIS_START <- 201850000L
NDUFB3_CIS_END   <- 203860000L

log_info("=== 07_7g retinal eQTL screening for NDUFB3 (v2.1) ===")

# ============================================================
# 1. Orozco eye-eQTL
# ============================================================
log_info("--- 1. Orozco eye-eQTL ---")

orozco_in <- "E:/retina_tools/pubdata/eQTL/GSE135092/eye_eqtl.csv"

orozco_out <- NULL
if (file.exists(orozco_in)) {
  orozco_raw <- fread(orozco_in)
  orozco_out <- orozco_raw[Gene == GENE_SYMBOL]
  orozco_out[, Usable_for_MR := ifelse(
    MAF >= 0.05 & (is.na(Notes) | !grepl("Pleiotropic", Notes)),
    "Yes", "No"
  )]
  fwrite(orozco_out, file.path(TAB, "NDUFB3_Orozco_screening.csv"))
  log_info(sprintf("Orozco: %d rows saved", nrow(orozco_out)))
} else {
  log_warn(sprintf("Orozco input not found: %s", orozco_in))
}

# ============================================================
# 2. Strunz
# ============================================================
log_info("--- 2. Strunz S1 Data ---")

strunz_in <- "E:/retina_tools/pubdata/eQTL/Strunz/pgen.1008934.s005.xlsx"

strunz_out <- NULL
if (file.exists(strunz_in)) {
  x_strunz <- read_excel(strunz_in, sheet = "S1", skip = 2)
  strunz_out <- as.data.table(x_strunz)[eGene == GENE_ENSEMBL]
  if (nrow(strunz_out) == 0) {
    strunz_out <- data.table(
      eVariant = NA_character_, eGene = GENE_ENSEMBL,
      `eGene Symbol` = GENE_SYMBOL, detected = "No")
    log_info("Strunz: NDUFB3 not among independent eQTL signals")
  } else {
    strunz_out[, detected := "Yes"]
    log_info(sprintf("Strunz: NDUFB3 detected in %d rows", nrow(strunz_out)))
  }
  fwrite(strunz_out, file.path(TAB, "NDUFB3_Strunz_screening.csv"))
} else {
  log_warn(sprintf("Strunz input not found: %s", strunz_in))
}

# ============================================================
# 3. EyeGEx cis-independent eQTLs
# ============================================================
log_info("--- 3. EyeGEx cis_independent_eQTLs ---")

eyeg_in <- "E:/retina_tools/pubdata/eQTL/GSE115828/NIHMS1518564-supplement-4.xlsx"

eyeg_cis_out <- NULL
if (file.exists(eyeg_in)) {
  x_eyeg <- read_excel(eyeg_in, sheet = "cis_independent_eQTLs", skip = 3)
  eyeg_cis_out <- as.data.table(x_eyeg)[external_gene_name == GENE_SYMBOL]
  if (nrow(eyeg_cis_out) == 0) {
    eyeg_cis_out <- data.table(
      gene_id = GENE_ENSEMBL, external_gene_name = GENE_SYMBOL,
      detected = "No")
    log_info("EyeGEx: NDUFB3 not among independent cis-eQTLs")
  } else {
    eyeg_cis_out[, detected := "Yes"]
    log_info(sprintf("EyeGEx: NDUFB3 detected in %d rows", nrow(eyeg_cis_out)))
  }
  fwrite(eyeg_cis_out, file.path(TAB, "NDUFB3_EyeGEx_cis_independent_screening.csv"))
} else {
  log_warn(sprintf("EyeGEx input not found: %s", eyeg_in))
}

# ============================================================
# 4. EyeGEx eGenes
# ============================================================
log_info("--- 4. EyeGEx eGenes ---")

eyeg_eg_out <- NULL
if (file.exists(eyeg_in)) {
  x_eg <- read_excel(eyeg_in, sheet = "eGenes", skip = 3)
  eyeg_eg_out <- as.data.table(x_eg)[external_gene_name == GENE_SYMBOL]
  if (nrow(eyeg_eg_out) == 0) {
    eyeg_eg_out <- data.table(
      gene = GENE_ENSEMBL, external_gene_name = GENE_SYMBOL,
      detected = "No")
    log_info("EyeGEx eGenes: NDUFB3 not present")
  } else {
    eyeg_eg_out[, detected := "Yes"]
    log_info(sprintf("EyeGEx eGenes: NDUFB3 detected in %d rows",
                     nrow(eyeg_eg_out)))
  }
  fwrite(eyeg_eg_out, file.path(TAB, "NDUFB3_EyeGEx_eGenes_screening.csv"))
} else {
  log_warn(sprintf("EyeGEx input not found: %s", eyeg_in))
}

# ============================================================
# 5. Wang et al. 2023 sc-multiomics
# ============================================================
log_info("--- 5. Wang 2023 sc-multiomics ---")

wang_xlsx <- "E:/retina_tools/pubdata/eQTL/sc_eQTL/13059_2023_3111_MOESM2_ESM.xlsx"

wang_s4_ndufb3 <- NULL
wang_s4_region <- NULL
wang_s8_region <- NULL

if (file.exists(wang_xlsx)) {

  # ---- 5a. Table S4: sc-eQTL (gene-level FDR < 0.1) ----
  s4_raw <- read_excel(wang_xlsx, sheet = "Table S4", col_names = FALSE)
  s4 <- as.data.table(s4_raw[-c(1, 2), ])
  setnames(s4, 1:15,
    c("celltype", "chr", "pos", "ref", "alt", "variant_id", "gene",
      "beta", "se", "pval", "tss_dis", "ma_samples", "ma_count",
      "maf", "FDR"))
  s4[, chr := as.integer(chr)]
  s4[, pos := as.numeric(pos)]
  s4[, beta := as.numeric(beta)]
  s4[, se := as.numeric(se)]
  s4[, maf := as.numeric(maf)]
  s4[, FDR := as.numeric(FDR)]

  wang_s4_ndufb3 <- s4[gene == GENE_SYMBOL]
  if (nrow(wang_s4_ndufb3) == 0) {
    log_info("Wang 2023 sc-eQTL: NDUFB3 not a sc-eGene in any cell type")
  } else {
    log_info(sprintf("Wang 2023 sc-eQTL: NDUFB3 detected in %d rows",
                     nrow(wang_s4_ndufb3)))
  }
  fwrite(wang_s4_ndufb3,
         file.path(TAB, "NDUFB3_Wang2023_sc-eQTL_NDUFB3.csv"))

  wang_s4_region <- s4[chr == NDUFB3_CIS_CHR &
                       pos >= NDUFB3_CIS_START &
                       pos <= NDUFB3_CIS_END]
  log_info(sprintf("Wang 2023 sc-eQTL: %d sc-eQTLs in NDUFB3 cis region",
                   nrow(wang_s4_region)))
  if (nrow(wang_s4_region) > 0) {
    log_info("  Cell types: %s",
             paste(sort(unique(wang_s4_region$celltype)), collapse = ", "))
    log_info("  Target genes: %s",
             paste(sort(unique(wang_s4_region$gene)), collapse = ", "))
  }
  fwrite(wang_s4_region,
         file.path(TAB, "NDUFB3_Wang2023_sc-eQTL_cis_region.csv"))

  # ---- 5b. Table S8: sc-caQTL (genome-level FDR < 0.1) ----
  s8_raw <- read_excel(wang_xlsx, sheet = "Table S8", col_names = FALSE)
  s8 <- as.data.table(s8_raw[-c(1, 2), ])
  setnames(s8, 1:5, c("celltype", "OCR_ID", "variant_id", "chr", "pos"))
  s8[, chr := as.integer(gsub("^chr", "", chr))]
  s8[, pos := as.numeric(pos)]

  wang_s8_region <- s8[chr == NDUFB3_CIS_CHR &
                       pos >= NDUFB3_CIS_START &
                       pos <= NDUFB3_CIS_END]
  log_info(sprintf("Wang 2023 sc-caQTL: %d sc-caQTLs in NDUFB3 cis region",
                   nrow(wang_s8_region)))
  if (nrow(wang_s8_region) > 0) {
    log_info("  Cell types: %s",
             paste(sort(unique(wang_s8_region$celltype)), collapse = ", "))
  }
  fwrite(wang_s8_region,
         file.path(TAB, "NDUFB3_Wang2023_sc-caQTL_cis_region.csv"))

} else {
  log_warn(sprintf("Wang 2023 input not found: %s", wang_xlsx))
}

# ============================================================
# 6. Unified summary table
# ============================================================
log_info("--- 6. Unified summary ---")

summary_rows <- list()

# ---- Orozco ----
if (!is.null(orozco_out) && nrow(orozco_out) > 0) {
  for (i in seq_len(nrow(orozco_out))) {
    r <- orozco_out[i]
    summary_rows[[length(summary_rows) + 1]] <- data.table(
      Source = "Orozco_eye_eqtl", Sample_size = 106L,
      Tissue = r$Tissue, SNP = r$SNP, FDR = r$FDR, Beta = r$Beta,
      MAF = r$MAF, eQTL_detected = "Yes",
      Usable_for_MR = r$Usable_for_MR, Notes = r$Notes)
  }
}

# ---- Strunz ----
if (!is.null(strunz_out)) {
  summary_rows[[length(summary_rows) + 1]] <- data.table(
    Source = "Strunz_mega_analysis", Sample_size = 311L,
    Tissue = "Retina", SNP = NA_character_, FDR = NA_real_,
    Beta = NA_real_, MAF = NA_real_,
    eQTL_detected = ifelse(
      nrow(strunz_out[strunz_out$detected == "Yes"]) > 0, "Yes", "No"),
    Usable_for_MR = "No",
    Notes = "Not among independent eQTL signals")
}

# ---- EyeGEx ----
if (!is.null(eyeg_cis_out)) {
  summary_rows[[length(summary_rows) + 1]] <- data.table(
    Source = "EyeGEx", Sample_size = 406L,
    Tissue = "Retina", SNP = NA_character_, FDR = NA_real_,
    Beta = NA_real_, MAF = NA_real_,
    eQTL_detected = ifelse(
      nrow(eyeg_cis_out[eyeg_cis_out$detected == "Yes"]) > 0, "Yes", "No"),
    Usable_for_MR = "No",
    Notes = "Not among independent cis-eQTLs")
}

# ---- Wang 2023 ----
if (!is.null(wang_s4_ndufb3)) {
  if (nrow(wang_s4_ndufb3) > 0) {
    for (i in seq_len(nrow(wang_s4_ndufb3))) {
      r <- wang_s4_ndufb3[i]
      summary_rows[[length(summary_rows) + 1]] <- data.table(
        Source = "Wang_2023_sc_multiomics", Sample_size = 20L,
        Tissue = sprintf("Retina_%s", r$celltype),
        SNP = r$variant_id, FDR = r$FDR, Beta = r$beta, MAF = r$maf,
        eQTL_detected = "Yes", Usable_for_MR = "No",
        Notes = "sc-eQTL in single-cell multiomics")
    }
  } else {
    n_region <- if (!is.null(wang_s4_region)) nrow(wang_s4_region) else 0L
    genes_region <- if (!is.null(wang_s4_region) && n_region > 0)
      paste(sort(unique(wang_s4_region$gene)), collapse = ",") else "none"
    n_caqtl <- if (!is.null(wang_s8_region)) nrow(wang_s8_region) else 0L

    summary_rows[[length(summary_rows) + 1]] <- data.table(
      Source = "Wang_2023_sc_multiomics", Sample_size = 20L,
      Tissue = "Retina (scRNA + scATAC)",
      SNP = NA_character_, FDR = NA_real_, Beta = NA_real_, MAF = NA_real_,
      eQTL_detected = "No", Usable_for_MR = "No",
      Notes = sprintf(
        "NDUFB3 not a sc-eGene in any major retinal cell type; %d sc-eQTLs in NDUFB3 cis region all target neighbours (%s); %d sc-caQTLs in cis region; PPIL3 sc-eQTLs replicated in AC/HC",
        n_region, genes_region, n_caqtl))
  }
}

if (length(summary_rows) > 0) {
  summary_tbl <- rbindlist(summary_rows, fill = TRUE)
  fwrite(summary_tbl,
         file.path(TAB, "NDUFB3_retinal_eQTL_screening_summary.csv"))
  log_info(sprintf("Summary saved: %d rows", nrow(summary_tbl)))
  print(summary_tbl)
} else {
  log_warn("No summary rows could be built")
}

log_info("=== 07_7g v2.1 done ===")
log_info(sprintf("Output dir: %s", OUT))

