# ============================================================
# 07_7_extract_eqtlgen_cis.r
# Purpose: extract cis-eQTLs for 13 OXPHOS genes + complement
#          positive-control genes from eQTLGen, and prepare
#          instrument variables (IVs) for MR
#
# v2 update：
#
# Windows-compatible: uses fread() to read .gz directly
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

if (!requireNamespace("TwoSampleMR", quietly = TRUE)) {
  message("TwoSampleMR 未安装。如果后续要做 MR，请先安装：")
  message("install.packages('TwoSampleMR', repos = c('https://mrcieu.r-universe.dev', 'https://cloud.r-project.org'))")
}

ROOT <- "E:/retina_tools/aging_rgc_axon/AMD"
OUT  <- file.path(ROOT, "results", "07_7_eqtlgen_cis")
TAB  <- file.path(OUT, "tables")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)

EQTLGEN_FILE <- "E:/retina_tools/pubdata/eQTL/eQTLGen/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"
stopifnot(file.exists(EQTLGEN_FILE))

log_info <- function(...) cat(sprintf("[%s] [INFO] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)))
log_warn <- function(...) cat(sprintf("[%s] [WARN] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...)))

write_tab <- function(df, name) {
  p <- file.path(TAB, name)
  fwrite(df, p)
  log_info(sprintf("saved: %s (%d rows)", name, nrow(df)))
  invisible(p)
}

log_info("=== 07_7 extract eQTLGen cis-eQTLs (v2: OXPHOS + complement) ===")

# ============================================================
# ============================================================
oxphos13 <- data.frame(
  symbol = c("NDUFAB1", "GABARAPL2", "NDUFB2", "NDUFC1",
             "COX7A2", "SDHB", "NDUFB3", "COX7C", "UQCRQ",
             "UQCRFS1", "UQCRH", "NDUFB1", "NDUFA12"),
  ensembl = c("ENSG00000004779", "ENSG00000034713", "ENSG00000090266",
              "ENSG00000109390", "ENSG00000131055", "ENSG00000117118",
              "ENSG00000119013", "ENSG00000127184", "ENSG00000164442",
              "ENSG00000169021", "ENSG00000173660", "ENSG00000183648",
              "ENSG00000184752"),
  gene_type = "OXPHOS",
  stringsAsFactors = FALSE
)

complement_genes <- data.frame(
  symbol = c("CFH", "CD59", "CFHR5", "C3", "CFB"),
  ensembl = c("ENSG00000000971",   # CFH
              "ENSG00000085063",   # CD59
              "ENSG00000134365",   # CFHR5
              "ENSG00000125730",   # C3
              "ENSG00000243649"),  # CFB
  gene_type = "complement_positive_control",
  stringsAsFactors = FALSE
)

target_genes <- rbind(oxphos13, complement_genes)
ens_target <- target_genes$ensembl

log_info(sprintf("Target genes: %d (OXPHOS: %d, complement: %d)",
                 nrow(target_genes),
                 sum(target_genes$gene_type == "OXPHOS"),
                 sum(target_genes$gene_type == "complement_positive_control")))
print(target_genes)

# ============================================================
# ============================================================
log_info("--- 3. Reading eQTLGen directly with fread ---")
log_info("(File size ~3.6 GB compressed; may take 2-5 minutes)")

eqtl_all <- fread(
  EQTLGEN_FILE,
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE,
  showProgress = TRUE
)

log_info(sprintf("Total rows read: %d", nrow(eqtl_all)))
log_info(sprintf("Columns: %s", paste(colnames(eqtl_all), collapse = ", ")))

required_cols <- c("Pvalue", "SNP", "SNPChr", "SNPPos", "AssessedAllele",
                   "OtherAllele", "Zscore", "Gene", "GeneSymbol",
                   "GeneChr", "GenePos", "NrCohorts", "NrSamples", "FDR")
missing_cols <- setdiff(required_cols, colnames(eqtl_all))
if (length(missing_cols) > 0) {
  stop(sprintf("缺少列: %s", paste(missing_cols, collapse = ", ")))
}

eqtl <- eqtl_all[Gene %in% ens_target, ]
rm(eqtl_all); gc()

log_info(sprintf("After gene filter: %d rows", nrow(eqtl)))
if (nrow(eqtl) == 0) {
  stop("没有提取到目标基因的行。请检查 Ensembl ID 是否正确。")
}

eqtl <- merge(eqtl, target_genes[, c("symbol", "ensembl", "gene_type")],
              by.x = "Gene", by.y = "ensembl", all.x = TRUE)
log_info(sprintf("After merge: %d rows", nrow(eqtl)))

missing_genes <- setdiff(target_genes$symbol, unique(eqtl$symbol))
if (length(missing_genes) > 0) {
  log_warn(sprintf("以下目标基因在 eQTLGen 中无 cis-eQTL 行: %s",
                   paste(missing_genes, collapse = ", ")))
}

write_tab(eqtl, "OXPHOS13_all_cis.csv")

# ============================================================
# 4. cis filter（±1 Mb）
# ============================================================
log_info("--- 4. cis filter (±1 Mb) ---")
eqtl$SNPPos <- as.numeric(eqtl$SNPPos)
eqtl$GenePos <- as.numeric(eqtl$GenePos)
eqtl$cis_dist <- abs(eqtl$SNPPos - eqtl$GenePos)
eqtl_cis <- eqtl[eqtl$cis_dist <= 1e6, ]
log_info(sprintf("After cis filter: %d", nrow(eqtl_cis)))

# ============================================================
# ============================================================
log_info("--- 5. Significance filter (P < 5e-8) ---")
eqtl_cis$Pvalue <- as.numeric(eqtl_cis$Pvalue)
eqtl_sig <- eqtl_cis[eqtl_cis$Pvalue < 5e-8, ]
log_info(sprintf("After significance filter: %d", nrow(eqtl_sig)))

gene_counts <- eqtl_sig %>%
  group_by(gene_type, symbol, Gene) %>%
  summarise(n_snp = n(), .groups = "drop") %>%
  arrange(gene_type, desc(n_snp))
cat("\n=== Per-gene cis-eQTL counts (P<5e-8) ===\n")
print(as.data.frame(gene_counts), row.names = FALSE)

# ============================================================
# ============================================================
log_info("--- 6. Build instrument variables ---")

eqtl_iv <- eqtl_sig %>%
  group_by(Gene, symbol, gene_type) %>%
  slice_min(Pvalue, n = 1, with_ties = FALSE) %>%
  ungroup()

# library(TwoSampleMR)
# eqtl_iv <- list()
# for (g in unique(eqtl_sig$Gene)) {
#   sub <- eqtl_sig[eqtl_sig$Gene == g, ]
#   if (nrow(sub) < 2) { eqtl_iv[[g]] <- sub; next }
#   clumped <- tryCatch({
#     clump_data(
#       data.frame(SNP = sub$SNP, pval = sub$Pvalue),
#       clump_kb = 10000, clump_r2 = 0.001, clump_p1 = 5e-8, pop = "EUR"
#     )
#   }, error = function(e) NULL)
#   if (is.null(clumped) || nrow(clumped) == 0) {
#     eqtl_iv[[g]] <- sub[0, ]
#   } else {
#     eqtl_iv[[g]] <- sub[sub$SNP %in% clumped$SNP, ]
#   }
# }
# eqtl_iv <- do.call(rbind, eqtl_iv)
# rownames(eqtl_iv) <- NULL

log_info(sprintf("Total IVs: %d", nrow(eqtl_iv)))

# ============================================================
# ============================================================
log_info("--- 7. Compute beta / se from Z-score ---")
eqtl_iv$Zscore <- as.numeric(eqtl_iv$Zscore)
eqtl_iv$NrSamples <- as.numeric(eqtl_iv$NrSamples)

eqtl_iv$beta <- eqtl_iv$Zscore / sqrt(eqtl_iv$NrSamples)
eqtl_iv$se   <- 1 / sqrt(eqtl_iv$NrSamples)
eqtl_iv$F_stat <- eqtl_iv$Zscore^2

write_tab(eqtl_iv, "OXPHOS13_IVs.csv")

# ============================================================
# ============================================================
summary_df <- eqtl_iv %>%
  group_by(gene_type, symbol, Gene) %>%
  summarise(
    n_iv = n(),
    max_F = max(F_stat, na.rm = TRUE),
    min_p = min(Pvalue, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(gene_type, desc(n_iv))

write_tab(summary_df, "OXPHOS13_summary.csv")
cat("\n=== IV summary per gene ===\n")
print(as.data.frame(summary_df), row.names = FALSE)

# ============================================================
# ============================================================
log_info("--- 9. Format for TwoSampleMR ---")
mr_format <- data.frame(
  SNP = eqtl_iv$SNP,
  exposure = paste0("eQTL_", eqtl_iv$symbol),
  effect_allele.exposure = eqtl_iv$AssessedAllele,
  other_allele.exposure = eqtl_iv$OtherAllele,
  eaf.exposure = NA,
  beta.exposure = eqtl_iv$beta,
  se.exposure = eqtl_iv$se,
  pval.exposure = eqtl_iv$Pvalue,
  samplesize.exposure = eqtl_iv$NrSamples,
  id.exposure = eqtl_iv$Gene,
  gene_type = eqtl_iv$gene_type,
  stringsAsFactors = FALSE
)

write_tab(mr_format, "OXPHOS13_IVs_TwoSampleMR.csv")

# ============================================================
# ============================================================
log_info("--- 10. Positive control integrity check ---")

comp_symbols <- complement_genes$symbol
comp_summary <- summary_df[summary_df$symbol %in% comp_symbols, ]

cat("\n=== Complement positive-control genes in eQTLGen ===\n")
if (nrow(comp_summary) == 0) {
  log_warn("没有任何补体阳性对照基因在 eQTLGen 中检出 cis-eQTL。")
  log_warn("建议：检查 Ensembl ID 是否正确，或改用 GTEx cis-eQTL。")
} else {
  print(as.data.frame(comp_summary), row.names = FALSE)
}

for (g in comp_symbols) {
  n_iv_g <- sum(summary_df$symbol == g)
  if (n_iv_g == 0) {
    log_warn(sprintf("  [MISSING] %s: 无 cis-eQTL (P<5e-8)", g))
  } else {
    log_info(sprintf("  [OK] %s: %d IV(s), max F = %.1f",
                     g, n_iv_g,
                     summary_df$max_F[summary_df$symbol == g]))
  }
}

# ============================================================
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("=== 07_7 summary (v2: OXPHOS + complement) ===\n")
cat(strrep("=", 70), "\n", sep = "")
cat("\nRaw cis-eQTLs:", nrow(eqtl), "\n")
cat("After cis filter:", nrow(eqtl_cis), "\n")
cat("After P<5e-8:", nrow(eqtl_sig), "\n")
cat("After IV selection:", nrow(eqtl_iv), "\n")

cat("\n--- OXPHOS genes ---\n")
ox_summary <- summary_df[summary_df$gene_type == "OXPHOS", ]
print(as.data.frame(ox_summary), row.names = FALSE)

cat("\n--- Complement positive controls ---\n")
print(as.data.frame(comp_summary), row.names = FALSE)

cat("\nOutput dir:", OUT, "\n")
log_info("=== 07_7 done ===")

