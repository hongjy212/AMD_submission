# ============================================================
# 07_7d_MR_robust_sensitivity.r (v3)
# Input：results/07_7b_MR_plink/NDUFB3_IV_details.csv
#       results/07_7_eqtlgen_cis/tables/OXPHOS13_all_cis.csv
#       E:/retina_tools/pubdata/GWAS/AMD/ebi-a-GCST010723.vcf.gz
# Output：results/07_7d_MR_robust_sensitivity/tables/
#
# ============================================================

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
  library(dplyr); library(readr)
})

ROOT     <- "E:/retina_tools/aging_rgc_axon/AMD"
GWAS_DIR <- "E:/retina_tools/pubdata/GWAS/AMD"
OUT      <- file.path(ROOT, "results/07_7d_MR_robust_sensitivity")
TAB      <- file.path(OUT, "tables")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
set.seed(42)
log_info("=== 07_7d MR robust sensitivity (v3) ===")

# ============================================================
# ============================================================
iv <- read_csv(file.path(ROOT,
  "results/07_7b_MR_plink/NDUFB3_IV_details.csv"),
  show_col_types = FALSE)
log_info(sprintf("NDUFB3 IV: %d SNP", nrow(iv)))

cis_all <- fread(file.path(ROOT,
  "results/07_7_eqtlgen_cis/tables/OXPHOS13_all_cis.csv"))

stopifnot("SNP" %in% names(cis_all))
stopifnot("NrSamples" %in% names(cis_all))

cis_ndufb3 <- cis_all[symbol == "NDUFB3",
                      .(SNP, NrSamples_cis = as.numeric(NrSamples))]

if ("NrSamples" %in% names(iv)) {
  iv <- iv %>% select(-NrSamples)
}

iv <- merge(iv, cis_ndufb3, by = "SNP", all.x = TRUE)
iv <- iv %>% rename(NrSamples = NrSamples_cis)

if (any(is.na(iv$NrSamples))) {
  stop("以下 SNP 在 OXPHOS13_all_cis.csv 中找不到 NrSamples: ",
       paste(iv$SNP[is.na(iv$NrSamples)], collapse = ", "))
}
log_info(sprintf("NrSamples per SNP: %s",
                 paste(sprintf("%s=%d", iv$SNP, iv$NrSamples),
                       collapse = ", ")))
print(iv)

exp_dat <- data.frame(
  SNP                    = iv$SNP,
  beta.exposure          = iv$Zscore / sqrt(iv$NrSamples),
  se.exposure            = 1 / sqrt(iv$NrSamples),
  pval.exposure          = iv$Pvalue,
  effect_allele.exposure = iv$effect_allele,
  other_allele.exposure  = iv$other_allele,
  eaf.exposure           = NA_real_,
  samplesize.exposure    = iv$NrSamples,
  exposure               = "NDUFB3_blood",
  id.exposure            = "NDUFB3_blood",
  stringsAsFactors       = FALSE
)

# ============================================================
# ============================================================
vcf_file <- file.path(GWAS_DIR, "ebi-a-GCST010723.vcf.gz")

con <- gzfile(vcf_file, "r")
n_header <- 0
repeat {
  l <- readLines(con, n = 1, warn = FALSE)
  if (length(l) == 0 || !grepl("^#", l)) break
  n_header <- n_header + 1
}
close(con)
log_info(sprintf("VCF header 行数: %d", n_header))

log_info("读取 VCF（约 11M 变体，需几分钟）...")
vcf <- fread(vcf_file, skip = n_header, header = FALSE, sep = "\t",
             col.names = c("CHROM","POS","ID","REF","ALT","QUAL","FILTER",
                           "INFO","FORMAT","SAMPLE"),
             showProgress = TRUE)
log_info(sprintf("VCF 总变体: %d", nrow(vcf)))

vcf_ndufb3 <- vcf[ID %in% iv$SNP]
log_info(sprintf("命中 NDUFB3 SNP: %d / %d", nrow(vcf_ndufb3), nrow(iv)))

parse_field <- function(x, idx) {
  parts <- strsplit(x, ":", fixed = TRUE)
  sapply(parts, function(p) suppressWarnings(as.numeric(p[idx])))
}
vcf_ndufb3[, ES := parse_field(SAMPLE, 1)]
vcf_ndufb3[, SE := parse_field(SAMPLE, 2)]
vcf_ndufb3[, LP := parse_field(SAMPLE, 3)]
vcf_ndufb3[, AF := parse_field(SAMPLE, 4)]

out_dat <- data.frame(
  SNP                   = vcf_ndufb3$ID,
  beta.outcome          = vcf_ndufb3$ES,
  se.outcome            = vcf_ndufb3$SE,
  effect_allele.outcome = vcf_ndufb3$ALT,
  other_allele.outcome  = vcf_ndufb3$REF,
  eaf.outcome           = vcf_ndufb3$AF,
  pval.outcome          = 10^(-vcf_ndufb3$LP),
  outcome               = "AMD_early_GCST010723",
  id.outcome            = "ebi-a-GCST010723",
  stringsAsFactors      = FALSE
)
print(out_dat)

# ============================================================
# 3. Harmonise
# ============================================================
dat <- harmonise_data(exp_dat, out_dat, action = 2)
log_info(sprintf("Harmonised SNP: %d", nrow(dat)))
print(dat[, c("SNP","beta.exposure","beta.outcome",
              "effect_allele.exposure","effect_allele.outcome","mr_keep")])

# ============================================================
# 4. robust MR
# ============================================================
methods <- c("mr_ivw", "mr_wald_ratio")
if (nrow(dat) >= 3) {
  methods <- c(methods, "mr_egger_regression",
               "mr_weighted_median", "mr_weighted_mode")
}
res <- mr(dat, method_list = methods)
cat("\n=== MR 稳健方法 ===\n")
print(res[, c("method","nsnp","b","se","pval")])
write_csv(res, file.path(TAB, "mr_robust_earlyAMD.csv"))

# ============================================================
# ============================================================
if (nrow(dat) >= 3) {
  het <- mr_heterogeneity(dat)
  cat("\n=== 异质性 ===\n")
  print(het)
  het_ivw <- het[het$method == "Inverse variance weighted", ]
  if (nrow(het_ivw) > 0 && !is.na(het_ivw$Q)) {
    I2 <- max(0, (het_ivw$Q - het_ivw$Q_df) / het_ivw$Q * 100)
    log_info(sprintf("I2 = %.1f%%", I2))
    het$I2 <- ifelse(het$method == "Inverse variance weighted", I2, NA)
  }
  write_csv(het, file.path(TAB, "heterogeneity_earlyAMD.csv"))
}

# ============================================================
# ============================================================
if (nrow(dat) >= 3) {
  pleio <- mr_pleiotropy_test(dat)
  cat("\n=== Egger 截距（多效性） ===\n")
  print(pleio)
  write_csv(pleio, file.path(TAB, "pleiotropy_earlyAMD.csv"))
}

# ============================================================
# ============================================================
if (nrow(dat) >= 3) {
  loo <- mr_leaveoneout(dat)
  cat("\n=== 留一法 ===\n")
  print(loo[, c("SNP","b","se","p")])
  write_csv(loo, file.path(TAB, "leaveoneout_earlyAMD.csv"))
}

# ============================================================
# ============================================================
single <- mr_singlesnp(dat)
cat("\n=== 单 SNP ===\n")
print(single)
write_csv(single, file.path(TAB, "single_snp_earlyAMD.csv"))

# ============================================================
# 9. IV strength: F-statistic  ← v4 Fix
# ============================================================
log_info("--- IV strength: F-statistic ---")

iv_strength <- exp_dat %>%
  mutate(
    F_stat = (beta.exposure / se.exposure)^2,
    R2_approx = F_stat / (F_stat + samplesize.exposure - 2)
  ) %>%
  select(SNP,
         beta = beta.exposure,
         NrSamples = samplesize.exposure,
         F_stat,
         R2_approx)

write_csv(iv_strength, file.path(TAB, "NDUFB3_IV_F_stat.csv"))

cat("\n=== NDUFB3 IV strength ===\n")
print(iv_strength)

cat(sprintf("\nMean F: %.1f\n", mean(iv_strength$F_stat, na.rm = TRUE)))
cat(sprintf("Min F: %.1f\n", min(iv_strength$F_stat, na.rm = TRUE)))
cat(sprintf("Total R² (approx): %.4f (%.2f%%)\n",
            sum(iv_strength$R2_approx, na.rm = TRUE),
            100 * sum(iv_strength$R2_approx, na.rm = TRUE)))

if (min(iv_strength$F_stat, na.rm = TRUE) > 10) {
  log_info("All IVs have F > 10: not weak instruments")
} else {
  log_warn("Some IVs have F <= 10: potential weak instrument bias")
}

log_info("=== 07_7d v3 done ===")

