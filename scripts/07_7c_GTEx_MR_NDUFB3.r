# ============================================================
# 07_7c_GTEx_MR_NDUFB3.r  (v5)
# Output: results/07_7c_GTEx_MR_NDUFB3/
# ============================================================

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(httr2); library(jsonlite); library(data.table)
  library(dplyr); library(TwoSampleMR)
})

ROOT     <- "E:/retina_tools/aging_rgc_axon/AMD"
OUT      <- file.path(ROOT, "results", "07_7c_GTEx_MR_NDUFB3")
DATA_DIR <- file.path(OUT, "data")
TAB_DIR  <- file.path(OUT, "tables")
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TAB_DIR,  recursive = TRUE, showWarnings = FALSE)
log_info("=== 07_7c GTEx MR NDUFB3 (v5) ===")

# ============================================================
# ============================================================
get_gtex_v10_samplesize <- function(tissue_id,
                                     cache_dir = NULL,
                                     verbose = TRUE) {
  candidates <- c(
    if (!is.null(cache_dir)) file.path(cache_dir,
      "GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt"),
    "E:/retina_tools/pubdata/eQTL/GTEx/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt",
    "E:/retina_tools/pubdata/eQTL/GTEx/metadata/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt",
    "E:/retina_tools/pubdata/eQTL/eQTLGen/GTEx_Analysis_v10_Annotations_SampleAttributesDS.txt"
  )
  for (f in candidates) {
    if (file.exists(f)) {
      meta <- tryCatch(data.table::fread(f, sep = "\t"),
                       error = function(e) NULL)
      if (!is.null(meta) && "SMTSD" %in% names(meta)) {
        n <- sum(meta$SMTSD == tissue_id, na.rm = TRUE)
        if (n > 0) {
          if (verbose) log_info(sprintf(
            "GTEx v10 samplesize for %s (local metadata): n=%d",
            tissue_id, n))
          return(as.integer(n))
        }
      }
    }
  }
  if (verbose) log_warn(sprintf(
    "GTEx v10 samplesize for %s not found locally; setting NA.",
    tissue_id))
  return(NA_integer_)
}

# ============================================================
# ============================================================
normalize_gtex_v10_eqtl <- function(dt) {
  dt <- as.data.table(dt)

  if ("snpId"  %in% names(dt)) setnames(dt, "snpId",  "SNP")
  if ("pValue" %in% names(dt)) setnames(dt, "pValue", "pval")

  required <- c("SNP", "pval", "nes", "variantId", "tissueSiteDetailId")
  missing_cols <- setdiff(required, names(dt))
  if (length(missing_cols) > 0) {
    stop("GTEx API 缺少必要列: ", paste(missing_cols, collapse = ", "))
  }

  parts <- tstrsplit(dt$variantId, "_", fixed = TRUE)
  # parts[[1]] = chr, parts[[2]] = pos, parts[[3]] = ref, parts[[4]] = alt, parts[[5]] = b38/b37
  dt[, ref_allele := parts[[3]]]
  dt[, alt_allele := parts[[4]]]
  dt[, effect_allele := alt_allele]
  dt[, other_allele  := ref_allele]
  dt[, c("ref_allele", "alt_allele") := NULL]

  # |z| = qnorm(pval/2, lower.tail = FALSE)
  # beta = nes
  # se   = |nes| / |z|
  dt[, beta := as.numeric(nes)]
  dt[, abs_z := qnorm(pval / 2, lower.tail = FALSE)]
  dt[, se := abs(beta) / abs_z]
  dt[, abs_z := NULL]

  keep_cols <- intersect(
    c("SNP", "beta", "se", "pval", "tissueSiteDetailId",
      "effect_allele", "other_allele",
      "geneSymbol", "gencodeId", "datasetId"),
    names(dt))
  dt <- dt[, ..keep_cols]

  return(dt)
}

# ============================================================
# ============================================================
gtex_file       <- file.path(DATA_DIR, "NDUFB3_eQTL_GTEx_Adipose.csv")
gtex_brain_file <- file.path(DATA_DIR, "NDUFB3_eQTL_GTEx_Brain.csv")

if (file.exists(gtex_file)) {
  log_info("Removing stale v4 output to regenerate with alleles...")
  file.remove(gtex_file)
}
if (file.exists(gtex_brain_file)) file.remove(gtex_brain_file)

log_info("Fetching GTEx v10 eQTL from API ...")
gencode_id <- "ENSG00000119013.10"

resp <- request("https://gtexportal.org/api/v2/association/singleTissueEqtl") |>
  req_url_query(
    gencodeId    = gencode_id,
    datasetId    = "gtex_v10",
    itemsPerPage = 1000
  ) |>
  req_perform()
result   <- resp_body_json(resp, simplifyVector = TRUE)
eqtl_raw <- as.data.table(result$data)

log_info(sprintf("API returned %d rows, columns: %s",
                 nrow(eqtl_raw),
                 paste(names(eqtl_raw), collapse = ", ")))

eqtl_all <- normalize_gtex_v10_eqtl(eqtl_raw)

ss_adipose <- get_gtex_v10_samplesize("Adipose_Subcutaneous",
                                       cache_dir = DATA_DIR)
ss_brain   <- get_gtex_v10_samplesize("Brain_Cortex",
                                       cache_dir = DATA_DIR)

# Adipose_Subcutaneous
adipose <- eqtl_all[tissueSiteDetailId == "Adipose_Subcutaneous"]
adipose[, samplesize := ss_adipose]
write.csv(adipose, gtex_file, row.names = FALSE)
log_info(sprintf("Saved GTEx adipose eQTL: %d SNPs, samplesize=%s",
                 nrow(adipose), as.character(ss_adipose)))

# Brain_Cortex
brain <- eqtl_all[tissueSiteDetailId == "Brain_Cortex"]
brain[, samplesize := ss_brain]
write.csv(brain, gtex_brain_file, row.names = FALSE)
log_info(sprintf("Saved GTEx brain eQTL: %d SNPs, samplesize=%s",
                 nrow(brain), as.character(ss_brain)))

# ============================================================
# ============================================================
adipose <- as.data.table(adipose)

if (!"pval" %in% names(adipose)) {
  stop("adipose 缺少 pval 列；列名：", paste(names(adipose), collapse = ", "))
}

top_snp <- adipose[order(pval)][1]
log_info(sprintf(
  "Top adipose SNP: %s (p = %.3e, beta = %.4f, se = %.4f, eff=%s, oth=%s)",
  top_snp$SNP, top_snp$pval, top_snp$beta, top_snp$se,
  top_snp$effect_allele, top_snp$other_allele))

# ============================================================
# 3. Read outcome
# ============================================================
VCF_EARLY <- "E:/retina_tools/pubdata/eQTL/gwas/ebi-a-GCST010723.vcf.gz"

read_outcome_vcf <- function(vcf_file, snp_ids) {
  con <- gzfile(vcf_file, "r")
  n_header <- 0
  repeat {
    l <- readLines(con, n = 1, warn = FALSE)
    if (length(l) == 0 || !grepl("^#", l)) break
    n_header <- n_header + 1
  }
  close(con)
  vcf <- fread(vcf_file, skip = n_header, header = FALSE, sep = "\t",
               col.names = c("CHROM","POS","ID","REF","ALT","QUAL","FILTER",
                             "INFO","FORMAT","SAMPLE"),
               showProgress = FALSE)
  vcf <- vcf[ID %in% snp_ids]
  if (nrow(vcf) == 0) return(NULL)
  parse_field <- function(x, idx) {
    sapply(strsplit(x, ":", fixed = TRUE),
           function(p) suppressWarnings(as.numeric(p[idx])))
  }
  data.table(
    SNP                   = vcf$ID,
    beta.outcome          = parse_field(vcf$SAMPLE, 1),
    se.outcome            = parse_field(vcf$SAMPLE, 2),
    effect_allele.outcome = vcf$ALT,
    other_allele.outcome  = vcf$REF,
    pval.outcome          = 10^(-parse_field(vcf$SAMPLE, 3))
  )
}

out_dat <- read_outcome_vcf(VCF_EARLY, top_snp$SNP)
if (is.null(out_dat) || nrow(out_dat) == 0) {
  stop("Top SNP not found in outcome VCF: ", top_snp$SNP)
}
log_info(sprintf("Outcome: %d SNP matched, eff=%s, oth=%s",
                 nrow(out_dat),
                 out_dat$effect_allele.outcome,
                 out_dat$other_allele.outcome))

# ============================================================
# 4. Wald ratio
# ============================================================
exp_dat <- data.frame(
  SNP                    = top_snp$SNP,
  beta.exposure          = as.numeric(top_snp$beta),
  se.exposure            = as.numeric(top_snp$se),
  pval.exposure          = as.numeric(top_snp$pval),
  effect_allele.exposure = top_snp$effect_allele,   # ★ v5 修复：从 variantId 解析
  other_allele.exposure  = top_snp$other_allele,    # ★ v5 修复
  eaf.exposure           = NA_real_,
  samplesize.exposure    = top_snp$samplesize,
  exposure               = "NDUFB3 (Adipose_Subcutaneous)",
  id.exposure            = "NDUFB3_adipose",
  stringsAsFactors       = FALSE
)
out_dat$outcome    <- "IAMDGC 2018 early AMD (GCST010723)"
out_dat$id.outcome <- "ebi-a-GCST010723"
out_dat$eaf.outcome <- NA_real_
out_dat$samplesize.outcome <- NA

dat <- suppressMessages(harmonise_data(exp_dat, out_dat, action = 2))
if (is.null(dat) || nrow(dat) == 0) {
  stop("Harmonisation failed: no SNP retained")
}
dat <- dat[dat$mr_keep, ]
if (nrow(dat) == 0) stop("Harmonisation failed: all SNPs removed")

log_info(sprintf("Harmonised: %d SNP, beta.exp=%.4f, beta.out=%.4f",
                 nrow(dat), dat$beta.exposure, dat$beta.outcome))

dat <- dat[1, ]

res <- suppressMessages(mr(dat, method_list = "mr_wald_ratio"))
if (is.na(res$b[1])) stop("Wald ratio failed")

log_info(sprintf("Wald ratio: OR = %.3f (%.3f-%.3f), p = %.3e",
                 exp(res$b), exp(res$b - 1.96*res$se),
                 exp(res$b + 1.96*res$se), res$pval))

# ============================================================
# 5. Output
# ============================================================
mr_out <- data.table(
  id.exposure     = res$id.exposure,
  id.outcome      = res$id.outcome,
  outcome         = out_dat$outcome[1],
  exposure        = exp_dat$exposure[1],
  method          = res$method,
  nsnp            = res$nsnp,
  b               = res$b,
  se              = res$se,
  pval            = res$pval,
  lo_ci           = res$b - 1.96 * res$se,
  up_ci           = res$b + 1.96 * res$se,
  or              = exp(res$b),
  or_lci95        = exp(res$b - 1.96 * res$se),
  or_uci95        = exp(res$b + 1.96 * res$se)
)
fwrite(mr_out, file.path(TAB_DIR, "MR_results.csv"))
log_info("Saved: MR_results.csv")
log_info("=== 07_7c done ===")

