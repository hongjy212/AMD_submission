# ============================================================
# 08_3_mr.r
# Output:  results/08_3_mr/
# ============================================================

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(TwoSampleMR)
})

ROOT     <- "E:/retina_tools/aging_rgc_axon/AMD"
DIR_THIS <- file.path(ROOT, "results", "08_3_mr")
DIR_MR   <- file.path(DIR_THIS, "01_mr")
DIR_ROB  <- file.path(DIR_THIS, "02_robust")
for (d in c(DIR_MR, DIR_ROB)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

IN_IV      <- file.path(ROOT, "results/08_2_eqtl_iv/02_iv/iv_all.csv")
GWAS_ALL   <- "E:/retina_tools/pubdata/eQTL/gwas/finn-b-H7_AMD.vcf.gz"
GWAS_EARLY <- "E:/retina_tools/pubdata/eQTL/gwas/ebi-a-GCST010723.vcf.gz"

log_info("=== 08_3 mr ===")
if (!file.exists(IN_IV)) stop("Run 08_2 first: ", IN_IV)
iv <- fread(IN_IV)
log_info("IVs: {nrow(iv)}")
log_info("IV columns: {paste(names(iv), collapse=', ')}")

read_outcome_vcf <- function(vcf_file) {
  con <- gzfile(vcf_file, "r"); n_header <- 0
  repeat {
    l <- readLines(con, n = 1, warn = FALSE)
    if (length(l) == 0 || !grepl("^#", l)) break
    n_header <- n_header + 1
  }
  close(con)
  vcf <- fread(vcf_file, skip = n_header, header = FALSE, sep = "\t",
               col.names = c("CHROM","POS","ID","REF","ALT","QUAL","FILTER",
                             "INFO","FORMAT","SAMPLE"), showProgress = FALSE)
  parse_field <- function(x, idx)
    sapply(strsplit(x, ":", fixed = TRUE),
           function(p) suppressWarnings(as.numeric(p[idx])))
  data.table(SNP = vcf$ID, EA_out = vcf$ALT, OA_out = vcf$REF,
             beta_out = parse_field(vcf$SAMPLE, 1),
             se_out   = parse_field(vcf$SAMPLE, 2),
             pval_out = 10^(-parse_field(vcf$SAMPLE, 3)))
}

build_dat <- function(iv_sub, outcome, label) {
  iv_sub <- copy(iv_sub)

  if ("AssessedAllele" %in% names(iv_sub))
    setnames(iv_sub, "AssessedAllele", "effect_allele")
  if ("OtherAllele" %in% names(iv_sub))
    setnames(iv_sub, "OtherAllele", "other_allele")

  merged <- merge(iv_sub, outcome, by = "SNP", allow.cartesian = TRUE)
  if (nrow(merged) == 0) return(NULL)
  g_sym <- unique(iv_sub$symbol)[1]

  need <- c("beta", "se", "effect_allele", "other_allele",
            "beta_out", "se_out", "EA_out", "OA_out", "pval_out")
  miss <- setdiff(need, names(merged))
  if (length(miss) > 0) {
    log_warn("build_dat [{g_sym}]: missing cols: {paste(miss, collapse=', ')}")
    return(NULL)
  }

  exp_dat <- data.frame(
    SNP = merged$SNP,
    beta.exposure = merged$beta,
    se.exposure   = merged$se,
    pval.exposure = merged$Pvalue,
    effect_allele.exposure = merged$effect_allele,
    other_allele.exposure  = merged$other_allele,
    eaf.exposure = NA_real_,
    samplesize.exposure = merged$NrSamples,
    exposure = paste0("eQTL_", g_sym),
    id.exposure = paste0("eQTL_", g_sym),
    stringsAsFactors = FALSE)

  out_dat <- data.frame(
    SNP = merged$SNP,
    beta.outcome = merged$beta_out,
    se.outcome   = merged$se_out,
    effect_allele.outcome = merged$EA_out,
    other_allele.outcome  = merged$OA_out,
    eaf.outcome = NA_real_,
    pval.outcome = merged$pval_out,
    outcome = label, id.outcome = label,
    samplesize.outcome = NA, stringsAsFactors = FALSE)

  dat <- tryCatch(suppressMessages(harmonise_data(exp_dat, out_dat, action = 2)),
                  error = function(e) NULL)
  if (is.null(dat) || nrow(dat) == 0) return(NULL)
  dat[dat$mr_keep, ]
}

log_info("--- main MR ---")
out_all   <- read_outcome_vcf(GWAS_ALL)
out_early <- read_outcome_vcf(GWAS_EARLY)

run_main_mr <- function(iv, outcome, label) {
  results_list <- list()
  for (g in unique(iv$Gene)) {
    iv_sub <- iv[Gene == g, ]
    dat <- build_dat(iv_sub, outcome, label)
    if (is.null(dat) || nrow(dat) == 0) next
    methods <- if (nrow(dat) == 1) "mr_wald_ratio" else "mr_ivw"
    res <- tryCatch(suppressMessages(mr(dat, method_list = methods)),
                    error = function(e) NULL)
    if (is.null(res) || nrow(res) == 0) next
    results_list[[length(results_list) + 1]] <- data.table(
      gene = g, symbol = unique(iv_sub$symbol)[1],
      gene_type = unique(iv_sub$gene_type)[1],
      outcome = label, n_iv = res$nsnp[1],
      method = ifelse(res$method[1] == "Inverse variance weighted",
                      "IVW", "Wald ratio"),
      beta_MR = res$b[1], se_MR = res$se[1], p_MR = res$pval[1],
      OR_MR = exp(res$b[1]),
      OR_lower = exp(res$b[1] - 1.96 * res$se[1]),
      OR_upper = exp(res$b[1] + 1.96 * res$se[1]))
  }
  if (length(results_list) == 0) return(NULL)
  rbindlist(results_list)
}

res_all   <- run_main_mr(iv, out_all,   "allstage")
res_early <- run_main_mr(iv, out_early, "early")

if (exists("res_all") && !is.null(res_all))
  fwrite(res_all, file.path(DIR_MR, "MR_allstage.csv"))
if (exists("res_early") && !is.null(res_early))
  fwrite(res_early, file.path(DIR_MR, "MR_early.csv"))

log_info("--- robust MR ---")
if (exists("res_early") && !is.null(res_early)) {
  genes_robust <- res_early[n_iv >= 3, unique(gene)]
  log_info("Genes with n_iv>=3: {length(genes_robust)}")

  robust_list <- list(); het_list <- list()
  pleio_list  <- list(); loo_list <- list(); single_list <- list()

  for (g in genes_robust) {
    iv_sub <- iv[Gene == g, ]
    dat <- build_dat(iv_sub, out_early, "early")
    if (is.null(dat) || nrow(dat) < 3) next
    g_sym <- unique(iv_sub$symbol)[1]

    res <- tryCatch(suppressMessages(mr(dat, method_list = c(
      "mr_ivw","mr_egger_regression","mr_weighted_median","mr_weighted_mode"))),
      error = function(e) NULL)
    if (!is.null(res)) { res$symbol <- g_sym; robust_list[[g_sym]] <- as.data.table(res) }

    het <- tryCatch(suppressMessages(mr_heterogeneity(dat)), error = function(e) NULL)
    if (!is.null(het)) { het$symbol <- g_sym; het_list[[g_sym]] <- as.data.table(het) }

    ple <- tryCatch(suppressMessages(mr_pleiotropy_test(dat)), error = function(e) NULL)
    if (!is.null(ple)) { ple$symbol <- g_sym; pleio_list[[g_sym]] <- as.data.table(ple) }

    loo <- tryCatch(suppressMessages(mr_leaveoneout(dat)), error = function(e) NULL)
    if (!is.null(loo)) { loo$symbol <- g_sym; loo_list[[g_sym]] <- as.data.table(loo) }

    sng <- tryCatch(suppressMessages(mr_singlesnp(dat)), error = function(e) NULL)
    if (!is.null(sng)) { sng$symbol <- g_sym; single_list[[g_sym]] <- as.data.table(sng) }
  }

  if (length(robust_list) > 0)  fwrite(rbindlist(robust_list, fill = TRUE),  file.path(DIR_ROB, "robust_methods.csv"))
  if (length(het_list) > 0)     fwrite(rbindlist(het_list, fill = TRUE),     file.path(DIR_ROB, "heterogeneity.csv"))
  if (length(pleio_list) > 0)   fwrite(rbindlist(pleio_list, fill = TRUE),   file.path(DIR_ROB, "pleiotropy.csv"))
  if (length(loo_list) > 0)     fwrite(rbindlist(loo_list, fill = TRUE),     file.path(DIR_ROB, "leaveoneout.csv"))
  if (length(single_list) > 0)  fwrite(rbindlist(single_list, fill = TRUE),  file.path(DIR_ROB, "single_snp.csv"))
  log_info("Robust genes: {length(robust_list)}")
}

log_info("=== 08_3 done ===")

