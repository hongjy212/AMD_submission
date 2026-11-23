# ============================================================
# Purpose: LD-clumped MR for OXPHOS genes (PLINK + 1kg.v3 EUR)
# Outcomes: finn-b-H7_AMD (all-stage) + ebi-a-GCST010723 (early)
# Method: ieugwasr::ld_clump() -> TwoSampleMR::harmonise_data() + mr()
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ieugwasr)
  library(TwoSampleMR)
})

ROOT <- "E:/retina_tools/aging_rgc_axon/AMD"
OUT  <- file.path(ROOT, "results", "07_7b_MR_plink")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

PLINK_BIN <- "D:/tools/plink_win64_20260923/plink.exe"
BFILE     <- "E:/retina_tools/pubdata/eQTL/LD_reference/1kg.v3/EUR"
if (!file.exists(PLINK_BIN)) stop("PLINK 未找到: ", PLINK_BIN)
if (!file.exists(paste0(BFILE, ".bed"))) stop("LD 参考面板未找到: ", BFILE)

# ============================================================
# ============================================================
EQTLGEN_TAB <- file.path(ROOT, "results", "07_7_eqtlgen_cis", "tables")
eqtl_all <- fread(file.path(EQTLGEN_TAB, "OXPHOS13_all_cis.csv"))
cat("Total cis-eQTL rows:", nrow(eqtl_all), "\n")

setnames(eqtl_all,
         old = c("AssessedAllele", "OtherAllele"),
         new = c("effect_allele", "other_allele"),
         skip_absent = TRUE)

eqtl_candidates <- eqtl_all[Pvalue < 1e-5]
cat("After P<1e-5:", nrow(eqtl_candidates), "\n")

cat("\n--- LD clumping per gene (main analysis: r2<0.001) ---\n")
iv_list <- list()
for (g in unique(eqtl_candidates$Gene)) {
  sub <- eqtl_candidates[Gene == g, ]
  g_sym <- unique(sub$symbol)[1]
  cat(sprintf("[%s] %s: %d candidate SNPs\n", g_sym, g, nrow(sub)))
  if (nrow(sub) < 2) { iv_list[[g]] <- sub; next }

  clumped <- tryCatch({
    ieugwasr::ld_clump(
      dat = dplyr::tibble(rsid = sub$SNP, pval = sub$Pvalue),
      plink_bin = PLINK_BIN, bfile = BFILE,
      clump_kb = 10000, clump_r2 = 0.001, clump_p = 1e-5
    )
  }, error = function(e) { cat("  Clumping failed:", conditionMessage(e), "\n"); NULL })

  if (is.null(clumped) || nrow(clumped) == 0) {
    iv_list[[g]] <- sub[0, ]
  } else {
    iv_list[[g]] <- sub[SNP %in% clumped$rsid, ]
    cat(sprintf("  After clumping: %d SNPs\n", nrow(iv_list[[g]])))
  }
}

eqtl_iv <- rbindlist(iv_list)
if (nrow(eqtl_iv) == 0) stop("No IVs survived clumping")

eqtl_iv[, Zscore := as.numeric(Zscore)]
eqtl_iv[, NrSamples := as.numeric(NrSamples)]
eqtl_iv[, beta := Zscore / sqrt(NrSamples)]
eqtl_iv[, se   := 1 / sqrt(NrSamples)]
eqtl_iv[, F_stat := Zscore^2]

cat("\nTotal IVs after clumping:", nrow(eqtl_iv), "\n")

read_outcome_vcf <- function(vcf_file) {
  cat("\nReading VCF:", vcf_file, "\n")
  vcf <- fread(vcf_file, sep = "\t", header = TRUE, skip = "#CHROM",
               stringsAsFactors = FALSE)
  fmt_col    <- colnames(vcf)[9]
  sample_col <- colnames(vcf)[10]
  fmt_keys   <- strsplit(vcf[[fmt_col]][1], ":", fixed = TRUE)[[1]]
  sample_vals <- tstrsplit(vcf[[sample_col]], ":", fixed = TRUE)
  names(sample_vals) <- fmt_keys

  get_field <- function(keys, candidates) {
    hit <- intersect(candidates, keys); if (length(hit) == 0) return(NULL); hit[1]
  }
  es_key <- get_field(fmt_keys, c("ES", "BETA", "Effect"))
  se_key <- get_field(fmt_keys, c("SE", "SEBETA", "StdErr"))
  lp_key <- get_field(fmt_keys, c("LP", "MLOG10P", "P"))

  outcome_sub <- data.table(
    SNP      = vcf$ID,
    EA_out   = vcf$ALT,
    OA_out   = vcf$REF,
    beta_out = as.numeric(sample_vals[[es_key]]),
    se_out   = as.numeric(sample_vals[[se_key]]),
    LP       = if (!is.null(lp_key)) as.numeric(sample_vals[[lp_key]]) else NA_real_
  )
  outcome_sub[, pval_out := ifelse(!is.na(LP) & LP > 1, 10^(-LP), LP)]
  outcome_sub <- outcome_sub[!is.na(beta_out) & !is.na(se_out) & SNP != "."]
  cat("Outcome SNPs:", nrow(outcome_sub), "\n")
  rm(vcf, sample_vals); gc()
  return(outcome_sub)
}

# ============================================================
# ============================================================
run_mr <- function(eqtl_iv_sub, outcome_sub, label, write_output = TRUE) {
  results_list <- list()

  for (g in unique(eqtl_iv_sub$Gene)) {
    iv_sub <- eqtl_iv_sub[Gene == g, ]
    g_sym  <- unique(iv_sub$symbol)[1]

    merged <- merge(iv_sub, outcome_sub, by = "SNP", allow.cartesian = TRUE)
    cat(sprintf("\n[%s|%s] Matched SNPs: %d\n", label, g_sym, nrow(merged)))
    if (nrow(merged) == 0) next

    if (!"NrSamples" %in% names(merged)) {
      merged[, NrSamples := (Zscore / beta)^2]
    }

    exp_dat <- data.frame(
      SNP                    = merged$SNP,
      beta.exposure          = merged$beta,
      se.exposure            = merged$se,
      pval.exposure          = merged$Pvalue,
      effect_allele.exposure = merged$effect_allele,
      other_allele.exposure  = merged$other_allele,
      eaf.exposure           = NA_real_,
      samplesize.exposure    = merged$NrSamples,
      exposure               = paste0("eQTL_", g_sym),
      id.exposure            = paste0("eQTL_", g_sym),
      stringsAsFactors       = FALSE
    )

    out_dat <- data.frame(
      SNP                   = merged$SNP,
      beta.outcome          = merged$beta_out,
      se.outcome            = merged$se_out,
      effect_allele.outcome = merged$EA_out,
      other_allele.outcome  = merged$OA_out,
      eaf.outcome           = NA_real_,
      pval.outcome          = merged$pval_out,
      outcome               = "AMD_outcome",
      id.outcome            = "AMD_outcome",
      samplesize.outcome    = NA,
      stringsAsFactors      = FALSE
    )

    dat <- tryCatch(
      suppressMessages(harmonise_data(exp_dat, out_dat, action = 2)),
      error = function(e) NULL
    )
    if (is.null(dat) || nrow(dat) == 0) next
    dat <- dat[dat$mr_keep, ]
    if (nrow(dat) == 0) next

    method_list <- if (nrow(dat) == 1) "mr_wald_ratio" else "mr_ivw"

    res <- tryCatch(
      suppressMessages(mr(dat, method_list = method_list)),
      error = function(e) NULL
    )
    if (is.null(res) || nrow(res) == 0) next

    results_list[[length(results_list) + 1]] <- data.table(
      gene      = g,
      symbol    = g_sym,
      exposure  = paste0("eQTL_", g_sym),
      n_iv      = res$nsnp[1],
      method    = ifelse(res$method[1] == "Inverse variance weighted",
                         "IVW", "Wald ratio"),
      beta_MR   = res$b[1],
      se_MR     = res$se[1],
      p_MR      = res$pval[1],
      OR_MR     = exp(res$b[1]),
      OR_lower  = exp(res$b[1] - 1.96 * res$se[1]),
      OR_upper  = exp(res$b[1] + 1.96 * res$se[1])
    )
  }

  if (length(results_list) == 0) return(NULL)
  results <- rbindlist(results_list)
  results[, p_FDR := p.adjust(p_MR, method = "BH")]
  setorder(results, p_MR)

  if (write_output) {
    fwrite(results, file.path(OUT, paste0("MR_", label, ".csv")))
  }
  return(results)
}

outcome_all   <- read_outcome_vcf("E:/retina_tools/pubdata/eQTL/gwas/finn-b-H7_AMD.vcf.gz")
outcome_early <- read_outcome_vcf("E:/retina_tools/pubdata/eQTL/gwas/ebi-a-GCST010723.vcf.gz")

res_all   <- run_mr(eqtl_iv, outcome_all,   "finnbH7_allstage")
res_early <- run_mr(eqtl_iv, outcome_early, "earlyAMD")

if (!is.null(res_all) && !is.null(res_early)) {
  cmp <- merge(
    res_all[,   .(gene, symbol, n_iv_all = n_iv, method_all = method,
                  OR_all = OR_MR, p_all = p_MR, FDR_all = p_FDR)],
    res_early[, .(gene, n_iv_early = n_iv, method_early = method,
                  OR_early = OR_MR, p_early = p_MR, FDR_early = p_FDR)],
    by = "gene", all = TRUE)
  fwrite(cmp, file.path(OUT, "MR_cross_cohort_comparison.csv"))
}

# ============================================================
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("Part 1.5: Build NDUFB3 main IV set (p<1e-5, r2<0.01)\n")
cat(strrep("=", 70), "\n", sep = "")

ndufb3_cis_all <- eqtl_all[symbol == "NDUFB3", ]
cat("NDUFB3 all cis rows:", nrow(ndufb3_cis_all), "\n")

cand_main <- ndufb3_cis_all[Pvalue < 1e-5, ]
cand_main[, Zscore := as.numeric(Zscore)]
cand_main[, NrSamples := as.numeric(NrSamples)]
cand_main[, beta := Zscore / sqrt(NrSamples)]
cand_main[, se   := 1 / sqrt(NrSamples)]
cand_main[, F_stat := Zscore^2]

clumped_main <- ieugwasr::ld_clump(
  dat = dplyr::tibble(rsid = cand_main$SNP, pval = cand_main$Pvalue),
  plink_bin = PLINK_BIN, bfile = BFILE,
  clump_kb = 10000, clump_r2 = 0.01, clump_p = 1e-5
)
ndufb3_iv_main <- cand_main[SNP %in% clumped_main$rsid, ]
cat("NDUFB3 main IV set:", nrow(ndufb3_iv_main), "SNPs\n")
print(ndufb3_iv_main[, .(SNP, Pvalue, Zscore, effect_allele, other_allele,
                          NrSamples, beta)])

# ============================================================
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("Part 2: NDUFB3 IV details + Steiger\n")
cat(strrep("=", 70), "\n", sep = "")

ndufb3_iv <- copy(ndufb3_iv_main)
N_OUT_AMD <- 105248L   # GCST010723: 14,034 cases + 91,214 controls

ndufb3_iv[, r2_exp := Zscore^2 / NrSamples]
ndufb3_iv[, r2_out_approx := {
  m <- merge(.SD, outcome_early[, .(SNP, beta_out, se_out)],
             by = "SNP", all.x = TRUE)
  m$beta_out^2 / (m$se_out^2 * N_OUT_AMD)
}]
ndufb3_iv[, steiger_correct := r2_exp > r2_out_approx]

iv_details <- ndufb3_iv[, .(
  SNP, SNPChr, SNPPos, Pvalue, Zscore,
  NrSamples,
  effect_allele, other_allele,
  in_known_locus = FALSE, known_locus_name = NA_character_
)]
fwrite(iv_details, file.path(OUT, "NDUFB3_IV_details.csv"))
cat("Saved: NDUFB3_IV_details.csv (", nrow(iv_details), " SNPs)\n")
print(iv_details)

steiger_out <- ndufb3_iv[, .(
  SNP,
  r2_exp = round(r2_exp, 6),
  r2_out = round(r2_out_approx, 6),
  steiger_correct_direction = steiger_correct
)]
fwrite(steiger_out, file.path(OUT, "NDUFB3_steiger_check.csv"))
cat("Saved: NDUFB3_steiger_check.csv\n")
print(steiger_out)

# ============================================================
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("Part 3: NDUFB3 sensitivity grid\n")
cat(strrep("=", 70), "\n", sep = "")

p_thrs  <- c(1e-8, 1e-6, 1e-5, 1e-4)
r2_thrs <- c(0.001, 0.01, 0.05, 0.1)

grid_rows <- list()

for (pt in p_thrs) {
  cand <- copy(ndufb3_cis_all[Pvalue < pt, ])
  n_cand <- nrow(cand)

  if (n_cand == 0) {
    for (r2t in r2_thrs) {
      for (ot in c("finnbH7_allstage", "earlyAMD")) {
        grid_rows[[length(grid_rows) + 1]] <- data.table(
          outcome = ot, p_thr = pt, r2_thr = r2t,
          n_candidates = 0L, n_iv = 0L, method = NA_character_,
          beta_MR = NA_real_, se_MR = NA_real_, p_MR = NA_real_,
          OR_MR = NA_real_, OR_lower = NA_real_, OR_upper = NA_real_,
          F_stat_max = NA_real_,
          p_thr_label = format(pt, scientific = TRUE))
      }
    }
    next
  }

  cand[, Zscore := as.numeric(Zscore)]
  cand[, NrSamples := as.numeric(NrSamples)]
  cand[, beta := Zscore / sqrt(NrSamples)]
  cand[, se   := 1 / sqrt(NrSamples)]
  cand[, F_stat := Zscore^2]

  for (r2t in r2_thrs) {
    if (n_cand < 2) {
      iv_sub <- cand
    } else {
      clumped <- tryCatch({
        ieugwasr::ld_clump(
          dat = dplyr::tibble(rsid = cand$SNP, pval = cand$Pvalue),
          plink_bin = PLINK_BIN, bfile = BFILE,
          clump_kb = 10000, clump_r2 = r2t, clump_p = pt
        )
      }, error = function(e) NULL)
      if (is.null(clumped) || nrow(clumped) == 0) {
        iv_sub <- cand[0, ]
      } else {
        iv_sub <- cand[SNP %in% clumped$rsid, ]
      }
    }

    for (ot in c("finnbH7_allstage", "earlyAMD")) {
      outc <- if (ot == "finnbH7_allstage") outcome_all else outcome_early

      if (nrow(iv_sub) == 0) {
        grid_rows[[length(grid_rows) + 1]] <- data.table(
          outcome = ot, p_thr = pt, r2_thr = r2t,
          n_candidates = n_cand, n_iv = 0L, method = NA_character_,
          beta_MR = NA_real_, se_MR = NA_real_, p_MR = NA_real_,
          OR_MR = NA_real_, OR_lower = NA_real_, OR_upper = NA_real_,
          F_stat_max = NA_real_,
          p_thr_label = format(pt, scientific = TRUE))
        next
      }

      mr_res <- tryCatch(
        run_mr(iv_sub, outc,
               paste0("NDUFB3_", ot, "_", pt, "_", r2t),
               write_output = FALSE),
        error = function(e) NULL)

      if (is.null(mr_res) || nrow(mr_res) == 0) {
        grid_rows[[length(grid_rows) + 1]] <- data.table(
          outcome = ot, p_thr = pt, r2_thr = r2t,
          n_candidates = n_cand, n_iv = nrow(iv_sub), method = "no_match",
          beta_MR = NA_real_, se_MR = NA_real_, p_MR = NA_real_,
          OR_MR = NA_real_, OR_lower = NA_real_, OR_upper = NA_real_,
          F_stat_max = max(iv_sub$F_stat, na.rm = TRUE),
          p_thr_label = format(pt, scientific = TRUE))
        next
      }

      grid_rows[[length(grid_rows) + 1]] <- data.table(
        outcome = ot, p_thr = pt, r2_thr = r2t,
        n_candidates = n_cand, n_iv = mr_res$n_iv[1], method = mr_res$method[1],
        beta_MR = mr_res$beta_MR[1], se_MR = mr_res$se_MR[1],
        p_MR = mr_res$p_MR[1], OR_MR = mr_res$OR_MR[1],
        OR_lower = mr_res$OR_lower[1], OR_upper = mr_res$OR_upper[1],
        F_stat_max = max(iv_sub$F_stat, na.rm = TRUE),
        p_thr_label = format(pt, scientific = TRUE))
    }
  }
}

grid <- rbindlist(grid_rows, fill = TRUE)
grid[, FDR := p.adjust(p_MR, method = "BH"), by = outcome]
fwrite(grid, file.path(OUT, "NDUFB3_sensitivity_grid.csv"))
cat("Saved: NDUFB3_sensitivity_grid.csv (", nrow(grid), " rows)\n")
print(grid[, .(outcome, p_thr, r2_thr, n_candidates, n_iv, method,
               OR_MR, p_MR, FDR)])

# ============================================================
# Part 4：NDUFB3 FinnGen IVW + leave-one-out
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("Part 4: NDUFB3 FinnGen IVW + leave-one-out\n")
cat(strrep("=", 70), "\n", sep = "")

# --- FinnGen IVW ---
ndufb3_mr_finn <- run_mr(ndufb3_iv_main, outcome_all, "NDUFB3_FinnGen_IVW")
if (!is.null(ndufb3_mr_finn)) {
  finn_ivw <- ndufb3_mr_finn[, .(outcome = "FinnGen_allstage_AMD", n_iv, method,
                                  beta_MR, se_MR, OR_MR, OR_lower, OR_upper, p_MR)]
  fwrite(finn_ivw, file.path(OUT, "NDUFB3_FinnGen_IVW.csv"))
  cat("Saved: NDUFB3_FinnGen_IVW.csv\n")
  print(finn_ivw)
}

# --- LOO（TwoSampleMR::mr_leaveoneout）---
if (nrow(ndufb3_iv_main) >= 2) {
  exp_main <- data.frame(
    SNP                    = ndufb3_iv_main$SNP,
    beta.exposure          = ndufb3_iv_main$beta,
    se.exposure            = ndufb3_iv_main$se,
    pval.exposure          = ndufb3_iv_main$Pvalue,
    effect_allele.exposure = ndufb3_iv_main$effect_allele,
    other_allele.exposure  = ndufb3_iv_main$other_allele,
    eaf.exposure           = NA_real_,
    samplesize.exposure    = ndufb3_iv_main$NrSamples,
    exposure               = "NDUFB3_blood",
    id.exposure            = "NDUFB3_blood",
    stringsAsFactors       = FALSE
  )
  merged_finn <- merge(ndufb3_iv_main, outcome_all, by = "SNP", allow.cartesian = TRUE)
  out_finn <- data.frame(
    SNP                   = merged_finn$SNP,
    beta.outcome          = merged_finn$beta_out,
    se.outcome            = merged_finn$se_out,
    effect_allele.outcome = merged_finn$EA_out,
    other_allele.outcome  = merged_finn$OA_out,
    eaf.outcome           = NA_real_,
    pval.outcome          = merged_finn$pval_out,
    outcome               = "FinnGen_allstage_AMD",
    id.outcome            = "finn-b-H7_AMD",
    samplesize.outcome    = NA,
    stringsAsFactors      = FALSE
  )
  dat_finn <- suppressMessages(harmonise_data(exp_main, out_finn, action = 2))
  dat_finn <- dat_finn[dat_finn$mr_keep, ]
  if (nrow(dat_finn) >= 2) {
    loo_finn <- as.data.table(mr_leaveoneout(dat_finn))   # ← 修复点
    loo_out <- loo_finn[, .(
      excluded_SNP = SNP,
      OR_loo   = exp(b),
      OR_lower = exp(b - 1.96 * se),
      OR_upper = exp(b + 1.96 * se),
      p_loo    = p
    )]
    fwrite(loo_out, file.path(OUT, "NDUFB3_FinnGen_LOO.csv"))
    cat("Saved: NDUFB3_FinnGen_LOO.csv\n")
    print(loo_out)
  }
} else {
  cat("[SKIP] FinnGen 匹配 SNP < 2，无法做 LOO\n")
}

cat("\nResults saved to:", OUT, "\n")
cat("=== 07_7b_MR_plink done ===\n")

