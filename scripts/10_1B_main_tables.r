# 10_1B_main_tables.R (v7.1)
# Main tables 1-4

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressWarnings(suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(knitr)
  library(flextable); library(officer)
}))

results_dir <- DIR_RESULT
table_dir   <- file.path(results_dir, "10_1B_main_tables")
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)
log_info("=== 10_1B Main tables (v7.1) ===")

save_table <- function(df, name, caption = NULL) {
  if (is.null(df)) {
    log_warn(sprintf("skip table: %s (NULL)", name))
    return(invisible())
  }
  write_csv(df, file.path(table_dir, paste0(name, ".csv")), progress = FALSE)
  if (!is.null(caption)) {
    html <- kable(df, format = "html", caption = caption, row.names = FALSE)
    writeLines(as.character(html), file.path(table_dir, paste0(name, ".html")))
  }
  log_info(sprintf("saved table: %s (%d x %d)", name, nrow(df), ncol(df)))
}

read_required <- function(path) {
  if (!file.exists(path)) stop(sprintf("REQUIRED FILE MISSING: %s", path))
  read_csv(path, show_col_types = FALSE, progress = FALSE)
}

check_cols <- function(df, required, name) {
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    stop(sprintf("%s: missing columns: %s", name, paste(missing, collapse = ", ")))
  }
}

normalize_complex <- function(x) {
  x <- toupper(trimws(x))
  x <- gsub("^COMPLEX[ _]*", "", x)
  x <- gsub("^C[ _]*", "", x)
  paste0("Complex ", x)
}

fmt_num <- function(x, digits = 3) {
  if (is.numeric(x)) {
    ifelse(is.na(x), NA_character_,
           ifelse(abs(x) < 0.001 & x != 0,
                  formatC(x, format = "e", digits = 2),
                  formatC(x, format = "f", digits = digits)))
  } else x
}

make_ft <- function(df, font_size = 8) {
  flextable(df) %>%
    fontsize(size = font_size, part = "all") %>%
    font(fontname = "Times New Roman", part = "all") %>%
    align(align = "center", part = "all") %>%
    align(j = 1, align = "left", part = "all") %>%
    valign(valign = "center", part = "all") %>%
    bold(part = "header") %>%
    border_remove() %>%
    hline_top(part = "header", border = fp_border(width = 1)) %>%
    hline_bottom(part = "header", border = fp_border(width = 0.5)) %>%
    hline_bottom(part = "body", border = fp_border(width = 1)) %>%
    set_table_properties(layout = "autofit")
}

# Table 1
sig16  <- read_required(file.path(results_dir,
  "02_1_453_discovery/tables/sig16_genes.csv"))
poly3f <- read_required(file.path(results_dir,
  "02_1_453_discovery/tables/genes16_poly3F_adjusted.csv"))

check_cols(poly3f,
  c("gene", "model", "F_stat", "p", "FDR",
    "p1_logFC", "p1_FDR", "p2_logFC", "p2_FDR",
    "p3_logFC", "p3_FDR"),
  "genes16_poly3F_adjusted")

model_use <- "M1_age_PMI_RIN"
poly3f_sel <- poly3f %>%
  filter(model == model_use) %>%
  select(gene, F_stat, F_p = p, F_FDR = FDR,
         p1_logFC, p1_FDR, p2_logFC, p2_FDR,
         p3_logFC, p3_FDR)

table1 <- sig16 %>%
  left_join(poly3f_sel, by = "gene") %>%
  arrange(F_p) %>%
  select(gene, F_stat, F_p, F_FDR,
         p1_logFC, p1_FDR, p2_logFC, p2_FDR,
         p3_logFC, p3_FDR)

save_table(table1, "Table1",
  paste0(
    "Table 1. Sixteen-gene candidate signature and poly3_F statistics. ",
    "Values are from the 453-sample discovery set (02_1_453_discovery). ",
    "F_stat, F_p, and F_FDR refer to the joint F-test across p1, p2, and p3 ",
    "using model ", model_use, " (age spline, sex, RIN, PMI as covariates). ",
    "F_FDR is the BH-corrected p-value across the 16 candidate genes. ",
    "p1_logFC, p2_logFC, and p3_logFC are the coefficients for the first-, ",
    "second-, and third-order polynomial terms. ",
    "p1_FDR, p2_FDR, and p3_FDR are the BH-corrected p-values across all ",
    "14,233 expressed genes. ",
    "The p3 term alone (p_poly3F) has nominal p-values between 1.69e-04 and ",
    "0.136, and local FDR q-values between 0.086 and 0.460; none of the 16 ",
    "genes reaches q < 0.05 (Supplementary Table S2). ",
    "C6orf120 lacks joint F-test values in genes16_poly3F_adjusted.csv ",
    "and is therefore NA for F_stat/F_p/F_FDR and for the p1/p2/p3 ",
    "coefficient columns."
  )
)

# Table 2
bulk_cplx <- read_required(file.path(results_dir,
  "07_5c_OXPHOS_analysis/tables/OXPHOS13_complex_score.csv"))
rgc_cplx  <- read_required(file.path(results_dir,
  "07_3_NDUFB8_coexpression/tables/RGC_complex_100iter_summary.csv"))

check_cols(bulk_cplx, c("complex", "n_subunits", "beta", "p", "FDR"),
           "OXPHOS13_complex_score")
check_cols(rgc_cplx, c("complex", "n_genes", "log2FC_median",
                       "log2FC_IQR_lo", "log2FC_IQR_hi",
                       "pct_p_lt_005", "pct_log2FC_pos"),
           "RGC_complex_100iter_summary")

bulk_cplx <- bulk_cplx %>% mutate(complex_std = normalize_complex(complex))
rgc_cplx  <- rgc_cplx  %>% mutate(complex_std = normalize_complex(complex))

table2 <- full_join(
  bulk_cplx %>% select(complex_std,
                       bulk_n_subunits = n_subunits,
                       bulk_beta = beta, bulk_p = p, bulk_FDR = FDR),
  rgc_cplx %>% select(complex_std,
                      rgc_n_genes = n_genes,
                      rgc_log2FC = log2FC_median,
                      rgc_IQR_lo = log2FC_IQR_lo,
                      rgc_IQR_hi = log2FC_IQR_hi,
                      pct_p_lt_005, pct_log2FC_pos),
  by = "complex_std"
) %>%
  mutate(
    complex_std = factor(complex_std,
                         levels = c("Complex I", "Complex II",
                                    "Complex III", "Complex IV", "Complex V")),
    bulk_direction = ifelse(bulk_beta < 0, "down", "up"),
    rgc_direction  = ifelse(rgc_log2FC < 0, "down", "up"),
    concordant     = bulk_direction == rgc_direction
  ) %>%
  arrange(complex_std)

perSD <- read_required(file.path(results_dir,
  "07_5c_OXPHOS_analysis/tables/OXPHOS13_composite_perSD.csv")) %>%
  filter(model == "adjusted") %>%
  mutate(
    complex_std      = "Composite (per 1 SD)",
    bulk_n_subunits  = 13L,
    bulk_beta        = beta_per_SD,
    bulk_p           = p_value,
    bulk_FDR         = NA_real_,
    rgc_n_genes      = NA_integer_,
    rgc_log2FC       = NA_real_,
    rgc_IQR_lo       = NA_real_,
    rgc_IQR_hi       = NA_real_,
    pct_p_lt_005     = NA_real_,
    pct_log2FC_pos   = NA_real_,
    bulk_direction   = ifelse(bulk_beta < 0, "down", "up"),
    rgc_direction    = NA_character_,
    concordant       = NA
  ) %>%
  select(complex_std, bulk_n_subunits, bulk_beta, bulk_p, bulk_FDR,
         rgc_n_genes, rgc_log2FC, rgc_IQR_lo, rgc_IQR_hi,
         pct_p_lt_005, pct_log2FC_pos,
         bulk_direction, rgc_direction, concordant)

pathway_perSD <- read_required(file.path(results_dir,
  "07_5c_OXPHOS_analysis/tables/OXPHOS_pathway_perSD.csv")) %>%
  filter(model == "adjusted") %>%
  mutate(
    complex_std      = "OXPHOS pathway (per 1 SD)",
    bulk_n_subunits  = 187L,
    bulk_beta        = beta_per_SD,
    bulk_p           = p_value,
    bulk_FDR         = NA_real_,
    rgc_n_genes      = NA_integer_,
    rgc_log2FC       = NA_real_,
    rgc_IQR_lo       = NA_real_,
    rgc_IQR_hi       = NA_real_,
    pct_p_lt_005     = NA_real_,
    pct_log2FC_pos   = NA_real_,
    bulk_direction   = ifelse(bulk_beta < 0, "down", "up"),
    rgc_direction    = NA_character_,
    concordant       = NA
  ) %>%
  select(complex_std, bulk_n_subunits, bulk_beta, bulk_p, bulk_FDR,
         rgc_n_genes, rgc_log2FC, rgc_IQR_lo, rgc_IQR_hi,
         pct_p_lt_005, pct_log2FC_pos,
         bulk_direction, rgc_direction, concordant)

table2$complex_std <- as.character(table2$complex_std)
table2 <- bind_rows(table2, perSD, pathway_perSD) %>%
  mutate(complex_std = factor(complex_std,
           levels = c("Complex I", "Complex II", "Complex III",
                      "Complex IV", "Complex V",
                      "Composite (per 1 SD)",
                      "OXPHOS pathway (per 1 SD)")))

log_info(sprintf("Table 2: %d rows", nrow(table2)))

save_table(table2, "Table2",
  paste0(
    "Table 2. OXPHOS complex-level bulk and single-cell comparison, ",
    "with composite and pathway scores per 1 SD. ",
    "Bulk complex scores were computed as the mean z-score across the genes ",
    "of each complex in the 13-gene eQTM-enriched OXPHOS subset ",
    "(Complex I: 6 genes; Complex II: 1 gene; Complex III: 3 genes; ",
    "Complex IV: 2 genes; Complex V: not represented in the 13-gene subset). ",
    "Bulk beta and FDR are from a linear model with age, sex, RIN, and PMI ",
    "as covariates. ",
    "The 'Composite (per 1 SD)' row reports the adjusted per-1-SD beta ",
    "for the composite score across all 13 OXPHOS genes ",
    "(beta = -0.025, 95% CI -0.129 to 0.078, p = 0.631; n = 418 after ",
    "excluding one technical outlier). ",
    "The 'OXPHOS pathway (per 1 SD)' row reports the adjusted per-1-SD beta ",
    "for the full OXPHOS pathway score (187 genes), from a standardized ",
    "model with the same covariates ",
    "(beta = -0.080, 95% CI -0.184 to 0.025, p = 0.135; n = 418). ",
    "Single-cell log2FC, IQR, and iteration percentages are ",
    "from 100 subsampling iterations (114 RGC per donor) comparing dAMD to ",
    "healthy donors. Within the 13-gene eQTM-enriched subset, Complex I had ",
    "the most negative beta (-0.0171), followed by Complex III (-0.0151) and ",
    "Complex IV (-0.0106); no complex reached FDR < 0.05. The 'Complex III ",
    "dominance' description applies only to this eQTM-enriched subset and ",
    "should not be generalised to the full OXPHOS set. Complex II has only ",
    "one gene (SDHB) and a bulk complex-level beta could not be computed; ",
    "Complex II is therefore NA for bulk_n_subunits, bulk_beta, bulk_p, and ",
    "bulk_FDR. ",
    "NDUFB3, a Complex I subunit within the 13-gene eQTM-enriched subset, ",
    "was further evaluated by two-sample MR and colocalization analysis ",
    "(Table 3). A systematic screening of four independent retinal eQTL ",
    "resources for NDUFB3 is presented in Table 4."
  )
)

# Table 3
coloc <- read_required(file.path(results_dir,
  "07_7f_coloc_NDUFB3/coloc_NDUFB3_both_tissues.csv"))
check_cols(coloc,
  c("analysis", "tissue", "n_samples", "PP_H1", "PP_H3", "PP_H4",
    "r2_top_variants", "direction_concordant"),
  "coloc_NDUFB3_both_tissues")

table3 <- coloc %>%
  mutate(
    analysis = factor(analysis,
      levels = c("NDUFB3_adipose_vs_earlyAMD",
                 "NDUFB3_blood_vs_earlyAMD")),
    tissue = factor(tissue,
      levels = c("Adipose_Subcutaneous", "Whole_Blood")),
    n_samples = as.integer(n_samples)
  ) %>%
  arrange(analysis)

adipose_row <- table3 %>% filter(tissue == "Adipose_Subcutaneous")
blood_row   <- table3 %>% filter(tissue == "Whole_Blood")

pp_h4_adi <- round(adipose_row$PP_H4, 4)
pp_h3_adi <- round(adipose_row$PP_H3, 4)
pp_h1_adi <- round(adipose_row$PP_H1, 4)
r2_adi    <- round(adipose_row$r2_top_variants, 2)

pp_h4_blo <- round(blood_row$PP_H4, 4)
pp_h3_blo <- round(blood_row$PP_H3, 4)
pp_h1_blo <- round(blood_row$PP_H1, 4)
r2_blo    <- round(blood_row$r2_top_variants, 2)

log_info(sprintf("Table 3: adipose PP.H4=%.4f PP.H3=%.4f PP.H1=%.4f r2=%.2f",
                 pp_h4_adi, pp_h3_adi, pp_h1_adi, r2_adi))
log_info(sprintf("Table 3: blood PP.H4=%.4f PP.H3=%.4f PP.H1=%.4f r2=%.2f",
                 pp_h4_blo, pp_h3_blo, pp_h1_blo, r2_blo))

save_table(table3, "Table3",
  paste0(
    "Table 3. Colocalization analysis between NDUFB3 cis-eQTL and early AMD GWAS. ",
    "Colocalization was performed with COLSTATS v2.3 ",
    "(https://colstats.irgb.cnr.it) in targeted mode with default priors ",
    "(p1 = 1e-4, p2 = 1e-4, p12 = 1e-5) and 1000 Genomes phase 3 EUR as ",
    "the LD reference. The analysis region was chr2:201,120,000-203,133,000 ",
    "(GRCh38). Summary statistics for NDUFB3 cis-eQTLs were obtained from ",
    "GTEx v6p (Adipose_Subcutaneous, n = 581; Whole_Blood, n = 670). ",
    "The COLSTATS raw output for adipose tissue and web-session screenshots ",
    "for both tissues are retained in results/07_7f_coloc_NDUFB3/ for ",
    "verification. ",
    "Early AMD GWAS summary statistics were obtained from ebi-a-GCST010723 ",
    "(IAMDGC 2018, 14,034 cases and 91,214 controls). ",
    "Only variants present in both datasets were retained. ",
    "PP.H4 is the posterior probability of a shared causal variant. ",
    "Subcutaneous adipose tissue showed moderate colocalization ",
    sprintf("(PP.H4 = %.4f, PP.H3 = %.4f, PP.H1 = %.4f), ",
            pp_h4_adi, pp_h3_adi, pp_h1_adi),
    "with the top colocalizing variants in high linkage disequilibrium ",
    sprintf("(r2 = %.2f) and concordant direction of effect. ", r2_adi),
    "Whole blood showed no colocalization ",
    sprintf("(PP.H4 = %.4f, PP.H3 = %.4f, PP.H1 = %.4f), ",
            pp_h4_blo, pp_h3_blo, pp_h1_blo),
    "and the top eQTL variants in the two tissues were completely unlinked ",
    sprintf("(r2 = %.2f). ", r2_blo),
    "This tissue-specific and partially discordant regulatory architecture ",
    "explains the instrument-fragile nature of the blood-derived NDUFB3 MR ",
    "estimate and supports the interpretation that blood eQTL instruments ",
    "are not valid proxies for NDUFB3 regulation in disease-relevant tissues."
  )
)

# Table 4
table4 <- read_required(file.path(results_dir,
  "07_7g_retinal_eqtl_screening_NDUFB3/tables/NDUFB3_retinal_eQTL_screening_summary.csv"))
check_cols(table4,
  c("Source", "Sample_size", "Tissue", "SNP", "FDR", "Beta",
    "MAF", "eQTL_detected", "Usable_for_MR", "Notes"),
  "NDUFB3_retinal_eQTL_screening_summary")

log_info(sprintf("Table 4: %d rows (incl. Wang 2023)", nrow(table4)))

wang_row <- table4 %>% filter(grepl("Wang", Source, ignore.case = TRUE))
wang_notes <- wang_row$Notes[1]
log_info(sprintf("Table 4: Wang row notes = %s", wang_notes))

save_table(table4, "Table4",
  paste0(
    "Table 4. Systematic screening of four independent retinal eQTL ",
    "resources for NDUFB3 cis-eQTL. ",
    "Four independent retinal eQTL resources were queried for NDUFB3 ",
    "cis-eQTL: the Orozco eye-eQTL database (GSE135092, macula retina and ",
    "non-macula RPE, n = 106), the Strunz retinal eQTL mega-analysis ",
    "(311 healthy retinas), EyeGEx (Ratnapriya et al. 2019; 406 retinas), ",
    "and the single-cell multiomics retinal resource of Wang et al. 2023 ",
    "(20 healthy donors; 192,792 snRNA-seq nuclei and 245,940 snATAC-seq ",
    "nuclei; GSE247157/GSE226108). ",
    "In the Orozco macula retina dataset, a cis-eQTL for NDUFB3 was ",
    "identified (rs13413479, FDR = 0.0248, beta = -0.23, MAF = 0.14), ",
    "but this SNP also strongly regulates the adjacent gene PPIL3 ",
    "(FDR = 8.83e-12), precluding its use as a clean instrument due to ",
    "potential horizontal pleiotropy. A second NDUFB3-specific SNP in ",
    "non-macula RPE (rs71406913, FDR = 0.0409, beta = -1.13) had ",
    "MAF = 0.01, insufficient for MR. ",
    "NDUFB3 was not among the 3,751 independent eQTL signals in the ",
    "Strunz mega-analysis (311 retinas) or the 10,474 eGenes in EyeGEx ",
    "(406 retinas). ",
    "In the Wang et al. 2023 single-cell multiomics resource, NDUFB3 ",
    "did not reach gene-level significance as a sc-eGene in any major ",
    "retinal cell type (gene-level FDR < 0.1), and 10 sc-caQTLs were ",
    "identified within the NDUFB3 cis region (chr2:201,850,000-203,860,000, ",
    "Hg19), but none targeted NDUFB3 itself; 40 significant sc-eQTLs were ",
    "identified in this region, and all targeted neighbouring genes, ",
    "including TRAK2, CASP8, CASP10, NOP58, CFLAR-AS1, PPIL3, STRADB, ",
    "ALS2, and BMPR2. Notably, PPIL3 sc-eQTLs were identified in ",
    "amacrine and horizontal cells, independently replicating the ",
    "horizontal pleiotropy observed in the Orozco dataset. ",
    "These findings indicate that NDUFB3 lacks a robust, retina-specific ",
    "genetic instrument, and its prioritisation remains provisional ",
    "pending larger retinal eQTL datasets with suitable instrument ",
    "characteristics. Full results are in Supplementary Results S17 and ",
    "Supplementary Table S10."
  )
)

# Word output
log_info("--- Writing Word (landscape) ---")

ft_list <- list()
ft_list$Table1 <- make_ft(table1 %>% mutate(across(where(is.numeric), fmt_num)),
                          font_size = 8)
ft_list$Table2 <- make_ft(table2 %>% mutate(
                            across(where(is.numeric), fmt_num),
                            across(where(is.character),
                                   ~ ifelse(is.na(.x), "\u2014", .x))),
                          font_size = 8)
ft_list$Table3 <- make_ft(table3 %>% mutate(
                            across(where(is.numeric) & !any_of("n_samples"), fmt_num),
                            n_samples = as.character(n_samples)),
                          font_size = 8)

table4_fmt <- table4
table4_fmt$Sample_size <- suppressWarnings(as.integer(table4_fmt$Sample_size))
ft_list$Table4 <- make_ft(table4_fmt, font_size = 8)

docx_file <- file.path(table_dir, "10_1B_main_tables.docx")

doc <- read_docx() %>%
  body_set_default_section(
    prop_section(
      page_margins = page_mar(bottom = 0.5, top = 0.5,
                              right = 0.5, left = 0.5),
      type = "continuous",
      page_size = page_size(orient = "landscape")
    )
  ) %>%
  body_add_par("Main Tables", style = "heading 1") %>%
  body_add_par("", style = "Normal")

for (nm in names(ft_list)) {
  doc <- doc %>%
    body_add_par(nm, style = "heading 2") %>%
    body_add_flextable(ft_list[[nm]]) %>%
    body_add_par("", style = "Normal")
}

print(doc, target = docx_file)
log_info(sprintf("saved Word: %s", docx_file))

log_info("=== 10_1B v7.1 done (landscape) ===")

