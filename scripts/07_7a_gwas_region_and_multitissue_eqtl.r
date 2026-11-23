# ============================================================
# 07_7a_gwas_region_and_multitissue_eqtl.r (v3.2)
# Purpose: (A) Extract gene region from AMD GWAS VCF
#          (B) Multi-tissue eQTL from GTEx
#
# v3.2 Fix：
# ============================================================

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(ggplot2)
  library(httr2); library(jsonlite)
})

# ============================================================
# ============================================================
save_table <- function(df, dir, filename) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  filepath <- file.path(dir, filename)
  write.csv(df, filepath, row.names = FALSE)
  log_info(sprintf("Saved: %s (%d rows)", filename, nrow(df)))
  invisible(filepath)
}

save_fig <- function(p, dir, filename, width = 8, height = 6) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  f_pdf <- file.path(dir, paste0(filename, ".pdf"))
  f_png <- file.path(dir, paste0(filename, ".png"))
  ggsave(f_pdf, p, width = width, height = height)
  ggsave(f_png, p, width = width, height = height, dpi = 300)
  log_info(sprintf("Saved: %s.pdf / .png", filename))
  invisible(f_pdf)
}

# ============================================================
# 1. Parameters
# ============================================================
args <- commandArgs(trailingOnly = TRUE)
STEP <- if (length(args) >= 1) args[1] else "all"
GENE <- if (length(args) >= 2) args[2] else "NDUFB3"

DIR_THIS <- setup_output_dir(sprintf("07_7a_%s", GENE))
log_info(sprintf("=== 07_7a %s (step=%s) v3.2 start ===", GENE, STEP))

# ---- Input paths ----
F_GWAS_VCF    <- "E:/retina_tools/pubdata/eQTL/gwas/ebi-a-GCST010723.vcf.gz"
F_GWAS_REGION <- sprintf("E:/retina_tools/pubdata/eQTL/gwas/GCST010723_%s_region.rds", GENE)
F_GTEX_LOCAL  <- sprintf("E:/retina_tools/pubdata/eQTL/GTEx/%s_eQTL_by_tissue.csv", GENE)
GTEX_V8_DIR   <- "E:/retina_tools/pubdata/eQTL/eQTLGen/GTEx_Analysis_v8_eQTL"

# ---- Gene coordinates (GRCh37) ----
GENE_COORD <- list(
  NDUFB3 = list(chr = "2",  start = 201934000L, end = 201947000L),
  UQCRQ  = list(chr = "5",  start = 132202000L, end = 132204000L),
  NDUFB8 = list(chr = "10", start = 102288000L, end = 102290000L),
  SDHB   = list(chr = "1",  start = 173452000L, end = 173477000L)
)

GENE_GENCODE <- list(
  NDUFB3 = "ENSG00000119013.10",
  UQCRQ  = "ENSG00000164442.13",
  NDUFB8 = "ENSG00000166136.15",
  SDHB   = "ENSG00000117118.11"
)

if (!GENE %in% names(GENE_COORD))
  stop(sprintf("Coordinates for %s not defined", GENE))

COORD <- GENE_COORD[[GENE]]
CIS_HALF_WINDOW <- 1e6

# ============================================================
# Step A: Extract GWAS region
# ============================================================
step_extract <- function() {
  log_info("--- Step A: Extract GWAS region ---")

  START <- max(0L, COORD$start - CIS_HALF_WINDOW)
  END   <- COORD$end + CIS_HALF_WINDOW

  log_info(sprintf("Region: chr%s:%d-%d", COORD$chr, START, END))

  if (!file.exists(F_GWAS_VCF)) {
    log_warn(sprintf("GWAS VCF not found: %s", F_GWAS_VCF))
    return(NULL)
  }

  con <- gzfile(F_GWAS_VCF, "r")

  repeat {
    l <- readLines(con, n = 1)
    if (length(l) == 0) break
    if (!startsWith(l, "#")) break
  }

  gwas_list <- list()
  n_read <- 0L; n_kept <- 0L

  repeat {
    lines <- readLines(con, n = 100000L)
    if (length(lines) == 0) break
    n_read <- n_read + length(lines)

    keep <- grepl(sprintf("^%s\t", COORD$chr), lines)
    if (any(keep)) {
      for (l in lines[keep]) {
        f <- strsplit(l, "\t", fixed = TRUE)[[1]]
        if (length(f) < 10) next
        pos <- as.integer(f[2])
        if (is.na(pos) || pos < START || pos > END) next

        fmt_keys <- strsplit(f[9], ":", fixed = TRUE)[[1]]
        fmt_vals <- strsplit(f[10], ":", fixed = TRUE)[[1]]

        d <- data.frame(snp = f[3], chr = f[1], pos = pos,
                        ref = f[4], alt = f[5],
                        stringsAsFactors = FALSE)
        for (k in c("ES", "SE", "LP", "AF", "SS", "NC", "EZ")) {
          idx <- which(fmt_keys == k)
          d[[k]] <- if (length(idx) > 0 && length(fmt_vals) >= idx) {
            suppressWarnings(as.numeric(fmt_vals[idx]))
          } else NA_real_
        }
        gwas_list[[length(gwas_list) + 1L]] <- d
        n_kept <- n_kept + 1L
      }
    }
    if (n_read %% 1e6 == 0)
      log_info(sprintf("  %d rows read, %d kept", n_read, n_kept))
  }
  close(con)

  if (length(gwas_list) == 0) {
    log_warn("No variants found in the region")
    return(NULL)
  }

  gwas_df <- do.call(rbind, gwas_list)
  gwas_df$pval <- 10^(-gwas_df$LP)
  gwas_df$z    <- gwas_df$ES / gwas_df$SE

  saveRDS(gwas_df, F_GWAS_REGION)

  summary_df <- data.frame(
    gene = GENE, chr = COORD$chr,
    region_start = START, region_end = END,
    n_variants = nrow(gwas_df),
    n_with_ES_SE = sum(!is.na(gwas_df$ES) & !is.na(gwas_df$SE)),
    n_p_lt_0.05 = sum(gwas_df$pval < 0.05, na.rm = TRUE),
    n_p_lt_5e8 = sum(gwas_df$pval < 5e-8, na.rm = TRUE),
    min_p = min(gwas_df$pval, na.rm = TRUE),
    stringsAsFactors = FALSE)
  save_table(summary_df, DIR_THIS, sprintf("%s_GWAS_region_summary.csv", GENE))

  cat("\n=== GWAS region summary ===\n")
  print(summary_df, row.names = FALSE, digits = 3)

  gwas_df
}

# ============================================================
# Step B: Multi-tissue eQTL (GTEx)
# ============================================================
step_eqtl <- function() {
  log_info("--- Step B: Multi-tissue eQTL (GTEx) ---")

  gtex <- NULL

  if (!GENE %in% names(GENE_GENCODE)) {
    log_warn(sprintf("%s 的 GENCODE ID 未定义，尝试本地文件", GENE))
  } else {
    gencode_id <- GENE_GENCODE[[GENE]]
    log_info(sprintf("GTEx API 拉取 %s (%s)", GENE, gencode_id))

    for (ds in c("gtex_v10", "gtex_v8")) {
      log_info(sprintf("  尝试 datasetId = %s", ds))

      gtex <- tryCatch({
        resp <- request("https://gtexportal.org/api/v2/association/singleTissueEqtl") |>
          req_url_query(gencodeId = gencode_id, datasetId = ds,
                        itemsPerPage = 1000) |>
          req_perform()
        result <- resp_body_json(resp, simplifyVector = TRUE)

        if (is.null(result$data) || length(result$data) == 0) {
          log_warn(sprintf("  %s 返回空数据", ds))
          NULL
        } else {
          df <- as.data.frame(result$data)
          log_info(sprintf("  API 返回列: %s", paste(colnames(df), collapse = ", ")))
          log_info(sprintf("  API 返回行数: %d", nrow(df)))

          get_col <- function(df, candidates, default = NA) {
            for (cn in candidates) {
              if (cn %in% colnames(df)) return(df[[cn]])
            }
            rep(default, nrow(df))
          }

          out <- data.frame(
            tissue     = get_col(df, c("tissueSiteDetailId", "tissueSiteDetail",
                                        "tissue", "tissueId")),
            nes        = as.numeric(get_col(df, c("nes", "Nes", "slope",
                                                   "slopeOfBeta"))),
            pval       = as.numeric(get_col(df, c("pValue", "pval", "p_value",
                                                   "pvalue"))),
            m_value    = as.numeric(get_col(df, c("mValue", "mvalue", "m_value"))),
            n_variants = as.integer(get_col(df, c("nVariants", "nvariants",
                                                   "n_variants"), default = NA_integer_)),
            stringsAsFactors = FALSE)

          out <- out[!is.na(out$nes), ]
          if (nrow(out) == 0) {
            log_warn(sprintf("  %s 提取后无有效数据", ds))
            NULL
          } else {
            out <- out %>%
              group_by(tissue) %>%
              slice_max(abs(nes), n = 1, with_ties = FALSE) %>%
              ungroup()

            out <- out[order(-abs(out$nes)), ]
            log_info(sprintf("  %s 成功: %d 个组织", ds, nrow(out)))
            out
          }
        }
      }, error = function(e) {
        log_warn(sprintf("  %s API 失败: %s", ds, conditionMessage(e)))
        NULL
      })

      if (!is.null(gtex) && nrow(gtex) > 0) break
    }
  }

  if (is.null(gtex) && file.exists(F_GTEX_LOCAL)) {
    log_info(sprintf("回退到本地 CSV: %s", F_GTEX_LOCAL))
    gtex <- fread(F_GTEX_LOCAL)
    setnames(gtex, tolower(colnames(gtex)), skip_absent = TRUE)
    if ("nes" %in% colnames(gtex) && "tissue" %in% colnames(gtex)) {
      gtex <- gtex %>%
        group_by(tissue) %>%
        slice_max(abs(nes), n = 1, with_ties = FALSE) %>%
        ungroup() %>%
        arrange(-abs(nes))
    }
  }

  if (is.null(gtex)) {
    if (dir.exists(GTEX_V8_DIR)) {
      log_info("从 GTEx v8 egenes 文件提取...")

      files <- list.files(GTEX_V8_DIR,
                          pattern = "\\.v8\\.egenes\\.txt\\.gz$",
                          full.names = TRUE)
      log_info(sprintf("找到 %d 个组织 egenes 文件", length(files)))

      ens_full  <- GENE_GENCODE[[GENE]]
      ens_nover <- sub("\\..*$", "", ens_full)
      gene_sym  <- GENE

      gtex_list <- list()
      for (f in files) {
        tissue <- sub("\\.v8\\.egenes\\.txt\\.gz$", "", basename(f))

        res <- tryCatch({
          dt <- fread(f)

          dt_sub <- dt[gene_id == ens_full]
          if (nrow(dt_sub) == 0) dt_sub <- dt[gene_id == ens_nover]
          if (nrow(dt_sub) == 0 && "gene_name" %in% colnames(dt))
            dt_sub <- dt[gene_name == gene_sym]

          if (nrow(dt_sub) > 0) {
            data.frame(
              tissue     = tissue,
              nes        = as.numeric(dt_sub$slope[1]),
              pval       = as.numeric(dt_sub$pval_nominal[1]),
              qval       = if ("qval" %in% colnames(dt_sub))
                             as.numeric(dt_sub$qval[1]) else NA_real_,
              m_value    = NA_real_,
              n_variants = if ("num_var" %in% colnames(dt_sub))
                             as.integer(dt_sub$num_var[1]) else NA_integer_,
              stringsAsFactors = FALSE)
          } else NULL
        }, error = function(e) {
          log_warn(sprintf("  %s 读取失败: %s", tissue, conditionMessage(e)))
          NULL
        })

        if (!is.null(res)) gtex_list[[tissue]] <- res
      }

      if (length(gtex_list) > 0) {
        gtex <- do.call(rbind, gtex_list)
        gtex <- gtex[order(-abs(gtex$nes)), ]
        rownames(gtex) <- NULL
        log_info(sprintf("从 egenes 提取到 %d 个组织是 eGene", nrow(gtex)))
      } else {
        log_warn(sprintf("%d 个组织文件中均未找到 %s（不是任何组织的 eGene）",
                         length(files), GENE))
      }
    } else {
      log_warn(sprintf("GTEx v8 目录不存在: %s", GTEX_V8_DIR))
    }
  }

  # ---------- Output ----------
  if (is.null(gtex) || nrow(gtex) == 0) {
    log_warn("GTEx 数据不可用")
    return(NULL)
  }

  save_table(gtex, DIR_THIS, sprintf("%s_multi_tissue_eQTL.csv", GENE))
  log_info(sprintf("Saved: %s_multi_tissue_eQTL.csv (%d tissues)",
                   GENE, nrow(gtex)))

  cat("\n=== Multi-tissue eQTL ===\n")
  print(gtex, row.names = FALSE, digits = 3)

  gtex <- gtex[order(gtex$nes), ]
  gtex$tissue <- factor(gtex$tissue, levels = unique(gtex$tissue))

  p <- ggplot(gtex, aes(x = tissue, y = nes,
                        fill = nes > 0)) +
    geom_col() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    coord_flip() + theme_bw(base_size = 11) +
    scale_fill_manual(values = c("TRUE" = "#C44E52", "FALSE" = "#4C72B0"),
                      name = "Positive slope", guide = "none") +
    labs(x = "", y = "eQTL effect size (slope)",
         title = sprintf("%s: multi-tissue eQTL (GTEx)", GENE))

  save_fig(p, DIR_THIS, sprintf("%s_multi_tissue_eQTL", GENE),
           width = 8, height = max(4, nrow(gtex) * 0.3))

  gtex
}

# ============================================================
# Main
# ============================================================
if (STEP %in% c("extract", "all")) step_extract()
if (STEP %in% c("eqtl", "all"))    step_eqtl()

log_info(sprintf("=== 07_7a %s (step=%s) v3.2 done ===", GENE, STEP))

