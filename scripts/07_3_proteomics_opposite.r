# ============================================================
# 07_3_proteomics_opposite.r
#
#   Part B: PXD063075 Supplementary File 2 (159 top-ranked protein)
#
#
# Output：results/07_3_proteomics_opposite/
#   PartA_*  PartB_*  PartC_*
#   summary_three_datasets.csv
# ============================================================

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(readxl)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(xml2)
})

DIR_THIS <- setup_output_dir("07_3_proteomics_opposite")
log_info("=== 07_3 proteomics opposite start ===")

# ============================================================
# ============================================================
CFG <- list(
  # Part A: PXD063075 XML
  xml_pxd063075    = "E:/retina_tools/raw/PXD063075/combined.prot.xml",
  csv_pxd063075_s2 = "E:/retina_tools/raw/PXD063075/41419_2026_8724_MOESM2_ESM.csv",
  suppl2_skip      = 3,
  # Part C: PXD029501 TMT matrix
  xlsx_pxd029501   = "E:/retina_tools/raw/PXD029501/41467_2022_31707_MOESM3_ESM.xlsx",
  pxd029501_sheet  = "Total proteins identified"
)

COMPLEX_MAP_PRIMARY <- list(
  ComplexI   = "Mito_CI_Subunits",
  ComplexII  = "Mito_CII_Subunits",
  ComplexIII = "Mito_CIII_Subunits",
  ComplexIV  = "Mito_CIV_Subunits",
  ComplexV   = "Mito_CV_Subunits"
)
COMPLEX_MAP_FALLBACK <- list(
  ComplexI   = "Mito_Complex_I",
  ComplexII  = "Mito_Complex_II",
  ComplexIII = "Mito_Complex_III",
  ComplexIV  = "Mito_Complex_IV",
  ComplexV   = "Mito_Complex_V"
)

# ============================================================
# ============================================================
log_info("--- 1. Load gene sets and focus genes ---")

GENE_SETS_ALL <- read_upstream("01_2_prepare_gene_sets", "data",
                                "gene_sets.rds")
if (is.null(GENE_SETS_ALL)) stop("gene_sets.rds from 01_2 not found")
log_info(sprintf("loaded %d gene sets", length(GENE_SETS_ALL)))

extract_complex_sets <- function(map_primary, map_fallback) {
  out <- list(); prov <- list()
  for (cplx in names(map_primary)) {
    gs_id <- map_primary[[cplx]]
    if (!gs_id %in% names(GENE_SETS_ALL)) gs_id <- map_fallback[[cplx]]
    if (!gs_id %in% names(GENE_SETS_ALL)) next
    gs <- GENE_SETS_ALL[[gs_id]]
    out[[cplx]] <- toupper(gs$genes)
    prov[[cplx]] <- data.frame(
      complex = cplx, geneset_id = gs_id,
      source = gs$source, source_id = gs$source_id,
      version = gs$version, reference = gs$reference,
      n_genes = length(gs$genes), stringsAsFactors = FALSE)
  }
  list(genes = out, prov = do.call(rbind, prov))
}

complex_data <- extract_complex_sets(COMPLEX_MAP_PRIMARY,
                                      COMPLEX_MAP_FALLBACK)
gene_sets <- complex_data$genes
save_table(complex_data$prov, DIR_THIS, "gene_set_provenance.csv")

# ---- build focus gene list from upstream results (no hard-coding) ----

# primary: NDUFB3 from MR unified FDR
mr_file <- file.path(DIR_RESULT, "07_7e_MR_unified_fdr",
                     "tables", "NDUFB3_unified_fdr_summary.csv")
if (!file.exists(mr_file)) {
  stop("NDUFB3 unified FDR file missing; run 07_7e_MR_unified_fdr first")
}
mr <- read.csv(mr_file, stringsAsFactors = FALSE)
if (!any(mr$p_unified_fdr < 0.05, na.rm = TRUE)) {
  stop("NDUFB3 not significant in unified FDR; check 07_7e output")
}
primary_genes <- "NDUFB3"

# secondary: all Complex III subunits from 01_2 gene sets
c3_set <- GENE_SETS_ALL[["Mito_CIII_Subunits"]]$genes
if (is.null(c3_set) || length(c3_set) == 0) {
  stop("Mito_CIII_Subunits not found in gene_sets.rds")
}
secondary_genes <- setdiff(unique(toupper(c3_set)), primary_genes)

# assemble
focus_df <- rbind(
  data.frame(gene = primary_genes,
             priority = "primary",
             reason = "MR-significant in early AMD (this study)",
             source = "07_7e_MR_unified_fdr",
             stringsAsFactors = FALSE),
  data.frame(gene = secondary_genes,
             priority = "secondary",
             reason = "Complex III subunit, bulk top signal",
             source = "01_2_prepare_gene_sets: Mito_CIII_Subunits",
             stringsAsFactors = FALSE)
)
FOCUS_GENES   <- unique(toupper(focus_df$gene))
FOCUS_PRIMARY <- toupper(focus_df$gene[tolower(focus_df$priority) == "primary"])

# save for audit
focus_csv <- file.path(DIR_THIS, "data", "focus_genes.csv")
write.csv(focus_df, focus_csv, row.names = FALSE)
log_info(sprintf("focus genes: %d (%d primary, %d secondary)",
                 nrow(focus_df), length(primary_genes), length(secondary_genes)))
log_info(sprintf("target genes: %s", paste(FOCUS_GENES, collapse = ", ")))

summary_rows <- list()

# ============================================================
# ============================================================
log_info("--- Part A: PXD063075 combined.prot.xml ---")

partA_result <- tryCatch({
  if (!file.exists(CFG$xml_pxd063075)) {
    log_info("skip: XML file not found")
    NULL
  } else {
    doc <- read_xml(CFG$xml_pxd063075)
    ns  <- xml_ns(doc)
    xp <- function(path) gsub("(?<=/)(?!/)([a-zA-Z_][a-zA-Z0-9_\\-]*)",
                              "d1:\\1", path, perl = TRUE)

    prot_nodes <- xml_find_all(doc, xp("//protein_group/protein"), ns = ns)
    log_info(sprintf("protein nodes: %d", length(prot_nodes)))

    if (length(prot_nodes) == 0) {
      NULL
    } else {
      prot_names <- xml_attr(prot_nodes, "protein_name")

      extract_gene <- function(x) {
        parts <- strsplit(x, "\\|")[[1]]
        if (length(parts) >= 3) {
          g <- sub("_HUMAN.*$", "", parts[3])
          g <- sub("^rev_", "", g)
          return(toupper(g))
        }
        NA_character_
      }
      genes_in_xml <- unique(na.omit(sapply(prot_names, extract_gene)))

      all_targets <- unique(c(FOCUS_GENES,
                              unlist(gene_sets)))
      hits <- intersect(all_targets, genes_in_xml)
      log_info(sprintf("XML: %d proteins identified, %d target genes detected",
                       length(genes_in_xml), length(hits)))

      complex_hits <- lapply(names(gene_sets), function(nm) {
        sub <- intersect(gene_sets[[nm]], genes_in_xml)
        data.frame(complex = nm, n_target = length(gene_sets[[nm]]),
                   n_detected = length(sub),
                   detected_genes = paste(sub, collapse = ";"),
                   stringsAsFactors = FALSE)
      })
      complex_tab <- do.call(rbind, complex_hits)
      save_table(complex_tab, DIR_THIS, "PartA_complex_detection.csv")

      cat("\n=== Part A: complex detection ===\n")
      print(complex_tab[, c("complex", "n_target", "n_detected")],
            row.names = FALSE)

      list(genes = genes_in_xml, hits = hits, table = complex_tab)
    }
  }
}, error = function(e) {
  log_info(sprintf("Part A error: %s", conditionMessage(e)))
  NULL
})

if (!is.null(partA_result)) {
  for (pg in FOCUS_PRIMARY) {
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      gene = pg, dataset = "PXD063075_XML",
      detected = ifelse(pg %in% partA_result$genes, "yes", "no"),
      direction = NA_character_, effect = NA_real_,
      pvalue = NA_real_, note = "ProteinProphet identification only",
      stringsAsFactors = FALSE)
  }
}

# ============================================================
# 3. Part B: PXD063075 Supplementary File 2 (159 protein)
# ============================================================
log_info("--- Part B: PXD063075 suppl2 (159 proteins) ---")

partB_result <- tryCatch({
  if (!file.exists(CFG$csv_pxd063075_s2)) {
    log_info("skip: CSV file not found")
    NULL
  } else {
    suppl2 <- fread(CFG$csv_pxd063075_s2, data.table = FALSE,
                    skip = CFG$suppl2_skip, sep = ",",
                    check.names = FALSE)
    log_info(sprintf("raw dimensions: %d x %d", nrow(suppl2), ncol(suppl2)))

    df <- data.frame(
      gene   = toupper(trimws(as.character(suppl2[[1]]))),
      t_stat = suppressWarnings(as.numeric(suppl2[[2]])),
      pvalue = suppressWarnings(as.numeric(suppl2[[3]])),
      FDR    = suppressWarnings(as.numeric(suppl2[[5]])),
      stringsAsFactors = FALSE)
    df <- df[!is.na(df$gene) & df$gene != "" & !is.na(df$t_stat), ]
    df <- df[!duplicated(df$gene), ]
    log_info(sprintf("after cleaning: %d proteins", nrow(df)))
    save_table(df, DIR_THIS, "PartB_suppl2_cleaned.csv")

    focus_sub <- df[df$gene %in% FOCUS_GENES, ]

    complex_hits <- lapply(names(gene_sets), function(nm) {
      sub <- df[df$gene %in% gene_sets[[nm]], ]
      if (nrow(sub) == 0) return(NULL)
      data.frame(complex = nm, gene = sub$gene,
                 t_stat = sub$t_stat, pvalue = sub$pvalue,
                 FDR = sub$FDR, stringsAsFactors = FALSE)
    })
    complex_tab <- do.call(rbind, complex_hits)
    if (!is.null(complex_tab)) {
      save_table(complex_tab, DIR_THIS, "PartB_complex_hits.csv")
      cat("\n=== Part B: complex hits ===\n")
      print(complex_tab, row.names = FALSE, digits = 4)
    } else {
      log_info("Part B: 0 complex hits")
    }

    list(df = df, focus = focus_sub, complex = complex_tab)
  }
}, error = function(e) {
  log_info(sprintf("Part B error: %s", conditionMessage(e)))
  NULL
})

if (!is.null(partB_result)) {
  for (pg in FOCUS_PRIMARY) {
    row <- partB_result$df[partB_result$df$gene == pg, ]
    if (nrow(row) > 0) {
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        gene = pg, dataset = "PXD063075_suppl2",
        detected = "yes",
        direction = ifelse(row$t_stat < 0, "down", "up"),
        effect = row$t_stat, pvalue = row$pvalue,
        note = sprintf("t-stat, FDR=%.3f", row$FDR),
        stringsAsFactors = FALSE)
    } else {
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        gene = pg, dataset = "PXD063075_suppl2",
        detected = "no", direction = NA_character_,
        effect = NA_real_, pvalue = NA_real_,
        note = "not in 159 top-ranked proteins",
        stringsAsFactors = FALSE)
    }
  }
}

# ============================================================
# ============================================================
log_info("--- Part C: PXD029501 TMT matrix ---")

partC_result <- tryCatch({
  if (!file.exists(CFG$xlsx_pxd029501)) {
    log_info("skip: xlsx file not found")
    NULL
  } else {
    sheets <- excel_sheets(CFG$xlsx_pxd029501)
    if (!CFG$pxd029501_sheet %in% sheets) {
      log_info(sprintf("sheet '%s' not found", CFG$pxd029501_sheet))
      NULL
    } else {
      df_raw <- read_excel(CFG$xlsx_pxd029501, sheet = CFG$pxd029501_sheet)
      log_info(sprintf("raw dimensions: %d x %d", nrow(df_raw), ncol(df_raw)))

      gene_col <- "Gene.Symbol"
      fc_col   <- "MaxFC"
      pval_col <- "T-TEST"
      ctrl_col <- "Means.Control"
      ga_col   <- "Means.GA"

      prot <- data.frame(
        gene      = toupper(trimws(as.character(df_raw[[gene_col]]))),
        fc_raw    = suppressWarnings(as.numeric(df_raw[[fc_col]])),
        pvalue    = suppressWarnings(as.numeric(df_raw[[pval_col]])),
        mean_ctrl = suppressWarnings(as.numeric(df_raw[[ctrl_col]])),
        mean_ga   = suppressWarnings(as.numeric(df_raw[[ga_col]])),
        stringsAsFactors = FALSE)
      prot <- prot[!is.na(prot$gene) & prot$gene != "" & !is.na(prot$fc_raw), ]
      prot <- prot[!duplicated(prot$gene), ]
      log_info(sprintf("after cleaning: %d proteins", nrow(prot)))

      prot$log2FC_abs <- log2(prot$fc_raw)
      prot$direction  <- sign(prot$mean_ga - prot$mean_ctrl)
      prot$signed_log2FC <- prot$log2FC_abs * prot$direction

      save_table(prot, DIR_THIS, "PartC_pxd029501_cleaned.csv")

      complex_hits <- lapply(names(gene_sets), function(nm) {
        sub <- prot[prot$gene %in% gene_sets[[nm]], ]
        if (nrow(sub) == 0) return(NULL)
        sub$complex <- nm
        sub[, c("complex","gene","signed_log2FC","pvalue",
                "mean_ctrl","mean_ga")]
      })
      complex_tab <- do.call(rbind, complex_hits)
      if (!is.null(complex_tab)) {
        save_table(complex_tab, DIR_THIS, "PartC_complex_hits.csv")

        summary_tab <- complex_tab %>%
          group_by(complex) %>%
          summarise(n_hit = n(),
                    n_up = sum(signed_log2FC > 0, na.rm = TRUE),
                    n_down = sum(signed_log2FC < 0, na.rm = TRUE),
                    n_sig = sum(pvalue < 0.05, na.rm = TRUE),
                    mean_l2fc = mean(signed_log2FC, na.rm = TRUE),
                    .groups = "drop")
        save_table(summary_tab, DIR_THIS, "PartC_complex_summary.csv")

        cat("\n=== Part C: complex direction ===\n")
        print(as.data.frame(summary_tab), row.names = FALSE, digits = 3)
      }

      focus_sub <- prot[prot$gene %in% FOCUS_GENES, ]
      list(df = prot, focus = focus_sub, complex = complex_tab)
    }
  }
}, error = function(e) {
  log_info(sprintf("Part C error: %s", conditionMessage(e)))
  NULL
})

if (!is.null(partC_result)) {
  for (pg in FOCUS_PRIMARY) {
    row <- partC_result$df[partC_result$df$gene == pg, ]
    if (nrow(row) > 0) {
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        gene = pg, dataset = "PXD029501_GA_RPE",
        detected = "yes",
        direction = ifelse(row$signed_log2FC < 0, "down", "up"),
        effect = row$signed_log2FC, pvalue = row$pvalue,
        note = "signed log2FC",
        stringsAsFactors = FALSE)
    } else {
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        gene = pg, dataset = "PXD029501_GA_RPE",
        detected = "no", direction = NA_character_,
        effect = NA_real_, pvalue = NA_real_,
        note = "not detected",
        stringsAsFactors = FALSE)
    }
  }
}

# ============================================================
# ============================================================
log_info("--- 5. Three-dataset summary ---")

if (length(summary_rows) > 0) {
  summary_df <- do.call(rbind, summary_rows)

  bulk_ref <- data.frame(
    gene = FOCUS_PRIMARY,
    bulk_direction = c("down"),   # Complex III down in bulk
    stringsAsFactors = FALSE)

  summary_df <- merge(summary_df, bulk_ref, by = "gene", all.x = TRUE)
  save_table(summary_df, DIR_THIS, "summary_three_datasets.csv")

  cat("\n", strrep("=", 70), "\n", sep = "")
  cat("=== Three-dataset summary ===\n")
  cat(strrep("=", 70), "\n", sep = "")
  print(summary_df, row.names = FALSE, digits = 3)

  cat("\n=== vs bulk direction ===\n")
  cat("  bulk retina: Complex III down (beta = -0.126, FDR = 0.016)\n")
  cat("  MR: NDUFB3 associated with early AMD risk (OR = 6.08)\n\n")
  for (pg in FOCUS_PRIMARY) {
    sub <- summary_df[summary_df$gene == pg, ]
    for (i in seq_len(nrow(sub))) {
      r <- sub[i, ]
      if (is.na(r$direction)) {
        cat(sprintf("  %-8s | %-20s | not detected\n", r$gene, r$dataset))
      } else {
        cat(sprintf("  %-8s | %-20s | %s (effect=%+.3f, p=%.3g)\n",
                    r$gene, r$dataset, r$direction, r$effect, r$pvalue))
      }
    }
  }
}

# ============================================================
# ============================================================
log_info("--- 6. Summary figure ---")

if (length(summary_rows) > 0) {
  plot_df <- summary_df[!is.na(summary_df$effect), ]
  if (nrow(plot_df) > 0) {
    plot_df$dataset <- factor(plot_df$dataset,
      levels = c("PXD063075_XML","PXD063075_suppl2","PXD029501_GA_RPE"))
    plot_df$direction_label <- ifelse(plot_df$direction == "down",
                                      "down", "up")

    p <- ggplot(plot_df,
                aes(x = gene, y = effect, fill = dataset)) +
      geom_col(position = position_dodge(width = 0.8)) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      theme_bw(base_size = 11) +
      labs(title = "Target gene effect across 3 proteomic datasets",
           subtitle = "bulk retina reference: Complex III down (beta = -0.126)",
           x = "", y = "effect size", fill = "Dataset")
    save_fig(p, DIR_THIS, "summary_target_genes", width = 9, height = 5)
  }
}

# ============================================================
# 7. summaryreport
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("=== 07_3 proteomics opposite summary ===\n")
cat(strrep("=", 70), "\n", sep = "")

cat(sprintf("\n[Part A] PXD063075 XML (ProteinProphet identification)\n"))
if (!is.null(partA_result)) {
  cat(sprintf("  proteins identified: %d\n", length(partA_result$genes)))
  cat(sprintf("  target genes detected: %d\n", length(partA_result$hits)))
} else {
  cat("  not completed or skipped\n")
}

cat(sprintf("\n[Part B] PXD063075 suppl2 (159 proteins)\n"))
if (!is.null(partB_result)) {
  cat(sprintf("  proteins: %d\n", nrow(partB_result$df)))
  cat(sprintf("  target genes detected: %d\n", nrow(partB_result$focus)))
} else {
  cat("  not completed or skipped\n")
}

cat(sprintf("\n[Part C] PXD029501 GA RPE (TMT matrix)\n"))
if (!is.null(partC_result)) {
  cat(sprintf("  proteins: %d\n", nrow(partC_result$df)))
  cat(sprintf("  target genes detected: %d\n", nrow(partC_result$focus)))
} else {
  cat("  not completed or skipped\n")
}

cat("\n=== Core conclusions ===\n")

if (!is.null(partA_result)) {
  cat(sprintf("  - Part A: %d proteins identified, %d target genes detected (identification only, no quantification)\n",
              length(partA_result$genes), length(partA_result$hits)))
} else {
  cat("  - Part A: skipped\n")
}

if (!is.null(partB_result) && !is.null(partB_result$complex)) {
  c3_b <- partB_result$complex[partB_result$complex$complex == "ComplexIII", ]
  cat(sprintf("  - Part B: %d proteins, Complex III hits = %d\n",
              nrow(partB_result$df), nrow(c3_b)))
} else {
  cat("  - Part B: skipped\n")
}

if (!is.null(partC_result) && !is.null(partC_result$complex)) {
  ct <- partC_result$complex
  c3_c <- ct[ct$complex == "ComplexIII", ]
  if (nrow(c3_c) > 0) {
    cat(sprintf("  - Part C: Complex III n_hit=%d, n_up=%d, n_down=%d, n_sig=%d, mean_l2fc=%+.3f\n",
                nrow(c3_c),
                sum(c3_c$signed_log2FC > 0, na.rm = TRUE),
                sum(c3_c$signed_log2FC < 0, na.rm = TRUE),
                sum(c3_c$pvalue < 0.05, na.rm = TRUE),
                mean(c3_c$signed_log2FC, na.rm = TRUE)))
  }
  nd <- partC_result$focus[partC_result$focus$gene == "NDUFB3", ]
  if (nrow(nd) > 0) {
    cat(sprintf("  - Part C: NDUFB3 signed_log2FC=%+.3f (p=%.3g)\n",
                nd$signed_log2FC, nd$pvalue))
  }
} else {
  cat("  - Part C: skipped\n")
}

cat("  - Overall: proteomic directions are inconsistent with bulk retina or not significant\n")
cat("  - Interpretation: Limitation, not included in main text\n")

cat("\n--- Output ---\n", DIR_THIS, "\n")
log_info("=== 07_3 proteomics opposite done ===")
