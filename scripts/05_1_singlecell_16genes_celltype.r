# ============================================================
# 05_1_singlecell_16genes_celltype.r
# ============================================================

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(Matrix); library(dplyr); library(tidyr); library(ggplot2)
  library(pheatmap)
})

DIR_THIS <- setup_output_dir("05_1_singlecell_16genes_celltype")
log_info("=== 05_1_singlecell_16genes_celltype.r start ===")

DIR_SC <- "E:/retina_tools/raw/GSE221042"
mtx <- readMM(file.path(DIR_SC, "GSE221042_data.mtx.gz"))
mtx <- as(mtx, "CsparseMatrix")

genes <- read.csv(gzfile(file.path(DIR_SC, "GSE221042_genes.csv.gz")),
                   stringsAsFactors = FALSE)
meta <- read.csv(gzfile(file.path(DIR_SC, "GSE221042_sample_annotations.csv.gz")),
                  row.names = 1, stringsAsFactors = FALSE)

log_info("mtx: {nrow(mtx)} × {ncol(mtx)}")
log_info("genes 列名: {paste(colnames(genes), collapse=', ')}")
log_info("meta 列名: {paste(colnames(meta), collapse=', ')}")

if (nrow(mtx) == nrow(genes)) {
  rownames(mtx) <- toupper(trimws(genes[[1]]))
  colnames(mtx) <- rownames(meta)
} else {
  mtx <- t(mtx)
  rownames(mtx) <- toupper(trimws(genes[[1]]))
  colnames(mtx) <- rownames(meta)
}

# ---------- 2. QC ----------
libsize <- Matrix::colSums(mtx)
n_gene <- Matrix::colSums(mtx > 0)
mt_idx <- grep("^MT-", rownames(mtx))
mt_pct <- if (length(mt_idx) > 0) {
  Matrix::colSums(mtx[mt_idx, , drop = FALSE]) / libsize * 100
} else rep(0, ncol(mtx))

keep <- libsize > 1000 & n_gene > 200 & n_gene < 8000 & mt_pct < 20
mtx <- mtx[, keep]
meta <- meta[keep, ]
log_info("QC 后: {ncol(mtx)} 细胞")

# ---------- 3. normalization ----------
libsize2 <- Matrix::colSums(mtx)
norm <- log1p(t(t(mtx) / libsize2 * 1e4))

celltype_in_meta <- grep("celltype|cell_type|CellType|author_cell_type",
                          colnames(meta), value = TRUE, ignore.case = TRUE)[1]

if (!is.na(celltype_in_meta)) {
  log_info("使用 meta 里的 celltype 列: {celltype_in_meta}")
  meta$cell_type <- meta[[celltype_in_meta]]
} else {
  log_info("meta 无 celltype 列，从 01_4 的 markers 打分")

  markers <- read_upstream("01_4_prepare_celltype_markers", "data",
                            "celltype_markers.rds")
  log_info("01_4 markers: {length(markers)} 细胞类型")

  scores <- do.call(cbind, lapply(names(markers), function(ct) {
    g <- intersect(markers[[ct]], rownames(norm))
    if (length(g) < 5) return(rep(0, ncol(norm)))
    Matrix::colMeans(norm[g, , drop = FALSE])
  }))
  colnames(scores) <- names(markers)

  max_score <- apply(scores, 1, max)
  best_type <- colnames(scores)[apply(scores, 1, which.max)]
  meta$cell_type <- ifelse(max_score > 0.1, best_type, "Unknown")
}

log_info("cell_type 分布:")
print(table(meta$cell_type))

sig16_file <- file.path(DIR_RESULT, "02_1_differential_expression_main",
                         "tables", "sig16_genes.csv")
sig16 <- read.csv(sig16_file, stringsAsFactors = FALSE)$gene
sig16_in <- intersect(sig16, rownames(norm))
log_info("16 基因匹配: {length(sig16_in)} / 16")

spec_tbl <- data.frame()
for (g in sig16_in) {
  v <- as.numeric(norm[g, ])
  by_ct <- data.frame(cell_type = meta$cell_type, expr = v) %>%
    group_by(cell_type) %>%
    summarise(mean_expr = mean(expr), .groups = "drop") %>%
    arrange(desc(mean_expr))
  spec_tbl <- rbind(spec_tbl, data.frame(
    gene = g,
    top_celltype = by_ct$cell_type[1],
    top_mean_expr = round(by_ct$mean_expr[1], 4),
    pct_cells_expressing = round(100 * mean(v > 0), 1),
    stringsAsFactors = FALSE))
}
save_table(spec_tbl, DIR_THIS, "16genes_celltype_specificity.csv")

muller_labels <- grep("muller|Müller|^MG$|Muller_glia|Muller_glial",
                       unique(meta$cell_type), ignore.case = TRUE, value = TRUE)
log_info("Müller 标签: {paste(muller_labels, collapse=', ')}")

if (length(muller_labels) > 0) {
  muller_idx <- which(meta$cell_type %in% muller_labels)

  donor_col <- grep("donor|sample|orig.ident|patient",
                     colnames(meta), value = TRUE, ignore.case = TRUE)[1]
  cond_col <- grep("condition|disease|group|status",
                    colnames(meta), value = TRUE, ignore.case = TRUE)[1]

  if (!is.na(donor_col) && !is.na(cond_col)) {
    log_info("donor: {donor_col}, condition: {cond_col}")

    # Pseudo-bulk
    muller_meta <- meta[muller_idx, ]
    donor_ct <- table(muller_meta[[donor_col]])
    valid_donors <- names(donor_ct)[donor_ct >= 10]

    pb <- data.frame()
    for (d in valid_donors) {
      idx <- muller_idx[meta[[donor_col]][muller_idx] == d]
      mean_expr <- rowMeans(norm[sig16_in, idx, drop = FALSE])
      pb <- rbind(pb, data.frame(
        donor = d,
        condition = meta[[cond_col]][idx[1]],
        t(mean_expr),
        stringsAsFactors = FALSE))
    }

    # dAMD vs Healthy
    cond_vals <- unique(pb$condition)
    amd_label <- grep("AMD|dAMD|disease", cond_vals, ignore.case = TRUE, value = TRUE)[1]
    ctrl_label <- grep("Healthy|control|normal|Ctrl", cond_vals, ignore.case = TRUE, value = TRUE)[1]

    if (!is.na(amd_label) && !is.na(ctrl_label)) {
      test_res <- data.frame()
      for (g in sig16_in) {
        x <- pb[[g]][pb$condition == amd_label]
        y <- pb[[g]][pb$condition == ctrl_label]
        if (length(x) < 3 || length(y) < 3) next
        tt <- tryCatch(t.test(x, y), error = function(e) NULL)
        if (is.null(tt)) next
        test_res <- rbind(test_res, data.frame(
          gene = g,
          amd_mean = round(mean(x), 4),
          ctrl_mean = round(mean(y), 4),
          p = tt$p.value,
          stringsAsFactors = FALSE))
      }
      if (nrow(test_res) > 0) {
        test_res$FDR <- p.adjust(test_res$p, "BH")
        test_res <- test_res[order(test_res$p), ]
        save_table(test_res, DIR_THIS, "muller_16genes_dAMD_vs_healthy.csv")
      }
    }
  }
}

if ("S100B" %in% sig16_in && "C3" %in% rownames(norm) &&
    length(muller_labels) > 0) {
  muller_idx <- which(meta$cell_type %in% muller_labels)
  s100b_mg <- as.numeric(norm["S100B", muller_idx])
  c3_mg <- as.numeric(norm["C3", muller_idx])
  r <- cor(s100b_mg, c3_mg)
  log_info("S100B vs C3 in Müller: r={round(r, 3)}")
  save_table(data.frame(
    pair = "S100B_C3", celltype = "MG",
    r_singlecell = round(r, 3), n_cells = length(muller_idx),
    stringsAsFactors = FALSE),
    DIR_THIS, "s100b_vs_c3_correlation.csv")
}

# ---------- 9. summary ----------
cat("\n", strrep("=", 70), "\n", sep = "")
cat("=== 05_1 单细胞 16 基因验证 汇总 ===\n")
cat(strrep("=", 70), "\n", sep = "")
cat("QC 后细胞:", ncol(norm), "\n")
cat("16 基因匹配:", length(sig16_in), "/ 16\n\n")
cat("--- 16 基因细胞来源 ---\n")
print(spec_tbl, row.names = FALSE)

log_info("=== 05_1 done ===")
