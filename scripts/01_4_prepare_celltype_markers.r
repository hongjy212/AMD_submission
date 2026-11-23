# ============================================================
# 01_4_prepare_celltype_markers.r
#   data/celltype_markers.rds
#   tables/celltype_marker_summary.csv
#   tables/celltype_marker_full.csv
# ============================================================

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages(library(readxl))

DIR_THIS <- setup_output_dir("01_4_prepare_celltype_markers")
log_info("=== 01_4_prepare_celltype_markers.r start ===")

# ============================================================
# 0. MOESM6 path
# ============================================================
MOESM6_FILE <- "E:/retina_tools/pubdata/retina_fingerprint/41467_2024_50853_MOESM6_ESM.xlsx"

if (!file.exists(MOESM6_FILE)) {
  alt <- list.files("E:/retina_tools/pubdata/retina_fingerprint",
                    pattern = "MOESM6.*\\.xlsx$",
                    full.names = TRUE, recursive = TRUE)
  if (length(alt) > 0) {
    MOESM6_FILE <- alt[1]
    log_info("自动找到 MOESM6: {MOESM6_FILE}")
  } else {
    log_err("MOESM6 不存在")
    stop("缺 MOESM6")
  }
}
log_info("MOESM6: {MOESM6_FILE}")

# ============================================================
# ============================================================
sheets_all <- readxl::excel_sheets(MOESM6_FILE)
log_info("sheet 数: {length(sheets_all)}")
for (s in sheets_all) log_info("  {s}")

CELLTYPE_CANDIDATES <- c("BC", "Cone", "Rod", "MG", "PRPC", "NRPC",
                          "RGC", "AC", "HC")

celltype_sheets <- intersect(CELLTYPE_CANDIDATES, sheets_all)
log_info("识别到 celltype sheet: {paste(celltype_sheets, collapse=', ')}")

if (length(celltype_sheets) == 0) {
  log_err("MOESM6 里没有识别到 celltype sheet")
  log_err("现有 sheet: {paste(sheets_all, collapse=', ')}")
  stop("无 celltype sheet")
}

# ============================================================
# ============================================================
log_info("--- 读 sheet ---")

find_gene_col <- function(nms) {
  for (cand in c("^gene$", "^symbol$", "^Gene$", "^Symbol$",
                 "^gene_name$", "gene", "symbol")) {
    hit <- grep(cand, nms, ignore.case = TRUE, value = TRUE)
    if (length(hit) > 0) return(hit[1])
  }
  nms[1]
}

marker_list <- list()
for (ct in celltype_sheets) {
  df <- tryCatch(readxl::read_excel(MOESM6_FILE, sheet = ct),
                 error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) {
    log_warn("  {ct}: 空 sheet，跳过")
    next
  }
  nms <- colnames(df)
  gene_col <- find_gene_col(nms)
  genes <- toupper(trimws(as.character(df[[gene_col]])))
  genes <- genes[!is.na(genes) & genes != ""]
  genes <- unique(genes)

  if (length(genes) < 5) {
    log_warn("  {ct}: {length(genes)} 基因 < 5，跳过")
    next
  }

  marker_list[[ct]] <- genes
  log_info("  {ct}: {length(genes)} 基因（列: {gene_col}）")
}

if (length(marker_list) == 0) {
  log_err("没有解析到任何 celltype marker")
  stop("无 marker")
}

# ============================================================
# ============================================================
save_rds(marker_list, DIR_THIS, "celltype_markers.rds")

summary_df <- data.frame(
  celltype  = names(marker_list),
  n_markers = sapply(marker_list, length),
  stringsAsFactors = FALSE)
summary_df <- summary_df[order(-summary_df$n_markers), ]
save_table(summary_df, DIR_THIS, "celltype_marker_summary.csv")

full_list <- do.call(rbind, lapply(names(marker_list), function(ct) {
  data.frame(celltype = ct, gene = marker_list[[ct]],
             stringsAsFactors = FALSE)
}))
save_table(full_list, DIR_THIS, "celltype_marker_full.csv")

# ============================================================
# 4. summary
# ============================================================
cat("\n========================================\n")
cat("=== 01_4_prepare_celltype_markers.r 汇总 ===\n")
cat("========================================\n")
cat("cell types:", length(marker_list), "\n")
cat("每个 celltype marker 数:\n")
print(summary_df, row.names = FALSE)

cat("\n--- 输出位置 ---\n")
cat("目录:", DIR_THIS, "\n")
cat("  data/celltype_markers.rds\n")
cat("  tables/celltype_marker_summary.csv\n")
cat("  tables/celltype_marker_full.csv\n")

log_info("=== 01_4_prepare_celltype_markers.r done ===")

