# ============================================================
# 01_3_prepare_grn_data.r
#   data/GRN_tissue.rds
#   tables/GRN_summary.csv
#   tables/GRN_celltype_matrix.csv
# ============================================================

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(readxl)
})

DIR_THIS <- setup_output_dir("01_3_prepare_grn_data")
log_info("=== 01_3_prepare_grn_data.r start ===")

# ============================================================
# 0. path
# ============================================================
FP_DIR <- "E:/retina_tools/pubdata/retina_fingerprint"
LAMBERT_FILE <- "E:/retina_tools/pubdata/TF_list/Lambert2018_TF_names_v_1.01.txt"

if (!dir.exists(FP_DIR)) {
  log_err("retina_fingerprint 目录不存在: {FP_DIR}")
  stop("目录不存在")
}

xlsx_files <- list.files(FP_DIR, pattern = "\\.xlsx$",
                          full.names = TRUE, recursive = TRUE)
log_info("retina_fingerprint 下 xlsx 数: {length(xlsx_files)}")
for (f in xlsx_files) log_info("  {basename(f)}")
if (length(xlsx_files) == 0) stop("没有 xlsx")

# ============================================================
# 1. Lambert TF
# ============================================================
if (file.exists(LAMBERT_FILE)) {
  tf_raw <- toupper(trimws(readLines(LAMBERT_FILE)))
  TF_REFERENCE <- unique(tf_raw[tf_raw != ""])
  log_info("Lambert TF: {length(TF_REFERENCE)}")
} else {
  TF_REFERENCE <- NULL
  log_warn("Lambert 不存在，跳过 TF 过滤")
}

# ============================================================
# ============================================================
log_info("--- 扫描 xlsx ---")

guess_col <- function(nms, candidates) {
  for (cand in candidates) {
    hit <- grep(cand, nms, ignore.case = TRUE, value = TRUE)
    if (length(hit) > 0) return(hit[1])
  }
  NA_character_
}

parse_grn_sheet <- function(df, sheet_name) {
  nms <- colnames(df)
  if (length(nms) == 0) return(NULL)

  col_tf     <- guess_col(nms, c("^source$", "^tf$", "^tf_name$",
                                  "transcription factor", "^regulator$"))
  col_target <- guess_col(nms, c("^target$", "^gene$", "^target_gene$",
                                  "^symbol$"))
  col_weight <- guess_col(nms, c("^weight$", "^w_tissue$", "^w$",
                                  "^importance$", "^score$"))

  if (is.na(col_tf) || is.na(col_target)) return(NULL)

  tf_v     <- toupper(trimws(as.character(df[[col_tf]])))
  target_v <- toupper(trimws(as.character(df[[col_target]])))
  weight_v <- if (!is.na(col_weight)) suppressWarnings(as.numeric(df[[col_weight]])) else 1

  keep <- !is.na(tf_v) & tf_v != "" & !is.na(target_v) & target_v != ""
  if (sum(keep) == 0) return(NULL)

  data.frame(
    tf       = tf_v[keep],
    target   = target_v[keep],
    weight   = weight_v[keep],
    celltype = sub("\\s*GRN.*$", "", sheet_name, ignore.case = TRUE),
    stringsAsFactors = FALSE)
}

grn_all <- list()
for (f in xlsx_files) {
  sheets <- tryCatch(readxl::excel_sheets(f), error = function(e) NULL)
  if (is.null(sheets)) {
    log_warn("  读 {basename(f)} sheet 失败，跳过")
    next
  }

  grn_sheets <- sheets[grepl("GRN", sheets, ignore.case = TRUE)]
  if (length(grn_sheets) == 0) {
    log_info("  {basename(f)}: 无 GRN sheet（{length(sheets)} 个），跳过")
    next
  }
  log_info("  {basename(f)}: {length(grn_sheets)} 个 GRN sheet")

  for (s in grn_sheets) {
    df <- tryCatch(readxl::read_excel(f, sheet = s),
                   error = function(e) NULL)
    if (is.null(df)) {
      log_warn("    读 sheet '{s}' 失败")
      next
    }
    parsed <- parse_grn_sheet(df, s)
    if (!is.null(parsed)) {
      grn_all[[length(grn_all) + 1]] <- parsed
      log_info("    {s}: {nrow(parsed)} 边")
    } else {
      log_warn("    {s}: 无法识别列，列名: {paste(colnames(df), collapse=' | ')}")
    }
  }
}

if (length(grn_all) == 0) {
  log_err("没有解析到任何 GRN 数据")
  stop("无 GRN")
}

grn <- do.call(rbind, grn_all)
log_info("合并后 GRN: {nrow(grn)} 边")

# ============================================================
# 3. filter
# ============================================================
log_info("--- 过滤 ---")

if (!is.null(TF_REFERENCE)) {
  n_before <- length(unique(grn$tf))
  grn <- grn[grn$tf %in% TF_REFERENCE, , drop = FALSE]
  log_info("TF 过滤: {n_before} → {length(unique(grn$tf))}")
}

grn <- grn[!duplicated(paste(grn$tf, grn$target, grn$celltype, sep = "|")), , drop = FALSE]
log_info("去重后: {nrow(grn)} 边")

ct_counts <- table(grn$celltype)
valid_ct <- names(ct_counts)[ct_counts >= 100]
grn <- grn[grn$celltype %in% valid_ct, , drop = FALSE]
log_info("celltype ≥100 边: {paste(valid_ct, collapse=', ')}")

# ============================================================
# ============================================================
save_rds(grn, DIR_THIS, "GRN_tissue.rds")

summary_df <- grn %>%
  group_by(celltype) %>%
  summarise(
    n_edges  = n(),
    n_tf     = n_distinct(tf),
    n_target = n_distinct(target),
    .groups  = "drop"
  ) %>%
  arrange(desc(n_edges))
save_table(summary_df, DIR_THIS, "GRN_summary.csv")

ct_tf_mat <- grn %>%
  count(celltype, tf) %>%
  tidyr::pivot_wider(names_from = tf, values_from = n, values_fill = 0)
save_table(as.data.frame(ct_tf_mat), DIR_THIS, "GRN_celltype_matrix.csv")

# ============================================================
# 5. summary
# ============================================================
cat("\n========================================\n")
cat("=== 01_3_prepare_grn_data.r 汇总 ===\n")
cat("========================================\n")
cat("GRN 边数:", nrow(grn), "\n")
cat("唯一 TF:", length(unique(grn$tf)), "\n")
cat("唯一靶基因:", length(unique(grn$target)), "\n")
cat("cell type 数:", length(unique(grn$celltype)), "\n\n")
cat("每个 celltype 统计:\n")
print(summary_df, row.names = FALSE)

cat("\n--- 输出位置 ---\n")
cat("目录:", DIR_THIS, "\n")
cat("  data/GRN_tissue.rds\n")
cat("  tables/GRN_summary.csv\n")
cat("  tables/GRN_celltype_matrix.csv\n")

log_info("=== 01_3_prepare_grn_data.r done ===")

