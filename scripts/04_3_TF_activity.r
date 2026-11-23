#!/usr/bin/env Rscript
# =============================================================================
# 04_3_TF_activity.R
#
#   - 01_1_prepare_bulk_data.r  → expr_corrected.rds, clinical.rds, gene_annotation.rds
#   - 01_2_prepare_gene_sets.r  → gene_sets.rds
# =============================================================================

rm(list = ls()); gc()
setwd("E:/retina_tools/aging_rgc_axon/AMD")
set.seed(2026)

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(data.table)
  library(ggplot2); library(logger); library(decoupleR)
  library(limma); library(splines); library(readxl)
  library(tibble)
})

output_dir   <- "results/04_3_TF_activity"
data_dir     <- file.path(output_dir, "data")
tables_dir   <- file.path(output_dir, "tables")
figures_dir  <- file.path(output_dir, "figures")
for_04_4_dir <- file.path(output_dir, "for_04_4")
for (d in c(data_dir, tables_dir, figures_dir, for_04_4_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

log_file <- file.path(output_dir, "04_3_TF_activity.log")
log_appender(appender_tee(log_file))
log_layout(layout_glue_generator(format = '{time} [{level}] {msg}'))
log_info("==================== 04_3_TF_activity 启动 ====================")
log_info("方法框架：limma multi (7 对比) + poly3_F + 年龄双模型分类（与 02_1 一致）")

# =============================================================================
# =============================================================================
DIR_01_1  <- "results/01_1_prepare_bulk_data"
DIR_01_2  <- "results/01_2_prepare_gene_sets"

FP_DIR    <- "E:/retina_tools/pubdata/retina_fingerprint"
PROP_XLSX <- "E:/retina_tools/raw/GSE265801/41588_2025_2454_MOESM4_ESM.xlsx"

# =============================================================================
# =============================================================================
GENE_SETS_PATH <- file.path(DIR_01_2, "data", "gene_sets.rds")
if (!file.exists(GENE_SETS_PATH)) {
  stop("未找到 gene_sets.rds: ", GENE_SETS_PATH)
}
GENE_SETS <- readRDS(GENE_SETS_PATH)
log_info("GENE_SETS: {length(GENE_SETS)} 个基因集")

# =============================================================================
# 2. loadexpressionmatrix
# =============================================================================
expr_path <- file.path(DIR_01_1, "data", "expr_corrected.rds")
anno_path <- file.path(DIR_01_1, "data", "gene_annotation.rds")
clin_path <- file.path(DIR_01_1, "data", "clinical.rds")
for (p in c(expr_path, anno_path, clin_path)) {
  if (!file.exists(p)) stop("未找到: ", p)
}
expr <- readRDS(expr_path)
anno <- readRDS(anno_path)
clin <- readRDS(clin_path)
log_info("expr_corrected: {nrow(expr)} × {ncol(expr)}")

n_na_expr <- sum(is.na(expr))
log_info("expr 总 NA 数: {n_na_expr}")
if (n_na_expr > 0) {
  keep_genes <- rowSums(is.na(expr)) == 0
  expr <- expr[keep_genes, , drop = FALSE]
  log_info("移除 NA 基因后: {nrow(expr)} × {ncol(expr)}")
}
if (any(!is.finite(expr))) {
  log_warn("expr 有非有限值，将被替换为 0")
  expr[!is.finite(expr)] <- 0
}

anno_sub <- anno[anno$ensembl %in% rownames(expr), ]
anno_sub <- anno_sub[!is.na(anno_sub$symbol), ]
expr <- expr[anno_sub$ensembl, , drop = FALSE]
expr <- expr[!duplicated(anno_sub$symbol), , drop = FALSE]
rownames(expr) <- anno_sub$symbol[!duplicated(anno_sub$symbol)]
log_info("expr 转 SYMBOL: {nrow(expr)} × {ncol(expr)}")

if (sum(is.na(expr)) > 0) {
  keep_genes <- rowSums(is.na(expr)) == 0
  expr <- expr[keep_genes, , drop = FALSE]
  log_info("二次清理 NA 后: {nrow(expr)} × {ncol(expr)}")
}

common_samples <- intersect(colnames(expr), rownames(clin))
expr <- expr[, common_samples, drop = FALSE]
clin <- clin[common_samples, , drop = FALSE]
stopifnot(identical(colnames(expr), rownames(clin)))

clin$grade <- as.integer(clin$grade)
log_info("样本数: {nrow(clin)}")
log_info("grade 分布: {paste(names(table(clin$grade)), table(clin$grade), sep='=', collapse=' ')}")

# =============================================================================
# =============================================================================
log_info("========== 3. 加载 celltype 特异 GRN（最原始输入） ==========")

if (!dir.exists(FP_DIR)) stop("retina_fingerprint 目录不存在: ", FP_DIR)
xlsx_files <- list.files(FP_DIR, pattern = "\\.xlsx$",
                          full.names = TRUE, recursive = TRUE)
log_info("retina_fingerprint xlsx: {length(xlsx_files)} 个")

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
  col_weight <- guess_col(nms, c("^estimate$", "^weight$", "^w_tissue$",
                                  "^w$", "^importance$", "^score$", "^cor$"))
  if (is.na(col_tf) || is.na(col_target)) return(NULL)
  tf_v     <- toupper(trimws(as.character(df[[col_tf]])))
  target_v <- toupper(trimws(as.character(df[[col_target]])))
  weight_v <- if (!is.na(col_weight)) suppressWarnings(as.numeric(df[[col_weight]])) else NA_real_
  keep <- !is.na(tf_v) & tf_v != "" & !is.na(target_v) & target_v != "" &
          !is.na(weight_v) & is.finite(weight_v)
  if (sum(keep) == 0) return(NULL)
  data.frame(
    tf       = tf_v[keep],
    target   = target_v[keep],
    weight   = weight_v[keep],
    celltype = trimws(sub("\\s*GRN.*$", "", sheet_name, ignore.case = TRUE)),
    stringsAsFactors = FALSE)
}

grn_all <- list()
for (f in xlsx_files) {
  sheets <- tryCatch(readxl::excel_sheets(f), error = function(e) NULL)
  if (is.null(sheets)) next
  grn_sheets <- sheets[grepl("GRN", sheets, ignore.case = TRUE)]
  if (length(grn_sheets) == 0) next
  log_info("  {basename(f)}: {length(grn_sheets)} 个 GRN sheet")
  for (s in grn_sheets) {
    df <- tryCatch(readxl::read_excel(f, sheet = s), error = function(e) NULL)
    if (is.null(df)) next
    parsed <- parse_grn_sheet(df, s)
    if (!is.null(parsed)) {
      grn_all[[length(grn_all) + 1]] <- parsed
      log_info("    {s}: {nrow(parsed)} 边")
    }
  }
}
if (length(grn_all) == 0) stop("没有解析到 GRN")

grn_ct <- do.call(rbind, grn_all)
grn_ct$tf       <- toupper(trimws(grn_ct$tf))
grn_ct$target   <- toupper(trimws(grn_ct$target))
grn_ct$weight   <- as.numeric(grn_ct$weight)
grn_ct$celltype <- trimws(as.character(grn_ct$celltype))

log_info("合并后 GRN: {nrow(grn_ct)} 边, {length(unique(grn_ct$tf))} TF")
log_info("weight 范围: [{round(min(grn_ct$weight, na.rm=TRUE), 4)}, {round(max(grn_ct$weight, na.rm=TRUE), 4)}]")
log_info("weight NA: {sum(is.na(grn_ct$weight))}, weight 非有限: {sum(!is.finite(grn_ct$weight))}")

if (!file.exists(PROP_XLSX)) {
  stop("未找到细胞比例表: ", PROP_XLSX)
}

read_cell_prop <- function(xlsx_path, sheet = "Table S4") {
  sheets_all <- readxl::excel_sheets(xlsx_path)
  if (!sheet %in% sheets_all) {
    cand <- grep("S4|Table\\s*4", sheets_all, ignore.case = TRUE, value = TRUE)
    if (length(cand) > 0) sheet <- cand[1] else return(NULL)
  }
  raw <- tryCatch(
    suppressMessages(readxl::read_excel(xlsx_path, sheet = sheet, col_names = FALSE)),
    error = function(e) NULL)
  if (is.null(raw) || nrow(raw) < 3) return(NULL)
  raw <- as.data.frame(raw, stringsAsFactors = FALSE)

  header_cand <- which(apply(raw, 1, function(r) {
    txt <- paste(as.character(r), collapse = " ")
    grepl("Major\\s*class", txt, ignore.case = TRUE) &&
      grepl("Proportion", txt, ignore.case = TRUE) &&
      grepl("#?\\s*[Cc]ells", txt, ignore.case = TRUE)
  }))
  if (length(header_cand) == 0) {
    log_warn("  [read_cell_prop] 未找到表头行")
    return(NULL)
  }
  header_row <- header_cand[length(header_cand)]
  log_info("  [read_cell_prop] 表头行 = {header_row}")

  df <- raw[(header_row + 1):nrow(raw), , drop = FALSE]
  new_names <- trimws(as.character(raw[header_row, ]))
  empty_idx <- is.na(new_names) | new_names == ""
  new_names[empty_idx] <- paste0("V", which(empty_idx))
  colnames(df) <- make.unique(new_names)
  df <- as.data.frame(df, stringsAsFactors = FALSE)

  is_empty <- apply(df, 1, function(r) all(is.na(r) | trimws(as.character(r)) == ""))
  df <- df[!is_empty, , drop = FALSE]

  ct_col   <- grep("Major\\s*class|cell.?type", colnames(df), ignore.case = TRUE, value = TRUE)[1]
  prop_col <- grep("^Proportion$|Proportion", colnames(df), ignore.case = TRUE, value = TRUE)[1]
  num_col  <- grep("^#?\\s*[Cc]ells$|^cells$|^count$|^number$", colnames(df),
                   ignore.case = TRUE, value = TRUE)[1]

  log_info("  [read_cell_prop] ct_col='{ct_col}', prop_col='{prop_col}', num_col='{num_col}'")
  if (is.na(ct_col) || is.na(prop_col)) {
    log_warn("  [read_cell_prop] 列名识别失败，列名: {paste(colnames(df), collapse=', ')}")
    return(NULL)
  }

  df <- df[!is.na(df[[ct_col]]) & trimws(df[[ct_col]]) != "", , drop = FALSE]
  df <- df[!grepl("^\\s*(major\\s*class|total|sum)\\s*$",
                  trimws(df[[ct_col]]), ignore.case = TRUE), , drop = FALSE]

  cell_type <- trimws(df[[ct_col]])

  n_cells <- if (!is.na(num_col)) {
    suppressWarnings(as.numeric(gsub("[^0-9]", "", as.character(df[[num_col]]))))
  } else NA_real_

  prop_raw <- df[[prop_col]]
  if (is.numeric(prop_raw)) {
    prop_num <- prop_raw
  } else {
    prop_raw <- trimws(as.character(prop_raw))
    is_pct <- grepl("%", prop_raw)
    prop_num <- suppressWarnings(as.numeric(gsub("[^0-9eE.+-]", "", prop_raw)))
    prop_num[is_pct] <- prop_num[is_pct] / 100
    if (!any(is_pct) && all(prop_num > 1, na.rm = TRUE)) prop_num <- prop_num / 100
  }

  out <- data.frame(cell_type = cell_type, n_cells = n_cells,
                    proportion = prop_num, stringsAsFactors = FALSE)

  if (all(is.na(out$proportion)) && !all(is.na(out$n_cells))) {
    out$proportion <- out$n_cells / sum(out$n_cells, na.rm = TRUE)
  }

  if (all(is.na(out$proportion))) return(NULL)
  out$proportion <- out$proportion / sum(out$proportion)
  out
}

cell_prop <- read_cell_prop(PROP_XLSX)

if (is.null(cell_prop) || nrow(cell_prop) == 0) {
  stop("Excel 读取失败，cell_prop 为 NULL。请检查 ", PROP_XLSX)
}
if (all(is.na(cell_prop$proportion))) {
  stop("Excel 读取成功但 proportion 全 NA，请检查表结构")
}
if (abs(sum(cell_prop$proportion, na.rm = TRUE) - 1) > 0.01) {
  stop("proportion 和 = ", sum(cell_prop$proportion, na.rm = TRUE), "，不等于 1")
}

log_info("细胞比例: {nrow(cell_prop)} 个细胞类型")
print(cell_prop, row.names = FALSE)

map_celltype <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub("[[:space:]]", "", x)
  out <- rep(NA_character_, length(x))
  out[grepl("^rod",          x)] <- "rod"
  out[grepl("^cone",         x)] <- "cone"
  out[grepl("^bc|bipolar",   x)] <- "bc"
  out[grepl("^ac|amacrine",  x)] <- "ac"
  out[grepl("^hc|horizontal", x)] <- "hc"
  out[grepl("^rgc|ganglion", x)] <- "rgc"
  out[grepl("^mg|muller|müller", x)] <- "mg"
  out[grepl("^astro",        x)] <- "astrocyte"
  out[grepl("^micro",        x)] <- "microglia"
  out[grepl("^rpe",          x)] <- "rpe"
  out[grepl("^prpc",         x)] <- "prpc"
  out[grepl("^nrpc",         x)] <- "nrpc"
  out
}

cell_prop$ct_norm <- map_celltype(cell_prop$cell_type)
grn_ct$ct_norm    <- map_celltype(grn_ct$celltype)

log_info("cell_prop ct_norm: {paste(cell_prop$ct_norm, collapse=', ')}")
log_info("grn_ct ct_norm 唯一值: {paste(sort(unique(grn_ct$ct_norm)), collapse=', ')}")

common_ct <- intersect(unique(grn_ct$ct_norm), cell_prop$ct_norm)
common_ct <- common_ct[!is.na(common_ct)]
log_info("匹配到的 cell type: {paste(common_ct, collapse=', ')}")
if (length(common_ct) == 0) {
  stop("GRN 的 celltype 与细胞比例表无法匹配")
}

grn_ct <- grn_ct[grn_ct$ct_norm %in% common_ct, , drop = FALSE]
log_info("过滤后 grn_ct: {nrow(grn_ct)} 边, {length(unique(grn_ct$tf))} TF")

prop_df <- data.frame(
  ct_norm    = cell_prop$ct_norm,
  proportion = cell_prop$proportion,
  stringsAsFactors = FALSE
)
prop_df <- prop_df[!duplicated(prop_df$ct_norm), , drop = FALSE]

grn_ct <- merge(grn_ct, prop_df, by = "ct_norm", all.x = TRUE)
log_info("merge 后 grn_ct: {nrow(grn_ct)} 边")
log_info("grn_ct$proportion: min={round(min(grn_ct$proportion, na.rm=TRUE), 4)}, max={round(max(grn_ct$proportion, na.rm=TRUE), 4)}, NA={sum(is.na(grn_ct$proportion))}")

if (any(is.na(grn_ct$proportion))) {
  stop("proportion 匹配失败")
}

denom <- sum(unique(prop_df$proportion[prop_df$ct_norm %in% common_ct]), na.rm = TRUE)
log_info("denom = {round(denom, 4)}")
if (denom <= 0) stop("denom <= 0")

grn_ct$prop_norm <- grn_ct$proportion / denom
log_info("grn_ct$prop_norm: min={round(min(grn_ct$prop_norm, na.rm=TRUE), 4)}, max={round(max(grn_ct$prop_norm, na.rm=TRUE), 4)}")

grn_ct$key <- paste(grn_ct$tf, grn_ct$target, sep = "|||")
split_idx <- split(seq_len(nrow(grn_ct)), grn_ct$key)

w_tissue_vec     <- numeric(length(split_idx))
n_cell_types_vec <- integer(length(split_idx))
cell_types_vec   <- character(length(split_idx))
tf_vec           <- character(length(split_idx))
target_vec       <- character(length(split_idx))

for (i in seq_along(split_idx)) {
  idx <- split_idx[[i]]
  est <- grn_ct$weight[idx]
  pn  <- grn_ct$prop_norm[idx]
  valid <- !is.na(est) & !is.na(pn) & is.finite(est) & is.finite(pn)
  if (sum(valid) == 0) {
    w_tissue_vec[i] <- 0
  } else {
    w_tissue_vec[i] <- sum(est[valid] * pn[valid])
  }
  n_cell_types_vec[i] <- sum(valid)
  cell_types_vec[i]   <- paste(sort(unique(grn_ct$celltype[idx][valid])), collapse = ";")
  tf_vec[i]           <- grn_ct$tf[idx[1]]
  target_vec[i]       <- grn_ct$target[idx[1]]
}

grn_tissue <- data.frame(
  tf = tf_vec, target = target_vec, w_tissue = w_tissue_vec,
  n_cell_types = n_cell_types_vec, cell_types = cell_types_vec,
  stringsAsFactors = FALSE
)
grn_tissue <- grn_tissue[order(-abs(grn_tissue$w_tissue)), ]

log_info("组织水平 GRN: {nrow(grn_tissue)} 边, {length(unique(grn_tissue$tf))} TF")
log_info("w_tissue 范围: [{round(min(grn_tissue$w_tissue), 4)}, {round(max(grn_tissue$w_tissue), 4)}]")
log_info("w_tissue 正边: {sum(grn_tissue$w_tissue > 0)}, 负边: {sum(grn_tissue$w_tissue < 0)}, 零边: {sum(grn_tissue$w_tissue == 0)}")

if (all(grn_tissue$w_tissue == 0)) {
  stop("w_tissue 全是 0")
}

net <- data.frame(
  source = grn_tissue$tf, target = grn_tissue$target,
  weight = grn_tissue$w_tissue, stringsAsFactors = FALSE
)

genes_in_mat <- rownames(expr)
n_before <- nrow(net)
net <- net[net$target %in% genes_in_mat, ]
log_info("组织水平 GRN 过滤后: {n_before} → {nrow(net)} 边, {length(unique(net$source))} TF")
log_info("net$weight 范围: [{round(min(net$weight), 4)}, {round(max(net$weight), 4)}]")

if (all(net$weight == 0) || sum(abs(net$weight)) < 1e-10) {
  stop("net$weight 全 0")
}

saveRDS(grn_tissue, file.path(data_dir, "GRN_tissue_weighted.rds"))
saveRDS(net,        file.path(data_dir, "net_filtered.rds"))

# =============================================================================
# =============================================================================
log_info("========== TF 活性（decoupleR::run_wmean） ==========")

w_max <- max(abs(net$weight), na.rm = TRUE)
log_info("net$weight max|w| = {round(w_max, 4)}")
if (w_max <= 0) stop("w_max <= 0")
net$weight <- net$weight / w_max
log_info("net$weight 缩放后范围: [{round(min(net$weight), 4)}, {round(max(net$weight), 4)}]")

tf_counts <- table(net$source)
log_info("TF 靶基因数: min={min(tf_counts)}, median={as.integer(median(tf_counts))}, max={max(tf_counts)}")
log_info("TF 靶基因数 >= 5: {sum(tf_counts >= 5)}")

tf_activity_long <- run_wmean(
  mat = expr, network = net,
  .source = "source", .target = "target", .mor = "weight",
  times = 100, minsize = 5)

log_info("run_wmean: {nrow(tf_activity_long)} 行, {length(unique(tf_activity_long$source))} TF")
log_info("run_wmean score 范围: [{round(min(tf_activity_long$score, na.rm=TRUE), 4)}, {round(max(tf_activity_long$score, na.rm=TRUE), 4)}]")
log_info("run_wmean score NA 数: {sum(is.na(tf_activity_long$score))}")

if (sum(!is.na(tf_activity_long$score)) == 0) {
  log_error("run_wmean score 全 NA！前 6 行：")
  print(head(tf_activity_long, 6))
  stop("run_wmean 失败")
}

if ("condition" %in% colnames(tf_activity_long) &&
    !"sample" %in% colnames(tf_activity_long)) {
  colnames(tf_activity_long)[colnames(tf_activity_long) == "condition"] <- "sample"
}

stat_val <- unique(tf_activity_long$statistic)
log_info("run_wmean statistic 值: {paste(stat_val, collapse=', ')}")
if ("wmean" %in% stat_val) {
  tf_activity_long <- tf_activity_long[tf_activity_long$statistic == "wmean", ]
  log_info("筛选 statistic == 'wmean' 后: {nrow(tf_activity_long)} 行")
}

tf_score_mat <- tf_activity_long %>%
  dplyr::select(sample, source, score) %>%
  pivot_wider(names_from = source, values_from = score,
              values_fn = list(score = mean)) %>%
  tibble::column_to_rownames("sample") %>% as.matrix()
storage.mode(tf_score_mat) <- "numeric"

log_info("tf_score_mat: {nrow(tf_score_mat)} × {ncol(tf_score_mat)}")
log_info("tf_score_mat NA: {sum(is.na(tf_score_mat))}, NaN: {sum(is.nan(tf_score_mat))}")

if (any(!is.finite(tf_score_mat))) {
  log_warn("tf_score_mat 有非有限值，替换为 0")
  tf_score_mat[!is.finite(tf_score_mat)] <- 0
}

common_samples <- intersect(rownames(tf_score_mat), rownames(clin))
tf_score_mat <- tf_score_mat[common_samples, , drop = FALSE]
clin <- clin[common_samples, , drop = FALSE]
stopifnot(identical(rownames(tf_score_mat), rownames(clin)))

tf_activity_z <- scale(tf_score_mat)
storage.mode(tf_activity_z) <- "numeric"
tf_list <- colnames(tf_activity_z)
log_info("有效 TF: {length(tf_list)}")
log_info("tf_activity_z NA: {sum(is.na(tf_activity_z))}")

if (any(is.na(tf_activity_z))) {
  bad_cols <- apply(tf_activity_z, 2, function(x) any(is.na(x)))
  tf_activity_z <- tf_activity_z[, !bad_cols, drop = FALSE]
  tf_score_mat  <- tf_score_mat[,  !bad_cols, drop = FALSE]
  tf_list <- colnames(tf_activity_z)
  log_info("移除 NA 列后 TF 数: {length(tf_list)}")
}

if (length(tf_list) == 0) stop("tf_activity_z 没有任何有效 TF")

saveRDS(tf_score_mat,  file.path(data_dir, "TF_activity_raw.rds"))
saveRDS(tf_activity_z, file.path(data_dir, "TF_activity_zscore.rds"))

# =============================================================================
# =============================================================================
log_info("========== 5. 协变量 ==========")

stopifnot(identical(rownames(tf_activity_z), rownames(clin)))

sp <- splines::ns(clin$age_n, df = 3)

design_df <- data.frame(
  grade  = factor(clin$grade),
  age_s1 = sp[,1],
  age_s2 = sp[,2],
  age_s3 = sp[,3],
  sex_M  = as.integer(clin$sex == "Male"),
  rin    = clin$rin_n,
  pmi    = clin$pmi_n,
  row.names = rownames(clin)
)

keep_idx <- complete.cases(design_df)
design_df <- design_df[keep_idx, , drop = FALSE]
tf_act_sub <- tf_activity_z[keep_idx, , drop = FALSE]

stopifnot(identical(rownames(design_df), rownames(tf_act_sub)))
log_info("样本数: {nrow(design_df)}, TFs: {ncol(tf_act_sub)}")

n_na_rows <- sum(apply(tf_act_sub, 1, function(x) all(is.na(x))))
if (n_na_rows > 0) {
  keep2 <- !apply(tf_act_sub, 1, function(x) all(is.na(x)))
  tf_act_sub <- tf_act_sub[keep2, , drop = FALSE]
  design_df  <- design_df[keep2, , drop = FALSE]
  log_info("移除 NA 行后样本数: {nrow(design_df)}")
}

if (nrow(design_df) < 20) stop("有效样本数 < 20")

# =============================================================================
# =============================================================================
log_info("========== 6. limma multi (7 对比) ==========")

design_multi <- model.matrix(
  ~ 0 + grade + age_s1 + age_s2 + age_s3 + sex_M + rin + pmi,
  data = design_df)
grade_idx <- grep("^grade", colnames(design_multi))
colnames(design_multi)[grade_idx] <- paste0("g", 1:4)

log_info("design_multi: {nrow(design_multi)} x {ncol(design_multi)}")
log_info("design_multi colnames: {paste(colnames(design_multi), collapse=', ')}")

cm <- limma::makeContrasts(
  g2_vs_g1      = g2 - g1,
  g3_vs_g2      = g3 - g2,
  g4_vs_g3      = g4 - g3,
  g4_vs_g1      = g4 - g1,
  anyAMD_vs_g1  = (g2 + g3 + g4)/3 - g1,
  late_vs_early = (g3 + g4)/2 - (g1 + g2)/2,
  g4_vs_rest    = g4 - (g1 + g2 + g3)/3,
  levels = design_multi
)
contrast_names <- colnames(cm)
log_info("7 对比: {paste(contrast_names, collapse=', ')}")

fit_multi <- eBayes(contrasts.fit(lmFit(t(tf_act_sub), design_multi), cm))

multi_res <- data.frame()
for (ct in contrast_names) {
  tt <- topTable(fit_multi, coef = ct, number = Inf, sort.by = "none")
  multi_res <- rbind(multi_res, data.frame(
    TF = rownames(tt), contrast = ct,
    logFC = tt$logFC, t = tt$t, p = tt$P.Value, FDR = tt$adj.P.Val,
    stringsAsFactors = FALSE))
}
write.csv(multi_res, file.path(tables_dir, "limma_multi_full.csv"), row.names = FALSE)
log_info("limma multi: {nrow(multi_res)} 行")

cat("\n=== limma multi 各对比 FDR<0.05 ===\n")
for (ct in contrast_names) {
  n_sig <- sum(multi_res$FDR[multi_res$contrast == ct] < 0.05, na.rm = TRUE)
  cat(sprintf("  %-15s  %d\n", ct, n_sig))
}

# =============================================================================
# 7. poly3_F
# =============================================================================
log_info("========== 7. poly3_F ==========")

poly_mat <- poly(clin$grade, 3)
colnames(poly_mat) <- c("grade_p1", "grade_p2", "grade_p3")
rownames(poly_mat) <- rownames(clin)

design_poly_df <- data.frame(
  grade_p1 = poly_mat[,1],
  grade_p2 = poly_mat[,2],
  grade_p3 = poly_mat[,3],
  age_s1   = sp[,1],
  age_s2   = sp[,2],
  age_s3   = sp[,3],
  sex_M    = as.integer(clin$sex == "Male"),
  rin      = clin$rin_n,
  pmi      = clin$pmi_n,
  row.names = rownames(clin)
)
design_poly_df <- design_poly_df[keep_idx, , drop = FALSE]
tf_act_poly <- tf_activity_z[keep_idx, , drop = FALSE]

stopifnot(identical(rownames(design_poly_df), rownames(tf_act_poly)))

design_poly <- model.matrix(
  ~ grade_p1 + grade_p2 + grade_p3 + age_s1 + age_s2 + age_s3 + sex_M + rin + pmi,
  data = design_poly_df)
log_info("design_poly: {nrow(design_poly)} x {ncol(design_poly)}")

fit_poly <- eBayes(lmFit(t(tf_act_poly), design_poly))

coef_poly <- grep("^grade_p", colnames(design_poly), value = TRUE)
tt_F <- topTable(fit_poly, coef = coef_poly, number = Inf, sort.by = "none")
polyF <- data.frame(
  TF = rownames(tt_F),
  F = tt_F$F, p = tt_F$P.Value, FDR = tt_F$adj.P.Val,
  stringsAsFactors = FALSE)
write.csv(polyF, file.path(tables_dir, "poly3_F_full.csv"), row.names = FALSE)

log_info("poly3_F: FDR<0.05 = {sum(polyF$FDR < 0.05, na.rm=TRUE)}")
log_info("poly3_F: p<0.001   = {sum(polyF$p < 0.001, na.rm=TRUE)}")
log_info("poly3_F: p<0.01    = {sum(polyF$p < 0.01, na.rm=TRUE)}")

# =============================================================================
# =============================================================================
log_info("========== 8. 候选 TF ==========")

sigTF <- polyF$TF[polyF$FDR < 0.05]
sigTF_detail <- polyF[polyF$TF %in% sigTF, ]
sigTF_detail <- sigTF_detail[order(sigTF_detail$FDR), ]
write.csv(sigTF_detail, file.path(tables_dir, "sigTF_poly3_detail.csv"), row.names = FALSE)
write.csv(data.frame(TF = sigTF), file.path(tables_dir, "sigTF_genes.csv"), row.names = FALSE)

log_info("poly3_F FDR<0.05: {length(sigTF)} TF")
if (length(sigTF) > 0) print(sigTF_detail, row.names = FALSE)

age_class_tbl <- data.frame()
for (g in sigTF) {
  df <- data.frame(y = tf_activity_z[, g], grade = clin$grade,
                   age = clin$age_n, sex = clin$sex,
                   rin = clin$rin_n, pmi = clin$pmi_n)
  df <- df[complete.cases(df), ]

  sp_g <- splines::ns(df$age, df = 3)

  fit1 <- lm(y ~ poly(grade, 3) + sp_g + sex + rin + pmi, data = df)
  p_g <- anova(fit1)["poly(grade, 3)", "Pr(>F)"]

  fit2 <- lm(y ~ sp_g + poly(grade, 3) + sex + rin + pmi, data = df)
  p_a <- anova(fit2)["sp_g", "Pr(>F)"]

  r_age <- cor(df$y, df$age, method = "pearson")
  r_age_s <- cor(df$y, df$age, method = "spearman")

  resid_y <- residuals(lm(y ~ factor(grade) + sex + rin + pmi, data = df))
  resid_age <- residuals(lm(age ~ factor(grade) + sex + rin + pmi, data = df))
  r_age_partial <- cor(resid_y, resid_age)

  age_class_tbl <- rbind(age_class_tbl, data.frame(
    TF = g,
    p_grade = p_g, p_age = p_a,
    r_age = round(r_age, 3),
    r_age_spearman = round(r_age_s, 3),
    r_age_partial_grade = round(r_age_partial, 3),
    stringsAsFactors = FALSE))
}

if (nrow(age_class_tbl) > 0) {
  age_class_tbl$type <- "无关"
  age_class_tbl$type[age_class_tbl$p_grade < 0.05 &
                     age_class_tbl$p_age > 0.10] <- "纯疾病"
  age_class_tbl$type[age_class_tbl$p_grade < 0.05 &
                     age_class_tbl$p_age < 0.05] <- "年龄+疾病"
  age_class_tbl$type[age_class_tbl$p_grade > 0.10 &
                     age_class_tbl$p_age < 0.05] <- "纯年龄"

  age_class_tbl <- age_class_tbl[order(age_class_tbl$type,
                                        -abs(age_class_tbl$r_age)), ]
  print(table(age_class_tbl$type))
  print(age_class_tbl, row.names = FALSE)
}

write.csv(age_class_tbl,
          file.path(tables_dir, "sigTF_age_classification.csv"),
          row.names = FALSE)

# =============================================================================
# 9. Top TF
# =============================================================================
top_tf_poly <- head(polyF[order(polyF$FDR), ], 30)
write.csv(top_tf_poly,
          file.path(tables_dir, "top_TF_poly3F.csv"), row.names = FALSE)

tf_targets <- net %>% group_by(source) %>%
  summarise(n_targets = n(), .groups = "drop")
colnames(tf_targets)[1] <- "tf"

tf_cell_map <- grn_ct %>%
  dplyr::select(tf, celltype) %>% distinct() %>%
  group_by(tf) %>%
  summarise(cell_types = paste(sort(unique(celltype)), collapse = "; "),
            n_cell_types = n(), .groups = "drop")

tf_annotation <- data.frame(TF = unique(net$source), stringsAsFactors = FALSE)
tf_annotation$cell_types   <- tf_cell_map$cell_types[match(tf_annotation$TF, tf_cell_map$tf)]
tf_annotation$n_cell_types <- tf_cell_map$n_cell_types[match(tf_annotation$TF, tf_cell_map$tf)]
tf_annotation$n_targets    <- tf_targets$n_targets[match(tf_annotation$TF, tf_targets$tf)]
write.csv(tf_annotation, file.path(tables_dir, "retinal_TF_annotation.csv"),
          row.names = FALSE)

# =============================================================================
# 10. figure
# =============================================================================
log_info("========== 10. 图 ==========")

library(ggplot2)

if (length(sigTF) > 0) {
  sigTF_detail$TF <- factor(
    sigTF_detail$TF,
    levels = rev(sigTF_detail$TF[order(sigTF_detail$F)]))
  p1 <- ggplot(sigTF_detail, aes(x = TF, y = F, fill = F > 0)) +
    geom_col() +
    coord_flip() +
    theme_bw(base_size = 11) +
    scale_fill_manual(values = c("TRUE" = "#C44E52", "FALSE" = "#4C72B0"),
                       guide = "none") +
    labs(x = "", y = "poly3_F F statistic",
         title = "TF activity poly3_F (age=spline)",
         subtitle = sprintf("FDR<0.05: %d TFs", length(sigTF)))
  ggsave(file.path(figures_dir, "sigTF_poly3_F.pdf"), p1,
         width = 6, height = max(4, length(sigTF) * 0.3))
  ggsave(file.path(figures_dir, "sigTF_poly3_F.png"), p1,
         width = 6, height = max(4, length(sigTF) * 0.3), dpi = 200)
}

if (nrow(top_tf_poly) > 0) {
  top_tf_poly$TF <- factor(top_tf_poly$TF,
                            levels = rev(top_tf_poly$TF[order(-top_tf_poly$FDR)]))
  p2 <- ggplot(top_tf_poly, aes(x = TF, y = -log10(p))) +
    geom_col(fill = "#4C72B0") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
    coord_flip() +
    theme_bw(base_size = 11) +
    labs(x = "", y = "-log10(p)",
         title = "Top 30 TF by poly3_F p-value",
         subtitle = sprintf("min FDR = %.3f", min(polyF$FDR, na.rm = TRUE)))
  ggsave(file.path(figures_dir, "top_TF_poly3F.pdf"), p2, width = 6, height = 8)
  ggsave(file.path(figures_dir, "top_TF_poly3F.png"), p2, width = 6, height = 8, dpi = 200)
}

if (nrow(age_class_tbl) > 0) {
  plot_df <- age_class_tbl
  plot_df$neglogp_grade <- -log10(pmax(plot_df$p_grade, 1e-300))

  if (requireNamespace("ggrepel", quietly = TRUE)) {
    library(ggrepel)
    p3 <- ggplot(plot_df, aes(x = r_age, y = neglogp_grade,
                               color = type, label = TF)) +
      geom_point(size = 3) +
      geom_text_repel(size = 3, max.overlaps = 20) +
      geom_vline(xintercept = 0, color = "gray70") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
      theme_bw(base_size = 12) +
      scale_color_manual(values = c(
        "纯疾病" = "#4C72B0", "年龄+疾病" = "#C44E52",
        "无关" = "#B0B0B0", "纯年龄" = "#FDB462")) +
      labs(x = "Pearson r with age", y = "-log10(p_grade)",
           color = "类型", title = "候选 TF：年龄相关 vs 疾病显著性")
    ggsave(file.path(figures_dir, "sigTF_age_vs_grade.pdf"), p3, width = 7, height = 5)
    ggsave(file.path(figures_dir, "sigTF_age_vs_grade.png"), p3, width = 7, height = 5, dpi = 200)
  }
}

# =============================================================================
# =============================================================================
log_info("========== 11. for_04_4 接口 ==========")

top50 <- head(polyF[order(polyF$FDR), ], 50)
top50$rank <- seq_len(nrow(top50))
write.csv(top50, file.path(for_04_4_dir, "TF_top50_poly3F.csv"), row.names = FALSE)

if (length(sigTF) > 0) {
  write.csv(sigTF_detail,
            file.path(for_04_4_dir, "candidate_TFs.csv"), row.names = FALSE)
} else {
  write.csv(data.frame(TF = character(0), F = numeric(0),
                       p = numeric(0), FDR = numeric(0)),
            file.path(for_04_4_dir, "candidate_TFs.csv"), row.names = FALSE)
}

write.csv(tf_annotation,
          file.path(for_04_4_dir, "TF_cell_type_annotation.csv"), row.names = FALSE)
write.csv(top_tf_poly,
          file.path(for_04_4_dir, "top_TF_poly3F.csv"), row.names = FALSE)

min_fdr_poly <- min(polyF$FDR, na.rm = TRUE)
min_fdr_limma <- min(multi_res$FDR, na.rm = TRUE)

run_summary <- data.frame(
  field = c("method", "input_expr", "input_clinical",
            "input_grn", "cell_prop_source",
            "n_samples", "n_genes_expr", "n_grn_edges", "n_TFs",
            "limma_min_FDR", "poly3F_min_FDR",
            "n_candidates_poly3F_FDR05",
            "n_pure_disease", "n_age_disease", "n_pure_age", "n_irrelevant",
            "top_TF_poly3F"),
  value = c("limma multi + poly3_F + age dual-model",
            expr_path, clin_path,
            "retina_fingerprint_xlsx", PROP_XLSX,
            nrow(clin), nrow(expr), nrow(net), length(tf_list),
            round(min_fdr_limma, 3),
            round(min_fdr_poly, 3),
            length(sigTF),
            if (nrow(age_class_tbl) > 0) sum(age_class_tbl$type == "纯疾病") else 0,
            if (nrow(age_class_tbl) > 0) sum(age_class_tbl$type == "年龄+疾病") else 0,
            if (nrow(age_class_tbl) > 0) sum(age_class_tbl$type == "纯年龄") else 0,
            if (nrow(age_class_tbl) > 0) sum(age_class_tbl$type == "无关") else 0,
            paste(head(top_tf_poly$TF, 10), collapse = ", ")),
  stringsAsFactors = FALSE)
write.csv(run_summary, file.path(for_04_4_dir, "run_summary.csv"), row.names = FALSE)

# =============================================================================
# 12. save + report
# =============================================================================
saveRDS(list(clinical = clin, net = net, tf_annotation = tf_annotation,
             grn_tissue = grn_tissue, grn_ct = grn_ct,
             tf_score_mat = tf_score_mat, tf_activity_z = tf_activity_z,
             design_multi = design_multi, contrasts_mat = cm,
             fit_multi = fit_multi, fit_poly = fit_poly,
             multi_res = multi_res, polyF = polyF,
             sigTF = sigTF, sigTF_detail = sigTF_detail,
             age_class_tbl = age_class_tbl,
             top_tf_poly = top_tf_poly,
             params = list(method = "limma multi + poly3_F + age dual-model",
                           n_samples = nrow(clin),
                           n_TFs = length(tf_list),
                           version = "04_3_v6_poly3F")),
        file.path(data_dir, "04_3_TF_activity_results.rds"))

report_lines <- c(
  "# 04_3_TF_activity 报告",
  paste("Date:", Sys.Date()),
  "",
  "## 方法",
  "- 与 02_1 一致：limma multi (7 对比) + poly3_F + 年龄双模型分类",
  "- 输入：01_1 的 expr_corrected.rds（批次校正后）",
  "- TF 活性：decoupleR::run_wmean + w_tissue（细胞比例加权）",
  "- 修复：清理 expr NA；run_wmean 只用 wmean statistic；细胞比例从 Excel 读取",
  "",
  "## 结果",
  paste0("- 有效 TF: ", length(tf_list)),
  paste0("- limma multi 最小 FDR: ", round(min_fdr_limma, 3)),
  paste0("- poly3_F 最小 FDR:     ", round(min_fdr_poly, 3)),
  paste0("- poly3_F FDR<0.05:     ", length(sigTF), " TF"),
  "",
  "## 输出文件",
  "",
  "### tables/",
  "- limma_multi_full.csv",
  "- poly3_F_full.csv",
  "- sigTF_poly3_detail.csv",
  "- sigTF_genes.csv",
  "- sigTF_age_classification.csv",
  "- top_TF_poly3F.csv",
  "- retinal_TF_annotation.csv",
  "",
  "### figures/",
  "- sigTF_poly3_F.pdf/.png",
  "- top_TF_poly3F.pdf/.png",
  "- sigTF_age_vs_grade.pdf/.png"
)
writeLines(report_lines, file.path(output_dir, "report.md"))

log_info("04_3_TF_activity 完成")
log_info("  samples: {nrow(clin)}")
log_info("  TFs: {length(tf_list)}")
log_info("  limma min FDR: {round(min_fdr_limma, 3)}")
log_info("  poly3_F min FDR: {round(min_fdr_poly, 3)}")
log_info("  poly3_F FDR<0.05: {length(sigTF)}")
log_info("End time: {Sys.time()}")
log_info("==================== 04_3_TF_activity 正常退出 ====================")

