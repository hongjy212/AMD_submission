# ============================================================
# Script name: 07_8a_GSE135092_prep
# Input:  E:/retina_tools/raw/GSE135092/
# Output: E:/retina_tools/aging_rgc_axon/AMD/results/07_8a_GSE135092_prep/data/
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tibble)
  library(readr)
})

raw_dir <- "E:/retina_tools/raw/GSE135092"
out_dir <- "E:/retina_tools/aging_rgc_axon/AMD/results/07_8a_GSE135092_prep"
data_dir <- file.path(out_dir, "data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

cat("==================== 07_8a GSE135092 prep 启动 ====================\n\n")

# ============================================================
# ============================================================
extract_dir <- file.path(raw_dir, "extracted")
if (!dir.exists(extract_dir) ||
    length(list.files(extract_dir, pattern = "\\.tsv\\.gz$")) == 0) {
  dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)
  cat("解压 RAW tar 到:", extract_dir, "\n")
  untar(file.path(raw_dir, "GSE135092_RAW.tar"), exdir = extract_dir)
} else {
  cat("已解压，跳过\n")
}

tsv_files <- list.files(extract_dir, pattern = "\\.tsv\\.gz$",
                        full.names = TRUE)
cat("找到", length(tsv_files), "个 .tsv.gz 文件\n\n")
if (length(tsv_files) == 0) stop("没有 .tsv.gz 文件")

# ============================================================
# ============================================================
cat("========== 2. 检查格式 ==========\n")
f1 <- tsv_files[1]
cat("样例:", basename(f1), "\n")
dt1 <- fread(f1, header = TRUE)
cat("维度:", nrow(dt1), "x", ncol(dt1), "\n")
cat("列名:", paste(colnames(dt1), collapse=", "), "\n")
print(head(dt1, 3))

# ============================================================
# ============================================================
cat("\n========== 3. 批量读入 ==========\n")

value_list <- list()
count_list <- list()

t0 <- Sys.time()
for (i in seq_along(tsv_files)) {
  f <- tsv_files[i]
  dt <- tryCatch(fread(f, header = TRUE, showProgress = FALSE),
                 error = function(e) NULL)
  if (is.null(dt)) next
  if (!all(c("ID_REF", "VALUE", "count") %in% colnames(dt))) next

  sm <- sub("\\.tsv\\.gz$", "", basename(f))
  sm <- sub("_.*$", "", sm)

  val_vec <- as.numeric(dt[["VALUE"]])
  names(val_vec) <- dt[["ID_REF"]]
  value_list[[sm]] <- val_vec

  cnt_int <- as.integer(dt[["count"]])
  names(cnt_int) <- dt[["ID_REF"]]
  count_list[[sm]] <- cnt_int

  if (i %% 50 == 0) cat("  已读", i, "/", length(tsv_files), "\n")
}
cat("读入", length(value_list), "个样本，耗时",
    round(as.numeric(difftime(Sys.time(), t0, units='mins')), 1), "分钟\n")

# ============================================================
# ============================================================
cat("\n========== 4. 合并矩阵 ==========\n")

all_genes <- Reduce(union, lapply(value_list, names))
cat("总基因数:", length(all_genes), "\n")

expr_mat <- matrix(NA_real_,
                   nrow = length(all_genes),
                   ncol = length(value_list),
                   dimnames = list(all_genes, names(value_list)))
for (sm in names(value_list)) {
  v <- value_list[[sm]]
  expr_mat[names(v), sm] <- v
}

all_genes_c <- Reduce(union, lapply(count_list, names))
count_mat <- matrix(0L,
                    nrow = length(all_genes_c),
                    ncol = length(count_list),
                    dimnames = list(all_genes_c, names(count_list)))
for (sm in names(count_list)) {
  v <- count_list[[sm]]
  count_mat[names(v), sm] <- v
}

cat("VALUE 矩阵:", nrow(expr_mat), "x", ncol(expr_mat), "\n")
cat("count 矩阵:", nrow(count_mat), "x", ncol(count_mat), "\n")

# ============================================================
# ============================================================
cat("\n========== 5. 读 series_matrix ==========\n")

sm_path <- file.path(raw_dir, "GSE135092_series_matrix.txt.gz")
sm_lines <- readLines(gzfile(sm_path))
cat("series_matrix 行数:", length(sm_lines), "\n")

parse_meta <- function(lines, prefix) {
  idx <- grep(paste0("^!Sample_", prefix), lines)
  if (length(idx) == 0) return(NULL)
  line <- lines[idx[1]]
  parts <- strsplit(line, "\t")[[1]]
  vals <- parts[-1]
  gsub('^"|"$', "", vals)
}

parse_characteristics <- function(lines) {
  idx <- grep("^!Sample_characteristics_ch1", lines)
  if (length(idx) == 0) return(NULL)
  result <- list()
  for (i in idx) {
    parts <- strsplit(lines[i], "\t")[[1]]
    vals <- gsub('^"|"$', "", parts[-1])
    if (length(vals) == 0) next
    first_val <- vals[1]
    if (grepl("^[a-zA-Z_]+:", first_val)) {
      field_name <- sub(":.*$", "", first_val)
      vals_clean <- sub("^[a-zA-Z_]+:\\s*", "", vals)
      result[[field_name]] <- vals_clean
    }
  }
  result
}

gsm_ids <- parse_meta(sm_lines, "geo_accession")
titles  <- parse_meta(sm_lines, "title")
char_list <- parse_characteristics(sm_lines)

cat("characteristics 字段:", paste(names(char_list), collapse=", "), "\n")

age      <- if ("age"        %in% names(char_list)) char_list[["age"]]        else rep(NA, length(gsm_ids))
amd_stat <- if ("amd_status" %in% names(char_list)) char_list[["amd_status"]] else rep(NA, length(gsm_ids))
samid    <- if ("samid"      %in% names(char_list)) char_list[["samid"]]      else rep(NA, length(gsm_ids))
tissue   <- if ("tissue"     %in% names(char_list)) char_list[["tissue"]]     else rep(NA, length(gsm_ids))

clinical <- data.frame(
  gsm        = gsm_ids,
  title      = titles,
  age        = suppressWarnings(as.numeric(age)),
  amd_status = amd_stat,
  samid      = samid,
  tissue     = tissue,
  stringsAsFactors = FALSE
)

clinical$tissue_type <- ifelse(grepl("^Retina", clinical$tissue), "Retina",
                        ifelse(grepl("^RPE",    clinical$tissue), "RPE", NA))

clinical$region <- ifelse(grepl("non-Macula", clinical$tissue), "non-Macula",
                   ifelse(grepl("Macula",     clinical$tissue), "Macula", NA))

clinical$amd_binary <- factor(
  ifelse(clinical$amd_status == "AMD", "AMD", "Control"),
  levels = c("Control", "AMD"))

cat("\n组织 x 疾病:\n")
print(table(clinical$tissue_type, clinical$amd_binary, useNA = "ifany"))
cat("\n区域 x 疾病:\n")
print(table(clinical$region, clinical$amd_binary, useNA = "ifany"))

# ============================================================
# 6. alignsamples
# ============================================================
cat("\n========== 6. 对齐样本 ==========\n")

common_samples <- base::intersect(colnames(expr_mat), clinical$gsm)
cat("表达矩阵样本:", ncol(expr_mat), "\n")
cat("metadata 样本:", nrow(clinical), "\n")
cat("交集:", length(common_samples), "\n")

if (length(common_samples) == 0) stop("样本名不匹配")

expr_mat  <- expr_mat[, common_samples, drop = FALSE]
count_mat <- count_mat[, common_samples, drop = FALSE]
clinical  <- clinical[match(common_samples, clinical$gsm), ]

stopifnot(identical(colnames(expr_mat), clinical$gsm))
stopifnot(identical(colnames(count_mat), clinical$gsm))

# ============================================================
# 7. Ensembl -> Symbol convert
# ============================================================
cat("\n========== 7. Ensembl -> Symbol ==========\n")

gene_ids <- rownames(expr_mat)
is_ensembl <- mean(grepl("^ENSG", gene_ids)) > 0.5
cat("是否 Ensembl:", is_ensembl, "\n")

if (is_ensembl) {
  local_anno_file <- "E:/retina_tools/aging_rgc_axon/AMD/results/01_1_prepare_bulk_data/data/gene_annotation.rds"
  if (file.exists(local_anno_file)) {
    anno <- readRDS(local_anno_file) |> as.data.frame()
    col_ens <- intersect(c("ensembl","ENSEMBL","ensembl_id"), names(anno))[1]
    col_sym <- intersect(c("symbol","SYMBOL","gene_symbol"), names(anno))[1]
    anno <- anno |>
      select(ensembl = all_of(col_ens), symbol = all_of(col_sym)) |>
      filter(!is.na(symbol)) |> distinct()
    
    gene_ids_clean <- sub("\\.\\d+$", "", gene_ids)
    sym_map <- anno$symbol[match(gene_ids_clean, anno$ensembl)]
    
    cat("通过本地 annotation 映射:", sum(!is.na(sym_map)), "/", length(sym_map), "\n")
  } else {
    cat("本地 annotation 不存在，无法转换\n")
    sym_map <- rep(NA, length(gene_ids))
  }
  
  # build Symbol expressionmatrix（deduplicate）
  keep <- !is.na(sym_map)
  expr_sym <- expr_mat[keep, , drop = FALSE]
  rownames(expr_sym) <- sym_map[keep]
  expr_sym <- expr_sym[!duplicated(rownames(expr_sym)), , drop = FALSE]
  
  cat("Symbol 矩阵:", nrow(expr_sym), "x", ncol(expr_sym), "\n")
} else {
  expr_sym <- expr_mat
}

# ============================================================
# 8. save
# ============================================================
cat("\n========== 8. 保存 ==========\n")

saveRDS(expr_sym,   file.path(data_dir, "expr_symbol.rds"))
saveRDS(expr_mat,   file.path(data_dir, "expr_ensembl.rds"))
saveRDS(count_mat,  file.path(data_dir, "expr_count.rds"))
saveRDS(clinical,   file.path(data_dir, "clinical.rds"))

cat("保存:\n")
cat("  ", file.path(data_dir, "expr_symbol.rds"), "\n")
cat("  ", file.path(data_dir, "expr_ensembl.rds"), "\n")
cat("  ", file.path(data_dir, "expr_count.rds"), "\n")
cat("  ", file.path(data_dir, "clinical.rds"), "\n")

# ============================================================
# ============================================================
cat("\n========== 9. 摘要 ==========\n")
cat("GSE135092 总样本:", nrow(clinical), "\n")
cat("  Retina:", sum(clinical$tissue_type == "Retina", na.rm = TRUE), "\n")
cat("  RPE:",    sum(clinical$tissue_type == "RPE", na.rm = TRUE), "\n")
cat("  Macula:", sum(clinical$region == "Macula", na.rm = TRUE), "\n")
cat("  non-Macula:", sum(clinical$region == "non-Macula", na.rm = TRUE), "\n\n")

cat("Macular RPE 样本数:\n")
macular_rpe <- clinical |> filter(region == "Macula", tissue_type == "RPE")
print(table(macular_rpe$amd_binary))

cat("\nMacular Retina 样本数:\n")
macular_retina <- clinical |> filter(region == "Macula", tissue_type == "Retina")
print(table(macular_retina$amd_binary))

cat("\n✅ 07_8a GSE135092 prep 完成\n")
cat("Output:", data_dir, "\n")
