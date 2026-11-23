# ============================================================
# 08_1_oxphos_gene_list.r
# Output:  results/08_1_oxphos_gene_list/
# ============================================================

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(data.table); library(dplyr)
})

DIR_THIS <- setup_output_dir("08_1_oxphos_gene_list")
ROOT     <- "E:/retina_tools/aging_rgc_axon/AMD"

log_info("=== 08_1 oxphos gene list ===")

gs_file <- file.path(ROOT, "results/01_2_prepare_gene_sets/data/gene_sets.rds")
if (!file.exists(gs_file)) stop("gene_sets.rds not found: ", gs_file)
gs <- readRDS(gs_file)
log_info("Loaded {length(gs)} gene sets")

oxphos_sets <- c(
  "Mito_Complex_I","Mito_Complex_II","Mito_Complex_III",
  "Mito_Complex_IV","Mito_Complex_V",
  "Mito_CI_Subunits","Mito_CII_Subunits","Mito_CIII_Subunits",
  "Mito_CIV_Subunits","Mito_CV_Subunits",
  "Mito_CI_Assembly","Mito_CIII_Assembly",
  "Mito_CIV_Assembly","Mito_CV_Assembly",
  "Mito_Respirasome"
)
available <- intersect(oxphos_sets, names(gs))
missing   <- setdiff(oxphos_sets, names(gs))
log_info("Available: {length(available)} / {length(oxphos_sets)}")
if (length(missing) > 0) log_warn("Missing: {paste(missing, collapse=', ')}")
if (length(available) == 0) stop("No OXPHOS gene sets found")

gene_dt <- rbindlist(lapply(available, function(nm) {
  data.table(
    geneset = nm,
    complex = fifelse(grepl("Complex_I$|CI_", nm), "CI",
              fifelse(grepl("Complex_II$|CII_", nm), "CII",
              fifelse(grepl("Complex_III$|CIII_", nm), "CIII",
              fifelse(grepl("Complex_IV$|CIV_", nm), "CIV",
              fifelse(grepl("Complex_V$|CV_", nm), "CV",
              fifelse(grepl("Respirasome", nm), "Respirasome", "Other")))))),
    gene = gs[[nm]]$genes
  )
}))

gene_unique <- gene_dt[, .(
  complexes = paste(sort(unique(complex)), collapse = "|"),
  n_sets = uniqueN(geneset),
  sets = paste(sort(unique(geneset)), collapse = "|")
), by = gene]

log_info("Total unique OXPHOS genes: {nrow(gene_unique)}")

id_map_file <- file.path(ROOT,
  "results/01_1_prepare_bulk_data/data/gene_annotation.rds")
if (file.exists(id_map_file)) {
  id_map <- as.data.frame(readRDS(id_map_file))
  col_ens <- intersect(c("ensembl","ENSEMBL","ensembl_id","gene_id"), names(id_map))[1]
  col_sym <- intersect(c("symbol","SYMBOL","gene_symbol","gene_name"), names(id_map))[1]
  if (!is.na(col_ens) && !is.na(col_sym)) {
    id_map <- id_map[, c(col_ens, col_sym)]
    colnames(id_map) <- c("ensembl", "symbol")
    id_map$symbol <- toupper(trimws(id_map$symbol))
    id_map <- id_map[!duplicated(id_map$symbol), ]
    gene_unique <- merge(gene_unique, id_map,
                         by.x = "gene", by.y = "symbol", all.x = TRUE)
    log_info("Mapped to Ensembl: {sum(!is.na(gene_unique$ensembl))} / {nrow(gene_unique)}")
  } else {
    gene_unique$ensembl <- NA_character_
  }
} else {
  gene_unique$ensembl <- NA_character_
  log_warn("gene_annotation.rds not found; ensembl = NA")
}

save_rds(gene_unique, DIR_THIS, "oxphos_target_genes.rds")
save_table(gene_unique, DIR_THIS, "oxphos_target_genes.csv")

summary_dt <- gene_unique[, .(
  n_genes = .N, n_with_ensembl = sum(!is.na(ensembl))
), by = complexes][order(-n_genes)]
save_table(summary_dt, DIR_THIS, "oxphos_gene_summary.csv")

by_complex <- gene_unique[, .(gene, ensembl), by = complexes]
save_table(by_complex, DIR_THIS, "oxphos_gene_by_complex.csv")

cat("\n=== OXPHOS gene list summary ===\n")
cat("Total unique genes:", nrow(gene_unique), "\n")
cat("With Ensembl ID:", sum(!is.na(gene_unique$ensembl)), "\n")
print(summary_dt)
log_info("=== 08_1 done ===")

