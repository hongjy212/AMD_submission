# ============================================================
# 07_6_prepare_sem_bnlearn_data.r (v2: 419)
# Purpose: shared data preparation for both
#   - 07_6a_bayesian_network.r  (bnlearn, uses df_disc)
#   - 07_6b_lavaan_mediation.r  (lavaan / SEM, uses df_cont)
#   - 07_6c_bn_sensitivity.r    (optional, uses df_disc)
#
#   - clinical.rds → clinical_wgcna.rds
#   - expr_corrected.rds → expr_wgcna.rds
#
# Input:
#   results/05_1_consensus_wgcna/data/consensusNet.rds
#   results/01_1b_batch_enhance/data/{clinical_wgcna,expr_wgcna}.rds   ★ 419
#   results/01_1_prepare_bulk_data/data/gene_annotation.rds
#   results/07_5c_OXPHOS_analysis/tables/OXPHOS13_gene_vs_grade.csv
#   results/01_2_prepare_gene_sets/data/gene_sets.rds
# Output:
#   results/07_6_prepare_sem_bnlearn_data/data/df_cont.rds
#   results/07_6_prepare_sem_bnlearn_data/data/df_disc.rds
#   results/07_6_prepare_sem_bnlearn_data/tables/
# ============================================================
source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(dplyr); library(WGCNA)
})

ROOT <- "E:/retina_tools/aging_rgc_axon/AMD"
OUT  <- file.path(ROOT, "results", "07_6_prepare_sem_bnlearn_data")
TAB  <- file.path(OUT, "tables")
DAT  <- file.path(OUT, "data")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(DAT, recursive = TRUE, showWarnings = FALSE)

write_tab <- function(df, name) {
  p <- file.path(TAB, name)
  write.csv(df, p, row.names = FALSE)
  log_info("saved: {name} ({nrow(df)} rows)")
  invisible(p)
}

log_info("=== 07_6_prepare_sem_bnlearn_data v2 (419) ===")

# ============================================================
# 0. Load  ★ v2：419
# ============================================================
log_info("--- 0. Load ---")

clin <- readRDS(file.path(ROOT, "results", "01_1b_batch_enhance",
                           "data", "clinical_wgcna.rds"))
expr <- readRDS(file.path(ROOT, "results", "01_1b_batch_enhance",
                           "data", "expr_wgcna.rds"))
anno <- readRDS(file.path(ROOT, "results", "01_1_prepare_bulk_data",
                           "data", "gene_annotation.rds"))
net  <- readRDS(file.path(ROOT, "results", "05_1_consensus_wgcna",
                           "data", "consensusNet.rds"))
GENE_SETS <- readRDS(file.path(ROOT, "results", "01_2_prepare_gene_sets",
                                "data", "gene_sets.rds"))

log_info("clinical: {nrow(clin)} x {ncol(clin)}")
log_info("expr:     {nrow(expr)} x {ncol(expr)}")
log_info("net$colors length: {length(net$colors)}")

# ============================================================
# 1. Module eigengenes
# ============================================================
log_info("--- 1. Module eigengenes ---")

if (!identical(names(net$colors), rownames(expr))) {
  log_warn("net$colors names != expr rownames; reordering...")
  if (all(rownames(expr) %in% names(net$colors))) {
    net$colors <- net$colors[rownames(expr)]
  } else {
    stop("Cannot match net$colors to expr rownames")
  }
}

MEs <- moduleEigengenes(t(expr), colors = net$colors)$eigengenes
log_info("MEs: {nrow(MEs)} x {ncol(MEs)}")

ME_candidates <- c("MEgreenyellow", "MEblack", "MEgreen")
ME_use <- intersect(ME_candidates, colnames(MEs))
log_info("ME used: {paste(ME_use, collapse=', ')}")

saveRDS(MEs, file.path(ROOT, "results", "05_1_consensus_wgcna",
                        "data", "MEs_global.rds"))
log_info("Also saved: MEs_global.rds")

# ============================================================
# 2. OXPHOS scores
# ============================================================
log_info("--- 2. OXPHOS scores ---")

ens2sym <- setNames(toupper(anno$symbol), sub("\\..*$", "", anno$ensembl))
expr_sym <- expr
rownames(expr_sym) <- ens2sym[rownames(expr)]

oxphos13_file <- file.path(ROOT, "results", "07_5c_OXPHOS_analysis",
                            "tables", "OXPHOS13_gene_vs_grade.csv")
if (!file.exists(oxphos13_file)) stop("07_5c not run")
oxphos13 <- toupper(read.csv(oxphos13_file)$gene)
ox13_in <- intersect(oxphos13, rownames(expr_sym))
log_info("OXPHOS13: {length(ox13_in)} / {length(oxphos13)} genes in expr")

oxphos_path <- toupper(GENE_SETS$OXPHOS$genes)
oxphos_in <- intersect(oxphos_path, rownames(expr_sym))
log_info("OXPHOS pathway: {length(oxphos_in)} / {length(oxphos_path)} genes in expr")

common <- Reduce(intersect, list(rownames(MEs), clin$sample_id, colnames(expr)))
log_info("Common samples: {length(common)}")

OXPHOS13 <- colMeans(t(scale(t(expr_sym[ox13_in, common, drop = FALSE]))),
                      na.rm = TRUE)
OXPHOS_pathway <- colMeans(t(scale(t(expr_sym[oxphos_in, common, drop = FALSE]))),
                            na.rm = TRUE)

# ============================================================
# 2b. Complex III scores
# ============================================================
log_info("--- 2b. Complex III scores ---")

cIII_genes <- c("UQCRQ", "UQCRFS1", "UQCRH")
cIII_in <- intersect(cIII_genes, rownames(expr_sym))
log_info("Complex III genes in expr: {paste(cIII_in, collapse=', ')} ({length(cIII_in)}/{length(cIII_genes)})")

if (length(cIII_in) < 3) {
  log_warn("Not all Complex III genes found; ComplexIII_mean will use available genes.")
}

expr_cIII_z <- t(scale(t(expr_sym[cIII_in, common, drop = FALSE])))
ComplexIII_mean <- colMeans(expr_cIII_z, na.rm = TRUE)

if ("UQCRQ" %in% cIII_in) {
  UQCRQ_z <- as.numeric(expr_cIII_z["UQCRQ", ])
} else {
  UQCRQ_z <- rep(NA_real_, length(common))
  log_warn("UQCRQ not found; UQCRQ_z set to NA.")
}

if (all(c("UQCRFS1","UQCRH") %in% cIII_in)) {
  UQCRFS1_gene_z <- as.numeric(expr_cIII_z["UQCRFS1", ])
  UQCRH_gene_z   <- as.numeric(expr_cIII_z["UQCRH",   ])
} else {
  UQCRFS1_gene_z <- UQCRH_gene_z <- rep(NA_real_, length(common))
  log_warn("UQCRFS1 or UQCRH not found; latent SEM Model B will be skipped.")
}

# ============================================================
# 3. Build df_cont (for SEM)
# ============================================================
log_info("--- 3. Build df_cont ---")

clin2 <- clin[match(common, clin$sample_id), ]

df_cont <- data.frame(
  sample_id         = common,
  MEgreenyellow     = as.numeric(MEs[common, "MEgreenyellow"]),
  MEblack           = as.numeric(MEs[common, "MEblack"]),
  MEgreen           = as.numeric(MEs[common, "MEgreen"]),
  OXPHOS13          = as.numeric(OXPHOS13),
  OXPHOS_pathway    = as.numeric(OXPHOS_pathway),
  ComplexIII_mean   = as.numeric(ComplexIII_mean),
  UQCRQ_z           = as.numeric(UQCRQ_z),
  UQCRFS1_gene_z    = as.numeric(UQCRFS1_gene_z),
  UQCRH_gene_z      = as.numeric(UQCRH_gene_z),
  age               = as.numeric(clin2$age_n),
  sex               = factor(clin2$sex, levels = c("Female", "Male")),
  RIN               = as.numeric(clin2$rin_n),
  PMI               = as.numeric(clin2$pmi_n),
  grade             = as.integer(clin2$grade),
  stringsAsFactors  = FALSE
)

df_cont <- df_cont[complete.cases(df_cont), ]
log_info("df_cont: {nrow(df_cont)} x {ncol(df_cont)}")

saveRDS(df_cont, file.path(DAT, "df_cont.rds"))
log_info("saved: df_cont.rds")

# ============================================================
# 4. Build df_disc (for bnlearn)
# ============================================================
log_info("--- 4. Build df_disc ---")

bn_exclude <- c("UQCRFS1_gene_z", "UQCRH_gene_z")
bn_nodes <- setdiff(colnames(df_cont), c("sample_id", bn_exclude))
df_disc <- df_cont[, bn_nodes, drop = FALSE]

for (v in setdiff(bn_nodes, c("sex", "grade"))) {
  qs <- unique(quantile(df_disc[[v]], probs = c(0, 1/3, 2/3, 1), na.rm = TRUE))
  if (length(qs) < 4) {
    qs <- unique(quantile(df_disc[[v]], probs = c(0, 0.5, 1), na.rm = TRUE))
    df_disc[[v]] <- cut(df_disc[[v]], breaks = qs,
                         labels = FALSE, include.lowest = TRUE)
  } else {
    df_disc[[v]] <- cut(df_disc[[v]], breaks = qs,
                         labels = c("low", "mid", "high"), include.lowest = TRUE)
  }
}

df_disc$grade <- factor(df_disc$grade, levels = 1:4,
                         labels = c("G1", "G2", "G3", "G4"))
df_disc <- df_disc[complete.cases(df_disc), ]
log_info("df_disc: {nrow(df_disc)} x {ncol(df_disc)}")

saveRDS(df_disc, file.path(DAT, "df_disc.rds"))
log_info("saved: df_disc.rds")

# ============================================================
# 5. Summaries
# ============================================================
log_info("--- 5. Summaries ---")

num_vars <- c("MEgreenyellow", "MEblack", "MEgreen",
               "OXPHOS13", "OXPHOS_pathway",
               "ComplexIII_mean", "UQCRQ_z",
               "age", "RIN", "PMI")

sum_cont <- data.frame(
  variable = num_vars,
  n        = nrow(df_cont),
  mean     = sapply(df_cont[, num_vars], mean, na.rm = TRUE),
  sd       = sapply(df_cont[, num_vars], sd, na.rm = TRUE),
  min      = sapply(df_cont[, num_vars], min, na.rm = TRUE),
  max      = sapply(df_cont[, num_vars], max, na.rm = TRUE)
)
write_tab(sum_cont, "continuous_summary.csv")

grade_tab <- as.data.frame(table(df_cont$grade))
colnames(grade_tab) <- c("grade", "n")
write_tab(grade_tab, "grade_distribution.csv")

sex_tab <- as.data.frame(table(df_cont$sex))
colnames(sex_tab) <- c("sex", "n")
write_tab(sex_tab, "sex_distribution.csv")

var_dict <- data.frame(
  variable = c("sample_id",
                "MEgreenyellow", "MEblack", "MEgreen",
                "OXPHOS13", "OXPHOS_pathway",
                "ComplexIII_mean", "UQCRQ_z",
                "UQCRFS1_gene_z", "UQCRH_gene_z",
                "age", "sex", "RIN", "PMI", "grade"),
  type = c("ID",
            "continuous", "continuous", "continuous",
            "continuous", "continuous",
            "continuous", "continuous",
            "continuous", "continuous",
            "continuous", "factor", "continuous", "continuous", "ordinal"),
  role = c("ID",
            "module", "module", "module",
            "score", "score",
            "score", "score",
            "indicator", "indicator",
            "covariate", "covariate", "covariate", "covariate", "outcome"),
  description = c(
    "Sample ID",
    "Greenyellow module eigengene (OXPHOS)",
    "Black module eigengene",
    "Green module eigengene",
    "Mean z-score of 13 OXPHOS complex genes",
    "Mean z-score of OXPHOS pathway (187 genes)",
    "Mean z-score of Complex III subunits (UQCRQ, UQCRFS1, UQCRH)",
    "UQCRQ single-gene z-score",
    "UQCRFS1 gene z-score (latent SEM indicator)",
    "UQCRH gene z-score (latent SEM indicator)",
    "Age at death",
    "Sex (Female/Male)",
    "RNA integrity number",
    "Post-mortem interval (hours)",
    "AMD grade (1-4, MGS classification)"
  ),
  stringsAsFactors = FALSE
)
write_tab(var_dict, "variable_dictionary.csv")

# ============================================================
# 6. Summary output
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("=== 07_6_prepare_sem_bnlearn_data v2 summary (419) ===\n")
cat(strrep("=", 70), "\n", sep = "")

cat("\nSamples (df_cont):", nrow(df_cont), "\n")
cat("Samples (df_disc):", nrow(df_disc), "\n")
cat("Variables:", ncol(df_cont) - 1, "\n\n")

cat("Grade distribution:\n")
print(grade_tab, row.names = FALSE)
cat("\nSex distribution:\n")
print(sex_tab, row.names = FALSE)

cat("\nContinuous summary:\n")
print(sum_cont, row.names = FALSE, digits = 3)

cat("\nOutput dir:", OUT, "\n")
cat("  data/df_cont.rds    (continuous, for SEM/lavaan)\n")
cat("  data/df_disc.rds    (discretized, for BN/bnlearn)\n")
cat("\nNext:\n")
cat("  source 07_6a_bayesian_network.r    (reads df_disc)\n")
cat("  source 07_6b_lavaan_mediation.r    (reads df_cont)\n")
log_info("=== 07_6_prepare_sem_bnlearn_data v2 (419) done ===")

