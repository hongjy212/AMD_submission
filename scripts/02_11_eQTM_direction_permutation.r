# ============================================================
# 02_11_eQTM_direction_permutation.r (v3)
# Input：
#   pathway_MGS1_vs_MGS3_eQTM.csv           (488 CpG，p3_pathway)
#   gene16_MGS1_vs_MGS3_eQTM.csv            (16 CpG，sig16)
#   02_10_OXPHOS_variance_decomposition/.../threeway_consistency_greenyellow.csv
#                                            (8 CpG，greenyellow ∩ eQTM)
# Output：results/02_11_eQTM_direction_permutation/tables/
# ============================================================

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tidyr)
})

ROOT <- DIR_ROOT
OUT  <- file.path(DIR_RESULT, "02_11_eQTM_direction_permutation")
TAB  <- file.path(OUT, "tables")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
set.seed(SEED)
log_info("=== 02_11 eQTM direction permutation (v3) ===")

# ------------------------------------------------------------
# ------------------------------------------------------------
run_block_permutation <- function(df, label, n_perm = 10000) {
  df <- df %>% filter(!is.na(consistent))
  if (nrow(df) < 3) {
    log_warn(sprintf("[%s] too few CpGs (%d), skip", label, nrow(df)))
    return(NULL)
  }
  if (!"gene_sym" %in% colnames(df)) {
    log_warn(sprintf("[%s] no gene_sym", label))
    return(NULL)
  }

  n_cpg <- nrow(df)
  n_consistent <- sum(df$consistent)
  obs_pct <- n_consistent / n_cpg
  binom_p <- binom.test(n_consistent, n_cpg, p = 0.5)$p.value

  genes <- unique(df$gene_sym)
  n_genes <- length(genes)

  log_info(sprintf("[%s] %d CpGs from %d genes, observed %d/%d (%.1f%%)",
                   label, n_cpg, n_genes, n_consistent, n_cpg, 100 * obs_pct))

  gene_idx_list <- split(seq_len(n_cpg), df$gene_sym)
  consistent_vec <- df$consistent

  perm_pct <- numeric(n_perm)
  for (i in seq_len(n_perm)) {
    perm_genes <- sample(genes, n_genes, replace = FALSE)
    perm_consistent <- logical(n_cpg)
    pos <- 1
    for (g in perm_genes) {
      idx <- gene_idx_list[[g]]
      n_g <- length(idx)
      perm_consistent[pos:(pos + n_g - 1)] <- consistent_vec[idx]
      pos <- pos + n_g
    }
    perm_pct[i] <- mean(perm_consistent)
  }

  perm_p <- mean(perm_pct >= obs_pct)

  data.frame(
    dataset        = label,
    n_cpg          = n_cpg,
    n_genes        = n_genes,
    n_consistent   = n_consistent,
    pct_consistent = round(100 * obs_pct, 1),
    binom_p        = signif(binom_p, 3),
    perm_p_block   = signif(perm_p, 4),
    perm_mean_pct  = round(100 * mean(perm_pct), 1),
    perm_q025_pct  = round(100 * quantile(perm_pct, 0.025), 1),
    perm_q975_pct  = round(100 * quantile(perm_pct, 0.975), 1),
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------
# ------------------------------------------------------------
log_info("--- pathway eQTM (488 CpG) ---")
pathway_file <- file.path(DIR_RESULT,
  "07_5_GSE231536_methylation/tables/pathway_MGS1_vs_MGS3_eQTM.csv")
pathway <- read_csv(pathway_file, show_col_types = FALSE) %>%
  mutate(consistent = as.logical(consistent))
res_pathway <- run_block_permutation(pathway, "p3_pathway_488CpG")
if (!is.null(res_pathway)) {
  cat("\n=== pathway 层区块置换 ===\n")
  print(res_pathway, row.names = FALSE)
  write_csv(res_pathway, file.path(TAB, "block_permutation_pathway.csv"))
}

# ------------------------------------------------------------
# ------------------------------------------------------------
log_info("--- sig16 eQTM (16 CpG) ---")
sig16_file <- file.path(DIR_RESULT,
  "07_5_GSE231536_methylation/tables/gene16_MGS1_vs_MGS3_eQTM.csv")
if (file.exists(sig16_file)) {
  sig16 <- read_csv(sig16_file, show_col_types = FALSE) %>%
    mutate(consistent = as.logical(consistent))
  res_sig16 <- run_block_permutation(sig16, "sig16_16CpG")
  if (!is.null(res_sig16)) {
    cat("\n=== sig16 层区块置换 ===\n")
    print(res_sig16, row.names = FALSE)
    write_csv(res_sig16, file.path(TAB, "block_permutation_sig16.csv"))
  }
} else {
  log_warn("sig16 eQTM file not found")
  res_sig16 <- NULL
}

# ------------------------------------------------------------
# ------------------------------------------------------------
log_info("--- greenyellow ∩ eQTM (8 CpG) ---")
gy_file <- file.path(DIR_RESULT,
  "02_10_OXPHOS_variance_decomposition/tables/threeway_consistency_greenyellow.csv")
if (file.exists(gy_file)) {
  gy <- read_csv(gy_file, show_col_types = FALSE)
  if (!"gene_sym" %in% colnames(gy) && "symbol" %in% colnames(gy)) {
    gy <- gy %>% rename(gene_sym = symbol)
  }
  gy <- gy %>%
    mutate(consistent = as.logical(consistent)) %>%
    filter(!is.na(consistent))
  res_gy <- run_block_permutation(gy, "greenyellow_8CpG")
  if (!is.null(res_gy)) {
    cat("\n=== greenyellow ∩ eQTM 区块置换 ===\n")
    print(res_gy, row.names = FALSE)
    write_csv(res_gy, file.path(TAB, "block_permutation_greenyellow.csv"))
  }
} else {
  log_warn("greenyellow threeway file not found")
  res_gy <- NULL
}

# ------------------------------------------------------------
# ------------------------------------------------------------
log_info("--- pathway independent CpG (1 per gene) ---")
if ("p_FDR" %in% colnames(pathway)) {
  pathway_rep <- pathway %>%
    group_by(gene_sym) %>%
    slice_min(p_FDR, n = 1, with_ties = FALSE) %>%
    ungroup()
} else {
  pathway_rep <- pathway %>%
    group_by(gene_sym) %>%
    slice(1) %>%
    ungroup()
}
log_info(sprintf("pathway independent CpG: %d genes", nrow(pathway_rep)))

res_pathway_rep <- run_block_permutation(pathway_rep, "p3_pathway_independent")
if (!is.null(res_pathway_rep)) {
  cat("\n=== pathway 独立 CpG 区块置换 ===\n")
  print(res_pathway_rep, row.names = FALSE)
  write_csv(res_pathway_rep, file.path(TAB, "block_permutation_pathway_independent.csv"))
}

# ------------------------------------------------------------
# 5. summary
# ------------------------------------------------------------
cat("\n", strrep("=", 60), "\n", sep = "")
cat("=== 02_11 汇总 ===\n")
cat(strrep("=", 60), "\n", sep = "")

all_res <- bind_rows(
  res_pathway, res_sig16, res_gy, res_pathway_rep
)

if (nrow(all_res) > 0) {
  cat("\n所有置换检验结果：\n")
  print(all_res[, c("dataset","n_cpg","n_genes","pct_consistent",
                    "binom_p","perm_p_block",
                    "perm_mean_pct","perm_q025_pct","perm_q975_pct")],
        row.names = FALSE)
  write_csv(all_res, file.path(TAB, "block_permutation_summary.csv"))
}

log_info("=== 02_11 done ===")

