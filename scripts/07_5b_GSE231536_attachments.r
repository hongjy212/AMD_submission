# ============================================================
# 07_5b_GSE231536_attachments.r
#   Part B：AMD causalgenescross（GSE135092 mmc5）
# Output：results/07_5b_GSE231536_attachments/
# ============================================================
source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({ library(readxl); library(dplyr) })

ROOT <- "E:/retina_tools/aging_rgc_axon/AMD"
OUT  <- file.path(ROOT, "results", "07_5b_GSE231536_attachments")
TAB  <- file.path(OUT, "tables")
FIG  <- file.path(OUT, "figures")
dir.create(TAB, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)

write_tab <- function(df, name) {
  p <- file.path(TAB, name)
  write.csv(df, p, row.names = FALSE)
  log_info("saved: {name} ({nrow(df)} rows)")
  invisible(p)
}

log_info("=== 07_5b attachments + AMD causal ===")

DIR_MET <- "E:/retina_tools/raw/GSE231536"
MMC5    <- "E:/retina_tools/pubdata/eQTL/GSE135092/mmc5.xlsx"

# ============================================================
# ============================================================
log_info("--- 0. 读目标基因 ---")

sig16 <- toupper(read.csv(file.path(
  ROOT, "results", "02_1_differential_expression_main",
  "tables", "sig16_genes.csv"))$gene)
log_info("sig16: {length(sig16)}")

polyF <- read.csv(file.path(
  ROOT, "results", "04_2_pathway_score",
  "tables", "pathway_poly3F_adjusted.csv"))
p3_paths <- polyF$pathway[polyF$pattern == "p3_only"]
GENE_SETS <- readRDS(file.path(
  ROOT, "results", "01_2_prepare_gene_sets", "data", "gene_sets.rds"))
p3_genes <- unique(toupper(unlist(lapply(
  intersect(p3_paths, names(GENE_SETS)),
  function(p) GENE_SETS[[p]]$genes))))
log_info("p3_only 通路基因: {length(p3_genes)}")

net <- readRDS(file.path(
  ROOT, "results", "05_1_consensus_wgcna", "data", "consensusNet.rds"))
anno <- readRDS(file.path(
  ROOT, "results", "01_1_prepare_bulk_data", "data", "gene_annotation.rds"))
ens2sym <- setNames(toupper(anno$symbol), sub("\\..*$", "", anno$ensembl))

mod_genes <- list(
  greenyellow = names(net$colors)[net$colors == "greenyellow"],
  black       = names(net$colors)[net$colors == "black"],
  green       = names(net$colors)[net$colors == "green"]
)
mod_sym <- lapply(mod_genes, function(ens) {
  s <- ens2sym[sub("\\..*$", "", ens)]
  s[!is.na(s)]
})
for (m in names(mod_sym)) log_info("{m}: {length(mod_sym[[m]])} 基因")

targets <- list(
  sig16           = sig16,
  p3_pathway      = p3_genes,
  mod_greenyellow = mod_sym$greenyellow,
  mod_black       = mod_sym$black,
  mod_green       = mod_sym$green
)

# ============================================================
# ============================================================
log_info("--- Part A: GSE231536 附件 ---")

extract_targets <- function(df) {
  col <- grep("Significant.*target genes|eQTM target genes",
              colnames(df), ignore.case = TRUE)[1]
  if (is.na(col)) {
    log_warn("extract_targets: no matching column, returning empty")
    return(character(0))
  }
  all_str <- paste(df[[col]], collapse = ",")
  genes <- trimws(unlist(strsplit(all_str, "[,;]")))
  genes <- genes[!is.na(genes) & genes != "" & genes != "-"]
  toupper(unique(genes))
}

m12 <- read_excel(file.path(DIR_MET, "41467_2024_46063_MOESM12_ESM.xlsx"),
                  sheet = "eQTM_target_genes_GO")
m12_targets <- extract_targets(m12)

m14 <- read_excel(file.path(DIR_MET, "41467_2024_46063_MOESM14_ESM.xlsx"))
m14_targets <- extract_targets(m14)

eqtm_targets <- unique(c(m12_targets, m14_targets))
log_info("eQTM targets: {length(eqtm_targets)}")

m6 <- read_excel(file.path(DIR_MET, "41467_2024_46063_MOESM6_ESM.xlsx"),
                 sheet = "eGenes")
egenes <- toupper(unique(m6$Gene_Symbol))
log_info("eGenes: {length(egenes)}")

refs_A <- list(
  eQTM_targets = eqtm_targets,
  eGenes       = egenes
)

cross_A <- expand.grid(
  target = names(targets),
  reference = names(refs_A),
  stringsAsFactors = FALSE)
cross_A$n_target <- sapply(cross_A$target,
                            function(t) length(targets[[t]]))
cross_A$n_ref    <- sapply(cross_A$reference,
                            function(r) length(refs_A[[r]]))
cross_A$n_hit    <- mapply(function(t, r)
  length(intersect(targets[[t]], refs_A[[r]])),
  cross_A$target, cross_A$reference)
cross_A$pct_hit  <- round(100 * cross_A$n_hit / cross_A$n_target, 1)
cross_A$pct_bg   <- round(100 * cross_A$n_ref / 20000, 1)
cross_A$fold     <- ifelse(cross_A$pct_bg > 0,
                            round(cross_A$pct_hit / cross_A$pct_bg, 2),
                            NA)

write_tab(cross_A, "cross_A_eqtm_egene.csv")

cat("\n=== Part A：eQTM + eGene 交叉 ===\n")
print(cross_A[, c("target","reference","n_target","n_ref",
                   "n_hit","pct_hit","pct_bg","fold")],
      row.names = FALSE)

detail_A <- list()
for (tgt in names(targets)) {
  for (ref in names(refs_A)) {
    hits <- intersect(targets[[tgt]], refs_A[[ref]])
    if (length(hits) > 0) {
      detail_A[[length(detail_A) + 1]] <- data.frame(
        target = tgt, reference = ref, gene = hits,
        stringsAsFactors = FALSE)
    }
  }
}
if (length(detail_A) > 0) {
  write_tab(do.call(rbind, detail_A), "cross_A_detail.csv")
}

gy_eqtm <- intersect(targets$mod_greenyellow, refs_A$eQTM_targets)
write_tab(data.frame(gene = sort(gy_eqtm)), "greenyellow_eqtm_genes.csv")
log_info("greenyellow ∩ eQTM targets: {length(gy_eqtm)} genes")

p3_names_norm <- toupper(gsub("[^A-Z0-9]+","_", p3_paths))
p3_clean <- sub("^GOBP_|^GOCC_|^KEGG_|^HALLMARK_", "", p3_names_norm)
go_all <- toupper(unique(c(m12$`Gene set`, m14$`Gene set`)))
hit_go <- intersect(p3_clean, go_all)
log_info("GO 直接匹配: {length(hit_go)}")

# ============================================================
# Part B: GSE135092 mmc5 — AMD causalgenes
# ============================================================
log_info("--- Part B: AMD 因果基因 ---")

if (!file.exists(MMC5)) {
  log_warn("mmc5 不存在，跳过 Part B")
} else {
  raw <- read_excel(MMC5, col_names = FALSE, .name_repair = "minimal")
  row1 <- as.character(raw[1, ]); row2 <- as.character(raw[2, ])
  row1[is.na(row1)] <- ""; row2[is.na(row2)] <- ""

  last <- ""
  for (i in seq_along(row1)) {
    if (row1[i] != "") last <- row1[i]
    else row1[i] <- last
  }

  col_names <- mapply(function(r1, r2) {
    if (r1 == "" && r2 == "") return(paste0("V", seq_along(row1)[1]))
    if (r1 == "") return(r2)
    if (r2 == "" || r1 == r2) return(r1)
    paste(r1, r2, sep = " | ")
  }, row1, row2, SIMPLIFY = TRUE)
  col_names <- make.unique(trimws(gsub("\\s+", " ", col_names)))

  dat <- raw[-c(1, 2), , drop = FALSE]
  colnames(dat) <- col_names
  log_info("mmc5: {nrow(dat)} 行 × {ncol(dat)} 列")

  find_col <- function(pat) {
    h <- grep(pat, colnames(dat), ignore.case = TRUE, value = TRUE)
    if (length(h) == 0) NA_character_ else h[1]
  }
  c_gene      <- find_col("Putative.*Causal")
  c_coloc     <- find_col("COLOC.*Probability|Probability$")
  c_eqtl_r_ma <- find_col("Retina, macula.*eQTL symbol")
  c_eqtl_r_nm <- find_col("Retina, non macula.*eQTL symbol")
  c_eqtl_p_ma <- find_col("RPE, macula.*eQTL symbol")
  c_eqtl_p_nm <- find_col("RPE, non-macula.*eQTL symbol")

  if (!is.na(c_coloc)) dat[[c_coloc]] <- suppressWarnings(as.numeric(dat[[c_coloc]]))
  dat$gene_causal <- toupper(trimws(as.character(dat[[c_gene]])))
  dat$coloc_p     <- if (!is.na(c_coloc)) dat[[c_coloc]] else NA_real_

  dat$tier <- NA_character_
  dat$tier[!is.na(dat$coloc_p) & dat$coloc_p >= 0.8] <- "Tier1"
  dat$tier[!is.na(dat$coloc_p) & dat$coloc_p >= 0.5 &
           dat$coloc_p < 0.8] <- "Tier2"

  collect_genes <- function(col) {
    if (is.na(col)) return(character(0))
    v <- toupper(trimws(as.character(dat[[col]])))
    v[!is.na(v) & v != "" & v != "NA"]
  }
  gene_causal    <- unique(na.omit(dat$gene_causal))
  all_amd_causal <- unique(c(gene_causal,
                              collect_genes(c_eqtl_r_ma),
                              collect_genes(c_eqtl_r_nm),
                              collect_genes(c_eqtl_p_ma),
                              collect_genes(c_eqtl_p_nm)))
  tier1_genes <- unique(na.omit(dat$gene_causal[dat$tier == "Tier1"]))
  tier2_genes <- unique(na.omit(dat$gene_causal[dat$tier == "Tier2"]))

  log_info("mmc5 位点: {nrow(dat)}, Tier1: {sum(dat$tier=='Tier1', na.rm=TRUE)}, Tier2: {sum(dat$tier=='Tier2', na.rm=TRUE)}")
  log_info("causal_gene: {length(gene_causal)}")
  log_info("all_amd_causal: {length(all_amd_causal)}")

  refs_B <- list(
    causal_gene   = gene_causal,
    causal_any    = all_amd_causal,
    tier1         = tier1_genes,
    tier2         = tier2_genes,
    retina_macula = collect_genes(c_eqtl_r_ma),
    retina_nonmac = collect_genes(c_eqtl_r_nm),
    rpe_macula    = collect_genes(c_eqtl_p_ma),
    rpe_nonmac    = collect_genes(c_eqtl_p_nm)
  )

  cross_B <- expand.grid(
    target = names(targets),
    reference = names(refs_B),
    stringsAsFactors = FALSE)
  cross_B$n_target <- sapply(cross_B$target,
                              function(t) length(targets[[t]]))
  cross_B$n_ref    <- sapply(cross_B$reference,
                              function(r) length(refs_B[[r]]))
  cross_B$n_hit    <- mapply(function(t, r)
    length(intersect(targets[[t]], refs_B[[r]])),
    cross_B$target, cross_B$reference)
  cross_B$pct_hit  <- round(100 * cross_B$n_hit / cross_B$n_target, 1)
  cross_B$pct_bg   <- round(100 * cross_B$n_ref / 20000, 2)
  cross_B$fold     <- ifelse(cross_B$pct_bg > 0,
                              round(cross_B$pct_hit / cross_B$pct_bg, 2),
                              NA)

  write_tab(cross_B, "cross_B_AMD_causal.csv")

  cat("\n=== Part B：AMD 因果基因交叉 ===\n")
  sig16_B <- cross_B[cross_B$target == "sig16", ]
  print(sig16_B[, c("reference","n_ref","n_hit","pct_hit","pct_bg","fold")],
        row.names = FALSE)

  detail_B <- list()
  for (tgt in names(targets)) {
    for (ref in names(refs_B)) {
      hits <- intersect(targets[[tgt]], refs_B[[ref]])
      if (length(hits) > 0) {
        detail_B[[length(detail_B) + 1]] <- data.frame(
          target = tgt, reference = ref, gene = hits,
          stringsAsFactors = FALSE)
      }
    }
  }
  if (length(detail_B) > 0) {
    write_tab(do.call(rbind, detail_B), "cross_B_detail.csv")
  }

  safe_cols <- c(c_gene, c_coloc, c_eqtl_r_ma, c_eqtl_p_ma, "tier")
  safe_cols <- safe_cols[!is.na(safe_cols)]
  write_tab(dat[, safe_cols, drop = FALSE], "mmc5_AMD_causal_full.csv")
}

# ============================================================
# summary
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("=== 07_5b 汇总 ===\n")
cat(strrep("=", 70), "\n", sep = "")

cat("\n【Part A】GSE231536 附件\n")
print(cross_A[, c("target","reference","n_target","n_ref",
                   "n_hit","pct_hit","fold")],
      row.names = FALSE)
cat("\nGO 直接匹配:", length(hit_go), "\n")
cat("\ngreenyellow ∩ eQTM targets:", length(gy_eqtm), "genes\n")

if (exists("cross_B")) {
  cat("\n【Part B】AMD 因果基因\n")
  print(cross_B[, c("target","reference","n_target","n_ref",
                     "n_hit","pct_hit","fold")],
        row.names = FALSE)
}

cat("\n关键结论：\n")
if (exists("cross_B")) {
  sig16_hit_A <- cross_A[cross_A$target == "sig16", ]
  sig16_hit_B <- cross_B[cross_B$target == "sig16", ]
  cat("  sig16 ∩ eQTM targets =",
      sig16_hit_A$n_hit[sig16_hit_A$reference == "eQTM_targets"], "\n")
  cat("  sig16 ∩ eGenes =",
      sig16_hit_A$n_hit[sig16_hit_A$reference == "eGenes"], "/", length(sig16), "\n")
  cat("  sig16 ∩ AMD causal =",
      max(sig16_hit_B$n_hit, na.rm = TRUE), "/", length(sig16), "\n")
}

cat("\nOUT =", OUT, "\n")
log_info("=== 07_5b done ===")

