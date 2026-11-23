#!/usr/bin/env Rscript
# ============================================================
# 05_2_module_trait_and_validation.r
#
#   - modulesummary
#
# validation（L0-L6）：
#
# Output：results/05_2_module_trait/
# ============================================================

if (exists(".AMD_SETUP_LOADED", envir = .GlobalEnv)) {
  rm(".AMD_SETUP_LOADED", envir = .GlobalEnv)
}
source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")

suppressPackageStartupMessages({
  library(WGCNA); library(limma); library(dplyr); library(tidyr)
  library(tibble); library(ggplot2); library(pheatmap); library(ggrepel)
})

allowWGCNAThreads()

DIR_THIS <- setup_output_dir("05_2_module_trait")
log_info("=== 05_2_module_trait_and_validation.r start ===")

PADJ_THRESHOLD <- 0.05

# ============================================================
# ============================================================
strip_version <- function(x) sub("\\..*$", "", x)

GENES_16 <- toupper(c("FCF1","UGGT1","CCDC3","S100B","BPNT1","BRI3BP","BTBD19",
                      "DPH7","AKR1B1","ZNF623","ZNF587","RSPH3","TAF3","NDUFB8",
                      "ZNF417","C6orf120"))

run_enrichment <- function(module_genes_symbol, gene_sets, background_symbols) {
  results <- data.frame()
  for (gs_name in names(gene_sets)) {
    gs <- gene_sets[[gs_name]]
    gs_genes <- toupper(gs$genes)
    gs_in_bg <- intersect(gs_genes, background_symbols)
    if (length(gs_in_bg) < 5) next

    hit <- intersect(module_genes_symbol, gs_in_bg)
    n_hit <- length(hit)
    if (n_hit < 3) next

    p <- phyper(n_hit - 1,
                length(gs_in_bg),
                length(background_symbols) - length(gs_in_bg),
                length(module_genes_symbol),
                lower.tail = FALSE)
    expected <- length(module_genes_symbol) * length(gs_in_bg) / length(background_symbols)

    results <- rbind(results, data.frame(
      Geneset = gs_name,
      Group = if (!is.null(gs$group)) gs$group else NA_character_,
      n_geneset = as.integer(length(gs_in_bg)),
      n_mod = as.integer(length(module_genes_symbol)),
      n_hit = as.integer(n_hit),
      expected = round(expected, 2),
      fold = round(n_hit / expected, 2),
      p = p,
      stringsAsFactors = FALSE))
  }
  if (nrow(results) > 0) {
    results$FDR <- p.adjust(results$p, method = "BH")
    results <- results[order(results$p), ]
  }
  results
}

# ============================================================
# ============================================================
log_info("--- 1. 载入 ---")

expr <- read_upstream("01_1b_batch_enhance", "data", "expr_wgcna.rds")
clin <- read_upstream("01_1b_batch_enhance", "data", "clinical_wgcna.rds")
consensusNet <- read_upstream("05_1_consensus_wgcna", "data", "consensusNet.rds")
module_df <- read_upstream("05_1_consensus_wgcna", "tables",
                            "consensus_module_genes_annotated.csv")
aging_scores <- read_upstream("04_1_aging_score", "data", "aging_scores.rds")
score_pathway <- read_upstream("04_2_pathway_score", "data",
                                "pathway_scores_ssgsea.rds")
GENE_SETS <- read_upstream("01_2_prepare_gene_sets", "data", "gene_sets.rds")
anno <- read_upstream("01_1_prepare_bulk_data", "data", "gene_annotation.rds")

common <- intersect(colnames(expr), rownames(clin))
expr <- expr[, common, drop = FALSE]
clin <- clin[common, , drop = FALSE]

log_info("expr: {nrow(expr)} × {ncol(expr)}")
log_info("模块数: {length(unique(module_df$Module))}")

sym_map <- setNames(toupper(anno$symbol), strip_version(anno$ensembl))
expr_row_sym <- sym_map[strip_version(rownames(expr))]

# ============================================================
# ============================================================
log_info("--- 2. 全局 MEs ---")

all_expr_t <- t(expr)
stopifnot(identical(names(consensusNet$colors), colnames(all_expr_t)))

MEs_global <- moduleEigengenes(
  all_expr_t,
  colors = consensusNet$colors,
  verbose = 0
)$eigengenes
MEs_global <- MEs_global[, !grepl("^MEgrey", colnames(MEs_global))]
log_info("全局 MEs: {nrow(MEs_global)} × {ncol(MEs_global)}")
save_rds(MEs_global, DIR_THIS, "MEs_global.rds")

all_modules <- setdiff(unique(module_df$Module), "grey")
log_info("真实模块: {length(all_modules)} 个")

# ============================================================
# ============================================================
log_info("--- 3. 构建性状表 ---")

clin$grade_int <- as.integer(clin$grade)
clin$age_v <- as.numeric(clin$age_n)
clin$rin_v <- as.numeric(clin$rin_n)
clin$pmi_v <- as.numeric(clin$pmi_n)
clin$batch_f <- factor(clin$batch)

sid <- as.character(clin$sample_id)
clin$layer_clean <- sub("^.*_(L[0-9]+)$", "\\1", sid)
clin$layer_clean <- ifelse(grepl("^L[0-9]+$", clin$layer_clean),
                            clin$layer_clean, NA_character_)
clin$layer_f <- factor(clin$layer_clean, levels = paste0("L", 1:8))

log_info("layer NA: {sum(is.na(clin$layer_f))}")

common_aging <- intersect(rownames(MEs_global), colnames(aging_scores))
hUSI_v     <- aging_scores["hUSI", common_aging]
ML_Aging_v <- aging_scores["ML_Aging", common_aging]

genes_16_in <- rownames(expr)[toupper(expr_row_sym) %in% GENES_16]
log_info("16 基因命中: {length(genes_16_in)} / {length(GENES_16)}")
if (length(genes_16_in) > 0) {
  log_info("命中基因: {paste(toupper(expr_row_sym[genes_16_in]), collapse=', ')}")
}
missing_16 <- setdiff(GENES_16, toupper(expr_row_sym[genes_16_in]))
if (length(missing_16) > 0) {
  log_warn("未命中基因: {paste(missing_16, collapse=', ')}")
}

if (length(genes_16_in) >= 5) {
  expr_z_16 <- t(scale(t(expr[genes_16_in, , drop = FALSE])))
  expr_z_16[is.na(expr_z_16)] <- 0
  score_16 <- colMeans(expr_z_16)
} else {
  score_16 <- setNames(rep(NA_real_, ncol(expr)), colnames(expr))
}

mito_pw <- c("Mito_CI_Subunits","Mito_CIII_Subunits","Mito_CIV_Subunits",
              "Mito_CV_Subunits","Mito_Complex_V","Mito_TCA","OXPHOS")
mito_pw <- intersect(mito_pw, rownames(score_pathway))
if (length(mito_pw) > 0) {
  mito_score <- colMeans(score_pathway[mito_pw, , drop = FALSE])
} else {
  mito_score <- setNames(rep(NA_real_, ncol(score_pathway)), colnames(score_pathway))
}

trait_df <- data.frame(
  Sample   = rownames(MEs_global),
  grade    = clin[rownames(MEs_global), "grade_int"],
  age      = clin[rownames(MEs_global), "age_v"],
  rin      = clin[rownames(MEs_global), "rin_v"],
  pmi      = clin[rownames(MEs_global), "pmi_v"],
  batch    = as.integer(clin[rownames(MEs_global), "batch_f"]),
  layer    = as.integer(clin[rownames(MEs_global), "layer_f"]),
  sex_M    = as.integer(clin[rownames(MEs_global), "sex"] == "Male"),
  hUSI     = hUSI_v[rownames(MEs_global)],
  ML_Aging = ML_Aging_v[rownames(MEs_global)],
  genes16  = score_16[rownames(MEs_global)],
  mito     = mito_score[rownames(MEs_global)],
  stringsAsFactors = FALSE
)
log_info("性状表: {nrow(trait_df)} × {ncol(trait_df)}")

# ============================================================
# ============================================================
log_info("--- 4. 模块-性状相关 ---")

trait_cols <- c("grade","age","rin","pmi","batch","layer","sex_M",
                 "hUSI","ML_Aging","genes16","mito")

cor_long <- data.frame()
for (mod_col in colnames(MEs_global)) {
  me_vec <- MEs_global[, mod_col]
  mod_name <- gsub("^ME", "", mod_col)
  for (tr in trait_cols) {
    tr_vec <- trait_df[[tr]]
    if (sum(is.finite(me_vec) & is.finite(tr_vec)) < 10) next
    ct <- tryCatch(cor.test(me_vec, tr_vec, use = "complete.obs"),
                    error = function(e) NULL)
    if (is.null(ct)) next
    cor_long <- rbind(cor_long, data.frame(
      Module = mod_name, Trait = tr,
      r = round(ct$estimate, 3), p = ct$p.value,
      stringsAsFactors = FALSE))
  }
}
cor_long$FDR <- p.adjust(cor_long$p, method = "BH")
cor_long$sig <- cor_long$FDR < PADJ_THRESHOLD
save_table(cor_long, DIR_THIS, "module_trait_correlation.csv")

cat("\n=== 模块-性状相关（FDR<0.05）===\n")
sig_cor <- cor_long[cor_long$FDR < PADJ_THRESHOLD, ]
sig_cor <- sig_cor[order(sig_cor$FDR), ]
print(head(sig_cor, 20), row.names = FALSE, digits = 4)

# ============================================================
# 5. module vs grade ANOVA
# ============================================================
log_info("--- 5. 模块 vs grade ANOVA ---")

anova_df <- data.frame()
for (mod_col in colnames(MEs_global)) {
  me_vec <- MEs_global[, mod_col]
  mod_name <- gsub("^ME", "", mod_col)
  grade_f <- factor(trait_df$grade)
  fit <- tryCatch(aov(me_vec ~ grade_f), error = function(e) NULL)
  if (is.null(fit)) next
  s <- summary(fit)[[1]]
  anova_df <- rbind(anova_df, data.frame(
    Module = mod_name,
    F_stat = round(s[1, "F value"], 2),
    p = s[1, "Pr(>F)"],
    R2 = round(s[1, "Sum Sq"] / sum(s[, "Sum Sq"]), 3),
    stringsAsFactors = FALSE))
}
anova_df$FDR <- p.adjust(anova_df$p, method = "BH")
anova_df <- anova_df[order(anova_df$p), ]
save_table(anova_df, DIR_THIS, "module_vs_grade_anova.csv")

# ============================================================
# ============================================================
log_info("--- 6. 偏相关 ---")

trait_df$batch_f <- factor(trait_df$batch)
trait_df$layer_f <- factor(trait_df$layer)
trait_df$sex_f   <- factor(trait_df$sex_M)

tech_terms <- c("rin + pmi + age")
if (nlevels(trait_df$batch_f) > 1) tech_terms <- paste(tech_terms, "+ batch_f")
if (nlevels(trait_df$layer_f) > 1) tech_terms <- paste(tech_terms, "+ layer_f")
if (nlevels(trait_df$sex_f) > 1)   tech_terms <- paste(tech_terms, "+ sex_f")

design_tech <- model.matrix(
  as.formula(paste("~", tech_terms)),
  data = trait_df)[, -1, drop = FALSE]
log_info("design_tech: {nrow(design_tech)} × {ncol(design_tech)}")

resid_grade <- residuals(lm(trait_df$grade ~ design_tech))

partial_df <- data.frame()
for (mod_col in colnames(MEs_global)) {
  mod_name <- gsub("^ME", "", mod_col)
  me_vec <- as.numeric(MEs_global[, mod_col])
  fit_me <- lm(me_vec ~ design_tech)
  resid_me <- residuals(fit_me)
  ct <- cor.test(resid_me, resid_grade, use = "complete.obs")
  partial_df <- rbind(partial_df, data.frame(
    Module = mod_name,
    partial_r = round(ct$estimate, 3),
    p = ct$p.value,
    stringsAsFactors = FALSE))
}
partial_df$FDR <- p.adjust(partial_df$p, method = "BH")
partial_df$sig <- partial_df$FDR < PADJ_THRESHOLD
partial_df <- partial_df[order(partial_df$p), ]
save_table(partial_df, DIR_THIS, "module_trait_partial.csv")

# ============================================================
# ============================================================
log_info("--- 7. 模块 vs 生物学轴 ---")

key_traits <- c("hUSI","ML_Aging","genes16","mito")
key_cor <- cor_long[cor_long$Trait %in% key_traits, ]

hUSI_cor <- cor_long[cor_long$Trait == "hUSI", c("Module","r","FDR")]
colnames(hUSI_cor)[2:3] <- c("r_hUSI", "FDR_hUSI")
save_table(hUSI_cor, DIR_THIS, "module_hUSI_correlation.csv")

g16_cor <- cor_long[cor_long$Trait == "genes16", c("Module","r","FDR")]
colnames(g16_cor)[2:3] <- c("r_16genes", "FDR_16genes")
save_table(g16_cor, DIR_THIS, "module_vs_16genes.csv")

mito_cor <- cor_long[cor_long$Trait == "mito", c("Module","r","FDR")]
colnames(mito_cor)[2:3] <- c("r_mito", "FDR_mito")
save_table(mito_cor, DIR_THIS, "module_vs_mito.csv")

biology_wide <- Reduce(function(x, y) merge(x, y, by = "Module", all = TRUE),
                       list(hUSI_cor, g16_cor, mito_cor))
save_table(biology_wide, DIR_THIS, "module_vs_biology.csv")

# ============================================================
# ============================================================
log_info("--- 8. 模块基因集富集 ---")

background <- rownames(expr)
bg_sym <- toupper(sym_map[strip_version(background)])
bg_sym <- bg_sym[!is.na(bg_sym) & bg_sym != ""]
bg_sym <- unique(bg_sym)

enrich_df <- data.frame()
for (mod in all_modules) {
  mod_ens <- module_df$ENSEMBL[module_df$Module == mod]
  mod_sym <- toupper(sym_map[strip_version(mod_ens)])
  mod_sym <- mod_sym[!is.na(mod_sym) & mod_sym != ""]
  mod_sym <- unique(mod_sym)

  if (length(mod_sym) < 10) next

  res <- run_enrichment(mod_sym, GENE_SETS, bg_sym)
  if (nrow(res) > 0) {
    res$Module <- mod
    enrich_df <- rbind(enrich_df, res)
  }
}
if (nrow(enrich_df) > 0) {
  enrich_df <- enrich_df[order(enrich_df$p), ]
}
save_table(enrich_df, DIR_THIS, "module_enrichment.csv")

# ============================================================
# ============================================================
log_info("--- 9. 基因模块归属 ---")

module_df$SYMBOL_FIXED <- sym_map[strip_version(module_df$ENSEMBL)]
module_df$SYMBOL_FIXED[is.na(module_df$SYMBOL_FIXED)] <-
  module_df$ENSEMBL[is.na(module_df$SYMBOL_FIXED)]

mito_key_genes <- c("NDUFA1","NDUFB8","NDUFS1","COX4I1","COX5A",
                     "ATP5F1A","SDHA","UQCRC1","NDUFA2","NDUFV2",
                     "UQCRFS1","COX7C","ATP5PB")

check_genes <- unique(c(GENES_16, toupper(mito_key_genes)))

hits_df <- data.frame()
for (g in check_genes) {
  mod_hit <- module_df$Module[toupper(module_df$SYMBOL_FIXED) == toupper(g)]
  if (length(mod_hit) > 0) {
    hits_df <- rbind(hits_df, data.frame(
      Gene = g, Module = mod_hit[1],
      Type = ifelse(g %in% GENES_16, "16gene", "mito"),
      stringsAsFactors = FALSE))
  }
}
save_table(hits_df, DIR_THIS, "module_gene_hits.csv")

# ============================================================
# 10. module-modulecorrelation
# ============================================================
log_info("--- 10. 模块-模块相关 ---")

mod_cor <- cor(MEs_global, use = "pairwise.complete.obs")
save_table(as.data.frame(round(mod_cor, 3)),
            DIR_THIS, "module_module_correlation.csv")

pdf(file.path(DIR_THIS, "figures", "module_module_heatmap.pdf"),
    width = 10, height = 10)
pheatmap(mod_cor,
         cluster_rows = TRUE, cluster_cols = TRUE,
         display_numbers = FALSE,
         color = colorRampPalette(c("#4C72B0","#FFFFFF","#C44E52"))(50),
         main = "Module-module correlation")
dev.off()

# ============================================================
# 11. modulesummary
# ============================================================
log_info("--- 11. 模块汇总 ---")

get_r <- function(mod, trait) {
  v <- cor_long$r[cor_long$Module == mod & cor_long$Trait == trait]
  if (length(v) > 0) as.numeric(v[1]) else NA_real_
}
get_fdr <- function(mod, trait) {
  v <- cor_long$FDR[cor_long$Module == mod & cor_long$Trait == trait]
  if (length(v) > 0) as.numeric(v[1]) else NA_real_
}

summary_df <- data.frame()
for (mod in all_modules) {
  mod_size <- sum(module_df$Module == mod)
  if (!paste0("ME", mod) %in% colnames(MEs_global)) next

  top_enrich <- enrich_df[enrich_df$Module == mod, ]
  top_enrich <- top_enrich[order(top_enrich$p), ]
  top1 <- if (nrow(top_enrich) > 0) as.character(top_enrich$Geneset[1]) else NA_character_
  top1_fdr <- if (nrow(top_enrich) > 0) as.numeric(top_enrich$FDR[1]) else NA_real_

  summary_df <- rbind(summary_df, data.frame(
    Module         = mod,
    N_genes        = as.integer(mod_size),
    r_grade        = get_r(mod, "grade"),
    FDR_grade      = get_fdr(mod, "grade"),
    r_hUSI         = get_r(mod, "hUSI"),
    FDR_hUSI       = get_fdr(mod, "hUSI"),
    r_ML_Aging     = get_r(mod, "ML_Aging"),
    FDR_ML_Aging   = get_fdr(mod, "ML_Aging"),
    r_mito         = get_r(mod, "mito"),
    FDR_mito       = get_fdr(mod, "mito"),
    r_genes16      = get_r(mod, "genes16"),
    FDR_genes16    = get_fdr(mod, "genes16"),
    r_batch        = get_r(mod, "batch"),
    r_layer        = get_r(mod, "layer"),
    r_rin          = get_r(mod, "rin"),
    top_enrich     = top1,
    top_enrich_FDR = top1_fdr,
    stringsAsFactors = FALSE))
}
summary_df <- summary_df[order(-abs(summary_df$r_grade)), ]
save_table(summary_df, DIR_THIS, "module_summary_annotated.csv")

# ============================================================
# 12. figure
# ============================================================
log_info("--- 12. 图 ---")

cor_wide <- reshape(cor_long[, c("Module","Trait","r")],
                     idvar = "Module", timevar = "Trait", direction = "wide")
rownames(cor_wide) <- cor_wide$Module
cor_wide$Module <- NULL
colnames(cor_wide) <- gsub("^r\\.", "", colnames(cor_wide))

keep_traits <- c("grade","hUSI","ML_Aging","genes16","mito",
                   "age","rin","pmi","batch","layer","sex_M")
cor_wide <- cor_wide[, intersect(keep_traits, colnames(cor_wide)), drop = FALSE]

pheatmap(as.matrix(cor_wide),
         cluster_rows = TRUE, cluster_cols = TRUE,
         display_numbers = TRUE, number_format = "%.2f",
         color = colorRampPalette(c("#4C72B0","#FFFFFF","#C44E52"))(50),
         main = "Module-trait correlation",
         filename = file.path(DIR_THIS, "figures",
                              "module_trait_heatmap.png"),
         width = 9, height = 6)

top8 <- head(summary_df$Module[order(-abs(summary_df$r_grade))], 8)

curve_df <- data.frame()
for (mod in top8) {
  me_col <- paste0("ME", mod)
  if (!me_col %in% colnames(MEs_global)) next
  for (g in sort(unique(trait_df$grade))) {
    idx <- trait_df$grade == g
    vals <- MEs_global[idx, me_col]
    curve_df <- rbind(curve_df, data.frame(
      Module = mod, Grade = g,
      mean_ME = mean(vals, na.rm = TRUE),
      se = sd(vals, na.rm = TRUE) / sqrt(sum(idx)),
      stringsAsFactors = FALSE))
  }
}

p <- ggplot(curve_df, aes(x = Grade, y = mean_ME,
                            color = Module, group = Module)) +
  geom_line(linewidth = 0.9) + geom_point(size = 2.5) +
  geom_errorbar(aes(ymin = mean_ME - se, ymax = mean_ME + se),
                 width = 0.15, alpha = 0.5) +
  theme_bw(base_size = 11) +
  scale_x_continuous(breaks = 1:4, labels = paste0("G", 1:4)) +
  labs(x = "AMD grade", y = "Module eigengene (mean ± SE)",
       title = "Top 8 modules by |r_grade|")
save_fig(p, DIR_THIS, "module_grade_curves", width = 10, height = 6)

# ============================================================
# 13. validation：selectsignificantmodule
# ============================================================
log_info("--- 13. 验证：选择显著模块 ---")

targets <- summary_df$Module[summary_df$FDR_grade < 0.10]
if (length(targets) == 0) {
  log_warn("FDR_grade < 0.10 无模块，改用 |r_grade| 前 5")
  targets <- head(summary_df$Module[order(-abs(summary_df$r_grade))], 5)
}
log_info("待验证模块: {paste(targets, collapse=', ')}")

if (length(targets) == 0) {
  log_warn("无验证模块，跳过 L0-L6")
} else {
  val_df <- data.frame(
    Sample = rownames(MEs_global),
    grade  = trait_df$grade,
    age    = trait_df$age,
    RIN    = trait_df$rin,
    PMI    = trait_df$pmi,
    stringsAsFactors = FALSE
  )

  # ============================================================
  # 14. L0: PC1-5 ~ covariates
  # ============================================================
  log_info("--- 14. L0: PC1-5 ~ covariates ---")
  pc_cols <- intersect(paste0("PC", 1:5), colnames(clin))
  if (length(pc_cols) > 0) {
    L0 <- data.frame()
    for (pc in pc_cols) {
      for (v in c("grade","age","batch","RIN","PMI","seq_date")) {
        if (!v %in% colnames(clin)) next
        r <- suppressWarnings(cor(as.numeric(clin[[pc]]),
                                  as.numeric(factor(clin[[v]])),
                                  use = "complete.obs"))
        L0 <- rbind(L0, data.frame(PC = pc, covariate = v, r = r))
      }
    }
    save_table(L0, DIR_THIS, "validation_L0_PC_vs_covariates.csv")
  } else {
    log_warn("clin 中无 PC1-5，跳过 L0")
  }

  # ============================================================
  # 15. L1: ME ~ seq_date
  # ============================================================
  if ("seq_date" %in% colnames(clin)) {
    log_info("--- 15. L1: ME ~ seq_date ---")
    L1 <- data.frame()
    for (m in targets) {
      me_col <- paste0("ME", m)
      if (!me_col %in% colnames(MEs_global)) next
      df <- data.frame(y = MEs_global[[me_col]], seq = clin$seq_date)
      a <- anova(lm(y ~ seq, data = df))
      L1 <- rbind(L1, data.frame(
        module = m,
        R2 = a$`Sum Sq`[1] / sum(a$`Sum Sq`),
        p = a$`Pr(>F)`[1]))
    }
    save_table(L1, DIR_THIS, "validation_L1_ME_vs_seqdate.csv")
  } else {
    log_warn("clin 中无 seq_date，跳过 L1-L3")
  }

  # ============================================================
  # ============================================================
  if ("seq_date" %in% colnames(clin)) {
    log_info("--- 16. L2: 分层 ---")
    L2 <- data.frame()
    for (d in unique(clin$seq_date)) {
      idx <- clin$seq_date == d
      if (sum(idx) < 30) next
      for (m in targets) {
        me_col <- paste0("ME", m)
        if (!me_col %in% colnames(MEs_global)) next
        df <- data.frame(
          y = MEs_global[idx, me_col],
          grade = clin$grade[idx],
          age = clin$age_n[idx],
          RIN = clin$rin_n[idx],
          PMI = clin$pmi_n[idx])
        fit <- tryCatch(lm(y ~ grade + age + RIN + PMI, data = df),
                        error = function(e) NULL)
        if (is.null(fit)) next
        co <- summary(fit)$coefficients
        if ("grade" %in% rownames(co)) {
          L2 <- rbind(L2, data.frame(
            module = m, seq_date = d,
            beta = co["grade", "Estimate"],
            p = co["grade", "Pr(>|t|)"]))
        }
      }
    }
    save_table(L2, DIR_THIS, "validation_L2_stratified.csv")
  }

  # ============================================================
  # 17. L3: LODO
  # ============================================================
  if ("seq_date" %in% colnames(clin)) {
    log_info("--- 17. L3: LODO ---")
    L3 <- data.frame()
    for (d in unique(clin$seq_date)) {
      idx <- clin$seq_date != d
      for (m in targets) {
        me_col <- paste0("ME", m)
        if (!me_col %in% colnames(MEs_global)) next
        df <- data.frame(
          y = MEs_global[idx, me_col],
          grade = clin$grade[idx],
          age = clin$age_n[idx],
          RIN = clin$rin_n[idx],
          PMI = clin$pmi_n[idx])
        fit <- lm(y ~ grade + age + RIN + PMI, data = df)
        co <- summary(fit)$coefficients
        L3 <- rbind(L3, data.frame(
          module = m, removed = d,
          beta = co["grade", "Estimate"],
          p = co["grade", "Pr(>|t|)"]))
      }
    }
    save_table(L3, DIR_THIS, "validation_L3_LODO.csv")
  }

  # ============================================================
  # 18. L4: permutation
  # ============================================================
  log_info("--- 18. L4: permutation ---")
  L4 <- data.frame()
  for (m in targets) {
    me_col <- paste0("ME", m)
    if (!me_col %in% colnames(MEs_global)) next
    df <- data.frame(
      y = MEs_global[[me_col]],
      grade = trait_df$grade,
      age = trait_df$age,
      RIN = trait_df$rin,
      PMI = trait_df$pmi)
    obs <- coef(lm(y ~ grade + age + RIN + PMI, data = df))["grade"]
    null <- numeric(N_PERM)
    for (i in seq_len(N_PERM)) {
      df$gp <- sample(df$grade)
      null[i] <- coef(lm(y ~ gp + age + RIN + PMI, data = df))["gp"]
    }
    L4 <- rbind(L4, data.frame(
      module = m, obs_beta = obs,
      p_perm = mean(abs(null) >= abs(obs)),
      n_perm = N_PERM))
  }
  save_table(L4, DIR_THIS, "validation_L4_permutation_summary.csv")
  print(L4, row.names = FALSE)

  # ============================================================
  # 19. L6: panel enrichment
  # ============================================================
  log_info("--- 19. L6: panel 富集 ---")
  L6 <- data.frame()
  for (m in targets) {
    mod_ens <- module_df$ENSEMBL[module_df$Module == m]
    mod_sym <- toupper(sym_map[strip_version(mod_ens)])
    mod_sym <- mod_sym[!is.na(mod_sym) & mod_sym != ""]
    mod_sym <- unique(mod_sym)
    if (length(mod_sym) < 5) next

    res <- run_enrichment(mod_sym, GENE_SETS, bg_sym)
    if (nrow(res) > 0) {
      res$module <- m
      L6 <- rbind(L6, res)
    }
  }
  if (nrow(L6) > 0) {
    L6 <- L6[order(L6$p), ]
  }
  save_table(L6, DIR_THIS, "validation_L6_panel_enrich.csv")
}

# ============================================================
# 20. summary
# ============================================================
log_info("--- 20. 汇总 ---")

cat("\n", strrep("=", 70), "\n", sep = "")
cat("=== 05_2 汇总（主分析 + 验证）===\n")
cat(strrep("=", 70), "\n", sep = "")

cat("\n【Q1】模块-性状相关\n")
cat("  总测试:", nrow(cor_long), "\n")
cat("  显著(FDR<0.05):", sum(cor_long$FDR < 0.05, na.rm = TRUE), "\n")

cat("\n【Q2】模块 vs grade ANOVA\n")
sig_anova <- anova_df[anova_df$FDR < 0.05, ]
cat("  FDR<0.05:", nrow(sig_anova), "\n")

cat("\n【Q3】偏相关\n")
cat("  FDR<0.05:", sum(partial_df$FDR < 0.05, na.rm = TRUE), "\n")

cat("\n【Q4】16 基因命中\n")
cat("  命中:", length(genes_16_in), " / ", length(GENES_16), "\n")
if (length(missing_16) > 0) {
  cat("  未命中:", paste(missing_16, collapse = ", "), "\n")
}

cat("\n【Q5】关键基因归属\n")
print(hits_df, row.names = FALSE)

if (exists("L4")) {
  cat("\n【Q6】验证 L4 permutation\n")
  print(L4, row.names = FALSE)
}

cat("\n--- 输出位置 ---\n")
cat("目录:", DIR_THIS, "\n")

log_info("=== 05_2_module_trait_and_validation.r done ===")

