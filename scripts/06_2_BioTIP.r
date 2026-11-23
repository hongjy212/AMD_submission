#!/usr/bin/env Rscript
# ============================================================
# 06_2_BioTIP.r
#
#
# ★ v2 Fix：
#
#   - covariate：age + sex + batch + layer + RIN + PMI
#
# Input：
#   results/01_1b_batch_enhance/data/{expr_wgcna,clinical_wgcna}.rds
#   results/01_1_prepare_bulk_data/data/gene_annotation.rds
#   results/01_2_prepare_gene_sets/data/gene_sets.rds
#   results/06_1_TRADE_vs_adjacent/tables/poly3F_gene_level.csv（optional）
#
#   tables/MCI_scores_<comp>.csv
#   tables/CTS_genes_<comp>.csv
#   tables/MCI_summary_across.csv
#   tables/CTS_artifact_check.csv
#   tables/IEG_mediation_analysis.csv
#   tables/CTS_poly3F_validation.csv
#   tables/sd_selection_diagnostics.csv
#   tables/correction_residual_diagnostics_subset.csv
#   figures/MCI_barplot_<comp>.pdf/.png
#   figures/MCI_comparison_across.pdf/.png
#   figures/PCA_subset_<comp>.pdf/.png
#   figures/CTS_poly3F_enrichment.pdf/.png
# ============================================================

if (exists(".AMD_SETUP_LOADED", envir = .GlobalEnv)) {
  rm(".AMD_SETUP_LOADED", envir = .GlobalEnv)
}
source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")

suppressPackageStartupMessages({
  library(BioTIP)
  library(dplyr); library(tidyr); library(ggplot2)
  library(limma); library(splines); library(pheatmap)
})

DIR_THIS <- setup_output_dir("06_2_BioTIP")
log_info("=== 06_2_BioTIP.r (v2) start ===")

# ============================================================
# ============================================================
SD_CUTOFF       <- 0.01
MIN_MODULE_SIZE <- 20
NETWORK_FDR     <- 1
ARTIFACT_R      <- 0.5
P_THRESHOLD     <- 0.05
BATCH_R2_WARN   <- 0.1
SEX_R2_WARN     <- 0.1
RESIDUALIZE_BEFORE_BIOTIP <- TRUE   # ★ v2：BioTIP 前是否 residualize

# ============================================================
# ============================================================
log_info("--- 1. 载入 ---")

expr <- read_upstream("01_1b_batch_enhance", "data", "expr_wgcna.rds")
clin <- read_upstream("01_1b_batch_enhance", "data", "clinical_wgcna.rds")
anno <- read_upstream("01_1_prepare_bulk_data", "data", "gene_annotation.rds")
GENE_SETS <- read_upstream("01_2_prepare_gene_sets", "data", "gene_sets.rds")

common <- intersect(colnames(expr), rownames(clin))
expr <- expr[, common, drop = FALSE]
clin <- clin[common, , drop = FALSE]

log_info("expr: {nrow(expr)} × {ncol(expr)}")

# ============================================================
# ============================================================
log_info("--- 2. 整理变量 ---")

clin$grade <- as.integer(clin$grade)
clin$age_v <- as.numeric(clin$age_n)
clin$rin_v <- as.numeric(clin$rin_n)
clin$pmi_v <- as.numeric(clin$pmi_n)
clin$batch_f <- factor(clin$batch)
sid <- as.character(clin$sample_id)
clin$layer_f <- factor(sub("^.*_(L[0-9]+)$", "\\1", sid),
                        levels = paste0("L", 1:8))
clin$sex_f <- factor(clin$sex)
rownames(clin) <- colnames(expr)
clin$sample_id <- rownames(clin)

log_info("grade: {paste(names(table(clin$grade)), table(clin$grade), sep='=', collapse=' ')}")

# ============================================================
# ============================================================
log_info("--- 3. 行名 → SYMBOL ---")

strip_version <- function(x) sub("\\..*$", "", x)
sym_map <- setNames(toupper(anno$symbol), strip_version(anno$ensembl))

expr_mat <- expr
rownames(expr_mat) <- sym_map[strip_version(rownames(expr_mat))]

if (any(duplicated(rownames(expr_mat)))) {
  vars <- apply(expr_mat, 1, var, na.rm = TRUE)
  vars[is.na(vars)] <- -Inf
  expr_mat <- expr_mat[order(rownames(expr_mat), -vars), ]
  expr_mat <- expr_mat[!duplicated(rownames(expr_mat)), ]
}
expr_mat <- as.matrix(expr_mat)
log_info("expr_mat: {nrow(expr_mat)} × {ncol(expr_mat)}")

# ============================================================
# ============================================================
log_info("--- 4. BioTIP 工具函数 ---")

extract_top_mci_per_state <- function(MCI_list) {
  empty_df <- data.frame(
    State = character(0), TopCommunity = character(0),
    MCI = numeric(0), n_communities = integer(0),
    stringsAsFactors = FALSE)
  if (is.null(MCI_list) || length(MCI_list) == 0) return(empty_df)

  rows <- list()
  for (state_name in names(MCI_list)) {
    state_mci <- MCI_list[[state_name]]
    if (is.null(state_mci) || length(state_mci) == 0) next

    if (is.list(state_mci)) {
      vals <- unlist(state_mci, use.names = TRUE)
    } else {
      vals <- as.numeric(state_mci); names(vals) <- names(state_mci)
    }
    if (length(vals) == 0 || all(is.na(vals))) next

    rows[[state_name]] <- data.frame(
      State = state_name,
      TopCommunity = names(vals)[which.max(vals)],
      MCI = max(vals, na.rm = TRUE),
      n_communities = length(state_mci),
      stringsAsFactors = FALSE)
  }
  if (length(rows) == 0) return(empty_df)
  MCI_df <- do.call(rbind, rows)
  MCI_df <- MCI_df[order(-MCI_df$MCI), ]
  rownames(MCI_df) <- NULL
  MCI_df
}

flatten_CTS <- function(CTS_res) {
  if (is.null(CTS_res)) return(NULL)
  if (is.character(CTS_res) && !is.list(CTS_res)) {
    return(CTS_res[!is.na(CTS_res) & CTS_res != ""])
  }
  if (is.list(CTS_res)) {
    flat <- unique(unlist(CTS_res, use.names = FALSE))
    return(flat[!is.na(flat) & flat != ""])
  }
  as.character(CTS_res)
}

residualize_for_biotip <- function(expr_sub, clin_sub) {
  sp_age <- splines::ns(clin_sub$age_v, df = 3)
  covar <- data.frame(
    age_s1 = sp_age[, 1],
    age_s2 = sp_age[, 2],
    age_s3 = sp_age[, 3],
    sex_M  = as.integer(clin_sub$sex_f == "Male"),
    rin    = clin_sub$rin_v,
    pmi    = clin_sub$pmi_v,
    stringsAsFactors = FALSE)

  design <- model.matrix(
    ~ age_s1 + age_s2 + age_s3 + sex_M + rin + pmi,
    data = covar)

  fit <- limma::lmFit(expr_sub, design)
  resid <- residuals(fit, expr_sub)
  list(resid = resid, design_ncol = ncol(design),
       n_samples = nrow(design))
}

run_biotip_core <- function(expr_input, samples_sub, state_vec,
                             sd_cutoff, min_module_size,
                             network_fdr = 1) {
  res <- list(MCI_df = NULL, MCI_res = NULL, CTS_res = NULL,
              igraphL = NULL, cluster = NULL,
              status = "ok", msg = NA_character_,
              n_selected = NA_integer_,
              n_selected_per_group = NA_character_)

  expr_sub <- expr_input[, samples_sub, drop = FALSE]

  # 1) sd_selection
  selected_matrix_list <- tryCatch({
    sd_selection(expr_sub,
                 samplesL = split(samples_sub, state_vec),
                 cutoff   = sd_cutoff)
  }, error = function(e) {
    res$status <<- "skip"
    res$msg    <<- paste0("sd_selection failed: ", e$message)
    NULL
  })
  if (is.null(selected_matrix_list)) return(res)

  n_per_group <- sapply(selected_matrix_list, function(x) {
    if (is.matrix(x) || is.data.frame(x)) return(nrow(x))
    length(x)
  })
  res$n_selected_per_group <- paste(names(selected_matrix_list), n_per_group,
                                     sep = "=", collapse = "; ")
  res$n_selected <- sum(n_per_group)

  if (sum(n_per_group) < 50) {
    res$status <- "skip"; res$msg <- "selected_genes < 50"; return(res)
  }

  # 2) getNetwork
  igraphL <- tryCatch(
    getNetwork(selected_matrix_list, fdr = network_fdr),
    error = function(e) {
      res$status <<- "skip"
      res$msg    <<- paste0("getNetwork failed: ", e$message)
      NULL
    })
  if (is.null(igraphL)) return(res)
  res$igraphL <- igraphL

  # 3) getCluster_methods
  cluster <- tryCatch(
    getCluster_methods(igraphL),
    error = function(e) {
      res$status <<- "skip"
      res$msg    <<- paste0("getCluster_methods failed: ", e$message)
      NULL
    })
  if (is.null(cluster)) return(res)
  res$cluster <- cluster

  # 4) getMCI
  membersL <- tryCatch(
    getMCI(cluster, selected_matrix_list),
    error = function(e) {
      res$status <<- "skip"
      res$msg    <<- paste0("getMCI failed: ", e$message)
      NULL
    })
  if (is.null(membersL)) return(res)
  res$MCI_res <- membersL
  res$MCI_df <- extract_top_mci_per_state(membersL[["MCI"]])

  # 5) getMaxMCImember
  maxMCIms <- tryCatch(
    getMaxMCImember(membersL[["members"]], membersL[["MCI"]],
                    min = min_module_size),
    error = function(e) {
      log_warn("getMaxMCImember failed: {e$message}，尝试回退")
      NULL
    })

  if (is.null(maxMCIms)) {
    log_warn("使用手动回退构造 maxMCIms")
    MCI_list <- membersL[["MCI"]]
    members_list <- membersL[["members"]]

    idx_list <- list(); members_top <- list()
    for (state_name in names(MCI_list)) {
      state_mci <- MCI_list[[state_name]]
      if (is.null(state_mci) || length(state_mci) == 0) next

      if (is.list(state_mci)) {
        vals <- unlist(state_mci, use.names = TRUE)
      } else {
        vals <- as.numeric(state_mci); names(vals) <- names(state_mci)
      }

      state_membership <- members_list[[state_name]]
      if (is.null(state_membership)) next

      membership_vec <- as.numeric(state_membership)
      names(membership_vec) <- names(state_membership)
      comm_sizes <- table(membership_vec)
      valid_comms <- names(comm_sizes)[comm_sizes >= min_module_size]
      if (length(valid_comms) == 0) next

      vals_valid <- vals[names(vals) %in% valid_comms]
      if (length(vals_valid) == 0) next

      top_comm <- names(vals_valid)[which.max(vals_valid)]
      idx_list[[state_name]] <- top_comm
      members_top[[state_name]] <- names(membership_vec)[
        membership_vec == as.numeric(top_comm)]
    }

    if (length(idx_list) > 0) {
      maxMCIms <- list(idx = idx_list, members = members_top)
    } else {
      res$status <- "partial"; res$msg <- "no valid community"
      return(res)
    }
  }

  # 6) getMaxStats
  maxMCI <- tryCatch(
    getMaxStats(membersL[["MCI"]], maxMCIms[["idx"]]),
    error = function(e) {
      log_warn("getMaxStats failed: {e$message}，尝试回退")
      NULL
    })

  if (is.null(maxMCI)) {
    MCI_list <- membersL[["MCI"]]
    maxMCI <- sapply(names(maxMCIms[["idx"]]), function(state_name) {
      state_mci <- MCI_list[[state_name]]
      if (is.list(state_mci)) state_mci <- unlist(state_mci, use.names = TRUE)
      as.numeric(state_mci[maxMCIms[["idx"]][[state_name]]])
    })
    names(maxMCI) <- names(maxMCIms[["idx"]])
  }

  # 7) getCTS
  CTS_raw <- tryCatch(
    getCTS(maxMCI, maxMCIms[["members"]]),
    error = function(e) {
      log_warn("getCTS failed: {e$message}，尝试回退")
      NULL
    })

  if (is.list(CTS_raw) && !is.character(CTS_raw)) {
    log_info("  getCTS 返回 list of {length(CTS_raw)} states，展平中...")
  }
  CTS_res <- flatten_CTS(CTS_raw)

  if (is.null(CTS_res) || length(CTS_res) == 0) {
    CTS_res <- unique(unlist(maxMCIms[["members"]]))
    CTS_res <- CTS_res[!is.na(CTS_res) & CTS_res != ""]
    if (length(CTS_res) == 0) CTS_res <- NULL
  }

  res$CTS_res <- CTS_res

  if (!is.null(res$MCI_df) && nrow(res$MCI_df) > 0) {
    for (k in seq_len(nrow(res$MCI_df))) {
      log_info("    {res$MCI_df$State[k]}: top comm={res$MCI_df$TopCommunity[k]}, MCI={round(res$MCI_df$MCI[k], 3)}")
    }
  }
  if (!is.null(res$CTS_res)) {
    log_info("  CTS 基因数: {length(res$CTS_res)}")
  }

  res
}

# ============================================================
# ============================================================
log_info("--- 5. BioTIP 分析 ---")

comparisons <- list(
  list(label = "2_vs_1", low = 1, high = 2),
  list(label = "3_vs_2", low = 2, high = 3),
  list(label = "4_vs_3", low = 3, high = 4)
)

all_MCI <- list()
all_CTS <- list()
diag_sub_df <- data.frame()
sd_diag_df  <- data.frame()
resid_diag_df <- data.frame()

for (comp in comparisons) {
  label <- comp$label
  low   <- comp$low
  high  <- comp$high

  set.seed(SEED)

  log_info("========== {label} (MGS{low} vs MGS{high}) ==========")

  keep_idx <- which(clin$grade %in% c(low, high))
  samples_sub <- rownames(clin)[keep_idx]
  expr_sub <- expr_mat[, samples_sub, drop = FALSE]
  clin_sub <- clin[keep_idx, ]

  state_vec <- ifelse(clin$grade[keep_idx] == high,
                       paste0("MGS", high), paste0("MGS", low))
  names(state_vec) <- samples_sub

  n_low  <- sum(state_vec == paste0("MGS", low))
  n_high <- sum(state_vec == paste0("MGS", high))
  log_info("  MGS{low}: {n_low}, MGS{high}: {n_high}")

  log_info("  [5a] 子集内校正残留诊断")
  pca_sub <- prcomp(t(expr_sub), center = TRUE, scale. = FALSE)
  var_pct_sub <- summary(pca_sub)$importance[2, 1:min(5, ncol(pca_sub))] * 100

  r2_grade_sub <- {
    fit <- aov(pca_sub$x[, 1] ~ factor(clin$grade[keep_idx]))
    summary(fit)[[1]][1, "Sum Sq"] / sum(summary(fit)[[1]][, "Sum Sq"])
  }

  r2_batch_sub <- sapply(c("batch_f", "layer_f"), function(bn) {
    fit <- tryCatch(
      aov(pca_sub$x[, 1] ~ clin[[bn]][keep_idx]),
      error = function(e) NULL)
    if (is.null(fit)) return(NA_real_)
    summary(fit)[[1]][1, "Sum Sq"] / sum(summary(fit)[[1]][, "Sum Sq"])
  })

  r2_sex_sub <- NA_real_
  fit_sex <- tryCatch(
    aov(pca_sub$x[, 1] ~ clin$sex_f[keep_idx]),
    error = function(e) NULL)
  if (!is.null(fit_sex)) {
    r2_sex_sub <- summary(fit_sex)[[1]][1, "Sum Sq"] /
                  sum(summary(fit_sex)[[1]][, "Sum Sq"])
  }

  log_info("  [{label}] PC1 方差: {round(var_pct_sub[1], 1)}%")
  log_info("  [{label}] PC1 ~ grade R²: {round(r2_grade_sub, 3)}")
  log_info("  [{label}] PC1 ~ batch R²: {round(r2_batch_sub['batch_f'], 3)}")
  log_info("  [{label}] PC1 ~ layer R²: {round(r2_batch_sub['layer_f'], 3)}")
  log_info("  [{label}] PC1 ~ sex R²: {round(r2_sex_sub, 3)}")

  diag_sub_df <- rbind(diag_sub_df, data.frame(
    Comparison    = label,
    n_low         = n_low,
    n_high        = n_high,
    PC1_var_pct   = round(var_pct_sub[1], 2),
    PC1_grade_R2  = round(r2_grade_sub, 3),
    PC1_batch_R2  = round(r2_batch_sub["batch_f"], 3),
    PC1_layer_R2  = round(r2_batch_sub["layer_f"], 3),
    PC1_sex_R2    = round(r2_sex_sub, 3),
    stringsAsFactors = FALSE))

  # PCA figure
  pca_df <- data.frame(
    PC1   = pca_sub$x[, 1],
    PC2   = if (ncol(pca_sub$x) >= 2) pca_sub$x[, 2] else 0,
    batch = factor(clin$batch_f[keep_idx]),
    grade = factor(clin$grade[keep_idx])
  )
  p_pca <- ggplot(pca_df, aes(PC1, PC2, color = batch, shape = grade)) +
    geom_point(size = 2.5, alpha = 0.8) +
    labs(title = paste0("Subset PCA - ", label),
         subtitle = sprintf("PC1~batch R2=%.3f, PC1~grade R2=%.3f, PC1~sex R2=%.3f",
                            r2_batch_sub["batch_f"],
                            r2_grade_sub, r2_sex_sub)) +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 9))
  save_fig(p_pca, DIR_THIS, paste0("PCA_subset_", label), width = 6, height = 4.5)

  # --- 5b-pre. ★ v2：residualize age + sex + RIN + PMI（keep grade） ---
  if (RESIDUALIZE_BEFORE_BIOTIP) {
    log_info("  [5b-pre] residualize age + sex + RIN + PMI")

    resid_res <- residualize_for_biotip(expr_sub, clin_sub)
    expr_sub_resid <- resid_res$resid

    sd_before <- apply(expr_sub, 1, sd, na.rm = TRUE)
    sd_after  <- apply(expr_sub_resid, 1, sd, na.rm = TRUE)

    log_info("  [5b-pre] 原始 SD 中位数: {round(median(sd_before, na.rm=TRUE), 3)}")
    log_info("  [5b-pre] 残差 SD 中位数: {round(median(sd_after, na.rm=TRUE), 3)}")
    log_info("  [5b-pre] design 列数: {resid_res$design_ncol}")

    resid_diag_df <- rbind(resid_diag_df, data.frame(
      Comparison      = label,
      SD_median_raw   = round(median(sd_before, na.rm = TRUE), 3),
      SD_median_resid = round(median(sd_after, na.rm = TRUE), 3),
      SD_ratio        = round(median(sd_after, na.rm = TRUE) /
                              median(sd_before, na.rm = TRUE), 3),
      design_ncol     = resid_res$design_ncol,
      stringsAsFactors = FALSE))

    expr_sub_for_biotip <- expr_sub_resid
  } else {
    log_warn("  [5b-pre] 跳过 residualize（RESIDUALIZE_BEFORE_BIOTIP = FALSE）")
    expr_sub_for_biotip <- expr_sub
  }

  log_info("  [5b] BioTIP 主流程")

  sd_all <- apply(expr_sub_for_biotip, 1, sd, na.rm = TRUE)
  log_info("  子集内基因 SD: median={round(median(sd_all, na.rm=TRUE),3)}")

  res_main <- run_biotip_core(expr_sub_for_biotip, samples_sub, state_vec,
                              SD_CUTOFF, MIN_MODULE_SIZE, NETWORK_FDR)

  sd_diag_df <- rbind(sd_diag_df, data.frame(
    Comparison = label, Version = "main",
    n_selected = res_main$n_selected,
    n_selected_per_group = res_main$n_selected_per_group,
    status = res_main$status,
    msg = res_main$msg,
    residualized = RESIDUALIZE_BEFORE_BIOTIP,
    stringsAsFactors = FALSE))

  if (res_main$status == "skip") {
    log_warn("  [{label}] 主版本跳过：{res_main$msg}")
    next
  }

  if (!is.null(res_main$MCI_df) && nrow(res_main$MCI_df) > 0) {
    write.csv(res_main$MCI_df,
              file.path(DIR_THIS, "tables", paste0("MCI_scores_", label, ".csv")),
              row.names = FALSE)
    all_MCI[[label]] <- res_main$MCI_res

    p_bar <- ggplot(res_main$MCI_df,
                    aes(x = reorder(State, -MCI), y = MCI, fill = State)) +
      geom_col(width = 0.6) +
      geom_text(aes(label = round(MCI, 3)), vjust = -0.3, size = 3.5) +
      scale_fill_manual(values = c("MGS2" = "#27ae60", "MGS3" = "#c0392b",
                                    "MGS1" = "#3498db", "MGS4" = "#8e44ad")) +
      labs(title = paste("MCI scores -", label),
           x = "State", y = "MCI score") +
      theme_bw(base_size = 12) +
      theme(legend.position = "none",
            plot.title = element_text(hjust = 0.5, face = "bold"))
    save_fig(p_bar, DIR_THIS, paste0("MCI_barplot_", label), width = 6, height = 4)
  }

  if (!is.null(res_main$CTS_res) && length(res_main$CTS_res) > 0) {
    write.csv(data.frame(Gene = res_main$CTS_res),
              file.path(DIR_THIS, "tables", paste0("CTS_genes_", label, ".csv")),
              row.names = FALSE)
    all_CTS[[label]] <- res_main$CTS_res
    log_info("  [{label}] CTS 基因数: {length(res_main$CTS_res)}")
  }
}

if (nrow(diag_sub_df) > 0) {
  save_table(diag_sub_df, DIR_THIS, "correction_residual_diagnostics_subset.csv")
}
if (nrow(sd_diag_df) > 0) {
  save_table(sd_diag_df, DIR_THIS, "sd_selection_diagnostics.csv")
}
if (nrow(resid_diag_df) > 0) {
  save_table(resid_diag_df, DIR_THIS, "residualize_diagnostics.csv")
  log_info("residualize 诊断:")
  print(resid_diag_df, row.names = FALSE)
}

# ============================================================
# ============================================================
log_info("--- 6. 跨对比 MCI 汇总 ---")

MCI_summary <- NULL
if (length(all_MCI) > 0) {
  MCI_summary <- data.frame()
  for (label in names(all_MCI)) {
    MCI_df <- extract_top_mci_per_state(all_MCI[[label]][["MCI"]])
    if (nrow(MCI_df) == 0) next
    MCI_df$Comparison <- label
    MCI_summary <- rbind(MCI_summary, MCI_df)
  }
  if (nrow(MCI_summary) > 0) {
    MCI_summary <- MCI_summary[order(MCI_summary$Comparison, -MCI_summary$MCI), ]
    save_table(MCI_summary, DIR_THIS, "MCI_summary_across.csv")

    p_cross <- ggplot(MCI_summary,
                      aes(x = Comparison, y = State, fill = MCI)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = round(MCI, 2)), size = 3.5, color = "black") +
      scale_fill_gradient(low = "white", high = "#c0392b", name = "MCI") +
      labs(title = "MCI scores across comparisons",
           subtitle = "Higher MCI = more likely tipping point",
           x = "Comparison", y = "State") +
      theme_bw(base_size = 12) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5, size = 9))
    save_fig(p_cross, DIR_THIS, "MCI_comparison_across", width = 7, height = 5)
  }
}

# ============================================================
# ============================================================
log_info("--- 7. poly3_F 交叉验证 ---")

CTS_poly3F_df <- NULL
poly3F_path <- file.path(DIR_RESULT, "06_1_TRADE_vs_adjacent", "tables",
                          "poly3F_gene_level.csv")

if (file.exists(poly3F_path) && length(all_CTS) > 0) {
  poly3F <- read.csv(poly3F_path, stringsAsFactors = FALSE)
  log_info("poly3F_gene_level: {nrow(poly3F)} 基因")

  bg_p3 <- mean(poly3F$p3_FDR < 0.05, na.rm = TRUE)
  bg_F  <- mean(poly3F$F_FDR < 0.05, na.rm = TRUE)
  log_info("背景 p3_FDR<0.05 比例: {round(bg_p3, 4)}")
  log_info("背景 F_FDR<0.05 比例: {round(bg_F, 4)}")

  CTS_poly3F_df <- data.frame()
  for (label in names(all_CTS)) {
    cts_genes <- all_CTS[[label]]
    if (length(cts_genes) < 5) next

    sub <- poly3F[poly3F$SYMBOL %in% cts_genes, ]
    n_total <- nrow(sub)
    if (n_total < 5) next

    n_p3 <- sum(sub$p3_FDR < 0.05, na.rm = TRUE)
    n_F  <- sum(sub$F_FDR < 0.05, na.rm = TRUE)

    p_hyper_p3 <- phyper(n_p3 - 1, ceiling(bg_p3 * nrow(poly3F)),
                          nrow(poly3F) - ceiling(bg_p3 * nrow(poly3F)),
                          n_total, lower.tail = FALSE)
    p_hyper_F <- phyper(n_F - 1, ceiling(bg_F * nrow(poly3F)),
                         nrow(poly3F) - ceiling(bg_F * nrow(poly3F)),
                         n_total, lower.tail = FALSE)

    CTS_poly3F_df <- rbind(CTS_poly3F_df, data.frame(
      Comparison       = label,
      CTS_n_genes      = length(cts_genes),
      CTS_n_in_poly3F  = n_total,
      n_p3_FDR05       = n_p3,
      n_F_FDR05        = n_F,
      pct_p3_FDR05     = round(100 * n_p3 / n_total, 2),
      pct_F_FDR05      = round(100 * n_F / n_total, 2),
      enrich_p3        = round((n_p3 / n_total) / bg_p3, 2),
      enrich_F         = round((n_F / n_total) / bg_F, 2),
      p_hyper_p3       = signif(p_hyper_p3, 3),
      p_hyper_F        = signif(p_hyper_F, 3),
      stringsAsFactors = FALSE))
  }

  if (nrow(CTS_poly3F_df) > 0) {
    CTS_poly3F_df$FDR_hyper_p3 <- round(p.adjust(CTS_poly3F_df$p_hyper_p3,
                                                   method = "BH"), 4)
    CTS_poly3F_df$FDR_hyper_F  <- round(p.adjust(CTS_poly3F_df$p_hyper_F,
                                                   method = "BH"), 4)
    save_table(CTS_poly3F_df, DIR_THIS, "CTS_poly3F_validation.csv")

    log_info("CTS poly3F 交叉验证:")
    print(CTS_poly3F_df, row.names = FALSE, digits = 4)

    p_enrich <- ggplot(CTS_poly3F_df,
                       aes(x = Comparison, y = enrich_p3,
                           fill = p_hyper_p3 < 0.05)) +
      geom_col(width = 0.6) +
      geom_hline(yintercept = 1, linetype = "dashed", color = "gray50") +
      geom_text(aes(label = sprintf("%.2f", enrich_p3)), vjust = -0.3, size = 3.5) +
      scale_fill_manual(values = c("TRUE" = "#c0392b", "FALSE" = "#B0B0B0"),
                         name = "p_hyper < 0.05") +
      labs(title = "CTS 基因的 p3_FDR<0.05 富集",
           x = "Comparison", y = "Enrichment fold (p3)") +
      theme_bw(base_size = 11) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    save_fig(p_enrich, DIR_THIS, "CTS_poly3F_enrichment", width = 6, height = 4)
  }
} else {
  log_warn("未找到 06_1 的 poly3F_gene_level.csv，跳过 poly3F 交叉验证")
}

# ============================================================
# ============================================================
log_info("--- 8. CTS 伪影检查 ---")

CTS_artifact_df <- NULL
if (length(all_CTS) > 0) {
  art_panels <- intersect(
    c("IEG", "IFN_I", "IFN_II", "ER_stress",
      "Oxidative_stress", "Inflammatory"),
    names(GENE_SETS))

  art_score_list <- list()
  for (pn in art_panels) {
    genes <- toupper(GENE_SETS[[pn]]$genes)
    genes_in <- intersect(genes, rownames(expr_mat))
    if (length(genes_in) < 5) next
    sub <- expr_mat[genes_in, , drop = FALSE]
    sub_z <- t(scale(t(sub)))
    art_score_list[[pn]] <- colMeans(sub_z, na.rm = TRUE)
  }

  if (length(art_score_list) > 0) {
    art_score_full <- as.data.frame(art_score_list)
    rownames(art_score_full) <- colnames(expr_mat)

    CTS_artifact_df <- data.frame()
    for (label in names(all_CTS)) {
      cts_genes <- all_CTS[[label]]
      cts_in <- intersect(cts_genes, rownames(expr_mat))
      if (length(cts_in) < 5) next

      lab_low  <- as.integer(sub("_vs_.*", "", label))
      lab_high <- as.integer(sub(".*_vs_", "", label))
      sub_samples <- rownames(clin)[clin$grade %in% c(lab_low, lab_high)]
      sub_samples <- intersect(sub_samples, colnames(expr_mat))
      if (length(sub_samples) < 5) next

      cts_expr <- colMeans(expr_mat[cts_in, sub_samples, drop = FALSE])
      art_sub  <- art_score_full[sub_samples, , drop = FALSE]

      for (pn in colnames(art_sub)) {
        if (sd(cts_expr, na.rm = TRUE) == 0 ||
            sd(art_sub[[pn]], na.rm = TRUE) == 0) next
        ct <- tryCatch(cor.test(cts_expr, art_sub[[pn]]),
                       error = function(e) NULL)
        if (is.null(ct)) next
        CTS_artifact_df <- rbind(CTS_artifact_df, data.frame(
          Comparison = label, Panel = pn,
          r = round(unname(ct$estimate), 3),
          p = ct$p.value,
          stringsAsFactors = FALSE))
      }

      sex_sub <- clin$sex_f[match(sub_samples, rownames(clin))]
      if (length(unique(sex_sub[!is.na(sex_sub)])) == 2) {
        tt <- tryCatch(t.test(cts_expr ~ sex_sub), error = function(e) NULL)
        if (!is.null(tt)) {
          CTS_artifact_df <- rbind(CTS_artifact_df, data.frame(
            Comparison = label, Panel = "Sex",
            r = NA_real_, p = tt$p.value,
            stringsAsFactors = FALSE))
        }
      }
    }

    if (nrow(CTS_artifact_df) > 0) {
      CTS_artifact_df$p_adj <- NA_real_
      for (lab in unique(CTS_artifact_df$Comparison)) {
        idx <- CTS_artifact_df$Comparison == lab
        CTS_artifact_df$p_adj[idx] <- p.adjust(CTS_artifact_df$p[idx],
                                                method = "BH")
      }
      CTS_artifact_df$flagged <- !is.na(CTS_artifact_df$r) &
                                 abs(CTS_artifact_df$r) > ARTIFACT_R &
                                 CTS_artifact_df$p_adj < P_THRESHOLD
      CTS_artifact_df$flagged[CTS_artifact_df$Panel == "Sex"] <-
        CTS_artifact_df$p_adj[CTS_artifact_df$Panel == "Sex"] < P_THRESHOLD

      save_table(CTS_artifact_df, DIR_THIS, "CTS_artifact_check.csv")
      log_info("CTS-伪影 flagged 数: {sum(CTS_artifact_df$flagged, na.rm=TRUE)}")
    }
  }
}

# ============================================================
# ============================================================
log_info("--- 9. IEG 中介分析 ---")

IEG_mediation_df <- NULL
if ("IEG" %in% names(GENE_SETS) && length(all_CTS) > 0) {
  ieg_genes <- toupper(GENE_SETS$IEG$genes)
  ieg_in <- intersect(ieg_genes, rownames(expr_mat))

  if (length(ieg_in) >= 5) {
    ieg_sub <- expr_mat[ieg_in, , drop = FALSE]
    ieg_z <- t(scale(t(ieg_sub)))
    ieg_score_full <- colMeans(ieg_z, na.rm = TRUE)

    IEG_mediation_df <- data.frame()
    for (label in names(all_CTS)) {
      cts_genes <- all_CTS[[label]]
      cts_in <- intersect(cts_genes, rownames(expr_mat))
      if (length(cts_in) < 5) next

      lab_low  <- as.integer(sub("_vs_.*", "", label))
      lab_high <- as.integer(sub(".*_vs_", "", label))
      sub_samples <- rownames(clin)[clin$grade %in% c(lab_low, lab_high)]
      sub_samples <- intersect(sub_samples, colnames(expr_mat))
      if (length(sub_samples) < 10) next

      cts_score <- colMeans(expr_mat[cts_in, sub_samples, drop = FALSE])
      ieg_score <- ieg_score_full[sub_samples]
      pmi_sub   <- clin$pmi_v[match(sub_samples, rownames(clin))]
      rin_sub   <- clin$rin_v[match(sub_samples, rownames(clin))]

      IEG_mediation_df <- rbind(IEG_mediation_df, data.frame(
        Comparison    = label,
        CTS_n         = length(cts_in),
        r_CTS_IEG     = round(cor(cts_score, ieg_score, use = "complete.obs"), 3),
        r_CTS_PMI     = round(cor(cts_score, pmi_sub,   use = "complete.obs"), 3),
        r_IEG_PMI     = round(cor(ieg_score, pmi_sub,   use = "complete.obs"), 3),
        r_CTS_RIN     = round(cor(cts_score, rin_sub,   use = "complete.obs"), 3),
        r_IEG_RIN     = round(cor(ieg_score, rin_sub,   use = "complete.obs"), 3),
        stringsAsFactors = FALSE))
    }

    if (nrow(IEG_mediation_df) > 0) {
      save_table(IEG_mediation_df, DIR_THIS, "IEG_mediation_analysis.csv")
    }
  }
}

# ============================================================
# 10. summary
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("=== 06_2 BioTIP 汇总 (v2) ===\n")
cat(strrep("=", 70), "\n", sep = "")

cat("\n【Q0】residualize 诊断\n")
if (nrow(resid_diag_df) > 0) {
  print(resid_diag_df, row.names = FALSE, digits = 4)
}

cat("\n【Q1】MCI 汇总\n")
if (!is.null(MCI_summary) && nrow(MCI_summary) > 0) {
  print(MCI_summary, row.names = FALSE, digits = 4)
}

cat("\n【Q2】CTS 基因数\n")
if (length(all_CTS) > 0) {
  for (lab in names(all_CTS)) {
    cat(sprintf("  %-10s  %d\n", lab, length(all_CTS[[lab]])))
  }
}

cat("\n【Q3】CTS poly3F 交叉验证\n")
if (!is.null(CTS_poly3F_df) && nrow(CTS_poly3F_df) > 0) {
  print(CTS_poly3F_df, row.names = FALSE, digits = 4)
}

cat("\n【Q4】CTS 伪影 flagged\n")
if (!is.null(CTS_artifact_df)) {
  cat(sprintf("  %d / %d\n",
              sum(CTS_artifact_df$flagged, na.rm = TRUE),
              nrow(CTS_artifact_df)))
}

cat("\n【Q5】sd_selection 诊断\n")
if (nrow(sd_diag_df) > 0) {
  print(sd_diag_df, row.names = FALSE)
}

cat("\n--- 输出位置 ---\n")
cat("目录:", DIR_THIS, "\n")

# ============================================================
# 11. save RDS
# ============================================================
save_rds(list(
  clinical         = clin,
  expr_mat         = expr_mat,
  all_MCI          = all_MCI,
  all_CTS          = all_CTS,
  MCI_summary      = MCI_summary,
  diag_sub_df      = diag_sub_df,
  sd_diag_df       = sd_diag_df,
  resid_diag_df    = resid_diag_df,
  CTS_artifact_df  = CTS_artifact_df,
  IEG_mediation_df = IEG_mediation_df,
  CTS_poly3F_df    = CTS_poly3F_df,
  params = list(
    sd_cutoff         = SD_CUTOFF,
    min_module_size   = MIN_MODULE_SIZE,
    network_fdr       = NETWORK_FDR,
    residualize       = RESIDUALIZE_BEFORE_BIOTIP,
    comparisons       = c("2_vs_1", "3_vs_2", "4_vs_3"),
    input_expr        = "01_1b_batch_enhance/expr_wgcna.rds",
    version           = "v2"
  )
), DIR_THIS, "mw06_2_biotip_results.rds")

log_info("=== 06_2_BioTIP.r (v2) done ===")

