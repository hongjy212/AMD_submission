# ============================================================
# 06_3_cellchat.r
#
# Output：results/06_3_cellchat/
#   cellchat.rds
#   rpe_subset_matrix.rds
#   tables/*.csv
#   figures/*.pdf
# ============================================================

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(CellChat); library(Matrix); library(SingleCellExperiment)
  library(dplyr); library(tibble); library(tidyr); library(ggplot2)
  library(pheatmap); library(circlize)
})

DIR_THIS <- setup_output_dir("06_3_cellchat")
log_info("=== 06_3_cellchat start ===")

DIR_SC  <- "E:/retina_tools/raw/GSE221042"
N_PERM  <- 1000
N_BOOT  <- 1000
SEED    <- 42

my_colors <- c(
  AC   = "#E41A1C", BC   = "#377EB8", Cone = "#4DAF4A",
  HC   = "#984EA3", MG   = "#FF7F00", RGC  = "#FFFF33",
  Rod  = "#A65628", RPE  = "#F781BF"
)

# ============================================================
# ============================================================

# ------------------------------------------------------------
# A1. read
# ------------------------------------------------------------
log_info("--- A1. 读入 ---")

mtx      <- readMM(file.path(DIR_SC, "GSE221042_data.mtx.gz"))
genes    <- read.csv(file.path(DIR_SC, "GSE221042_genes.csv.gz"),
                     header = FALSE, stringsAsFactors = FALSE)
meta_raw <- read.csv(file.path(DIR_SC, "GSE221042_sample_annotations.csv.gz"),
                     check.names = FALSE, stringsAsFactors = FALSE)

log_info(sprintf("矩阵维度: %d × %d", nrow(mtx), ncol(mtx)))

# ------------------------------------------------------------
# ------------------------------------------------------------
if (nrow(mtx) == nrow(meta_raw)) {
  log_info("矩阵方向: cells × genes")
  mtx_cells_genes <- mtx
} else if (ncol(mtx) == nrow(meta_raw)) {
  log_info("矩阵方向: genes × cells（转置）")
  mtx_cells_genes <- t(mtx)
} else {
  stop("矩阵维度与 meta_raw 行数不匹配")
}

gene_names <- toupper(trimws(sub(" \\(.*", "", as.character(genes$V2))))
colnames(mtx_cells_genes) <- make.unique(gene_names)
rownames(mtx_cells_genes) <- make.unique(as.character(meta_raw[[1]]))

log_info(sprintf("QC 前矩阵: %d cells × %d genes",
                 nrow(mtx_cells_genes), ncol(mtx_cells_genes)))

# ------------------------------------------------------------
# A3. QC
# ------------------------------------------------------------
log_info("--- A2. QC（200–7,500 基因，< 10% 线粒体） ---")

nFeature <- Matrix::rowSums(mtx_cells_genes > 0)
mito_genes <- grep("^MT-", colnames(mtx_cells_genes), ignore.case = TRUE)
if (length(mito_genes) > 0) {
  percent_mt <- Matrix::rowSums(mtx_cells_genes[, mito_genes, drop = FALSE]) /
                Matrix::rowSums(mtx_cells_genes) * 100
} else {
  percent_mt <- rep(0, nrow(mtx_cells_genes))
}

keep <- nFeature > 200 & nFeature < 7500 & percent_mt < 10
log_info(sprintf("QC 前细胞数: %d", nrow(mtx_cells_genes)))
log_info(sprintf("QC 后细胞数: %d", sum(keep)))

mtx_cells_genes <- mtx_cells_genes[keep, ]
meta_raw <- meta_raw[keep, ]
rownames(mtx_cells_genes) <- make.unique(as.character(meta_raw[[1]]))

log_info(sprintf("QC 后矩阵: %d cells × %d genes",
                 nrow(mtx_cells_genes), ncol(mtx_cells_genes)))

# ------------------------------------------------------------
# ------------------------------------------------------------
log_info("--- A3. 细胞类型注释 ---")

MARKERS <- list(
  Rod  = c("RHO","GNAT1","PDE6B","NRL"),
  Cone = c("ARR3","GNAT2","PDE6H"),
  HC   = c("ONECUT1","ONECUT2","LHX1","TFAP2B"),
  BC   = c("VSX1","VSX2","GRM6"),
  AC   = c("GAD1","GAD2"),
  MG   = c("RLBP1","GLUL"),
  RGC  = c("RBPMS","NEFL","NEFM"),
  RPE  = c("RPE65","BEST1","TTR"))

gene_ids <- colnames(mtx_cells_genes)
libsize_all <- Matrix::rowSums(mtx_cells_genes)
norm_all    <- log1p(mtx_cells_genes / libsize_all * 1e4)

scores <- do.call(cbind, lapply(MARKERS, function(g) {
  r <- which(gene_ids %in% g)
  if (length(r) == 0) return(rep(0, nrow(norm_all)))
  Matrix::rowMeans(norm_all[, r, drop = FALSE])
}))
colnames(scores) <- names(MARKERS)
cell_type <- colnames(scores)[apply(scores, 1, which.max)]

log_info("细胞类型分布:")
print(table(cell_type))

if (length(unique(cell_type)) < 2) {
  stop("只注释到 1 种细胞类型，无法做 CellChat。")
}

# ------------------------------------------------------------
# A5. build CellChat object
# ------------------------------------------------------------
log_info("--- A4. 构建 CellChat 对象 ---")

mtx_genes_cells <- t(mtx_cells_genes)
mtx_sparse <- as(mtx_genes_cells, "CsparseMatrix")
rownames(mtx_sparse) <- colnames(mtx_cells_genes)
colnames(mtx_sparse) <- rownames(mtx_cells_genes)

cell_id <- colnames(mtx_sparse)
stopifnot(length(cell_id) == length(cell_type))

cellchat <- createCellChat(
  object = mtx_sparse,
  meta = data.frame(
    labels = cell_type,
    samples = factor(meta_raw$sample_labels),
    row.names = cell_id
  ),
  group.by = "labels"
)

cellchat@DB <- CellChatDB.human

rpe_cells_all <- cell_id[cell_type == "RPE"]
if (length(rpe_cells_all) > 0) {
  rpe_subset_matrix <- mtx_sparse[, rpe_cells_all, drop = FALSE]
  saveRDS(rpe_subset_matrix, file.path(DIR_THIS, "rpe_subset_matrix.rds"))
  log_info(sprintf("RPE 子集矩阵已保存: %d genes × %d cells",
                   nrow(rpe_subset_matrix), ncol(rpe_subset_matrix)))
}

# ------------------------------------------------------------
# ------------------------------------------------------------
log_info("--- A5. 预处理 ---")

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

log_info(sprintf("L-R pairs: %d", nrow(cellchat@LR$LRsig)))
if (nrow(cellchat@LR$LRsig) == 0) stop("L-R pairs = 0")

# ------------------------------------------------------------
# ------------------------------------------------------------
log_info("--- A6. 计算通讯概率（约 25 分钟） ---")

cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
log_info(sprintf("通路数: %d", length(cellchat@netP$pathways)))

cellchat <- aggregateNet(cellchat)
log_info(sprintf("net$count: %s", paste(dim(cellchat@net$count), collapse = " × ")))

# ------------------------------------------------------------
# A8. extractResults
# ------------------------------------------------------------
log_info("--- A7. 提取通讯结果 ---")

comm_count <- subsetCommunication(cellchat)
write.csv(comm_count, file.path(DIR_THIS, "tables", "communication_count.csv"),
          row.names = FALSE)

comm_prob <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(comm_prob, file.path(DIR_THIS, "tables", "communication_prob.csv"),
          row.names = FALSE)

pathway_prob <- cellchat@netP$pathways
write.csv(data.frame(pathway = pathway_prob),
          file.path(DIR_THIS, "tables", "pathway_prob.csv"), row.names = FALSE)

top_axis <- comm_count %>%
  group_by(source, target) %>%
  summarise(n_interactions = n(), .groups = "drop") %>%
  arrange(desc(n_interactions))
write.csv(top_axis, file.path(DIR_THIS, "tables", "top_axis.csv"), row.names = FALSE)

saveRDS(cellchat, file.path(DIR_THIS, "cellchat.rds"))

# ------------------------------------------------------------
# ------------------------------------------------------------
log_info("--- A8. 主图 ---")

stopifnot(all(levels(cellchat@idents) %in% names(my_colors)))

mat <- as.matrix(cellchat@net$count)
storage.mode(mat) <- "numeric"
mat <- mat[names(my_colors), names(my_colors)]

pdf(file.path(DIR_THIS, "figures", "Fig1_circle.pdf"), width = 8, height = 8)
circos.clear()
circos.par(gap.after = 2)
chordDiagram(as.data.frame(as.table(mat)),
             grid.col = my_colors,
             transparency = 0.3,
             annotationTrack = c("name", "grid"))
title("Number of interactions")
dev.off()

pdf(file.path(DIR_THIS, "figures", "Fig2_heatmap.pdf"), width = 8, height = 8)
pheatmap(mat,
         color = colorRampPalette(c("white", "red"))(100),
         cluster_rows = FALSE, cluster_cols = FALSE,
         display_numbers = TRUE, number_format = "%.0f",
         fontsize_number = 8,
         main = "Number of interactions")
dev.off()

pdf(file.path(DIR_THIS, "figures", "Fig3_top_pathway_network.pdf"),
    width = 10, height = 8)
netVisual_aggregate(cellchat, signaling = pathway_prob[1],
                    layout = "circle", color.use = my_colors)
dev.off()

pdf(file.path(DIR_THIS, "figures", "Fig4_top_pathways_bubble.pdf"),
    width = 10, height = 8)
rankNet(cellchat, mode = "single", stacked = TRUE,
        color.use = c("grey", "firebrick"))
dev.off()

# ============================================================
# ============================================================

# ------------------------------------------------------------
# ------------------------------------------------------------
log_info("--- B1. 通路重要性排名 ---")

netP_prob <- cellchat@netP$prob
pathways_all <- dimnames(netP_prob)[[3]]
cell_types_all <- dimnames(netP_prob)[[1]]

role_rows <- list()
for (pw in pathways_all) {
  m <- netP_prob[, , pw]
  if (is.null(m) || !is.matrix(m)) next
  for (ct in cell_types_all) {
    outdeg <- sum(m[ct, ])
    indeg  <- sum(m[, ct])
    role_rows[[paste(pw, ct)]] <- data.frame(
      pathway = pw, cell_type = ct,
      outdeg = outdeg, indeg = indeg,
      stringsAsFactors = FALSE
    )
  }
}
role_all <- bind_rows(role_rows)

if (nrow(role_all) > 0) {
  write.csv(role_all, file.path(DIR_THIS, "tables", "signaling_role_all.csv"),
            row.names = FALSE)

  rpe_sender <- role_all %>%
    filter(cell_type == "RPE") %>%
    arrange(desc(outdeg))
  write.csv(rpe_sender, file.path(DIR_THIS, "tables", "RPE_sender_top.csv"),
            row.names = FALSE)

  cat("\n=== RPE 作为 sender 的 top 10 通路 ===\n")
  print(head(rpe_sender, 10), row.names = FALSE, digits = 3)

  top_sender <- role_all %>%
    group_by(pathway) %>%
    slice_max(outdeg, n = 1) %>%
    ungroup() %>%
    arrange(desc(outdeg))
  write.csv(top_sender, file.path(DIR_THIS, "tables", "top_sender_per_pathway.csv"),
            row.names = FALSE)

  cat("\n=== 每条通路的最强 sender（top 10） ===\n")
  print(head(top_sender, 10), row.names = FALSE, digits = 3)
}

# ------------------------------------------------------------
# ------------------------------------------------------------
log_info("--- B2. RPE 相关通讯对筛选 ---")

rpe_comm <- comm_count %>%
  filter(source == "RPE" | target == "RPE") %>%
  arrange(desc(prob))
write.csv(rpe_comm, file.path(DIR_THIS, "tables", "RPE_related_communications.csv"),
          row.names = FALSE)
log_info(sprintf("RPE 相关通讯对: %d", nrow(rpe_comm)))

rpe_as_source <- comm_count %>% filter(source == "RPE")
rpe_as_target <- comm_count %>% filter(target == "RPE")
cat(sprintf("RPE 作为 source: %d 对\n", nrow(rpe_as_source)))
cat(sprintf("RPE 作为 target: %d 对\n", nrow(rpe_as_target)))

if (nrow(rpe_comm) > 0) {
  rpe_pathway_count <- rpe_comm %>%
    group_by(pathway_name) %>%
    summarise(n = n(), .groups = "drop") %>%
    arrange(desc(n))
  write.csv(rpe_pathway_count,
            file.path(DIR_THIS, "tables", "RPE_pathway_count.csv"),
            row.names = FALSE)

  cat("\n=== RPE 相关通讯对富集的 top 15 通路 ===\n")
  print(head(rpe_pathway_count, 15), row.names = FALSE)

  top_n <- head(rpe_pathway_count, 20)
  top_n$pathway_name <- factor(top_n$pathway_name,
                               levels = top_n$pathway_name[order(top_n$n)])
  p <- ggplot(top_n, aes(x = pathway_name, y = n, fill = n)) +
    geom_col() + coord_flip() + theme_bw(base_size = 11) +
    scale_fill_gradient(low = "#FDDBC7", high = "#B2182B") +
    labs(x = "", y = "Number of RPE-related communications",
         title = "Pathways enriched in RPE-related communications") +
    theme(legend.position = "none")
  save_fig(p, DIR_THIS, "RPE_pathway_bar", width = 8, height = 6)
}

# ============================================================
# ============================================================

# ------------------------------------------------------------
# ------------------------------------------------------------
log_info("--- C1. RPE vs Rod/Cone 分组 ---")

rodcone_comm <- comm_count %>%
  filter((source %in% c("Rod", "Cone") | target %in% c("Rod", "Cone")) &
         source != "RPE" & target != "RPE")

log_info(sprintf("RPE 组通讯对: %d", nrow(rpe_comm)))
log_info(sprintf("Rod/Cone 组通讯对: %d", nrow(rodcone_comm)))

all_pathways <- sort(unique(comm_count$pathway_name))
all_pathways <- all_pathways[!is.na(all_pathways) & all_pathways != ""]

compute_log2_ratio <- function(rpe_vec, rc_vec, all_pathways) {
  n_rpe <- table(factor(rpe_vec, levels = all_pathways))
  n_rc  <- table(factor(rc_vec, levels = all_pathways))
  prop_rpe <- as.numeric(n_rpe) / sum(n_rpe)
  prop_rc  <- as.numeric(n_rc) / sum(n_rc)
  log2((prop_rpe + 1e-6) / (prop_rc + 1e-6))
}

observed_log2 <- compute_log2_ratio(rpe_comm$pathway_name,
                                     rodcone_comm$pathway_name,
                                     all_pathways)

observed_df <- data.frame(
  pathway = all_pathways,
  n_RPE = as.numeric(table(factor(rpe_comm$pathway_name, levels = all_pathways))),
  n_RodCone = as.numeric(table(factor(rodcone_comm$pathway_name, levels = all_pathways))),
  log2_ratio = observed_log2,
  stringsAsFactors = FALSE
)

write.csv(observed_df,
          file.path(DIR_THIS, "tables", "RPE_vs_RodCone_observed.csv"),
          row.names = FALSE)

# ------------------------------------------------------------
# ------------------------------------------------------------
log_info("--- C2. 置换检验（1000 次） ---")

set.seed(SEED)

all_comm_pathways <- c(rpe_comm$pathway_name, rodcone_comm$pathway_name)
n_rpe_total <- length(rpe_comm$pathway_name)

perm_matrix <- matrix(NA, nrow = N_PERM, ncol = length(all_pathways))
colnames(perm_matrix) <- all_pathways

for (i in seq_len(N_PERM)) {
  if (i %% 200 == 0) log_info(sprintf("  置换进度: %d / %d", i, N_PERM))
  shuffled <- sample(all_comm_pathways, replace = FALSE)
  fake_rpe <- shuffled[1:n_rpe_total]
  fake_rc  <- shuffled[(n_rpe_total + 1):length(shuffled)]
  perm_matrix[i, ] <- compute_log2_ratio(fake_rpe, fake_rc, all_pathways)
}

emp_p <- sapply(seq_along(all_pathways), function(j) {
  obs <- observed_log2[j]
  mean(abs(perm_matrix[, j]) >= abs(obs))
})
emp_p <- pmax(emp_p, 1 / N_PERM)
fdr <- p.adjust(emp_p, method = "BH")

null_ci_low  <- apply(perm_matrix, 2, quantile, probs = 0.025, na.rm = TRUE)
null_ci_high <- apply(perm_matrix, 2, quantile, probs = 0.975, na.rm = TRUE)

perm_df <- data.frame(
  pathway = all_pathways,
  n_RPE = observed_df$n_RPE,
  n_RodCone = observed_df$n_RodCone,
  log2_ratio = observed_log2,
  emp_p = emp_p,
  FDR = fdr,
  null_ci_low = null_ci_low,
  null_ci_high = null_ci_high,
  stringsAsFactors = FALSE
) %>%
  arrange(FDR, desc(abs(log2_ratio)))

write.csv(perm_df,
          file.path(DIR_THIS, "tables", "RPE_vs_RodCone_permutation.csv"),
          row.names = FALSE)

cat("\n=== 置换检验 + FDR 校正后 top 20 ===\n")
print(head(perm_df %>%
             select(pathway, n_RPE, n_RodCone, log2_ratio, emp_p, FDR) %>%
             arrange(FDR, desc(abs(log2_ratio))), 20),
      row.names = FALSE, digits = 4)

n_sig_005 <- sum(perm_df$FDR < 0.05, na.rm = TRUE)
n_sig_010 <- sum(perm_df$FDR < 0.10, na.rm = TRUE)
log_info(sprintf("FDR < 0.05 通路数: %d", n_sig_005))
log_info(sprintf("FDR < 0.10 通路数: %d", n_sig_010))

# ------------------------------------------------------------
# C3. Bootstrap 95% CI
# ------------------------------------------------------------
log_info("--- C3. Bootstrap 95% CI（1000 次） ---")

set.seed(SEED)

top_enriched <- observed_df %>% arrange(desc(log2_ratio)) %>% head(15) %>% pull(pathway)
top_depleted <- observed_df %>% arrange(log2_ratio) %>% head(15) %>% pull(pathway)
top_pathways <- unique(c(top_enriched, top_depleted))
top_pathways <- top_pathways[!is.na(top_pathways) & top_pathways != ""]

bootstrap_ci <- function(rpe_vec, rc_vec, target_pathway, n_boot = N_BOOT) {
  log2_boot <- numeric(n_boot)
  for (b in seq_len(n_boot)) {
    rpe_boot <- sample(rpe_vec, replace = TRUE)
    rc_boot  <- sample(rc_vec, replace = TRUE)
    n_rpe_target <- sum(rpe_boot == target_pathway, na.rm = TRUE)
    n_rc_target  <- sum(rc_boot == target_pathway, na.rm = TRUE)
    prop_rpe <- n_rpe_target / length(rpe_boot)
    prop_rc  <- n_rc_target / length(rc_boot)
    log2_boot[b] <- log2((prop_rpe + 1e-6) / (prop_rc + 1e-6))
  }
  ci <- quantile(log2_boot, probs = c(0.025, 0.975), na.rm = TRUE)
  c(ci_low = as.numeric(ci[1]),
    ci_high = as.numeric(ci[2]),
    boot_sd = as.numeric(sd(log2_boot)))
}

boot_results <- list()
for (pw in top_pathways) {
  boot_results[[pw]] <- bootstrap_ci(rpe_comm$pathway_name,
                                      rodcone_comm$pathway_name, pw)
}

boot_list <- lapply(names(boot_results), function(pw) {
  idx <- match(pw, observed_df$pathway)
  r <- boot_results[[pw]]
  list(
    pathway = pw,
    log2_ratio = as.numeric(observed_df$log2_ratio[idx]),
    boot_ci_low = unname(as.numeric(r["ci_low"])),
    boot_ci_high = unname(as.numeric(r["ci_high"])),
    boot_sd = unname(as.numeric(r["boot_sd"]))
  )
})

boot_df <- bind_rows(lapply(boot_list, as.data.frame,
                             stringsAsFactors = FALSE)) %>%
  left_join(perm_df %>% select(pathway, n_RPE, n_RodCone, emp_p, FDR),
            by = "pathway") %>%
  arrange(desc(log2_ratio))

cat(sprintf("\nboot_df: %d 行\n", nrow(boot_df)))
cat(sprintf("  boot_ci_low  NA 数量: %d\n", sum(is.na(boot_df$boot_ci_low))))
cat(sprintf("  boot_ci_high NA 数量: %d\n", sum(is.na(boot_df$boot_ci_high))))

write.csv(boot_df,
          file.path(DIR_THIS, "tables", "RPE_vs_RodCone_bootstrap_CI.csv"),
          row.names = FALSE)

cat("\n=== Bootstrap 95% CI（top 富集通路） ===\n")
print(head(boot_df %>% arrange(desc(log2_ratio)), 15),
      row.names = FALSE, digits = 3)

cat("\n=== Bootstrap 95% CI（top 缺失通路） ===\n")
print(head(boot_df %>% arrange(log2_ratio), 15),
      row.names = FALSE, digits = 3)

# ------------------------------------------------------------
# ------------------------------------------------------------
final_df <- perm_df %>%
  left_join(boot_df %>% select(pathway, boot_ci_low, boot_ci_high, boot_sd),
            by = "pathway") %>%
  arrange(FDR, desc(abs(log2_ratio)))

write.csv(final_df,
          file.path(DIR_THIS, "tables", "RPE_vs_RodCone_final.csv"),
          row.names = FALSE)

cat(sprintf("\nfinal_df: %d 行, 有 CI 的通路: %d\n",
            nrow(final_df), sum(!is.na(final_df$boot_ci_low))))

# ============================================================
# ============================================================

# ------------------------------------------------------------
# D1. Volcano plot
# ------------------------------------------------------------
log_info("--- D1. 统计可视化 ---")

plot_df <- perm_df %>%
  mutate(
    sig = case_when(
      FDR < 0.05 ~ "FDR < 0.05",
      FDR < 0.10 ~ "FDR < 0.10",
      TRUE ~ "NS"
    ),
    neglog10fdr = -log10(pmax(FDR, 1e-10))
  )

p_volcano <- ggplot(plot_df, aes(x = log2_ratio, y = neglog10fdr, color = sig)) +
  geom_point(aes(size = n_RPE + n_RodCone), alpha = 0.7) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
  scale_color_manual(values = c("FDR < 0.05" = "#C44E52",
                                 "FDR < 0.10" = "#F0AD4E",
                                 "NS" = "grey70")) +
  scale_size_continuous(name = "Total pairs", range = c(2, 10)) +
  theme_bw(base_size = 11) +
  labs(x = "log2(RPE / RodCone)", y = "-log10(FDR)",
       title = "RPE vs Rod/Cone pathway usage (permutation test)",
       color = "")
save_fig(p_volcano, DIR_THIS, "RPE_vs_RodCone_volcano", width = 8, height = 6)

# ------------------------------------------------------------
# D2. Forest plot
# ------------------------------------------------------------
forest_df <- final_df %>%
  filter(!is.na(boot_ci_low) & !is.na(boot_ci_high)) %>%
  arrange(desc(log2_ratio)) %>%
  head(20)

if (nrow(forest_df) > 0) {
  forest_df$pathway <- factor(forest_df$pathway,
                               levels = forest_df$pathway[order(forest_df$log2_ratio)])
  p_forest <- ggplot(forest_df, aes(x = log2_ratio, y = pathway)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = boot_ci_low, xmax = boot_ci_high),
                   height = 0.2, color = "grey40") +
    geom_point(aes(color = FDR < 0.05), size = 3) +
    scale_color_manual(values = c("TRUE" = "#C44E52", "FALSE" = "grey60"),
                       name = "FDR < 0.05") +
    theme_bw(base_size = 11) +
    labs(x = "log2(RPE / RodCone)", y = "",
         title = "RPE vs Rod/Cone: top pathways (bootstrap 95% CI)")
  save_fig(p_forest, DIR_THIS, "RPE_vs_RodCone_CI_forest", width = 8, height = 7)
}

# ------------------------------------------------------------
# D3. RPE marker QC
# ------------------------------------------------------------
log_info("--- D3. RPE 细胞质量检查 ---")

rpe_subset_file <- file.path(DIR_THIS, "rpe_subset_matrix.rds")
if (file.exists(rpe_subset_file)) {
  mtx_rpe <- readRDS(rpe_subset_file)
  cat(sprintf("RPE 子集矩阵: %d genes × %d cells\n",
              nrow(mtx_rpe), ncol(mtx_rpe)))

  check_markers <- list(
    RPE_marker  = c("RPE65", "BEST1", "TTR"),
    Rod_marker  = c("RHO", "GNAT1", "PDE6B", "NRL"),
    Cone_marker = c("ARR3", "GNAT2", "PDE6H"),
    MG_marker   = c("RLBP1", "GLUL"),
    RGC_marker  = c("RBPMS", "NEFL", "NEFM"),
    Immune      = c("PTPRC", "CD68", "CD14"),
    Endo        = c("PECAM1", "VWF", "CLDN5")
  )

  qc_rows <- list()
  for (ct in names(check_markers)) {
    for (g in check_markers[[ct]]) {
      if (!g %in% rownames(mtx_rpe)) next
      expr <- mtx_rpe[g, ]
      qc_rows[[paste(ct, g)]] <- data.frame(
        category = ct, gene = g,
        n_expr = sum(expr > 0),
        pct_expr = 100 * sum(expr > 0) / length(expr),
        mean_expr = mean(expr),
        stringsAsFactors = FALSE
      )
    }
  }
  qc_table <- bind_rows(qc_rows)
  write.csv(qc_table, file.path(DIR_THIS, "tables", "RPE_marker_QC.csv"),
            row.names = FALSE)

  cat("\n=== RPE 细胞 marker 表达检查 ===\n")
  print(qc_table, row.names = FALSE, digits = 3)

  cat("\n=== 按类别汇总 ===\n")
  qc_summary <- qc_table %>%
    group_by(category) %>%
    summarise(mean_pct = mean(pct_expr),
              max_pct = max(pct_expr),
              .groups = "drop") %>%
    arrange(desc(mean_pct))
  print(qc_summary, row.names = FALSE, digits = 2)

  if ("RPE_marker" %in% qc_summary$category &&
      "Rod_marker" %in% qc_summary$category) {
    rpe_pct <- qc_summary$mean_pct[qc_summary$category == "RPE_marker"]
    rod_pct <- qc_summary$mean_pct[qc_summary$category == "Rod_marker"]
    cat(sprintf("\nRPE marker 平均表达率: %.1f%%\n", rpe_pct))
    cat(sprintf("Rod marker 平均表达率: %.1f%%\n", rod_pct))
    if (rod_pct > rpe_pct) {
      cat("⚠ RPE 细胞中 Rod marker 表达率高于 RPE marker，可能存在污染\n")
    } else {
      cat("✅ RPE 细胞 marker 表达正常，注释可信\n")
    }
  }
}

# ------------------------------------------------------------
# D4. OXPHOS correlationpathwaycheck
# ------------------------------------------------------------
log_info("--- D4. OXPHOS 相关通路检查 ---")

oxphos_keywords <- c("OXPHOS", "mTOR", "AMPK", "VEGF", "MHC", "CD39",
                     "CD73", " adenosine", "TGFb", "BMP", "FGF",
                     "IGF", "WNT", "NOTCH", "HEDGEHOG", "PPAR",
                     "lipid", "metabol", "mitochond")

all_pathways_comm <- unique(comm_count$pathway_name)
oxphos_hit <- all_pathways_comm[grepl(paste(oxphos_keywords, collapse = "|"),
                                       all_pathways_comm, ignore.case = TRUE)]

cat(sprintf("\n数据库共 %d 条通路，匹配代谢关键词 %d 条\n",
            length(all_pathways_comm), length(oxphos_hit)))
if (length(oxphos_hit) > 0) print(oxphos_hit)

if (length(oxphos_hit) > 0) {
  oxphos_compare <- final_df %>%
    filter(pathway %in% oxphos_hit) %>%
    arrange(desc(log2_ratio))
  write.csv(oxphos_compare,
            file.path(DIR_THIS, "tables", "OXPHOS_pathway_RPE_vs_RodCone.csv"),
            row.names = FALSE)

  cat("\n=== 代谢/线粒体相关通路在 RPE vs Rod/Cone 里的对比 ===\n")
  print(oxphos_compare %>%
          select(pathway, n_RPE, n_RodCone, log2_ratio, emp_p, FDR,
                 boot_ci_low, boot_ci_high),
        row.names = FALSE, digits = 3)
}

# ============================================================
# summary
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("=== 06_3_cellchat 总结 ===\n")
cat(strrep("=", 70), "\n", sep = "")
cat(sprintf("QC 后细胞数: %d\n", ncol(mtx_sparse)))
cat(sprintf("基因数: %d\n", nrow(mtx_sparse)))
cat(sprintf("细胞类型: %s\n", paste(unique(cell_type), collapse = ", ")))
cat(sprintf("通讯对数: %d\n", nrow(comm_count)))
cat(sprintf("通路数: %d\n", length(pathway_prob)))
cat(sprintf("RPE 相关通讯对: %d\n", nrow(rpe_comm)))
cat(sprintf("RPE 作为 source: %d\n", nrow(rpe_as_source)))
cat(sprintf("RPE 作为 target: %d\n", nrow(rpe_as_target)))
cat(sprintf("RPE 细胞数: %d\n", length(rpe_cells_all)))
cat(sprintf("置换检验次数: %d\n", N_PERM))
cat(sprintf("Bootstrap 次数: %d\n", N_BOOT))
cat(sprintf("FDR < 0.05 通路数: %d\n", n_sig_005))
cat(sprintf("FDR < 0.10 通路数: %d\n", n_sig_010))
cat(sprintf("有 CI 的通路数: %d\n", sum(!is.na(final_df$boot_ci_low))))
cat(sprintf("输出目录: %s\n", DIR_THIS))

log_info("=== 06_3_cellchat done ===")

