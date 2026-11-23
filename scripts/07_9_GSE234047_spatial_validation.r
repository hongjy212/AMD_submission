# ============================================================
# 07_9_GSE234047_spatial_validation.r (v3.6)
# GSE234047 spatialtranscriptomevalidation NDUFB3 / Complex III
# Outputdirectory：results/07_9_GSE234047_spatial_validation/
#
#     （C6orf120 vs C6ORF120）
# ============================================================

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(data.table)
  library(dplyr); library(tidyr); library(ggplot2)
  library(patchwork); library(EBImage); library(R.utils)
})

DIR_THIS <- setup_output_dir("07_9_GSE234047_spatial_validation")
log_info("=== 07_9 GSE234047 spatial validation v3.6 start ===")

# ============================================================
# 0. Parameters
# ============================================================
DIR_SC     <- "E:/retina_tools/raw/GSE234047"
GENE_FOCUS <- "NDUFB3"
SEED       <- 42
MIN_SPOTS  <- 0

META_COLS <- c(
  "barcode","prediction_celltype",
  "Bipolar","Cone","Endothelial","Fibroblast","Immune","Interneuron",
  "Melanocyte","Muller.Astrocyte","Pericyte.SMC","RGC","Rod","RPE.x","Schwann",
  "res_ss","region","tissue","percent_CNV","image")

CELLTYPE_COLS <- c(
  "Bipolar","Cone","Endothelial","Fibroblast","Immune","Interneuron",
  "Melanocyte","Muller.Astrocyte","Pericyte.SMC","RGC","Rod","RPE.x","Schwann")

CIII_GENES <- c("UQCRQ","UQCRFS1","UQCRH","UQCRC1","UQCRC2","UQCRB","UQCR10")

CI_GENES <- c(
  "NDUFA1","NDUFA2","NDUFA3","NDUFA4","NDUFA5","NDUFA6","NDUFA7","NDUFA8",
  "NDUFA9","NDUFA10","NDUFA11","NDUFA12","NDUFA13",
  "NDUFB1","NDUFB2","NDUFB3","NDUFB4","NDUFB5","NDUFB6","NDUFB7","NDUFB8",
  "NDUFB9","NDUFB10","NDUFB11",
  "NDUFC1","NDUFC2",
  "NDUFS1","NDUFS2","NDUFS3","NDUFS4","NDUFS5","NDUFS6","NDUFS7","NDUFS8",
  "NDUFV1","NDUFV2","NDUFV3")

GENES_16 <- c("FCF1","UGGT1","CCDC3","S100B","BPNT1","BRI3BP","BTBD19",
              "DPH7","AKR1B1","ZNF623","ZNF587","RSPH3","TAF3","NDUFB8",
              "ZNF417","C6orf120")

normalize_tag <- function(x) {
  x <- tolower(x)
  x <- sub("^gsm[0-9]+[_-]?", "", x)
  x <- sub("\\.tif$", "", x)
  x <- sub("_normalized\\.csv\\.gz$", "", x)
  x <- gsub("[[:space:]\\-]", "_", x)
  x
}

set_assay_data_compat <- function(seu, data) {
  if ("layer" %in% names(formals(SetAssayData))) {
    SetAssayData(seu, layer = "data", new.data = data)
  } else {
    SetAssayData(seu, slot  = "data", new.data = data)
  }
}

# ============================================================
# 1. design
# ============================================================
log_info("--- 1. 读取 design ---")

design <- read.delim(file.path(DIR_SC, "GSE234047_design.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
design$group     <- ifelse(grepl("control", design$sample, ignore.case = TRUE),
                           "Control", "CNV")
design$image_tag <- normalize_tag(design$image)
print(design)

# ============================================================
# 2. expressionmatrix
# ============================================================
log_info("--- 2. 读入表达矩阵 ---")

norm_files <- list.files(DIR_SC, pattern = "_normalized\\.csv\\.gz$",
                         recursive = TRUE, full.names = TRUE)
log_info(sprintf("找到 %d 个 normalized 文件", length(norm_files)))

seu_list         <- list()
sample_meta_list <- list()

for (f in norm_files) {
  bn  <- basename(f)
  tag <- normalize_tag(bn)
  idx <- match(tag, design$image_tag)

  if (is.na(idx)) {
    log_info(sprintf("★ 无法匹配: tag=%s", tag)); next
  }

  sample_name <- design$sample[idx]
  grp         <- design$group[idx]
  slide_id    <- design$slide[idx]
  area_id     <- design$area[idx]

  log_info(sprintf("读取 %s → sample=%s, group=%s", bn, sample_name, grp))

  ok <- tryCatch({
    dt <- fread(f, data.table = FALSE)

    missing_cols <- setdiff(META_COLS, colnames(dt))
    if (length(missing_cols) > 0) {
      log_info(sprintf("★ 缺列: %s", paste(missing_cols, collapse = ", ")))
      stop("missing cols")
    }

    meta      <- dt[, META_COLS, drop = FALSE]
    gene_cols <- setdiff(colnames(dt), META_COLS)

    expr_mat <- as.matrix(dt[, gene_cols, drop = FALSE])
    expr_mat <- t(expr_mat)
    colnames(expr_mat) <- meta$barcode
    rownames(expr_mat) <- make.unique(gene_cols)

    seu <- CreateSeuratObject(counts = expr_mat, project = sample_name,
                              min.cells = 0, min.features = MIN_SPOTS)

    keep_cells <- colnames(seu)
    meta <- meta[match(keep_cells, meta$barcode), , drop = FALSE]
    if (any(is.na(meta$barcode))) {
      log_info("★ meta 对齐失败，跳过该样本"); stop("meta align failed")
    }

    seu <- set_assay_data_compat(seu, expr_mat[, keep_cells, drop = FALSE])

    seu$sample   <- sample_name
    seu$group    <- grp
    seu$slide    <- slide_id
    seu$area     <- area_id
    seu$image    <- tag
    seu$prediction_celltype <- meta$prediction_celltype
    seu$region   <- meta$region
    seu$tissue   <- meta$tissue
    seu$percent_CNV <- meta$percent_CNV

    for (ct in CELLTYPE_COLS) {
      seu@meta.data[[paste0("prop_", ct)]] <- meta[[ct]]
    }

    seu_list[[sample_name]] <- seu
    sample_meta_list[[sample_name]] <- data.frame(
      sample = sample_name, group = grp, slide = slide_id,
      area = area_id, tag = tag, n_spots = ncol(seu),
      stringsAsFactors = FALSE)

    log_info(sprintf("  ✓ %s: %d spots", sample_name, ncol(seu)))
    TRUE
  }, error = function(e) {
    log_info(sprintf("★ 处理失败 %s: %s", bn, conditionMessage(e)))
    FALSE
  })

  if (!ok) next
}

if (length(seu_list) == 0) stop("所有样本处理失败")

sample_info <- do.call(rbind, sample_meta_list)
print(sample_info)

if (length(seu_list) == 1) {
  seu <- seu_list[[1]]
} else {
  seu <- merge(seu_list[[1]], y = seu_list[-1], add.cell.ids = names(seu_list))
}
seu <- tryCatch(JoinLayers(seu), error = function(e) seu)

log_info(sprintf("合并后: %d spots × %d genes", ncol(seu), nrow(seu)))

# ============================================================
# ============================================================
log_info("--- 3. 基因检查 ---")

present <- function(gs, ref = rownames(seu)) {
  ref_upper <- toupper(ref)
  gs_upper  <- toupper(gs)
  idx <- match(gs_upper, ref_upper)
  ref[idx[!is.na(idx)]]
}

ci_present  <- present(CI_GENES)
c3_present  <- present(CIII_GENES)
g16_present <- present(GENES_16)

log_info(sprintf("Complex I: %d/%d", length(ci_present),  length(CI_GENES)))
log_info(sprintf("Complex III: %d/%d", length(c3_present), length(CIII_GENES)))
log_info(sprintf("16-gene: %d/%d", length(g16_present), length(GENES_16)))
log_info(sprintf("NDUFB3: %s",
                 ifelse(any(toupper(rownames(seu)) == toupper(GENE_FOCUS)),
                        "存在", "缺失")))

# ============================================================
# 4. Module Score
# ============================================================
log_info("--- 4. Module Score ---")

score_sets <- list()
if (length(c3_present)  >= 2) score_sets[["ComplexIII"]] <- c3_present
if (length(ci_present)  >= 2) score_sets[["ComplexI"]]   <- ci_present
if (length(g16_present) >= 3) score_sets[["Sig16"]]      <- g16_present

if (length(score_sets) > 0) {
  seu <- AddModuleScore(seu, features = score_sets, name = "MS_", ctrl = 50)

  for (i in seq_along(score_sets)) {
    old_col <- paste0("MS_", i)
    new_col <- paste0(names(score_sets)[i], "_score")
    if (old_col %in% colnames(seu@meta.data)) {
      colnames(seu@meta.data)[match(old_col, colnames(seu@meta.data))] <- new_col
      log_info(sprintf("重命名: %s -> %s", old_col, new_col))
    } else {
      log_info(sprintf("★ 未找到列: %s", old_col))
    }
  }

  log_info(sprintf("score 列: %s",
                   paste(grep("_score$", colnames(seu@meta.data), value = TRUE),
                         collapse = ", ")))
} else {
  log_info("★ 无基因集满足条件，跳过 AddModuleScore")
}

# ============================================================
# 5. summary
# ============================================================
log_info("--- 5. 汇总 ---")

meta <- seu@meta.data
expr_vars <- unique(c(GENE_FOCUS, c3_present, ci_present, g16_present))
expr_vars <- expr_vars[expr_vars %in% rownames(seu)]
score_vars <- grep("_score$", colnames(meta), value = TRUE)
all_vars <- unique(c(expr_vars, score_vars))

meta_all <- cbind(meta, FetchData(seu, vars = expr_vars))

grp_tab <- meta_all %>%
  group_by(group) %>%
  summarise(across(all_of(all_vars), ~mean(.x, na.rm = TRUE)),
            n_spots = n(), .groups = "drop")
save_table(grp_tab, DIR_THIS, "GSE234047_by_group.csv")

tissue_tab <- meta_all %>%
  group_by(tissue, group) %>%
  summarise(across(all_of(all_vars), ~mean(.x, na.rm = TRUE)),
            n_spots = n(), .groups = "drop")
save_table(tissue_tab, DIR_THIS, "GSE234047_by_tissue.csv")

region_tab <- meta_all %>%
  group_by(region, group) %>%
  summarise(across(all_of(all_vars), ~mean(.x, na.rm = TRUE)),
            n_spots = n(), .groups = "drop")
save_table(region_tab, DIR_THIS, "GSE234047_by_region.csv")

sample_tab <- meta_all %>%
  group_by(sample, group, slide, area) %>%
  summarise(across(all_of(all_vars), ~mean(.x, na.rm = TRUE)),
            n_spots = n(), .groups = "drop")
save_table(sample_tab, DIR_THIS, "GSE234047_by_sample.csv")

retina_tab <- meta_all %>%
  filter(tissue %in% c("retina","rpe_choroid")) %>%
  group_by(tissue, group) %>%
  summarise(across(all_of(all_vars), ~mean(.x, na.rm = TRUE)),
            n_spots = n(), .groups = "drop")
save_table(retina_tab, DIR_THIS, "GSE234047_retina_rpe_by_group.csv")

# ============================================================
# ============================================================
log_info("--- 6. 与 07_3 交叉 ---")

DIR_KO <- file.path(DIR_RESULT, "07_3_scTenifoldKnk_NDUFB3_RGC", "tables")

# 6.1 Complex I cross
ci_file <- file.path(DIR_KO, "NDUFB3_vs_ComplexI_RGC.csv")
ci_cross <- NULL
if (file.exists(ci_file)) {
  ci_cross <- read.csv(ci_file)
  save_table(ci_cross, DIR_THIS, "NDUFB3_vs_ComplexI_KO.csv")
  log_info(sprintf("Complex I 交叉: %d 基因, 显著 %d",
                   nrow(ci_cross),
                   sum(ci_cross[["p.adj"]] < 0.05, na.rm = TRUE)))
} else {
  log_info(sprintf("★ 未找到: %s", ci_file))
}

# 6.2 Complex III cross
dr_file <- file.path(DIR_KO, "NDUFB3_diffRegulation_RGC.csv")
ciii_dr <- NULL
if (file.exists(dr_file)) {
  dr <- read.csv(dr_file)
  log_info(sprintf("diffRegulation: %d 基因", nrow(dr)))

  ciii_dr <- dr[toupper(dr$gene) %in% toupper(CIII_GENES), ]
  ciii_dr$sig <- ifelse(ciii_dr[["p.adj"]] < 0.05, "yes", "no")

  save_table(ciii_dr, DIR_THIS, "NDUFB3_vs_ComplexIII_KO.csv")
  log_info(sprintf("Complex III 在 NDUFB3 KO: 命中 %d, 显著 %d",
                   nrow(ciii_dr),
                   sum(ciii_dr[["p.adj"]] < 0.05, na.rm = TRUE)))
  print(ciii_dr[, c("gene", "Z", "p.value", "p.adj", "sig")])

  n_pos <- sum(ciii_dr$Z > 0, na.rm = TRUE)
  n_tot <- sum(!is.na(ciii_dr$Z))
  if (n_tot > 0) {
    bt <- binom.test(n_pos, n_tot, p = 0.5)
    log_info(sprintf("Complex III 方向: %d/%d 正 Z, binomial p = %.4f",
                     n_pos, n_tot, bt$p.value))
    dir_tab <- data.frame(
      metric = c("n_pos", "n_total", "pct_pos", "binom_p"),
      value  = c(n_pos, n_tot, round(100 * n_pos / n_tot, 1),
                 signif(bt$p.value, 4)))
    save_table(dir_tab, DIR_THIS, "ComplexIII_directionality_test.csv")
  }
} else {
  log_info(sprintf("★ 未找到: %s", dr_file))
}

g16_file <- file.path(DIR_KO, "NDUFB3_vs_16genes_RGC.csv")
g16_cross <- NULL
missing_16 <- character(0)
if (file.exists(g16_file)) {
  g16_cross <- read.csv(g16_file)
  save_table(g16_cross, DIR_THIS, "NDUFB3_vs_16genes_KO.csv")
  log_info(sprintf("16 基因交叉: %d 基因, 显著 %d",
                   nrow(g16_cross),
                   sum(g16_cross[["p.adj"]] < 0.05, na.rm = TRUE)))

  n_pos_16 <- sum(g16_cross$Z > 0, na.rm = TRUE)
  n_tot_16 <- sum(!is.na(g16_cross$Z))
  if (n_tot_16 > 0) {
    bt16 <- binom.test(n_pos_16, n_tot_16, p = 0.5)
    log_info(sprintf("16 基因方向: %d/%d 正 Z, binomial p = %.4f",
                     n_pos_16, n_tot_16, bt16$p.value))
  }

  g16_upper     <- toupper(GENES_16)
  present_upper <- toupper(g16_cross$gene)
  missing_16    <- GENES_16[!(g16_upper %in% present_upper)]

  log_info(sprintf("16 基因集在 scTenifoldKnk 中可评估: %d/%d",
                   nrow(g16_cross), length(GENES_16)))
  if (length(missing_16) > 0) {
    log_info(sprintf("未出现在 scTenifoldKnk 网络 (%d): %s",
                     length(missing_16),
                     paste(missing_16, collapse = ", ")))
  } else {
    log_info("所有 16 基因都在网络中")
  }
} else {
  log_info(sprintf("★ 未找到: %s", g16_file))
}

# ============================================================
# 7. visualization
# ============================================================
log_info("--- 7. 可视化 ---")

plot_vars <- unique(c(GENE_FOCUS, c3_present, score_vars))
plot_vars <- plot_vars[plot_vars %in% colnames(seu@meta.data)]

if (length(plot_vars) > 0) {
  p1 <- VlnPlot(seu, features = plot_vars, group.by = "group",
                ncol = 4, pt.size = 0) +
    plot_annotation(title = "GSE234047: Control vs CNV")
  save_fig(p1, DIR_THIS, "violin_by_group", width = 16, height = 10)
}

if (length(plot_vars) > 0) {
  p2 <- VlnPlot(seu, features = plot_vars, group.by = "tissue",
                ncol = 4, pt.size = 0) +
    plot_annotation(title = "GSE234047: by tissue")
  save_fig(p2, DIR_THIS, "violin_by_tissue", width = 16, height = 10)
}

seu_sub <- subset(seu, tissue %in% c("retina","rpe_choroid"))
if (ncol(seu_sub) > 0 && length(plot_vars) > 0) {
  p3 <- VlnPlot(seu_sub, features = plot_vars, group.by = "tissue",
                split.by = "group", ncol = 4, pt.size = 0) +
    plot_annotation(title = "Retina / RPE-choroid: tissue x group")
  save_fig(p3, DIR_THIS, "violin_retina_rpe_tissue_x_group",
           width = 16, height = 10)
}

if (GENE_FOCUS %in% colnames(meta_all) &&
    "ComplexIII_score" %in% colnames(meta_all)) {
  p4 <- ggplot(meta_all,
               aes(x = .data[[GENE_FOCUS]], y = ComplexIII_score,
                   color = group, shape = tissue)) +
    geom_point(alpha = 0.4, size = 0.8) +
    geom_smooth(method = "lm", se = FALSE, aes(group = group)) +
    theme_bw(base_size = 12) +
    labs(title = "NDUFB3 vs Complex III score",
         x = "NDUFB3", y = "Complex III score")
  save_fig(p4, DIR_THIS, "NDUFB3_vs_ComplexIII_scatter", width = 9, height = 6)
}

if (length(c3_present) > 0) {
  tissue_mean <- meta_all %>%
    filter(!is.na(tissue)) %>%
    group_by(tissue) %>%
    summarise(across(all_of(c3_present), ~mean(.x, na.rm = TRUE)),
              .groups = "drop") %>%
    tibble::column_to_rownames("tissue")
  tissue_long <- tissue_mean %>%
    tibble::rownames_to_column("tissue") %>%
    pivot_longer(-tissue, names_to = "gene", values_to = "expr")
  p5 <- ggplot(tissue_long, aes(x = gene, y = tissue, fill = expr)) +
    geom_tile() + scale_fill_viridis_c(option = "magma") +
    theme_minimal(base_size = 11) +
    labs(title = "Complex III genes across tissues", x = "", y = "")
  save_fig(p5, DIR_THIS, "ComplexIII_tissue_heatmap", width = 9, height = 5)
}

if (all(paste0("prop_", CELLTYPE_COLS) %in% colnames(meta_all))) {
  ct_long <- meta_all %>%
    select(tissue, group, all_of(paste0("prop_", CELLTYPE_COLS))) %>%
    pivot_longer(-c(tissue, group), names_to = "celltype", values_to = "prop") %>%
    mutate(celltype = sub("^prop_", "", celltype)) %>%
    group_by(tissue, celltype) %>%
    summarise(prop = mean(prop, na.rm = TRUE), .groups = "drop")
  p6 <- ggplot(ct_long, aes(x = tissue, y = prop, fill = celltype)) +
    geom_col(position = "stack") + theme_bw(base_size = 11) +
    labs(title = "Cell type composition by tissue",
         x = "", y = "Mean proportion") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_fig(p6, DIR_THIS, "celltype_composition_by_tissue",
           width = 10, height = 6)
}

if (!is.null(ciii_dr) && nrow(ciii_dr) > 0) {
  z_df <- data.frame(gene = ciii_dr$gene, Z = ciii_dr$Z,
                     set = "Complex III")

  if (!is.null(ci_cross) && "Z" %in% colnames(ci_cross)) {
    ci_top <- ci_cross[order(-abs(ci_cross$Z)), ]
    ci_top <- head(ci_top, 20)
    z_df <- rbind(z_df, data.frame(gene = ci_top$gene, Z = ci_top$Z,
                                   set = "Complex I (top 20)"))
  }
  if (!is.null(g16_cross) && "Z" %in% colnames(g16_cross)) {
    z_df <- rbind(z_df, data.frame(gene = g16_cross$gene, Z = g16_cross$Z,
                                   set = "16-gene"))
  }

  p7 <- ggplot(z_df, aes(x = reorder(gene, Z), y = Z, fill = set)) +
    geom_col() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    coord_flip() +
    facet_wrap(~ set, scales = "free_y") +
    theme_bw(base_size = 11) +
    labs(title = "NDUFB3 KO: Z-scores in Complex I / III / 16-gene",
         x = "", y = "Z")
  save_fig(p7, DIR_THIS, "KO_Zscore_comparison", width = 12, height = 8)
}

# ============================================================
# ============================================================
log_info("--- 8. TIF 图像 ---")

tif_gz <- list.files(DIR_SC, pattern = "\\.tif\\.gz$",
                     recursive = TRUE, full.names = TRUE)
img_list <- list()

for (tf in tif_gz) {
  bn  <- basename(tf)
  tag <- normalize_tag(bn)
  idx <- match(tag, design$image_tag)
  if (is.na(idx)) next
  sample_name <- design$sample[idx]

  img <- tryCatch({
    tmp <- tempfile(fileext = ".tif")
    R.utils::gunzip(tf, tmp, remove = FALSE, overwrite = TRUE)
    img <- EBImage::readImage(tmp)
    unlink(tmp)
    img
  }, error = function(e) {
    log_info(sprintf("★ 失败: %s — %s", bn, conditionMessage(e)))
    NULL
  })

  if (!is.null(img)) {
    img_list[[sample_name]] <- img
    thumb_file <- file.path(DIR_THIS, "figures",
                            paste0("thumb_", sample_name, ".png"))
    try(EBImage::writeImage(img, thumb_file, quality = 85), silent = TRUE)
    log_info(sprintf("  ✓ %s", sample_name))
  }
}

if (length(img_list) > 0) {
  pdf(file.path(DIR_THIS, "figures", "tissue_overview.pdf"),
      width = 12, height = 10)
  n <- length(img_list); nc <- min(3, n); nr <- ceiling(n / nc)
  par(mfrow = c(nr, nc), mar = c(1,1,3,1))
  for (sn in names(img_list)) {
    img <- img_list[[sn]]
    gray <- tryCatch({
      if (length(dim(img)) >= 3 && dim(img)[3] >= 3)
        EBImage::channel(img, "gray") else img
    }, error = function(e) img)
    plot(EBImage::resize(gray, w = 400), main = sn)
  }
  dev.off()
}

# ============================================================
# 9. save
# ============================================================
log_info("--- 9. 保存 ---")

saveRDS(seu, file.path(DIR_THIS, "GSE234047_seurat.rds"))
if (length(img_list) > 0)
  saveRDS(img_list, file.path(DIR_THIS, "GSE234047_images.rds"))

write.csv(meta_all,
          file.path(DIR_THIS, "tables", "GSE234047_metadata_full.csv"),
          row.names = FALSE)

# ============================================================
# 10. summary
# ============================================================
cat("\n", strrep("=", 70), "\n", sep = "")
cat("=== 07_9 GSE234047 spatial validation v3.6 汇总 ===\n")
cat(strrep("=", 70), "\n", sep = "")

cat(sprintf("\n【1】样本: %d 切片\n", nrow(sample_info)))
print(sample_info)

cat(sprintf("\n【2】spot 总数: %d\n", ncol(seu)))

cat("\n【3】tissue 分布:\n"); print(table(seu$tissue))
cat("\n【4】region 分布:\n");  print(table(seu$region))
cat("\n【5】group 分布:\n");   print(table(seu$group))

cat(sprintf("\n【6】NDUFB3 检出 spot: %d / %d (%.2f%%)\n",
            sum(FetchData(seu, vars = GENE_FOCUS)[, 1] > 0),
            ncol(seu),
            100 * mean(FetchData(seu, vars = GENE_FOCUS)[, 1] > 0)))

cat("\n【7】score 列:\n")
print(grep("_score$", colnames(seu@meta.data), value = TRUE))

cat("\n【8】按 group 均值:\n")
print(as.data.frame(grp_tab)[, c("group",
  intersect(c("NDUFB3","UQCRQ","ComplexIII_score","ComplexI_score","Sig16_score"),
            colnames(grp_tab)))],
  row.names = FALSE, digits = 4)

cat("\n【9】retina + rpe_choroid 按 tissue x group 均值:\n")
print(as.data.frame(retina_tab)[, c("tissue","group",
  intersect(c("NDUFB3","ComplexIII_score","ComplexI_score","Sig16_score"),
            colnames(retina_tab)))],
  row.names = FALSE, digits = 4)

cat("\n", strrep("-", 70), "\n", sep = "")
cat("【10】与 07_3 scTenifoldKnk NDUFB3 KO 交叉\n")
cat(strrep("-", 70), "\n", sep = "")

if (!is.null(ci_cross)) {
  cat(sprintf("\nComplex I 交叉: %d 基因, 显著 %d\n",
              nrow(ci_cross),
              sum(ci_cross[["p.adj"]] < 0.05, na.rm = TRUE)))
}
if (!is.null(ciii_dr)) {
  cat(sprintf("\nComplex III 交叉: %d 基因, 显著 %d\n",
              nrow(ciii_dr),
              sum(ciii_dr[["p.adj"]] < 0.05, na.rm = TRUE)))
  print(ciii_dr[, c("gene","Z","p.adj")], row.names = FALSE, digits = 4)
  n_pos <- sum(ciii_dr$Z > 0, na.rm = TRUE)
  n_tot <- sum(!is.na(ciii_dr$Z))
  if (n_tot > 0) {
    cat(sprintf("方向: %d/%d 正 Z, binomial p = %.4f\n",
                n_pos, n_tot, binom.test(n_pos, n_tot)$p.value))
  }
}
if (!is.null(g16_cross)) {
  cat(sprintf("\n16 基因交叉: %d / %d 可评估, 显著 %d\n",
              nrow(g16_cross), length(GENES_16),
              sum(g16_cross[["p.adj"]] < 0.05, na.rm = TRUE)))
  n_pos16 <- sum(g16_cross$Z > 0, na.rm = TRUE)
  n_tot16 <- sum(!is.na(g16_cross$Z))
  if (n_tot16 > 0) {
    cat(sprintf("方向: %d/%d 正 Z, binomial p = %.4f\n",
                n_pos16, n_tot16, binom.test(n_pos16, n_tot16)$p.value))
  }
  if (length(missing_16) > 0) {
    cat(sprintf("未出现在 scTenifoldKnk 网络 (%d): %s\n",
                length(missing_16),
                paste(missing_16, collapse = ", ")))
  } else {
    cat("所有 16 基因都在网络中\n")
  }
}

cat(sprintf("\n【11】图像读取: %d 张\n", length(img_list)))

cat("\n--- 输出位置 ---\n")
cat("目录:", DIR_THIS, "\n")
log_info("=== 07_9 GSE234047 v3.6 done ===")

