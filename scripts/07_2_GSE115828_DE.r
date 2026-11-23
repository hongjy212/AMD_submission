#!/usr/bin/env Rscript
# ============================================================
# ============================================================

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(ggplot2)
})

DIR_THIS <- setup_output_dir("07_2_GSE115828_DE")
log_info("=== 07_2_GSE115828_DE.r v2 start ===")

# ---- 1. read ----
xlsx_path <- "E:/retina_tools/raw/GSE115828/NIHMS1518564-supplement-6.xlsx"
raw_all <- read_excel(xlsx_path, sheet = "DifferentialExpression", col_names = FALSE)
hdr_row <- which(apply(raw_all, 1, function(r) "gene_id" %in% as.character(unlist(r))))[1]

de_raw <- read_excel(xlsx_path, sheet = "DifferentialExpression", skip = hdr_row - 1)
colnames(de_raw)[1] <- "gene_id"
colnames(de_raw)[2] <- "gene_name"

logFC_all <- grep("^logFC", colnames(de_raw), value = TRUE)
adjP_all  <- grep("^adj\\.P\\.Val", colnames(de_raw), value = TRUE)

traj <- data.frame(
  gene_id   = de_raw$gene_id,
  gene_name = toupper(de_raw$gene_name),
  G1 = 0,
  G2 = de_raw[[logFC_all[1]]],
  G3 = de_raw[[logFC_all[2]]],
  G4 = de_raw[[logFC_all[3]]],
  adjP_2_1 = de_raw[[adjP_all[1]]],
  adjP_3_1 = de_raw[[adjP_all[2]]],
  adjP_4_1 = de_raw[[adjP_all[3]]],
  stringsAsFactors = FALSE)
traj <- traj[complete.cases(traj[, c("G2", "G3", "G4")]), ]
log_info(paste0("有效基因: ", nrow(traj)))

classify_traj <- function(g1, g2, g3, g4, min_amp = 0.05) {
  d12 <- g2 - g1
  d23 <- g3 - g2
  d34 <- g4 - g3

  g3_max <- g3 > g2 && g3 > g4
  g3_min <- g3 < g2 && g3 < g4
  amp   <- min(abs(g3 - g2), abs(g3 - g4))

  if ((g3_max || g3_min) && amp >= min_amp) {
    return(ifelse(g3_max, "G3_up", "G3_down"))
  }

  if ((d12 >= 0 && d23 >= 0 && d34 >= 0) ||
      (d12 <= 0 && d23 <= 0 && d34 <= 0)) {
    return("monotonic")
  }

  mid_avg <- (g2 + g3) / 2
  if ((g1 - mid_avg) * (g4 - mid_avg) > 0) {
    return("U_shape")
  }

  "other"
}

traj$shape <- mapply(classify_traj, traj$G1, traj$G2, traj$G3, traj$G4,
                     min_amp = 0.05)

traj$amp_23 <- abs(traj$G3 - traj$G2) + abs(traj$G3 - traj$G4)
traj$amp_12 <- abs(traj$G2 - traj$G1)
traj$amp_34 <- abs(traj$G4 - traj$G3)
traj$is_g3_dominant <- traj$amp_23 > (traj$amp_12 + traj$amp_34)

traj$any_sig <- (traj$adjP_2_1 < 0.05) | (traj$adjP_3_1 < 0.05) | (traj$adjP_4_1 < 0.05)

save_table(traj, DIR_THIS, "shape_from_paper_DE.csv")

# ---- 4. summary ----
cat("\n=== 原论文 DE 数据的轨迹形状分布（全部基因）===\n")
print(table(traj$shape, useNA = "ifany"))

cat("\n=== 限制到至少一个对比 adj.P.Val<0.05 的基因 ===\n")
sub_sig <- traj[traj$any_sig, ]
cat("显著基因数:", nrow(sub_sig), " / ", nrow(traj), "\n")
print(table(sub_sig$shape, useNA = "ifany"))

GENES_16 <- c("FCF1","UGGT1","CCDC3","S100B","BPNT1","BRI3BP","BTBD19",
              "DPH7","AKR1B1","ZNF623","ZNF587","RSPH3","TAF3","NDUFB8",
              "ZNF417","C6orf120")
g16 <- traj[traj$gene_name %in% GENES_16, ]
save_table(g16, DIR_THIS, "gene16_in_paper_DE.csv")

cat("\n=== 16 基因在原论文 DE 数据里的形状 ===\n")
print(g16[, c("gene_name", "G2", "G3", "G4", "shape",
               "adjP_2_1", "adjP_3_1", "adjP_4_1")],
      row.names = FALSE, digits = 4)

GENE_SETS <- read_upstream("01_2_prepare_gene_sets", "data", "gene_sets.rds")
polyF <- read.csv(file.path(DIR_RESULT, "04_2_pathway_score", "tables",
                             "pathway_poly3F_adjusted.csv"))
p3_only_paths <- intersect(polyF$pathway[polyF$pattern == "p3_only"],
                            names(GENE_SETS))

pm <- data.frame()
for (pw in p3_only_paths) {
  gs_sym <- toupper(GENE_SETS[[pw]]$genes)
  sub <- traj[traj$gene_name %in% gs_sym, ]
  if (nrow(sub) < 5) next

  pm <- rbind(pm, data.frame(
    pathway = pw, n_members = nrow(sub),
    n_G3_up = sum(sub$shape == "G3_up", na.rm = TRUE),
    n_G3_down = sum(sub$shape == "G3_down", na.rm = TRUE),
    n_monotonic = sum(sub$shape == "monotonic", na.rm = TRUE),
    n_U_shape = sum(sub$shape == "U_shape", na.rm = TRUE),
    n_other = sum(sub$shape == "other", na.rm = TRUE),
    pct_G3 = round(100 * mean(sub$shape %in% c("G3_up", "G3_down"), na.rm = TRUE), 1),
    stringsAsFactors = FALSE))
}
pm <- pm[order(-pm$pct_G3), ]
save_table(pm, DIR_THIS, "pathway_members_paper_DE.csv")

cat("\n=== 19 通路成员在原论文 DE 里的 G3 特异比例 ===\n")
print(pm, row.names = FALSE)

# ---- 7. figure ----
p1 <- ggplot(data.frame(shape = traj$shape), aes(x = shape, fill = shape)) +
  geom_bar() + theme_bw(base_size = 11) +
  labs(title = "Original DE data: trajectory shape distribution",
       x = "Shape", y = "Count")
save_fig(p1, DIR_THIS, "shape_summary", width = 6, height = 4)

# ---- 8. summary ----
cat("\n", strrep("=", 70), "\n", sep = "")
cat("=== 07_2_GSE115828_DE v2 汇总 ===\n")
cat(strrep("=", 70), "\n", sep = "")
cat("基因数:", nrow(traj), "\n")
cat("全部 pattern 分布:\n"); print(table(traj$shape))
cat("\n显著基因 pattern 分布:\n"); print(table(sub_sig$shape))
cat("\n16 基因 pattern:\n"); print(table(g16$shape))
cat("\n目录:", DIR_THIS, "\n")
log_info("=== 07_2_GSE115828_DE.r v2 done ===")

