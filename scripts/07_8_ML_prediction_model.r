# ============================================================
# Script name: 07_8_ML_prediction_model
# Purpose: ML classification based on 16-gene signature + greenyellow 
#          + black modules, with two external validation cohorts 
#          (GSE29801 + GSE135092) and batch sensitivity assessment
# Features: 330 (16-gene + greenyellow + black)
# Mode: G1 vs G4
# External cohorts:
#   - GSE29801 (retina, from 07_2_gse29801_replication)
#   - GSE135092 (retina + RPE, from 07_8a_GSE135092_prep)
#   - GSE135092: log2-transformed to align with discovery cohort
# Data policy: only non-mw26 directories
# Output: E:/retina_tools/aging_rgc_axon/AMD/results/07_8_ML_prediction_model/
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(glmnet)
  library(ranger)
  library(xgboost)
  library(pROC)
})

root_dir <- "E:/retina_tools/aging_rgc_axon/AMD/results"
out_dir  <- file.path(root_dir, "07_8_ML_prediction_model")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

genes_16 <- c("NDUFB8", "UGGT1", "C6orf120", "DPH7", "AKR1B1",
              "RSPH3", "FCF1", "BRI3BP", "ZNF623", "ZNF417",
              "ZNF587", "TAF3", "BPNT1", "BTBD19", "S100B", "CCDC3")

# ============================================================
# 1. Load training data
# ============================================================
cat("========== 1. Load data ==========\n")

expr_train <- readRDS(file.path(root_dir,
  "01_1_prepare_bulk_data/data/expr_corrected.rds"))
sample_tbl <- readRDS(file.path(root_dir,
  "05_1_consensus_wgcna/data/sample_table.rds"))

id_map_raw <- readRDS(file.path(root_dir,
  "01_1_prepare_bulk_data/data/gene_annotation.rds")) |> as.data.frame()

col_ens <- intersect(c("ensembl","ENSEMBL","ensembl_id","Ensembl_ID","gene_id"),
                     names(id_map_raw))[1]
col_sym <- intersect(c("symbol","SYMBOL","gene_symbol","gene_name","Symbol"),
                     names(id_map_raw))[1]

if (is.na(col_ens) || is.na(col_sym)) {
  stop("Cannot find Ensembl/Symbol columns in gene_annotation.rds")
}

id_map <- id_map_raw |>
  select(ensembl = all_of(col_ens), symbol = all_of(col_sym)) |>
  filter(!is.na(symbol)) |> distinct()

cat("expr dim:", paste(dim(expr_train), collapse=" x "), "\n")
cat("sample_table dim:", paste(dim(sample_tbl), collapse=" x "), "\n")

# ============================================================
# 2. Load module genes
# ============================================================
cat("\n========== 2. Load module genes ==========\n")

gy_ens <- read_csv(file.path(root_dir,
  "05_1_consensus_wgcna/tables/module_genes_greenyellow.csv"),
  show_col_types = FALSE) |> pull(ENSEMBL) |> unique()

bl_ens <- read_csv(file.path(root_dir,
  "05_1_consensus_wgcna/tables/module_genes_black.csv"),
  show_col_types = FALSE) |> pull(ENSEMBL) |> unique()

cat("greenyellow:", length(gy_ens), "\n")
cat("black:", length(bl_ens), "\n")

ens_16 <- id_map |> filter(symbol %in% genes_16) |>
  pull(ensembl) |> unique() |> intersect(rownames(expr_train))
cat("16-gene Ensembl:", length(ens_16), "/ 16\n")

feature_ens <- unique(c(ens_16, gy_ens, bl_ens))
feature_ens <- intersect(feature_ens, rownames(expr_train))

cat("\nTotal features:", length(feature_ens), "\n")
cat("  - 16-gene:", length(ens_16), "\n")
cat("  - greenyellow in expr:", length(intersect(gy_ens, rownames(expr_train))), "\n")
cat("  - black in expr:", length(intersect(bl_ens, rownames(expr_train))), "\n")

# ============================================================
# 3. Labels G1 vs G4
# ============================================================
cat("\n========== 3. Labels: G1 vs G4 ==========\n")

sample_tbl_use <- sample_tbl |>
  filter(Condition %in% c("G1", "G4")) |>
  mutate(label = ifelse(Condition == "G1", 0L, 1L))

print(table(sample_tbl_use$Condition, sample_tbl_use$label))

# ============================================================
# 4. Build training set
# ============================================================
common_samples <- intersect(colnames(expr_train), sample_tbl_use$Sample)

X_train <- t(expr_train[feature_ens, common_samples, drop = FALSE])
X_train <- as.data.frame(X_train)

sym_lookup <- id_map$symbol[match(colnames(X_train), id_map$ensembl)]
sym_lookup[is.na(sym_lookup)] <- colnames(X_train)[is.na(sym_lookup)]
colnames(X_train) <- make.unique(sym_lookup)

df_train <- X_train |>
  rownames_to_column("Sample") |>
  inner_join(sample_tbl_use |> select(Sample, label), by = "Sample") |>
  select(-Sample) |>
  mutate(label = factor(label, levels = c(0, 1)))

genes_avail <- setdiff(colnames(df_train), "label")

cat("\nTraining set:", nrow(df_train), "samples,",
    sum(df_train$label == 1), "cases,",
    sum(df_train$label == 0), "controls\n")
cat("Features:", length(genes_avail), "\n\n")

# ============================================================
# 5. Internal 5-fold CV
# ============================================================
cat("========== 4. Internal 5-fold CV ==========\n")

set.seed(42)
k <- 5
fold_id <- sample(rep(1:k, length.out = nrow(df_train)))
folds <- split(seq_len(nrow(df_train)), fold_id)

cv_auc <- function(method) {
  aucs <- c()
  for (i in seq_along(folds)) {
    test_idx  <- folds[[i]]
    train_idx <- setdiff(seq_len(nrow(df_train)), test_idx)

    X_tr <- as.matrix(df_train[train_idx, genes_avail])
    y_tr <- as.integer(as.numeric(df_train$label[train_idx]) - 1)
    X_te <- as.matrix(df_train[test_idx, genes_avail])
    y_te <- as.integer(as.numeric(df_train$label[test_idx]) - 1)

    pred <- tryCatch({
      if (method == "lasso") {
        fit <- glmnet::cv.glmnet(X_tr, y_tr, family = "binomial", alpha = 1)
        as.numeric(predict(fit, X_te, s = "lambda.min", type = "response"))
      } else if (method == "rf") {
        fit <- ranger::ranger(x = X_tr, y = factor(y_tr),
                              num.trees = 500, probability = TRUE)
        predict(fit, data = as.data.frame(X_te))$predictions[, "1"]
      } else {
        dtrain <- xgboost::xgb.DMatrix(data = X_tr, label = y_tr)
        dtest  <- xgboost::xgb.DMatrix(data = X_te)
        fit <- xgboost::xgb.train(
          data = dtrain, nrounds = 200,
          params = list(objective = "binary:logistic",
                        eval_metric = "auc"),
          verbose = 0
        )
        as.numeric(predict(fit, dtest))
      }
    }, error = function(e) {
      cat("  Error in", method, "fold", i, ":", conditionMessage(e), "\n")
      rep(NA_real_, length(y_te))
    })

    if (all(is.na(pred))) { aucs[i] <- NA; next }
    aucs[i] <- as.numeric(pROC::auc(pROC::roc(y_te, pred, quiet = TRUE)))
  }
  mean(aucs, na.rm = TRUE)
}

internal_auc <- tibble(
  feature_set = "16+greenyellow+black",
  n_features = length(genes_avail),
  model = c("LASSO", "RandomForest", "XGBoost"),
  internal_cv_auc = c(cv_auc("lasso"), cv_auc("rf"), cv_auc("xgb"))
)
print(internal_auc)
write_csv(internal_auc, file.path(out_dir,
  "Table1_internal_CV_AUC_module.csv"))

# ============================================================
# 6. Final fit on all data
# ============================================================
cat("\n========== 5. Final fit ==========\n")

X_all <- as.matrix(df_train[, genes_avail])
y_all <- as.integer(as.numeric(df_train$label) - 1)

# ---- 6.1 LASSO ----
lasso_fit <- glmnet::cv.glmnet(X_all, y_all, family = "binomial", alpha = 1)

lasso_coef_mat <- as.matrix(coef(lasso_fit, s = "lambda.min"))
lasso_coef <- data.frame(
  gene = rownames(lasso_coef_mat),
  coefficient = lasso_coef_mat[, 1],
  row.names = NULL
) |>
  filter(gene != "(Intercept)", coefficient != 0) |>
  arrange(desc(abs(coefficient)))

cat("\nLASSO selected:", nrow(lasso_coef), "genes\n")
print(head(lasso_coef, 20))
write_csv(lasso_coef, file.path(out_dir,
  "Table2_LASSO_coefficients_module.csv"))

# ---- 6.2 XGBoost ----
dtrain_all <- xgboost::xgb.DMatrix(data = X_all, label = y_all)
xgb_fit <- xgboost::xgb.train(
  data = dtrain_all, nrounds = 200,
  params = list(objective = "binary:logistic", eval_metric = "auc"),
  verbose = 0
)
xgb_imp <- xgboost::xgb.importance(model = xgb_fit)
cat("\nTop 20 XGBoost features:\n")
print(head(as_tibble(xgb_imp), 20))
write_csv(as_tibble(xgb_imp), file.path(out_dir,
  "Table4_XGB_importance_module.csv"))

# ---- 6.3 Random Forest ----
rf_fit <- ranger::ranger(
  x = X_all, y = factor(y_all),
  num.trees = 500, probability = TRUE,
  importance = "impurity"
)
rf_imp <- tibble(
  Feature = names(rf_fit$variable.importance),
  Importance = as.numeric(rf_fit$variable.importance)
) |> arrange(desc(Importance))
cat("\nTop 20 RandomForest features:\n")
print(head(rf_imp, 20))
write_csv(rf_imp, file.path(out_dir,
  "Table3_RF_importance_module.csv"))

# ============================================================
# 7. Generic external evaluation function
# ============================================================
evaluate_external <- function(model, model_name, meta_df, expr_mat,
                              feature_genes, cohort_name,
                              id_col = "sample") {
  common_genes <- intersect(feature_genes, rownames(expr_mat))
  cat(" ", cohort_name, "- matching genes:", length(common_genes), "/",
      length(feature_genes), "\n")

  if (length(common_genes) < 5) return(NULL)

  X_ext <- t(expr_mat[common_genes, , drop = FALSE])
  X_ext <- as.data.frame(X_ext)
  missing <- setdiff(feature_genes, colnames(X_ext))
  for (g in missing) X_ext[[g]] <- 0
  X_ext <- X_ext[, feature_genes, drop = FALSE]

  meta_use <- meta_df |>
    rename(sample = all_of(id_col)) |>
    select(sample, label)

  df_ext <- X_ext |>
    rownames_to_column("sample") |>
    inner_join(meta_use, by = "sample") |>
    select(-sample) |>
    mutate(label = factor(label, levels = c(0, 1)))

  if (nrow(df_ext) < 10 ||
      sum(df_ext$label == 1) < 3 ||
      sum(df_ext$label == 0) < 3) return(NULL)

  X_te <- as.matrix(df_ext[, feature_genes])
  y_te <- as.integer(as.numeric(df_ext$label) - 1)

  pred <- tryCatch({
    if (model_name == "LASSO") {
      as.numeric(predict(model, X_te, s = "lambda.min", type = "response"))
    } else if (model_name == "XGBoost") {
      as.numeric(predict(model, xgboost::xgb.DMatrix(data = X_te)))
    } else if (model_name == "RandomForest") {
      predict(model, as.data.frame(X_te))$predictions[, "1"]
    }
  }, error = function(e) {
    cat("  Error:", conditionMessage(e), "\n")
    rep(NA_real_, length(y_te))
  })

  if (all(is.na(pred))) return(NULL)

  roc_obj <- pROC::roc(y_te, pred, quiet = TRUE)
  ci_obj  <- pROC::ci.auc(roc_obj)

  tibble(
    cohort    = cohort_name,
    model     = model_name,
    auc       = as.numeric(pROC::auc(roc_obj)),
    ci_low    = ci_obj[1],
    ci_high   = ci_obj[3],
    n_total   = length(y_te),
    n_case    = sum(y_te == 1),
    n_control = sum(y_te == 0)
  )
}

# ============================================================
# 8. External cohort 1: GSE29801 (retina, G1 vs G4)
# ============================================================
cat("\n========== 6. External GSE29801 ==========\n")

expr_ext <- readRDS(file.path(root_dir,
  "07_2_gse29801_replication/data/expr_log2.rds"))
meta_ext <- readRDS(file.path(root_dir,
  "07_2_gse29801_replication/data/meta.rds")) |> as.data.frame()

meta_ext <- meta_ext |>
  filter(grade %in% c("G1", "G4")) |>
  mutate(label = ifelse(grade == "G1", 0L, 1L))

cat("GSE29801 G1+G4 samples:", nrow(meta_ext), "\n")
print(table(meta_ext$grade, meta_ext$label))

gse29801_res <- bind_rows(
  evaluate_external(lasso_fit, "LASSO", meta_ext, expr_ext,
                    genes_avail, "GSE29801_G1G4", id_col = "sample"),
  evaluate_external(xgb_fit, "XGBoost", meta_ext, expr_ext,
                    genes_avail, "GSE29801_G1G4", id_col = "sample"),
  evaluate_external(rf_fit, "RandomForest", meta_ext, expr_ext,
                    genes_avail, "GSE29801_G1G4", id_col = "sample")
)
print(gse29801_res)
write_csv(gse29801_res, file.path(out_dir,
  "Table5_external_AUC_GSE29801.csv"))

# ============================================================
# 9. External cohort 2: GSE135092 (retina + RPE, log2-transformed)
# ============================================================
cat("\n========== 7. External GSE135092 ==========\n")

gse135092_data_dir <- file.path(root_dir, "07_8a_GSE135092_prep/data")

if (!dir.exists(gse135092_data_dir)) {
  cat("07_8a_GSE135092_prep not found. Skipping GSE135092.\n")
  gse135092_res <- tibble()
} else {
  expr_135092 <- readRDS(file.path(gse135092_data_dir, "expr_symbol.rds"))
  meta_135092 <- readRDS(file.path(gse135092_data_dir, "clinical.rds"))

  cat("Before log2: range =", paste(round(range(expr_135092, na.rm=TRUE), 2),
                                     collapse=" - "), "\n")
  expr_135092 <- log2(expr_135092 + 1)
  cat("After  log2: range =", paste(round(range(expr_135092, na.rm=TRUE), 2),
                                     collapse=" - "), "\n")

  cat("expr dim:", paste(dim(expr_135092), collapse=" x "), "\n")
  cat("meta dim:", paste(dim(meta_135092), collapse=" x "), "\n")

  # Macular RPE
  macular_rpe <- meta_135092 |>
    filter(region == "Macula", tissue_type == "RPE") |>
    mutate(label = ifelse(amd_binary == "Control", 0L, 1L))

  # Macular Retina
  macular_retina <- meta_135092 |>
    filter(region == "Macula", tissue_type == "Retina") |>
    mutate(label = ifelse(amd_binary == "Control", 0L, 1L))

  cat("\nMacular RPE:", nrow(macular_rpe),
      "(Control", sum(macular_rpe$label==0),
      ", AMD", sum(macular_rpe$label==1), ")\n")
  cat("Macular Retina:", nrow(macular_retina),
      "(Control", sum(macular_retina$label==0),
      ", AMD", sum(macular_retina$label==1), ")\n\n")

  gse135092_res <- bind_rows(
    evaluate_external(lasso_fit, "LASSO", macular_rpe,
                      expr_135092, genes_avail, "GSE135092_macular_RPE",
                      id_col = "gsm"),
    evaluate_external(xgb_fit, "XGBoost", macular_rpe,
                      expr_135092, genes_avail, "GSE135092_macular_RPE",
                      id_col = "gsm"),
    evaluate_external(rf_fit, "RandomForest", macular_rpe,
                      expr_135092, genes_avail, "GSE135092_macular_RPE",
                      id_col = "gsm"),
    evaluate_external(lasso_fit, "LASSO", macular_retina,
                      expr_135092, genes_avail, "GSE135092_macular_Retina",
                      id_col = "gsm"),
    evaluate_external(xgb_fit, "XGBoost", macular_retina,
                      expr_135092, genes_avail, "GSE135092_macular_Retina",
                      id_col = "gsm"),
    evaluate_external(rf_fit, "RandomForest", macular_retina,
                      expr_135092, genes_avail, "GSE135092_macular_Retina",
                      id_col = "gsm")
  )
  print(gse135092_res)
  write_csv(gse135092_res, file.path(out_dir,
    "Table5b_external_AUC_GSE135092.csv"))
}

external_res_all <- bind_rows(gse29801_res, gse135092_res)
print(external_res_all)
write_csv(external_res_all, file.path(out_dir,
  "Table5c_external_AUC_all_cohorts.csv"))

# ============================================================
# 10. Batch sensitivity of LASSO-selected genes
# ============================================================
cat("\n========== 8. Batch sensitivity of LASSO genes ==========\n")

batch_16 <- read_csv(file.path(root_dir,
  "diag/diag_batch_sensitivity/tables/poly3F_16genes_batch.csv"),
  show_col_types = FALSE)

lasso_genes <- lasso_coef$gene

overlap_16 <- intersect(lasso_genes, batch_16$feature)
cat("LASSO genes overlapping with 16-gene signature:",
    length(overlap_16), "\n")
print(overlap_16)

overlap_detail <- batch_16 |>
  filter(feature %in% overlap_16) |>
  select(feature, FDR_nobatch, FDR_batch) |>
  mutate(
    batch_robust = FDR_batch < 0.05,
    status = ifelse(batch_robust, "Robust after batch", "Batch-sensitive")
  )
cat("\nOverlap genes batch status:\n")
print(overlap_detail)
write_csv(overlap_detail, file.path(out_dir,
  "Table6_LASSO_16gene_batch_status.csv"))

gy_sym <- id_map$symbol[match(gy_ens, id_map$ensembl)]
bl_sym <- id_map$symbol[match(bl_ens, id_map$ensembl)]

lasso_in_gy <- intersect(lasso_genes, gy_sym)
lasso_in_bl <- intersect(lasso_genes, bl_sym)
lasso_in_16 <- overlap_16

breakdown <- tibble(
  source = c("16-gene signature", "greenyellow", "black", "Total unique"),
  n_genes = c(length(lasso_in_16), length(lasso_in_gy),
              length(lasso_in_bl),
              length(unique(c(lasso_in_16, lasso_in_gy, lasso_in_bl))))
)
cat("\n===== LASSO genes breakdown =====\n")
print(breakdown)
cat("\nIn greenyellow:\n"); print(lasso_in_gy)
cat("\nIn black:\n"); print(lasso_in_bl)
write_csv(breakdown, file.path(out_dir,
  "Table7_LASSO_gene_source_breakdown.csv"))

# ---- completeannotation ----
lasso_annotated <- lasso_coef |>
  mutate(
    in_16gene = gene %in% lasso_in_16,
    in_greenyellow = gene %in% lasso_in_gy,
    in_black = gene %in% lasso_in_bl,
    source = case_when(
      in_16gene & in_greenyellow ~ "16gene+greenyellow",
      in_16gene & in_black ~ "16gene+black",
      in_16gene ~ "16gene",
      in_greenyellow ~ "greenyellow",
      in_black ~ "black",
      TRUE ~ "other"
    )
  )
write_csv(lasso_annotated, file.path(out_dir,
  "Table8_LASSO_annotated.csv"))

# ============================================================
# 11. Output summary
# ============================================================
cat("\n========== Output files ==========\n")
print(list.files(out_dir, full.names = FALSE))
cat("\nDirectory:", out_dir, "\n")

cat("\n===== Key results =====\n")
cat("Internal CV AUC (LASSO/RF/XGB):",
    paste(round(internal_auc$internal_cv_auc, 3), collapse=" / "), "\n\n")

cat("External AUC:\n")
print(external_res_all |>
  select(cohort, model, auc, ci_low, ci_high, n_total, n_case, n_control))

cat("\nLASSO selected:", nrow(lasso_coef), "genes\n")
cat("  - from 16-gene signature:", length(lasso_in_16), "\n")
cat("  - from greenyellow:", length(lasso_in_gy), "\n")
cat("  - from black:", length(lasso_in_bl), "\n")
cat("\nAll 16-gene signature LASSO members batch robust:",
    all(overlap_detail$batch_robust), "\n")
	
	
