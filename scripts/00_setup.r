# ============================================================
# 00_setup.r
# ============================================================

suppressPackageStartupMessages({
  library(glue)
})

DIR_ROOT   <- "E:/retina_tools/aging_rgc_axon/AMD"
DIR_RESULT <- file.path(DIR_ROOT, "results")
DIR_RAW    <- "E:/retina_tools/raw"

SEED    <- 20260925
N_PERM  <- 1000
COVARS  <- c("age", "sex", "batch", "RIN", "PMI",
             "PC1", "PC2", "PC3", "PC4", "PC5")

.log_ts <- function() format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
.log_msg <- function(level, msg) {
  txt <- tryCatch(as.character(glue::glue(msg)), error = function(e) msg)
  cat(sprintf("%s [%s] %s\n", .log_ts(), level, txt))
  invisible(txt)
}
log_info  <- function(msg) .log_msg("INFO", msg)
log_warn  <- function(msg) .log_msg("WARN", msg)
log_error <- function(msg) .log_msg("ERROR", msg)

setup_output_dir <- function(name, subdir = NULL) {
  path <- if (is.null(subdir) || subdir == "") file.path(DIR_RESULT, name)
          else file.path(DIR_RESULT, subdir, name)
  for (sub in c("data", "tables", "figures"))
    dir.create(file.path(path, sub), recursive = TRUE, showWarnings = FALSE)
  path
}

save_rds <- function(obj, dir_this, name) {
  if (!grepl("\\.rds$", name)) name <- paste0(name, ".rds")
  path <- file.path(dir_this, "data", name)
  saveRDS(obj, path)
  log_info("saved: {name}")
  invisible(path)
}

save_table <- function(df, dir_this, name) {
  if (!grepl("\\.csv$", name)) name <- paste0(name, ".csv")
  path <- file.path(dir_this, "tables", name)
  utils::write.csv(df, path, row.names = FALSE)
  log_info("saved: {name} ({nrow(df)} rows)")
  invisible(path)
}

save_fig <- function(p, dir_this, name, width = 8, height = 6, dpi = 300) {
  if (grepl("\\.(pdf|png)$", name)) name <- sub("\\.(pdf|png)$", "", name)
  pdf_path <- file.path(dir_this, "figures", paste0(name, ".pdf"))
  png_path <- file.path(dir_this, "figures", paste0(name, ".png"))
  if (inherits(p, "ggplot") || inherits(p, "gg")) {
    ggplot2::ggsave(pdf_path, plot = p, width = width, height = height, device = "pdf")
    ggplot2::ggsave(png_path, plot = p, width = width, height = height, dpi = dpi, device = "png")
  } else if (is.function(p)) {
    grDevices::pdf(pdf_path, width = width, height = height); p(); grDevices::dev.off()
    grDevices::png(png_path, width = width, height = height, units = "in", res = dpi); p(); grDevices::dev.off()
  } else {
    log_warn("save_fig: 不支持的对象类型")
  }
  log_info("saved: {name}.pdf/.png")
  invisible(c(pdf = pdf_path, png = png_path))
}

read_upstream <- function(script_name, type = "data", file_name, subdir = NULL) {
  if (!type %in% c("data", "tables", "figures"))
    stop("type 必须是 data / tables / figures 之一")
  base <- if (is.null(subdir) || subdir == "") file.path(DIR_RESULT, script_name, type)
          else file.path(DIR_RESULT, subdir, script_name, type)
  path <- file.path(base, file_name)
  if (!file.exists(path)) stop(sprintf("文件不存在: %s", path))
  if (grepl("\\.rds$", file_name)) readRDS(path)
  else if (grepl("\\.csv$", file_name)) utils::read.csv(path, stringsAsFactors = FALSE)
  else stop(sprintf("不支持的文件类型: %s", file_name))
}

log_info("=== 00_setup.r loaded ===")
log_info("DIR_ROOT   = {DIR_ROOT}")
log_info("DIR_RESULT = {DIR_RESULT}")
log_info("DIR_RAW    = {DIR_RAW}")
log_info("SEED = {SEED}, N_PERM = {N_PERM}")
log_info("COVARS = {paste(COVARS, collapse=', ')}")
