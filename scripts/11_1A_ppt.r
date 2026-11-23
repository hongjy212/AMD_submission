# 11_1A_ppt.R (v2.5)
# Purpose: assemble manuscript.docx + Figures 1-7 + Graphical Abstract
#          + supplementary figures (optional) into a single PPTX
# Output:  results/11_1A_ppt/

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(officer)
  library(dplyr); library(stringr); library(readr)
  library(pdftools)
})

has_magick <- requireNamespace("magick", quietly = TRUE)
if (!has_magick) {
  log_warn("magick not installed; fallback to pdf_pages_size (may distort aspect)")
}

if (!exists("log_info")) {
  log_info <- function(...) message(sprintf("[%s] [INFO] %s", Sys.time(), paste0(...)))
  log_warn <- function(...) message(sprintf("[%s] [WARN] %s", Sys.time(), paste0(...)))
}

results_dir <- DIR_RESULT
out_dir     <- file.path(results_dir, "11_1A_ppt")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
log_info("=== 11_1A PPT generator (v2.5) ===")

# ---- 1. Input paths ----
manuscript_file <- file.path(DIR_ROOT, "manuscript", "manuscript.docx")
fig_main_dir    <- file.path(results_dir, "10_1A_main_figures")
fig_supp_dir    <- file.path(results_dir, "10_1C_supp_figures")
fig_ga_dir      <- file.path(results_dir, "10_1E_graphical_abstract")

if (!file.exists(manuscript_file)) {
  stop(sprintf("manuscript.docx not found: %s", manuscript_file))
}
log_info(sprintf("manuscript file: %s", manuscript_file))

# ---- 2. Extract Abstract from manuscript.docx (Purpose/Results/Conclusions) ----
log_info("--- reading manuscript.docx ---")
doc_src <- read_docx(manuscript_file)
paras <- docx_summary(doc_src)
paras <- paras[paras$content_type == "paragraph", c("style_name", "text")]
paras <- paras[nzchar(paras$text), ]
log_info(sprintf("paragraphs: %d", nrow(paras)))

get_abstract_section <- function(keyword, max_chars = 1500) {
  idx <- which(grepl(keyword, paras$text, ignore.case = TRUE))
  if (length(idx) == 0) return("")
  next_idx <- idx[1] + 1
  txt <- gsub("\\*\\*", "", paras$text[idx[1]])
  txt <- sub("^(Purpose|Background|Methods|Results|Conclusions)\\.\\s*", "", txt)
  while (next_idx <= nrow(paras)) {
    line <- gsub("\\*\\*", "", paras$text[next_idx])
    if (grepl("^(Purpose|Background|Methods|Results|Conclusions)\\.", line)) break
    txt <- paste(txt, line, sep = " ")
    next_idx <- next_idx + 1
  }
  out <- substr(txt, 1, max_chars)
  out <- sub("\\.[^.]*$", ".", out)
  out
}

abstract_pur <- get_abstract_section("Purpose\\.", 600)
abstract_res <- get_abstract_section("Results\\.", 1200)
abstract_con <- get_abstract_section("Conclusions\\.", 800)

log_info(sprintf("Abstract Purpose: %d chars", nchar(abstract_pur)))
log_info(sprintf("Abstract Results: %d chars", nchar(abstract_res)))
log_info(sprintf("Abstract Conclusions: %d chars", nchar(abstract_con)))

# ---- 3. Main figures (titles aligned with legend.docx) ----
main_figs <- c(
  "Figure1" = "Figure 1. A multi-scale causal triangulation framework validated by complement positive controls",
  "Figure2" = "Figure 2. A non-linear transcriptomic framework nominates G2-G3 as a candidate inflection point",
  "Figure3" = "Figure 3. OXPHOS downregulation is directionally consistent but not significant; causal analyses do not support a strong causal role",
  "Figure4" = "Figure 4. Dual epigenetic and genetic axes of OXPHOS regulation",
  "Figure5" = "Figure 5. MR nominates NDUFB3 as a provisional, instrument-fragile, stage-specific hypothesis for early AMD",
  "Figure6" = "Figure 6. Single-cell complex redistribution in retinal ganglion cells",
  "Figure7" = "Figure 7. NDUFB3 virtual knockout reveals cell-type-specific network perturbation"
)

# ---- 4. Supplementary figures: keys = actual file stems; S10/S11 skipped ----
supp_titles <- c(
  "FigureS1"                      = "Figure S1. OXPHOS13 single-gene breakdown",
  "FigureS2"                      = "Figure S2. BN sensitivity analysis after excluding age",
  "FigureS3"                      = "Figure S3. SEM mediation models",
  "FigureS4"                      = "Figure S4. MR forest: OXPHOS genes vs AMD risk (two outcomes)",
  "FigureS5"                      = "Figure S5. ML performance across cohorts",
  "FigureS6"                      = "Figure S6. Cell-type proportions (dAMD vs Healthy)",
  "FigureS7"                      = "Figure S7. TRADE panel-level perturbation at G3 vs G2",
  "FigureS8"                      = "Figure S8. BioTIP MCI scores",
  "FigureS9"                      = "Figure S9. Permutation test histograms",
  "FigureS12"                     = "Figure S12. NDUFB3 multi-tissue eQTL landscape (GTEx V10)",
  "FigureS13_G2G3_robustness"     = "Figure S13. G2-G3 signal robustness assessment",
  "FigureS14_OXPHOS_variance"     = "Figure S14. OXPHOS variance decomposition and complex-level direction",
  "FigureS15_eQTM_background"     = "Figure S15. eQTM background sensitivity",
  "FigureS16_NDUFB3_KO_RGC_ROD"   = "Figure S16. NDUFB3 virtual knockout in RGC and ROD",
  "FigureS17a_ATAC_cohen_d"       = "Figure S17a. ATAC-seq promoter accessibility (Cohen's d)",
  "FigureS17b_ATAC_vs_RNA"        = "Figure S17b. ATAC vs RNA consistency (tissue-level)",
  "FigureS18"                     = "Figure S18. CellChat RPE communications",
  "FigureS19_greenyellow_ssGSEA"  = "Figure S19. greenyellow module ssGSEA by AMD grade",
  "FigureS20_CellChat_RPE"        = "Figure S20. CellChat RPE vs Rod/Cone pathway usage",
  "FigureS21_NDUFB3_KO_RPE"       = "Figure S21. NDUFB3 virtual knockout in RPE",
  "FigureS22_CpG_island_and_motif_hit" = "Figure S22. CpG island distribution and descriptive motif hit rates",
  "FigureS23_greenyellow_motif_enrichment" = "Figure S23. Motif enrichment of greenyellow-closed promoters (monaLisa)"
)

supp_figs <- intersect(paste0(names(supp_titles), ".pdf"),
                       list.files(fig_supp_dir, pattern = "^FigureS.*\\.pdf$"))
supp_figs <- supp_figs[order(as.integer(gsub("\\D", "", supp_figs)))]
log_info(sprintf("supplementary figures (kept): %d", length(supp_figs)))
log_info(sprintf("supplementary list: %s", paste(supp_figs, collapse = ", ")))

# ---- 5. Graphical Abstract ----
ga_fig <- file.path(fig_ga_dir, "GraphicalAbstract.pdf")
ga_exists <- file.exists(ga_fig)
log_info(sprintf("Graphical Abstract exists: %s", ga_exists))

# ---- Helper: PDF -> PNG -> scale to pixel aspect, centre on slide ----
add_pdf_image_slide <- function(doc, title, pdf_path, slide_w = 10.5, slide_h = 7.4) {
  if (has_magick) {
    tmp_png <- tempfile(fileext = ".png")
    pdftools::pdf_convert(pdf_path, format = "png", dpi = 150,
                          filenames = tmp_png, verbose = FALSE)
    info <- magick::image_info(magick::image_read(tmp_png))
    aspect <- info$width / info$height
    src <- tmp_png
  } else {
    ps <- tryCatch(pdftools::pdf_pages_size(pdf_path), error = function(e) NULL)
    if (is.null(ps) || nrow(ps) == 0 || is.na(ps$width[1]) || is.na(ps$height[1])) {
      log_warn(sprintf("cannot get size of %s, fallback aspect 1.414", basename(pdf_path)))
      aspect <- 1.414
    } else {
      aspect <- ps$width[1] / ps$height[1]
    }
    src <- pdf_path
  }

  if (aspect > slide_w / slide_h) {
    w <- slide_w; h <- slide_w / aspect
  } else {
    h <- slide_h; w <- slide_h * aspect
  }
  left <- 0.5 + (slide_w - w) / 2
  top  <- 1.2 + (slide_h - h) / 2

  add_slide(doc, layout = "Title and Content", master = "Office Theme") %>%
    ph_with(value = title, location = ph_location_type(type = "title")) %>%
    ph_with(
      value = external_img(src, width = w, height = h),
      location = ph_location(left = left, top = top, width = w, height = h)
    )
}

# ---- 6. Build pptx ----
log_info("--- building pptx ---")
doc <- read_pptx()

# ---- Slide 1: Title slide (no journal name) ----
doc <- add_slide(doc, layout = "Title Slide", master = "Office Theme") %>%
  ph_with(
    value = fpar(ftext(
      "A multi-scale causal triangulation framework validated by positive controls distinguishes causal drivers from concomitant features in age-related macular degeneration",
      fp_text(font.size = 24, bold = TRUE)
    )),
    location = ph_location_type(type = "ctrTitle")
  ) %>%
  ph_with(
    value = fpar(ftext(
      "Jiayi Hong, Jiafei Xi, Liqiang Wang",
      fp_text(font.size = 16)
    )),
    location = ph_location_type(type = "subTitle")
  )
log_info("added slide: Title")

# ---- Slide 2: Graphical Abstract (fixed 9.4 x 7.3 in) ----
if (ga_exists) {
  if (has_magick) {
    tmp_ga <- tempfile(fileext = ".png")
    pdftools::pdf_convert(ga_fig, format = "png", dpi = 150,
                          filenames = tmp_ga, verbose = FALSE)
    src_ga <- tmp_ga
  } else {
    src_ga <- ga_fig
  }
  doc <- doc %>%
    add_slide(layout = "Title and Content", master = "Office Theme") %>%
    ph_with(value = "Graphical Abstract", location = ph_location_type(type = "title")) %>%
    ph_with(
      value = external_img(src_ga, width = 9.4, height = 7.3),
      location = ph_location(left = 0.3, top = 1.1, width = 9.4, height = 7.3)
    )
  log_info("added slide: Graphical Abstract")
} else {
  log_warn("Graphical Abstract PDF not found, skipped")
}

# ---- Slide 3: Abstract (Purpose / Results / Conclusions) ----
abstract_text <- paste0(
  "Purpose: ", abstract_pur, "\n\n",
  "Results: ", abstract_res, "\n\n",
  "Conclusions: ", abstract_con
)
doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme") %>%
  ph_with(value = "Abstract", location = ph_location_type(type = "title")) %>%
  ph_with(
    value = fpar(ftext(abstract_text, fp_text(font.size = 11))),
    location = ph_location_type(type = "body")
  )
log_info("added slide: Abstract")

# ---- Slides 4-10: Figure 1-7 ----
for (nm in names(main_figs)) {
  f_pdf <- file.path(fig_main_dir, paste0(nm, ".pdf"))
  if (!file.exists(f_pdf)) {
    log_warn(sprintf("missing main figure: %s", f_pdf))
    next
  }
  doc <- add_pdf_image_slide(doc, main_figs[[nm]], f_pdf)
  log_info(sprintf("added slide: %s", nm))
}

# ---- Slide 11: Key conclusions ----
conclusions_text <- paste0(
  "1. OXPHOS downregulation is a reproducible correlate, not a strong causal driver.\n\n",
  "2. Two-axis regulation: genetic (16-gene, eGene-enriched) and epigenetic (greenyellow, eQTM-enriched) layers.\n\n",
  "3. NDUFB3 is the only gene-level candidate (early AMD, instrument-fragile, provisional).\n\n",
  "4. The framework recovered established causal signals for C3 (all-stage) and CFH (early AMD); CD59 was non-significant; no OXPHOS gene reached unified FDR except NDUFB3 in early AMD (instrument-fragile)."
)
doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme") %>%
  ph_with(value = "Key conclusions", location = ph_location_type(type = "title")) %>%
  ph_with(
    value = fpar(ftext(conclusions_text, fp_text(font.size = 16))),
    location = ph_location_type(type = "body")
  )
log_info("added slide: Key conclusions")

# ---- Slide 12+: Supplementary figures (descriptive titles; S10/S11 skipped) ----
for (i in seq_along(supp_figs)) {
  f_pdf <- file.path(fig_supp_dir, supp_figs[i])
  key <- gsub("\\.pdf$", "", supp_figs[i])
  title <- supp_titles[[key]]
  if (is.null(title)) title <- key
  doc <- add_pdf_image_slide(doc, title, f_pdf)
  log_info(sprintf("added slide: %s", key))
}

# ---- Save ----
out_pptx <- file.path(out_dir, "11_1A_slides.pptx")
print(doc, target = out_pptx)
log_info(sprintf("saved pptx: %s", out_pptx))

# ---- Summary ----
log_info(sprintf("total slides: %d", length(doc)))
log_info("=== 11_1A done (v2.5) ===")