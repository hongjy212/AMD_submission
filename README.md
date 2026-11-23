# AMD multi-omics analysis — code repository

Analysis code for the manuscript "Retinal Oxidative Phosphorylation Gene Expression Across Age-Related Macular Degeneration Grades and Genetic Prioritisation of NDUFB3 in Early Disease".

## Contents

- `scripts/` — 97 analysis scripts (96 R + 1 Python)
  - `00_setup.r` — entry point (paths, parameters, I/O helpers)
  - `01_*` to `11_*` — analysis stages (see Run order)
  - `diag/` — diagnostic and robustness scripts
- `reproducibility/07_7f_coloc_NDUFB3/` — COLSTATS input, output and screenshots
- `DATA_SOURCES.md` — complete data inventory (accessions, versions, URLs)
- `Supplementary_Code_1.txt` — COLSTATS v2.3 colocalization workflow

## Requirements

- R >= 4.6.1
- Python >= 3.9 (packages: pandas, liftover)
- R packages: see sessionInfo() after running 00_setup.r

## Entry point

    source("scripts/00_setup.r")

00_setup.r defines:

- DIR_ROOT, DIR_RESULT, DIR_RAW
- SEED = 20260925, N_PERM = 1000
- COVARS = age, sex, batch, RIN, PMI, PC1-PC5
- helpers: setup_output_dir, save_rds, save_table, save_fig, read_upstream

## Run order

    01_*  data preparation          (5 scripts)
    02_*  bulk discovery            (15 scripts)
    03_*  deconvolution/mediation   (4 scripts)
    04_*  aging/pathway scoring     (6 scripts)
    05_*  WGCNA modules             (4 scripts)
    06_*  single-cell/TRADE/BioTIP  (5 scripts)
    07_*  validation layer          (41 scripts)
    08_*  full OXPHOS MR/coloc      (5 scripts)
    10_*  figure/table generation   (5 scripts)
    11_*  presentation              (1 script)
    diag/ diagnostic scripts        (5 scripts)

## Data

All raw data are public. See DATA_SOURCES.md for:

- GEO accessions: GSE115828, GSE221042, GSE231536, GSE265801,
  GSE29801, GSE99287, GSE135092, GSE234047, GSE262151
- ProteomeXchange: PXD029501, PXD063075
- GWAS: ebi-a-GCST010723 (GWAS Catalog), finn-b-H7_AMD (FinnGen)
- eQTL: eQTLGen (2019-12-11), GTEx v6p, GTEx v8, GTEx v10
- LD reference: 1000 Genomes phase 3 EUR
- Motifs: JASPAR2024 CORE (Bioconductor)

External data paths (local, not included):

- E:/retina_tools/raw/ — 11 GEO/ProteomeXchange datasets
- E:/retina_tools/pubdata/eQTL/ — eQTLGen, GTEx, GWAS, LD reference

## Colocalization (Table 3 / Figure 5D)

Colocalization was performed on the COLSTATS v2.3 web platform
(https://colstats.irgb.cnr.it), which is not scriptable. The complete
workflow is documented in Supplementary_Code_1.txt.

Input, raw output and web-session screenshots are retained in
reproducibility/07_7f_coloc_NDUFB3/:

- NDUFB3_earlyAMD_COLSTATS_hg38.tsv — GWAS input (GRCh38)
- NDUFB3_earlyAMD_COLSTATS_input.tsv — intermediate upload
- NDUFB3_colstats.txt — COLSTATS raw output (adipose)
- figure/*.png — web-session screenshots (both tissues)

The tidy result table is produced by scripts/07_7f2_coloc_NDUFB3.r.

## Excluded from this repository

- AMD260626 to AMD261021 directories — historical snapshots
  (17 directories; early mw* naming, later migrated to 01_* to 08_*)
- -02_5_causal_mediation_analysis.r — deprecated
  (superseded by 03_3_causal_mediation.r)
- pubdata/script/ — 7 unreferenced legacy R files
- pubdata/eQTL/sc_eQTL/single_cell_QTL-main/ — third-party code
  accompanying Wang et al. 2023, not used by this study

## Reproducibility

All analyses use a fixed seed (SEED = 20260925) and fixed permutation
count (N_PERM = 1000). Paths in 00_setup.r assume a Windows layout
(E:/retina_tools/...); adjust DIR_ROOT, DIR_RESULT, and DIR_RAW when
running on another system.

## License

GPL-3.0.
