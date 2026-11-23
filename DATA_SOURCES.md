# Data Sources

Complete inventory of datasets, versions, accessions, and code
dependencies for the AMD multi-omics analysis.

## Part 1: raw data (11 datasets)

| Dataset | Type | Accession | Used by | Figures/Tables |
|---|---|---|---|---|
| GSE115828 | Bulk RNA-seq (n=453) | GEO | 01_1_prepare_bulk_data.r, 01_1b_batch_enhance.r, 07_2_GSE115828_DE.r | Table1, Fig2/3, S2/S3 |
| GSE221042 | scRNA-seq (10x) | GEO | 06_1_singlecell.r, 07_3_* (11 scripts), 06_3_cellchat.r | Fig6D, Fig7, S6/S11/S16/S21 |
| GSE231536 | Methylation and eQTM | GEO | 07_5_GSE231536_methylation.r, 07_5b, 07_5c, 07_10, 07_10c | Fig4B/C, S5/S15/S22 |
| GSE265801 | snRNA-seq reference | GEO | 03_2_deconvolution.r, 04_3_TF_activity.r | Table2, S11 |
| GSE29801 | Bulk microarray (n=62) | GEO | 07_2_gse29801_replication.r | S1, S12 |
| GSE99287 | ATAC-seq (RPE and Retina) | GEO | 07_1_atac_validation.r | Fig4D, S16/S17 |
| GSE135092 | Orozco eye-eQTL | GEO | 07_8a GSE135092 prep.r | Table4, S5 |
| GSE234047 | Spatial transcriptomics | GEO | 07_9_GSE234047_spatial_validation.r | Spatial S figures |
| GSE262151 | RPE snRNA-seq | GEO | 07_3_scTenifoldKnk_NDUFB3_RPE.r | Fig7, S21 |
| PXD029501 | Proteomics | ProteomeXchange | 07_3_proteomics_opposite.r | 07_3_proteomics_opposite |
| PXD063075 | Proteomics | ProteomeXchange | 07_3_proteomics_opposite.r | 07_3_proteomics_opposite |

Download URLs:

- GEO datasets: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=<GSE_ID>
- ProteomeXchange: https://www.ebi.ac.uk/pride/archive/projects/<PXD_ID>

## Part 2: external summary statistics

| Dataset | Version | Accession | Used by | Figures/Tables |
|---|---|---|---|---|
| eQTLGen cis-eQTL | 2019-12-11 | eQTLGen | 07_7_extract_eqtlgen_cis.r, 07_7b_MR_plink.r | Fig5A, Table3 |
| GTEx eQTL | v8 (98 tissues) | GTEx Portal | 08_5_coloc.r | S14 |
| GTEx eQTL | v10 (API) | GTEx Portal API v2 | 07_7c_GTEx_MR_NDUFB3.r | Fig5A/C |
| GTEx eQTL | v6p | COLSTATS cs_4612_c / cs_1067364_c | 07_7f_coloc_NDUFB3.r | Table3, Fig5D |
| FinnGen AMD GWAS | finn-b-H7_AMD | FinnGen | 07_7b_MR_plink.r, 08_3_mr.r | Fig1B, Fig3C |
| IAMDGC 2018 early AMD | ebi-a-GCST010723 | GWAS Catalog | 07_7b, 07_7c, 07_7f, 08_3, 08_5 | Fig5D, Table3 |
| 1000 Genomes EUR | phase 3 | 1000G | 07_7b_MR_plink.r, 08_3_mr.r | LD clumping |
| Orozco eye-eQTL | GSE135092 | GEO | 07_7g_retinal_eqtl_screening_NDUFB3.r | Table4 |
| Strunz retinal eQTL | pgen.1008934 | PLoS Genet | 07_7g_retinal_eqtl_screening_NDUFB3.r | Table4 |
| EyeGEx | NIHMS1518564 | GEO | 07_7g_retinal_eqtl_screening_NDUFB3.r | Table4 |
| Wang 2023 sc-eQTL | 13059_2023_3111 | Genome Biol | 07_7g_retinal_eqtl_screening_NDUFB3.r | Table4 |

Download URLs:

- eQTLGen: https://www.eqtlgen.org/cis-eqtls.html
- GTEx v8: https://gtexportal.org/home/downloads/adult-gtex/qtl
- GTEx v10 API: https://gtexportal.org/api/v2/
- FinnGen: https://r9.finngen.fi/pheno/H7_AMD
- IAMDGC GWAS: https://www.ebi.ac.uk/gwas/studies/GCST010723
- 1000 Genomes EUR: https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eur.zip
- Strunz: https://doi.org/10.1371/journal.pgen.1008934
- Wang 2023: https://doi.org/10.1186/s13059-023-03111-8

## Part 3: external resources (Bioconductor and online platforms)

| Resource | Version | Access | Used by | Figures/Tables |
|---|---|---|---|---|
| JASPAR2024 CORE | 2024 | Bioconductor | 07_10b_motif_enrichment.r, 07_1b_atac_motif_enrichment.r | S22/S23 |
| COLSTATS | v2.3 | https://colstats.irgb.cnr.it | 07_7f_coloc_NDUFB3.r (manual upload) | Table3, Fig5D |
| MitoCarta | 3.0 | Broad Institute | 01_2_prepare_gene_sets.r | OXPHOS gene sets |
| MSigDB | v7.x | Broad Institute | 01_2_prepare_gene_sets.r | Pathway scores |
| CollecTRI | v2.x | decoupleR | 04_3_TF_activity.r | TF activity |

Download URLs:

- JASPAR2024: BiocManager::install("JASPAR2024")
- MitoCarta 3.0: https://personal.broadinstitute.org/scalvo/MitoCarta3.0/Human.MitoCarta3.0.xls
- MSigDB: https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp

## Colocalization workflow (Table 3 and Figure 5D)

1. R: 07_7f_coloc_NDUFB3.r produces ndufb3_gwas_region.csv
   (6,787 SNPs, GRCh37).
2. Python: 07_7f_coloc_NDUFB3.py produces
   NDUFB3_earlyAMD_COLSTATS_hg38.tsv (GRCh38).
3. COLSTATS v2.3 web platform (targeted mode; GTEx v6p):
   - Adipose_Subcutaneous (cs_4612_c, n = 581):
     raw output saved as NDUFB3_colstats.txt (session 2026-10-11).
   - Whole_Blood (cs_1067364_c, n = 670):
     raw web output not saved as a separate file; results transcribed
     from session 2026-10-12, with screenshots in figure/.
4. R: 07_7f2_coloc_NDUFB3.r produces coloc_NDUFB3_both_tissues.csv
   (tidy table).

Files retained in reproducibility/07_7f_coloc_NDUFB3/ for verification:

- NDUFB3_colstats.txt (COLSTATS raw output, adipose)
- NDUFB3_earlyAMD_COLSTATS_hg38.tsv (upload input, GRCh38)
- NDUFB3_earlyAMD_COLSTATS_input.tsv (upload input, intermediate)
- figure/*.png (web-session screenshots, both tissues)

## Data availability statement

All raw data analysed in this study are publicly available. Bulk retinal
transcriptomes (GSE115828), single-cell retinal transcriptomes (GSE221042,
GSE262151), retinal DNA methylation and eQTM data (GSE231536), snRNA-seq
deconvolution reference (GSE265801), external bulk replication (GSE29801),
retinal ATAC-seq (GSE99287), Orozco eye-eQTL (GSE135092), spatial
transcriptomics (GSE234047), and proteomics data (PXD029501, PXD063075)
were obtained from GEO and ProteomeXchange. Summary statistics for early
AMD (ebi-a-GCST010723, IAMDGC 2018) and all-stage AMD (finn-b-H7_AMD,
FinnGen) were downloaded from the GWAS Catalog and FinnGen, respectively.
eQTLGen cis-eQTL (2019-12-11 release), GTEx v8 eQTL and GTEx v10 eQTL
(via the GTEx Portal API v2) were used for Mendelian randomization and
sensitivity analyses. GTEx v6p NDUFB3 cis-eQTL and the COLSTATS v2.3
platform (https://colstats.irgb.cnr.it) were used for colocalization
(Table 3, Figure 5D). The COLSTATS raw output for adipose tissue, the
GRCh38 upload input, and web-session screenshots for both adipose and
whole blood are retained in results/07_7f_coloc_NDUFB3/ for verification.
LD reference panels were from the 1000 Genomes Project phase 3 EUR
population. JASPAR2024 CORE motifs were obtained via the JASPAR2024
Bioconductor package. A complete inventory of datasets, versions, and
accession numbers is provided in DATA_SOURCES.md.
