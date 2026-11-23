# ============================================================
# 01_2_prepare_gene_sets.r
#
#
#   data/gene_sets.rds
#   tables/gene_sets_metadata.csv
#   tables/gene_sets_full.csv
# ============================================================

source("E:/retina_tools/aging_rgc_axon/AMD/scripts/00_setup.r")
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(readxl)
})

DIR_THIS <- setup_output_dir("01_2_prepare_gene_sets")
log_info("=== 01_2_prepare_gene_sets.r start ===")

MSIGDB_DIR      <- "E:/retina_tools/pubdata/msigdb"
MITOCARTA_FILE  <- "E:/retina_tools/pubdata/MitoCarta/Human.MitoCarta3.0.xls"
TF_LIST_DIR     <- "E:/retina_tools/pubdata/TF_list"
LAMBERT_FILE    <- file.path(TF_LIST_DIR, "Lambert2018_TF_names_v_1.01.txt")
HPA_RETINA_TSV  <- file.path(TF_LIST_DIR, "retina.tsv")
AMD_ANCHOR_FILE <- "E:/retina_tools/pubdata/eQTL/GSE135092/mmc5.xlsx"

MIN_GENES_GLOBAL <- 5
MAX_GENES_GLOBAL <- 500      # ★ 从 150 提高到 500

# ============================================================
# ============================================================
log_info("--- 1. 来源注册表 ---")

SOURCE_REGISTRY <- list(
  IEG = list(group = "artifact",
    primary = list(type = "IEG_TFS", id = "scCustomize_intersect_Lambert"),
    reference = "Tyssowski 2018; Lambert 2018", min_genes = 20),
  IFN_I = list(group = "artifact",
    primary = list(type = "Hallmark", id = "HALLMARK_INTERFERON_ALPHA_RESPONSE"),
    reference = "Liberzon 2015", min_genes = 5),
  IFN_II = list(group = "artifact",
    primary = list(type = "Hallmark", id = "HALLMARK_INTERFERON_GAMMA_RESPONSE"),
    reference = "Liberzon 2015", min_genes = 5),
  Inflammatory = list(group = "artifact",
    primary = list(type = "Hallmark", id = "HALLMARK_INFLAMMATORY_RESPONSE"),
    reference = "Liberzon 2015", min_genes = 5),
  Oxidative_stress = list(group = "artifact",
    primary = list(type = "Hallmark", id = "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY"),
    reference = "Liberzon 2015", min_genes = 5),
  Hypoxia = list(group = "artifact",
    primary = list(type = "Hallmark", id = "HALLMARK_HYPOXIA"),
    reference = "Liberzon 2015", min_genes = 5),
  ER_stress = list(group = "artifact",
    primary = list(type = "Hallmark", id = "HALLMARK_UNFOLDED_PROTEIN_RESPONSE"),
    reference = "Liberzon 2015", min_genes = 5),
  Apoptosis = list(group = "artifact",
    primary = list(type = "Hallmark", id = "HALLMARK_APOPTOSIS"),
    reference = "Liberzon 2015", min_genes = 5),
  TNFA_NFKB = list(group = "artifact",
    primary = list(type = "Hallmark", id = "HALLMARK_TNFA_SIGNALING_VIA_NFKB"),
    reference = "Liberzon 2015", min_genes = 5),

  OXPHOS = list(group = "metabolism",
    primary = list(type = "Hallmark", id = "HALLMARK_OXIDATIVE_PHOSPHORYLATION"),
    reference = "Liberzon 2015", min_genes = 10),
  Glycolysis = list(group = "metabolism",
    primary = list(type = "Hallmark", id = "HALLMARK_GLYCOLYSIS"),
    reference = "Liberzon 2015", min_genes = 5),
  Fatty_acid_oxidation = list(group = "metabolism",
    primary = list(type = "Hallmark", id = "HALLMARK_FATTY_ACID_METABOLISM"),
    reference = "Liberzon 2015", min_genes = 5),
  Mito_biogenesis = list(group = "metabolism",
    primary = list(type = "MitoCarta", id = "Mitochondrial central dogma"),
    reference = "Rath 2021", min_genes = 5),
  Mito_dynamics = list(group = "metabolism",
    primary = list(type = "MitoCarta", id = "Fusion"),
    extra   = list(type = "MitoCarta", id = "Fission"),
    reference = "Rath 2021", min_genes = 5),
  Mito_mitophagy = list(group = "metabolism",
    primary = list(type = "MitoCarta", id = "Mitophagy"),
    reference = "Rath 2021", min_genes = 3),
  Mito_mtDNA_rep = list(group = "metabolism",
    primary = list(type = "MitoCarta", id = "mtDNA replication"),
    reference = "Rath 2021", min_genes = 3),
  Mito_antioxidant = list(group = "metabolism",
    primary = list(type = "MitoCarta", id = "ROS and glutathione metabolism"),
    reference = "Rath 2021", min_genes = 3),

  # ---- Group C: panel ----
  Notch = list(group = "panel",
    primary = list(type = "KEGG", id = "KEGG_NOTCH_SIGNALING_PATHWAY"),
    fallback = list(type = "Hallmark", id = "HALLMARK_NOTCH_SIGNALING"),
    reference = "Kanehisa 2000", min_genes = 5),
  Complement = list(group = "panel",
    primary = list(type = "KEGG", id = "KEGG_COMPLEMENT_AND_COAGULATION_CASCADES"),
    reference = "Kanehisa 2000", min_genes = 5),
  Spliceosome = list(group = "panel",
    primary = list(type = "KEGG", id = "KEGG_SPLICEOSOME"),
    reference = "Kanehisa 2000", min_genes = 5),
  Proteasome = list(group = "panel",
    primary = list(type = "KEGG", id = "KEGG_PROTEASOME"),
    reference = "Kanehisa 2000", min_genes = 5),
  Autophagy = list(group = "panel",
    primary = list(type = "KEGG", id = "KEGG_REGULATION_OF_AUTOPHAGY"),
    reference = "Kanehisa 2000", min_genes = 5),
  Phototransduction = list(group = "panel",
    primary = list(type = "GOBP", id = "GOBP_PHOTOTRANSDUCTION"),
    reference = "Liberzon 2015", min_genes = 5),
  Photoreceptor_OuterSegment = list(group = "panel",
    primary = list(type = "GOCC", id = "GOCC_PHOTORECEPTOR_OUTER_SEGMENT"),
    reference = "GO Consortium", min_genes = 5),
  Photoreceptor_InnerSegment = list(group = "panel",
    primary = list(type = "GOCC", id = "GOCC_PHOTORECEPTOR_INNER_SEGMENT"),
    reference = "GO Consortium", min_genes = 5),
  Photoreceptor_ConnectingCilium = list(group = "panel",
    primary = list(type = "GOCC", id = "GOCC_PHOTORECEPTOR_CONNECTING_CILIUM"),
    reference = "GO Consortium", min_genes = 5),
  AMPK = list(group = "panel",
    primary = list(type = "GOBP", id = "GOBP_CAMKK_AMPK_SIGNALING_CASCADE"),
    reference = "Liberzon 2015", min_genes = 5),
  mTORC1 = list(group = "panel",
    primary = list(type = "Hallmark", id = "HALLMARK_MTORC1_SIGNALING"),
    reference = "Liberzon 2015", min_genes = 5),
  p53 = list(group = "panel",
    primary = list(type = "Hallmark", id = "HALLMARK_P53_PATHWAY"),
    reference = "Liberzon 2015", min_genes = 5),
  DNA_repair = list(group = "panel",
    primary = list(type = "Hallmark", id = "HALLMARK_DNA_REPAIR"),
    reference = "Liberzon 2015", min_genes = 5),
  Apical_junction = list(group = "panel",
    primary = list(type = "Hallmark", id = "HALLMARK_APICAL_JUNCTION"),
    reference = "Liberzon 2015", min_genes = 5)
)

VISUAL_REGISTRY <- list(
  Phototransduction_KEGG = list(group = "visual",
    primary = list(type = "KEGG", id = "KEGG_PHOTOTRANSDUCTION"),
    alt_ids = c("KEGG_PHOTOTRANSDUCTION"),
    reference = "Kanehisa 2000", min_genes = 5),
  Phototransduction_Activation = list(group = "visual",
    primary = list(type = "Reactome",
                    id = "REACTOME_ACTIVATION_OF_THE_PHOTOTRANSDUCTION_CASCADE"),
    alt_ids = c("REACTOME_ACTIVATION_OF_THE_PHOTOTRANSDUCTION_CASCADE",
                "REACTOME_THE_ACTIVATION_OF_THE_PHOTOTRANSDUCTION_CASCADE"),
    reference = "Reactome", min_genes = 5),
  Phototransduction_Inactivation = list(group = "visual",
    primary = list(type = "Reactome",
                    id = "REACTOME_INACTIVATION_RECOVERY_AND_REGULATION_OF_THE_PHOTOTRANSDUCTION_CASCADE"),
    alt_ids = c("REACTOME_INACTIVATION_RECOVERY_AND_REGULATION_OF_THE_PHOTOTRANSDUCTION_CASCADE",
                "REACTOME_INACTIVATION_RECOVERY_AND_REGULATION_OF_THE_PHOTOTRANSDUCTION"),
    reference = "Reactome", min_genes = 5),

  Visual_Perception = list(group = "visual",
    primary = list(type = "GOBP", id = "GOBP_VISUAL_PERCEPTION"),
    alt_ids = c("GOBP_VISUAL_PERCEPTION",
                "GOBP_VISUAL_PERCEPTION_LIGHT_STIMULUS"),
    reference = "GO Consortium", min_genes = 5),
  Detection_of_Light = list(group = "visual",
    primary = list(type = "GOBP",
                    id = "GOBP_DETECTION_OF_LIGHT_STIMULUS_INVOLVED_IN_VISUAL_PERCEPTION"),
    alt_ids = c("GOBP_DETECTION_OF_LIGHT_STIMULUS_INVOLVED_IN_VISUAL_PERCEPTION",
                "GOBP_DETECTION_OF_LIGHT_STIMULUS",
                "GOBP_DETECTION_OF_VISIBLE_LIGHT"),
    reference = "GO Consortium", min_genes = 5),

  Photoreceptor_Differentiation = list(group = "visual",
    primary = list(type = "GOBP", id = "GOBP_PHOTORECEPTOR_CELL_DIFFERENTIATION"),
    reference = "GO Consortium", min_genes = 5),
  Photoreceptor_Maintenance = list(group = "visual",
    primary = list(type = "GOBP", id = "GOBP_PHOTORECEPTOR_CELL_MAINTENANCE"),
    reference = "GO Consortium", min_genes = 5),

  Retina_Development = list(group = "visual",
    primary = list(type = "GOBP", id = "GOBP_RETINA_DEVELOPMENT"),
    alt_ids = c("GOBP_RETINA_DEVELOPMENT",
                "GOBP_RETINA_MORPHOGENESIS",
                "GOBP_RETINA_DEVELOPMENT_IN_CAMERA_TYPE_EYE"),
    reference = "GO Consortium", min_genes = 5),
  Retina_Homeostasis = list(group = "visual",
    primary = list(type = "GOBP", id = "GOBP_RETINA_HOMEOSTASIS"),
    reference = "GO Consortium", min_genes = 5),

  Retinoid_Cycle = list(group = "visual",
    primary = list(type = "Reactome", id = "REACTOME_THE_RETINOID_CYCLE"),
    alt_ids = c("REACTOME_THE_RETINOID_CYCLE",
                "REACTOME_RETINOID_CYCLE",
                "REACTOME_VISUAL_PHOTOTRANSDUCTION"),
    reference = "Reactome", min_genes = 3),

  Sensory_Cilium_Assembly = list(group = "visual",
    primary = list(type = "GOBP", id = "GOBP_SENSORY_CILIUM_ASSEMBLY"),
    alt_ids = c("GOBP_SENSORY_CILIUM_ASSEMBLY",
                "GOBP_SENSORY_CILIUM_ORGANIZATION"),
    reference = "GO Consortium", min_genes = 5),
  Cilium_Assembly = list(group = "visual",
    primary = list(type = "GOBP", id = "GOBP_CILIUM_ASSEMBLY"),
    alt_ids = c("GOBP_CILIUM_ASSEMBLY",
                "GOBP_CILIUM_ORGANIZATION",
                "GOBP_CILIUM_MORPHOGENESIS"),
    reference = "GO Consortium", min_genes = 10),

  Axonal_Transport = list(group = "axon",
    primary = list(type = "GOBP", id = "GOBP_AXONAL_TRANSPORT"),
    reference = "GO Consortium", min_genes = 10),
  Axon_Regeneration = list(group = "axon",
    primary = list(type = "GOBP", id = "GOBP_AXON_REGENERATION"),
    reference = "GO Consortium", min_genes = 5),
  Axon_Ensheathment = list(group = "axon",
    primary = list(type = "GOBP", id = "GOBP_AXON_ENSHEATHMENT"),
    reference = "GO Consortium", min_genes = 5),
  Axonal_Growth_Cone = list(group = "axon",
    primary = list(type = "GOCC", id = "GOCC_AXONAL_GROWTH_CONE"),
    reference = "GO Consortium", min_genes = 5),
  Axon_Development = list(group = "axon",
    primary = list(type = "GOBP", id = "GOBP_AXON_DEVELOPMENT"),
    reference = "GO Consortium", min_genes = 10),
  Axonogenesis = list(group = "axon",
    primary = list(type = "GOBP", id = "GOBP_AXONOGENESIS"),
    reference = "GO Consortium", min_genes = 10),
  Axon_Guidance_GOBP = list(group = "axon",
    primary = list(type = "GOBP", id = "GOBP_AXON_GUIDANCE"),
    reference = "GO Consortium", min_genes = 10),
  Growth_Cone = list(group = "axon",
    primary = list(type = "GOCC", id = "GOCC_GROWTH_CONE"),
    alt_ids = c("GOCC_GROWTH_CONE",
                "GOCC_GROWTH_CONE_MEMBRANE"),
    reference = "GO Consortium", min_genes = 5),
  Optic_Nerve_Development = list(group = "axon",
    primary = list(type = "GOBP", id = "GOBP_OPTIC_NERVE_DEVELOPMENT"),
    reference = "GO Consortium", min_genes = 3)
)

AMD_PATHOLOGY_REGISTRY <- list(
  Complement_Cascade_Reactome = list(group = "amd_pathology",
    primary = list(type = "Reactome", id = "REACTOME_COMPLEMENT_CASCADE"),
    reference = "Reactome", min_genes = 5),
  Complement_Activation_GOBP = list(group = "amd_pathology",
    primary = list(type = "GOBP", id = "GOBP_COMPLEMENT_ACTIVATION"),
    reference = "GO Consortium", min_genes = 5),
  Complement_Regulation = list(group = "amd_pathology",
    primary = list(type = "GOBP", id = "GOBP_REGULATION_OF_COMPLEMENT_ACTIVATION"),
    reference = "GO Consortium", min_genes = 5),
  Complement_Classical = list(group = "amd_pathology",
    primary = list(type = "GOBP",
                    id = "GOBP_COMPLEMENT_ACTIVATION_CLASSICAL_PATHWAY"),
    reference = "GO Consortium", min_genes = 3),
  Complement_Alternative = list(group = "amd_pathology",
    primary = list(type = "GOBP",
                    id = "GOBP_COMPLEMENT_ACTIVATION_ALTERNATIVE_PATHWAY"),
    reference = "GO Consortium", min_genes = 3),

  # ---- cGAS-STING ----
  cGAS_STING_STING_Mediated = list(group = "cgast_amp",
    primary = list(type = "Reactome",
                    id = "REACTOME_STING_MEDIATED_INDUCTION_OF_HOST_IMMUNE_RESPONSES"),
    reference = "Reactome", min_genes = 5),
  cGAS_STING_Cytosolic_DNA_Sensors = list(group = "cgast_amp",
    primary = list(type = "Reactome",
                    id = "REACTOME_CYTOSOLIC_SENSORS_OF_PATHOGEN_ASSOCIATED_DNA"),
    reference = "Reactome", min_genes = 5),
  cGAS_STING_GOBP = list(group = "cgast_amp",
    primary = list(type = "GOBP",
                    id = "GOBP_CGAS_STING_SIGNALING_PATHWAY"),
    alt_ids = c("GOBP_CGAS_STING_SIGNALING_PATHWAY",
                "GOBP_CGAS_STING"),
    reference = "GO Consortium", min_genes = 5),

  Inflammasome = list(group = "amd_pathology",
    primary = list(type = "Reactome", id = "REACTOME_INFLAMMASOMES"),
    reference = "Reactome", min_genes = 5),
  Innate_Immune = list(group = "amd_pathology",
    primary = list(type = "Reactome", id = "REACTOME_INNATE_IMMUNE_SYSTEM"),
    reference = "Reactome", min_genes = 10),

  Cellular_Senescence = list(group = "amd_pathology",
    primary = list(type = "Reactome", id = "REACTOME_CELLULAR_SENESCENCE"),
    reference = "Reactome", min_genes = 5),
  Senescence_Hallmark = list(group = "amd_pathology",
    primary = list(type = "Hallmark", id = "HALLMARK_CELLULAR_SENESCENCE"),
    alt_ids = c("HALLMARK_CELLULAR_SENESCENCE",
                "HALLMARK_SENESCENCE"),
    reference = "Liberzon 2015", min_genes = 10),
  Oxidative_Senescence = list(group = "amd_pathology",
    primary = list(type = "Reactome",
                    id = "REACTOME_OXIDATIVE_STRESS_INDUCED_SENESCENCE"),
    reference = "Reactome", min_genes = 3),

  ECM_Organization = list(group = "amd_pathology",
    primary = list(type = "Reactome",
                    id = "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION"),
    reference = "Reactome", min_genes = 10),
  Collagen_Formation = list(group = "amd_pathology",
    primary = list(type = "Reactome", id = "REACTOME_COLLAGEN_FORMATION"),
    reference = "Reactome", min_genes = 10),

  VEGF_Signaling = list(group = "amd_pathology",
    primary = list(type = "KEGG", id = "KEGG_VEGF_SIGNALING_PATHWAY"),
    reference = "Kanehisa 2000", min_genes = 10)
)

MITO_SUBFUNC_REGISTRY <- list(

  Mito_Complex_I = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "Complex I"),
    reference = "Rath 2021", min_genes = 10),
  Mito_Complex_II = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "Complex II"),
    reference = "Rath 2021", min_genes = 3),
  Mito_Complex_III = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "Complex III"),
    reference = "Rath 2021", min_genes = 3),
  Mito_Complex_IV = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "Complex IV"),
    reference = "Rath 2021", min_genes = 5),
  Mito_Complex_V = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "Complex V"),
    reference = "Rath 2021", min_genes = 5),

  Mito_CI_Subunits = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "CI subunits"),
    reference = "Rath 2021", min_genes = 10),
  Mito_CI_Assembly = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "CI assembly factors"),
    reference = "Rath 2021", min_genes = 5),
  Mito_CII_Subunits = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "CII subunits"),
    reference = "Rath 2021", min_genes = 3),
  Mito_CIII_Subunits = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "CIII subunits"),
    reference = "Rath 2021", min_genes = 3),
  Mito_CIII_Assembly = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "CIII assembly factors"),
    reference = "Rath 2021", min_genes = 3),
  Mito_CIV_Subunits = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "CIV subunits"),
    reference = "Rath 2021", min_genes = 5),
  Mito_CIV_Assembly = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "CIV assembly factors"),
    reference = "Rath 2021", min_genes = 5),
  Mito_CV_Subunits = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "CV subunits"),
    reference = "Rath 2021", min_genes = 5),
  Mito_CV_Assembly = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "CV assembly factors"),
    reference = "Rath 2021", min_genes = 3),
  Mito_Respirasome = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "Respirasome assembly"),
    reference = "Rath 2021", min_genes = 3),

  # ---- TCA ----
  Mito_TCA = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "TCA cycle"),
    reference = "Rath 2021", min_genes = 5),
  Mito_TCA_associated = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "TCA-associated"),
    reference = "Rath 2021", min_genes = 5),

  Mito_mtDNA_maintenance = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "mtDNA maintenance"),
    reference = "Rath 2021", min_genes = 5),
  Mito_mtDNA_replication = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "mtDNA replication"),
    reference = "Rath 2021", min_genes = 3),
  Mito_mtDNA_repair = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "mtDNA repair"),
    reference = "Rath 2021", min_genes = 3),
  Mito_mtRNA_metabolism = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "mtRNA metabolism"),
    reference = "Rath 2021", min_genes = 10),
  Mito_Transcription = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "Transcription"),
    reference = "Rath 2021", min_genes = 3),
  Mito_Translation = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "Translation"),
    reference = "Rath 2021", min_genes = 10),
  Mito_Ribosome = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "Mitochondrial ribosome"),
    reference = "Rath 2021", min_genes = 10),
  Mito_Ribosome_Assembly = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "Mitochondrial ribosome assembly"),
    reference = "Rath 2021", min_genes = 5),
  Mito_Translation_Factors = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "Translation factors"),
    reference = "Rath 2021", min_genes = 5),
  Mito_tRNA_synthetases = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "mt-tRNA synthetases"),
    reference = "Rath 2021", min_genes = 5),
  Mito_central_dogma = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "Mitochondrial central dogma"),
    reference = "Rath 2021", min_genes = 10),

  Mito_Proteostasis = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "Protein homeostasis"),
    reference = "Rath 2021", min_genes = 10),
  Mito_Proteases = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "Proteases"),
    reference = "Rath 2021", min_genes = 5),
  Mito_Chaperones = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "Chaperones"),
    reference = "Rath 2021", min_genes = 5),
  Mito_Import = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "Protein import and sorting"),
    reference = "Rath 2021", min_genes = 10),
  Mito_TIM23 = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "TIM23 presequence pathway"),
    reference = "Rath 2021", min_genes = 3),
  Mito_TIM22 = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "TIM22 carrier pathway"),
    reference = "Rath 2021", min_genes = 3),
  Mito_TOM = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "TOM"),
    reference = "Rath 2021", min_genes = 3),

  Mito_Fusion = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "Fusion"),
    reference = "Rath 2021", min_genes = 3),
  Mito_Fission = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "Fission"),
    reference = "Rath 2021", min_genes = 3),
  Mito_Mitophagy = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "Mitophagy"),
    reference = "Rath 2021", min_genes = 3),
  Mito_Autophagy = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "Autophagy"),
    reference = "Rath 2021", min_genes = 5),
  Mito_Apoptosis = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "Apoptosis"),
    reference = "Rath 2021", min_genes = 5),
  Mito_Cristae = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "Cristae formation"),
    reference = "Rath 2021", min_genes = 5),
  Mito_Dynamics_Surveillance = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "Mitochondrial dynamics and surveillance"),
    reference = "Rath 2021", min_genes = 10),

  Mito_Carriers = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "SLC25A family"),
    reference = "Rath 2021", min_genes = 10),
  Mito_SmallMoleculeTransport = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "Small molecule transport"),
    reference = "Rath 2021", min_genes = 10),
  Mito_Calcium = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "Calcium homeostasis"),
    reference = "Rath 2021", min_genes = 5),
  Mito_CalciumCycle = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "Calcium cycle"),
    reference = "Rath 2021", min_genes = 3),
  Mito_mPTP = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "Mitochondrial permeability transition pore"),
    reference = "Rath 2021", min_genes = 3),

  Mito_FattyAcidOxidation = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "Fatty acid oxidation"),
    reference = "Rath 2021", min_genes = 10),
  Mito_Cardiolipin = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "Cardiolipin synthesis"),
    reference = "Rath 2021", min_genes = 5),
  Mito_ROS_Glutathione = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "ROS and glutathione metabolism"),
    reference = "Rath 2021", min_genes = 5),
  Mito_FeS_Cluster = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "Fe-S cluster biosynthesis"),
    reference = "Rath 2021", min_genes = 5),
  Mito_Heme = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "Heme synthesis and processing"),
    reference = "Rath 2021", min_genes = 5),
  Mito_CoQ = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "Coenzyme Q metabolism"),
    reference = "Rath 2021", min_genes = 5),

  Mito_Immune = list(group = "mito_subfunc",
    primary = list(type = "MitoCarta", id = "Immune response"),
    reference = "Rath 2021", min_genes = 3)
)

SOURCE_REGISTRY <- c(SOURCE_REGISTRY, VISUAL_REGISTRY,
                      AMD_PATHOLOGY_REGISTRY, MITO_SUBFUNC_REGISTRY)
log_info("注册基因集（合并后）: {length(SOURCE_REGISTRY)} 个")

# ============================================================
# ============================================================
log_info("--- 2. 加载外部数据库 ---")

# 2a. MSigDB GMT
find_gmt <- function(prefix) {
  patterns <- c(
    paste0("^", prefix, "\\.all\\.v.*symbols\\.gmt$"),
    paste0("^", prefix, "\\..*symbols\\.gmt$"),
    paste0("^", prefix, ".*\\.gmt$"))
  for (pat in patterns) {
    f <- list.files(MSIGDB_DIR, pattern = pat, full.names = TRUE,
                    ignore.case = TRUE)
    if (length(f) > 0) return(f[1])
  }
  NULL
}
gmt_files <- list(
  hallmark = find_gmt("h\\.all"),
  c5       = find_gmt("c5"),
  c2       = find_gmt("c2"))
for (nm in names(gmt_files)) {
  if (!is.null(gmt_files[[nm]])) log_info("  {nm}: {basename(gmt_files[[nm]])}")
}

read_gmt_safe <- function(f) {
  if (is.null(f) || !file.exists(f)) return(NULL)
  tryCatch(clusterProfiler::read.gmt(f),
           error = function(e) { log_warn("读 GMT 失败: {f}"); NULL })
}

log_info("加载 MSigDB...")
hallmark_df <- read_gmt_safe(gmt_files$hallmark)
c2_df       <- read_gmt_safe(gmt_files$c2)
c5_df       <- read_gmt_safe(gmt_files$c5)

msigdb <- list()
if (!is.null(hallmark_df)) msigdb$hallmark <- hallmark_df
if (!is.null(c2_df)) {
  msigdb$kegg     <- c2_df[grepl("^KEGG_",     c2_df$term), ]
  msigdb$reactome <- c2_df[grepl("^REACTOME_", c2_df$term), ]
}
if (!is.null(c5_df)) {
  msigdb$gobp <- c5_df[grepl("^GOBP_", c5_df$term), ]
  msigdb$gocc <- c5_df[grepl("^GOCC_", c5_df$term), ]
}
for (nm in names(msigdb)) {
  log_info("  {nm}: {length(unique(msigdb[[nm]]$term))} terms")
}

extract_gmt <- function(df, term_id) {
  if (is.null(df)) return(NULL)
  g <- df$gene[df$term == term_id]
  if (length(g) == 0) return(NULL)
  unique(g)
}

extract_gmt_multi <- function(df, ids) {
  if (is.null(df) || length(ids) == 0) return(NULL)
  for (id in ids) {
    g <- extract_gmt(df, id)
    if (!is.null(g) && length(g) > 0) {
      attr(g, "matched_id") <- id
      return(g)
    }
  }
  NULL
}

extract_gmt_version <- function(f) {
  if (is.null(f)) return(NA_character_)
  m <- regmatches(basename(f), regexpr("v[0-9]+\\.[0-9]+", basename(f)))
  if (length(m) == 0) return(NA_character_)
  m
}
gmt_versions <- sapply(gmt_files, extract_gmt_version)

# 2b. MitoCarta
log_info("加载 MitoCarta")
mito_pathways <- NULL
if (file.exists(MITOCARTA_FILE)) {
  sheets <- readxl::excel_sheets(MITOCARTA_FILE)
  target <- sheets[grepl("MitoPathway", sheets, ignore.case = TRUE)][1]
  if (is.na(target)) target <- sheets[1]
  mito_df <- readxl::read_excel(MITOCARTA_FILE, sheet = target)
  colnames(mito_df) <- trimws(colnames(mito_df))
  pc <- grep("^MitoPathway$", colnames(mito_df), value = TRUE)[1]
  hc <- grep("MitoPathways Hierarchy", colnames(mito_df), value = TRUE)[1]
  gc <- grep("^Genes$", colnames(mito_df), value = TRUE)[1]
  if (!is.na(pc) && !is.na(gc)) {
    mito_pathways <- data.frame(
      pathway = mito_df[[pc]],
      hierarchy = if (!is.na(hc)) mito_df[[hc]] else NA_character_,
      genes = mito_df[[gc]], stringsAsFactors = FALSE)
    log_info("  MitoCarta: {nrow(mito_pathways)} rows")
  }
}

extract_mito <- function(df, name) {
  if (is.null(df)) return(NULL)

  idx <- integer(0)

  if ("pathway" %in% colnames(df)) {
    idx <- which(df$pathway == name)
  }

  if (length(idx) == 0 && "hierarchy" %in% colnames(df) &&
      any(!is.na(df$hierarchy))) {
    idx <- which(df$hierarchy == name)
  }

  if (length(idx) == 0 && "hierarchy" %in% colnames(df) &&
      any(!is.na(df$hierarchy))) {
    idx <- which(endsWith(df$hierarchy, paste0(" > ", name)))
  }

  if (length(idx) == 0 && "pathway" %in% colnames(df)) {
    idx <- which(tolower(df$pathway) == tolower(name))
  }

  if (length(idx) == 0 && "hierarchy" %in% colnames(df) &&
      any(!is.na(df$hierarchy))) {
    idx <- which(grepl(name, df$hierarchy, fixed = TRUE))
  }

  if (length(idx) == 0) return(NULL)

  if (length(idx) > 1 && "hierarchy" %in% colnames(df)) {
    h <- df$hierarchy[idx]
    depth <- sapply(h, function(x) {
      if (is.na(x) || x == "") return(Inf)
      length(strsplit(x, " > ", fixed = TRUE)[[1]])
    })
    idx <- idx[which.min(depth)]
  }

  gs <- df$genes[idx[1]]
  if (is.na(gs)) return(NULL)
  g <- trimws(unlist(strsplit(gs, ",\\s*")))
  unique(g[g != ""])
}

# 2c. Lambert
log_info("加载 Lambert 2018")
TF_REFERENCE <- NULL
if (file.exists(LAMBERT_FILE)) {
  tf <- toupper(trimws(readLines(LAMBERT_FILE)))
  tf <- tf[tf != ""]
  TF_REFERENCE <- list(genes = unique(tf),
                       source = "Lambert 2018", version = "v1.01",
                       reference = "Lambert 2018")
  log_info("  Lambert TF: {length(TF_REFERENCE$genes)}")
}

# 2d. HPA retina
log_info("加载 HPA retina")
HPA_RETINA <- NULL
if (file.exists(HPA_RETINA_TSV)) {
  hpa <- read.delim(HPA_RETINA_TSV, stringsAsFactors = FALSE,
                    check.names = FALSE, quote = "")
  names(hpa) <- gsub('^"|"$', "", names(hpa))
  hpa$gene_clean <- toupper(trimws(gsub('^"|"$', "", hpa$Gene)))
  hpa$spec_clean <- trimws(gsub('^"|"$', "", hpa$`RNA tissue specificity`))
  HPA_RETINA <- hpa
  log_info("  HPA retina: {nrow(hpa)}")
}

# 2e. scCustomize IEG
log_info("加载 scCustomize IEG")
IEG_ALL <- NULL
if (requireNamespace("scCustomize", quietly = TRUE)) {
  ieg_all <- toupper(scCustomize::ieg_gene_list$Homo_sapiens_IEG)
  IEG_ALL <- list(genes = unique(ieg_all), source = "scCustomize",
                  version = as.character(packageVersion("scCustomize")),
                  reference = "Tyssowski 2018")
  log_info("  全 IEG: {length(IEG_ALL$genes)}")
}

IEG_TFS <- NULL
if (!is.null(IEG_ALL) && !is.null(TF_REFERENCE)) {
  ieg_tfs <- intersect(IEG_ALL$genes, TF_REFERENCE$genes)
  IEG_TFS <- list(genes = sort(ieg_tfs),
                  source = "scCustomize ∩ Lambert2018",
                  version = paste0("scCustomize v", IEG_ALL$version),
                  reference = "Tyssowski 2018; Lambert 2018")
  log_info("  IEG ∩ TF: {length(ieg_tfs)}")
}

# ============================================================
# 2g. AMD positiveanchor
# ============================================================
log_info("加载 AMD 阳性锚点")
AMD_TIER1 <- character(0); AMD_TIER2 <- character(0)

read_anchor <- function(f) {
  if (!file.exists(f)) return(NULL)
  raw <- readxl::read_excel(f, col_names = FALSE, .name_repair = "minimal")
  if (nrow(raw) < 2) return(NULL)
  row1 <- as.character(raw[1, ]); row2 <- as.character(raw[2, ])
  row1[is.na(row1)] <- ""; row2[is.na(row2)] <- ""
  last <- ""
  for (i in seq_along(row1)) {
    if (row1[i] != "") last <- row1[i]
    else if (row2[i] != "") row1[i] <- last
  }
  two_hdr <- any(grepl("Locus name", row1, ignore.case = TRUE)) &&
             any(grepl("Probability|eQTL symbol", row2, ignore.case = TRUE))
  if (two_hdr) {
    cn <- ifelse(row1 == "" & row2 == "", paste0("V", seq_along(row1)),
          ifelse(row1 == "", row2,
          ifelse(row2 == "" | row1 == row2, row1,
                 paste(row1, row2, sep = " | "))))
    dat <- raw[-c(1, 2), , drop = FALSE]
  } else {
    cn <- ifelse(row1 == "", paste0("V", seq_along(row1)), row1)
    dat <- raw[-1, , drop = FALSE]
  }
  colnames(dat) <- make.unique(trimws(gsub("\\s+", " ", cn)))
  dat
}
find_col <- function(nms, pat) {
  h <- grep(pat, nms, ignore.case = TRUE, value = TRUE)
  if (length(h) == 0) return(NA_character_)
  h[1]
}

if (file.exists(AMD_ANCHOR_FILE)) {
  dat <- read_anchor(AMD_ANCHOR_FILE)
  if (!is.null(dat)) {
    nms <- colnames(dat)
    c_gene    <- find_col(nms, "Putative.*Causal.*gene")
    c_gwas_p  <- find_col(nms, "^GWAS P-value$|GWAS.*P-value$")
    c_coloc   <- find_col(nms, "COLOC.*Probability|^Probability$")
    c_tissue  <- find_col(nms, "COLOC.*Tissue|^Tissue$")
    c_de_r    <- find_col(nms, "Control vs AMD in Retina.*Adjusted")
    c_de_rpe  <- find_col(nms, "Control vs AMD in RPE.*Adjusted")
    c_eqtl_rs <- find_col(nms, "Retina.*macula.*eQTL symbol")
    c_eqtl_rf <- find_col(nms, "Retina.*macula.*eQTL FDR")
    c_eqtl_ps <- find_col(nms, "RPE.*macula.*eQTL symbol")
    c_eqtl_pf <- find_col(nms, "RPE.*macula.*eQTL FDR")
    need <- c(c_gene, c_gwas_p, c_coloc, c_tissue,
              c_de_r, c_de_rpe, c_eqtl_rs, c_eqtl_rf, c_eqtl_ps, c_eqtl_pf)
    if (all(!is.na(need))) {
      for (cc in c(c_gwas_p, c_coloc, c_de_r, c_de_rpe,
                   c_eqtl_rf, c_eqtl_pf)) {
        dat[[cc]] <- suppressWarnings(as.numeric(dat[[cc]]))
      }
      tissue_chr <- as.character(dat[[c_tissue]])
      retina_ok <- !is.na(dat[[c_eqtl_rf]]) & dat[[c_eqtl_rf]] < 0.05 &
                   !is.na(dat[[c_de_r]])   & dat[[c_de_r]]   < 0.05
      rpe_ok    <- !is.na(dat[[c_eqtl_pf]]) & dat[[c_eqtl_pf]] < 0.05 &
                   !is.na(dat[[c_de_rpe]]) & dat[[c_de_rpe]] < 0.05
      macula_coloc <- !is.na(tissue_chr) &
                      grepl("macula", tissue_chr, ignore.case = TRUE) &
                      !grepl("non-macula", tissue_chr, ignore.case = TRUE)
      causal_gene <- ifelse(
        retina_ok, toupper(trimws(as.character(dat[[c_eqtl_rs]]))),
        ifelse(rpe_ok, toupper(trimws(as.character(dat[[c_eqtl_ps]]))),
               NA_character_))
      base_ok <- !is.na(causal_gene) & causal_gene != "" &
                 !is.na(dat[[c_gwas_p]]) & dat[[c_gwas_p]] < 5e-8 &
                 macula_coloc & (retina_ok | rpe_ok)
      t1 <- base_ok & !is.na(dat[[c_coloc]]) & dat[[c_coloc]] >= 0.8
      t2 <- base_ok & !is.na(dat[[c_coloc]]) & dat[[c_coloc]] >= 0.5 &
            dat[[c_coloc]] < 0.8
      AMD_TIER1 <- sort(unique(causal_gene[t1]))
      AMD_TIER2 <- sort(unique(causal_gene[t2]))
      log_info("  AMD Tier1: {length(AMD_TIER1)} 个")
      log_info("  AMD Tier2: {length(AMD_TIER2)} 个")
    } else {
      log_warn("AMD 锚点：列名匹配不全，跳过")
      miss <- names(need)[is.na(need)]
      log_warn("  未匹配: {paste(miss, collapse=', ')}")
    }
  }
}

# ============================================================
# ============================================================
log_info("--- 3. 组装 gene_sets ---")

GENE_SETS <- list(); META <- data.frame()
UNMATCHED <- data.frame()

for (gs in names(SOURCE_REGISTRY)) {
  spec <- SOURCE_REGISTRY[[gs]]
  genes <- NULL; src_used <- NA_character_; src_id <- NA_character_
  ver <- NA_character_; is_fb <- FALSE
  matched_id <- NA_character_

  p <- spec$primary
  candidates <- unique(c(p$id, spec$alt_ids))

  if (p$type == "IEG_TFS") {
    if (!is.null(IEG_TFS)) {
      genes <- IEG_TFS$genes; src_used <- IEG_TFS$source
      src_id <- "scCustomize_intersect_Lambert"; ver <- IEG_TFS$version
      matched_id <- src_id
    }
  } else if (p$type == "Hallmark") {
    g <- extract_gmt_multi(msigdb$hallmark, candidates)
    if (!is.null(g)) {
      genes <- g; src_used <- "MSigDB Hallmark"
      matched_id <- attr(g, "matched_id"); src_id <- matched_id
      ver <- gmt_versions[["hallmark"]]
    }
  } else if (p$type == "GOBP") {
    g <- extract_gmt_multi(msigdb$gobp, candidates)
    if (!is.null(g)) {
      genes <- g; src_used <- "MSigDB GOBP"
      matched_id <- attr(g, "matched_id"); src_id <- matched_id
      ver <- gmt_versions[["c5"]]
    }
  } else if (p$type == "GOCC") {
    g <- extract_gmt_multi(msigdb$gocc, candidates)
    if (!is.null(g)) {
      genes <- g; src_used <- "MSigDB GOCC"
      matched_id <- attr(g, "matched_id"); src_id <- matched_id
      ver <- gmt_versions[["c5"]]
    }
  } else if (p$type == "KEGG") {
    g <- extract_gmt_multi(msigdb$kegg, candidates)
    if (!is.null(g)) {
      genes <- g; src_used <- "MSigDB KEGG"
      matched_id <- attr(g, "matched_id"); src_id <- matched_id
      ver <- gmt_versions[["c2"]]
    }
  } else if (p$type == "Reactome") {
    g <- extract_gmt_multi(msigdb$reactome, candidates)
    if (!is.null(g)) {
      genes <- g; src_used <- "MSigDB Reactome"
      matched_id <- attr(g, "matched_id"); src_id <- matched_id
      ver <- gmt_versions[["c2"]]
    }
  } else if (p$type == "MitoCarta") {
    g <- extract_mito(mito_pathways, p$id)
    if (!is.null(g)) {
      genes <- g; src_used <- "MitoCarta 3.0"
      matched_id <- p$id; src_id <- p$id; ver <- "3.0"
    }
  }

  # min filter
  if (!is.null(genes) && length(genes) < spec$min_genes) genes <- NULL

  # max filter
  if (!is.null(genes) && length(genes) > MAX_GENES_GLOBAL) {
    log_warn("{gs}: {length(genes)} > {MAX_GENES_GLOBAL}，跳过")
    genes <- NULL
  }

  # extra
  if (!is.null(spec$extra) && !is.null(genes)) {
    e <- spec$extra; eg <- NULL
    if (e$type == "MitoCarta") eg <- extract_mito(mito_pathways, e$id)
    if (!is.null(eg)) {
      genes <- unique(c(genes, eg))
      src_id <- paste0(src_id, " + ", e$id)
    }
  }

  # fallback
  if (is.null(genes) && !is.null(spec$fallback)) {
    f <- spec$fallback; fg <- NULL
    if (f$type == "Hallmark") fg <- extract_gmt(msigdb$hallmark, f$id)
    else if (f$type == "KEGG") fg <- extract_gmt(msigdb$kegg, f$id)
    if (!is.null(fg) && length(fg) >= spec$min_genes) {
      genes <- fg; src_id <- f$id; is_fb <- TRUE
      log_warn("{gs}: fallback 到 '{f$id}'")
    }
  }

  if (is.null(genes) || length(genes) < spec$min_genes) {
    log_warn("{gs}: 未命中或过小，跳过")
    UNMATCHED <- rbind(UNMATCHED, data.frame(
      geneset = gs, group = spec$group, type = p$type,
      id = p$id, alt_ids = paste(spec$alt_ids, collapse = "|"),
      min_genes = spec$min_genes, stringsAsFactors = FALSE))
    next
  }

  GENE_SETS[[gs]] <- list(
    genes = sort(unique(toupper(genes))), group = spec$group,
    source = src_used, source_id = src_id, version = ver,
    reference = spec$reference, n_genes = length(unique(genes)),
    is_fallback = is_fb)
  META <- rbind(META, data.frame(
    geneset = gs, group = spec$group,
    source = src_used, source_id = src_id, version = ver,
    n_genes = length(unique(genes)),
    reference = spec$reference, is_fallback = is_fb,
    stringsAsFactors = FALSE))
  log_info("  {gs}: {length(unique(genes))} genes")
}

# TF_REFERENCE
if (!is.null(TF_REFERENCE)) {
  GENE_SETS$TF_REFERENCE <- list(
    genes = TF_REFERENCE$genes, group = "tf_reference",
    source = TF_REFERENCE$source, source_id = "Lambert2018_full",
    version = TF_REFERENCE$version, reference = TF_REFERENCE$reference,
    n_genes = length(TF_REFERENCE$genes), is_fallback = FALSE)
  META <- rbind(META, data.frame(
    geneset = "TF_REFERENCE", group = "tf_reference",
    source = TF_REFERENCE$source, source_id = "Lambert2018_full",
    version = TF_REFERENCE$version,
    n_genes = length(TF_REFERENCE$genes),
    reference = TF_REFERENCE$reference, is_fallback = FALSE,
    stringsAsFactors = FALSE))
}

# known_retinal_tfs
if (!is.null(HPA_RETINA) && !is.null(TF_REFERENCE)) {
  enriched <- c("Tissue enriched", "Group enriched", "Tissue enhanced")
  hpa_enriched <- HPA_RETINA$gene_clean[HPA_RETINA$spec_clean %in% enriched]
  krt <- sort(unique(intersect(hpa_enriched, TF_REFERENCE$genes)))
  if (length(krt) >= 5) {
    GENE_SETS$known_retinal_tfs <- list(
      genes = krt, group = "tf_retinal",
      source = "HPA retina ∩ Lambert",
      source_id = "HPA_retina_enriched_x_Lambert",
      version = "HPA_Lambert", reference = "HPA + Lambert 2018",
      n_genes = length(krt), is_fallback = FALSE)
    META <- rbind(META, data.frame(
      geneset = "known_retinal_tfs", group = "tf_retinal",
      source = "HPA retina ∩ Lambert",
      source_id = "HPA_retina_enriched_x_Lambert",
      version = "HPA_Lambert", n_genes = length(krt),
      reference = "HPA + Lambert 2018", is_fallback = FALSE,
      stringsAsFactors = FALSE))
  }
}

# IEG_ALL
if (!is.null(IEG_ALL)) {
  GENE_SETS$IEG_ALL <- list(
    genes = IEG_ALL$genes, group = "artifact_reference",
    source = IEG_ALL$source, source_id = "scCustomize_Homo_sapiens_IEG",
    version = IEG_ALL$version, reference = IEG_ALL$reference,
    n_genes = length(IEG_ALL$genes), is_fallback = FALSE)
  META <- rbind(META, data.frame(
    geneset = "IEG_ALL", group = "artifact_reference",
    source = IEG_ALL$source, source_id = "scCustomize_Homo_sapiens_IEG",
    version = IEG_ALL$version, n_genes = length(IEG_ALL$genes),
    reference = IEG_ALL$reference, is_fallback = FALSE,
    stringsAsFactors = FALSE))
}

# AMD anchor
if (length(AMD_TIER1) > 0) {
  GENE_SETS$AMD_positive_anchors_Tier1 <- list(
    genes = AMD_TIER1, group = "amd_anchor",
    source = "GSE135092 mmc5.xlsx", source_id = "Tier1",
    version = "file", reference = "Orozco 2020",
    n_genes = length(AMD_TIER1), is_fallback = FALSE)
  META <- rbind(META, data.frame(
    geneset = "AMD_positive_anchors_Tier1", group = "amd_anchor",
    source = "GSE135092 mmc5.xlsx", source_id = "Tier1",
    version = "file", n_genes = length(AMD_TIER1),
    reference = "Orozco 2020", is_fallback = FALSE,
    stringsAsFactors = FALSE))
}
if (length(AMD_TIER2) > 0) {
  GENE_SETS$AMD_positive_anchors_Tier2 <- list(
    genes = AMD_TIER2, group = "amd_anchor",
    source = "GSE135092 mmc5.xlsx", source_id = "Tier2",
    version = "file", reference = "Orozco 2020",
    n_genes = length(AMD_TIER2), is_fallback = FALSE)
  META <- rbind(META, data.frame(
    geneset = "AMD_positive_anchors_Tier2", group = "amd_anchor",
    source = "GSE135092 mmc5.xlsx", source_id = "Tier2",
    version = "file", n_genes = length(AMD_TIER2),
    reference = "Orozco 2020", is_fallback = FALSE,
    stringsAsFactors = FALSE))
}
AMD_ALL <- sort(unique(c(AMD_TIER1, AMD_TIER2)))
if (length(AMD_ALL) > 0) {
  GENE_SETS$AMD_positive_anchors <- list(
    genes = AMD_ALL, group = "amd_anchor",
    source = "GSE135092 mmc5.xlsx", source_id = "Tier1+Tier2",
    version = "file", reference = "Orozco 2020",
    n_genes = length(AMD_ALL), is_fallback = FALSE)
  META <- rbind(META, data.frame(
    geneset = "AMD_positive_anchors", group = "amd_anchor",
    source = "GSE135092 mmc5.xlsx", source_id = "Tier1+Tier2",
    version = "file", n_genes = length(AMD_ALL),
    reference = "Orozco 2020", is_fallback = FALSE,
    stringsAsFactors = FALSE))
}

log_info("总基因集数: {length(GENE_SETS)}")

# ============================================================
# 4. save
# ============================================================
save_rds(GENE_SETS, DIR_THIS, "gene_sets.rds")
save_table(META, DIR_THIS, "gene_sets_metadata.csv")
if (nrow(UNMATCHED) > 0) {
  save_table(UNMATCHED, DIR_THIS, "gene_sets_unmatched.csv")
}

full_list <- do.call(rbind, lapply(names(GENE_SETS), function(nm) {
  gs <- GENE_SETS[[nm]]
  data.frame(geneset = nm, group = gs$group,
             source = gs$source, gene = gs$genes,
             stringsAsFactors = FALSE)
}))
save_table(full_list, DIR_THIS, "gene_sets_full.csv")

# ============================================================
# 5. summary
# ============================================================
cat("\n========================================\n")
cat("=== 01_2_prepare_gene_sets.r 汇总 ===\n")
cat("========================================\n")
cat("基因集总数:", length(GENE_SETS), "\n")
cat("fallback:", sum(META$is_fallback %in% TRUE), "\n\n")

cat("按 group:\n"); print(table(META$group))

if (nrow(UNMATCHED) > 0) {
  cat("\n未命中 id 数:", nrow(UNMATCHED), "\n")
  cat("未命中清单:\n")
  print(UNMATCHED, row.names = FALSE)
}

cat("\n基因集清单:\n")
print(META[, c("geneset", "group", "source", "n_genes")], row.names = FALSE)

cat("\n--- 输出位置 ---\n")
cat("目录:", DIR_THIS, "\n")
cat("  data/gene_sets.rds\n")
cat("  tables/gene_sets_metadata.csv\n")
cat("  tables/gene_sets_full.csv\n")
if (nrow(UNMATCHED) > 0) cat("  tables/gene_sets_unmatched.csv\n")

log_info("=== 01_2_prepare_gene_sets.r done ===")

