# Settings for the analysis. CANCER ("BRCA" or "SKCM") is set in main.R.

# ---- Cancer-specific settings ----

if (CANCER == "BRCA") {
  RNA_RDS_FILE    <- "BRCA-rna_data.RDS"
  MIRNA_RDS_FILE  <- "BRCA-mirna_data.RDS"
  CELL_LINES_FILE <- "Breast_Cancer_Cell_Lines.csv"
} else if (CANCER == "SKCM") {
  RNA_RDS_FILE    <- "SKCM_rna_data.RDS"
  MIRNA_RDS_FILE  <- "SKCM_mirna_data.RDS"
  CELL_LINES_FILE <- "Melanoma_cell_lines_hg38_FULL.csv"
} else {
  stop("CANCER must be \"BRCA\" or \"SKCM\"")
}

STUDY     <- paste0("TCGA-", CANCER)
DATA_DIR  <- here("data_files", STUDY)   # inputs and saved outputs
PLOTS_DIR <- here("plots", CANCER)

# ---- Files in DATA_DIR ----

GENE_MAPPING_FILE          <- "gene_mapping.csv"
RNA_FILTERED_FILE          <- "rNA_filtered_data.csv"
MIRNA_FILTERED_FILE        <- "miRNA_filtered_data.csv"
CORRELATION_FILE           <- paste0("TCGA_", CANCER, "_mRNA_miRNA_sig_cor.csv")
MIRNA_TF_INTERACTIONS_FILE <- "interactions_TFs.csv"
TF_TARGET_DF_FILE          <- "TF_target_df.csv"
TF_TARGET_COUNT_FILE       <- "TF_target_gene_count.csv"
MIRNA_DIRECT_TARGETS_FILE  <- "miRNA_direct_targets.csv"
ORA_POSCOR_FILE            <- "overlap_analysis_poscor_sig.csv"
ORA_NEGCOR_FILE            <- "overlap_analysis_negcor_sig.csv"

# ---- Shared files ----

HUMAN_TFS_FILE <- here("data_files", "HumanTFs.csv")  # http://humantfs.ccbr.utoronto.ca
TF_TARGETS_DIR <- here("data_files", "TF_targets")    # ChIP-Atlas target gene files

# ---- Analysis parameters ----

EXPRESSION_THRESHOLD   <- 60     # drop genes/miRNAs that are zero in more than this % of samples
CORRELATION_METHOD     <- "pearson"
CORRELATION_FDR        <- 0.05
ORA_FDR                <- 0.05
CELL_LINE_THRESHOLD    <- 0.5    # keep TF targets with peaks in at least this fraction of cell lines
CHIP_ATLAS_GENOME      <- "hg38"
CHIP_ATLAS_TSS_KB      <- 5      # peak distance from TSS (kb)
VALIDATION_EXPERIMENTS <- c("luciferase", "western blot")  # for the interaction plots
ORA_CORES              <- 1      # >1 runs the ORA loop in parallel
