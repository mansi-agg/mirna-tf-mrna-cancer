# Settings for the analysis. CANCER ("BRCA" or "SKCM") is set in main.R.

# ---- Cancer-specific settings ----

if (CANCER == "BRCA") {
  RNA_RDS_FILE    <- "BRCA-rna_data.RDS"
  MIRNA_RDS_FILE  <- "BRCA-mirna_data.RDS"
} else if (CANCER == "SKCM") {
  RNA_RDS_FILE    <- "SKCM_rna_data.RDS"
  MIRNA_RDS_FILE  <- "SKCM_mirna_data.RDS"
} else {
  stop("CANCER must be \"BRCA\" or \"SKCM\"")
}

STUDY       <- paste0("TCGA-", CANCER)
DATA_DIR    <- here("data_files", STUDY)      # inputs and saved outputs of steps 1-5 and 7
RESULTS_DIR <- here("data_files_v2", STUDY)   # saved outputs of steps 6 and 8 (TF targets, ORA)
PLOTS_DIR   <- here("plots_v2", CANCER)

# Cell lines made by cell_lines/build_cell_lines.R. The previous versions were
# data_files/TCGA-BRCA/Breast_Cancer_Cell_Lines.csv and
# data_files/TCGA-SKCM/Melanoma_cell_lines_hg38_FULL.csv (used with RESULTS_DIR = DATA_DIR
# and PLOTS_DIR = here("plots", CANCER)).
CELL_LINES_FILE <- here("cell_lines", paste0(CANCER, "_cell_lines.csv"))

# ---- Files in DATA_DIR / RESULTS_DIR ----

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
MIN_CELL_LINES         <- 1      # keep TF targets bound in at least this many different cancer cell lines
                                 # (was: bound in >= 50% of the cell-line experiments)
CHIP_ATLAS_GENOME      <- "hg38"
CHIP_ATLAS_TSS_KB      <- 5      # peak distance from TSS (kb)
VALIDATION_EXPERIMENTS <- c("luciferase", "western blot")  # for the interaction plots
ORA_CORES              <- 1      # >1 runs the ORA loop in parallel
