# =============================================================================
# config_brca.R
# Variables for the TCGA-BRCA (breast cancer) analysis.
# Anything not set here falls back to default_config() in config_common.R.
# =============================================================================

brca_config <- function(cfg = default_config()) {
  cfg$cancer            <- "BRCA"
  cfg$study             <- "TCGA-BRCA"
  cfg$plot_title_prefix <- "BRCA"

  study_dir <- file.path(cfg$paths$common_dir, cfg$study)
  cfg$paths$study_dir   <- study_dir
  cfg$paths$results_dir <- study_dir
  cfg$paths$plots_dir   <- file.path(here::here(), "plots_main", cfg$cancer)

  # File names inside study_dir (inputs and saved outputs)
  cfg$files <- utils::modifyList(cfg$files, list(
    rna_rds              = "BRCA-rna_data.RDS",
    mirna_rds            = "BRCA-mirna_data.RDS",
    cell_lines           = "Breast_Cancer_Cell_Lines.csv",
    gene_mapping         = "gene_mapping.csv",
    rna_filtered         = "rNA_filtered_data.csv",
    mirna_filtered       = "miRNA_filtered_data.csv",
    correlation          = "TCGA_BRCA_mRNA_miRNA_sig_cor.csv",
    mirna_tf_interactions = "interactions_TFs.csv",
    tf_target_df         = "TF_target_df.csv",
    tf_target_count      = "TF_target_gene_count.csv",
    mirna_direct_targets = "miRNA_direct_targets.csv",
    ora_poscor           = "overlap_analysis_poscor_sig.csv",
    ora_negcor           = "overlap_analysis_negcor_sig.csv"
  ))

  # BRCA-specific parameter overrides go here, e.g.
  # cfg$filters$expression_threshold <- 60
  # cfg$plots$y_limit <- 10.5

  cfg
}
