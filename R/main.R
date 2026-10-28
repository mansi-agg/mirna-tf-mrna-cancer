# Exploring microRNA-mRNA relationships mediated by transcription factors in cancer.
#
# Set the options below and source this file. Each step either reruns its
# analysis (TRUE) or reads its saved output from data_files/TCGA-<CANCER>/ (FALSE).
# If you rerun a step, also rerun the steps after it that use its output.

CANCER <- "SKCM"   # "BRCA" or "SKCM"

RUN_DOWNLOAD_TCGA         <- FALSE  # 1. download TCGA RNA / miRNA data (several GB)
RUN_EXPRESSION            <- FALSE  # 2. filter and log-transform expression data
RUN_CORRELATION           <- FALSE  # 3. mRNA-miRNA correlation
RUN_MIRNA_TF_INTERACTIONS <- FALSE  # 4. miRNA -> TF interactions from multiMiR (needs internet)
RUN_DOWNLOAD_CHIP_ATLAS   <- FALSE  # 5. download ChIP-Atlas TF target files (needs internet)
RUN_TF_TARGETS            <- FALSE  # 6. TF target genes in the cancer cell lines
RUN_MIRNA_DIRECT_TARGETS  <- FALSE  # 7. miRNA direct targets from multiMiR (needs internet)
RUN_ORA                   <- FALSE  # 8. over-representation analysis
MAKE_PLOTS                <- TRUE   # 9. make and save the figures


library(here)
library(tidyverse)
library(ggrepel)

source(here("R", "config.R"))
source(here("R", "utils.R"))
source(here("R", "data_preparation.R"))
source(here("R", "correlation_analysis.R"))
source(here("R", "tf_interactions.R"))
source(here("R", "tf_targets.R"))
source(here("R", "ora_analysis.R"))
source(here("R", "plots.R"))
source(here("R", "visualizations.R"))

log_msg("Running ", STUDY)


# ---- 1. TCGA data ----

if (RUN_DOWNLOAD_TCGA) {
  download_tcga_data(STUDY, file.path(DATA_DIR, RNA_RDS_FILE), file.path(DATA_DIR, MIRNA_RDS_FILE))
}


# ---- 2. Expression data ----

if (RUN_EXPRESSION) {
  rna_data <- readRDS(file.path(DATA_DIR, RNA_RDS_FILE))
  mirna_data <- readRDS(file.path(DATA_DIR, MIRNA_RDS_FILE))

  gene_mapping <- get_gene_mapping(rna_data)
  expression <- prepare_expression_data(rna_data, mirna_data, EXPRESSION_THRESHOLD)
  rna_expr <- expression$rna
  mirna_expr <- expression$mirna

  save_table(gene_mapping, file.path(DATA_DIR, GENE_MAPPING_FILE))
  save_matrix(rna_expr, file.path(DATA_DIR, RNA_FILTERED_FILE))
  save_matrix(mirna_expr, file.path(DATA_DIR, MIRNA_FILTERED_FILE))
} else {
  rna_expr <- read_matrix(file.path(DATA_DIR, RNA_FILTERED_FILE))
  mirna_expr <- read_matrix(file.path(DATA_DIR, MIRNA_FILTERED_FILE))
  # The gene ID -> symbol table is only needed to rerun steps 3 and 8
  if (RUN_CORRELATION || RUN_ORA) {
    gene_mapping <- load_gene_mapping(file.path(DATA_DIR, GENE_MAPPING_FILE),
                                      file.path(DATA_DIR, RNA_RDS_FILE))
  }
}


# ---- 3. mRNA-miRNA correlation ----

if (RUN_CORRELATION) {
  cor_result <- compute_mrna_mirna_correlation(rna_expr, mirna_expr, gene_mapping,
                                               method = CORRELATION_METHOD, fdr = CORRELATION_FDR)
  save_table(cor_result, file.path(DATA_DIR, CORRELATION_FILE))
} else {
  cor_result <- read_table(file.path(DATA_DIR, CORRELATION_FILE))
}


# ---- 4. miRNA -> TF interactions ----

known_tfs <- load_known_tfs(HUMAN_TFS_FILE)

if (RUN_MIRNA_TF_INTERACTIONS) {
  interactions_tfs <- query_mirna_tf_interactions(rownames(mirna_expr), known_tfs)
  save_table(interactions_tfs, file.path(DATA_DIR, MIRNA_TF_INTERACTIONS_FILE))
} else {
  interactions_tfs <- read_interactions(file.path(DATA_DIR, MIRNA_TF_INTERACTIONS_FILE))
}


# ---- 5. ChIP-Atlas TF target files ----

if (RUN_DOWNLOAD_CHIP_ATLAS) {
  download_chip_atlas_files(unique(interactions_tfs$target_symbol), TF_TARGETS_DIR,
                            CHIP_ATLAS_GENOME, CHIP_ATLAS_TSS_KB)
}
chip_atlas_files <- list_chip_atlas_files(TF_TARGETS_DIR, CHIP_ATLAS_GENOME, CHIP_ATLAS_TSS_KB)

# miRNA-TF pairs where the miRNA has correlated genes and the TF has ChIP-Atlas data
overlap <- build_mirna_tf_overlap(cor_result, interactions_tfs, names(chip_atlas_files))


# ---- 6. TF target genes ----

if (RUN_TF_TARGETS) {
  cell_lines <- read_table(file.path(DATA_DIR, CELL_LINES_FILE))
  tf_targets <- build_tf_targets(chip_atlas_files, cell_lines, CELL_LINE_THRESHOLD)
  tf_target_df <- tf_targets$target_df
  tf_target_count <- tf_targets$target_count

  save_table(tf_target_df, file.path(DATA_DIR, TF_TARGET_DF_FILE))
  save_table(tf_target_count, file.path(DATA_DIR, TF_TARGET_COUNT_FILE))
} else {
  tf_target_df <- read_table(file.path(DATA_DIR, TF_TARGET_DF_FILE))
  tf_target_count <- read_table(file.path(DATA_DIR, TF_TARGET_COUNT_FILE))
}


# ---- 7. miRNA direct targets ----

if (RUN_MIRNA_DIRECT_TARGETS) {
  mirna_direct_targets <- query_mirna_targets(unique(overlap$miRNA))
  save_table(mirna_direct_targets, file.path(DATA_DIR, MIRNA_DIRECT_TARGETS_FILE))
} else {
  mirna_direct_targets <- read_interactions(file.path(DATA_DIR, MIRNA_DIRECT_TARGETS_FILE))
}


# ---- 8. Over-representation analysis ----

if (RUN_ORA) {
  background_genes <- gene_mapping$gene_symbol[gene_mapping$gene_id %in% rownames(rna_expr)]
  ora <- run_ora_analysis(overlap, cor_result, tf_target_df, mirna_direct_targets,
                          background_genes, fdr = ORA_FDR, cores = ORA_CORES)
  ora_poscor <- ora$poscor
  ora_negcor <- ora$negcor

  save_table(ora_poscor, file.path(DATA_DIR, ORA_POSCOR_FILE))
  save_table(ora_negcor, file.path(DATA_DIR, ORA_NEGCOR_FILE))
} else {
  ora_poscor <- read_table(file.path(DATA_DIR, ORA_POSCOR_FILE))
  ora_negcor <- read_table(file.path(DATA_DIR, ORA_NEGCOR_FILE))
}

log_msg("Significant miRNA-TF pairs: ", nrow(ora_poscor), " positive, ", nrow(ora_negcor), " negative")


# ---- 9. Plots ----

if (MAKE_PLOTS) {
  plots <- make_all_plots(CANCER, PLOTS_DIR, rna_expr, mirna_expr, cor_result,
                          interactions_tfs, overlap, tf_target_count,
                          ora_poscor, ora_negcor, VALIDATION_EXPERIMENTS)
}

log_msg("Done: ", STUDY)
