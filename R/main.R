# =============================================================================
# main.R
# Exploring microRNA-mRNA Relationships Mediated by Transcription Factors
# in Cancer.
#
# One-button run: set the options below and click "Source" in RStudio
# (or run `Rscript R/main.R` from the project root).
#
# Variables per cancer type live in R/config_brca.R and R/config_skcm.R;
# shared settings live in R/config_common.R.
# =============================================================================

# ---- 1. Choose what to run --------------------------------------------------

# "BRCA", "SKCM", or c("BRCA", "SKCM") to run both
CANCER <- "BRCA"

# Steps to recompute instead of reading saved outputs from data_files/<study>/.
#   character(0)              -> read everything that already exists (fast)
#   "all"                     -> recompute every analysis step
#   c("tf_targets", "ora")    -> recompute just these (and, by default,
#                                everything downstream of them)
# Analysis steps: gene_mapping, expression, correlation, mirna_tf_interactions,
#                 tf_targets, mirna_direct_targets, ora
# Downloads (never part of "all"): tcga_data, chip_atlas
REGENERATE <- character(0)

# Build and save the figures
MAKE_PLOTS <- TRUE

# Optional overrides of any config value, e.g.
#   list(plots = list(save = FALSE), analysis = list(ora_cores = 4L))
OVERRIDES <- list()


# ---- 2. Load libraries and modules -------------------------------------------

suppressPackageStartupMessages({
  library(here)
  library(R6)
  library(tidyverse)
  library(SummarizedExperiment)
  library(clusterProfiler)
  library(multiMiR)
  library(ggrepel)
})

module_files <- c(
  "utils.R",
  "config_common.R", "config_brca.R", "config_skcm.R",
  "data_preparation.R", "correlation_analysis.R", "tf_interactions.R",
  "tf_targets.R", "ora_analysis.R",
  "plots.R", "visualizations.R",
  "analysis_pipeline.R"
)
for (f in module_files) source(here("R", f))


# ---- 3. Run ------------------------------------------------------------------

pipelines <- list()
results <- list()

for (cancer in CANCER) {
  cfg <- load_config(cancer, regenerate = REGENERATE, overrides = OVERRIDES)
  # With more than one cancer, keep each one's results in `results$<CANCER>`
  # instead of overwriting the same global variables.
  if (length(CANCER) > 1) cfg$export_to_global <- FALSE

  pipelines[[cancer]] <- MirTfPipeline$new(cfg)
  results[[cancer]] <- pipelines[[cancer]]$run(make_plots = MAKE_PLOTS)
}

# Handy afterwards:
#   pipelines$BRCA$get("correlation")    # any single step
#   pipelines$BRCA$status()              # what is read vs regenerated
#   results$BRCA$plots$top_pairs_positive
