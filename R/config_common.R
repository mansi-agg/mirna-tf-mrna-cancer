# =============================================================================
# config_common.R
# Settings shared by every cancer type. Cancer-specific files (config_brca.R,
# config_skcm.R) start from these defaults and override what differs.
# =============================================================================

#' Names of every pipeline step that can be regenerated.
#'
#' "all" in `regenerate` expands to `REGENERABLE_STEPS`. The two download steps
#' are deliberately excluded from "all" and must be named explicitly.
REGENERABLE_STEPS <- c(
  "gene_mapping",          # gene_mapping.csv              (needs the TCGA RDS files)
  "expression",            # rNA/miRNA_filtered_data.csv   (needs the TCGA RDS files)
  "correlation",           # TCGA_<X>_mRNA_miRNA_sig_cor.csv
  "mirna_tf_interactions", # interactions_TFs.csv          (multiMiR query)
  "tf_targets",            # TF_target_df.csv + TF_target_gene_count.csv
  "mirna_direct_targets",  # miRNA_direct_targets.csv      (multiMiR query)
  "ora"                    # overlap_analysis_{pos,neg}cor_sig.csv
)
DOWNLOAD_STEPS <- c(
  "tcga_data",             # GDC download of the raw RNA / miRNA data
  "chip_atlas"             # ChIP-Atlas TF target files
)

default_config <- function() {
  root <- here::here()

  list(
    # ---- Paths ---------------------------------------------------------------
    paths = list(
      common_dir     = file.path(root, "data_files"),
      tf_targets_dir = file.path(root, "data_files", "TF_targets"),
      # study_dir / results_dir / plots_dir are filled in by the cancer config
      study_dir      = NULL,  # where inputs and previously saved outputs live
      results_dir    = NULL,  # where (re)generated outputs are written
      plots_dir      = NULL
    ),

    # ---- Common input files --------------------------------------------------
    files = list(
      human_tfs = "HumanTFs.csv"  # http://humantfs.ccbr.utoronto.ca (Lambert et al. 2018)
    ),

    # ---- Data filters --------------------------------------------------------
    filters = list(
      # Drop transcripts that are zero in more than this % of samples
      expression_threshold = 60,
      gene_type            = "protein_coding",
      rna_assay            = "tpm_unstrand",
      mirna_value_prefix   = "reads_per_million_miRNA_mapped_"
    ),

    # ---- Analysis parameters -------------------------------------------------
    analysis = list(
      correlation_method      = "pearson",
      correlation_fdr         = 0.05,
      ora_fdr                 = 0.05,
      p_adjust_method         = "BH",
      multimir_table          = "validated",
      # Keep TF target genes with peaks in at least this fraction of cell lines
      gene_peaks_in_cell_lines_threshold = 0.5,
      # Experiments counted as strong validation for the interaction plots
      validation_experiments  = c("luciferase", "western blot"),
      # Cores for the per-miRNA ORA loop (>1 uses forked parallelism)
      ora_cores               = 1L
    ),

    # ---- ChIP-Atlas TF target genes -----------------------------------------
    chip_atlas = list(
      genome                 = "hg38",
      peak_distance_from_tss = 5,  # kb
      base_url               = "https://chip-atlas.dbcls.jp/data"
    ),

    # ---- Plots ---------------------------------------------------------------
    plots = list(
      save            = TRUE,
      width           = 8,
      height          = 5,
      point_size      = 7,
      alpha           = 0.7,
      annot_text_size = 5,
      top_n_labels    = 5,
      y_limit         = 10.5,
      colors = list(
        positive         = "#0D84CD",
        negative         = "#9C354B",
        interaction_pos  = "#8939FD",
        interaction_neg  = "#9F0A0C",
        mirna            = "#89b7d7",
        tf               = "#ffbb78",
        distribution     = c("#1a9850", "#4575b4", "#999999")
      )
    ),

    # ---- Run behaviour -------------------------------------------------------
    # Steps to recompute instead of reading from disk. Use "all" for every
    # analysis step; downloads ("tcga_data", "chip_atlas") must be named.
    regenerate = character(0),
    # Also regenerate every step downstream of a regenerated one
    cascade_regeneration = TRUE,
    # Copy all results into the global environment after run()
    export_to_global = TRUE,
    verbose = TRUE
  )
}

#' Build the full config for a cancer type.
#'
#' @param cancer "BRCA" or "SKCM".
#' @param regenerate Character vector of step names (see REGENERABLE_STEPS),
#'   "all", or character(0) to read everything from disk.
#' @param overrides Optional nested list merged on top of the result.
load_config <- function(cancer, regenerate = character(0), overrides = list()) {
  cancer <- toupper(cancer)
  builder <- switch(cancer,
    BRCA = brca_config,
    SKCM = skcm_config,
    stop("Unknown cancer type '", cancer, "'. Use \"BRCA\" or \"SKCM\".", call. = FALSE)
  )

  cfg <- builder(default_config())
  cfg$regenerate <- expand_regenerate(regenerate)
  cfg <- utils::modifyList(cfg, overrides)
  validate_config(cfg)
  cfg
}

expand_regenerate <- function(regenerate) {
  regenerate <- unique(as.character(regenerate))
  if ("all" %in% regenerate) {
    regenerate <- union(setdiff(regenerate, "all"), REGENERABLE_STEPS)
  }
  unknown <- setdiff(regenerate, c(REGENERABLE_STEPS, DOWNLOAD_STEPS))
  if (length(unknown)) {
    stop("Unknown step(s) in `regenerate`: ", paste(unknown, collapse = ", "),
         "\nValid steps: ", paste(c(REGENERABLE_STEPS, DOWNLOAD_STEPS), collapse = ", "),
         call. = FALSE)
  }
  regenerate
}

validate_config <- function(cfg) {
  required <- c("cancer", "study", "plot_title_prefix")
  missing <- required[!required %in% names(cfg)]
  if (length(missing)) stop("Config is missing: ", paste(missing, collapse = ", "), call. = FALSE)

  for (p in c("study_dir", "results_dir", "plots_dir")) {
    if (is.null(cfg$paths[[p]])) stop("Config path `", p, "` is not set.", call. = FALSE)
  }
  if (!dir.exists(cfg$paths$study_dir)) {
    stop("Study directory not found: ", cfg$paths$study_dir, call. = FALSE)
  }
  invisible(cfg)
}
