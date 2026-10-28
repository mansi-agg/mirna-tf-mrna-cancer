# =============================================================================
# visualizations.R
# Prepares the data behind every figure (runAnalysis.Rmd section 6) and
# generates / saves all figures for one cancer type with the builders in
# plots.R.
# =============================================================================

# -----------------------------------------------------------------------------
# Plot data
# -----------------------------------------------------------------------------

#' Counts of miRNAs per correlation direction and TF-overlap category.
build_mirna_distribution_data <- function(cor_result, overlap, ora_poscor, ora_negcor) {
  conditions <- c("Targets TF and significant overlap",
                  "Targets TF but no significant overlap",
                  "Not targeting TFs")

  n_targeting_tf <- dplyr::n_distinct(overlap$miRNA)
  n_not_targeting <- length(setdiff(unique(cor_result$miRNA), unique(overlap$miRNA)))
  n_sig_pos <- dplyr::n_distinct(ora_poscor$miRNA)
  n_sig_neg <- dplyr::n_distinct(ora_negcor$miRNA)

  tibble::tibble(
    correlation_type = rep(c("miRNAs with positively correlated genes",
                             "miRNAs with negatively correlated genes"), each = 3),
    condition = factor(rep(conditions, 2), levels = conditions),
    value = c(n_sig_pos, n_targeting_tf - n_sig_pos, n_not_targeting,
              n_sig_neg, n_targeting_tf - n_sig_neg, n_not_targeting)
  ) %>%
    dplyr::filter(value > 0)
}

#' miRNA-TF pairs validated by any of `experiments` (e.g. luciferase, western blot).
build_validated_mir_tf <- function(interactions_tfs, experiments) {
  pattern <- paste(experiments, collapse = "|")
  interactions_tfs %>%
    dplyr::mutate(experiment_lower = tolower(experiment)) %>%
    dplyr::filter(!is.na(experiment_lower),
                  stringr::str_detect(experiment_lower, pattern)) %>%
    dplyr::mutate(miR_TF = make_mir_tf_label(miRNA, target_symbol)) %>%
    dplyr::distinct(miR_TF, .keep_all = TRUE) %>%
    dplyr::select(miR_TF, experiment_lower)
}

#' Z-scored miRNA and TF expression for every validated miRNA-TF pair.
build_interaction_data <- function(overlap, mirna_zscores, gene_zscores, validated) {
  overlap %>%
    dplyr::left_join(dplyr::rename(mirna_zscores, miRNA_exp = z_transform_exp), by = "miRNA") %>%
    dplyr::left_join(gene_zscores %>%
                       dplyr::select(RNA_symbol, TF_exp = z_transform_exp),
                     by = c("TF" = "RNA_symbol"), relationship = "many-to-many") %>%
    dplyr::select(miRNA, TF, miRNA_exp, TF_exp) %>%
    dplyr::arrange(miRNA, TF) %>%
    dplyr::mutate(miR_TF = make_mir_tf_label(miRNA, TF)) %>%
    dplyr::inner_join(validated, by = "miR_TF")
}

#' Long format (one row per miRNA / TF value) restricted to significant pairs.
build_interaction_long_data <- function(interaction_data, ora_sig) {
  interaction_data %>%
    dplyr::filter(miR_TF %in% unique(ora_sig$miR_TF)) %>%
    tidyr::pivot_longer(c(miRNA_exp, TF_exp), names_to = "variable",
                        values_to = "value", values_drop_na = TRUE) %>%
    dplyr::mutate(variable = factor(dplyr::recode(variable, miRNA_exp = "miRNA", TF_exp = "TF"),
                                    levels = c("miRNA", "TF")))
}

# -----------------------------------------------------------------------------
# Generate all figures
# -----------------------------------------------------------------------------

#' Build (and save, if cfg$plots$save) every figure.
#'
#' @param res Named list of pipeline results (MirTfPipeline$results()).
#' @param cfg Config list.
#' @return Named list of ggplot objects.
generate_all_plots <- function(res, cfg) {
  pc <- cfg$plots
  dir <- cfg$paths$plots_dir
  prefix <- cfg$plot_title_prefix
  plots <- list()

  add_plot <- function(key, p, title) {
    save_plot(p, paste0(prefix, " - ", title), pc, dir)
    plots[[key]] <<- p
  }

  add_plot("targets_vs_replicates",
           plot_targets_vs_replicates(res$tf_target_count, pc),
           "TF target counts are not correlated with the number of biological replicates")

  add_plot("top_pairs_positive",
           plot_top_mir_tf_pairs(res$ora_poscor, pc, "pos"),
           "Top miR-TF pairs - Positive Correlation")

  add_plot("top_pairs_negative",
           plot_top_mir_tf_pairs(res$ora_negcor, pc, "neg"),
           "Top miR-TF pairs - Negative Correlation")

  distribution <- build_mirna_distribution_data(res$cor_result, res$overlap,
                                                res$ora_poscor, res$ora_negcor)
  add_plot("mirna_distribution",
           plot_mirna_distribution(distribution, pc),
           "Distribution of miRNAs Based on Correlation with Genes and Overlap with TF Targets")

  validated <- build_validated_mir_tf(res$interactions_tfs, cfg$analysis$validation_experiments)
  interaction_data <- build_interaction_data(res$overlap, res$mirna_zscores,
                                             res$gene_zscores, validated)

  add_plot("interactions_positive",
           plot_interactions(build_interaction_long_data(interaction_data, res$ora_poscor),
                             pc, pc$colors$interaction_pos),
           "Interaction Plot - miRs with positively correlated genes")

  add_plot("interactions_negative",
           plot_interactions(build_interaction_long_data(interaction_data, res$ora_negcor),
                             pc, pc$colors$interaction_neg),
           "Interaction Plot - miRs with negatively correlated genes")

  if (isTRUE(pc$save)) log_msg("Saved ", length(plots), " plots to ", dir, verbose = cfg$verbose)
  plots
}
