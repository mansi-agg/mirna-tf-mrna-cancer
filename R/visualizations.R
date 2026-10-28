# Data behind each figure, and make_all_plots() to build and save every figure.

# Number of miRNAs per correlation direction and TF-overlap category
build_mirna_distribution_data <- function(cor_result, overlap, ora_poscor, ora_negcor) {
  conditions <- c("Targets TF and significant overlap",
                  "Targets TF but no significant overlap",
                  "Not targeting TFs")

  n_targeting_tf <- dplyr::n_distinct(overlap$miRNA)
  n_not_targeting <- length(setdiff(unique(cor_result$miRNA), unique(overlap$miRNA)))
  n_sig_pos <- dplyr::n_distinct(ora_poscor$miRNA)
  n_sig_neg <- dplyr::n_distinct(ora_negcor$miRNA)

  tibble::tibble(correlation_type = rep(c("miRNAs with positively correlated genes",
                                          "miRNAs with negatively correlated genes"), each = 3),
         condition = factor(rep(conditions, 2), levels = conditions),
         value = c(n_sig_pos, n_targeting_tf - n_sig_pos, n_not_targeting,
                   n_sig_neg, n_targeting_tf - n_sig_neg, n_not_targeting)) %>%
    dplyr::filter(value > 0)
}

# miRNA-TF pairs validated by any of `experiments` (e.g. luciferase, western blot)
build_validated_mir_tf <- function(interactions_tfs, experiments) {
  interactions_tfs %>%
    dplyr::mutate(experiment_lower = tolower(experiment)) %>%
    dplyr::filter(!is.na(experiment_lower),
                  stringr::str_detect(experiment_lower, paste(experiments, collapse = "|"))) %>%
    dplyr::mutate(miR_TF = make_mir_tf_label(miRNA, target_symbol)) %>%
    dplyr::distinct(miR_TF, .keep_all = TRUE) %>%
    dplyr::select(miR_TF, experiment_lower)
}

# Z-scored miRNA and TF expression for every validated miRNA-TF pair
build_interaction_data <- function(overlap, mirna_zscores, gene_zscores, validated) {
  overlap %>%
    dplyr::left_join(dplyr::rename(mirna_zscores, miRNA_exp = z_transform_exp), by = "miRNA") %>%
    dplyr::left_join(dplyr::select(gene_zscores, RNA_symbol, TF_exp = z_transform_exp),
                     by = c("TF" = "RNA_symbol"), relationship = "many-to-many") %>%
    dplyr::select(miRNA, TF, miRNA_exp, TF_exp) %>%
    dplyr::arrange(miRNA, TF) %>%
    dplyr::mutate(miR_TF = make_mir_tf_label(miRNA, TF)) %>%
    dplyr::inner_join(validated, by = "miR_TF")
}

# Long format (one row per miRNA / TF value), only the significant pairs
build_interaction_long_data <- function(interaction_data, ora_sig) {
  interaction_data %>%
    dplyr::filter(miR_TF %in% ora_sig$miR_TF) %>%
    tidyr::pivot_longer(c(miRNA_exp, TF_exp), names_to = "variable",
                        values_to = "value", values_drop_na = TRUE) %>%
    dplyr::mutate(variable = factor(dplyr::recode(variable, miRNA_exp = "miRNA", TF_exp = "TF"),
                                    levels = c("miRNA", "TF")))
}

# Build every figure, save them to plots_dir, and return them as a list
make_all_plots <- function(cancer, plots_dir, rna_expr, mirna_expr, cor_result,
                           interactions_tfs, overlap, tf_target_count,
                           ora_poscor, ora_negcor, validation_experiments) {
  plots <- list()

  plots$targets_vs_replicates <- plot_targets_vs_replicates(tf_target_count)
  plots$top_pairs_positive <- plot_top_mir_tf_pairs(ora_poscor, "pos")
  plots$top_pairs_negative <- plot_top_mir_tf_pairs(ora_negcor, "neg")

  distribution <- build_mirna_distribution_data(cor_result, overlap, ora_poscor, ora_negcor)
  plots$mirna_distribution <- plot_mirna_distribution(distribution)

  validated <- build_validated_mir_tf(interactions_tfs, validation_experiments)
  interaction_data <- build_interaction_data(overlap,
                                             expression_zscores(mirna_expr, "miRNA"),
                                             gene_expression_zscores(rna_expr, cor_result),
                                             validated)
  plots$interactions_positive <- plot_interactions(
    build_interaction_long_data(interaction_data, ora_poscor), line_color = "#8939FD")
  plots$interactions_negative <- plot_interactions(
    build_interaction_long_data(interaction_data, ora_negcor), line_color = "#9F0A0C")

  titles <- c(
    targets_vs_replicates = "TF target counts are not correlated with the number of biological replicates",
    top_pairs_positive    = "Top miR-TF pairs - Positive Correlation",
    top_pairs_negative    = "Top miR-TF pairs - Negative Correlation",
    mirna_distribution    = "Distribution of miRNAs Based on Correlation with Genes and Overlap with TF Targets",
    interactions_positive = "Interaction Plot - miRs with positively correlated genes",
    interactions_negative = "Interaction Plot - miRs with negatively correlated genes"
  )
  for (name in names(plots)) {
    save_plot(plots[[name]], paste0(cancer, " - ", titles[[name]]), plots_dir)
  }
  log_msg("Saved ", length(plots), " plots to ", plots_dir)

  plots
}
