# =============================================================================
# plots.R
# ggplot builders. Each function takes ready-made plot data plus the `plots`
# section of the config and returns a ggplot object; saving is done by
# save_plot(). Plot data is prepared in visualizations.R.
# =============================================================================

save_plot <- function(p, title, plot_cfg, plots_dir) {
  if (!isTRUE(plot_cfg$save)) return(invisible(NULL))
  ensure_dir(plots_dir)
  path <- file.path(plots_dir, title_to_filename(title))
  ggplot2::ggsave(plot = p, filename = path, width = plot_cfg$width,
                  height = plot_cfg$height, bg = "white")
  invisible(path)
}

#' Supplementary: TF target counts vs number of biological replicates.
plot_targets_vs_replicates <- function(target_count, plot_cfg) {
  ggplot2::ggplot(target_count,
                  ggplot2::aes(x = log2(replicate.count), y = log2(target.count))) +
    ggplot2::geom_point() +
    ggplot2::labs(x = "Log2 Replicate Count",
                  y = "Log2 Target Gene Count",
                  title = "TF target counts are not correlated with the number of biological replicates") +
    ggplot2::theme_minimal()
}

#' Scatter of significant miRNA-TF pairs: top correlation vs -log10(q).
#'
#' @param direction "pos" or "neg".
plot_top_mir_tf_pairs <- function(df, plot_cfg, direction = c("pos", "neg")) {
  direction <- match.arg(direction)
  df <- as.data.frame(df)

  if (direction == "neg") {
    df$top_cor_coef <- abs(df$top_cor_coef)
    point_color <- plot_cfg$colors$negative
    x_label <- "Absolute Correlation Coefficient of \nTop Correlated Gene"
  } else {
    point_color <- plot_cfg$colors$positive
    x_label <- "Correlation Coefficient of \nTop Correlated Gene"
  }

  top_pairs <- df %>%
    dplyr::arrange(dplyr::desc(top_cor_coef), qval_plot) %>%
    utils::head(plot_cfg$top_n_labels)

  ggplot2::ggplot(df, ggplot2::aes(x = top_cor_coef, y = qval_plot)) +
    ggplot2::geom_point(color = point_color, size = plot_cfg$point_size, alpha = plot_cfg$alpha) +
    ggplot2::scale_y_continuous(limits = c(0, plot_cfg$y_limit)) +
    ggrepel::geom_text_repel(data = top_pairs,
                             ggplot2::aes(label = miR_TF),
                             force = 1,
                             direction = "both",
                             segment.size = 0.2,
                             color = "black",
                             nudge_x = 0.1,
                             nudge_y = 0.1,
                             max.overlaps = 1000,
                             size = plot_cfg$annot_text_size) +
    ggplot2::labs(x = x_label, y = expression(-log[10](qvalue))) +
    ggplot2::theme_minimal()
}

#' Stacked bar: miRNAs by correlation direction and TF-target overlap.
plot_mirna_distribution <- function(distribution_data, plot_cfg) {
  ggplot2::ggplot(distribution_data,
                  ggplot2::aes(x = correlation_type, y = value, fill = condition)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::geom_text(ggplot2::aes(label = value), color = "white",
                       position = ggplot2::position_stack(vjust = 0.5),
                       size = plot_cfg$annot_text_size) +
    ggplot2::theme_minimal() +
    ggplot2::scale_fill_manual(values = plot_cfg$colors$distribution) +
    ggplot2::scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 25)) +
    ggplot2::labs(x = "", y = "Number of miRs", fill = "")
}

#' Paired miRNA / TF z-scored expression, one line per significant pair.
#'
#' @param interactions_long Output of build_interaction_long_data().
plot_interactions <- function(interactions_long, plot_cfg, line_color) {
  mirna_labels <- interactions_long %>%
    dplyr::filter(variable == "miRNA") %>%
    dplyr::distinct(miRNA, .keep_all = TRUE)

  tf_labels <- interactions_long %>%
    dplyr::filter(variable == "TF") %>%
    dplyr::distinct(TF, .keep_all = TRUE)

  pd <- ggplot2::position_dodge(0.2)
  label_layer <- function(data, label_col) {
    ggrepel::geom_text_repel(data = data,
                             ggplot2::aes(label = .data[[label_col]]),
                             position = pd,
                             force = 1.5,
                             direction = "both",
                             segment.size = 0.2,
                             color = "black",
                             max.overlaps = 1000,
                             size = plot_cfg$annot_text_size,
                             min.segment.length = 0)
  }

  ggplot2::ggplot(interactions_long,
                  ggplot2::aes(x = variable, y = value, group = miR_TF)) +
    ggplot2::geom_point(ggplot2::aes(color = variable), position = pd, show.legend = FALSE,
                        size = plot_cfg$point_size, alpha = plot_cfg$alpha) +
    ggplot2::scale_color_manual(values = c(miRNA = plot_cfg$colors$mirna,
                                           TF = plot_cfg$colors$tf),
                                guide = "none") +
    label_layer(mirna_labels, "miRNA") +
    label_layer(tf_labels, "TF") +
    ggplot2::geom_line(color = line_color, alpha = 0.5, position = pd) +
    ggplot2::labs(x = NULL, y = "Expression Level \n(z-transformed)") +
    ggplot2::theme_minimal() +
    ggplot2::guides(alpha = "none")
}
