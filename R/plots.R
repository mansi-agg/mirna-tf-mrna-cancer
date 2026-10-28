# ggplot functions. Each takes ready-made plot data (see visualizations.R) and
# returns a ggplot.

POINT_SIZE <- 7
POINT_ALPHA <- 0.7
LABEL_SIZE <- 5

save_plot <- function(p, title, plots_dir, width = 8, height = 5) {
  dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)
  file_name <- paste0(gsub("\\s+", " ", gsub("[\r\n]", " ", title)), ".png")
  ggsave(plot = p, filename = file.path(plots_dir, file_name),
         width = width, height = height, bg = "white")
}

# TF target counts vs number of biological replicates (supplementary)
plot_targets_vs_replicates <- function(tf_target_count) {
  ggplot(tf_target_count, aes(x = log2(replicate.count), y = log2(target.count))) +
    geom_point() +
    labs(x = "Log2 Replicate Count",
         y = "Log2 Target Gene Count",
         title = "TF target counts are not correlated with the number of biological replicates") +
    theme_minimal()
}

# Significant miRNA-TF pairs: correlation of the top gene vs -log10(q-value).
# direction is "pos" or "neg".
plot_top_mir_tf_pairs <- function(ora_df, direction, n_labels = 5, y_limit = 10.5) {
  df <- as.data.frame(ora_df)
  if (direction == "neg") {
    df$top_cor_coef <- abs(df$top_cor_coef)
    point_color <- "#9C354B"
    x_label <- "Absolute Correlation Coefficient of \nTop Correlated Gene"
  } else {
    point_color <- "#0D84CD"
    x_label <- "Correlation Coefficient of \nTop Correlated Gene"
  }

  top_pairs <- df %>%
    dplyr::arrange(dplyr::desc(top_cor_coef), qval_plot) %>%
    head(n_labels)

  ggplot(df, aes(x = top_cor_coef, y = qval_plot)) +
    geom_point(color = point_color, size = POINT_SIZE, alpha = POINT_ALPHA) +
    scale_y_continuous(limits = c(0, y_limit)) +
    geom_text_repel(data = top_pairs, aes(label = miR_TF),
                    force = 1, direction = "both", segment.size = 0.2, color = "black",
                    nudge_x = 0.1, nudge_y = 0.1, max.overlaps = 1000, size = LABEL_SIZE) +
    labs(x = x_label, y = expression(-log[10](qvalue))) +
    theme_minimal()
}

# Stacked bar of miRNAs by correlation direction and TF-target overlap
plot_mirna_distribution <- function(distribution_data) {
  ggplot(distribution_data, aes(x = correlation_type, y = value, fill = condition)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = value), color = "white",
              position = position_stack(vjust = 0.5), size = LABEL_SIZE) +
    theme_minimal() +
    scale_fill_manual(values = c("#1a9850", "#4575b4", "#999999")) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 25)) +
    labs(x = "", y = "Number of miRs", fill = "")
}

# Z-scored miRNA and TF expression, one line per significant pair
plot_interactions <- function(interactions_long, line_color) {
  mirna_labels <- interactions_long %>%
    dplyr::filter(variable == "miRNA") %>%
    dplyr::distinct(miRNA, .keep_all = TRUE)
  tf_labels <- interactions_long %>%
    dplyr::filter(variable == "TF") %>%
    dplyr::distinct(TF, .keep_all = TRUE)

  pd <- position_dodge(0.2)
  label_layer <- function(data, label_col) {
    geom_text_repel(data = data, aes(label = .data[[label_col]]), position = pd,
                    force = 1.5, direction = "both", segment.size = 0.2, color = "black",
                    max.overlaps = 1000, size = LABEL_SIZE, min.segment.length = 0)
  }

  ggplot(interactions_long, aes(x = variable, y = value, group = miR_TF)) +
    geom_point(aes(color = variable), position = pd, show.legend = FALSE,
               size = POINT_SIZE, alpha = POINT_ALPHA) +
    scale_color_manual(values = c(miRNA = "#89b7d7", TF = "#ffbb78"), guide = "none") +
    label_layer(mirna_labels, "miRNA") +
    label_layer(tf_labels, "TF") +
    geom_line(color = line_color, alpha = 0.5, position = pd) +
    labs(x = NULL, y = "Expression Level \n(z-transformed)") +
    theme_minimal() +
    guides(alpha = "none")
}
