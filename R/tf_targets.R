# TF target genes from ChIP-seq experiments in cancer cell lines.

# Target genes of one TF, using a two-level rule so that each cell line counts once:
#   1. Within a cell line: the gene is "bound in that line" if it has a peak in at
#      least `line_threshold` of that line's experiments (replicates agree).
#   2. Across cell lines: the gene is a target if it is bound in at least
#      `across_threshold` of the TF's cell lines.
# Example: MCF-7 bound in 460/529 experiments (yes), T-47D in 30/108 (no),
# ZR-75-1 in 5/5 (yes) -> bound in 2 of 3 lines (67%) -> target.
# Returns NULL if the TF has no experiments in these cell lines.
consensus_tf_targets <- function(tf_file, experiment_line, line_threshold, across_threshold) {
  tab <- read.table(tf_file, header = TRUE, sep = "\t")
  names(tab) <- gsub("\\..*", "", names(tab))  # "SRX186798.Breast_cancer" -> "SRX186798"

  experiments <- intersect(names(tab), names(experiment_line))
  if (length(experiments) == 0) return(NULL)

  bound <- as.matrix(tab[, experiments, drop = FALSE]) != 0
  lines <- experiment_line[experiments]

  # Level 1: fraction of each line's experiments with a peak (genes x lines)
  line_fraction <- sapply(unique(lines), function(l) rowMeans(bound[, lines == l, drop = FALSE]))
  bound_in_line <- matrix(line_fraction, nrow = nrow(tab)) >= line_threshold

  # Level 2: fraction of lines in which the gene is bound
  tab$Target_genes[rowMeans(bound_in_line) >= across_threshold]
}

# Consensus targets for every TF.
# Returns list(target_df = TF, TF_targets; target_count = TF, replicate.count, target.count)
build_tf_targets <- function(tf_files, cell_lines, line_threshold, across_threshold) {
  log_msg("Reading ", length(tf_files), " ChIP-Atlas files")
  # Cell line of each experiment (Cellosaurus ID; older lists only have the cell type name)
  line_col <- if ("cellosaurus_id" %in% names(cell_lines)) "cellosaurus_id" else "Cell type"
  experiment_line <- setNames(cell_lines[[line_col]], cell_lines$`SRX ID`)
  targets <- lapply(tf_files, consensus_tf_targets,
                    experiment_line = experiment_line,
                    line_threshold = line_threshold, across_threshold = across_threshold)
  targets <- Filter(Negate(is.null), targets)

  target_df <- tibble::tibble(TF         = rep(names(targets), lengths(targets)),
                              TF_targets = unlist(targets, use.names = FALSE))

  # Per TF: experiments (replicate.count) and distinct cell lines in the list
  replicates <- cell_lines %>%
    dplyr::group_by(TF = `Track type`) %>%
    dplyr::summarise(replicate.count = dplyr::n(),
                     cell_line.count = dplyr::n_distinct(.data[[line_col]]), .groups = "drop")
  target_count <- tibble::tibble(TF = names(targets), target.count = lengths(targets)) %>%
    dplyr::left_join(replicates, by = "TF") %>%
    dplyr::select(TF, replicate.count, cell_line.count, target.count)

  list(target_df = target_df, target_count = target_count)
}
