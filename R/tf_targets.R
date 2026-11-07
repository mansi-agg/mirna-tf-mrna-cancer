# TF target genes from ChIP-seq experiments in cancer cell lines.

# Target genes of one TF that are bound in at least `min_cell_lines` different
# cancer cell lines. Experiments are first merged per cell line: a gene counts
# as bound in a line if any of that line's experiments has a peak near it, so a
# heavily studied line (e.g. MCF-7) counts once. Returns NULL if the TF has no
# data in these cell lines.
consensus_tf_targets <- function(tf_file, experiment_line, min_cell_lines) {
  tab <- read.table(tf_file, header = TRUE, sep = "\t")
  names(tab) <- gsub("\\..*", "", names(tab))  # "SRX186798.Breast_cancer" -> "SRX186798"

  experiments <- intersect(names(tab), names(experiment_line))
  if (length(experiments) == 0) return(NULL)

  bound <- as.matrix(tab[, experiments, drop = FALSE]) != 0
  lines <- experiment_line[experiments]
  bound_per_line <- sapply(unique(lines), function(l) rowSums(bound[, lines == l, drop = FALSE]) > 0)
  n_lines_bound <- rowSums(matrix(bound_per_line, nrow = nrow(tab)))

  tab$Target_genes[n_lines_bound >= min_cell_lines]
}

# Consensus targets for every TF.
# Returns list(target_df = TF, TF_targets; target_count = TF, replicate.count, target.count)
build_tf_targets <- function(tf_files, cell_lines, min_cell_lines) {
  log_msg("Reading ", length(tf_files), " ChIP-Atlas files")
  # Cell line of each experiment (Cellosaurus ID; older lists only have the cell type name)
  line_col <- if ("cellosaurus_id" %in% names(cell_lines)) "cellosaurus_id" else "Cell type"
  experiment_line <- setNames(cell_lines[[line_col]], cell_lines$`SRX ID`)
  targets <- lapply(tf_files, consensus_tf_targets,
                    experiment_line = experiment_line, min_cell_lines = min_cell_lines)
  targets <- Filter(Negate(is.null), targets)

  target_df <- tibble::tibble(TF         = rep(names(targets), lengths(targets)),
                              TF_targets = unlist(targets, use.names = FALSE))

  replicates <- dplyr::count(cell_lines, TF = `Track type`, name = "replicate.count")
  target_count <- tibble::tibble(TF = names(targets), target.count = lengths(targets)) %>%
    dplyr::left_join(replicates, by = "TF") %>%
    dplyr::select(TF, replicate.count, target.count)

  list(target_df = target_df, target_count = target_count)
}
