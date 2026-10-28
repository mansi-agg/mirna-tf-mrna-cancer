# TF target genes from ChIP-seq experiments in cancer cell lines.

# Target genes of one TF with a peak in at least `threshold` of the cancer
# cell lines. Returns NULL if the TF has no data in these cell lines.
consensus_tf_targets <- function(tf_file, cell_line_ids, threshold) {
  tab <- read.table(tf_file, header = TRUE, sep = "\t")
  names(tab) <- gsub("\\..*", "", names(tab))  # "SRX186798.Breast_cancer" -> "SRX186798"

  cell_cols <- which(names(tab) %in% cell_line_ids)
  if (length(cell_cols) == 0) return(NULL)

  consensus <- rowMeans(as.matrix(tab[, cell_cols, drop = FALSE]) != 0)
  tab$Target_genes[consensus >= threshold]
}

# Consensus targets for every TF.
# Returns list(target_df = TF, TF_targets; target_count = TF, replicate.count, target.count)
build_tf_targets <- function(tf_files, cell_lines, threshold) {
  log_msg("Reading ", length(tf_files), " ChIP-Atlas files")
  targets <- lapply(tf_files, consensus_tf_targets,
                    cell_line_ids = cell_lines$`SRX ID`, threshold = threshold)
  targets <- Filter(Negate(is.null), targets)

  target_df <- tibble::tibble(TF         = rep(names(targets), lengths(targets)),
                              TF_targets = unlist(targets, use.names = FALSE))

  replicates <- dplyr::count(cell_lines, TF = `Track type`, name = "replicate.count")
  target_count <- tibble::tibble(TF = names(targets), target.count = lengths(targets)) %>%
    dplyr::left_join(replicates, by = "TF") %>%
    dplyr::select(TF, replicate.count, target.count)

  list(target_df = target_df, target_count = target_count)
}
