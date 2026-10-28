# =============================================================================
# tf_targets.R
# "True" TF target genes from ChIP-seq in cancer cell lines
# (runAnalysis.Rmd section 4.5).
# =============================================================================

#' Cancer cell-line ChIP-seq experiments (one row per SRX experiment).
load_cell_lines <- function(path) {
  read_table_csv(path)
}

#' Number of cell-line experiments per TF ("biological replicates").
count_replicates <- function(cell_lines) {
  cell_lines %>%
    dplyr::count(TF = `Track type`, name = "replicate.count")
}

#' Target genes of one TF that have a peak in at least `threshold` of the
#' cancer cell lines. Returns NULL when the TF has no cancer cell-line data.
consensus_tf_targets <- function(tf_file, cancer_cell_lines, threshold) {
  tab <- utils::read.table(tf_file, header = TRUE, sep = "\t")
  # Column names look like "SRX186798.Breast_cancer" -> keep the SRX ID
  names(tab) <- gsub("\\..*", "", names(tab))

  cell_cols <- which(names(tab) %in% cancer_cell_lines)
  if (!length(cell_cols)) return(NULL)

  consensus <- rowMeans(as.matrix(tab[, cell_cols, drop = FALSE]) != 0)
  tab$Target_genes[consensus >= threshold]
}

#' Consensus targets for every TF with a ChIP-Atlas file.
#'
#' @param tf_files Named vector TF -> file path (from list_tf_target_files()).
#' @return list(target_df = tibble(TF, TF_targets),
#'              target_count = tibble(TF, replicate.count, target.count))
build_tf_targets <- function(tf_files, cell_lines, threshold, verbose = TRUE) {
  cancer_cell_lines <- cell_lines$`SRX ID`

  log_msg("Reading ", length(tf_files), " ChIP-Atlas target files", verbose = verbose)
  targets <- lapply(tf_files, consensus_tf_targets,
                    cancer_cell_lines = cancer_cell_lines, threshold = threshold)
  targets <- targets[!vapply(targets, is.null, logical(1))]

  target_df <- tibble::tibble(
    TF         = rep(names(targets), lengths(targets)),
    TF_targets = unlist(targets, use.names = FALSE)
  )

  target_count <- tibble::tibble(TF = names(targets), target.count = lengths(targets)) %>%
    dplyr::left_join(count_replicates(cell_lines), by = "TF") %>%
    dplyr::select(TF, replicate.count, target.count)

  list(target_df = target_df, target_count = target_count)
}

#' Spearman test: are target counts driven by the number of replicates?
test_replicate_bias <- function(target_count) {
  suppressWarnings(stats::cor.test(target_count$replicate.count,
                                   target_count$target.count,
                                   method = "spearman"))
}

#' Console summary of TF / cell-line coverage.
summarise_tf_coverage <- function(overlap, cell_lines, study) {
  tfs_with_cell_lines <- unique(cell_lines$`Track type`)
  list(
    study                         = study,
    n_tfs_with_chip_and_mirna     = dplyr::n_distinct(overlap$TF),
    n_tfs_with_cancer_cell_lines  = length(tfs_with_cell_lines),
    n_tfs_in_both                 = length(intersect(overlap$TF, tfs_with_cell_lines))
  )
}
