# =============================================================================
# tf_interactions.R
# Which miRNAs target known transcription factors (runAnalysis.Rmd 4.2 - 4.4).
# =============================================================================

#' Known human TFs (Lambert et al. 2018, http://humantfs.ccbr.utoronto.ca).
load_known_tfs <- function(path) {
  utils::read.csv(path, row.names = 1) %>%
    dplyr::filter(Is.TF. == "Yes")
}

#' Validated miRNA targets from multiMiR, with a lower-case `miRNA` column
#' matching the TCGA miRNA IDs.
query_mirna_targets <- function(mirnas, table = "validated") {
  res <- multiMiR::get_multimir(org = "hsa",
                                mirna = mirnas,
                                table = table,
                                summary = TRUE,
                                use.tibble = TRUE)
  as.data.frame(res@data) %>%
    dplyr::mutate(miRNA = tolower(mature_mirna_id))
}

#' Validated miRNA -> known TF interactions.
query_mirna_tf_interactions <- function(mirnas, known_tfs, table = "validated") {
  query_mirna_targets(mirnas, table) %>%
    dplyr::filter(target_symbol %in% known_tfs$HGNC.symbol)
}

#' Add the lower-case `miRNA` column if a saved file predates it.
standardise_interactions <- function(df) {
  if (!"miRNA" %in% names(df)) df <- dplyr::mutate(df, miRNA = tolower(mature_mirna_id))
  df
}

# -----------------------------------------------------------------------------
# ChIP-Atlas TF target files
# -----------------------------------------------------------------------------

tf_target_file_suffix <- function(chip_atlas) {
  paste0("_", chip_atlas$genome, "_TSS", chip_atlas$peak_distance_from_tss, "kb.tsv")
}

#' Download ChIP-Atlas target-gene tables (https://chip-atlas.org/target_genes).
#' Only runs when "chip_atlas" is in `regenerate`.
download_tf_target_files <- function(tfs, dest_dir, chip_atlas, overwrite = FALSE, verbose = TRUE) {
  ensure_dir(dest_dir)
  suffix <- tf_target_file_suffix(chip_atlas)
  failed <- character(0)

  for (tf in unique(tfs)) {
    dest <- file.path(dest_dir, paste0(tf, suffix))
    if (file.exists(dest) && !overwrite) next

    url <- sprintf("%s/%s/target/%s.%s.tsv", chip_atlas$base_url, chip_atlas$genome,
                   tf, chip_atlas$peak_distance_from_tss)
    ok <- tryCatch({
      utils::download.file(url, destfile = dest, quiet = TRUE)
      TRUE
    }, error = function(e) FALSE, warning = function(w) FALSE)

    if (!ok) {
      if (file.exists(dest)) unlink(dest)
      failed <- c(failed, tf)
    }
  }
  if (length(failed)) log_msg("No ChIP-Atlas data for ", length(failed), " TFs: ",
                              paste(failed, collapse = ", "), verbose = verbose)
  invisible(failed)
}

#' Named vector TF -> path of every available ChIP-Atlas target file.
list_tf_target_files <- function(dir, chip_atlas) {
  suffix <- tf_target_file_suffix(chip_atlas)
  files <- list.files(dir, pattern = paste0(gsub(".", "\\.", suffix, fixed = TRUE), "$"),
                      full.names = TRUE)
  stats::setNames(files, sub(suffix, "", basename(files), fixed = TRUE))
}

# -----------------------------------------------------------------------------
# miRNA x TF overlap table
# -----------------------------------------------------------------------------

#' One row per (miRNA, TF) pair where the miRNA has significant gene
#' correlations, targets the TF, and the TF has ChIP-Atlas target data.
build_mirna_tf_overlap <- function(cor_result, interactions_tfs, tf_available) {
  mir_tf <- interactions_tfs %>%
    dplyr::select(miRNA, target_symbol) %>%
    dplyr::filter(target_symbol %in% tf_available) %>%
    dplyr::distinct()

  cor_counts <- cor_result %>%
    dplyr::group_by(miRNA) %>%
    dplyr::summarise(positively_cor_genes = sum(coef > 0),
                     negatively_cor_genes = sum(coef < 0),
                     .groups = "drop")

  dplyr::inner_join(cor_counts, mir_tf, by = "miRNA") %>%
    dplyr::rename(TF = target_symbol) %>%
    dplyr::add_count(miRNA, name = "num_TF_interactions") %>%
    dplyr::relocate(num_TF_interactions, .before = TF)
}
