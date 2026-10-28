# Which miRNAs target known transcription factors (TFs), and ChIP-Atlas TF target files.

# Known human TFs (Lambert et al. 2018)
load_known_tfs <- function(path) {
  read.csv(path, row.names = 1) %>%
    dplyr::filter(Is.TF. == "Yes")
}

# Validated targets of the given miRNAs from multiMiR (needs internet)
query_mirna_targets <- function(mirnas) {
  res <- multiMiR::get_multimir(org = "hsa", mirna = mirnas, table = "validated",
                                summary = TRUE, use.tibble = TRUE)
  as.data.frame(res@data) %>%
    dplyr::mutate(miRNA = tolower(mature_mirna_id))
}

# Validated miRNA -> known TF interactions
query_mirna_tf_interactions <- function(mirnas, known_tfs) {
  query_mirna_targets(mirnas) %>%
    dplyr::filter(target_symbol %in% known_tfs$HGNC.symbol)
}

# Read a saved multiMiR table (older files have no lower-case miRNA column)
read_interactions <- function(path) {
  df <- read_table(path)
  if (!"miRNA" %in% names(df)) df <- dplyr::mutate(df, miRNA = tolower(mature_mirna_id))
  df
}

# ---- ChIP-Atlas TF target files ----

chip_atlas_suffix <- function(genome, tss_kb) {
  paste0("_", genome, "_TSS", tss_kb, "kb.tsv")
}

# Download ChIP-Atlas target gene tables (https://chip-atlas.org/target_genes)
download_chip_atlas_files <- function(tfs, dest_dir, genome = "hg38", tss_kb = 5) {
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)
  failed <- character(0)

  for (tf in unique(tfs)) {
    dest <- file.path(dest_dir, paste0(tf, chip_atlas_suffix(genome, tss_kb)))
    url <- sprintf("https://chip-atlas.dbcls.jp/data/%s/target/%s.%s.tsv", genome, tf, tss_kb)
    ok <- tryCatch({
      download.file(url, destfile = dest, quiet = TRUE)
      TRUE
    }, error = function(e) FALSE, warning = function(w) FALSE)

    if (!ok) {
      unlink(dest)
      failed <- c(failed, tf)
    }
  }
  if (length(failed)) log_msg("No ChIP-Atlas data for ", length(failed), " TFs: ",
                              paste(failed, collapse = ", "))
}

# Named vector TF -> path of each ChIP-Atlas file in `dir`
list_chip_atlas_files <- function(dir, genome = "hg38", tss_kb = 5) {
  suffix <- chip_atlas_suffix(genome, tss_kb)
  files <- list.files(dir, full.names = TRUE)
  files <- files[endsWith(files, suffix)]
  setNames(files, sub(suffix, "", basename(files), fixed = TRUE))
}

# ---- miRNA x TF overlap ----

# One row per miRNA-TF pair where the miRNA has significant gene correlations,
# targets the TF, and the TF has ChIP-Atlas data
build_mirna_tf_overlap <- function(cor_result, interactions_tfs, tfs_with_chip_data) {
  mir_tf <- interactions_tfs %>%
    dplyr::select(miRNA, target_symbol) %>%
    dplyr::filter(target_symbol %in% tfs_with_chip_data) %>%
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
