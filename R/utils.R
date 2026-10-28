# Small helpers used by the other scripts: logging, reading and saving files.

log_msg <- function(...) {
  message(format(Sys.time(), "[%H:%M:%S] "), ...)
}

read_table <- function(path) {
  suppressMessages(readr::read_csv(path, show_col_types = FALSE, progress = FALSE))
}

save_table <- function(df, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(as.data.frame(df), path, na = "NA")
  log_msg("Saved ", path)
}

# Expression matrices are saved with the gene / miRNA IDs in the first column
read_matrix <- function(path) {
  df <- as.data.frame(read_table(path))
  rownames(df) <- df[[1]]
  df[[1]] <- NULL
  df
}

save_matrix <- function(df, path) {
  save_table(tibble::rownames_to_column(as.data.frame(df), var = "id"), path)
}

# Label for a miRNA-TF pair, e.g. "miR-21-5p  -  STAT3"
make_mir_tf_label <- function(mirna, tf) {
  paste(sub("^hsa-", "", mirna), " - ", tf)
}
