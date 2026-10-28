# =============================================================================
# utils.R
# Small helpers shared by all modules: logging, file I/O and labels.
# =============================================================================

log_msg <- function(..., verbose = TRUE) {
  if (isTRUE(verbose)) {
    message(sprintf("[%s] %s", format(Sys.time(), "%H:%M:%S"), paste0(...)))
  }
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

read_table_csv <- function(path) {
  suppressMessages(readr::read_csv(path, show_col_types = FALSE, progress = FALSE))
}

write_table_csv <- function(df, path) {
  ensure_dir(dirname(path))
  readr::write_csv(as.data.frame(df), path, na = "NA")
  invisible(path)
}

#' Read a CSV whose first column holds row names (e.g. expression matrices).
read_matrix_csv <- function(path) {
  df <- suppressMessages(readr::read_csv(path, show_col_types = FALSE, progress = FALSE))
  df <- as.data.frame(df)
  rownames(df) <- df[[1]]
  df[[1]] <- NULL
  df
}

write_matrix_csv <- function(df, path, id_col = "id") {
  ensure_dir(dirname(path))
  out <- tibble::rownames_to_column(as.data.frame(df), var = id_col)
  readr::write_csv(out, path)
  invisible(path)
}

#' Label used for every miRNA-TF pair (kept identical across all tables so
#' joins between ORA results and validated interactions line up).
make_mir_tf_label <- function(mirna, tf) {
  paste(sub("^hsa-", "", mirna), " - ", tf)
}

#' Turn a plot title into a safe file name (no new lines / repeated spaces).
title_to_filename <- function(title, ext = ".png") {
  name <- gsub("\\s+", " ", gsub("[\r\n]", " ", title))
  name <- gsub("[/\\\\:*?\"<>|]", "_", trimws(name))
  paste0(name, ext)
}

`%||%` <- function(a, b) if (is.null(a)) b else a
