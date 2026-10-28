# =============================================================================
# data_preparation.R
# Download / load TCGA RNA-seq and miRNA-seq data and turn them into filtered,
# log2-transformed, sample-matched expression matrices.
# (runAnalysis.Rmd sections 2 and 3)
# =============================================================================

#' Download TCGA RNA-seq and miRNA-seq data from the GDC and save as RDS.
#' Only runs when "tcga_data" is in `regenerate`; this is a large download.
download_tcga_data <- function(study, rna_rds_path, mirna_rds_path,
                               download_dir = file.path(dirname(rna_rds_path), "GDCdata")) {
  if (!requireNamespace("TCGAbiolinks", quietly = TRUE)) {
    stop("TCGAbiolinks is required to download TCGA data.", call. = FALSE)
  }

  query_and_prepare <- function(data_type, strategy, sub_dir) {
    query <- TCGAbiolinks::GDCquery(project = study,
                                    data.category = "Transcriptome Profiling",
                                    data.type = data_type,
                                    experimental.strategy = strategy)
    dir <- file.path(download_dir, sub_dir)
    TCGAbiolinks::GDCdownload(query = query, files.per.chunk = 10, directory = dir)
    TCGAbiolinks::GDCprepare(query = query, directory = dir)
  }

  rna_data <- query_and_prepare("Gene Expression Quantification", "RNA-Seq", "RNA")
  saveRDS(rna_data, rna_rds_path)

  mirna_data <- query_and_prepare("miRNA Expression Quantification", "miRNA-Seq", "miRNA")
  saveRDS(mirna_data, mirna_rds_path)

  list(rna = rna_data, mirna = mirna_data)
}

#' Protein-coding TPM matrix (genes x samples) from the TCGA SummarizedExperiment.
extract_rna_matrix <- function(rna_se, gene_type = "protein_coding", assay_name = "tpm_unstrand") {
  keep <- which(SummarizedExperiment::rowData(rna_se)$gene_type == gene_type)
  rna_se <- rna_se[keep, ]

  mat <- SummarizedExperiment::assay(rna_se, assay_name)
  colnames(mat) <- rna_se$sample
  rownames(mat) <- rownames(rna_se)
  mat
}

#' Gene ID -> gene symbol lookup. Kept separately because several gene IDs
#' share a symbol, so symbols cannot be used as row names.
extract_gene_mapping <- function(rna_se, gene_type = "protein_coding") {
  row_data <- SummarizedExperiment::rowData(rna_se)
  keep <- row_data$gene_type == gene_type
  tibble::tibble(gene_id     = rownames(rna_se)[keep],
                 gene_symbol = row_data$gene_name[keep])
}

#' miRNA RPM matrix (miRNAs x samples) with sample names matching the RNA data.
extract_mirna_matrix <- function(mirna_df, value_prefix = "reads_per_million_miRNA_mapped_") {
  value_cols <- grep(value_prefix, colnames(mirna_df), fixed = TRUE, value = TRUE)

  mat <- as.matrix(mirna_df[, value_cols, drop = FALSE])
  rownames(mat) <- mirna_df$miRNA_ID

  # Drop the last three "-"-separated fields of the barcode so miRNA sample
  # IDs match the RNA sample IDs.
  colnames(mat) <- sub("(-[^-]*){3}$", "", sub(value_prefix, "", value_cols, fixed = TRUE))
  mat
}

#' Remove all-zero rows and rows with more than `max_zero_pct` % zeros.
filter_low_expression <- function(mat, max_zero_pct) {
  zero_pct <- rowMeans(mat == 0, na.rm = TRUE) * 100
  keep <- rowSums(mat) > 0 & zero_pct <= max_zero_pct
  mat[keep, , drop = FALSE]
}

#' Full preparation: match samples, filter, log2(x + 1).
#'
#' @return list(rna = data.frame, mirna = data.frame), genes/miRNAs in rows.
prepare_expression_data <- function(rna_se, mirna_df, filters, verbose = TRUE) {
  rna <- extract_rna_matrix(rna_se, filters$gene_type, filters$rna_assay)
  mirna <- extract_mirna_matrix(mirna_df, filters$mirna_value_prefix)
  log_msg("Raw data: ", nrow(rna), " genes x ", ncol(rna), " samples; ",
          nrow(mirna), " miRNAs x ", ncol(mirna), " samples", verbose = verbose)

  common <- intersect(colnames(rna), colnames(mirna))
  rna <- rna[, common, drop = FALSE]
  mirna <- mirna[, common, drop = FALSE]
  log_msg(length(common), " samples have matched mRNA & miRNA data", verbose = verbose)

  rna <- log2(filter_low_expression(rna, filters$expression_threshold) + 1)
  mirna <- log2(filter_low_expression(mirna, filters$expression_threshold) + 1)
  log_msg("After filtering: ", nrow(rna), " genes, ", nrow(mirna), " miRNAs", verbose = verbose)

  list(rna = as.data.frame(rna), mirna = as.data.frame(mirna))
}

#' Mean expression per row, z-transformed across rows.
expression_zscores <- function(expr, id_col) {
  tibble::tibble(!!id_col := rownames(expr),
                 z_transform_exp = as.numeric(scale(rowMeans(as.matrix(expr)))))
}
