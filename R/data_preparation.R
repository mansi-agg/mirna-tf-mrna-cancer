# Download TCGA data and turn it into filtered, log2 expression matrices.

# Download RNA-seq and miRNA-seq data from the GDC and save them as RDS files.
# This is a large download (several GB).
download_tcga_data <- function(study, rna_rds_path, mirna_rds_path) {
  download_dir <- file.path(dirname(rna_rds_path), "GDCdata")

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
  log_msg("Saved ", rna_rds_path)

  mirna_data <- query_and_prepare("miRNA Expression Quantification", "miRNA-Seq", "miRNA")
  saveRDS(mirna_data, mirna_rds_path)
  log_msg("Saved ", mirna_rds_path)
}

# Gene ID -> gene symbol for protein-coding genes
get_gene_mapping <- function(rna_data) {
  row_data <- SummarizedExperiment::rowData(rna_data)
  keep <- row_data$gene_type == "protein_coding"
  tibble::tibble(gene_id     = rownames(rna_data)[keep],
                 gene_symbol = row_data$gene_name[keep])
}

# Read the saved gene mapping, or make it from the RNA data if it isn't saved yet
load_gene_mapping <- function(gene_mapping_path, rna_rds_path) {
  if (!file.exists(gene_mapping_path)) {
    save_table(get_gene_mapping(readRDS(rna_rds_path)), gene_mapping_path)
  }
  read_table(gene_mapping_path)
}

# Protein-coding TPM matrix (genes x samples)
get_rna_matrix <- function(rna_data) {
  rna_data <- rna_data[SummarizedExperiment::rowData(rna_data)$gene_type == "protein_coding", ]
  mat <- SummarizedExperiment::assay(rna_data, "tpm_unstrand")
  colnames(mat) <- rna_data$sample
  mat
}

# miRNA reads-per-million matrix (miRNAs x samples), with sample IDs shortened
# to match the RNA sample IDs
get_mirna_matrix <- function(mirna_data) {
  prefix <- "reads_per_million_miRNA_mapped_"
  value_cols <- grep(prefix, colnames(mirna_data), fixed = TRUE, value = TRUE)

  mat <- as.matrix(mirna_data[, value_cols])
  rownames(mat) <- mirna_data$miRNA_ID
  colnames(mat) <- sub("(-[^-]*){3}$", "", sub(prefix, "", value_cols, fixed = TRUE))
  mat
}

# Remove rows that are all zero or zero in more than `max_zero_pct` % of samples
filter_low_expression <- function(mat, max_zero_pct) {
  zero_pct <- rowMeans(mat == 0, na.rm = TRUE) * 100
  mat[rowSums(mat) > 0 & zero_pct <= max_zero_pct, , drop = FALSE]
}

# Keep samples with both RNA and miRNA data, filter, and log2(x + 1) transform.
# Returns list(rna, mirna) with genes / miRNAs in rows.
prepare_expression_data <- function(rna_data, mirna_data, max_zero_pct) {
  rna <- get_rna_matrix(rna_data)
  mirna <- get_mirna_matrix(mirna_data)

  samples <- intersect(colnames(rna), colnames(mirna))
  log_msg(length(samples), " samples have both mRNA and miRNA data")

  rna <- log2(filter_low_expression(rna[, samples], max_zero_pct) + 1)
  mirna <- log2(filter_low_expression(mirna[, samples], max_zero_pct) + 1)
  log_msg("After filtering: ", nrow(rna), " genes, ", nrow(mirna), " miRNAs")

  list(rna = as.data.frame(rna), mirna = as.data.frame(mirna))
}

# Mean expression of each row, z-transformed across rows
expression_zscores <- function(expr, id_col) {
  tibble::tibble(!!id_col := rownames(expr),
                 z_transform_exp = as.numeric(scale(rowMeans(as.matrix(expr)))))
}
