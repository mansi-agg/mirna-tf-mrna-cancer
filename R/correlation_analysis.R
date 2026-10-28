# mRNA-miRNA correlation.

# Correlate every gene with every miRNA and keep the FDR-significant pairs.
# Gives the same Pearson r and t-test p-value as Hmisc::rcorr(), but only
# computes the gene x miRNA pairs, so it is much faster.
compute_mrna_mirna_correlation <- function(rna, mirna, gene_mapping,
                                           method = "pearson", fdr = 0.05) {
  rna <- as.matrix(rna)
  mirna <- as.matrix(mirna)[, colnames(rna)]

  r <- cor(t(rna), t(mirna), method = method)  # genes x miRNAs
  n <- ncol(rna)
  t_stat <- r * sqrt((n - 2) / pmax(1 - r^2, .Machine$double.eps))
  p <- 2 * pt(-abs(t_stat), df = n - 2)
  p_adj <- p.adjust(p, method = "BH")

  sig <- which(p_adj < fdr)
  rows_cols <- arrayInd(sig, dim(r))

  tibble::tibble(RNA       = rownames(r)[rows_cols[, 1]],
                 miRNA     = colnames(r)[rows_cols[, 2]],
                 coef      = r[sig],
                 p_val     = p[sig],
                 p_val_adj = p_adj[sig]) %>%
    dplyr::arrange(p_val) %>%
    dplyr::mutate(RNA_ENSEMBL = sub("\\..*", "", RNA)) %>%
    dplyr::left_join(gene_mapping, by = c("RNA" = "gene_id")) %>%
    dplyr::rename(RNA_symbol = gene_symbol)
}

# Z-scored mean gene expression, with the gene symbols from the correlation table
gene_expression_zscores <- function(rna, cor_result) {
  symbols <- dplyr::distinct(cor_result, RNA, RNA_symbol) %>%
    dplyr::distinct(RNA, .keep_all = TRUE)

  expression_zscores(rna, "RNA") %>%
    dplyr::left_join(symbols, by = "RNA")
}
