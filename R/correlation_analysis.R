# =============================================================================
# correlation_analysis.R
# mRNA-miRNA correlation (runAnalysis.Rmd section 4.1).
#
# The Rmd ran Hmisc::rcorr() on the stacked (genes + miRNAs) matrix, which
# computes every gene-gene and miRNA-miRNA pair too and then throws them away.
# Here only the gene x miRNA block is computed, with the same Pearson r and
# the same t-test p-value, so it is far faster and uses far less memory.
# =============================================================================

#' Correlate every gene with every miRNA and keep FDR-significant pairs.
#'
#' @param rna,mirna Expression data frames (features x samples), same samples.
#' @param gene_mapping tibble(gene_id, gene_symbol).
#' @return tibble(RNA, miRNA, coef, p_val, p_val_adj, RNA_ENSEMBL, RNA_symbol)
compute_mrna_mirna_correlation <- function(rna, mirna, gene_mapping,
                                           method = "pearson", fdr = 0.05,
                                           p_adjust_method = "BH") {
  stopifnot(setequal(colnames(rna), colnames(mirna)))
  rna <- as.matrix(rna)
  mirna <- as.matrix(mirna)[, colnames(rna), drop = FALSE]

  r <- stats::cor(t(rna), t(mirna), method = method)  # genes x miRNAs
  n <- ncol(rna)
  t_stat <- r * sqrt((n - 2) / pmax(1 - r^2, .Machine$double.eps))
  p <- 2 * stats::pt(-abs(t_stat), df = n - 2)
  p_adj <- stats::p.adjust(p, method = p_adjust_method)

  sig <- which(p_adj < fdr)          # p.adjust() drops dims: use linear indices
  rc <- arrayInd(sig, dim(r))

  tibble::tibble(
    RNA       = rownames(r)[rc[, 1]],
    miRNA     = colnames(r)[rc[, 2]],
    coef      = r[sig],
    p_val     = p[sig],
    p_val_adj = p_adj[sig]
  ) %>%
    dplyr::arrange(p_val) %>%
    dplyr::mutate(RNA_ENSEMBL = sub("\\..*", "", RNA)) %>%
    dplyr::left_join(gene_mapping, by = c("RNA" = "gene_id")) %>%
    dplyr::rename(RNA_symbol = gene_symbol)
}

#' Console summary of the significant correlations.
summarise_correlation <- function(cor_result) {
  list(
    n_pairs          = nrow(cor_result),
    n_mirnas         = dplyr::n_distinct(cor_result$miRNA),
    n_genes          = dplyr::n_distinct(cor_result$RNA_symbol),
    n_positive       = sum(cor_result$coef > 0),
    n_negative       = sum(cor_result$coef < 0),
    top_positive     = dplyr::slice_max(cor_result, coef, n = 1, with_ties = TRUE),
    top_negative     = dplyr::slice_min(cor_result, coef, n = 1, with_ties = TRUE)
  )
}

#' Z-scored mean expression of genes, labelled with the symbols used in the
#' correlation table (genes without significant correlations get NA symbols,
#' exactly as in the Rmd).
gene_expression_zscores <- function(rna, cor_result) {
  symbols <- dplyr::distinct(cor_result, RNA, RNA_symbol) %>%
    dplyr::distinct(RNA, .keep_all = TRUE)

  expression_zscores(rna, "RNA") %>%
    dplyr::left_join(symbols, by = "RNA")
}
