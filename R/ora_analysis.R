# =============================================================================
# ora_analysis.R
# Over-representation of TF targets among miRNA-correlated genes
# (runAnalysis.Rmd section 5 and get_top_cor_stats() from section 6).
#
# For each miRNA, the miRNA's own validated direct targets are removed from its
# positively / negatively correlated genes, and each set is tested against the
# TF target gene sets with clusterProfiler::enricher().
# =============================================================================

#' TERM2GENE table: TF target sets restricted to the expressed background,
#' plus the background itself as a term (as in the Rmd).
build_term2gene <- function(tf_target_df, background_genes) {
  tf_target_df %>%
    dplyr::filter(TF_targets %in% background_genes) %>%
    dplyr::bind_rows(tibble::tibble(TF = "background", TF_targets = background_genes))
}

#' Run enricher() on one gene set and return a tidy result with prefixed columns.
run_enricher <- function(genes, mirna, prefix, term2gene, universe) {
  empty <- tibble::tibble(miRNA = character(), ID = character(),
                          TF_target_BgRatio = character(),
                          !!paste0(prefix, "_GeneRatio") := character(),
                          !!paste0(prefix, "_overlap_geneID") := character(),
                          !!paste0(prefix, "_ora_pval") := numeric(),
                          !!paste0(prefix, "_ora_qval") := numeric())
  if (!length(genes)) return(empty)

  res <- suppressMessages(clusterProfiler::enricher(
    gene          = genes,
    pvalueCutoff  = 1,
    pAdjustMethod = "BH",
    qvalueCutoff  = 1,
    universe      = universe,
    TERM2GENE     = term2gene,
    minGSSize     = 0,
    maxGSSize     = length(universe) - 2
  ))
  if (is.null(res) || !nrow(res@result)) return(empty)

  tibble::as_tibble(res@result) %>%
    dplyr::transmute(miRNA = mirna, ID,
                     TF_target_BgRatio = BgRatio,
                     !!paste0(prefix, "_GeneRatio") := GeneRatio,
                     !!paste0(prefix, "_overlap_geneID") := geneID,
                     !!paste0(prefix, "_ora_pval") := pvalue,
                     !!paste0(prefix, "_ora_qval") := qvalue)
}

#' ORA for one miRNA, positive and negative correlations.
ora_for_mirna <- function(mirna, cor_result, direct_targets, term2gene, universe) {
  excluded <- direct_targets$target_symbol[direct_targets$miRNA == mirna]
  genes <- cor_result %>%
    dplyr::filter(miRNA == mirna, !(RNA_symbol %in% excluded))

  list(
    poscor = run_enricher(genes$RNA_symbol[genes$coef > 0], mirna, "poscor", term2gene, universe),
    negcor = run_enricher(genes$RNA_symbol[genes$coef < 0], mirna, "negcor", term2gene, universe)
  )
}

#' Per miRNA-TF pair: the most extreme and mean correlation of the overlap genes.
add_top_correlation_stats <- function(ora_df, cor_result, prefix) {
  gene_col <- paste0(prefix, "_overlap_geneID")
  pick_idx <- if (prefix == "poscor") which.max else which.min

  stats <- ora_df %>%
    dplyr::distinct(miRNA, TF, gene = .data[[gene_col]]) %>%
    tidyr::separate_rows(gene, sep = "/") %>%
    dplyr::inner_join(dplyr::select(cor_result, miRNA, RNA_symbol, coef),
                      by = c("miRNA", "gene" = "RNA_symbol"),
                      relationship = "many-to-many") %>%
    dplyr::group_by(miRNA, TF) %>%
    dplyr::summarise(top_cor_gene = gene[pick_idx(coef)],
                     top_cor_coef = coef[pick_idx(coef)],
                     avg_cor_coef = mean(coef),
                     .groups = "drop")

  dplyr::left_join(ora_df, stats, by = c("miRNA", "TF"))
}

#' Full ORA for all miRNA-TF pairs.
#'
#' @return list(poscor = tibble, negcor = tibble) of significant pairs, with
#'   top-correlation stats, -log10 q-value and miR_TF label.
run_ora_analysis <- function(overlap, cor_result, tf_target_df, direct_targets,
                             background_genes, fdr = 0.05, p_adjust_method = "BH",
                             cores = 1L, verbose = TRUE) {
  # Evaluate lazy arguments once here, not separately inside each worker
  force(cor_result); force(direct_targets)
  term2gene <- build_term2gene(tf_target_df, background_genes)
  overlap <- dplyr::filter(overlap, TF %in% tf_target_df$TF)
  mirnas <- unique(overlap$miRNA)

  log_msg("Running ORA for ", length(mirnas), " miRNAs on ", cores, " core(s)", verbose = verbose)
  cor_result <- dplyr::filter(cor_result, miRNA %in% mirnas)
  run_one <- function(m) ora_for_mirna(m, cor_result, direct_targets, term2gene, background_genes)
  per_mirna <- if (cores > 1) {
    parallel::mclapply(mirnas, run_one, mc.cores = cores)
  } else {
    lapply(mirnas, run_one)
  }

  # Attach ORA results to the miRNA-TF pairs
  results <- lapply(c(poscor = "poscor", negcor = "negcor"), function(prefix) {
    ora <- dplyr::bind_rows(lapply(per_mirna, `[[`, prefix))
    overlap %>%
      dplyr::left_join(ora, by = c("miRNA", "TF" = "ID")) %>%
      dplyr::filter(!is.na(.data[[paste0(prefix, "_ora_pval")]]))
  })

  # Multiple-testing correction across both directions together
  p_all <- c(results$poscor$poscor_ora_pval, results$negcor$negcor_ora_pval)
  p_adj <- stats::p.adjust(p_all, method = p_adjust_method)
  n_pos <- nrow(results$poscor)
  results$poscor$poscor_ora_pval.adj <- p_adj[seq_len(n_pos)]
  results$negcor$negcor_ora_pval.adj <- p_adj[n_pos + seq_len(nrow(results$negcor))]

  epsilon <- 1e-10  # avoid log(0)
  lapply(stats::setNames(names(results), names(results)), function(prefix) {
    adj_col <- paste0(prefix, "_ora_pval.adj")
    results[[prefix]] %>%
      dplyr::filter(.data[[adj_col]] < fdr) %>%
      add_top_correlation_stats(cor_result, prefix) %>%
      dplyr::mutate(qval_plot = -log10(.data[[adj_col]] + epsilon),
                    miR_TF = make_mir_tf_label(miRNA, TF))
  })
}
