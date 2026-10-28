# Over-representation analysis (ORA): are a miRNA's correlated genes enriched
# for the targets of the TFs it regulates?
#
# For each miRNA, its own validated direct targets are removed from its
# positively / negatively correlated genes, and each gene set is tested against
# the TF target gene sets with clusterProfiler::enricher().

# Run enricher() on one gene set; columns are prefixed with "poscor" or "negcor"
run_enricher <- function(genes, mirna, prefix, term2gene, universe) {
  empty <- tibble::tibble(miRNA = character(), ID = character(), TF_target_BgRatio = character(),
                          !!paste0(prefix, "_GeneRatio") := character(),
                          !!paste0(prefix, "_overlap_geneID") := character(),
                          !!paste0(prefix, "_ora_pval") := numeric(),
                          !!paste0(prefix, "_ora_qval") := numeric())
  if (length(genes) == 0) return(empty)

  res <- suppressMessages(clusterProfiler::enricher(
    gene = genes, universe = universe, TERM2GENE = term2gene,
    pvalueCutoff = 1, qvalueCutoff = 1, pAdjustMethod = "BH",
    minGSSize = 0, maxGSSize = length(universe) - 2
  ))
  if (is.null(res) || nrow(res@result) == 0) return(empty)

  tibble::as_tibble(res@result) %>%
    dplyr::transmute(miRNA = mirna, ID,
                     TF_target_BgRatio = BgRatio,
                     !!paste0(prefix, "_GeneRatio") := GeneRatio,
                     !!paste0(prefix, "_overlap_geneID") := geneID,
                     !!paste0(prefix, "_ora_pval") := pvalue,
                     !!paste0(prefix, "_ora_qval") := qvalue)
}

# ORA for one miRNA, for its positively and negatively correlated genes
ora_for_mirna <- function(mirna, cor_result, direct_targets, term2gene, universe) {
  excluded <- direct_targets$target_symbol[direct_targets$miRNA == mirna]
  genes <- dplyr::filter(cor_result, miRNA == mirna, !(RNA_symbol %in% excluded))

  list(poscor = run_enricher(genes$RNA_symbol[genes$coef > 0], mirna, "poscor", term2gene, universe),
       negcor = run_enricher(genes$RNA_symbol[genes$coef < 0], mirna, "negcor", term2gene, universe))
}

# For each miRNA-TF pair: the most extreme and the mean correlation of the overlapping genes
add_top_correlation_stats <- function(ora_df, cor_result, prefix) {
  gene_col <- paste0(prefix, "_overlap_geneID")
  pick_idx <- if (prefix == "poscor") which.max else which.min

  stats <- ora_df %>%
    dplyr::distinct(miRNA, TF, gene = .data[[gene_col]]) %>%
    tidyr::separate_rows(gene, sep = "/") %>%
    dplyr::inner_join(dplyr::select(cor_result, miRNA, RNA_symbol, coef),
                      by = c("miRNA", "gene" = "RNA_symbol"), relationship = "many-to-many") %>%
    dplyr::group_by(miRNA, TF) %>%
    dplyr::summarise(top_cor_gene = gene[pick_idx(coef)],
                     top_cor_coef = coef[pick_idx(coef)],
                     avg_cor_coef = mean(coef),
                     .groups = "drop")

  dplyr::left_join(ora_df, stats, by = c("miRNA", "TF"))
}

# ORA for every miRNA-TF pair in `overlap`.
# Returns list(poscor, negcor): the significant pairs for each direction.
run_ora_analysis <- function(overlap, cor_result, tf_target_df, direct_targets,
                             background_genes, fdr = 0.05, cores = 1) {
  # TF target sets limited to expressed genes, plus the background as its own term
  term2gene <- tf_target_df %>%
    dplyr::filter(TF_targets %in% background_genes) %>%
    dplyr::bind_rows(tibble::tibble(TF = "background", TF_targets = background_genes))

  overlap <- dplyr::filter(overlap, TF %in% tf_target_df$TF)
  mirnas <- unique(overlap$miRNA)
  cor_result <- dplyr::filter(cor_result, miRNA %in% mirnas)

  log_msg("Running ORA for ", length(mirnas), " miRNAs")
  run_one <- function(m) ora_for_mirna(m, cor_result, direct_targets, term2gene, background_genes)
  per_mirna <- if (cores > 1) parallel::mclapply(mirnas, run_one, mc.cores = cores) else lapply(mirnas, run_one)

  poscor <- overlap %>%
    dplyr::left_join(dplyr::bind_rows(lapply(per_mirna, `[[`, "poscor")), by = c("miRNA", "TF" = "ID")) %>%
    dplyr::filter(!is.na(poscor_ora_pval))
  negcor <- overlap %>%
    dplyr::left_join(dplyr::bind_rows(lapply(per_mirna, `[[`, "negcor")), by = c("miRNA", "TF" = "ID")) %>%
    dplyr::filter(!is.na(negcor_ora_pval))

  # Adjust p-values across both directions together
  p_adj <- p.adjust(c(poscor$poscor_ora_pval, negcor$negcor_ora_pval), method = "BH")
  poscor$poscor_ora_pval.adj <- p_adj[seq_len(nrow(poscor))]
  negcor$negcor_ora_pval.adj <- p_adj[nrow(poscor) + seq_len(nrow(negcor))]

  epsilon <- 1e-10  # avoid log(0)
  poscor <- poscor %>%
    dplyr::filter(poscor_ora_pval.adj < fdr) %>%
    add_top_correlation_stats(cor_result, "poscor") %>%
    dplyr::mutate(qval_plot = -log10(poscor_ora_pval.adj + epsilon),
                  miR_TF = make_mir_tf_label(miRNA, TF))
  negcor <- negcor %>%
    dplyr::filter(negcor_ora_pval.adj < fdr) %>%
    add_top_correlation_stats(cor_result, "negcor") %>%
    dplyr::mutate(qval_plot = -log10(negcor_ora_pval.adj + epsilon),
                  miR_TF = make_mir_tf_label(miRNA, TF))

  list(poscor = poscor, negcor = negcor)
}
