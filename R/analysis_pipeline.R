# =============================================================================
# analysis_pipeline.R
# MirTfPipeline: an R6 object that runs the whole miRNA-TF analysis for one
# cancer type.
#
# Every intermediate result is a named "step". Asking for a step with
# pipeline$get("<step>") will, in order:
#   1. return it from memory if it was already loaded this session;
#   2. recompute and save it if it is listed in config$regenerate
#      (or is downstream of a listed step when cascade_regeneration = TRUE);
#   3. read it from disk if its output files exist;
#   4. otherwise compute it (pulling in whatever it depends on) and save it.
# So only the work that is actually needed is ever done.
# =============================================================================

#' Step definitions. Each step is a list with:
#'   deps          names of steps it depends on (used for cascading)
#'   files         config$files keys of its output (NULL = memory only)
#'   dir           "study" (default) or "common" for where files live
#'   input_only    TRUE if the file is a raw input that cannot be computed
#'   explicit_only TRUE if compute only runs when named in `regenerate`
#'   compute       function(p) -> value  (p is the pipeline)
#'   read          function(paths) -> value
#'   write         function(value, paths)
pipeline_steps <- function() {
  cfg_of <- function(p) p$config
  verbose <- function(p) p$config$verbose

  list(
    tcga_data = list(
      files = c("rna_rds", "mirna_rds"), explicit_only = TRUE,
      compute = function(p) {
        paths <- p$output_paths("tcga_data")
        download_tcga_data(cfg_of(p)$study, paths[["rna_rds"]], paths[["mirna_rds"]])
      },
      read = function(paths) list(rna = readRDS(paths[["rna_rds"]]),
                                  mirna = readRDS(paths[["mirna_rds"]])),
      write = NULL  # download_tcga_data() saves the RDS files itself
    ),

    gene_mapping = list(
      deps = "tcga_data", files = "gene_mapping",
      compute = function(p) extract_gene_mapping(p$get("tcga_data")$rna, cfg_of(p)$filters$gene_type),
      read = function(paths) read_table_csv(paths[[1]]),
      write = function(x, paths) write_table_csv(x, paths[[1]])
    ),

    expression = list(
      deps = "tcga_data", files = c("rna_filtered", "mirna_filtered"),
      compute = function(p) {
        raw <- p$get("tcga_data")
        prepare_expression_data(raw$rna, raw$mirna, cfg_of(p)$filters, verbose(p))
      },
      read = function(paths) list(rna = read_matrix_csv(paths[["rna_filtered"]]),
                                  mirna = read_matrix_csv(paths[["mirna_filtered"]])),
      write = function(x, paths) {
        write_matrix_csv(x$rna, paths[["rna_filtered"]])
        write_matrix_csv(x$mirna, paths[["mirna_filtered"]])
      }
    ),

    correlation = list(
      deps = c("expression", "gene_mapping"), files = "correlation",
      compute = function(p) {
        expr <- p$get("expression")
        a <- cfg_of(p)$analysis
        compute_mrna_mirna_correlation(expr$rna, expr$mirna, p$get("gene_mapping"),
                                       method = a$correlation_method,
                                       fdr = a$correlation_fdr,
                                       p_adjust_method = a$p_adjust_method)
      },
      read = function(paths) read_table_csv(paths[[1]]),
      write = function(x, paths) write_table_csv(x, paths[[1]])
    ),

    known_tfs = list(
      files = "human_tfs", dir = "common", input_only = TRUE,
      read = function(paths) load_known_tfs(paths[[1]])
    ),

    mirna_tf_interactions = list(
      deps = c("expression", "known_tfs"), files = "mirna_tf_interactions",
      compute = function(p) {
        query_mirna_tf_interactions(rownames(p$get("expression")$mirna),
                                    p$get("known_tfs"),
                                    cfg_of(p)$analysis$multimir_table)
      },
      read = function(paths) standardise_interactions(read_table_csv(paths[[1]])),
      write = function(x, paths) write_table_csv(x, paths[[1]])
    ),

    # Named vector TF -> ChIP-Atlas file. Downloads only if "chip_atlas" is
    # explicitly in `regenerate`; otherwise uses the files already on disk.
    chip_atlas = list(
      deps = "mirna_tf_interactions",
      compute = function(p) {
        cfg <- cfg_of(p)
        if ("chip_atlas" %in% cfg$regenerate) {
          tfs <- unique(p$get("mirna_tf_interactions")$target_symbol)
          log_msg("Downloading ChIP-Atlas targets for ", length(tfs), " TFs", verbose = cfg$verbose)
          download_tf_target_files(tfs, cfg$paths$tf_targets_dir, cfg$chip_atlas,
                                   overwrite = TRUE, verbose = cfg$verbose)
        }
        list_tf_target_files(cfg$paths$tf_targets_dir, cfg$chip_atlas)
      }
    ),

    cell_lines = list(
      files = "cell_lines", input_only = TRUE,
      read = function(paths) load_cell_lines(paths[[1]])
    ),

    tf_targets = list(
      deps = c("chip_atlas", "cell_lines"), files = c("tf_target_df", "tf_target_count"),
      compute = function(p) {
        build_tf_targets(p$get("chip_atlas"), p$get("cell_lines"),
                         cfg_of(p)$analysis$gene_peaks_in_cell_lines_threshold, verbose(p))
      },
      read = function(paths) list(target_df = read_table_csv(paths[["tf_target_df"]]),
                                  target_count = read_table_csv(paths[["tf_target_count"]])),
      write = function(x, paths) {
        write_table_csv(x$target_df, paths[["tf_target_df"]])
        write_table_csv(x$target_count, paths[["tf_target_count"]])
      }
    ),

    overlap = list(
      deps = c("correlation", "mirna_tf_interactions", "chip_atlas"),
      compute = function(p) {
        build_mirna_tf_overlap(p$get("correlation"), p$get("mirna_tf_interactions"),
                               names(p$get("chip_atlas")))
      }
    ),

    mirna_direct_targets = list(
      deps = "overlap", files = "mirna_direct_targets",
      compute = function(p) {
        query_mirna_targets(unique(p$get("overlap")$miRNA), cfg_of(p)$analysis$multimir_table)
      },
      read = function(paths) standardise_interactions(read_table_csv(paths[[1]])),
      write = function(x, paths) write_table_csv(x, paths[[1]])
    ),

    ora = list(
      deps = c("overlap", "correlation", "tf_targets", "mirna_direct_targets",
               "expression", "gene_mapping"),
      files = c("ora_poscor", "ora_negcor"),
      compute = function(p) {
        a <- cfg_of(p)$analysis
        gm <- p$get("gene_mapping")
        background <- gm$gene_symbol[gm$gene_id %in% rownames(p$get("expression")$rna)]
        run_ora_analysis(overlap          = p$get("overlap"),
                         cor_result       = p$get("correlation"),
                         tf_target_df     = p$get("tf_targets")$target_df,
                         direct_targets   = p$get("mirna_direct_targets"),
                         background_genes = background,
                         fdr              = a$ora_fdr,
                         p_adjust_method  = a$p_adjust_method,
                         cores            = a$ora_cores,
                         verbose          = verbose(p))
      },
      read = function(paths) list(poscor = read_table_csv(paths[["ora_poscor"]]),
                                  negcor = read_table_csv(paths[["ora_negcor"]])),
      write = function(x, paths) {
        write_table_csv(x$poscor, paths[["ora_poscor"]])
        write_table_csv(x$negcor, paths[["ora_negcor"]])
      }
    ),

    gene_zscores = list(
      deps = c("expression", "correlation"),
      compute = function(p) gene_expression_zscores(p$get("expression")$rna, p$get("correlation"))
    ),

    mirna_zscores = list(
      deps = "expression",
      compute = function(p) expression_zscores(p$get("expression")$mirna, "miRNA")
    )
  )
}

MirTfPipeline <- R6::R6Class("MirTfPipeline",
  public = list(
    config = NULL,

    initialize = function(config) {
      self$config <- config
      private$steps <- pipeline_steps()
      private$cache <- new.env(parent = emptyenv())
      private$stale <- private$resolve_stale(config$regenerate, config$cascade_regeneration)
      ensure_dir(config$paths$results_dir)
    },

    #' Get a step's result (from memory, disk, or by computing it).
    get = function(name) {
      if (!name %in% names(private$steps)) {
        stop("Unknown step '", name, "'. Steps: ", paste(names(private$steps), collapse = ", "),
             call. = FALSE)
      }
      if (!exists(name, envir = private$cache, inherits = FALSE)) {
        assign(name, private$load_or_compute(name), envir = private$cache)
      }
      base::get(name, envir = private$cache, inherits = FALSE)
    },

    #' Paths a step reads from (results_dir first, then study_dir).
    input_paths = function(name) {
      step <- private$steps[[name]]
      vapply(step$files, function(key) {
        file <- private$file_name(key)
        if (identical(step$dir, "common")) return(file.path(self$config$paths$common_dir, file))
        out <- file.path(self$config$paths$results_dir, file)
        if (file.exists(out)) out else file.path(self$config$paths$study_dir, file)
      }, character(1)) %>% stats::setNames(step$files)
    },

    #' Paths a step writes to.
    output_paths = function(name) {
      step <- private$steps[[name]]
      vapply(step$files, function(key) file.path(self$config$paths$results_dir, private$file_name(key)),
             character(1)) %>% stats::setNames(step$files)
    },

    #' Which steps will be read from disk vs regenerated.
    status = function() {
      file_steps <- Filter(function(s) !is.null(s$files), private$steps)
      tibble::tibble(
        step   = names(file_steps),
        action = vapply(names(file_steps), function(n) {
          if (n %in% private$stale && !isTRUE(file_steps[[n]]$input_only)) return("regenerate")
          if (all(file.exists(self$input_paths(n)))) "read from disk" else "compute (missing)"
        }, character(1)),
        files  = vapply(names(file_steps), function(n) paste(basename(self$input_paths(n)), collapse = ", "),
                        character(1))
      )
    },

    #' All results as a named list (loads / computes whatever is needed).
    results = function() {
      expr <- self$get("expression")
      tf_targets <- self$get("tf_targets")
      ora <- self$get("ora")
      list(
        config               = self$config,
        rna_expression       = expr$rna,
        mirna_expression     = expr$mirna,
        cor_result           = self$get("correlation"),
        known_tfs            = self$get("known_tfs"),
        interactions_tfs     = self$get("mirna_tf_interactions"),
        tf_available         = names(self$get("chip_atlas")),
        cell_lines           = self$get("cell_lines"),
        tf_target_df         = tf_targets$target_df,
        tf_target_count      = tf_targets$target_count,
        overlap              = self$get("overlap"),
        mirna_direct_targets = self$get("mirna_direct_targets"),
        ora_poscor           = ora$poscor,
        ora_negcor           = ora$negcor,
        gene_zscores         = self$get("gene_zscores"),
        mirna_zscores        = self$get("mirna_zscores")
      )
    },

    #' Key numbers reported in the Rmd, as a list.
    summary = function(res = self$results()) {
      bias <- test_replicate_bias(res$tf_target_count)
      list(
        samples              = ncol(res$rna_expression),
        genes                = nrow(res$rna_expression),
        mirnas               = nrow(res$mirna_expression),
        correlation          = summarise_correlation(res$cor_result),
        known_tfs_correlated = length(intersect(res$known_tfs$HGNC.symbol, res$cor_result$RNA_symbol)),
        mirna_tf_interactions = list(
          n_interactions = nrow(res$interactions_tfs),
          n_mirnas       = dplyr::n_distinct(res$interactions_tfs$miRNA),
          n_tfs          = dplyr::n_distinct(res$interactions_tfs$target_symbol),
          n_tfs_with_chip_atlas = length(res$tf_available)
        ),
        tf_coverage          = summarise_tf_coverage(res$overlap, res$cell_lines, self$config$study),
        replicate_bias       = list(rho = unname(bias$estimate), p_value = bias$p.value),
        ora = list(
          n_sig_pairs_pos  = nrow(res$ora_poscor),
          n_sig_pairs_neg  = nrow(res$ora_negcor),
          n_sig_mirnas_pos = dplyr::n_distinct(res$ora_poscor$miRNA),
          n_sig_mirnas_neg = dplyr::n_distinct(res$ora_negcor$miRNA)
        )
      )
    },

    #' Run everything: load/compute all steps, summarise, plot.
    run = function(make_plots = TRUE) {
      cfg <- self$config
      log_msg("==== ", cfg$study, " ====", verbose = cfg$verbose)
      if (cfg$verbose) print(self$status(), n = Inf)

      res <- self$results()
      res$summary <- self$summary(res)
      if (cfg$verbose) print_pipeline_summary(res$summary)

      if (make_plots) res$plots <- generate_all_plots(res, cfg)

      if (isTRUE(cfg$export_to_global)) {
        list2env(res, envir = globalenv())
        log_msg("Results copied to the global environment", verbose = cfg$verbose)
      }
      log_msg("Done: ", cfg$study, verbose = cfg$verbose)
      invisible(res)
    },

    print = function(...) {
      cat("<MirTfPipeline>", self$config$study, "\n")
      cat("  regenerate:", if (length(self$config$regenerate)) paste(self$config$regenerate, collapse = ", ") else "nothing", "\n")
      cat("  loaded:    ", paste(ls(private$cache), collapse = ", "), "\n")
      invisible(self)
    }
  ),

  private = list(
    steps = NULL,
    cache = NULL,
    stale = character(0),

    file_name = function(key) {
      f <- self$config$files[[key]]
      if (is.null(f)) stop("config$files$", key, " is not set.", call. = FALSE)
      f
    },

    # Requested steps plus, optionally, everything downstream of them.
    resolve_stale = function(requested, cascade) {
      stale <- intersect(requested, names(private$steps))
      if (!isTRUE(cascade)) return(stale)
      repeat {
        downstream <- names(Filter(function(s) any(s$deps %in% stale), private$steps))
        new <- setdiff(downstream, stale)
        if (!length(new)) return(stale)
        stale <- c(stale, new)
      }
    },

    load_or_compute = function(name) {
      step <- private$steps[[name]]
      verbose <- self$config$verbose
      if (is.null(step$files)) return(step$compute(self))

      paths <- self$input_paths(name)
      requested <- name %in% private$stale && !isTRUE(step$input_only)

      if (!requested && all(file.exists(paths))) {
        log_msg("Reading   ", name, " <- ", paste(basename(paths), collapse = ", "), verbose = verbose)
        return(step$read(paths))
      }

      if (isTRUE(step$input_only)) {
        stop("Input file(s) for '", name, "' not found: ", paste(paths, collapse = ", "), call. = FALSE)
      }
      if (isTRUE(step$explicit_only) && !name %in% self$config$regenerate) {
        stop("'", name, "' files not found: ", paste(paths, collapse = ", "),
             "\nAdd \"", name, "\" to REGENERATE to download them.", call. = FALSE)
      }

      log_msg(if (requested) "Regenerating " else "Computing (no saved output) ", name,
              verbose = verbose)
      value <- step$compute(self)
      if (!is.null(step$write)) {
        out <- self$output_paths(name)
        step$write(value, out)
        log_msg("Saved     ", name, " -> ", paste(basename(out), collapse = ", "), verbose = verbose)
      }
      value
    }
  )
)

print_pipeline_summary <- function(s) {
  cat("\n---------------- Summary ----------------\n")
  cat(sprintf("Samples: %d | genes: %d | miRNAs: %d\n", s$samples, s$genes, s$mirnas))
  cat(sprintf("Significant mRNA-miRNA correlations: %d (%d positive, %d negative)\n",
              s$correlation$n_pairs, s$correlation$n_positive, s$correlation$n_negative))
  cat(sprintf("  involving %d miRNAs and %d genes\n", s$correlation$n_mirnas, s$correlation$n_genes))
  cat(sprintf("Known TFs among correlated genes: %d\n", s$known_tfs_correlated))
  cat(sprintf("Validated miRNA-TF interactions: %d (%d miRNAs, %d TFs; ChIP-Atlas data for %d TFs)\n",
              s$mirna_tf_interactions$n_interactions, s$mirna_tf_interactions$n_mirnas,
              s$mirna_tf_interactions$n_tfs, s$mirna_tf_interactions$n_tfs_with_chip_atlas))
  cat(sprintf("TFs with ChIP data & correlated miRNA: %d | with cancer cell lines: %d | both: %d\n",
              s$tf_coverage$n_tfs_with_chip_and_mirna, s$tf_coverage$n_tfs_with_cancer_cell_lines,
              s$tf_coverage$n_tfs_in_both))
  cat(sprintf("Target count vs replicate count (Spearman): rho = %.3f, p = %.3g\n",
              s$replicate_bias$rho, s$replicate_bias$p_value))
  cat(sprintf("Significant miRNA-TF pairs (ORA): %d positive (%d miRNAs), %d negative (%d miRNAs)\n",
              s$ora$n_sig_pairs_pos, s$ora$n_sig_mirnas_pos,
              s$ora$n_sig_pairs_neg, s$ora$n_sig_mirnas_neg))
  cat("-----------------------------------------\n\n")
}
