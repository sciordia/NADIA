# =============================================================================
# Benchmarking_Multiple.R
# Multiple-method OpDEA benchmarking: import, rank, and visualize
#
# Compares opdea_metrics from multiple normalization/imputation combinations.
# Ranking follows OpDEA methodology (Peng et al., Nature Comms 2024):
#   rank_final = mean(rank_nMCC, rank_G_mean, rank_pAUC_001, rank_pAUC_005, rank_pAUC_010)
#
# MIT License | Copyright (c) 2025 Sergio Ciordia
# =============================================================================

# --- Conditional operator ---
if (!exists("%||%", mode = "function")) `%||%` <- function(a, b) if (!is.null(a)) a else b

# =============================================================================
# SECTION 1: INTERNAL HELPERS
# =============================================================================

#' Required columns in the combined opdea data.frame
#' @keywords internal
.BM_REQUIRED_COLS <- c("Assay", "Comparison", "nMCC", "G_mean",
                       "pAUC_001", "pAUC_005", "pAUC_010")

#' Default metrics for ranking
#' @keywords internal
.BM_METRICS <- c("nMCC", "G_mean", "pAUC_001", "pAUC_005", "pAUC_010")


#' Read a TSV file with readr fallback
#' @param path File path
#' @return data.frame
#' @keywords internal
.bm_read_tsv <- function(path) {
  if (requireNamespace("readr", quietly = TRUE)) {
    as.data.frame(readr::read_tsv(path, show_col_types = FALSE),
                  stringsAsFactors = FALSE)
  } else {
    read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  }
}


#' Validate combined opdea data.frame
#' @param df data.frame to validate
#' @keywords internal
.bm_validate_opdea <- function(df) {
  if (!is.data.frame(df))
    stop("opdea_combined must be a data.frame.")

  missing <- setdiff(.BM_REQUIRED_COLS, colnames(df))
  if (length(missing) > 0)
    stop("Missing required columns in opdea_combined: ",
         paste(missing, collapse = ", "),
         "\nRequired: ", paste(.BM_REQUIRED_COLS, collapse = ", "))

  if (nrow(df) == 0)
    stop("opdea_combined has 0 rows.")

  invisible(TRUE)
}


#' Aggregate metrics by Assay
#'
#' @param df Combined opdea data.frame (with Assay column)
#' @param metrics Character vector of metric column names
#' @param agg_fun Aggregation function (mean or median)
#' @return data.frame with one row per Assay
#' @keywords internal
.bm_aggregate_metrics <- function(df, metrics, agg_fun = mean) {
  assays <- unique(df$Assay)
  agg_list <- lapply(assays, function(a) {
    sub <- df[df$Assay == a, , drop = FALSE]
    vals <- vapply(metrics, function(m) {
      agg_fun(sub[[m]], na.rm = TRUE)
    }, numeric(1))
    row <- data.frame(Assay = a, stringsAsFactors = FALSE)
    for (i in seq_along(metrics)) row[[metrics[i]]] <- round(vals[i], 4)
    row
  })
  do.call(rbind, agg_list)
}


#' Rank methods by descending metric value
#'
#' @param agg_df Aggregated data.frame (one row per Assay)
#' @param metrics Character vector of metric column names
#' @return data.frame with rank columns and rank_final, sorted by rank_final
#' @keywords internal
.bm_rank_methods <- function(agg_df, metrics) {
  rank_df <- data.frame(Assay = agg_df$Assay, stringsAsFactors = FALSE)

  for (m in metrics) {
    vals <- agg_df[[m]]
    vals[is.na(vals)] <- -Inf
    rank_df[[paste0("rank_", m)]] <- rank(-vals, ties.method = "average")
  }

  rank_cols <- paste0("rank_", metrics)
  rank_df$rank_final <- round(rowMeans(rank_df[, rank_cols, drop = FALSE]), 2)
  rank_df <- rank_df[order(rank_df$rank_final), , drop = FALSE]
  rownames(rank_df) <- NULL
  rank_df
}


#' Export data.frame to TSV
#' @param data data.frame
#' @param filepath Output path
#' @keywords internal
.bm_export_data <- function(data, filepath) {
  if (requireNamespace("readr", quietly = TRUE)) {
    readr::write_tsv(data, filepath)
  } else {
    write.table(data, filepath, sep = "\t", quote = FALSE, row.names = FALSE)
  }
  invisible(filepath)
}


#' Export ggplot2 to PNG
#' @param gg ggplot2 object
#' @param filepath Output path
#' @param width Width in inches
#' @param height Height in inches
#' @param dpi Resolution
#' @keywords internal
.bm_export_gg_plot <- function(gg, filepath, width = 12, height = 8, dpi = 150) {
  if (is.null(gg)) return(invisible(NULL))

  tryCatch({
    ggplot2::ggsave(filepath, plot = gg, width = width, height = height,
                    dpi = dpi, bg = "white")
  }, error = function(e) {
    warning("Error exporting plot to ", filepath, ": ", e$message)
  })

  invisible(filepath)
}


# =============================================================================
# SECTION 2: IMPORT FUNCTION
# =============================================================================

#' Import OpDEA metrics from multiple benchmark result folders
#'
#' Reads \code{benchmark_opdea_metrics.tsv} files from subdirectories,
#' adding an \code{Assay} column derived from the folder name.
#'
#' @param results_dir Parent directory containing benchmark subfolders
#' @param pattern Regex pattern for the opdea metrics filename
#'   (default: \code{"benchmark_opdea_metrics\\.tsv$"})
#' @param method_names Optional character vector of Assay names.
#'   If NULL, names are derived from the parent folder of each file.
#' @param recursive Logical. Search subdirectories recursively? (default: TRUE)
#'
#' @return data.frame in long format with columns:
#'   Assay, Comparison, nMCC, G_mean, pAUC_001, pAUC_005, pAUC_010
#'
#' @examples
#' \dontrun{
#' # Folder structure:
#' #   results/benchmark_cycloess_Impseq_min/benchmark_opdea_metrics.tsv
#' #   results/benchmark_quantile_knn_min/benchmark_opdea_metrics.tsv
#'
#' opdea_all <- import_opdea_results("results")
#' }
#' @export
import_opdea_results <- function(results_dir,
                                 pattern      = "benchmark_opdea_metrics\\.tsv$",
                                 method_names = NULL,
                                 recursive    = TRUE) {
  if (!dir.exists(results_dir))
    stop("Directory not found: ", results_dir)

  files <- list.files(results_dir, pattern = pattern,
                      full.names = TRUE, recursive = recursive)
  if (length(files) == 0)
    stop("No files matching pattern '", pattern, "' found in: ", results_dir)

  # --- Determine Assay names ---
  if (is.null(method_names)) {
    method_names <- basename(dirname(files))
  }
  if (length(method_names) != length(files))
    stop("Length of method_names (", length(method_names),
         ") must match number of files (", length(files), ").")

  # --- Read and combine ---
  df_list <- vector("list", length(files))
  for (i in seq_along(files)) {
    df <- .bm_read_tsv(files[i])
    df$Assay <- method_names[i]
    df_list[[i]] <- df
  }

  combined <- do.call(rbind, df_list)
  rownames(combined) <- NULL

  .bm_validate_opdea(combined)

  n_methods <- length(unique(combined$Assay))
  n_comps   <- length(unique(combined$Comparison))
  message("import_opdea_results: loaded ", n_methods, " method(s), ",
          n_comps, " comparison(s), ", nrow(combined), " rows.")

  combined
}


# =============================================================================
# SECTION 3: RANKING COMPUTATION
# =============================================================================

#' Compute OpDEA ranking across multiple methods
#'
#' For each metric, aggregates across comparisons (mean and median),
#' then ranks methods in descending order (rank 1 = best).
#' Final rank = average of individual metric ranks.
#'
#' @param opdea_combined data.frame with columns: Assay, Comparison,
#'   nMCC, G_mean, pAUC_001, pAUC_005, pAUC_010
#' @param metrics Character vector of metrics to include in ranking
#'   (default: all 5 OpDEA metrics)
#'
#' @return Named list with:
#'   \describe{
#'     \item{mean_aggregated}{data.frame of mean metric values per Assay}
#'     \item{median_aggregated}{data.frame of median metric values per Assay}
#'     \item{mean_ranking}{data.frame of ranks based on mean aggregation}
#'     \item{median_ranking}{data.frame of ranks based on median aggregation}
#'     \item{opdea_combined}{Input combined data.frame}
#'     \item{n_methods}{Number of methods}
#'     \item{n_comparisons}{Number of comparisons}
#'     \item{metrics_used}{Metrics included in ranking}
#'   }
#'
#' @examples
#' \dontrun{
#' ranking <- bm_compute_ranking(opdea_all)
#' ranking$mean_ranking
#' }
#' @export
bm_compute_ranking <- function(opdea_combined,
                               metrics = .BM_METRICS) {
  .bm_validate_opdea(opdea_combined)

  # Validate requested metrics exist
  avail <- intersect(metrics, colnames(opdea_combined))
  if (length(avail) == 0)
    stop("None of the requested metrics found in opdea_combined.")
  if (length(avail) < length(metrics)) {
    missing <- setdiff(metrics, avail)
    warning("Metrics not found (skipped): ", paste(missing, collapse = ", "))
    metrics <- avail
  }

  # Aggregate
  mean_agg   <- .bm_aggregate_metrics(opdea_combined, metrics, mean)
  median_agg <- .bm_aggregate_metrics(opdea_combined, metrics, median)

  # Rank
  mean_rank   <- .bm_rank_methods(mean_agg, metrics)
  median_rank <- .bm_rank_methods(median_agg, metrics)

  list(
    mean_aggregated   = mean_agg,
    median_aggregated = median_agg,
    mean_ranking      = mean_rank,
    median_ranking    = median_rank,
    opdea_combined    = opdea_combined,
    n_methods         = length(unique(opdea_combined$Assay)),
    n_comparisons     = length(unique(opdea_combined$Comparison)),
    metrics_used      = metrics
  )
}


# =============================================================================
# SECTION 4: VISUALIZATIONS
# =============================================================================

#' Heatmap of metric ranks by method
#'
#' Rows = Assay (ordered by rank_final), columns = metrics.
#' Fill color = rank value (green = best, red = worst).
#'
#' @param ranking_result List returned by \code{bm_compute_ranking()}
#' @param type "mean" or "median" aggregation (default: "mean")
#' @param title Optional plot title
#'
#' @return ggplot2 object
#' @export
bm_plot_ranking_heatmap <- function(ranking_result,
                                    type  = "mean",
                                    title = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required.")
  if (!requireNamespace("tidyr", quietly = TRUE))
    stop("Package 'tidyr' is required.")

  rank_df <- if (type == "median") {
    ranking_result$median_ranking
  } else {
    ranking_result$mean_ranking
  }

  metrics <- ranking_result$metrics_used
  rank_cols <- paste0("rank_", metrics)

  # Pivot to long
  plot_data <- rank_df[, c("Assay", rank_cols, "rank_final"), drop = FALSE]

  # Order Assay by rank_final
  plot_data$Assay <- factor(plot_data$Assay,
                            levels = rev(rank_df$Assay))

  long <- tidyr::pivot_longer(plot_data,
                              cols      = c(rank_cols, "rank_final"),
                              names_to  = "Metric",
                              values_to = "Rank")

  # Clean metric names for display
  long$Metric <- sub("^rank_", "", long$Metric)
  metric_order <- c(metrics, "final")
  long$Metric <- factor(long$Metric, levels = metric_order)

  n_methods <- ranking_result$n_methods

  title <- title %||% paste0("OpDEA Ranking Heatmap (", type, "-based)")
  subtitle <- paste0(ranking_result$n_methods, " methods, ",
                     ranking_result$n_comparisons, " comparisons")

  gg <- ggplot2::ggplot(long, ggplot2::aes(x = Metric, y = Assay, fill = Rank)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = round(Rank, 1)),
                       size = 3.5, color = "black") +
    ggplot2::scale_fill_gradient(low = "#2ca02c", high = "#d62728",
                                 limits = c(1, n_methods),
                                 name = "Rank") +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid        = ggplot2::element_blank(),
      plot.title        = ggplot2::element_text(face = "bold"),
      legend.position   = "right"
    )

  gg
}


#' Bar chart of final ranks
#'
#' Horizontal bars of rank_final by Assay, ordered best to worst.
#'
#' @param ranking_result List returned by \code{bm_compute_ranking()}
#' @param type "mean" or "median" aggregation (default: "mean")
#' @param title Optional plot title
#'
#' @return ggplot2 object
#' @export
bm_plot_ranking_bars <- function(ranking_result,
                                 type  = "mean",
                                 title = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required.")

  rank_df <- if (type == "median") {
    ranking_result$median_ranking
  } else {
    ranking_result$mean_ranking
  }

  n_methods <- ranking_result$n_methods

  rank_df$Assay <- factor(rank_df$Assay,
                          levels = rev(rank_df$Assay))

  title <- title %||% paste0("OpDEA Final Ranking (", type, "-based)")
  subtitle <- paste0(ranking_result$n_methods, " methods, ",
                     ranking_result$n_comparisons, " comparisons")

  gg <- ggplot2::ggplot(rank_df,
                        ggplot2::aes(x = Assay, y = rank_final, fill = rank_final)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = round(rank_final, 1)),
                       hjust = -0.2, size = 3.5) +
    ggplot2::scale_fill_gradient(low = "#2ca02c", high = "#d62728",
                                 limits = c(1, n_methods),
                                 name = "Rank") +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
    ggplot2::labs(title = title, subtitle = subtitle,
                  x = NULL, y = "Final Rank (lower = better)") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(face = "bold"),
      legend.position = "none"
    )

  gg
}


#' Heatmap of aggregated metric values (not ranks)
#'
#' Shows the actual metric values per method, with per-column color scaling.
#'
#' @param ranking_result List returned by \code{bm_compute_ranking()}
#' @param type "mean" or "median" aggregation (default: "mean")
#' @param title Optional plot title
#'
#' @return ggplot2 object
#' @export
bm_plot_metrics_heatmap <- function(ranking_result,
                                    type  = "mean",
                                    title = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required.")
  if (!requireNamespace("tidyr", quietly = TRUE))
    stop("Package 'tidyr' is required.")
  if (!requireNamespace("dplyr", quietly = TRUE))
    stop("Package 'dplyr' is required.")

  agg_df <- if (type == "median") {
    ranking_result$median_aggregated
  } else {
    ranking_result$mean_aggregated
  }

  # Order by mean ranking
  rank_df <- if (type == "median") {
    ranking_result$median_ranking
  } else {
    ranking_result$mean_ranking
  }

  metrics <- ranking_result$metrics_used

  # Pivot to long
  long <- tidyr::pivot_longer(agg_df,
                              cols      = dplyr::all_of(metrics),
                              names_to  = "Metric",
                              values_to = "Value")

  long$Assay  <- factor(long$Assay, levels = rev(rank_df$Assay))
  long$Metric <- factor(long$Metric, levels = metrics)

  # Normalize values per metric to [0, 1] for consistent color scale
  long <- dplyr::group_by(long, Metric)
  long <- dplyr::mutate(long,
    Value_scaled = if (max(Value, na.rm = TRUE) == min(Value, na.rm = TRUE)) {
      0.5
    } else {
      (Value - min(Value, na.rm = TRUE)) /
        (max(Value, na.rm = TRUE) - min(Value, na.rm = TRUE))
    }
  )
  long <- dplyr::ungroup(long)

  title <- title %||% paste0("OpDEA Metric Values (", type, "-based)")
  subtitle <- paste0(ranking_result$n_methods, " methods, ",
                     ranking_result$n_comparisons, " comparisons | ",
                     "color scaled per metric")

  gg <- ggplot2::ggplot(long,
                        ggplot2::aes(x = Metric, y = Assay, fill = Value_scaled)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.4f", Value)),
                       size = 3, color = "black") +
    ggplot2::scale_fill_gradient(low = "#fee0d2", high = "#2ca02c",
                                 name = "Relative\nPerformance",
                                 labels = c("Worst", "Best"),
                                 breaks = c(0, 1)) +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid        = ggplot2::element_blank(),
      plot.title        = ggplot2::element_text(face = "bold"),
      legend.position   = "right"
    )

  gg
}


#' Boxplot comparison of metric distributions across methods
#'
#' Faceted boxplot showing the distribution of each metric across
#' comparisons for every method.
#'
#' @param opdea_combined data.frame with Assay, Comparison, and metric columns
#' @param metrics Character vector of metrics to plot (default: all 5)
#' @param title Optional plot title
#'
#' @return ggplot2 object
#' @export
bm_plot_metrics_comparison <- function(opdea_combined,
                                       metrics = .BM_METRICS,
                                       title   = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required.")
  if (!requireNamespace("tidyr", quietly = TRUE))
    stop("Package 'tidyr' is required.")
  if (!requireNamespace("dplyr", quietly = TRUE))
    stop("Package 'dplyr' is required.")

  avail <- intersect(metrics, colnames(opdea_combined))
  if (length(avail) == 0) stop("No valid metrics found in opdea_combined.")

  long <- tidyr::pivot_longer(opdea_combined,
                              cols      = dplyr::all_of(avail),
                              names_to  = "Metric",
                              values_to = "Value")

  long$Metric <- factor(long$Metric, levels = avail)

  n_methods <- length(unique(long$Assay))
  n_comps   <- length(unique(long$Comparison))

  title <- title %||% "OpDEA Metrics Distribution by Method"
  subtitle <- paste0(n_methods, " methods, ", n_comps, " comparisons")

  gg <- ggplot2::ggplot(long,
                        ggplot2::aes(x = Assay, y = Value, fill = Assay)) +
    ggplot2::geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    ggplot2::geom_jitter(width = 0.15, size = 1.2, alpha = 0.6) +
    ggplot2::facet_wrap(~ Metric, scales = "free_y", ncol = 3) +
    ggplot2::labs(title = title, subtitle = subtitle,
                  x = NULL, y = "Value") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_text(angle = 60, hjust = 1, size = 8),
      plot.title        = ggplot2::element_text(face = "bold"),
      legend.position   = "none",
      strip.text        = ggplot2::element_text(face = "bold")
    )

  gg
}


# =============================================================================
# SECTION 5: ORCHESTRATOR
# =============================================================================

#' Multiple-method OpDEA Benchmarking
#'
#' Orchestrator that imports, ranks, visualizes, and exports results from
#' multiple normalization/imputation method combinations.
#'
#' @param opdea_combined Optional pre-built data.frame with columns:
#'   Assay, Comparison, nMCC, G_mean, pAUC_001, pAUC_005, pAUC_010.
#'   If NULL, results are imported from \code{results_dir}.
#' @param results_dir Parent directory containing benchmark subfolders
#'   (each with \code{benchmark_opdea_metrics.tsv}). Used only if
#'   \code{opdea_combined} is NULL.
#' @param pattern Regex for opdea metrics filename (default:
#'   \code{"benchmark_opdea_metrics\\.tsv$"})
#' @param method_names Optional Assay names for file import
#' @param recursive Search subdirectories? (default: TRUE)
#' @param metrics Character vector of metrics for ranking
#'   (default: all 5 OpDEA metrics)
#' @param plots Which plots to generate: "all" or character vector of names.
#'   Valid names: "ranking_heatmap_mean", "ranking_heatmap_median",
#'   "ranking_bars_mean", "ranking_bars_median",
#'   "metrics_heatmap_mean", "metrics_heatmap_median",
#'   "metrics_comparison"
#' @param verbose Print progress messages (default: TRUE)
#' @param output_dir Directory for exporting results (NULL = no export)
#' @param export_plots Export plots as PNG (default: TRUE)
#' @param export_tables Export tables as TSV (default: TRUE)
#' @param plot_width Plot width in inches (default: 12)
#' @param plot_height Plot height in inches (default: 8)
#' @param plot_dpi Plot resolution (default: 150)
#'
#' @return Named list with:
#'   \describe{
#'     \item{opdea_combined}{Combined input data.frame}
#'     \item{mean_aggregated}{Mean of metrics per Assay}
#'     \item{median_aggregated}{Median of metrics per Assay}
#'     \item{mean_ranking}{Ranking table (mean-based)}
#'     \item{median_ranking}{Ranking table (median-based)}
#'     \item{gg_ranking_heatmap_mean}{Heatmap of ranks (mean)}
#'     \item{gg_ranking_heatmap_median}{Heatmap of ranks (median)}
#'     \item{gg_ranking_bars_mean}{Bar chart of final rank (mean)}
#'     \item{gg_ranking_bars_median}{Bar chart of final rank (median)}
#'     \item{gg_metrics_heatmap_mean}{Heatmap of metric values (mean)}
#'     \item{gg_metrics_heatmap_median}{Heatmap of metric values (median)}
#'     \item{gg_metrics_comparison}{Boxplot distributions}
#'     \item{parameters}{List of parameters used}
#'   }
#'
#' @examples
#' \dontrun{
#' # --- Option A: From files ---
#' result <- benchmarking_multiple(
#'   results_dir = "results",
#'   output_dir  = "results/bm_multiple"
#' )
#'
#' # --- Option B: In-memory ---
#' opdea_all <- do.call(rbind, list(
#'   cbind(bench1$opdea_metrics, Assay = "cycloess_Impseq_min"),
#'   cbind(bench2$opdea_metrics, Assay = "quantile_knn_min")
#' ))
#' result <- benchmarking_multiple(
#'   opdea_combined = opdea_all,
#'   output_dir     = "results/bm_multiple"
#' )
#'
#' result$mean_ranking
#' result$gg_ranking_heatmap_mean
#' }
#' @export
benchmarking_multiple <- function(opdea_combined = NULL,
                                  results_dir    = NULL,
                                  pattern        = "benchmark_opdea_metrics\\.tsv$",
                                  method_names   = NULL,
                                  recursive      = TRUE,
                                  metrics        = .BM_METRICS,
                                  plots          = "all",
                                  verbose        = TRUE,
                                  output_dir     = NULL,
                                  export_plots   = TRUE,
                                  export_tables  = TRUE,
                                  plot_width     = 12,
                                  plot_height    = 8,
                                  plot_dpi       = 150) {

  # --- Required packages ---
  for (pkg in c("ggplot2", "dplyr", "tidyr")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop("Package '", pkg, "' is required for benchmarking_multiple().")
  }

  # ==========================================================================
  # STEP 1: Input resolution
  # ==========================================================================
  if (!is.null(opdea_combined) && !is.null(results_dir)) {
    warning("Both opdea_combined and results_dir provided. Using opdea_combined.")
  }

  if (is.null(opdea_combined)) {
    if (is.null(results_dir))
      stop("Provide either opdea_combined or results_dir.")
    if (verbose) message("=== IMPORTING OPDEA RESULTS ===")
    opdea_combined <- import_opdea_results(results_dir, pattern,
                                           method_names, recursive)
  }

  .bm_validate_opdea(opdea_combined)

  n_methods <- length(unique(opdea_combined$Assay))
  n_comps   <- length(unique(opdea_combined$Comparison))

  if (verbose) {
    message("=== MULTIPLE-METHOD OPDEA BENCHMARKING ===")
    message("  Methods: ", n_methods,
            " | Comparisons: ", n_comps,
            " | Rows: ", nrow(opdea_combined))
  }

  # ==========================================================================
  # STEP 2: Compute ranking
  # ==========================================================================
  if (verbose) message("  Computing ranking ...")
  ranking <- bm_compute_ranking(opdea_combined, metrics)

  if (verbose) {
    message("  Best method (mean):   ", ranking$mean_ranking$Assay[1],
            " (rank_final = ", ranking$mean_ranking$rank_final[1], ")")
    message("  Best method (median): ", ranking$median_ranking$Assay[1],
            " (rank_final = ", ranking$median_ranking$rank_final[1], ")")
  }

  # ==========================================================================
  # STEP 3: Plot registry
  # ==========================================================================
  all_plot_names <- c("ranking_heatmap_mean", "ranking_heatmap_median",
                      "ranking_bars_mean", "ranking_bars_median",
                      "metrics_heatmap_mean", "metrics_heatmap_median",
                      "metrics_comparison")

  plot_fns <- list(
    ranking_heatmap_mean   = function() bm_plot_ranking_heatmap(ranking, "mean"),
    ranking_heatmap_median = function() bm_plot_ranking_heatmap(ranking, "median"),
    ranking_bars_mean      = function() bm_plot_ranking_bars(ranking, "mean"),
    ranking_bars_median    = function() bm_plot_ranking_bars(ranking, "median"),
    metrics_heatmap_mean   = function() bm_plot_metrics_heatmap(ranking, "mean"),
    metrics_heatmap_median = function() bm_plot_metrics_heatmap(ranking, "median"),
    metrics_comparison     = function() bm_plot_metrics_comparison(opdea_combined, metrics)
  )

  # Determine which plots to run
  if (identical(plots, "all")) {
    selected <- all_plot_names
  } else {
    unknown <- setdiff(plots, all_plot_names)
    if (length(unknown) > 0)
      warning("Unknown plot name(s) ignored: ", paste(unknown, collapse = ", "))
    selected <- intersect(plots, all_plot_names)
    if (length(selected) == 0)
      stop("No valid plot names provided.")
  }

  # Execute plots
  plot_result <- vector("list", length(selected))
  names(plot_result) <- selected

  for (nm in selected) {
    if (verbose) message("  Generating plot: ", nm, " ...")
    plot_result[[nm]] <- tryCatch(
      plot_fns[[nm]](),
      error = function(e) {
        warning("bm_plot_", nm, "() failed: ", conditionMessage(e))
        NULL
      }
    )
  }

  n_ok   <- sum(!vapply(plot_result, is.null, logical(1)))
  n_fail <- length(selected) - n_ok
  if (verbose) {
    message("  ", n_ok, " plot(s) generated",
            if (n_fail > 0) paste0(", ", n_fail, " failed") else ".")
  }

  # ==========================================================================
  # STEP 4: Export
  # ==========================================================================
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    if (verbose) message("  Exporting to: ", output_dir)

    if (export_tables) {
      .bm_export_data(opdea_combined,
                      file.path(output_dir, "bm_multiple_opdea_combined.tsv"))
      .bm_export_data(ranking$mean_aggregated,
                      file.path(output_dir, "bm_multiple_mean_aggregated.tsv"))
      .bm_export_data(ranking$median_aggregated,
                      file.path(output_dir, "bm_multiple_median_aggregated.tsv"))
      .bm_export_data(ranking$mean_ranking,
                      file.path(output_dir, "bm_multiple_mean_ranking.tsv"))
      .bm_export_data(ranking$median_ranking,
                      file.path(output_dir, "bm_multiple_median_ranking.tsv"))
      if (verbose) message("  Exported 5 TSV files.")
    }

    if (export_plots) {
      for (nm in selected) {
        if (!is.null(plot_result[[nm]])) {
          .bm_export_gg_plot(plot_result[[nm]],
                             file.path(output_dir, paste0("bm_multiple_", nm, ".png")),
                             plot_width, plot_height, plot_dpi)
        }
      }
      if (verbose) message("  Exported ", n_ok, " PNG files.")
    }
  }

  # ==========================================================================
  # STEP 5: Build result
  # ==========================================================================
  result <- list(
    opdea_combined    = opdea_combined,
    mean_aggregated   = ranking$mean_aggregated,
    median_aggregated = ranking$median_aggregated,
    mean_ranking      = ranking$mean_ranking,
    median_ranking    = ranking$median_ranking
  )

  # Add plots with gg_ prefix
  for (nm in selected) {
    result[[paste0("gg_", nm)]] <- plot_result[[nm]]
  }

  result$parameters <- list(
    metrics       = metrics,
    n_methods     = n_methods,
    n_comparisons = n_comps,
    plots         = selected
  )

  if (verbose) message("=== DONE ===")

  result
}


# =============================================================================
# SECTION 6: EXAMPLE WORKFLOW
# =============================================================================

# --- Option A: In-memory loop ---
#
# source("R/Processing.R")
# source("R/Benchmarking_Single.R")
# source("R/Benchmarking_Multiple.R")
#
# library(dplyr)
#
# # Define method combinations to test
# combos <- list(
#   list(norm = "cycloess", imp = "combo", mar = "Impseqrob", mnar = "min"),
#   list(norm = "quantile", imp = "combo", mar = "knn",       mnar = "min"),
#   list(norm = "log2Norm", imp = "combo", mar = "Impseqrob", mnar = "MinDet")
# )
#
# # Run processing + benchmarking for each combo
# opdea_list <- list()
# for (combo in combos) {
#   result <- process_proteomics(
#     file_path       = "data/input.tsv",
#     metadata_path   = "data/metadata.tsv",
#     norm_method     = combo$norm,
#     imp_method      = combo$imp,
#     mar_method      = combo$mar,
#     mnar_method     = combo$mnar
#   )
#
#   bench <- benchmarking_proteomics(
#     de_res          = result$DEPs_results,
#     species_df      = species_df,
#     expected_values = expected,
#     alpha           = 0.05,
#     output_dir      = paste0("results/benchmark_", combo$norm, "_",
#                              combo$mar, "_", combo$mnar)
#   )
#
#   opdea <- bench$opdea_metrics
#   opdea$Assay <- paste(combo$norm, combo$mar, combo$mnar, sep = "_")
#   opdea_list[[length(opdea_list) + 1]] <- opdea
# }
#
# opdea_all <- do.call(rbind, opdea_list)
#
# bm_result <- benchmarking_multiple(
#   opdea_combined = opdea_all,
#   output_dir     = "results/bm_multiple",
#   verbose        = TRUE
# )
#
# bm_result$mean_ranking
# bm_result$gg_ranking_heatmap_mean
#
#
# --- Option B: Import from files ---
#
# source("R/Benchmarking_Multiple.R")
#
# # Assumes folders like:
# #   results/benchmark_cycloess_Impseqrob_min/benchmark_opdea_metrics.tsv
# #   results/benchmark_quantile_knn_min/benchmark_opdea_metrics.tsv
# #   results/benchmark_log2Norm_Impseqrob_MinDet/benchmark_opdea_metrics.tsv
#
# bm_result <- benchmarking_multiple(
#   results_dir = "results",
#   output_dir  = "results/bm_multiple",
#   verbose     = TRUE
# )
#
# bm_result$mean_ranking
# bm_result$gg_ranking_bars_mean
