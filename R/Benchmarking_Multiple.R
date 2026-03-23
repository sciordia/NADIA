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

#' Required columns in the combined confusion data.frame
#' @keywords internal
.BM_CONFUSION_COLS <- c("Assay", "Comparison", "TP", "FP", "TN", "FN")

#' Minimum required columns in the combined classified data.frame
#' @keywords internal
.BM_CLASSIFIED_COLS <- c("Assay", "Comparison", "truth")

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


#' Validate combined confusion data.frame
#' @param df data.frame to validate
#' @keywords internal
.bm_validate_confusion <- function(df) {
  if (!is.data.frame(df))
    stop("confusion_combined must be a data.frame.")

  missing <- setdiff(.BM_CONFUSION_COLS, colnames(df))
  if (length(missing) > 0)
    stop("Missing required columns in confusion_combined: ",
         paste(missing, collapse = ", "),
         "\nRequired: ", paste(.BM_CONFUSION_COLS, collapse = ", "))

  if (nrow(df) == 0)
    stop("confusion_combined has 0 rows.")

  invisible(TRUE)
}


#' Validate combined classified data.frame
#' @param df data.frame to validate
#' @param p_col Name of the p-value column
#' @keywords internal
.bm_validate_classified <- function(df, p_col = "adj.P.Val") {
  if (!is.data.frame(df))
    stop("classified_combined must be a data.frame.")

  missing <- setdiff(.BM_CLASSIFIED_COLS, colnames(df))
  if (length(missing) > 0)
    stop("Missing required columns in classified_combined: ",
         paste(missing, collapse = ", "),
         "\nRequired: ", paste(.BM_CLASSIFIED_COLS, collapse = ", "))

  if (!p_col %in% colnames(df))
    stop("P-value column '", p_col, "' not found in classified_combined.")

  if (nrow(df) == 0)
    stop("classified_combined has 0 rows.")

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



#' Import confusion matrices from multiple benchmark result folders
#'
#' Reads \code{benchmark_confusion_overall.tsv} files from subdirectories,
#' adding an \code{Assay} column derived from the folder name.
#'
#' @param results_dir Parent directory containing benchmark subfolders
#' @param pattern Regex pattern for the confusion filename
#' @param method_names Optional character vector of Assay names
#' @param recursive Search subdirectories recursively? (default: TRUE)
#'
#' @return data.frame with columns:
#'   Assay, Comparison, N, TP, FP, TN, FN, plus percentage columns
#'
#' @export
import_confusion_results <- function(results_dir,
                                     pattern      = "benchmark_confusion_overall\\.tsv$",
                                     method_names = NULL,
                                     recursive    = TRUE) {
  if (!dir.exists(results_dir))
    stop("Directory not found: ", results_dir)

  files <- list.files(results_dir, pattern = pattern,
                      full.names = TRUE, recursive = recursive)
  if (length(files) == 0)
    stop("No files matching pattern '", pattern, "' found in: ", results_dir)

  if (is.null(method_names)) {
    method_names <- basename(dirname(files))
  }
  if (length(method_names) != length(files))
    stop("Length of method_names (", length(method_names),
         ") must match number of files (", length(files), ").")

  df_list <- vector("list", length(files))
  for (i in seq_along(files)) {
    df <- .bm_read_tsv(files[i])
    df$Assay <- method_names[i]
    df_list[[i]] <- df
  }

  combined <- do.call(rbind, df_list)
  rownames(combined) <- NULL

  .bm_validate_confusion(combined)

  n_methods <- length(unique(combined$Assay))
  n_comps   <- length(unique(combined$Comparison))
  message("import_confusion_results: loaded ", n_methods, " method(s), ",
          n_comps, " comparison(s), ", nrow(combined), " rows.")

  combined
}


#' Import classified results from multiple benchmark result folders
#'
#' Reads \code{benchmark_classified.tsv} files from subdirectories,
#' adding an \code{Assay} column derived from the folder name.
#'
#' @param results_dir Parent directory containing benchmark subfolders
#' @param pattern Regex pattern for the classified filename
#' @param method_names Optional character vector of Assay names
#' @param recursive Search subdirectories recursively? (default: TRUE)
#'
#' @return data.frame with columns from classified_df plus Assay
#'
#' @export
import_classified_results <- function(results_dir,
                                      pattern      = "benchmark_classified\\.tsv$",
                                      method_names = NULL,
                                      recursive    = TRUE) {
  if (!dir.exists(results_dir))
    stop("Directory not found: ", results_dir)

  files <- list.files(results_dir, pattern = pattern,
                      full.names = TRUE, recursive = recursive)
  if (length(files) == 0)
    stop("No files matching pattern '", pattern, "' found in: ", results_dir)

  if (is.null(method_names)) {
    method_names <- basename(dirname(files))
  }
  if (length(method_names) != length(files))
    stop("Length of method_names (", length(method_names),
         ") must match number of files (", length(files), ").")

  df_list <- vector("list", length(files))
  for (i in seq_along(files)) {
    df <- .bm_read_tsv(files[i])
    df$Assay <- method_names[i]
    df_list[[i]] <- df
  }

  combined <- do.call(rbind, df_list)
  rownames(combined) <- NULL

  n_methods <- length(unique(combined$Assay))
  n_comps   <- length(unique(combined$Comparison))
  message("import_classified_results: loaded ", n_methods, " method(s), ",
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


#' Compute OpDEA ranking separately for each comparison
#'
#' For each comparison, ranks methods directly on their metric values
#' (no aggregation needed since there is one value per method per comparison).
#'
#' @param opdea_combined data.frame with columns: Assay, Comparison,
#'   nMCC, G_mean, pAUC_001, pAUC_005, pAUC_010
#' @param metrics Character vector of metrics to include in ranking
#'
#' @return Named list with:
#'   \describe{
#'     \item{by_comparison}{Named list where each key is a comparison and
#'       each value is a data.frame with metric values, ranks, and rank_final}
#'     \item{ranking_combined}{data.frame with all comparisons in long format
#'       (includes Comparison column)}
#'   }
#'
#' @examples
#' \dontrun{
#' by_comp <- bm_compute_ranking_by_comparison(opdea_all)
#' by_comp$by_comparison[["B-A"]]
#' by_comp$ranking_combined
#' }
#' @export
bm_compute_ranking_by_comparison <- function(opdea_combined,
                                             metrics = .BM_METRICS) {
  .bm_validate_opdea(opdea_combined)

  avail <- intersect(metrics, colnames(opdea_combined))
  if (length(avail) == 0)
    stop("None of the requested metrics found in opdea_combined.")
  metrics <- avail

  comps <- sort(unique(opdea_combined$Comparison))
  by_comp <- vector("list", length(comps))
  names(by_comp) <- comps

  combined_list <- vector("list", length(comps))

  for (i in seq_along(comps)) {
    comp <- comps[i]
    sub <- opdea_combined[opdea_combined$Comparison == comp,
                          c("Assay", metrics), drop = FALSE]

    # Rank directly (no aggregation — 1 row per method)
    rank_df <- .bm_rank_methods(sub, metrics)

    # Merge values + ranks
    merged <- merge(sub, rank_df, by = "Assay", sort = FALSE)
    merged <- merged[order(merged$rank_final), , drop = FALSE]
    rownames(merged) <- NULL

    by_comp[[comp]] <- merged

    combined_entry <- merged
    combined_entry$Comparison <- comp
    combined_list[[i]] <- combined_entry
  }

  ranking_combined <- do.call(rbind, combined_list)
  rownames(ranking_combined) <- NULL

  list(
    by_comparison    = by_comp,
    ranking_combined = ranking_combined
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

  # Use distinct shapes + colors for comparisons on the jittered points
  comp_levels <- sort(unique(long$Comparison))
  long$Comparison <- factor(long$Comparison, levels = comp_levels)

  # Select a colorblind-friendly palette
  if (n_comps <= 8) {
    comp_colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                     "#0072B2", "#D55E00", "#CC79A7", "#999999")[seq_len(n_comps)]
  } else {
    comp_colors <- grDevices::hcl.colors(n_comps, palette = "Dynamic")
  }
  names(comp_colors) <- comp_levels

  # Shapes: cycle through distinguishable filled shapes
  shape_pool <- c(16, 17, 15, 18, 8, 4, 3, 7, 9, 10, 12, 13, 14)
  comp_shapes <- shape_pool[((seq_len(n_comps) - 1) %% length(shape_pool)) + 1]
  names(comp_shapes) <- comp_levels

  # Pastel fill palette for boxplots (one per Assay), with matching darker borders
  assay_levels <- unique(long$Assay)
  n_assays <- length(assay_levels)
  assay_fills <- grDevices::hcl.colors(n_assays, palette = "Pastel 1")
  names(assay_fills) <- assay_levels
  # Derive darker border colors by reducing luminance
  assay_borders <- vapply(assay_fills, function(hex) {
    rgb_vals <- grDevices::col2rgb(hex)[, 1] / 255
    darker <- pmax(rgb_vals * 0.55, 0)
    grDevices::rgb(darker[1], darker[2], darker[3])
  }, character(1))
  names(assay_borders) <- assay_levels

  # Map border color per row for geom_boxplot
  long$assay_border <- assay_borders[as.character(long$Assay)]

  gg <- ggplot2::ggplot(long, ggplot2::aes(x = Assay, y = Value)) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = Assay),
                          alpha = 0.5, outlier.shape = NA, linewidth = 0.4,
                          show.legend = FALSE) +
    ggplot2::geom_point(
      ggplot2::aes(color = Comparison, shape = Comparison),
      position = ggplot2::position_jitter(width = 0.15, seed = 42),
      size = 2.2, alpha = 0.85
    ) +
    ggplot2::scale_fill_manual(values = assay_fills) +
    ggplot2::scale_color_manual(values = comp_colors) +
    ggplot2::scale_shape_manual(values = comp_shapes) +
    ggplot2::facet_wrap(~ Metric, scales = "free_y", ncol = 3) +
    ggplot2::labs(title = title, subtitle = subtitle,
                  x = NULL, y = "Value",
                  color = "Comparison", shape = "Comparison") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      axis.text.x      = ggplot2::element_text(angle = 60, hjust = 1, size = 8),
      plot.title        = ggplot2::element_text(face = "bold"),
      legend.position   = "bottom",
      legend.title      = ggplot2::element_text(face = "bold", size = 9),
      legend.text       = ggplot2::element_text(size = 8),
      strip.text        = ggplot2::element_text(face = "bold")
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(nrow = 1, override.aes = list(size = 3)),
      shape = ggplot2::guide_legend(nrow = 1)
    )

  gg
}


#' Heatmap of metric ranks for a single comparison
#'
#' @param ranking_by_comp List returned by \code{bm_compute_ranking_by_comparison()}
#' @param comparison Character string: which comparison to plot
#' @param title Optional plot title
#'
#' @return ggplot2 object
#' @export
bm_plot_ranking_heatmap_by_comp <- function(ranking_by_comp,
                                            comparison,
                                            title = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required.")
  if (!requireNamespace("tidyr", quietly = TRUE))
    stop("Package 'tidyr' is required.")

  if (!comparison %in% names(ranking_by_comp$by_comparison))
    stop("Comparison '", comparison, "' not found. Available: ",
         paste(names(ranking_by_comp$by_comparison), collapse = ", "))

  comp_df <- ranking_by_comp$by_comparison[[comparison]]

  rank_cols <- grep("^rank_", colnames(comp_df), value = TRUE)
  plot_data <- comp_df[, c("Assay", rank_cols), drop = FALSE]

  plot_data$Assay <- factor(plot_data$Assay, levels = rev(comp_df$Assay))

  long <- tidyr::pivot_longer(plot_data,
                              cols      = rank_cols,
                              names_to  = "Metric",
                              values_to = "Rank")

  long$Metric <- sub("^rank_", "", long$Metric)
  metric_order <- sub("^rank_", "", rank_cols)
  long$Metric <- factor(long$Metric, levels = metric_order)

  n_methods <- nrow(comp_df)

  title <- title %||% paste0("OpDEA Ranking: ", comparison)

  gg <- ggplot2::ggplot(long, ggplot2::aes(x = Metric, y = Assay, fill = Rank)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = round(Rank, 1)),
                       size = 3.5, color = "black") +
    ggplot2::scale_fill_gradient(low = "#2ca02c", high = "#d62728",
                                 limits = c(1, n_methods),
                                 name = "Rank") +
    ggplot2::labs(title = title, x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x    = ggplot2::element_text(angle = 45, hjust = 1),
      panel.grid     = ggplot2::element_blank(),
      plot.title     = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    )

  gg
}


#' Bar chart of final ranks for a single comparison
#'
#' @param ranking_by_comp List returned by \code{bm_compute_ranking_by_comparison()}
#' @param comparison Character string: which comparison to plot
#' @param title Optional plot title
#'
#' @return ggplot2 object
#' @export
bm_plot_ranking_bars_by_comp <- function(ranking_by_comp,
                                         comparison,
                                         title = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required.")

  if (!comparison %in% names(ranking_by_comp$by_comparison))
    stop("Comparison '", comparison, "' not found. Available: ",
         paste(names(ranking_by_comp$by_comparison), collapse = ", "))

  comp_df <- ranking_by_comp$by_comparison[[comparison]]
  n_methods <- nrow(comp_df)

  comp_df$Assay <- factor(comp_df$Assay, levels = rev(comp_df$Assay))

  title <- title %||% paste0("OpDEA Final Ranking: ", comparison)

  gg <- ggplot2::ggplot(comp_df,
                        ggplot2::aes(x = Assay, y = rank_final, fill = rank_final)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = round(rank_final, 1)),
                       hjust = -0.2, size = 3.5) +
    ggplot2::scale_fill_gradient(low = "#2ca02c", high = "#d62728",
                                 limits = c(1, n_methods),
                                 name = "Rank") +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
    ggplot2::labs(title = title, x = NULL, y = "Final Rank (lower = better)") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(face = "bold"),
      legend.position = "none"
    )

  gg
}


#' Stacked bar chart of TP/FP/FN/TN (OpDEA Figure 5 style)
#'
#' Horizontal stacked bars showing confusion matrix counts per method.
#' If \code{comparison} is NULL, generates a faceted plot with all comparisons.
#'
#' @param confusion_combined data.frame with columns:
#'   Assay, Comparison, TP, FP, TN, FN
#' @param comparison Optional: single comparison to plot.
#'   If NULL, faceted plot with all comparisons.
#' @param title Optional plot title
#'
#' @return ggplot2 object
#' @export
bm_plot_confusion_stacked <- function(confusion_combined,
                                      comparison = NULL,
                                      title      = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package 'ggplot2' is required.")
  if (!requireNamespace("tidyr", quietly = TRUE))
    stop("Package 'tidyr' is required.")

  .bm_validate_confusion(confusion_combined)

  # Filter to single comparison if requested
  if (!is.null(comparison)) {
    confusion_combined <- confusion_combined[
      confusion_combined$Comparison == comparison, , drop = FALSE]
    if (nrow(confusion_combined) == 0)
      stop("Comparison '", comparison, "' not found.")
  }

  # OpDEA colors: TP=blue, TN=coral, FP=teal, FN=gold
  conf_colors <- c(TP = "#5494cc", TN = "#e18283", FP = "#0d898a", FN = "#f9cc52")
  categories <- c("TP", "TN", "FP", "FN")

  # Pivot to long
  long <- tidyr::pivot_longer(
    confusion_combined[, c("Assay", "Comparison", categories), drop = FALSE],
    cols      = categories,
    names_to  = "Category",
    values_to = "Count"
  )
  long$Category <- factor(long$Category, levels = categories)

  # Order methods by TP count (descending) within each comparison
  # Use mean TP across comparisons for consistent ordering
  tp_order <- tapply(
    confusion_combined$TP,
    confusion_combined$Assay,
    mean, na.rm = TRUE
  )
  assay_order <- names(sort(tp_order, decreasing = FALSE))
  long$Assay <- factor(long$Assay, levels = assay_order)

  # Only show label if segment is large enough (>2% of total per bar)
  total_per_bar <- tapply(long$Count, list(long$Assay, long$Comparison),
                          sum, default = 0)
  long$label_text <- vapply(seq_len(nrow(long)), function(i) {
    total <- total_per_bar[as.character(long$Assay[i]),
                           as.character(long$Comparison[i])]
    if (!is.na(long$Count[i]) && long$Count[i] > total * 0.02) {
      as.character(long$Count[i])
    } else {
      ""
    }
  }, character(1))

  title <- if (!is.null(title)) {
    title
  } else if (!is.null(comparison)) {
    paste0("Confusion Matrix: ", comparison)
  } else {
    "Confusion Matrix by Method"
  }

  gg <- ggplot2::ggplot(long,
                        ggplot2::aes(x = Assay, y = Count,
                                     fill = Category)) +
    ggplot2::geom_col(width = 0.75, color = "white", linewidth = 0.3) +
    ggplot2::geom_text(
      ggplot2::aes(label = label_text),
      position = ggplot2::position_stack(vjust = 0.5),
      size = 3, color = "black", fontface = "bold"
    ) +
    ggplot2::scale_fill_manual(values = conf_colors,
                                name = NULL) +
    ggplot2::coord_flip() +
    ggplot2::labs(title = title, x = NULL, y = "Number of proteins") +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(face = "bold"),
      axis.line.y     = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.text     = ggplot2::element_text(size = 9),
      panel.grid.major.y = ggplot2::element_blank()
    ) +
    ggplot2::guides(fill = ggplot2::guide_legend(nrow = 1))

  # Facet if multiple comparisons
  if (is.null(comparison)) {
    gg <- gg + ggplot2::facet_wrap(~ Comparison, scales = "free_x")
  }

  gg
}


#' ROC Curves by Method for a Single Comparison
#'
#' Generates ROC (or ROC zoom) plot where each curve represents
#' a different method (Assay) for one comparison.
#'
#' @param classified_combined data.frame with columns: Assay, Comparison,
#'   truth, and a p-value column
#' @param comparison Character. Single comparison to plot (e.g. "B-A")
#' @param p_col P-value column name (default: "adj.P.Val")
#' @param zoom Logical. If TRUE, zoom to FPR 0-10% and show pAUC in legend
#' @param palette RColorBrewer palette name (default: "Set2")
#'
#' @return ggplot2 object or NULL if pROC not available
#'
#' @export
bm_plot_roc <- function(classified_combined,
                        comparison,
                        p_col   = "adj.P.Val",
                        zoom    = FALSE,
                        palette = "Set2") {

  if (!requireNamespace("pROC", quietly = TRUE)) {
    warning("Package 'pROC' not installed. Cannot generate ROC curves.\n",
            "Install with: install.packages('pROC')")
    return(NULL)
  }

  # Filter to the requested comparison
  df <- classified_combined[classified_combined$Comparison == comparison, , drop = FALSE]

  if (nrow(df) == 0) {
    warning("No data for comparison '", comparison, "'")
    return(NULL)
  }

  # Compute score: -log10(p-value)
  df$score <- -log10(pmax(as.numeric(df[[p_col]]), 1e-300))

  # Filter valid rows
  df <- df[!is.na(df$truth) & !is.na(df$score), , drop = FALSE]

  if (nrow(df) == 0) {
    warning("No valid data for ROC in comparison '", comparison, "'")
    return(NULL)
  }

  # Build ROC objects per Assay
  assay_order <- unique(df$Assay)
  df$Assay <- factor(df$Assay, levels = assay_order)
  assay_list <- split(df, df$Assay)

  roc_list <- lapply(assay_list, function(d) {
    if (length(unique(d$truth)) < 2 || nrow(d) < 10) return(NULL)
    tryCatch(
      pROC::roc(response = d$truth, predictor = d$score, quiet = TRUE),
      error = function(e) NULL
    )
  })
  roc_list <- Filter(Negate(is.null), roc_list)

  if (length(roc_list) == 0) {
    warning("Could not generate ROC curves for comparison '", comparison, "'")
    return(NULL)
  }

  # AUC / pAUC labels
  if (zoom) {
    aucs <- vapply(roc_list, function(r) {
      tryCatch(
        as.numeric(pROC::auc(r,
                              partial.auc = c(1, 0.9),
                              partial.auc.correct = TRUE)),
        error = function(e) NA_real_
      )
    }, numeric(1))
    labels <- paste0(names(roc_list), " (pAUC=", sprintf("%.3f", aucs), ")")
  } else {
    aucs <- vapply(roc_list, function(r) as.numeric(pROC::auc(r)), numeric(1))
    labels <- paste0(names(roc_list), " (AUC=", sprintf("%.3f", aucs), ")")
  }

  # Build plot
  title_text <- if (zoom) {
    paste0("ROC Zoom \u2014 ", comparison)
  } else {
    paste0("ROC \u2014 ", comparison)
  }

  gg <- pROC::ggroc(roc_list, legacy.axes = TRUE, linewidth = 1) +
    ggplot2::geom_abline(
      slope = 1, intercept = 0,
      linetype = "dashed", color = "gray50", linewidth = 0.5
    ) +
    ggplot2::labs(
      title    = title_text,
      subtitle = if (zoom) {
        "Low False Positive Rate region (0-10%)"
      } else {
        "Diagonal = random classifier"
      },
      x     = "False Positive Rate (1 - Specificity)",
      y     = "True Positive Rate (Sensitivity)",
      color = "Method"
    ) +
    ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(hjust = 0.5, face = "bold",
                                            size = 15, color = "#1D3557"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 11,
                                            color = "#6C757D"),
      legend.position = "right",
      legend.text     = ggplot2::element_text(size = 10),
      legend.title    = ggplot2::element_text(face = "bold"),
      panel.grid.minor = ggplot2::element_blank()
    )

  # Apply palette — use hcl.colors to support any number of methods
  n_methods <- length(roc_list)
  pal_colors <- if (requireNamespace("RColorBrewer", quietly = TRUE) &&
                    n_methods <= RColorBrewer::brewer.pal.info[palette, "maxcolors"]) {
    RColorBrewer::brewer.pal(max(n_methods, 3), palette)[seq_len(n_methods)]
  } else {
    grDevices::hcl.colors(n_methods, palette = "Dynamic")
  }
  gg <- gg + ggplot2::scale_color_manual(values = pal_colors, labels = labels)

  # Zoom or full view
  if (zoom) {
    gg <- gg + ggplot2::coord_cartesian(xlim = c(0, 0.1), ylim = c(0, 1))
  } else {
    gg <- gg + ggplot2::coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
  }

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
#' @param confusion_combined Optional pre-built data.frame with columns:
#'   Assay, Comparison, TP, FP, TN, FN. If NULL and \code{results_dir}
#'   is provided, imported from \code{benchmark_confusion_overall.tsv} files.
#' @param classified_combined Optional pre-built data.frame with columns:
#'   Assay, Comparison, truth, and a p-value column. If NULL and
#'   \code{results_dir} is provided, imported from
#'   \code{benchmark_classified.tsv} files. Used for multi-method ROC curves.
#' @param p_col P-value column name for ROC curves (default: "adj.P.Val")
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
#'     \item{ranking_by_comparison}{Named list of ranking data.frames per comparison}
#'     \item{ranking_by_comparison_combined}{data.frame with all per-comparison rankings}
#'     \item{gg_ranking_heatmap_by_comp}{Named list of heatmaps per comparison}
#'     \item{gg_ranking_bars_by_comp}{Named list of bar charts per comparison}
#'     \item{confusion_combined}{Combined confusion data.frame (if available)}
#'     \item{gg_confusion_stacked}{Faceted confusion stacked bars (if available)}
#'     \item{gg_confusion_stacked_by_comp}{Named list of confusion plots per comparison}
#'     \item{classified_combined}{Combined classified data.frame (if available)}
#'     \item{gg_roc_by_comp}{Named list of ROC plots per comparison (if available)}
#'     \item{gg_roc_zoom_by_comp}{Named list of ROC zoom plots per comparison (if available)}
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
benchmarking_multiple <- function(opdea_combined       = NULL,
                                  confusion_combined   = NULL,
                                  classified_combined  = NULL,
                                  results_dir          = NULL,
                                  pattern              = "benchmark_opdea_metrics\\.tsv$",
                                  method_names         = NULL,
                                  recursive            = TRUE,
                                  metrics              = .BM_METRICS,
                                  p_col                = "adj.P.Val",
                                  plots                = "all",
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

  # --- Import confusion data (optional) ---
  has_confusion <- FALSE
  if (is.null(confusion_combined) && !is.null(results_dir)) {
    confusion_combined <- tryCatch({
      if (verbose) message("  Importing confusion data ...")
      import_confusion_results(results_dir,
                               method_names = method_names,
                               recursive    = recursive)
    }, error = function(e) {
      if (verbose) message("  No confusion data found (skipping): ", e$message)
      NULL
    })
  }
  if (!is.null(confusion_combined)) {
    tryCatch({
      .bm_validate_confusion(confusion_combined)
      has_confusion <- TRUE
    }, error = function(e) {
      warning("Invalid confusion_combined (skipping): ", e$message)
      confusion_combined <- NULL
    })
  }

  # --- Import classified data (optional, for ROC curves) ---
  has_classified <- FALSE
  if (is.null(classified_combined) && !is.null(results_dir)) {
    classified_combined <- tryCatch({
      if (verbose) message("  Importing classified data ...")
      import_classified_results(results_dir,
                                method_names = method_names,
                                recursive    = recursive)
    }, error = function(e) {
      if (verbose) message("  No classified data found (skipping): ", e$message)
      NULL
    })
  }
  if (!is.null(classified_combined)) {
    tryCatch({
      .bm_validate_classified(classified_combined, p_col)
      has_classified <- TRUE
    }, error = function(e) {
      warning("Invalid classified_combined (skipping): ", e$message)
      classified_combined <- NULL
    })
  }

  n_methods <- length(unique(opdea_combined$Assay))
  n_comps   <- length(unique(opdea_combined$Comparison))

  if (verbose) {
    message("=== MULTIPLE-METHOD OPDEA BENCHMARKING ===")
    message("  Methods: ", n_methods,
            " | Comparisons: ", n_comps,
            " | Rows: ", nrow(opdea_combined))
    if (has_confusion) message("  Confusion data: available")
    if (has_classified) message("  Classified data: available (ROC curves)")
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
  # STEP 2b: Compute ranking by comparison
  # ==========================================================================
  if (verbose) message("  Computing ranking by comparison ...")
  ranking_by_comp <- bm_compute_ranking_by_comparison(opdea_combined, metrics)
  comps <- names(ranking_by_comp$by_comparison)

  if (verbose) {
    for (comp in comps) {
      best <- ranking_by_comp$by_comparison[[comp]]$Assay[1]
      best_rf <- ranking_by_comp$by_comparison[[comp]]$rank_final[1]
      message("    ", comp, ": best = ", best, " (rank_final = ", best_rf, ")")
    }
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
    message("  ", n_ok, " global plot(s) generated",
            if (n_fail > 0) paste0(", ", n_fail, " failed") else ".")
  }

  # --- Per-comparison plots ---
  gg_heatmap_by_comp <- vector("list", length(comps))
  names(gg_heatmap_by_comp) <- comps
  gg_bars_by_comp <- vector("list", length(comps))
  names(gg_bars_by_comp) <- comps

  for (comp in comps) {
    if (verbose) message("  Generating plots for comparison: ", comp, " ...")
    gg_heatmap_by_comp[[comp]] <- tryCatch(
      bm_plot_ranking_heatmap_by_comp(ranking_by_comp, comp),
      error = function(e) {
        warning("bm_plot_ranking_heatmap_by_comp(", comp, ") failed: ",
                conditionMessage(e))
        NULL
      }
    )
    gg_bars_by_comp[[comp]] <- tryCatch(
      bm_plot_ranking_bars_by_comp(ranking_by_comp, comp),
      error = function(e) {
        warning("bm_plot_ranking_bars_by_comp(", comp, ") failed: ",
                conditionMessage(e))
        NULL
      }
    )
  }

  n_comp_plots <- sum(!vapply(gg_heatmap_by_comp, is.null, logical(1))) +
                  sum(!vapply(gg_bars_by_comp, is.null, logical(1)))
  if (verbose) message("  ", n_comp_plots, " per-comparison plot(s) generated.")

  # --- Confusion stacked bar plots ---
  gg_confusion_stacked <- NULL
  gg_confusion_by_comp <- list()

  if (has_confusion) {
    if (verbose) message("  Generating confusion stacked bar plots ...")

    # Faceted (all comparisons)
    gg_confusion_stacked <- tryCatch(
      bm_plot_confusion_stacked(confusion_combined),
      error = function(e) {
        warning("bm_plot_confusion_stacked() failed: ", conditionMessage(e))
        NULL
      }
    )

    # Per comparison
    conf_comps <- unique(confusion_combined$Comparison)
    gg_confusion_by_comp <- vector("list", length(conf_comps))
    names(gg_confusion_by_comp) <- conf_comps

    for (comp in conf_comps) {
      gg_confusion_by_comp[[comp]] <- tryCatch(
        bm_plot_confusion_stacked(confusion_combined, comparison = comp),
        error = function(e) {
          warning("bm_plot_confusion_stacked(", comp, ") failed: ",
                  conditionMessage(e))
          NULL
        }
      )
    }

    n_conf_ok <- (!is.null(gg_confusion_stacked)) +
                 sum(!vapply(gg_confusion_by_comp, is.null, logical(1)))
    if (verbose) message("  ", n_conf_ok, " confusion plot(s) generated.")
  }

  # --- ROC curves per comparison (one curve per method) ---
  gg_roc_by_comp      <- list()
  gg_roc_zoom_by_comp <- list()

  if (has_classified) {
    if (verbose) message("  Generating ROC curves per comparison ...")

    roc_comps <- unique(classified_combined$Comparison)
    gg_roc_by_comp      <- vector("list", length(roc_comps))
    names(gg_roc_by_comp) <- roc_comps
    gg_roc_zoom_by_comp <- vector("list", length(roc_comps))
    names(gg_roc_zoom_by_comp) <- roc_comps

    for (comp in roc_comps) {
      gg_roc_by_comp[[comp]] <- tryCatch(
        bm_plot_roc(classified_combined, comp, p_col = p_col, zoom = FALSE),
        error = function(e) {
          warning("bm_plot_roc(", comp, ") failed: ", conditionMessage(e))
          NULL
        }
      )
      gg_roc_zoom_by_comp[[comp]] <- tryCatch(
        bm_plot_roc(classified_combined, comp, p_col = p_col, zoom = TRUE),
        error = function(e) {
          warning("bm_plot_roc(", comp, ", zoom) failed: ", conditionMessage(e))
          NULL
        }
      )
    }

    n_roc_ok <- sum(!vapply(gg_roc_by_comp, is.null, logical(1))) +
                sum(!vapply(gg_roc_zoom_by_comp, is.null, logical(1)))
    if (verbose) message("  ", n_roc_ok, " ROC plot(s) generated.")
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

      # Per-comparison tables
      comp_dir <- file.path(output_dir, "by_comparison")
      if (!dir.exists(comp_dir)) dir.create(comp_dir, recursive = TRUE)

      .bm_export_data(ranking_by_comp$ranking_combined,
                      file.path(output_dir, "bm_multiple_ranking_by_comparison.tsv"))

      for (comp in comps) {
        safe_comp <- gsub("[^A-Za-z0-9_-]", "_", comp)
        .bm_export_data(ranking_by_comp$by_comparison[[comp]],
                        file.path(comp_dir, paste0("bm_ranking_", safe_comp, ".tsv")))
      }

      n_tables <- 5 + 1 + length(comps)
      if (verbose) message("  Exported ", n_tables, " TSV files.")
    }

    if (export_plots) {
      for (nm in selected) {
        if (!is.null(plot_result[[nm]])) {
          .bm_export_gg_plot(plot_result[[nm]],
                             file.path(output_dir, paste0("bm_multiple_", nm, ".png")),
                             plot_width, plot_height, plot_dpi)
        }
      }

      # Per-comparison plots
      comp_dir <- file.path(output_dir, "by_comparison")
      if (!dir.exists(comp_dir)) dir.create(comp_dir, recursive = TRUE)

      for (comp in comps) {
        safe_comp <- gsub("[^A-Za-z0-9_-]", "_", comp)
        if (!is.null(gg_heatmap_by_comp[[comp]])) {
          .bm_export_gg_plot(gg_heatmap_by_comp[[comp]],
                             file.path(comp_dir, paste0("bm_ranking_heatmap_", safe_comp, ".png")),
                             plot_width, plot_height, plot_dpi)
        }
        if (!is.null(gg_bars_by_comp[[comp]])) {
          .bm_export_gg_plot(gg_bars_by_comp[[comp]],
                             file.path(comp_dir, paste0("bm_ranking_bars_", safe_comp, ".png")),
                             plot_width, plot_height, plot_dpi)
        }
      }

      if (verbose) message("  Exported ", n_ok, " global + ",
                           n_comp_plots, " per-comparison PNG files.")
    }

    # Confusion exports
    if (has_confusion) {
      if (export_tables) {
        .bm_export_data(confusion_combined,
                        file.path(output_dir, "bm_multiple_confusion_combined.tsv"))
      }
      if (export_plots) {
        if (!is.null(gg_confusion_stacked)) {
          .bm_export_gg_plot(gg_confusion_stacked,
                             file.path(output_dir, "bm_multiple_confusion_stacked.png"),
                             plot_width, plot_height, plot_dpi)
        }
        comp_dir <- file.path(output_dir, "by_comparison")
        if (!dir.exists(comp_dir)) dir.create(comp_dir, recursive = TRUE)
        for (comp in names(gg_confusion_by_comp)) {
          if (!is.null(gg_confusion_by_comp[[comp]])) {
            safe_comp <- gsub("[^A-Za-z0-9_-]", "_", comp)
            .bm_export_gg_plot(gg_confusion_by_comp[[comp]],
                               file.path(comp_dir,
                                         paste0("bm_confusion_stacked_", safe_comp, ".png")),
                               plot_width, plot_height, plot_dpi)
          }
        }
      }
      if (verbose) message("  Exported confusion data and plots.")
    }

    # ROC exports
    if (has_classified) {
      if (export_tables) {
        .bm_export_data(classified_combined,
                        file.path(output_dir, "bm_multiple_classified_combined.tsv"))
      }
      if (export_plots) {
        comp_dir <- file.path(output_dir, "by_comparison")
        if (!dir.exists(comp_dir)) dir.create(comp_dir, recursive = TRUE)
        for (comp in names(gg_roc_by_comp)) {
          safe_comp <- gsub("[^A-Za-z0-9_-]", "_", comp)
          if (!is.null(gg_roc_by_comp[[comp]])) {
            .bm_export_gg_plot(gg_roc_by_comp[[comp]],
                               file.path(comp_dir, paste0("bm_roc_", safe_comp, ".png")),
                               plot_width, plot_height, plot_dpi)
          }
          if (!is.null(gg_roc_zoom_by_comp[[comp]])) {
            .bm_export_gg_plot(gg_roc_zoom_by_comp[[comp]],
                               file.path(comp_dir, paste0("bm_roc_zoom_", safe_comp, ".png")),
                               plot_width, plot_height, plot_dpi)
          }
        }
      }
      if (verbose) message("  Exported classified data and ROC plots.")
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
    median_ranking    = ranking$median_ranking,
    ranking_by_comparison          = ranking_by_comp$by_comparison,
    ranking_by_comparison_combined = ranking_by_comp$ranking_combined
  )

  # Add global plots with gg_ prefix
  for (nm in selected) {
    result[[paste0("gg_", nm)]] <- plot_result[[nm]]
  }

  # Add per-comparison plots
  result$gg_ranking_heatmap_by_comp <- gg_heatmap_by_comp
  result$gg_ranking_bars_by_comp    <- gg_bars_by_comp

  # Confusion data and plots
  if (has_confusion) {
    result$confusion_combined              <- confusion_combined
    result$gg_confusion_stacked            <- gg_confusion_stacked
    result$gg_confusion_stacked_by_comp    <- gg_confusion_by_comp
  }

  # Classified data and ROC plots
  if (has_classified) {
    result$classified_combined     <- classified_combined
    result$gg_roc_by_comp          <- gg_roc_by_comp
    result$gg_roc_zoom_by_comp     <- gg_roc_zoom_by_comp
  }

  result$parameters <- list(
    metrics       = metrics,
    n_methods     = n_methods,
    n_comparisons = n_comps,
    comparisons   = comps,
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
#   list(norm = "cycloess", imp = "none"),
#   list(norm = "log2Norm", imp = "combo", mar = "Impseqrob", mnar = "MinDet")
# )
#
# # Helper: build assay name from combo
# .make_assay_name <- function(combo) {
#   if (combo$imp == "none") {
#     paste0(combo$norm, "_none")
#   } else {
#     paste(combo$norm, combo$mar, combo$mnar, sep = "_")
#   }
# }
#
# # Run processing + benchmarking for each combo
# opdea_list      <- list()
# confusion_list  <- list()
# classified_list <- list()
# for (combo in combos) {
#   # Build process_proteomics args (omit mar/mnar when imp = "none")
#   proc_args <- list(
#     preprocessing = preprocessing,
#     norm_method   = combo$norm,
#     imp_method    = combo$imp
#   )
#   if (combo$imp != "none") {
#     proc_args$mar_method  <- combo$mar
#     proc_args$mnar_method <- combo$mnar
#   }
#   result <- do.call(process_proteomics, proc_args)
#
#   assay_name <- .make_assay_name(combo)
#
#   bench <- benchmarking_proteomics(
#     de_res          = result$DEPs_results,
#     species_df      = species_df,
#     expected_values = expected,
#     alpha           = 0.05,
#     output_dir      = paste0("results/benchmark_", assay_name)
#   )
#
#   opdea <- bench$opdea_metrics
#   opdea$Assay <- assay_name
#   opdea_list[[length(opdea_list) + 1]] <- opdea
#
#   conf <- bench$confusion_overall
#   conf$Assay <- assay_name
#   confusion_list[[length(confusion_list) + 1]] <- conf
#
#   cls <- bench$classified_df
#   cls$Assay <- assay_name
#   classified_list[[length(classified_list) + 1]] <- cls
# }
#
# opdea_all      <- do.call(rbind, opdea_list)
# confusion_all  <- do.call(rbind, confusion_list)
# classified_all <- do.call(rbind, classified_list)
#
# bm_result <- benchmarking_multiple(
#   opdea_combined      = opdea_all,
#   confusion_combined  = confusion_all,
#   classified_combined = classified_all,
#   output_dir          = "results/bm_multiple",
#   verbose             = TRUE
# )
#
# bm_result$mean_ranking
# bm_result$gg_ranking_heatmap_mean
# bm_result$gg_roc_by_comp[["B-A"]]
# bm_result$gg_roc_zoom_by_comp[["B-A"]]
#
#
# --- Option B: Import from files ---
#
# source("R/Benchmarking_Multiple.R")
#
# # Assumes folders like:
# #   results/benchmark_cycloess_Impseqrob_min/benchmark_opdea_metrics.tsv
# #   results/benchmark_cycloess_Impseqrob_min/benchmark_confusion_overall.tsv
# #   results/benchmark_cycloess_Impseqrob_min/benchmark_classified.tsv
# #   results/benchmark_quantile_knn_min/...
# #   results/benchmark_log2Norm_Impseqrob_MinDet/...
#
# bm_result <- benchmarking_multiple(
#   results_dir = "results",
#   output_dir  = "results/bm_multiple",
#   verbose     = TRUE
# )
#
# bm_result$mean_ranking
# bm_result$gg_ranking_bars_mean
# bm_result$gg_roc_by_comp[["B-A"]]
# bm_result$gg_roc_zoom_by_comp[["B-A"]]
