# =============================================================================
# Pattern Profiler Highcharts: cluster visualisation
# =============================================================================
#
# This script builds interactive Highcharts visualisations for the clustering
# results produced by Pattern_Profiler_Analysis.R.
#
# Input:
#   - Pattern_Profiler_Input.parquet: table in LONG format with:
#     - FeatureID, Cluster, Membership, z-scores per condition
#
# Output:
#   - Interactive Highcharts plots (profiles, centroids)
#
# Author: Sergio Ciordia
# License: GPL-3
# =============================================================================

# -----------------------------------------------------------------------------
# DEPENDENCIES
# -----------------------------------------------------------------------------



# =============================================================================
# COLOUR HELPER FUNCTIONS
# =============================================================================

#' Configure the colour palette for the clusters
#'
#' @param n_clusters Number of clusters
#' @param palette Palette to use (NULL, "ggsci::name", "brewer:name")
#'
#' @return Vector of colours named by cluster
configure_cluster_palette <- function(n_clusters, palette = NULL) {

  # Default palette: distinctive and accessible colours
  default_colors <- c(
    "#E63946", "#457B9D", "#2A9D8F", "#E9C46A", "#F4A261",
    "#264653", "#A8DADC", "#1D3557", "#F77F00", "#D62828",
    "#023E8A", "#0077B6", "#00B4D8", "#90E0EF", "#CAF0F8"
  )

  if (is.null(palette)) {
    # If more clusters are requested than there are base colours, interpolate
    # instead of recycling (this avoids the length mismatch when names() is
    # assigned further down).
    if (n_clusters > length(default_colors)) {
      colors <- colorRampPalette(default_colors)(n_clusters)
    } else {
      colors <- default_colors[seq_len(n_clusters)]
    }

  } else if (grepl("^ggsci::", palette)) {
    # ggsci palettes via paletteer
    if (requireNamespace("paletteer", quietly = TRUE)) {
      pal_name <- sub("^ggsci::", "", palette)
      colors <- tryCatch({
        as.character(paletteer::paletteer_d(paste0("ggsci::", pal_name), n_clusters))
      }, error = function(e) {
        warning("ggsci palette not found, using the default one")
        default_colors[seq_len(n_clusters)]
      })
    } else {
      warning("Package 'paletteer' not available, using the default palette")
      colors <- default_colors[seq_len(n_clusters)]
    }

  } else if (grepl("^brewer:", palette)) {
    # RColorBrewer palettes
    pal_name <- sub("^brewer:", "", palette)
    if (requireNamespace("RColorBrewer", quietly = TRUE)) {
      max_colors <- RColorBrewer::brewer.pal.info[pal_name, "maxcolors"]
      if (is.na(max_colors)) {
        warning("brewer palette not found, using the default one")
        colors <- default_colors[seq_len(n_clusters)]
      } else {
        colors <- RColorBrewer::brewer.pal(min(n_clusters, max_colors), pal_name)
        if (n_clusters > max_colors) {
          colors <- colorRampPalette(colors)(n_clusters)
        }
      }
    } else {
      warning("Package 'RColorBrewer' not available, using the default palette")
      colors <- default_colors[seq_len(n_clusters)]
    }

  } else {
    # Assume a vector of colours
    if (length(palette) >= n_clusters) {
      colors <- palette[seq_len(n_clusters)]
    } else {
      colors <- colorRampPalette(palette)(n_clusters)
    }
  }

  # Name them by cluster
  names(colors) <- seq_len(n_clusters)
  colors
}


# =============================================================================
# DATA READING
# =============================================================================

#' Read Pattern Profiler data from a parquet file
#'
#' @param file_path Path to the parquet file
#' @param min_membership Optional filter by minimum membership
#'
#' @return DataFrame with the clustering data
#'
#' @examples
#' if (requireNamespace("arrow", quietly = TRUE) &&
#'     requireNamespace("Mfuzz", quietly = TRUE) &&
#'     requireNamespace("e1071", quietly = TRUE)) {
#'   data(nadia_dia)
#'   res <- process_proteomics(nadia_dia, verbose = FALSE)
#'
#'   # Write the clustering to parquet, then read it back
#'   f <- file.path(tempdir(), "pattern_profiler.parquet")
#'   pattern_profiler_analysis(res$se_proc, res$DEPs_results,
#'                             assay_name = "Impseqrob_min",
#'                             auto_select_c = FALSE, c = 3,
#'                             output_file = f, verbose = FALSE)
#'
#'   pp_data <- read_pattern_profiler_data(f, min_membership = 0.5)
#'   unlink(f)
#'   print(head(pp_data))
#' }
#'
#' @export
read_pattern_profiler_data <- function(file_path, min_membership = NULL) {

  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }

  data <- arrow::read_parquet(file_path)
  data <- as.data.frame(data)

  # Validate the required columns
  required_cols <- c("FeatureID", "Cluster", "Membership")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing columns: ", paste(missing_cols, collapse = ", "))
  }

  # Filter by membership if requested
  if (!is.null(min_membership)) {
    n_before <- nrow(data)
    data <- data[data$Membership >= min_membership, ]
    n_after <- nrow(data)
    message(sprintf("Filtered by membership >= %.2f: %d -> %d rows",
                    min_membership, n_before, n_after))
  }

  data
}


#' Detect the condition columns (z-scores)
#'
#' @param data Pattern Profiler DataFrame
#' @return Vector with the names of the condition columns
detect_condition_columns <- function(data) {
  # Exclude the known columns
  exclude_cols <- c("FeatureID", "Cluster", "Membership")
  all_cols <- names(data)

  condition_cols <- setdiff(all_cols, exclude_cols)

  if (length(condition_cols) == 0) {
    stop("No condition columns (z-scores) found")
  }

  condition_cols
}


# =============================================================================
# VISUALISATION FUNCTIONS
# =============================================================================

#' Cluster profile plot with Highcharts
#'
#' Builds an interactive plot showing the expression profiles of the proteins in
#' a given cluster.
#'
#' @param data Pattern Profiler DataFrame (LONG format)
#' @param cluster Number of the cluster to plot
#' @param conditions Vector of condition names (order for the X axis)
#' @param min_membership Additional membership filter (NULL = no filter)
#' @param show_centroid Show the centroid line (default: TRUE)
#' @param centroid_summary Method for the centroid: "mean" or "median"
#' @param cluster_color Cluster colour (NULL = automatic)
#' @param line_width Width of the profile lines (default: 1)
#' @param line_opacity Opacity of the lines (default: 0.4)
#' @param centroid_width Width of the centroid line (default: 3)
#' @param title Custom title (optional)
#' @param height Plot height in pixels
#'
#' @return highchart object
#'
#' @examples
#' if (requireNamespace("Mfuzz", quietly = TRUE) &&
#'     requireNamespace("e1071", quietly = TRUE)) {
#'   data(nadia_dia)
#'   res <- process_proteomics(nadia_dia, verbose = FALSE)
#'   pp <- pattern_profiler_analysis(res$se_proc, res$DEPs_results,
#'                                   assay_name = "Impseqrob_min",
#'                                   auto_select_c = FALSE, c = 3,
#'                                   verbose = FALSE)
#'
#'   hc <- cluster_profile_highchart(pp$long_output, cluster = 1,
#'                                   conditions = c("A", "B", "D"))
#'   print(class(hc))
#' }
#'
#' @export
cluster_profile_highchart <- function(data,
                                       cluster,
                                       conditions = NULL,
                                       min_membership = NULL,
                                       show_centroid = TRUE,
                                       centroid_summary = c("mean", "median"),
                                       cluster_color = NULL,
                                       line_width = 1,
                                       line_opacity = 0.4,
                                       centroid_width = 3,
                                       title = NULL,
                                       height = NULL) {

  centroid_summary <- match.arg(centroid_summary)

  # Detect the conditions if they are not supplied
  if (is.null(conditions)) {
    conditions <- detect_condition_columns(data)
  }

  # Filter by cluster
  cluster_data <- data[data$Cluster == cluster, ]

  if (nrow(cluster_data) == 0) {
    warning(sprintf("No data for cluster %d", cluster))
    return(NULL)
  }

  # Apply the additional membership filter
  if (!is.null(min_membership)) {
    cluster_data <- cluster_data[cluster_data$Membership >= min_membership, ]
  }

  if (nrow(cluster_data) == 0) {
    warning(sprintf("Cluster %d: no proteins with membership >= %.2f",
                    cluster, min_membership))
    return(NULL)
  }

  # Extract the z-scores and the metadata
  zscores <- as.matrix(cluster_data[, conditions, drop = FALSE])
  feature_ids <- cluster_data$FeatureID
  n_proteins <- nrow(cluster_data)

  # Compute the centroid
  if (centroid_summary == "mean") {
    centroid <- colMeans(zscores, na.rm = TRUE)
  } else {
    centroid <- apply(zscores, 2, median, na.rm = TRUE)
  }
  centroid <- unname(as.numeric(centroid))

  # Configure the cluster colour
  if (is.null(cluster_color)) {
    n_clusters <- max(data$Cluster)
    palette <- configure_cluster_palette(n_clusters)
    cluster_color <- unname(palette[as.character(cluster)])
  }

  # Line colour with opacity
  line_color <- .hex_to_rgba(cluster_color, line_opacity)

  # Title
  if (is.null(title)) {
    min_mem <- min(cluster_data$Membership)
    title <- sprintf("Cluster %d (n = %d, membership >= %.2f)",
                     cluster, n_proteins, min_mem)
  }

  # ---------------------------------------------------------------------------
  # Build the profile line series (non-interactive)
  # ---------------------------------------------------------------------------
  profile_series <- lapply(seq_len(n_proteins), function(i) {
    zscore_row <- unname(as.numeric(zscores[i, ]))
    feature_id <- as.character(feature_ids[i])

    points <- lapply(seq_along(conditions), function(j) {
      list(
        x = as.integer(j - 1),
        y = round(zscore_row[j], 4)
      )
    })

    list(
      name = feature_id,
      type = "line",
      data = points,
      color = line_color,
      lineWidth = line_width,
      marker = list(enabled = FALSE),
      enableMouseTracking = FALSE,
      showInLegend = FALSE
    )
  })

  # ---------------------------------------------------------------------------
  # Centroid series (interactive)
  # ---------------------------------------------------------------------------
  centroid_series <- NULL
  if (show_centroid) {
    centroid_points <- lapply(seq_along(conditions), function(j) {
      list(
        x = as.integer(j - 1),
        y = round(centroid[j], 4),
        condition = unname(conditions[j])
      )
    })

    # Clean colours without names
    centroid_line_color <- .darken_hex(cluster_color, 0.2)
    centroid_marker_line <- .darken_hex(cluster_color, 0.3)

    centroid_series <- list(
      name = paste0("Centroid (", centroid_summary, ")"),
      type = "line",
      data = centroid_points,
      color = centroid_line_color,
      lineWidth = centroid_width,
      marker = list(
        enabled = TRUE,
        symbol = "circle",
        radius = 5,
        fillColor = cluster_color,
        lineColor = centroid_marker_line,
        lineWidth = 2
      ),
      zIndex = 10,
      showInLegend = TRUE,
      enableMouseTracking = TRUE
    )
  }

  # ---------------------------------------------------------------------------
  # Build the highchart
  # ---------------------------------------------------------------------------
  hc <- highchart() |>
    hc_chart(
      type = "line",
      backgroundColor = "#FFFFFF",
      style = list(fontFamily = "Inter, -apple-system, sans-serif"),
      height = height,
      zoomType = "xy"
    ) |>
    hc_title(
      text = title,
      style = list(
        fontSize = "18px",
        fontWeight = "600",
        color = "#1D3557"
      )
    ) |>
    hc_xAxis(
      categories = as.list(unname(conditions)),
      title = list(
        text = "Condition",
        style = list(
          fontSize = "13px",
          fontWeight = "bold",
          color = "#212529"
        )
      ),
      labels = list(
        style = list(
          fontSize = "12px",
          color = "#495057"
        )
      ),
      lineColor = "#DEE2E6",
      tickColor = "#DEE2E6",
      gridLineWidth = 0
    ) |>
    hc_yAxis(
      title = list(
        text = "z-score",
        style = list(
          fontSize = "13px",
          fontWeight = "bold",
          color = "#212529"
        )
      ),
      labels = list(
        style = list(
          fontSize = "11px",
          color = "#495057"
        )
      ),
      lineColor = "#DEE2E6",
      lineWidth = 1,
      gridLineColor = "#F1F3F4",
      gridLineDashStyle = "Dot",
      plotLines = list(
        list(
          value = 0,
          color = "#ADB5BD",
          width = 1,
          dashStyle = "Dash",
          zIndex = 1
        )
      )
    ) |>
    hc_tooltip(
      useHTML = TRUE,
      backgroundColor = "rgba(255, 255, 255, 0.95)",
      borderColor = "#DEE2E6",
      borderRadius = 8,
      shadow = TRUE,
      style = list(fontSize = "12px"),
      headerFormat = "",
      pointFormat = paste0(
        "<div style='padding: 4px;'>",
        "<b style='font-size: 13px; color: #1D3557;'>{series.name}</b><br/>",
        "<span style='color: #6C757D;'>Condition:</span> <b>{point.condition}</b><br/>",
        "<span style='color: #6C757D;'>z-score:</span> <b>{point.y:.3f}</b>",
        "</div>"
      )
    ) |>
    hc_legend(
      enabled = show_centroid,
      layout = "horizontal",
      align = "center",
      verticalAlign = "bottom",
      itemStyle = list(
        fontSize = "12px",
        fontWeight = "normal",
        color = "#495057"
      )
    ) |>
    hc_exporting(
      enabled = TRUE,
      buttons = list(
        contextButton = list(
          menuItems = c("downloadPNG", "downloadSVG", "downloadPDF")
        )
      )
    )

  # Add the profile series
  for (series in profile_series) {
    hc <- hc |> hc_add_series(
      name = series$name,
      type = series$type,
      data = series$data,
      color = series$color,
      lineWidth = series$lineWidth,
      marker = series$marker,
      enableMouseTracking = series$enableMouseTracking,
      showInLegend = series$showInLegend
    )
  }

  # Add the centroid
  if (!is.null(centroid_series)) {
    hc <- hc |> hc_add_series(
      name = centroid_series$name,
      type = centroid_series$type,
      data = centroid_series$data,
      color = centroid_series$color,
      lineWidth = centroid_series$lineWidth,
      marker = centroid_series$marker,
      zIndex = centroid_series$zIndex,
      showInLegend = centroid_series$showInLegend,
      enableMouseTracking = centroid_series$enableMouseTracking
    )
  }

  hc
}


#' List of cluster profile plots with Highcharts
#'
#' Builds a list of interactive plots for all the clusters.
#'
#' @param data Pattern Profiler DataFrame (LONG format)
#' @param conditions Vector of condition names
#' @param clusters Vector of clusters to plot (NULL = all)
#' @param min_membership Minimum membership filter
#' @param show_centroid Show the centroid line (default: TRUE)
#' @param centroid_summary Method for the centroid: "mean" or "median"
#' @param palette Palette for the cluster colours
#' @param line_width Width of the profile lines (default: 1)
#' @param line_opacity Opacity of the lines (default: 0.4)
#' @param centroid_width Width of the centroid line (default: 3)
#' @param height Height of each plot in pixels
#'
#' @return Named list of highchart objects
#'
#' @examples
#' if (requireNamespace("Mfuzz", quietly = TRUE) &&
#'     requireNamespace("e1071", quietly = TRUE)) {
#'   data(nadia_dia)
#'   res <- process_proteomics(nadia_dia, verbose = FALSE)
#'   pp <- pattern_profiler_analysis(res$se_proc, res$DEPs_results,
#'                                   assay_name = "Impseqrob_min",
#'                                   auto_select_c = FALSE, c = 3,
#'                                   verbose = FALSE)
#'
#'   hc_profiles <- cluster_profile_highchart_list(pp$long_output,
#'                                                 conditions = c("A", "B", "D"))
#'   print(names(hc_profiles))
#' }
#'
#' @export
cluster_profile_highchart_list <- function(data,
                                            conditions = NULL,
                                            clusters = NULL,
                                            min_membership = NULL,
                                            show_centroid = TRUE,
                                            centroid_summary = c("mean", "median"),
                                            palette = NULL,
                                            line_width = 1,
                                            line_opacity = 0.4,
                                            centroid_width = 3,
                                            height = NULL) {

  centroid_summary <- match.arg(centroid_summary)

  # Detect the conditions if they are not supplied
  if (is.null(conditions)) {
    conditions <- detect_condition_columns(data)
  }

  # Determine which clusters to plot
  if (is.null(clusters)) {
    clusters <- sort(unique(data$Cluster))
  }

  n_clusters <- max(data$Cluster)

  # Configure the colour palette
  cluster_colors <- configure_cluster_palette(n_clusters, palette)

  # Build the plots
  hc_list <- lapply(clusters, function(k) {
    hc <- tryCatch({
      cluster_profile_highchart(
        data = data,
        cluster = k,
        conditions = conditions,
        min_membership = min_membership,
        show_centroid = show_centroid,
        centroid_summary = centroid_summary,
        cluster_color = unname(cluster_colors[as.character(k)]),
        line_width = line_width,
        line_opacity = line_opacity,
        centroid_width = centroid_width,
        height = height
      )
    }, error = function(e) {
      warning(sprintf("Could not build the plot for Cluster %d: %s", k, e$message))
      return(NULL)
    })

    hc
  })

  names(hc_list) <- paste0("Cluster_", clusters)
  Filter(Negate(is.null), hc_list)
}


#' Centroid plot for all the clusters
#'
#' Builds a comparative plot with the centroids of all the clusters.
#'
#' @param data Pattern Profiler DataFrame (LONG format)
#' @param conditions Vector of condition names
#' @param clusters Clusters to include (NULL = all)
#' @param min_membership Membership filter applied before computing the centroids
#' @param centroid_summary Method: "mean" or "median"
#' @param palette Colour palette
#' @param line_width Line width (default: 2.5)
#' @param show_markers Show markers on the points (default: TRUE)
#' @param title Custom title
#' @param height Plot height
#'
#' @return highchart object
#'
#' @examples
#' if (requireNamespace("Mfuzz", quietly = TRUE) &&
#'     requireNamespace("e1071", quietly = TRUE)) {
#'   data(nadia_dia)
#'   res <- process_proteomics(nadia_dia, verbose = FALSE)
#'   pp <- pattern_profiler_analysis(res$se_proc, res$DEPs_results,
#'                                   assay_name = "Impseqrob_min",
#'                                   auto_select_c = FALSE, c = 3,
#'                                   verbose = FALSE)
#'
#'   # All the cluster centroids on one chart
#'   hc <- cluster_centroids_highchart(pp$long_output,
#'                                     conditions = c("A", "B", "D"),
#'                                     centroid_summary = "median")
#'   print(class(hc))
#' }
#'
#' @export
cluster_centroids_highchart <- function(data,
                                         conditions = NULL,
                                         clusters = NULL,
                                         min_membership = NULL,
                                         centroid_summary = c("mean", "median"),
                                         palette = NULL,
                                         line_width = 2.5,
                                         show_markers = TRUE,
                                         title = NULL,
                                         height = NULL) {

  centroid_summary <- match.arg(centroid_summary)

  # Detect the conditions
  if (is.null(conditions)) {
    conditions <- detect_condition_columns(data)
  }

  # Filter by membership if requested
  if (!is.null(min_membership)) {
    data <- data[data$Membership >= min_membership, ]
  }

  # Determine the clusters
  if (is.null(clusters)) {
    clusters <- sort(unique(data$Cluster))
  }

  n_clusters <- max(data$Cluster)

  # Configure the palette
  cluster_colors <- configure_cluster_palette(n_clusters, palette)

  # Compute the centroids per cluster
  agg_fun <- if (centroid_summary == "mean") mean else median

  centroids_list <- lapply(clusters, function(k) {
    cluster_data <- data[data$Cluster == k, conditions, drop = FALSE]
    n_proteins <- nrow(cluster_data)

    if (n_proteins == 0) return(NULL)

    centroid <- apply(cluster_data, 2, agg_fun, na.rm = TRUE)

    list(
      cluster = as.integer(k),
      centroid = unname(as.numeric(centroid)),
      n_proteins = as.integer(n_proteins),
      color = unname(cluster_colors[as.character(k)])
    )
  })

  centroids_list <- Filter(Negate(is.null), centroids_list)

  if (length(centroids_list) == 0) {
    warning("No data available to build the centroid plot")
    return(NULL)
  }

  # Title
  if (is.null(title)) {
    title <- sprintf("Cluster Centroids (%s)", centroid_summary)
  }

  # Build the highchart
  hc <- highchart() |>
    hc_chart(
      type = "line",
      backgroundColor = "#FFFFFF",
      style = list(fontFamily = "Inter, -apple-system, sans-serif"),
      height = height,
      zoomType = "xy"
    ) |>
    hc_title(
      text = title,
      style = list(
        fontSize = "18px",
        fontWeight = "600",
        color = "#1D3557"
      )
    ) |>
    hc_xAxis(
      categories = as.list(unname(conditions)),
      title = list(
        text = "Condition",
        style = list(
          fontSize = "13px",
          fontWeight = "bold",
          color = "#212529"
        )
      ),
      labels = list(
        style = list(
          fontSize = "12px",
          color = "#495057"
        )
      ),
      lineColor = "#DEE2E6",
      tickColor = "#DEE2E6"
    ) |>
    hc_yAxis(
      title = list(
        text = "z-score",
        style = list(
          fontSize = "13px",
          fontWeight = "bold",
          color = "#212529"
        )
      ),
      labels = list(
        style = list(
          fontSize = "11px",
          color = "#495057"
        )
      ),
      gridLineColor = "#F1F3F4",
      gridLineDashStyle = "Dot",
      plotLines = list(
        list(
          value = 0,
          color = "#ADB5BD",
          width = 1,
          dashStyle = "Dash",
          zIndex = 1
        )
      )
    ) |>
    hc_tooltip(
      useHTML = TRUE,
      backgroundColor = "rgba(255, 255, 255, 0.95)",
      borderColor = "#DEE2E6",
      borderRadius = 8,
      headerFormat = "",
      pointFormat = paste0(
        "<div style='padding: 4px;'>",
        "<b style='color: {series.color};'>{series.name}</b><br/>",
        "<span style='color: #6C757D;'>Condition:</span> <b>{point.condition}</b><br/>",
        "<span style='color: #6C757D;'>z-score:</span> <b>{point.y:.3f}</b>",
        "</div>"
      )
    ) |>
    hc_legend(
      enabled = TRUE,
      layout = "horizontal",
      align = "center",
      verticalAlign = "bottom",
      itemStyle = list(
        fontSize = "12px",
        fontWeight = "normal"
      )
    ) |>
    hc_exporting(
      enabled = TRUE,
      buttons = list(
        contextButton = list(
          menuItems = c("downloadPNG", "downloadSVG", "downloadPDF")
        )
      )
    )

  # Add the centroid series
  for (item in centroids_list) {
    # Build the points with clean values (without names)
    points <- lapply(seq_along(conditions), function(j) {
      list(
        x = as.integer(j - 1),
        y = round(item$centroid[j], 4),
        condition = unname(conditions[j])
      )
    })

    # Clean colour without names
    series_color <- unname(item$color)
    darker_color <- .darken_hex(series_color, 0.2)

    hc <- hc |> hc_add_series(
      name = sprintf("Cluster %d (n=%d)", item$cluster, item$n_proteins),
      type = "line",
      data = points,
      color = series_color,
      lineWidth = line_width,
      marker = list(
        enabled = show_markers,
        symbol = "circle",
        radius = 4,
        fillColor = series_color,
        lineColor = darker_color,
        lineWidth = 1
      ),
      connectNulls = TRUE
    )
  }

  hc
}


# =============================================================================
# SUMMARY HELPER FUNCTIONS
# =============================================================================

#' Summary of the Pattern Profiler data
#'
#' @param data Pattern Profiler DataFrame
#' @return List with summary statistics
#'
#' @examples
#' if (requireNamespace("Mfuzz", quietly = TRUE) &&
#'     requireNamespace("e1071", quietly = TRUE)) {
#'   data(nadia_dia)
#'   res <- process_proteomics(nadia_dia, verbose = FALSE)
#'   pp <- pattern_profiler_analysis(res$se_proc, res$DEPs_results,
#'                                   assay_name = "Impseqrob_min",
#'                                   auto_select_c = FALSE, c = 3,
#'                                   verbose = FALSE)
#'
#'   info <- summarize_pattern_profiler(pp$long_output)
#'   print(info$cluster_summary)
#' }
#'
#' @export
summarize_pattern_profiler <- function(data) {

  conditions <- detect_condition_columns(data)

  # Count per cluster
  cluster_summary <- data %>%
    group_by(Cluster) %>%
    summarise(
      n_entries = n(),
      n_unique_features = n_distinct(FeatureID),
      mean_membership = mean(Membership),
      min_membership = min(Membership),
      max_membership = max(Membership),
      .groups = "drop"
    )

  # Features in several clusters
  multi_cluster <- data %>%
    group_by(FeatureID) %>%
    summarise(n_clusters = n_distinct(Cluster), .groups = "drop") %>%
    filter(n_clusters > 1)

  list(
    n_total_entries = nrow(data),
    n_unique_features = n_distinct(data$FeatureID),
    n_clusters = n_distinct(data$Cluster),
    conditions = conditions,
    cluster_summary = cluster_summary,
    n_multi_cluster_features = nrow(multi_cluster),
    multi_cluster_features = multi_cluster
  )
}


# =============================================================================
# USAGE EXAMPLES
# =============================================================================

# --- Typical usage ---
#
# # Read the data
# data <- read_pattern_profiler_data(file.path(out_dir, "Pattern_Profiler_Input.parquet"))
#
# # Inspect the summary
# summary <- summarize_pattern_profiler(data)
# print(summary$cluster_summary)
#
# # Define the condition order
# conditions <- c("A", "B", "C", "D")
#
# # Plot for one specific cluster
# hc_c1 <- cluster_profile_highchart(data, cluster = 1, conditions = conditions)
# hc_c1
#
# # List of plots for all the clusters
# hc_profiles <- cluster_profile_highchart_list(data, conditions)
# hc_profiles[["Cluster_1"]]
# hc_profiles[["Cluster_2"]]
#
# # With a stricter membership filter
# hc_profiles <- cluster_profile_highchart_list(
#   data, conditions,
#   min_membership = 0.5
# )
#
# # Comparative centroid plot
# hc_centroids <- cluster_centroids_highchart(data, conditions)
# hc_centroids
#
# # With a custom palette
# hc_profiles <- cluster_profile_highchart_list(
#   data, conditions,
#   palette = "ggsci::nrc_npg"
# )
