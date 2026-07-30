# =============================================================================
# Interactive PCA Plot with Highcharts for Proteomics Data
# =============================================================================



# -----------------------------------------------------------------------------
# Function to build PCA scores
# -----------------------------------------------------------------------------

#' Build a PCA scores data frame from data in long format
#'
#' @param pca_input Data frame in long format with columns:
#'   - SampleID: Sample identifier
#'   - FeatureID: Protein/feature identifier
#'   - Intensity: Intensity value (log2)
#'   - Condition: Experimental condition
#'   - Replicate: Replicate number (optional)
#'   - sig_any: Logical flagging significance in any comparison (for mode="any")
#'   - adjP_*: Adjusted p-value columns, one per comparison (for mode="specific")
#' @param mode Protein filtering mode: "all", "any", or "specific"
#' @param alpha Significance threshold for mode "specific" (default: 0.05)
#' @param comparison Name of the comparison for mode "specific" (e.g. "B-A")
#' @param subset_label Custom label for the subset (optional)
#' @param center Center the data before PCA (default: TRUE)
#' @param scale. Scale the data before PCA (default: TRUE)
#' @param filter_samples_to_comparison Restrict samples to the conditions
#'   involved in the specific comparison (default: FALSE)
#' @param cond_col Name of the condition column (default: "Condition")
#'
#' @return Data frame with columns: SampleID, PC1, PC2, PC1_Perc, PC2_Perc,
#'   Subset, Condition, Replicate
build_pca_scores <- function(pca_input,
                             mode = c("all", "any", "specific"),
                             alpha = 0.05,
                             comparison = NULL,
                             subset_label = NULL,
                             center = TRUE,
                             scale. = TRUE,
                             filter_samples_to_comparison = FALSE,
                             cond_col = "Condition") {

  mode <- match.arg(mode)

  # Default label depending on the mode
  if (is.null(subset_label)) {
    subset_label <- switch(mode,
      all = "All proteins",
      any = "DEPs (any comparison)",
      specific = paste0("DEPs (", comparison, ")")
    )
  }

  # Convert to data.frame to avoid issues with tibbles
  pca_input <- as.data.frame(pca_input)

  # Per-sample metadata
  md_cols <- intersect(c("SampleID", "Condition", "Replicate"), names(pca_input))
  md <- pca_input[!duplicated(pca_input$SampleID), md_cols, drop = FALSE]
  rownames(md) <- md$SampleID

  # Filter samples if requested (only for mode = "specific")
  keep_samples <- md$SampleID
  if (isTRUE(filter_samples_to_comparison) && mode == "specific" && !is.null(comparison)) {
    conds <- unique(trimws(strsplit(comparison, "[-|:]")[[1]]))
    if (!(cond_col %in% names(md))) {
      stop("Column '", cond_col, "' does not exist in pca_input.")
    }
    keep_samples <- rownames(md)[md[[cond_col]] %in% conds]
    if (length(keep_samples) < 2) {
      stop("Fewer than 2 samples left after filtering by comparison: ", comparison)
    }
  }

  # Get the feature IDs according to the mode
  ids <- .get_feature_ids(pca_input, mode = mode, alpha = alpha, comparison = comparison)
  if (length(ids) < 2) {
    stop("Subset '", mode, "' has too few proteins for PCA (minimum 2).")
  }

  # Filter the data
  dt <- pca_input[pca_input$SampleID %in% keep_samples & pca_input$FeatureID %in% ids,
                  c("SampleID", "FeatureID", "Intensity"), drop = FALSE]
  dt <- dt[is.finite(dt$Intensity) & !is.na(dt$Intensity), , drop = FALSE]

  # Pivot to a matrix (samples x features)
  Xt <- with(dt, tapply(Intensity, list(SampleID, FeatureID), mean))
  Xt <- as.matrix(Xt)

  # Drop features with zero variance
  v <- apply(Xt, 2, var, na.rm = TRUE)
  Xt <- Xt[, is.finite(v) & v > 0, drop = FALSE]

  # prcomp does not accept NAs: tapply leaves NA for the sample x feature
  # combinations that are absent. Keep only complete features (no NA in any sample).
  complete_feats <- colSums(is.na(Xt)) == 0
  n_dropped <- sum(!complete_feats)
  if (n_dropped > 0) {
    warning(sprintf("PCA '%s': %d features with NAs discarded before prcomp.",
                    subset_label, n_dropped))
    Xt <- Xt[, complete_feats, drop = FALSE]
  }

  if (ncol(Xt) < 2) {
    stop("Too few proteins with variance > 0 for PCA in '", subset_label, "'.")
  }

  # Run the PCA
  pc <- stats::prcomp(Xt, center = center, scale. = scale.)
  var_exp <- (pc$sdev^2) / sum(pc$sdev^2)
  scores <- pc$x[, 1:2, drop = FALSE]

  # Build the output data frame
  out <- data.frame(
    SampleID = rownames(scores),
    PC1 = as.numeric(scores[, 1]),
    PC2 = as.numeric(scores[, 2]),
    PC1_Perc = round(100 * var_exp[1], 2),
    PC2_Perc = round(100 * var_exp[2], 2),
    Subset = subset_label,
    stringsAsFactors = FALSE
  )

  # Add the metadata
  if ("Condition" %in% names(md)) {
    out$Condition <- md[out$SampleID, "Condition"]
  }
  if ("Replicate" %in% names(md)) {
    out$Replicate <- md[out$SampleID, "Replicate"]
  }

  out
}


# -----------------------------------------------------------------------------
# Function to compute the convex hull per group
# -----------------------------------------------------------------------------

#' Compute the convex hull (enclosing polygon) per group
#'
#' @param scores_df Data frame with columns PC1, PC2 and the grouping column
#' @param group_col Name of the grouping column (default: "Condition")
#'
#' @return List of data frames, each one with columns: group, x, y
compute_hulls <- function(scores_df, group_col = "Condition") {

  # Convert to data.frame to avoid issues with tibbles
  scores_df <- as.data.frame(scores_df)

  required <- c("PC1", "PC2", group_col)
  missing <- setdiff(required, names(scores_df))
  if (length(missing) > 0) {
    stop("Missing required columns for hulls: ", paste(missing, collapse = ", "))
  }

  split_list <- split(scores_df, scores_df[[group_col]], drop = TRUE)

  hulls <- lapply(names(split_list), function(g) {
    d <- split_list[[g]]
    d <- d[is.finite(d$PC1) & is.finite(d$PC2), , drop = FALSE]

    # At least 3 points are needed for a hull
    if (nrow(d) < 3) return(NULL)

    # Compute the hull (indices of the outline)
    h <- grDevices::chull(d$PC1, d$PC2)
    # Close the polygon by repeating the first point
    h <- c(h, h[1])

    data.frame(
      group = g,
      x = d$PC1[h],
      y = d$PC2[h],
      stringsAsFactors = FALSE
    )
  })

  Filter(Negate(is.null), hulls)
}


# -----------------------------------------------------------------------------
# Function to compute the confidence ellipse per group
# -----------------------------------------------------------------------------

#' Compute the confidence ellipse per group
#'
#' Computes the coordinates of a confidence ellipse based on the chi-squared
#' distribution, similar to FactoMineR::coord.ellipse and ggplot2::stat_ellipse.
#'
#' @param scores_df Data frame with columns PC1, PC2 and the grouping column
#' @param group_col Name of the grouping column (default: "Condition")
#' @param level Confidence level (default: 0.95)
#' @param npoints Number of points used to draw the ellipse (default: 100)
#'
#' @return List of data frames, each one with columns: group, x, y
compute_confidence_ellipse <- function(scores_df,
                                       group_col = "Condition",
                                       level = 0.95,
                                       npoints = 100) {

  # Convert to data.frame
  scores_df <- as.data.frame(scores_df)

  required <- c("PC1", "PC2", group_col)
  missing <- setdiff(required, names(scores_df))
  if (length(missing) > 0) {
    stop("Missing required columns for ellipse: ", paste(missing, collapse = ", "))
  }

  split_list <- split(scores_df, scores_df[[group_col]], drop = TRUE)

  ellipses <- lapply(names(split_list), function(g) {
    d <- split_list[[g]]
    d <- d[is.finite(d$PC1) & is.finite(d$PC2), , drop = FALSE]

    # At least 3 points are needed for an ellipse
    if (nrow(d) < 3) return(NULL)

    # Coordinates
    x <- d$PC1
    y <- d$PC2

    # Center (mean)
    center_x <- mean(x)
    center_y <- mean(y)

    # Covariance matrix
    cov_mat <- cov(cbind(x, y))

    # Radius based on the chi-squared distribution with 2 degrees of freedom
    # Similar to FactoMineR: sqrt(qchisq(level, df = 2))
    radius <- sqrt(stats::qchisq(level, df = 2))

    # Eigen decomposition to obtain the axes of the ellipse
    eigen_decomp <- eigen(cov_mat)
    eigenvalues <- eigen_decomp$values
    eigenvectors <- eigen_decomp$vectors

    # Check that the eigenvalues are positive
    if (any(eigenvalues <= 0)) return(NULL)

    # Angles used to parameterize the ellipse
    theta <- seq(0, 2 * pi, length.out = npoints + 1)

    # Semi-axes of the ellipse
    a <- radius * sqrt(eigenvalues[1])
    b <- radius * sqrt(eigenvalues[2])

    # Rotation angle
    angle <- atan2(eigenvectors[2, 1], eigenvectors[1, 1])

    # Coordinates of the ellipse (parametric form)
    ellipse_x <- center_x + a * cos(theta) * cos(angle) - b * sin(theta) * sin(angle)
    ellipse_y <- center_y + a * cos(theta) * sin(angle) + b * sin(theta) * cos(angle)

    data.frame(
      group = g,
      x = ellipse_x,
      y = ellipse_y,
      stringsAsFactors = FALSE
    )
  })

  Filter(Negate(is.null), ellipses)
}


# -----------------------------------------------------------------------------
# Main function: PCA Highchart
# -----------------------------------------------------------------------------

#' Interactive PCA Plot with Highcharts
#'
#' Builds a PCA scatter plot with the option of showing per-group ellipses or
#' convex hulls, similar to factoextra::fviz_pca_ind.
#'
#' @param scores_df Data frame produced by build_pca_scores() with columns:
#'   PC1, PC2, PC1_Perc, PC2_Perc, Subset, Condition, Replicate, SampleID
#' @param color_by Column used to colour the points (default: "Condition")
#' @param group_order Vector with the order of the groups/conditions (optional)
#' @param palette Named vector of colours, or NULL for an automatic palette
#' @param title Chart title (optional, defaults to Subset)
#' @param addEllipses Show ellipses/hulls around the groups (default: TRUE)
#' @param ellipse_type Ellipse type: "convex" for a convex hull or "confidence"
#'   for a confidence ellipse based on the normal distribution (default: "convex")
#' @param ellipse_level Confidence level for ellipse_type = "confidence"
#'   (default: 0.95). Typical values: 0.95, 0.90, 0.68
#' @param ellipse_fill_opacity Fill opacity (0-1, default: 0.12)
#' @param ellipse_line_width Outline line width (default: 1)
#' @param ellipse_npoints Number of points used to draw the confidence ellipse
#'   (default: 100). Only applies to ellipse_type = "confidence"
#' @param point_size Point radius (default: 5)
#' @param show_labels Show the point labels (SampleID) (default: FALSE)
#' @param label_size Label font size in px (default: 10)
#'
#' @return A highchart object
#' @export
pca_highchart <- function(scores_df,
                          color_by = "Condition",
                          group_order = NULL,
                          palette = NULL,
                          title = NULL,
                          addEllipses = TRUE,
                          ellipse_type = c("convex", "confidence"),
                          ellipse_level = 0.95,
                          ellipse_fill_opacity = 0.12,
                          ellipse_line_width = 1,
                          ellipse_npoints = 100,
                          point_size = 5,
                          show_labels = FALSE,
                          label_size = 10) {

  ellipse_type <- match.arg(ellipse_type)

  # ---------------------------------------------------------------------------
  # 1) Input validation and conversion to data.frame
  # ---------------------------------------------------------------------------

  # Convert to data.frame to avoid issues with tibbles
  scores_df <- as.data.frame(scores_df)

  required <- c("PC1", "PC2", "SampleID")
  missing <- setdiff(required, names(scores_df))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "))
  }

  if (!is.null(color_by) && !(color_by %in% names(scores_df))) {
    stop("Column '", color_by, "' does not exist in scores_df.")
  }

  # ---------------------------------------------------------------------------
  # 2) Set up the axis labels with the % of variance
  # ---------------------------------------------------------------------------
  pc1p <- unique(scores_df$PC1_Perc)
  pc2p <- unique(scores_df$PC2_Perc)
  x_lab <- if (length(pc1p) == 1) paste0("PC1 (", pc1p, "%)") else "PC1"
  y_lab <- if (length(pc2p) == 1) paste0("PC2 (", pc2p, "%)") else "PC2"

  # ---------------------------------------------------------------------------
  # 3) Set up the group levels
  # ---------------------------------------------------------------------------
  if (!is.null(group_order)) {
    scores_df[[color_by]] <- factor(scores_df[[color_by]], levels = group_order)
  } else {
    scores_df[[color_by]] <- factor(scores_df[[color_by]])
  }
  lvls <- levels(scores_df[[color_by]])

  # ---------------------------------------------------------------------------
  # 4) Set up the colour palette
  # ---------------------------------------------------------------------------
  default_palette <- c(
    "#457B9D", "#E63946", "#2A9D8F", "#E9C46A",
    "#9B5DE5", "#F4A261", "#264653", "#00BBF9"
  )

  if (is.null(palette)) {
    # Default palette
    pal <- grDevices::hcl.colors(length(lvls), "Dark 3")
    palette <- stats::setNames(pal, lvls)
  } else if (is.character(palette) && length(palette) == 1 && grepl("::", palette)) {
    # paletteer palette (format "ggsci::category10_d3")
    if (!requireNamespace("paletteer", quietly = TRUE)) {
      stop("To use paletteer, install it with: install.packages('paletteer')")
    }
    pal <- tryCatch({
      raw_pal <- as.character(paletteer::paletteer_d(palette))
      # Normalize the colours (drop the alpha channel if present)
      vapply(raw_pal, .normalize_hex, character(1), USE.NAMES = FALSE)
    }, error = function(e) {
      stop("Error loading palette '", palette, "': ", e$message)
    })
    if (length(pal) < length(lvls)) {
      pal <- rep(pal, length.out = length(lvls))
    }
    palette <- stats::setNames(pal[seq_along(lvls)], lvls)
  } else if (is.character(palette) && length(palette) == 1 && startsWith(palette, "brewer:")) {
    # RColorBrewer palette (format "brewer:Set1")
    if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
      stop("To use RColorBrewer, install it with: install.packages('RColorBrewer')")
    }
    nm <- sub("^brewer:", "", palette)
    pal <- tryCatch({
      maxc <- RColorBrewer::brewer.pal.info[nm, "maxcolors"]
      raw_pal <- RColorBrewer::brewer.pal(maxc, nm)
      # Normalize the colours for consistency
      vapply(raw_pal, .normalize_hex, character(1), USE.NAMES = FALSE)
    }, error = function(e) {
      stop("Error loading brewer palette '", nm, "': ", e$message)
    })
    if (length(pal) < length(lvls)) {
      pal <- rep(pal, length.out = length(lvls))
    }
    palette <- stats::setNames(pal[seq_along(lvls)], lvls)
  } else if (is.character(palette) && length(palette) == 1) {
    # Single colour string - repeat it for every level
    palette <- stats::setNames(rep(palette, length(lvls)), lvls)
  } else if (is.character(palette) && length(palette) > 1 && is.null(names(palette))) {
    # Unnamed vector of colours
    if (length(palette) < length(lvls)) {
      palette <- rep(palette, length.out = length(lvls))
    }
    palette <- stats::setNames(palette[seq_along(lvls)], lvls)
  } else if (is.character(palette) && !is.null(names(palette))) {
    # Named vector of colours
    miss <- setdiff(lvls, names(palette))
    if (length(miss) > 0) {
      stop("Missing colours for levels: ", paste(miss, collapse = ", "))
    }
    palette <- palette[lvls]
  }

  # ---------------------------------------------------------------------------
  # 5) Chart title
  # ---------------------------------------------------------------------------
  # The subset is used as the title not only when `title` is NULL, but also when
  # it comes in empty or NA. This used to be handled by a local variant of
  # `%||%`; now that the operator is the canonical one (NULL only), the check is
  # made explicit.
  chart_title <- if (!is.null(title) && length(title) && !is.na(title[1])) {
    title
  } else {
    unique(scores_df$Subset)[1]
  }

  # ---------------------------------------------------------------------------
  # 6) Build the base highchart
  # ---------------------------------------------------------------------------
  hc <- highchart() |>
    hc_chart(
      type = "scatter",
      zoomType = "xy",
      backgroundColor = "#FFFFFF",
      style = list(fontFamily = "Inter, -apple-system, sans-serif")
    ) |>
    hc_title(
      text = chart_title,
      style = list(
        fontSize = "18px",
        fontWeight = "600",
        color = "#1D3557"
      )
    ) |>
    hc_xAxis(
      title = list(
        text = x_lab,
        style = list(
          fontSize = "13px",
          fontWeight = "bold",
          color = "#212529"
        )
      ),
      gridLineWidth = 1,
      gridLineColor = "#F1F3F4",
      gridLineDashStyle = "Dot",
      lineColor = "#DEE2E6",
      tickColor = "#DEE2E6",
      plotLines = list(
        list(value = 0, color = "#ADB5BD", width = 1, dashStyle = "Dash", zIndex = 1)
      )
    ) |>
    hc_yAxis(
      title = list(
        text = y_lab,
        style = list(
          fontSize = "13px",
          fontWeight = "bold",
          color = "#212529"
        )
      ),
      gridLineWidth = 1,
      gridLineColor = "#F1F3F4",
      gridLineDashStyle = "Dot",
      lineColor = "#DEE2E6",
      lineWidth = 1,
      plotLines = list(
        list(value = 0, color = "#ADB5BD", width = 1, dashStyle = "Dash", zIndex = 1)
      )
    ) |>
    hc_legend(
      enabled = TRUE,
      layout = "horizontal",
      align = "center",
      verticalAlign = "bottom",
      itemStyle = list(
        fontSize = "12px",
        fontWeight = "normal",
        color = "#495057"
      ),
      itemHoverStyle = list(color = "#1D3557")
    ) |>
    hc_plotOptions(
      scatter = list(
        marker = list(
          radius = point_size,
          symbol = "circle",
          states = list(
            hover = list(
              radiusPlus = 2,
              lineWidthPlus = 1
            )
          )
        )
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

  # ---------------------------------------------------------------------------
  # 7) Add one scatter series per group (with an id to link the ellipses)
  # ---------------------------------------------------------------------------
  split_list <- split(scores_df, scores_df[[color_by]], drop = TRUE)

  for (g in lvls) {
    if (!(g %in% names(split_list))) next

    d <- as.data.frame(split_list[[g]])
    group_id <- paste0("scatter_", gsub("[^a-zA-Z0-9]", "_", as.character(g)))

    # Build the list of points with explicit scalar values
    pts <- vector("list", nrow(d))
    for (i in seq_len(nrow(d))) {
      pts[[i]] <- list(
        x = as.numeric(d[i, "PC1"]),
        y = as.numeric(d[i, "PC2"]),
        SampleID = as.character(d[i, "SampleID"]),
        Condition = if ("Condition" %in% names(d)) as.character(d[i, "Condition"]) else "",
        Replicate = if ("Replicate" %in% names(d)) as.integer(d[i, "Replicate"]) else NA_integer_,
        Subset = as.character(d[i, "Subset"])
      )
    }

    # Get the group colour and a darkened version for the labels
    group_color <- unname(palette[as.character(g)])
    label_color <- .darken_hex(group_color, factor = 0.3)

    # Set up dataLabels if show_labels = TRUE
    data_labels_config <- if (isTRUE(show_labels)) {
      list(
        enabled = TRUE,
        format = "{point.SampleID}",
        style = list(
          fontSize = paste0(label_size, "px"),
          fontWeight = "bold",
          color = label_color,
          textOutline = "none"
        ),
        y = -10,
        allowOverlap = FALSE
      )
    } else {
      list(enabled = FALSE)
    }

    hc <- hc |>
      hc_add_series(
        data = pts,
        type = "scatter",
        id = group_id,
        name = as.character(g),
        color = group_color,
        zIndex = 5,
        marker = list(
          lineColor = label_color,
          lineWidth = 1
        ),
        dataLabels = data_labels_config,
        tooltip = list(
          headerFormat = "",
          pointFormat = paste0(
            "<div style='padding: 6px;'>",
            "<b style='font-size: 14px; color: #1D3557;'>{point.SampleID}</b><br/>",
            "<span style='color: #6C757D;'>Condition:</span> <b>{point.Condition}</b><br/>",
            "<span style='color: #6C757D;'>Replicate:</span> <b>{point.Replicate}</b><br/>",
            "<span style='color: #6C757D;'>PC1:</span> <b>{point.x:.3f}</b><br/>",
            "<span style='color: #6C757D;'>PC2:</span> <b>{point.y:.3f}</b>",
            "</div>"
          )
        )
      )
  }

  # ---------------------------------------------------------------------------
  # 8) Add the ellipses/hulls linked to the scatter series
  # ---------------------------------------------------------------------------
  if (isTRUE(addEllipses)) {

    # Compute the coordinates according to the ellipse type
    if (ellipse_type == "convex") {
      ellipse_coords <- compute_hulls(scores_df, group_col = color_by)
    } else {
      # ellipse_type == "confidence"
      ellipse_coords <- compute_confidence_ellipse(
        scores_df,
        group_col = color_by,
        level = ellipse_level,
        npoints = ellipse_npoints
      )
    }

    for (poly in ellipse_coords) {
      g <- unique(poly$group)
      group_id <- paste0("scatter_", gsub("[^a-zA-Z0-9]", "_", as.character(g)))
      base_color <- unname(palette[as.character(g)])
      rgba_color <- .hex_to_rgba(base_color, ellipse_fill_opacity)

      pts <- lapply(seq_len(nrow(poly)), function(k) {
        list(x = poly$x[k], y = poly$y[k])
      })

      hc <- hc |>
        hc_add_series(
          data = pts,
          type = "polygon",
          name = as.character(g),
          linkedTo = group_id,
          color = rgba_color,
          fillColor = rgba_color,
          fillOpacity = ellipse_fill_opacity,
          lineWidth = ellipse_line_width,
          lineColor = base_color,
          marker = list(enabled = FALSE),
          enableMouseTracking = FALSE,
          showInLegend = FALSE,
          zIndex = 0
        )
    }
  }

  hc
}


# -----------------------------------------------------------------------------
# Wrapper function: build a list of PCA plots
# -----------------------------------------------------------------------------

#' Build a List of PCA Plots for Multiple Subsets
#'
#' Automatically builds PCA plots for "all", "any" and/or specific comparisons.
#' Similar to factoextra::fviz_pca_ind with ellipse options.
#'
#' @param pca_input Data frame in long format (see build_pca_scores for the structure)
#' @param modes Vector of modes to generate: "all", "any", and/or comparison names
#'   (default: c("all", "any"))
#' @param alpha Significance threshold for DEPs (default: 0.05)
#' @param color_by Column used to colour the points (default: "Condition")
#' @param group_order Order of the groups/conditions (optional)
#' @param palette Named vector of colours, or NULL for an automatic palette
#' @param addEllipses Show ellipses/hulls around the groups (default: TRUE)
#' @param ellipse_type Ellipse type: "convex" or "confidence" (default: "convex")
#' @param ellipse_level Confidence level for ellipse_type = "confidence"
#'   (default: 0.95)
#' @param ellipse_fill_opacity Ellipse fill opacity (0-1, default: 0.12)
#' @param ellipse_line_width Outline line width (default: 1)
#' @param ellipse_npoints Number of points for the confidence ellipse (default: 100)
#' @param point_size Point radius (default: 5)
#' @param show_labels Show the point labels (SampleID) (default: FALSE)
#' @param label_size Label font size in px (default: 10)
#' @param center Center the data before PCA (default: TRUE)
#' @param scale. Scale the data before PCA (default: TRUE)
#' @param filter_samples_to_comparison For specific comparisons, restrict the
#'   samples to the conditions involved (default: FALSE)
#'
#' @return Named list of highchart objects
#'
#' @examples
#' \dontrun{
#' # Load the data
#' pca_input <- arrow::read_parquet("PCA_Input.parquet")
#'
#' # PCA with convex hull (default)
#' hc_pcas <- pca_highchart_list(
#'   pca_input   = pca_input,
#'   modes       = c("all", "any"),
#'   group_order = c("A", "B", "C", "D")
#' )
#'
#' # PCA with a 95% confidence ellipse
#' hc_pcas <- pca_highchart_list(
#'   pca_input     = pca_input,
#'   modes         = c("all"),
#'   group_order   = c("A", "B", "C", "D"),
#'   ellipse_type  = "confidence",
#'   ellipse_level = 0.95
#' )
#'
#' # PCA with visible labels
#' hc_pcas <- pca_highchart_list(
#'   pca_input   = pca_input,
#'   modes       = c("all"),
#'   show_labels = TRUE,
#'   label_size  = 9
#' )
#'
#' # Display
#' hc_pcas[["all"]]
#' }
#' @export
pca_highchart_list <- function(pca_input,
                               modes = c("all", "any"),
                               alpha = 0.05,
                               color_by = "Condition",
                               group_order = NULL,
                               palette = NULL,
                               addEllipses = TRUE,
                               ellipse_type = c("convex", "confidence"),
                               ellipse_level = 0.95,
                               ellipse_fill_opacity = 0.12,
                               ellipse_line_width = 1,
                               ellipse_npoints = 100,
                               point_size = 5,
                               show_labels = FALSE,
                               label_size = 10,
                               center = TRUE,
                               scale. = TRUE,
                               filter_samples_to_comparison = FALSE) {

  ellipse_type <- match.arg(ellipse_type)

  # ---------------------------------------------------------------------------
  # 1) Input validation
  # ---------------------------------------------------------------------------
  required_cols <- c("SampleID", "FeatureID", "Intensity")
  missing_cols <- setdiff(required_cols, names(pca_input))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (!is.null(color_by) && !(color_by %in% names(pca_input))) {
    warning("'", color_by, "' is not in the data frame. It will be ignored.")
    color_by <- NULL
  }

  # ---------------------------------------------------------------------------
  # 2) Detect the available comparisons (adjP_* columns)
  # ---------------------------------------------------------------------------
  adjp_cols <- grep("^adjP_", names(pca_input), value = TRUE)
  available_comparisons <- sub("^adjP_", "", adjp_cols)

  # ---------------------------------------------------------------------------
  # 3) Build the plots for each mode
  # ---------------------------------------------------------------------------
  hc_list <- list()

  for (m in modes) {

    # Determine the internal mode and the comparison (where applicable)
    if (m == "all") {
      internal_mode <- "all"
      comparison <- NULL
      subset_label <- "All proteins"
      plot_title <- "PCA (All proteins)"
    } else if (m == "any") {
      internal_mode <- "any"
      comparison <- NULL
      subset_label <- "DEPs (any comparison)"
      plot_title <- "PCA (DEPs in any comparison)"
    } else {
      # It is a specific comparison
      if (!(m %in% available_comparisons)) {
        warning("Comparison '", m, "' not found. It will be skipped.")
        next
      }
      internal_mode <- "specific"
      comparison <- m
      subset_label <- paste0("DEPs (", m, ")")
      plot_title <- paste0("PCA (DEPs ", m, ")")
    }

    # Try to build the scores (it may fail if there are too few proteins)
    scores_df <- tryCatch({
      build_pca_scores(
        pca_input = pca_input,
        mode = internal_mode,
        alpha = alpha,
        comparison = comparison,
        subset_label = subset_label,
        center = center,
        scale. = scale.,
        filter_samples_to_comparison = filter_samples_to_comparison
      )
    }, error = function(e) {
      warning("Error building the PCA for '", m, "': ", e$message)
      return(NULL)
    })

    if (is.null(scores_df)) next

    # Build the plot
    hc <- pca_highchart(
      scores_df = scores_df,
      color_by = color_by,
      group_order = group_order,
      palette = palette,
      title = plot_title,
      addEllipses = addEllipses,
      ellipse_type = ellipse_type,
      ellipse_level = ellipse_level,
      ellipse_fill_opacity = ellipse_fill_opacity,
      ellipse_line_width = ellipse_line_width,
      ellipse_npoints = ellipse_npoints,
      point_size = point_size,
      show_labels = show_labels,
      label_size = label_size
    )

    hc_list[[m]] <- hc
  }

  hc_list
}


# =============================================================================
# USAGE EXAMPLES
# =============================================================================

# --- Load the data ---
# pca_input <- arrow::read_parquet("PCA_Input.parquet")
# pca_input <- readr::read_tsv("PCA_Input.tsv")

# --- Basic example: PCA with convex hull (default) ---
# sc_all <- build_pca_scores(pca_input, mode = "all")
# p_all <- pca_highchart(sc_all, color_by = "Condition", group_order = c("A","B","C","D"))
# p_all

# --- PCA with a 95% confidence ellipse ---
# sc_all <- build_pca_scores(pca_input, mode = "all")
# p_all <- pca_highchart(
#   sc_all,
#   color_by = "Condition",
#   group_order = c("A","B","C","D"),
#   ellipse_type = "confidence",
#   ellipse_level = 0.95
# )
# p_all

# --- PCA with a 68% confidence ellipse (1 standard deviation) ---
# p_68 <- pca_highchart(
#   sc_all,
#   color_by = "Condition",
#   ellipse_type = "confidence",
#   ellipse_level = 0.68,
#   ellipse_fill_opacity = 0.2
# )
# p_68

# --- PCA without ellipses (points only) ---
# p_noellipse <- pca_highchart(
#   sc_all,
#   color_by = "Condition",
#   addEllipses = FALSE,
#   point_size = 6
# )
# p_noellipse

# --- Build several PCA plots with convex hulls ---
# hc_pcas <- pca_highchart_list(
#   pca_input   = pca_input,
#   modes       = c("all", "any", "B-A", "C-A", "D-A"),
#   alpha       = 0.05,
#   group_order = c("A", "B", "C", "D"),
#   ellipse_type = "convex"
# )
# hc_pcas[["all"]]
# hc_pcas[["any"]]
# hc_pcas[["B-A"]]

# --- Build PCA plots with confidence ellipses ---
# hc_pcas <- pca_highchart_list(
#   pca_input     = pca_input,
#   modes         = c("all", "any"),
#   group_order   = c("A", "B", "C", "D"),
#   ellipse_type  = "confidence",
#   ellipse_level = 0.95,
#   ellipse_fill_opacity = 0.15
# )
# hc_pcas[["all"]]

# --- With a custom palette ---
# my_palette <- c(A = "#457B9D", B = "#E63946", C = "#2A9D8F", D = "#E9C46A")
# hc_pcas <- pca_highchart_list(
#   pca_input   = pca_input,
#   modes       = c("all", "any"),
#   group_order = c("A", "B", "C", "D"),
#   palette     = my_palette
# )

# --- Comparing different confidence levels ---
# sc_all <- build_pca_scores(pca_input, mode = "all")
#
# # 68% (1 SD)
# p_68 <- pca_highchart(sc_all, ellipse_type = "confidence", ellipse_level = 0.68,
#                       title = "PCA - 68% CI")
# # 95% (approx. 2 SD)
# p_95 <- pca_highchart(sc_all, ellipse_type = "confidence", ellipse_level = 0.95,
#                       title = "PCA - 95% CI")
# # 99%
# p_99 <- pca_highchart(sc_all, ellipse_type = "confidence", ellipse_level = 0.99,
#                       title = "PCA - 99% CI")

# --- PCA with visible labels (useful for exporting) ---
# p_labels <- pca_highchart(
#   sc_all,
#   color_by = "Condition",
#   group_order = c("A","B","C","D"),
#   show_labels = TRUE,
#   label_size = 9
# )
# p_labels

# --- PCA with labels and a confidence ellipse ---
# p_labels_ellipse <- pca_highchart(
#   sc_all,
#   color_by = "Condition",
#   group_order = c("A","B","C","D"),
#   ellipse_type = "confidence",
#   ellipse_level = 0.95,
#   show_labels = TRUE,
#   label_size = 10
# )
# p_labels_ellipse

# --- Using pca_highchart_list with labels ---
# hc_pcas <- pca_highchart_list(
#   pca_input   = pca_input,
#   modes       = c("all", "any"),
#   group_order = c("A", "B", "C", "D"),
#   show_labels = TRUE,
#   label_size  = 9
# )
# hc_pcas[["all"]]

# --- With a paletteer palette ---
# hc_pcas <- pca_highchart_list(
#   pca_input   = pca_input,
#   modes       = c("all"),
#   group_order = c("A", "B", "C", "D"),
#   palette     = "ggsci::category10_d3"
# )
# hc_pcas[["all"]]

# --- With an RColorBrewer palette ---
# hc_pcas <- pca_highchart_list(
#   pca_input   = pca_input,
#   modes       = c("all"),
#   group_order = c("A", "B", "C", "D"),
#   palette     = "brewer:Set1"
# )
# hc_pcas[["all"]]
