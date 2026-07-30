# =============================================================================
# Interactive Highcharts boxplots for proteomics data
# =============================================================================


# -----------------------------------------------------------------------------
# Helper to compute the boxplot statistics
# -----------------------------------------------------------------------------

calc_boxplot_stats <- function(x, coef = 1.5) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NULL)

  q <- quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  iqr <- q[3] - q[1]

  lower_fence <- q[1] - coef * iqr
  upper_fence <- q[3] + coef * iqr

  # Values inside the whiskers
  whisker_low <- min(x[x >= lower_fence], na.rm = TRUE)
  whisker_high <- max(x[x <= upper_fence], na.rm = TRUE)

  # Outliers
  outliers <- x[x < lower_fence | x > upper_fence]

  list(
    low = whisker_low,
    q1 = unname(q[1]),
    median = unname(q[2]),
    q3 = unname(q[3]),
    high = whisker_high,
    outliers = outliers,
    mean = mean(x, na.rm = TRUE),
    n = length(x)
  )
}


# -----------------------------------------------------------------------------
# Main function: interactive Highcharts boxplot
# -----------------------------------------------------------------------------

#' Interactive Highcharts Boxplot for Proteomics
#'
#' @param data Data frame with columns: Column, Assay, Intensity, Condition
#' @param assays Vector of assays to include (NULL = all of them)
#' @param color_by Column used for colouring (default: "Condition")
#' @param group_order Order of the groups/conditions
#' @param palette Colour palette: "ggsci::palette", "brewer:Name", or a vector
#' @param title Chart title (optional).
#'   Use \code{\{assay\}} as a placeholder (e.g. "Boxplot: \{assay\}" -> "Boxplot: ImpSeqRob_Min")
#' @param subtitle Chart subtitle (optional).
#'   Use \code{\{assay\}} as a placeholder
#' @param show_outliers Show outliers lying outside the whiskers (default: TRUE)
#' @param outlier_jitter Amount of horizontal jitter applied to the outliers (default: 0.15)
#' @param outlier_size Radius of the outlier points (default: 3)
#' @param box_width Width of the boxplot boxes in pixels (default: 20)
#' @param horizontal Horizontal orientation (default: TRUE)
#' @param height Chart height in pixels
#'
#' @return List of highchart objects (one per assay)
#'
#' @examples
#' data(nadia_dia)
#' res <- process_proteomics(nadia_dia, verbose = FALSE)
#'
#' # One boxplot per assay; keep only the imputed one
#' hc_list <- boxplot_highchart_list(
#'   res$BoxPlot_Input,
#'   assays      = "Impseqrob_min",
#'   group_order = c("A", "B", "D"),
#'   title       = "Boxplot: {assay}"
#' )
#' names(hc_list)
#'
#' @export
boxplot_highchart_list <- function(
    data,
    assays = NULL,
    color_by = "Condition",
    group_order = NULL,
    palette = NULL,
    title = NULL,
    subtitle = NULL,
    show_outliers = TRUE,
    outlier_jitter = 0.15,
    outlier_size = 3,
    box_width = 20,
    horizontal = TRUE,
    height = NULL
) {

  # ---------------------------------------------------------------------------
  # 1) Input validation
  # ---------------------------------------------------------------------------
  required_cols <- c("Column", "Assay", "Intensity")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Check the colouring column
  if (!is.null(color_by) && !(color_by %in% names(data))) {
    warning("'", color_by, "' is not in the data frame. It will be ignored.")
    color_by <- NULL
  }

  # Filter assays
  if (!is.null(assays)) {
    data <- data[data$Assay %in% assays, , drop = FALSE]
  }

  # Drop NA and Inf
  data <- data[is.finite(data$Intensity), , drop = FALSE]

  available_assays <- unique(data$Assay)

  # ---------------------------------------------------------------------------
  # 2) Set up the group order
  # ---------------------------------------------------------------------------
  if (!is.null(color_by)) {
    if (!is.null(group_order)) {
      data[[color_by]] <- factor(data[[color_by]], levels = group_order)
    } else {
      data[[color_by]] <- factor(data[[color_by]], levels = unique(data[[color_by]]))
    }
    group_levels <- levels(data[[color_by]])
  } else {
    group_levels <- "All"
    data$`.group` <- "All"
    color_by <- ".group"
  }

  # ---------------------------------------------------------------------------
  # 3) Set up the colour palette
  # ---------------------------------------------------------------------------
  default_palette <- c(
    "#457B9D",
    "#E63946",
    "#2A9D8F",
    "#E9C46A",
    "#9B5DE5",
    "#F4A261",
    "#264653",
    "#00BBF9"
  )

  if (is.null(palette)) {
    pal <- default_palette
  } else if (is.character(palette) && length(palette) == 1 && grepl("::", palette)) {
    if (!requireNamespace("paletteer", quietly = TRUE)) {
      stop("To use paletteer, install it with: install.packages('paletteer')")
    }
    pal <- as.character(paletteer::paletteer_d(palette))
  } else if (is.character(palette) && length(palette) == 1 && startsWith(palette, "brewer:")) {
    nm <- sub("^brewer:", "", palette)
    maxc <- RColorBrewer::brewer.pal.info[nm, "maxcolors"]
    pal <- RColorBrewer::brewer.pal(maxc, nm)
  } else if (is.character(palette)) {
    pal <- palette
  } else {
    pal <- default_palette
  }

  if (length(pal) < length(group_levels)) {
    pal <- rep(pal, length.out = length(group_levels))
  }
  col_values <- stats::setNames(pal[seq_along(group_levels)], group_levels)

  # Assay labels
  assay_labels <- c(
    "log2" = "Log\u2082 Intensity",
    "LoessCyc" = "LOESS Cyclic Normalization",
    "raw" = "Raw Intensity",
    "vsn" = "VSN Normalized",
    "quantile" = "Quantile Normalized"
  )

  # ---------------------------------------------------------------------------
  # 4) Build one chart per assay
  # ---------------------------------------------------------------------------
  hc_list <- lapply(available_assays, function(current_assay) {

    dt <- data[data$Assay == current_assay, , drop = FALSE]
    if (nrow(dt) == 0) return(NULL)

    # Get the unique samples ordered by group
    samples_df <- dt %>%
      select(Column, all_of(color_by)) %>%
      distinct() %>%
      arrange(.data[[color_by]], Column)

    samples <- samples_df$Column
    sample_groups <- samples_df[[color_by]]

    # Compute the statistics per sample
    box_data <- lapply(seq_along(samples), function(i) {
      sample_name <- samples[i]
      group <- as.character(sample_groups[i])
      values <- dt$Intensity[dt$Column == sample_name]
      stats <- calc_boxplot_stats(values)

      if (is.null(stats)) return(NULL)

      list(
        sample = sample_name,
        group = group,
        color = unname(col_values[group]),
        stats = stats,
        index = i - 1  # 0-based index into the categories
      )
    })

    box_data <- Filter(Negate(is.null), box_data)

    # -------------------------------------------------------------------------
    # Build one boxplot series PER GROUP (so the legend stays interactive)
    # -------------------------------------------------------------------------
    boxplot_series_list <- lapply(group_levels, function(grp) {
      grp_data <- Filter(function(x) x$group == grp, box_data)
      if (length(grp_data) == 0) return(NULL)

      grp_color <- unname(col_values[grp])

      points <- lapply(grp_data, function(bd) {
        list(
          x = bd$index,
          low = bd$stats$low,
          q1 = bd$stats$q1,
          median = bd$stats$median,
          q3 = bd$stats$q3,
          high = bd$stats$high,
          name = bd$sample,
          n = bd$stats$n,
          mean = round(bd$stats$mean, 3)
        )
      })

      list(
        id = paste0("boxplot_", grp),
        name = grp,
        type = "boxplot",
        data = points,
        color = grp_color,
        fillColor = .hex_to_rgba(grp_color, 0.7),
        lineWidth = 1.5,
        whiskerLength = "50%",
        whiskerWidth = 2,
        medianColor = "#1D3557",
        medianWidth = 2,
        showInLegend = TRUE
      )
    })

    boxplot_series_list <- Filter(Negate(is.null), boxplot_series_list)

    # -------------------------------------------------------------------------
    # Build the outliers (per group, linked to the boxplots)
    # -------------------------------------------------------------------------
    outlier_series_list <- list()

    if (show_outliers) {
      outlier_series_list <- lapply(group_levels, function(grp) {
        grp_data <- Filter(function(x) x$group == grp, box_data)
        if (length(grp_data) == 0) return(NULL)

        grp_color <- unname(col_values[grp])

        points <- do.call(c, lapply(grp_data, function(bd) {
          if (length(bd$stats$outliers) == 0) return(NULL)
          lapply(bd$stats$outliers, function(o) {
            jittered_x <- bd$index + runif(1, -outlier_jitter, outlier_jitter)
            list(x = jittered_x, y = o, name = bd$sample)
          })
        }))

        if (length(points) == 0) return(NULL)

        list(
          name = grp,
          type = "scatter",
          data = points,
          color = grp_color,
          linkedTo = paste0("boxplot_", grp),  # Link to the boxplot of this group
          marker = list(
            symbol = "circle",
            radius = outlier_size,
            fillColor = grp_color,
            lineWidth = 1,
            lineColor = "#FFFFFF"
          ),
          tooltip = list(
            pointFormat = "<b>{point.name}</b><br/>Outlier: {point.y:.3f}"
          ),
          showInLegend = FALSE  # Not shown in the legend (linked to the boxplot)
        )
      })

      outlier_series_list <- Filter(Negate(is.null), outlier_series_list)
    }

    # Chart title (make sure it is a string)
    if (!is.null(title)) {
      chart_title <- gsub("{assay}", current_assay, as.character(title), fixed = TRUE)
    } else {
      chart_title <- if (current_assay %in% names(assay_labels)) {
        as.character(assay_labels[current_assay])
      } else {
        as.character(current_assay)
      }
    }

    # Build the highchart
    hc <- highchart() %>%
      hc_chart(
        type = "boxplot",
        inverted = horizontal,
        backgroundColor = "#FFFFFF",
        style = list(fontFamily = "Inter, -apple-system, sans-serif"),
        height = height
      ) %>%
      hc_title(
        text = chart_title,
        style = list(
          fontSize = "18px",
          fontWeight = "600",
          color = "#1D3557"
        )
      ) %>%
      hc_xAxis(
        categories = samples,
        title = list(
          text = "Samples",
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
        tickColor = "#DEE2E6",
        gridLineWidth = 0
      ) %>%
      hc_yAxis(
        title = list(
          text = "log<sub>2</sub> Intensity",
          useHTML = TRUE,
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
        gridLineDashStyle = "Dot"
      ) %>%
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
      ) %>%
      hc_plotOptions(
        boxplot = list(
          grouping = FALSE,  # IMPORTANT: keeps the series from being placed side by side
          groupPadding = 0.1,
          pointPadding = 0.05,
          borderRadius = 2,
          pointWidth = box_width
        ),
        scatter = list(
          jitter = list(x = 0, y = 0)  # The jitter is already applied manually
        )
      ) %>%
      hc_exporting(
        enabled = TRUE,
        buttons = list(
          contextButton = list(
            menuItems = c("downloadPNG", "downloadSVG", "downloadPDF", "separator", "downloadCSV")
          )
        )
      )

    # -------------------------------------------------------------------------
    # Add the boxplot series (one per group so the legend stays interactive)
    # -------------------------------------------------------------------------
    for (box_series in boxplot_series_list) {
      hc <- hc %>% hc_add_series(
        id = box_series$id,
        name = box_series$name,
        type = "boxplot",
        data = box_series$data,
        color = box_series$color,
        fillColor = box_series$fillColor,
        lineWidth = box_series$lineWidth,
        whiskerLength = box_series$whiskerLength,
        whiskerWidth = box_series$whiskerWidth,
        medianColor = box_series$medianColor,
        medianWidth = box_series$medianWidth,
        showInLegend = box_series$showInLegend,
        tooltip = list(
          headerFormat = "",
          pointFormat = paste0(
            "<div style='padding: 6px;'>",
            "<b style='font-size: 14px; color: #1D3557;'>{point.name}</b><br/>",
            "<span style='color: #6C757D;'>Max:</span> <b>{point.high:.3f}</b><br/>",
            "<span style='color: #6C757D;'>Q3:</span> <b>{point.q3:.3f}</b><br/>",
            "<span style='color: #6C757D;'>Median:</span> <b>{point.median:.3f}</b><br/>",
            "<span style='color: #6C757D;'>Q1:</span> <b>{point.q1:.3f}</b><br/>",
            "<span style='color: #6C757D;'>Min:</span> <b>{point.low:.3f}</b><br/>",
            "<span style='color: #ADB5BD; font-size: 11px;'>n = {point.n} | mean = {point.mean}</span>",
            "</div>"
          )
        )
      )
    }

    # -------------------------------------------------------------------------
    # Add the outlier series (linked to the boxplots)
    # -------------------------------------------------------------------------
    for (outlier_series in outlier_series_list) {
      hc <- hc %>% hc_add_series(
        name = outlier_series$name,
        type = "scatter",
        data = outlier_series$data,
        color = outlier_series$color,
        linkedTo = outlier_series$linkedTo,
        marker = outlier_series$marker,
        tooltip = outlier_series$tooltip,
        showInLegend = outlier_series$showInLegend
      )
    }

    # Subtitle
    if (!is.null(subtitle)) {
      sub_text <- gsub("{assay}", current_assay, as.character(subtitle), fixed = TRUE)
      hc <- hc %>% hc_subtitle(
        text = sub_text,
        style = list(
          fontSize = "13px",
          color = "#6C757D"
        )
      )
    }

    hc
  })

  names(hc_list) <- available_assays
  Filter(Negate(is.null), hc_list)
}


# =============================================================================
# USAGE EXAMPLES
# =============================================================================

# --- Basic example (boxplots + outliers by default) ---
# hc_boxplots <- boxplot_highchart_list(
#   data        = my_dataframe,
#   assays      = c("log2", "LoessCyc"),
#   color_by    = "Condition",
#   group_order = c("A", "B", "C", "D")
# )
# hc_boxplots[["log2"]]
# hc_boxplots[["LoessCyc"]]

# --- With a custom title ---
# hc_boxplots <- boxplot_highchart_list(
#   data     = se_proc,
#   assays   = "LoessCyc",
#   color_by = "Condition",
#   title    = "Intensity Distribution per Sample",
#   subtitle = "Cyclic LOESS normalization"
# )
# hc_boxplots[["LoessCyc"]]

# --- Vertical orientation ---
# hc_boxplots <- boxplot_highchart_list(
#   data       = se_proc,
#   color_by   = "Condition",
#   horizontal = FALSE
# )

# --- With a custom palette ---
# hc_boxplots <- boxplot_highchart_list(
#   data     = my_dataframe,
#   color_by = "Condition",
#   palette  = "ggsci::nrc_npg"
# )

# --- Without outliers (boxplots only) ---
# hc_boxplots <- boxplot_highchart_list(
#   data          = my_dataframe,
#   color_by      = "Condition",
#   show_outliers = FALSE
# )

# --- Customizing outlier size and jitter ---
# hc_boxplots <- boxplot_highchart_list(
#   data           = se_proc,
#   color_by       = "Condition",
#   outlier_jitter = 0,
#   outlier_size   = 3
# )

# --- Customizing the box width (useful depending on the number of samples) ---
# hc_boxplots <- boxplot_highchart_list(
#   data      = my_dataframe,
#   color_by  = "Condition",
#   box_width = 15  # Narrower when there are many samples
# )
