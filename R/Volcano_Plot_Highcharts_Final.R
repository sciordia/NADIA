
# Load the libraries


#' Interactive Volcano Plot with Highcharter
#'
#' @param de_res Data frame with differential expression results
#' @param ain Vector of assays to filter on (optional)
#' @param comparisons Vector of comparisons to include (optional)
#' @param lfc_thr log2 fold-change threshold (default: 0)
#' @param alpha Significance threshold (default: 0.05)
#' @param p_col Column of p-values to use
#' @param point_size Point size (default: 4)
#' @param colors List of colours for "up", "down", "ns" (ignored when palette is used)
#' @param palette Name of a paletteer palette (e.g. "ggsci::default_jco"). Uses 3 colours: up, down, ns
#' @param show_top_genes Number of top genes to label by significance (default: 0)
#' @param highlight_genes Vector of gene names to highlight manually (default: NULL)
#' @param title Custom chart title (default: NULL, uses "Comparison (Assay)").
#'   Use \code{\{comparison\}} as a placeholder (e.g. "Volcano Plot: \{comparison\}" -> "Volcano Plot: B-A")
#'
#' @return List of highchart objects
#' @export
volcano_highchart_list <- function(
    de_res,
    ain = NULL,
    comparisons = NULL,
    lfc_thr = 0,
    alpha = 0.05,
    p_col = "adj.P.Val",
    point_size = 4,
    colors = NULL,
    palette = NULL,
    show_top_genes = 0,
    highlight_genes = NULL,
    title = NULL
) {
  
  # --- Colour handling ---
  if (!is.null(palette)) {
    if (!requireNamespace("paletteer", quietly = TRUE)) {
      stop("The 'paletteer' package is required to use palettes. Install it with install.packages('paletteer')")
    }
    pal_colors <- as.character(paletteer::paletteer_d(palette, n = 3))
    colors <- list(
      up   = pal_colors[1],
      down = pal_colors[2],
      ns   = pal_colors[3]
    )
  } else {
    default_colors <- list(
      up   = "#457B9D",
      down = "#E63946",
      ns   = "#ADB5BD"
    )
    colors <- modifyList(default_colors, colors %||% list())
  }
  
  # --- Validation of the required columns ---
  validate_columns <- function(df) {
    required <- c("logFC", "Comparison")
    missing <- setdiff(required, names(df))
    if (length(missing) > 0) {
      stop("Missing required columns: ", paste(missing, collapse = ", "))
    }
  }
  
  validate_columns(de_res)
  
  # --- Filtering by Assay and Comparison ---
  if (!is.null(ain) && "Assay" %in% names(de_res)) {
    de_res <- de_res[de_res$Assay %in% ain, , drop = FALSE]
  }
  
  if (!is.null(comparisons) && "Comparison" %in% names(de_res)) {
    de_res <- de_res[de_res$Comparison %in% comparisons, , drop = FALSE]
  } else {
    comparisons <- unique(de_res$Comparison)
  }
  
  if (nrow(de_res) == 0) {
    stop("The data frame is empty after applying the filters.")
  }
  
  p_col <- detect_pvalue_col(de_res, p_col)
  
  # --- Helper that assigns the change category ---
  assign_change <- function(dt) {
    pvals <- suppressWarnings(as.numeric(dt[[p_col]]))
    lfc <- as.numeric(dt$logFC)

    # Use a standalone vector to avoid problems with factors/tibbles
    change_vec <- rep("Not Significant", nrow(dt))
    significant <- !is.na(pvals) & pvals < alpha & abs(lfc) >= lfc_thr
    change_vec[significant & lfc > 0] <- "Up"
    change_vec[significant & lfc < 0] <- "Down"

    dt$Change <- factor(change_vec, levels = c("Not Significant", "Up", "Down"))
    dt
  }
  
  color_map <- c(
    "Not Significant" = colors$ns,
    "Up"              = colors$up,
    "Down"            = colors$down
  )
  
  # --- Build the list of plots ---
  hclist <- lapply(comparisons, function(comp) {
    
    dt <- de_res[de_res$Comparison == comp, , drop = FALSE]
    if (nrow(dt) == 0) return(NULL)

    # Drop rows with an NA p-value: replacing -log10(NA) would push them to the
    # top of the volcano, which would be misleading. p = 0 is kept (Inf is
    # mapped to the real ceiling below).
    dt <- dt[!is.na(suppressWarnings(as.numeric(dt[[p_col]]))), , drop = FALSE]
    if (nrow(dt) == 0) return(NULL)

    dt <- assign_change(dt)
    dt$minusLog10P <- -log10(as.numeric(dt[[p_col]]))
    dt$pval_fmt <- sprintf("%.3g", as.numeric(dt[[p_col]]))

    finite_ml <- dt$minusLog10P[is.finite(dt$minusLog10P)]
    max_finite <- if (length(finite_ml) > 0) max(finite_ml) else 1
    dt$minusLog10P[!is.finite(dt$minusLog10P)] <- max_finite * 1.1
    
    # --- Identify the genes to highlight ---
    dt$Highlight <- FALSE
    
    if (show_top_genes > 0) {
      sig_genes <- dt[dt$Change != "Not Significant", ]
      if (nrow(sig_genes) > 0) {
        sig_genes <- sig_genes[order(as.numeric(sig_genes[[p_col]])), ]
        top_ids <- head(sig_genes$Gene.Names, show_top_genes)
        dt$Highlight[dt$Gene.Names %in% top_ids] <- TRUE
      }
    }
    
    if (!is.null(highlight_genes) && length(highlight_genes) > 0) {
      dt$Highlight[dt$Gene.Names %in% highlight_genes] <- TRUE
    }
    
    # --- Split the data: regular vs highlighted points ---
    dt_normal <- dt[!dt$Highlight, ]
    dt_highlighted <- dt[dt$Highlight, ]
    
    # --- Plot title ---
    if (!is.null(title)) {
      title_txt <- gsub("{comparison}", comp, title, fixed = TRUE)
    } else {
      assay_label <- get_assay_label(dt)
      title_txt <- if (assay_label != "") paste0(comp, " (", assay_label, ")") else comp
    }
    
    used_changes <- levels(droplevels(dt_normal$Change))
    used_colors <- unname(color_map[used_changes])
    
    x_plotlines <- list(
      list(value = 0, color = "#1D3557", width = 1, zIndex = 3)
    )
    
    if (lfc_thr > 0) {
      threshold_line <- list(color = "#6C757D", width = 1, dashStyle = "Dash", zIndex = 2)
      x_plotlines <- c(x_plotlines, list(
        modifyList(threshold_line, list(value = -lfc_thr)),
        modifyList(threshold_line, list(value = lfc_thr))
      ))
    }
    
    # --- Build the base highchart with the regular points ---
    hc <- highcharter::hchart(
      dt_normal,
      type = "scatter",
      highcharter::hcaes(
        x       = logFC,
        y       = minusLog10P,
        group   = Change,
        protein = Protein.IDs,
        gene    = Gene.Names,
        pval    = pval_fmt
      )
    ) |>
      highcharter::hc_chart(
        zoomType = "xy",
        backgroundColor = "#FFFFFF",
        style = list(fontFamily = "Inter, -apple-system, sans-serif")
      ) |>
      highcharter::hc_title(
        text = title_txt,
        style = list(fontSize = "16px", fontWeight = "600", color = "#1D3557")
      ) |>
      highcharter::hc_xAxis(
        title = list(
          text = "log<sub>2</sub> Fold Change",
          useHTML = TRUE,
          style = list(fontSize = "13px", color = "#495057")
        ),
        gridLineWidth = 0,
        lineColor = "#DEE2E6",
        tickColor = "#DEE2E6",
        plotLines = x_plotlines
      ) |>
      highcharter::hc_yAxis(
        title = list(
          text = "-log<sub>10</sub>(p-value)",
          useHTML = TRUE,
          style = list(fontSize = "13px", color = "#495057")
        ),
        lineColor = "#DEE2E6",
        lineWidth = 1,
        tickColor = "#DEE2E6",
        gridLineColor = "#F1F3F4",
        gridLineDashStyle = "Dot",
        plotLines = list(
          list(
            value = -log10(alpha),
            color = "#6C757D",
            width = 1,
            dashStyle = "Dash",
            zIndex = 2,
            label = list(
              text = paste0("\u03b1 = ", alpha),
              style = list(color = "#6C757D", fontSize = "10px"),
              align = "right",
              x = -10,
              y = 12
            )
          )
        )
      ) |>
      highcharter::hc_colors(used_colors) |>
      highcharter::hc_tooltip(
        useHTML = TRUE,
        backgroundColor = "rgba(255, 255, 255, 0.95)",
        borderColor = "#DEE2E6",
        borderRadius = 8,
        shadow = TRUE,
        style = list(fontSize = "12px"),
        headerFormat = "",
        pointFormat = paste0(
          "<div style='padding: 4px;'>",
          "<b style='font-size: 13px; color: #1D3557;'>{point.gene}</b><br/>",
          "<span style='color: #6C757D;'>Protein:</span> {point.protein}<br/>",
          "<span style='color: #6C757D;'>log\u2082FC:</span> <b>{point.x:.3f}</b><br/>",
          "<span style='color: #6C757D;'>", p_col, ":</span> <b>{point.pval}</b>",
          "</div>"
        )
      ) |>
      highcharter::hc_legend(
        title = list(text = "", style = list(fontStyle = "normal")),
        layout = "horizontal",
        align = "center",
        verticalAlign = "bottom",
        itemStyle = list(fontSize = "12px", fontWeight = "normal")
      ) |>
      highcharter::hc_plotOptions(
        scatter = list(
          marker = list(
            radius = point_size,
            symbol = "circle",
            lineWidth = 0,
            states = list(
              hover = list(
                radiusPlus = 2,
                lineWidthPlus = 1,
                lineColor = "#1D3557"
              )
            )
          ),
          turboThreshold = nrow(dt) + 100
        )
      ) |>
      highcharter::hc_exporting(
        enabled = TRUE,
        buttons = list(contextButton = list(menuItems = c("downloadPNG", "downloadSVG", "downloadPDF")))
      )
    
    # --- Add the series of highlighted points ---
    if (nrow(dt_highlighted) > 0) {

      # Pick the border colour according to the change category
      dt_highlighted$borderColor <- vapply(dt_highlighted$Change, function(ch) {
        if (ch == "Up") colors$up
        else if (ch == "Down") colors$down
        else colors$ns
      }, character(1))
      
      # Build the list of points with per-point formatting
      highlighted_points <- lapply(seq_len(nrow(dt_highlighted)), function(i) {
        row <- dt_highlighted[i, ]
        list(
          x        = row$logFC,
          y        = row$minusLog10P,
          gene     = row$Gene.Names,
          protein  = row$Protein.IDs,
          pval     = row$pval_fmt,
          marker   = list(
            fillColor   = "#FFFFFF",
            lineColor   = row$borderColor,
            lineWidth   = 2,
            radius      = point_size + 1,
            symbol      = "circle"
          ),
          dataLabels = list(
            enabled = TRUE,
            format  = row$Gene.Names,
            style   = list(
              fontSize     = "11px",
              fontWeight   = "600",
              color        = "#1D3557",
              textOutline  = "2px #FFFFFF"
            ),
            y            = -12,
            x            = 0,
            align        = "center",
            verticalAlign = "bottom",
            allowOverlap = FALSE
          )
        )
      })
      
      hc <- hc |>
        highcharter::hc_add_series(
          name = "Highlighted",
          type = "scatter",
          data = highlighted_points,
          showInLegend = FALSE,
          enableMouseTracking = TRUE,
          tooltip = list(
            headerFormat = "",
            pointFormat = paste0(
              "<div style='padding: 6px;'>",
              "<b style='font-size: 14px; color: #1D3557;'>{point.gene}</b>",
              "<span style='background: #E63946; color: white; padding: 2px 6px; border-radius: 3px; margin-left: 8px; font-size: 10px; position: relative; top: -2px;'>\u2605 Highlighted</span>",
              "<br/>",
              "<span style='color: #6C757D;'>Protein:</span> {point.protein}<br/>",
              "<span style='color: #6C757D;'>log\u2082FC:</span> <b>{point.x:.3f}</b><br/>",
              "<span style='color: #6C757D;'>", p_col, ":</span> <b>{point.pval}</b>",
              "</div>"
            )
          ),
          zIndex = 10
        )
    }
    
    hc
  })
  
  names(hclist) <- comparisons
  Filter(Negate(is.null), hclist)
}

# --- Helper functions ---

detect_pvalue_col <- function(df, preferred) {
  if (preferred %in% names(df)) return(preferred)
  candidates <- c("adj.P.Val", "P.Value", "pvalue", "p.value", "padj")
  found <- intersect(candidates, names(df))
  if (length(found) == 0) stop("No p-value column found.")
  found[1]
}

get_assay_label <- function(dt) {
  if (!"Assay" %in% names(dt)) return("")
  assays <- unique(dt$Assay)
  if (length(assays) == 1) assays else ""
}


# =============================================================================
# USAGE EXAMPLES
# =============================================================================

# DEPs_results <- read_tsv("DEPs_results.tsv")
# DEPs_results <- arrow::read_parquet("./data-raw/VolcanoPlot_Input.parquet")

# --- Basic example ---
# hc_volcanos <- volcano_highchart_list(
#   de_res      = my_dataframe,
#   ain         = "LoessCyc",
#   comparisons = c("B-A"),
#   alpha       = 0.05,
#   p_col       = "adj.P.Val",
#   point_size  = 3
# )

# --- With top genes ---
# hc_volcanos <- volcano_highchart_list(
#   de_res         = my_dataframe,
#   ain            = "LoessCyc",
#   comparisons    = c("B-A"),
#   alpha          = 0.05,
#   point_size     = 3,
#   show_top_genes = 10
# )

# --- With hand-picked genes ---
# hc_volcanos <- volcano_highchart_list(
#   de_res          = my_dataframe,
#   ain             = "LoessCyc",
#   comparisons     = c("B-A"),
#   alpha           = 0.05,
#   point_size      = 3,
#   highlight_genes = c("EGFR", "TP53", "BRCA1")
# )

# --- With a custom title ---
# hc_volcanos <- volcano_highchart_list(
#   de_res      = my_dataframe,
#   ain         = "LoessCyc",
#   comparisons = c("B-A"),
#   alpha       = 0.05,
#   point_size  = 3,
#   title       = "Treatment vs Control"
# )

# --- With a paletteer palette ---
# hc_volcanos <- volcano_highchart_list(
#   de_res      = my_dataframe,
#   ain         = "LoessCyc",
#   comparisons = c("B-A"),
#   alpha       = 0.05,
#   point_size  = 3,
#   palette     = "ggsci::default_jco"
# )

# --- Combining every option ---
# hc_volcanos <- volcano_highchart_list(
#   de_res          = my_dataframe,
#   ain             = "LoessCyc",
#   comparisons     = c("B-A", "C-A", "D-A"),
#   lfc_thr         = 0,
#   alpha           = 0.05,
#   point_size      = 3,
#   show_top_genes  = 5,
#   highlight_genes = c("EGFR", "plaP", "SEC6"),
#   title           = "Differential Analysis",
#   palette         = "ggsci::nrc_npg"
# )

# --- Display ---
# hc_volcanos[["B-A"]]
