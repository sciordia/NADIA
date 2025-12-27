# =============================================================================
# Boxplot Interactivo con Highcharts para Datos de Proteómica
# =============================================================================

library(highcharter)
library(dplyr)
library(tidyr)
library(RColorBrewer)

# -----------------------------------------------------------------------------
# Función para calcular estadísticas del boxplot
# -----------------------------------------------------------------------------

calc_boxplot_stats <- function(x, coef = 1.5) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NULL)

  q <- quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE)
  iqr <- q[3] - q[1]

  lower_fence <- q[1] - coef * iqr
  upper_fence <- q[3] + coef * iqr

  # Valores dentro de los bigotes
  whisker_low <- min(x[x >= lower_fence], na.rm = TRUE)
  whisker_high <- max(x[x <= upper_fence], na.rm = TRUE)

  # Outliers
  outliers <- x[x < lower_fence | x > upper_fence]

  # Valores no-outliers (dentro de los bigotes)
  non_outliers <- x[x >= lower_fence & x <= upper_fence]

  list(
    low = whisker_low,
    q1 = unname(q[1]),
    median = unname(q[2]),
    q3 = unname(q[3]),
    high = whisker_high,
    outliers = outliers,
    non_outliers = non_outliers,
    all_values = x,
    mean = mean(x, na.rm = TRUE),
    n = length(x)
  )
}


# -----------------------------------------------------------------------------
# Función principal: Boxplot interactivo con Highcharts
# -----------------------------------------------------------------------------

#' Boxplot Interactivo con Highcharts para Proteómica
#'
#' @param data Data frame con columnas: Column, Assay, Intensity, Condition
#' @param assays Vector de assays a incluir (NULL = todos)
#' @param color_by Columna para colorear (default: "Condition")
#' @param group_order Orden de los grupos/condiciones
#' @param palette Paleta de colores: "ggsci::palette", "brewer:Name", o vector
#' @param title Título del gráfico (opcional)
#' @param subtitle Subtítulo del gráfico (opcional)
#' @param show_points Mostrar todos los puntos individuales (default: TRUE)
#' @param show_outliers Mostrar outliers fuera de los bigotes (default: TRUE)
#' @param point_jitter Cantidad de jitter horizontal para los puntos (default: 0.15)
#' @param point_size Radio de los puntos (default: 3)
#' @param horizontal Orientación horizontal (default: TRUE)
#' @param height Altura del gráfico en píxeles
#'
#' @return Lista de objetos highchart (uno por assay)
#'
boxplot_highchart_list <- function(
    data,
    assays = NULL,
    color_by = "Condition",
    group_order = NULL,
    palette = NULL,
    title = NULL,
    subtitle = NULL,
    show_points = TRUE,
    show_outliers = TRUE,
    point_jitter = 0.15,
    point_size = 3,
    horizontal = TRUE,
    height = NULL
) {

  # ---------------------------------------------------------------------------
  # 1) Validación de inputs
  # ---------------------------------------------------------------------------
  required_cols <- c("Column", "Assay", "Intensity")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Columnas requeridas faltantes: ", paste(missing_cols, collapse = ", "))
  }

  # Verificar columna de color
  if (!is.null(color_by) && !(color_by %in% names(data))) {
    warning("'", color_by, "' no está en el data frame. Se ignorará.")
    color_by <- NULL
  }

  # Filtrar assays
  if (!is.null(assays)) {
    data <- data[data$Assay %in% assays, , drop = FALSE]
  }

  # Eliminar NA e Inf
  data <- data[is.finite(data$Intensity), , drop = FALSE]

  available_assays <- unique(data$Assay)

  # ---------------------------------------------------------------------------
  # 2) Configurar orden de grupos
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
  # 3) Configurar paleta de colores
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
      stop("Para usar paletteer, instala con: install.packages('paletteer')")
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

  # Etiquetas de assay
  assay_labels <- c(
    "log2" = "Log\u2082 Intensity",
    "LoessCyc" = "LOESS Cyclic Normalization",
    "raw" = "Raw Intensity",
    "vsn" = "VSN Normalized",
    "quantile" = "Quantile Normalized"
  )

  # ---------------------------------------------------------------------------
  # 4) Generar un gráfico por cada assay
  # ---------------------------------------------------------------------------
  hc_list <- lapply(available_assays, function(current_assay) {

    dt <- data[data$Assay == current_assay, , drop = FALSE]
    if (nrow(dt) == 0) return(NULL)

    # Obtener samples únicos ordenados por grupo
    samples_df <- dt %>%
      select(Column, all_of(color_by)) %>%
      distinct() %>%
      arrange(.data[[color_by]], Column)

    samples <- samples_df$Column
    sample_groups <- samples_df[[color_by]]

    # Calcular estadísticas por sample
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
        index = i - 1  # índice 0-based para las categorías
      )
    })

    box_data <- Filter(Negate(is.null), box_data)

    # -------------------------------------------------------------------------
    # Preparar datos para highcharts boxplot (UNA sola serie con colores individuales)
    # Formato: lista de objetos con low, q1, median, q3, high y color individual
    # -------------------------------------------------------------------------
    boxplot_points <- lapply(box_data, function(bd) {
      list(
        low = bd$stats$low,
        q1 = bd$stats$q1,
        median = bd$stats$median,
        q3 = bd$stats$q3,
        high = bd$stats$high,
        name = bd$sample,
        n = bd$stats$n,
        mean = round(bd$stats$mean, 3),
        color = bd$color,
        fillColor = paste0(bd$color, "B3")  # 70% opacity
      )
    })

    # -------------------------------------------------------------------------
    # Preparar puntos scatter (por grupo para la leyenda)
    # Formato: [x, y] donde x es el índice de la categoría
    # -------------------------------------------------------------------------
    scatter_series_list <- list()

    if (show_points) {
      # Crear una serie scatter por cada grupo (para leyenda con colores)
      scatter_series_list <- lapply(group_levels, function(grp) {
        grp_data <- Filter(function(x) x$group == grp, box_data)
        if (length(grp_data) == 0) return(NULL)

        # Generar puntos para todos los valores (con jitter)
        points <- do.call(c, lapply(grp_data, function(bd) {
          values <- bd$stats$all_values
          lapply(values, function(v) {
            # Añadir jitter horizontal
            jittered_x <- bd$index + runif(1, -point_jitter, point_jitter)
            list(x = jittered_x, y = v, name = bd$sample)
          })
        }))

        if (length(points) == 0) return(NULL)

        list(
          name = grp,
          type = "scatter",
          data = points,
          color = unname(col_values[grp]),
          marker = list(
            symbol = "circle",
            radius = point_size,
            fillColor = unname(col_values[grp]),
            lineWidth = 0.5,
            lineColor = "#FFFFFF"
          ),
          tooltip = list(
            pointFormat = "<b>{point.name}</b><br/>Valor: {point.y:.3f}"
          ),
          showInLegend = TRUE
        )
      })

      scatter_series_list <- Filter(Negate(is.null), scatter_series_list)

    } else if (show_outliers) {
      # Solo mostrar outliers si show_points = FALSE pero show_outliers = TRUE
      outlier_series_by_group <- lapply(group_levels, function(grp) {
        grp_data <- Filter(function(x) x$group == grp, box_data)
        if (length(grp_data) == 0) return(NULL)

        points <- do.call(c, lapply(grp_data, function(bd) {
          if (length(bd$stats$outliers) == 0) return(NULL)
          lapply(bd$stats$outliers, function(o) {
            jittered_x <- bd$index + runif(1, -point_jitter, point_jitter)
            list(x = jittered_x, y = o, name = bd$sample)
          })
        }))

        if (length(points) == 0) return(NULL)

        list(
          name = paste0(grp, " (outliers)"),
          type = "scatter",
          data = points,
          color = unname(col_values[grp]),
          marker = list(
            symbol = "circle",
            radius = point_size,
            fillColor = unname(col_values[grp]),
            lineWidth = 1,
            lineColor = "#FFFFFF"
          ),
          tooltip = list(
            pointFormat = "<b>{point.name}</b><br/>Outlier: {point.y:.3f}"
          ),
          showInLegend = TRUE
        )
      })

      scatter_series_list <- Filter(Negate(is.null), outlier_series_by_group)
    }

    # Título del gráfico
    if (!is.null(title)) {
      chart_title <- title
    } else {
      chart_title <- if (current_assay %in% names(assay_labels)) {
        assay_labels[current_assay]
      } else {
        current_assay
      }
    }

    # Construir highchart
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
          groupPadding = 0.1,
          pointPadding = 0.05,
          borderRadius = 2,
          lineWidth = 1.5,
          whiskerLength = "50%",
          whiskerWidth = 2,
          medianColor = "#1D3557",
          medianWidth = 2,
          colorByPoint = TRUE  # Permite colores individuales por boxplot
        ),
        scatter = list(
          jitter = list(x = 0, y = 0)  # El jitter ya se aplica manualmente
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
    # Añadir serie boxplot (una sola serie con colores individuales)
    # -------------------------------------------------------------------------
    hc <- hc %>% hc_add_series(
      name = "Boxplot",
      type = "boxplot",
      data = boxplot_points,
      showInLegend = FALSE,
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

    # -------------------------------------------------------------------------
    # Añadir series scatter (puntos por grupo)
    # -------------------------------------------------------------------------
    for (scatter_series in scatter_series_list) {
      hc <- hc %>% hc_add_series(
        name = scatter_series$name,
        type = "scatter",
        data = scatter_series$data,
        color = scatter_series$color,
        marker = scatter_series$marker,
        tooltip = scatter_series$tooltip,
        showInLegend = scatter_series$showInLegend
      )
    }

    # Subtítulo
    if (!is.null(subtitle)) {
      hc <- hc %>% hc_subtitle(
        text = subtitle,
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
# EJEMPLOS DE USO
# =============================================================================

# --- Ejemplo básico (con puntos por defecto) ---
# hc_boxplots <- boxplot_highchart_list(
#   data        = mi_dataframe,
#   assays      = c("log2", "LoessCyc"),
#   color_by    = "Condition",
#   group_order = c("A", "B", "C", "D")
# )
# hc_boxplots[["log2"]]
# hc_boxplots[["LoessCyc"]]

# --- Con título personalizado ---
# hc_boxplots <- boxplot_highchart_list(
#   data     = mi_dataframe,
#   assays   = "LoessCyc",
#   color_by = "Condition",
#   title    = "Distribución de Intensidades por Muestra",
#   subtitle = "Normalización LOESS cíclica"
# )
# hc_boxplots[["LoessCyc"]]

# --- Orientación vertical ---
# hc_boxplots <- boxplot_highchart_list(
#   data       = mi_dataframe,
#   color_by   = "Condition",
#   horizontal = FALSE
# )

# --- Con paleta personalizada ---
# hc_boxplots <- boxplot_highchart_list(
#   data     = mi_dataframe,
#   color_by = "Condition",
#   palette  = "ggsci::nrc_npg"
# )

# --- Sin puntos (solo boxplots) ---
# hc_boxplots <- boxplot_highchart_list(
#   data        = mi_dataframe,
#   color_by    = "Condition",
#   show_points = FALSE
# )

# --- Solo outliers (sin todos los puntos) ---
# hc_boxplots <- boxplot_highchart_list(
#   data          = mi_dataframe,
#   color_by      = "Condition",
#   show_points   = FALSE,
#   show_outliers = TRUE
# )

# --- Personalizar tamaño y jitter de puntos ---
# hc_boxplots <- boxplot_highchart_list(
#   data         = mi_dataframe,
#   color_by     = "Condition",
#   point_jitter = 0.2,
#   point_size   = 4
# )
