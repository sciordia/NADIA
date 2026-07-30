# =============================================================================
# Pattern Profiler Highcharts: Visualización de Clusters
# =============================================================================
#
# Este script genera visualizaciones interactivas con Highcharts para los
# resultados del clustering generados por Pattern_Profiler_Analysis.R.
#
# Input:
#   - Pattern_Profiler_Input.parquet: Tabla en formato LONG con:
#     - FeatureID, Cluster, Membership, z-scores por condición
#
# Output:
#   - Gráficos Highcharts interactivos (perfiles, centroides)
#
# Autor: Sergio Ciordia
# Licencia: MIT
# =============================================================================

# -----------------------------------------------------------------------------
# DEPENDENCIAS
# -----------------------------------------------------------------------------

library(highcharter)
library(dplyr)
library(arrow)


# =============================================================================
# FUNCIONES AUXILIARES DE COLOR
# =============================================================================

#' Configurar paleta de colores para clusters
#'
#' @param n_clusters Número de clusters
#' @param palette Paleta a usar (NULL, "ggsci::nombre", "brewer:nombre")
#'
#' @return Vector de colores nombrado por cluster
configure_cluster_palette <- function(n_clusters, palette = NULL) {

  # Paleta por defecto: colores distintivos y accesibles
  default_colors <- c(
    "#E63946", "#457B9D", "#2A9D8F", "#E9C46A", "#F4A261",
    "#264653", "#A8DADC", "#1D3557", "#F77F00", "#D62828",
    "#023E8A", "#0077B6", "#00B4D8", "#90E0EF", "#CAF0F8"
  )

  if (is.null(palette)) {
    # Si se piden mas clusters que colores base, interpolar en vez de reciclar
    # (evita el desajuste de longitud al asignar names() mas abajo).
    if (n_clusters > length(default_colors)) {
      colors <- colorRampPalette(default_colors)(n_clusters)
    } else {
      colors <- default_colors[seq_len(n_clusters)]
    }

  } else if (grepl("^ggsci::", palette)) {
    # Paletas de ggsci via paletteer
    if (requireNamespace("paletteer", quietly = TRUE)) {
      pal_name <- sub("^ggsci::", "", palette)
      colors <- tryCatch({
        as.character(paletteer::paletteer_d(paste0("ggsci::", pal_name), n_clusters))
      }, error = function(e) {
        warning("Paleta ggsci no encontrada, usando default")
        default_colors[seq_len(n_clusters)]
      })
    } else {
      warning("Paquete 'paletteer' no disponible, usando paleta default")
      colors <- default_colors[seq_len(n_clusters)]
    }

  } else if (grepl("^brewer:", palette)) {
    # Paletas de RColorBrewer
    pal_name <- sub("^brewer:", "", palette)
    if (requireNamespace("RColorBrewer", quietly = TRUE)) {
      max_colors <- RColorBrewer::brewer.pal.info[pal_name, "maxcolors"]
      if (is.na(max_colors)) {
        warning("Paleta brewer no encontrada, usando default")
        colors <- default_colors[seq_len(n_clusters)]
      } else {
        colors <- RColorBrewer::brewer.pal(min(n_clusters, max_colors), pal_name)
        if (n_clusters > max_colors) {
          colors <- colorRampPalette(colors)(n_clusters)
        }
      }
    } else {
      warning("Paquete 'RColorBrewer' no disponible, usando paleta default")
      colors <- default_colors[seq_len(n_clusters)]
    }

  } else {
    # Asumir vector de colores
    if (length(palette) >= n_clusters) {
      colors <- palette[seq_len(n_clusters)]
    } else {
      colors <- colorRampPalette(palette)(n_clusters)
    }
  }

  # Nombrar por cluster
  names(colors) <- seq_len(n_clusters)
  colors
}


# =============================================================================
# LECTURA DE DATOS
# =============================================================================

#' Leer datos de Pattern Profiler desde parquet
#'
#' @param file_path Ruta al archivo parquet
#' @param min_membership Filtro opcional por membership mínima
#'
#' @return DataFrame con datos de clustering
#'
#' @examples
#' \dontrun{
#' data <- read_pattern_profiler_data("data/Pattern_Profiler_Input.parquet")
#' data <- read_pattern_profiler_data("data/Pattern_Profiler_Input.parquet",
#'                                     min_membership = 0.5)
#' }
read_pattern_profiler_data <- function(file_path, min_membership = NULL) {

  if (!file.exists(file_path)) {
    stop("Archivo no encontrado: ", file_path)
  }

  data <- arrow::read_parquet(file_path)
  data <- as.data.frame(data)

  # Validar columnas requeridas
  required_cols <- c("FeatureID", "Cluster", "Membership")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Columnas faltantes: ", paste(missing_cols, collapse = ", "))
  }

  # Filtrar por membership si se especifica
  if (!is.null(min_membership)) {
    n_before <- nrow(data)
    data <- data[data$Membership >= min_membership, ]
    n_after <- nrow(data)
    message(sprintf("Filtrado por membership >= %.2f: %d -> %d filas",
                    min_membership, n_before, n_after))
  }

  data
}


#' Detectar columnas de condiciones (z-scores)
#'
#' @param data DataFrame de Pattern Profiler
#' @return Vector de nombres de columnas de condiciones
detect_condition_columns <- function(data) {
  # Excluir columnas conocidas
  exclude_cols <- c("FeatureID", "Cluster", "Membership")
  all_cols <- names(data)

  condition_cols <- setdiff(all_cols, exclude_cols)

  if (length(condition_cols) == 0) {
    stop("No se encontraron columnas de condiciones (z-scores)")
  }

  condition_cols
}


# =============================================================================
# FUNCIONES DE VISUALIZACIÓN
# =============================================================================

#' Gráfico de Perfil de Cluster con Highcharts
#'
#' Genera un gráfico interactivo mostrando los perfiles de expresión
#' de las proteínas en un cluster específico.
#'
#' @param data DataFrame de Pattern Profiler (formato LONG)
#' @param cluster Número de cluster a visualizar
#' @param conditions Vector de nombres de condiciones (orden para eje X)
#' @param min_membership Filtro adicional de membership (NULL = sin filtro)
#' @param show_centroid Mostrar línea del centroide (default: TRUE)
#' @param centroid_summary Método para centroide: "mean" o "median"
#' @param cluster_color Color del cluster (NULL = automático)
#' @param line_width Ancho de líneas de perfil (default: 1)
#' @param line_opacity Opacidad de líneas (default: 0.4)
#' @param centroid_width Ancho de línea del centroide (default: 3)
#' @param title Título personalizado (opcional)
#' @param height Altura del gráfico en píxeles
#'
#' @return Objeto highchart
#'
#' @examples
#' \dontrun{
#' data <- read_pattern_profiler_data("Pattern_Profiler_Input.parquet")
#' conditions <- c("A", "B", "C", "D")
#'
#' hc <- cluster_profile_highchart(data, cluster = 1, conditions = conditions)
#' hc
#' }
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

  # Detectar condiciones si no se especifican
  if (is.null(conditions)) {
    conditions <- detect_condition_columns(data)
  }

  # Filtrar por cluster
  cluster_data <- data[data$Cluster == cluster, ]

  if (nrow(cluster_data) == 0) {
    warning(sprintf("No hay datos para el cluster %d", cluster))
    return(NULL)
  }

  # Filtrar por membership adicional
  if (!is.null(min_membership)) {
    cluster_data <- cluster_data[cluster_data$Membership >= min_membership, ]
  }

  if (nrow(cluster_data) == 0) {
    warning(sprintf("Cluster %d: sin proteínas con membership >= %.2f",
                    cluster, min_membership))
    return(NULL)
  }

  # Extraer z-scores y metadata
  zscores <- as.matrix(cluster_data[, conditions, drop = FALSE])
  feature_ids <- cluster_data$FeatureID
  n_proteins <- nrow(cluster_data)

  # Calcular centroide
  if (centroid_summary == "mean") {
    centroid <- colMeans(zscores, na.rm = TRUE)
  } else {
    centroid <- apply(zscores, 2, median, na.rm = TRUE)
  }
  centroid <- unname(as.numeric(centroid))

  # Configurar color del cluster
  if (is.null(cluster_color)) {
    n_clusters <- max(data$Cluster)
    palette <- configure_cluster_palette(n_clusters)
    cluster_color <- unname(palette[as.character(cluster)])
  }

  # Color de líneas con opacidad
  line_color <- .hex_to_rgba(cluster_color, line_opacity)

  # Título
  if (is.null(title)) {
    min_mem <- min(cluster_data$Membership)
    title <- sprintf("Cluster %d (n = %d, membership >= %.2f)",
                     cluster, n_proteins, min_mem)
  }

  # ---------------------------------------------------------------------------
  # Construir series de líneas de perfil (sin interactividad)
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
  # Serie del centroide (interactiva)
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

    # Colores limpios sin nombres
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
  # Construir highchart
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

  # Añadir series de perfiles
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

  # Añadir centroide
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


#' Lista de Gráficos de Perfil de Clusters con Highcharts
#'
#' Genera una lista de gráficos interactivos para todos los clusters.
#'
#' @param data DataFrame de Pattern Profiler (formato LONG)
#' @param conditions Vector de nombres de condiciones
#' @param clusters Vector de clusters a visualizar (NULL = todos)
#' @param min_membership Filtro de membership mínima
#' @param show_centroid Mostrar línea central (default: TRUE)
#' @param centroid_summary Método para centroide: "mean" o "median"
#' @param palette Paleta para colores de clusters
#' @param line_width Ancho de líneas de perfil (default: 1)
#' @param line_opacity Opacidad de líneas (default: 0.4)
#' @param centroid_width Ancho de línea del centroide (default: 3)
#' @param height Altura de cada gráfico en píxeles
#'
#' @return Lista nombrada de objetos highchart
#'
#' @examples
#' \dontrun{
#' data <- read_pattern_profiler_data("Pattern_Profiler_Input.parquet")
#' conditions <- c("A", "B", "C", "D")
#'
#' hc_profiles <- cluster_profile_highchart_list(data, conditions)
#' hc_profiles[["Cluster_1"]]
#' hc_profiles[["Cluster_2"]]
#' }
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

  # Detectar condiciones si no se especifican
  if (is.null(conditions)) {
    conditions <- detect_condition_columns(data)
  }

  # Determinar clusters a visualizar
  if (is.null(clusters)) {
    clusters <- sort(unique(data$Cluster))
  }

  n_clusters <- max(data$Cluster)

  # Configurar paleta de colores
  cluster_colors <- configure_cluster_palette(n_clusters, palette)

  # Generar gráficos
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
      warning(sprintf("Error generando gráfico para Cluster %d: %s", k, e$message))
      return(NULL)
    })

    hc
  })

  names(hc_list) <- paste0("Cluster_", clusters)
  Filter(Negate(is.null), hc_list)
}


#' Gráfico de Centroides de Todos los Clusters
#'
#' Genera un gráfico comparativo con los centroides de todos los clusters.
#'
#' @param data DataFrame de Pattern Profiler (formato LONG)
#' @param conditions Vector de nombres de condiciones
#' @param clusters Clusters a incluir (NULL = todos)
#' @param min_membership Filtro de membership para calcular centroides
#' @param centroid_summary Método: "mean" o "median"
#' @param palette Paleta de colores
#' @param line_width Ancho de líneas (default: 2.5)
#' @param show_markers Mostrar marcadores en puntos (default: TRUE)
#' @param title Título personalizado
#' @param height Altura del gráfico
#'
#' @return Objeto highchart
#'
#' @examples
#' \dontrun{
#' data <- read_pattern_profiler_data("Pattern_Profiler_Input.parquet")
#' conditions <- c("A", "B", "C", "D")
#'
#' hc_centroids <- cluster_centroids_highchart(data, conditions)
#' hc_centroids
#' }
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

  # Detectar condiciones
  if (is.null(conditions)) {
    conditions <- detect_condition_columns(data)
  }

  # Filtrar por membership si se especifica
  if (!is.null(min_membership)) {
    data <- data[data$Membership >= min_membership, ]
  }

  # Determinar clusters
  if (is.null(clusters)) {
    clusters <- sort(unique(data$Cluster))
  }

  n_clusters <- max(data$Cluster)

  # Configurar paleta
  cluster_colors <- configure_cluster_palette(n_clusters, palette)

  # Calcular centroides por cluster
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
    warning("No hay datos para generar gráfico de centroides")
    return(NULL)
  }

  # Título
  if (is.null(title)) {
    title <- sprintf("Cluster Centroids (%s)", centroid_summary)
  }

  # Construir highchart
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

  # Añadir series de centroides
  for (item in centroids_list) {
    # Construir puntos con valores limpios (sin nombres)
    points <- lapply(seq_along(conditions), function(j) {
      list(
        x = as.integer(j - 1),
        y = round(item$centroid[j], 4),
        condition = unname(conditions[j])
      )
    })

    # Color limpio sin nombres
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
# FUNCIONES AUXILIARES DE RESUMEN
# =============================================================================

#' Resumen de datos de Pattern Profiler
#'
#' @param data DataFrame de Pattern Profiler
#' @return Lista con estadísticas resumidas
summarize_pattern_profiler <- function(data) {

  conditions <- detect_condition_columns(data)

  # Contar por cluster
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

  # Features en múltiples clusters
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
# EJEMPLOS DE USO
# =============================================================================

# --- Uso típico ---
# source("R/Pattern_Profiler_Highcharts.R")
#
# # Leer datos
# data <- read_pattern_profiler_data("data/Pattern_Profiler_Input.parquet")
#
# # Ver resumen
# summary <- summarize_pattern_profiler(data)
# print(summary$cluster_summary)
#
# # Definir orden de condiciones
# conditions <- c("A", "B", "C", "D")
#
# # Gráfico de un cluster específico
# hc_c1 <- cluster_profile_highchart(data, cluster = 1, conditions = conditions)
# hc_c1
#
# # Lista de gráficos para todos los clusters
# hc_profiles <- cluster_profile_highchart_list(data, conditions)
# hc_profiles[["Cluster_1"]]
# hc_profiles[["Cluster_2"]]
#
# # Con filtro de membership más estricto
# hc_profiles <- cluster_profile_highchart_list(
#   data, conditions,
#   min_membership = 0.5
# )
#
# # Gráfico de centroides comparativo
# hc_centroids <- cluster_centroids_highchart(data, conditions)
# hc_centroids
#
# # Con paleta personalizada
# hc_profiles <- cluster_profile_highchart_list(
#   data, conditions,
#   palette = "ggsci::nrc_npg"
# )
