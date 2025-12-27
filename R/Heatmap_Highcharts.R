# =============================================================================
# Heatmap Interactivo con Highcharts para Datos de Proteómica
# =============================================================================

library(highcharter)
library(dplyr)
library(tidyr)

# -----------------------------------------------------------------------------
# Función principal: Heatmap interactivo con Highcharts
# -----------------------------------------------------------------------------

#' Heatmap Interactivo con Highcharts para Proteómica
#'
#' @param data Data frame en formato largo con columnas: Row (proteína/gen),
#'             Column (muestra), Value (intensidad). O matriz/data.frame ancho.
#' @param row_col Nombre de la columna con identificadores de fila (proteínas)
#' @param col_col Nombre de la columna con identificadores de columna (muestras)
#' @param value_col Nombre de la columna con valores
#' @param annotation_col Data frame con anotaciones de columnas (opcional)
#' @param annotation_var Variable de anotación para colorear columnas
#' @param cluster_rows Aplicar clustering a filas (default: FALSE)
#' @param cluster_cols Aplicar clustering a columnas (default: FALSE)
#' @param scale_data Escalar datos: "none", "row", "column" (default: "none")
#' @param color_palette Paleta de colores: "RdBu", "viridis", "plasma", o vector
#' @param reverse_palette Invertir paleta (default: TRUE para divergente)
#' @param n_colors Número de colores en el gradiente (default: 9)
#' @param title Título del gráfico
#' @param subtitle Subtítulo del gráfico
#' @param show_rownames Mostrar nombres de filas (default: TRUE si <= 50)
#' @param show_colnames Mostrar nombres de columnas (default: TRUE)
#' @param row_label Etiqueta para el eje Y
#' @param col_label Etiqueta para el eje X
#' @param height Altura del gráfico en píxeles
#' @param width Anchura del gráfico en píxeles
#'
#' @return Objeto highchart
#'
heatmap_highchart <- function(
    data,
    row_col = NULL,
    col_col = NULL,
    value_col = NULL,
    annotation_col = NULL,
    annotation_var = NULL,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    scale_data = "none",
    color_palette = "RdBu",
    reverse_palette = TRUE,
    n_colors = 9,
    title = NULL,
    subtitle = NULL,
    show_rownames = NULL,
    show_colnames = TRUE,
    row_label = "Proteins",
    col_label = "Samples",
    height = NULL,
    width = NULL
) {

  # ---------------------------------------------------------------------------
  # 1) Detectar formato de datos y convertir a matriz
  # ---------------------------------------------------------------------------

  if (is.matrix(data)) {
    # Ya es matriz
    mat <- data
  } else if (!is.null(row_col) && !is.null(col_col) && !is.null(value_col)) {
    # Formato largo -> convertir a matriz
    data_wide <- data %>%
      select(all_of(c(row_col, col_col, value_col))) %>%
      pivot_wider(
        names_from = all_of(col_col),
        values_from = all_of(value_col)
      )

    row_names <- data_wide[[row_col]]
    mat <- as.matrix(data_wide[, -1])
    rownames(mat) <- row_names
  } else {
    # Asumir data.frame ancho (primera columna = rownames)
    if (is.character(data[[1]]) || is.factor(data[[1]])) {
      row_names <- data[[1]]
      mat <- as.matrix(data[, -1])
      rownames(mat) <- row_names
    } else {
      mat <- as.matrix(data)
    }
  }

  # Asegurar que es numérico
  storage.mode(mat) <- "double"

  # ---------------------------------------------------------------------------
  # 2) Escalar datos si se solicita
  # ---------------------------------------------------------------------------

  if (scale_data == "row") {
    mat <- t(scale(t(mat)))
  } else if (scale_data == "column") {
    mat <- scale(mat)
  }

  # Reemplazar NaN por NA
  mat[is.nan(mat)] <- NA

  # ---------------------------------------------------------------------------
  # 3) Clustering jerárquico
  # ---------------------------------------------------------------------------

  if (cluster_rows && nrow(mat) > 2) {
    # Eliminar filas con todos NA para clustering
    valid_rows <- rowSums(!is.na(mat)) > 0
    if (sum(valid_rows) > 2) {
      mat_for_clust <- mat[valid_rows, ]
      mat_for_clust[is.na(mat_for_clust)] <- 0

      row_dist <- dist(mat_for_clust)
      row_hclust <- hclust(row_dist, method = "complete")
      row_order <- row_hclust$order

      # Reordenar matriz completa
      mat <- mat[rownames(mat_for_clust)[row_order], , drop = FALSE]
    }
  }

  if (cluster_cols && ncol(mat) > 2) {
    valid_cols <- colSums(!is.na(mat)) > 0
    if (sum(valid_cols) > 2) {
      mat_for_clust <- mat[, valid_cols]
      mat_for_clust[is.na(mat_for_clust)] <- 0

      col_dist <- dist(t(mat_for_clust))
      col_hclust <- hclust(col_dist, method = "complete")
      col_order <- col_hclust$order

      mat <- mat[, colnames(mat_for_clust)[col_order], drop = FALSE]
    }
  }

  # ---------------------------------------------------------------------------
  # 4) Preparar datos para Highcharts
  # ---------------------------------------------------------------------------

  rows <- rownames(mat)
  cols <- colnames(mat)

  if (is.null(rows)) rows <- paste0("Row_", seq_len(nrow(mat)))
  if (is.null(cols)) cols <- paste0("Col_", seq_len(ncol(mat)))

  # Decidir si mostrar nombres de filas
  if (is.null(show_rownames)) {
    show_rownames <- length(rows) <= 50
  }

  # Convertir a formato largo para Highcharts
  heatmap_data <- lapply(seq_along(cols), function(j) {
    lapply(seq_along(rows), function(i) {
      val <- mat[i, j]
      list(
        x = j - 1,
        y = i - 1,
        value = if (is.na(val)) NULL else round(val, 4),
        row_name = rows[i],
        col_name = cols[j]
      )
    })
  })

  heatmap_data <- unlist(heatmap_data, recursive = FALSE)

  # Filtrar NAs
  heatmap_data <- Filter(function(p) !is.null(p$value), heatmap_data)

  # ---------------------------------------------------------------------------
  # 5) Configurar paleta de colores
  # ---------------------------------------------------------------------------

  get_color_stops <- function(palette_name, n, reverse = FALSE) {

    palettes <- list(
      # Divergentes
      RdBu = c("#67001F", "#B2182B", "#D6604D", "#F4A582", "#FDDBC7",
               "#F7F7F7", "#D1E5F0", "#92C5DE", "#4393C3", "#2166AC", "#053061"),
      RdYlBu = c("#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE090",
                 "#FFFFBF", "#E0F3F8", "#ABD9E9", "#74ADD1", "#4575B4", "#313695"),
      PiYG = c("#8E0152", "#C51B7D", "#DE77AE", "#F1B6DA", "#FDE0EF",
               "#F7F7F7", "#E6F5D0", "#B8E186", "#7FBC41", "#4D9221", "#276419"),
      BrBG = c("#543005", "#8C510A", "#BF812D", "#DFC27D", "#F6E8C3",
               "#F5F5F5", "#C7EAE5", "#80CDC1", "#35978F", "#01665E", "#003C30"),

      # Secuenciales
      viridis = c("#440154", "#482878", "#3E4A89", "#31688E", "#26828E",
                  "#1F9E89", "#35B779", "#6DCD59", "#B4DE2C", "#FDE725"),
      plasma = c("#0D0887", "#46039F", "#7201A8", "#9C179E", "#BD3786",
                 "#D8576B", "#ED7953", "#FB9F3A", "#FDC328", "#F0F921"),
      inferno = c("#000004", "#1B0C41", "#4A0C6B", "#781C6D", "#A52C60",
                  "#CF4446", "#ED6925", "#FB9A06", "#F7D03C", "#FCFFA4"),

      # Personalizada para proteómica
      proteomics = c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7",
                     "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F")
    )

    if (palette_name %in% names(palettes)) {
      pal <- palettes[[palette_name]]
    } else {
      # Usar RColorBrewer si está disponible
      if (requireNamespace("RColorBrewer", quietly = TRUE) &&
          palette_name %in% rownames(RColorBrewer::brewer.pal.info)) {
        max_n <- RColorBrewer::brewer.pal.info[palette_name, "maxcolors"]
        pal <- RColorBrewer::brewer.pal(max_n, palette_name)
      } else {
        pal <- palettes$RdBu
      }
    }

    if (reverse) pal <- rev(pal)

    # Interpolar a n colores
    color_func <- colorRampPalette(pal)
    colors <- color_func(n)

    # Crear stops para Highcharts (0 a 1)
    stops <- lapply(seq_along(colors), function(i) {
      list((i - 1) / (length(colors) - 1), colors[i])
    })

    stops
  }

  color_stops <- get_color_stops(color_palette, n_colors, reverse_palette)

  # ---------------------------------------------------------------------------
  # 6) Calcular rango de valores
  # ---------------------------------------------------------------------------

  values <- sapply(heatmap_data, function(x) x$value)
  min_val <- min(values, na.rm = TRUE)
  max_val <- max(values, na.rm = TRUE)

  # Para escalas divergentes, centrar en 0
  if (scale_data == "row" || scale_data == "column") {
    abs_max <- max(abs(c(min_val, max_val)), na.rm = TRUE)
    min_val <- -abs_max
    max_val <- abs_max
  }

  # ---------------------------------------------------------------------------
  # 7) Construir Highchart
  # ---------------------------------------------------------------------------

  # Calcular altura automática si no se especifica
  if (is.null(height)) {
    height <- max(400, min(800, 50 + length(rows) * 15))
  }

  hc <- highchart() %>%
    hc_chart(
      type = "heatmap",
      backgroundColor = "#FFFFFF",
      style = list(fontFamily = "Inter, -apple-system, sans-serif"),
      height = height,
      width = width,
      marginTop = 80,
      marginBottom = if (show_colnames) 120 else 60,
      marginLeft = if (show_rownames) 150 else 60,
      marginRight = 80
    ) %>%
    hc_title(
      text = title,
      style = list(
        fontSize = "18px",
        fontWeight = "600",
        color = "#1D3557"
      )
    ) %>%
    hc_subtitle(
      text = subtitle,
      style = list(
        fontSize = "13px",
        color = "#6C757D"
      )
    ) %>%
    hc_xAxis(
      categories = cols,
      title = list(
        text = col_label,
        style = list(
          fontSize = "13px",
          fontWeight = "bold",
          color = "#212529"
        ),
        margin = 15
      ),
      labels = list(
        enabled = show_colnames,
        rotation = -45,
        style = list(
          fontSize = "10px",
          color = "#495057"
        )
      ),
      lineColor = "#DEE2E6",
      tickColor = "#DEE2E6",
      opposite = FALSE
    ) %>%
    hc_yAxis(
      categories = rows,
      title = list(
        text = row_label,
        style = list(
          fontSize = "13px",
          fontWeight = "bold",
          color = "#212529"
        )
      ),
      labels = list(
        enabled = show_rownames,
        style = list(
          fontSize = "10px",
          color = "#495057"
        )
      ),
      lineColor = "#DEE2E6",
      reversed = TRUE,
      gridLineWidth = 0
    ) %>%
    hc_colorAxis(
      min = min_val,
      max = max_val,
      stops = color_stops,
      labels = list(
        style = list(
          fontSize = "11px",
          color = "#495057"
        )
      )
    ) %>%
    hc_legend(
      align = "right",
      layout = "vertical",
      verticalAlign = "middle",
      symbolHeight = 200
    ) %>%
    hc_tooltip(
      useHTML = TRUE,
      backgroundColor = "rgba(255, 255, 255, 0.95)",
      borderColor = "#DEE2E6",
      borderRadius = 8,
      shadow = TRUE,
      style = list(fontSize = "12px"),
      formatter = JS(paste0(
        "function() {",
        "  return '<div style=\"padding: 8px;\">' +",
        "    '<b style=\"font-size: 13px; color: #1D3557;\">' + this.point.row_name + '</b><br/>' +",
        "    '<span style=\"color: #6C757D;\">Sample:</span> <b>' + this.point.col_name + '</b><br/>' +",
        "    '<span style=\"color: #6C757D;\">Value:</span> <b>' + this.point.value.toFixed(3) + '</b>' +",
        "    '</div>';",
        "}"
      ))
    ) %>%
    hc_plotOptions(
      heatmap = list(
        borderWidth = 0.5,
        borderColor = "#FFFFFF",
        nullColor = "#F8F9FA"
      )
    ) %>%
    hc_add_series(
      name = "Intensity",
      data = heatmap_data,
      boostThreshold = 100,
      turboThreshold = 0
    ) %>%
    hc_exporting(
      enabled = TRUE,
      buttons = list(
        contextButton = list(
          menuItems = c("downloadPNG", "downloadSVG", "downloadPDF")
        )
      )
    )

  # ---------------------------------------------------------------------------
  # 8) Añadir anotaciones de columna si se proporcionan
  # ---------------------------------------------------------------------------

  if (!is.null(annotation_col) && !is.null(annotation_var)) {
    # Esto requeriría una implementación más compleja con múltiples ejes
    # Por ahora, añadimos info al tooltip
    message("Nota: Las anotaciones de columna se incluirán en versiones futuras.")
  }

  hc
}


# -----------------------------------------------------------------------------
# Función auxiliar: Crear heatmap desde data frame largo de proteómica
# -----------------------------------------------------------------------------

#' Heatmap desde datos de proteómica en formato largo
#'
#' @param data Data frame con columnas similares a boxplot (Column, Assay, Intensity, etc.)
#' @param protein_col Columna con identificadores de proteína/gen
#' @param sample_col Columna con identificadores de muestra
#' @param value_col Columna con valores de intensidad
#' @param assay Nombre del assay a usar (si hay múltiples)
#' @param top_n Número de proteínas más variables a mostrar (default: 50)
#' @param ... Argumentos adicionales para heatmap_highchart
#'
#' @return Objeto highchart
#'
heatmap_proteomics <- function(
    data,
    protein_col = "Protein",
    sample_col = "Column",
    value_col = "Intensity",
    assay = NULL,
    top_n = 50,
    ...
) {

  # Filtrar por assay si se especifica
  if (!is.null(assay) && "Assay" %in% names(data)) {
    data <- data[data$Assay == assay, , drop = FALSE]
  }

  # Verificar columnas
  if (!(protein_col %in% names(data))) {
    stop("Columna '", protein_col, "' no encontrada en los datos.")
  }
  if (!(sample_col %in% names(data))) {
    stop("Columna '", sample_col, "' no encontrada en los datos.")
  }
  if (!(value_col %in% names(data))) {
    stop("Columna '", value_col, "' no encontrada en los datos.")
  }

  # Pivotar a formato ancho
  mat_df <- data %>%
    select(all_of(c(protein_col, sample_col, value_col))) %>%
    pivot_wider(
      names_from = all_of(sample_col),
      values_from = all_of(value_col),
      values_fn = mean
    )

  proteins <- mat_df[[protein_col]]
  mat <- as.matrix(mat_df[, -1])
  rownames(mat) <- proteins

  # Seleccionar top_n proteínas más variables
  if (!is.null(top_n) && nrow(mat) > top_n) {
    row_vars <- apply(mat, 1, var, na.rm = TRUE)
    top_idx <- order(row_vars, decreasing = TRUE)[1:top_n]
    mat <- mat[top_idx, , drop = FALSE]
  }

  # Llamar a la función principal
  heatmap_highchart(
    data = mat,
    ...
  )
}


# =============================================================================
# EJEMPLOS DE USO
# =============================================================================

# --- Ejemplo con matriz ---
# mat <- matrix(rnorm(200), nrow = 20, ncol = 10)
# rownames(mat) <- paste0("Protein_", 1:20)
# colnames(mat) <- paste0("Sample_", 1:10)
#
# hc <- heatmap_highchart(
#   data = mat,
#   title = "Heatmap de Intensidades",
#   subtitle = "Top 20 proteínas más variables",
#   scale_data = "row",
#   cluster_rows = TRUE,
#   cluster_cols = TRUE
# )
# hc

# --- Ejemplo con data frame largo ---
# hc <- heatmap_highchart(
#   data = mi_dataframe_largo,
#   row_col = "Protein",
#   col_col = "Sample",
#   value_col = "Intensity",
#   scale_data = "row",
#   cluster_rows = TRUE,
#   color_palette = "viridis",
#   reverse_palette = FALSE
# )
# hc

# --- Ejemplo con función de proteómica ---
# hc <- heatmap_proteomics(
#   data = mi_dataframe,
#   protein_col = "Protein.ID",
#   sample_col = "Column",
#   value_col = "Intensity",
#   assay = "LoessCyc",
#   top_n = 30,
#   scale_data = "row",
#   cluster_rows = TRUE,
#   cluster_cols = TRUE,
#   title = "Top 30 Proteínas Variables",
#   color_palette = "RdBu"
# )
# hc

# --- Paletas disponibles ---
# Divergentes: "RdBu", "RdYlBu", "PiYG", "BrBG", "proteomics"
# Secuenciales: "viridis", "plasma", "inferno"
# También cualquier paleta de RColorBrewer
