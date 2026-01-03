# =============================================================================
# PCA Plot Interactivo con Highcharts para Datos de Proteómica
# =============================================================================

library(highcharter)
library(dplyr)


# -----------------------------------------------------------------------------
# Operador null-coalesce
# -----------------------------------------------------------------------------

`%||%` <- function(a, b) if (!is.null(a) && length(a) && !is.na(a[1])) a else b


# -----------------------------------------------------------------------------
# Función para convertir color hex a rgba
# -----------------------------------------------------------------------------

hex_to_rgba <- function(hex, alpha = 0.12) {
  hex <- gsub("^#", "", hex)
  r <- strtoi(substr(hex, 1, 2), base = 16)
  g <- strtoi(substr(hex, 3, 4), base = 16)
  b <- strtoi(substr(hex, 5, 6), base = 16)
  sprintf("rgba(%d, %d, %d, %.2f)", r, g, b, alpha)
}


# -----------------------------------------------------------------------------
# Funciones auxiliares para filtrado de proteínas
# -----------------------------------------------------------------------------

#' Construir nombre de columna adjP para una comparación
#'
#' @param comparison Nombre de la comparación (ej: "B-A")
#' @return Nombre de la columna (ej: "adjP_B-A")
adjp_col <- function(comparison) {
  paste0("adjP_", comparison)
}


#' Obtener IDs de proteínas según el modo de filtrado
#'
#' @param pca_input Data frame en formato long con columnas FeatureID, sig_any, adjP_*
#' @param mode Modo de filtrado: "all", "any", o "specific"
#' @param alpha Umbral de significancia para modo "specific" (default: 0.05)
#' @param comparison Nombre de la comparación para modo "specific"
#'
#' @return Vector de FeatureIDs que cumplen el criterio
get_feature_ids <- function(pca_input,
                            mode = c("all", "any", "specific"),
                            alpha = 0.05,
                            comparison = NULL) {

  mode <- match.arg(mode)

  # Convertir a data.frame para evitar problemas con tibbles
  pca_input <- as.data.frame(pca_input)

  # Obtener features únicos
  feat <- pca_input[!duplicated(pca_input$FeatureID), , drop = FALSE]

  if (mode == "all") {
    return(feat$FeatureID)
  }

  if (mode == "any") {
    if (!("sig_any" %in% names(feat))) {
      stop("La columna 'sig_any' es requerida para mode = 'any'")
    }
    return(feat$FeatureID[feat$sig_any == TRUE])
  }

  # mode == "specific"
  if (is.null(comparison)) {
    stop("El argumento 'comparison' es requerido para mode = 'specific'")
  }

  col <- adjp_col(comparison)
  if (!(col %in% names(feat))) {
    stop("No existe la columna: ", col)
  }

  feat$FeatureID[feat[[col]] <= alpha]
}


# -----------------------------------------------------------------------------
# Función para construir scores de PCA
# -----------------------------------------------------------------------------

#' Construir dataframe de scores de PCA desde datos en formato long
#'
#' @param pca_input Data frame en formato long con columnas:
#'   - SampleID: Identificador de muestra
#'   - FeatureID: Identificador de proteína/feature
#'   - Intensity: Valor de intensidad (log2)
#'   - Condition: Condición experimental
#'   - Replicate: Número de réplica (opcional)
#'   - sig_any: Lógico indicando significancia en cualquier comparación (para mode="any")
#'   - adjP_*: Columnas de p-valores ajustados por comparación (para mode="specific")
#' @param mode Modo de filtrado de proteínas: "all", "any", o "specific"
#' @param alpha Umbral de significancia para modo "specific" (default: 0.05)
#' @param comparison Nombre de la comparación para modo "specific" (ej: "B-A")
#' @param subset_label Etiqueta personalizada para el subset (opcional)
#' @param center Centrar datos antes de PCA (default: TRUE)
#' @param scale. Escalar datos antes de PCA (default: TRUE)
#' @param filter_samples_to_comparison Filtrar muestras solo a las condiciones
#'   de la comparación específica (default: FALSE)
#' @param cond_col Nombre de la columna de condición (default: "Condition")
#'
#' @return Data frame con columnas: SampleID, PC1, PC2, PC1_Perc, PC2_Perc,
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

  # Etiqueta por defecto según el modo
  if (is.null(subset_label)) {
    subset_label <- switch(mode,
      all = "All proteins",
      any = "DEPs (any comparison)",
      specific = paste0("DEPs (", comparison, ")")
    )
  }

  # Convertir a data.frame para evitar problemas con tibbles
  pca_input <- as.data.frame(pca_input)

  # Metadata por muestra
  md_cols <- intersect(c("SampleID", "Condition", "Replicate"), names(pca_input))
  md <- pca_input[!duplicated(pca_input$SampleID), md_cols, drop = FALSE]
  rownames(md) <- md$SampleID

  # Filtrar muestras si se solicita (solo para mode = "specific")
  keep_samples <- md$SampleID
  if (isTRUE(filter_samples_to_comparison) && mode == "specific" && !is.null(comparison)) {
    conds <- unique(trimws(strsplit(comparison, "[-|:]")[[1]]))
    if (!(cond_col %in% names(md))) {
      stop("No existe la columna '", cond_col, "' en pca_input.")
    }
    keep_samples <- rownames(md)[md[[cond_col]] %in% conds]
    if (length(keep_samples) < 2) {
      stop("Menos de 2 muestras tras filtrar por comparación: ", comparison)
    }
  }

  # Obtener IDs de features según el modo
  ids <- get_feature_ids(pca_input, mode = mode, alpha = alpha, comparison = comparison)
  if (length(ids) < 2) {
    stop("Subset '", mode, "' sin suficientes proteínas para PCA (mínimo 2).")
  }

  # Filtrar datos
  dt <- pca_input[pca_input$SampleID %in% keep_samples & pca_input$FeatureID %in% ids,
                  c("SampleID", "FeatureID", "Intensity"), drop = FALSE]
  dt <- dt[is.finite(dt$Intensity) & !is.na(dt$Intensity), , drop = FALSE]

  # Pivotar a matriz (muestras x features)
  Xt <- with(dt, tapply(Intensity, list(SampleID, FeatureID), mean))
  Xt <- as.matrix(Xt)

  # Eliminar features con varianza 0
  v <- apply(Xt, 2, var, na.rm = TRUE)
  Xt <- Xt[, is.finite(v) & v > 0, drop = FALSE]

  if (ncol(Xt) < 2) {
    stop("Demasiado pocas proteínas con varianza > 0 para PCA en '", subset_label, "'.")
  }

  # Ejecutar PCA
  pc <- stats::prcomp(Xt, center = center, scale. = scale.)
  var_exp <- (pc$sdev^2) / sum(pc$sdev^2)
  scores <- pc$x[, 1:2, drop = FALSE]

  # Construir dataframe de salida
  out <- data.frame(
    SampleID = rownames(scores),
    PC1 = as.numeric(scores[, 1]),
    PC2 = as.numeric(scores[, 2]),
    PC1_Perc = round(100 * var_exp[1], 2),
    PC2_Perc = round(100 * var_exp[2], 2),
    Subset = subset_label,
    stringsAsFactors = FALSE
  )

  # Añadir metadata
  if ("Condition" %in% names(md)) {
    out$Condition <- md[out$SampleID, "Condition"]
  }
  if ("Replicate" %in% names(md)) {
    out$Replicate <- md[out$SampleID, "Replicate"]
  }

  out
}


# -----------------------------------------------------------------------------
# Función para calcular convex hull por grupo
# -----------------------------------------------------------------------------

#' Calcular convex hull (polígono envolvente) por grupo
#'
#' @param scores_df Data frame con columnas PC1, PC2 y la columna de grupo
#' @param group_col Nombre de la columna de agrupación (default: "Condition")
#'
#' @return Lista de data frames, cada uno con columnas: group, x, y
compute_hulls <- function(scores_df, group_col = "Condition") {

  # Convertir a data.frame para evitar problemas con tibbles
  scores_df <- as.data.frame(scores_df)

  required <- c("PC1", "PC2", group_col)
  missing <- setdiff(required, names(scores_df))
  if (length(missing) > 0) {
    stop("Columnas requeridas faltantes para hulls: ", paste(missing, collapse = ", "))
  }

  split_list <- split(scores_df, scores_df[[group_col]], drop = TRUE)

  hulls <- lapply(names(split_list), function(g) {
    d <- split_list[[g]]
    d <- d[is.finite(d$PC1) & is.finite(d$PC2), , drop = FALSE]

    # Se necesitan al menos 3 puntos para un hull
    if (nrow(d) < 3) return(NULL)

    # Calcular hull (índices del contorno)
    h <- grDevices::chull(d$PC1, d$PC2)
    # Cerrar el polígono repitiendo el primer punto
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
# Función para calcular elipse de confianza por grupo
# -----------------------------------------------------------------------------

#' Calcular elipse de confianza por grupo
#'
#' Calcula las coordenadas de una elipse de confianza basada en la distribución
#' chi-cuadrado, similar a FactoMineR::coord.ellipse y ggplot2::stat_ellipse.
#'
#' @param scores_df Data frame con columnas PC1, PC2 y la columna de grupo
#' @param group_col Nombre de la columna de agrupación (default: "Condition")
#' @param level Nivel de confianza (default: 0.95)
#' @param npoints Número de puntos para dibujar la elipse (default: 100)
#'
#' @return Lista de data frames, cada uno con columnas: group, x, y
compute_confidence_ellipse <- function(scores_df,
                                       group_col = "Condition",
                                       level = 0.95,
                                       npoints = 100) {

  # Convertir a data.frame
  scores_df <- as.data.frame(scores_df)

  required <- c("PC1", "PC2", group_col)
  missing <- setdiff(required, names(scores_df))
  if (length(missing) > 0) {
    stop("Columnas requeridas faltantes para ellipse: ", paste(missing, collapse = ", "))
  }

  split_list <- split(scores_df, scores_df[[group_col]], drop = TRUE)

  ellipses <- lapply(names(split_list), function(g) {
    d <- split_list[[g]]
    d <- d[is.finite(d$PC1) & is.finite(d$PC2), , drop = FALSE]

    # Se necesitan al menos 3 puntos para una elipse
    if (nrow(d) < 3) return(NULL)

    # Coordenadas
    x <- d$PC1
    y <- d$PC2

    # Centro (media)
    center_x <- mean(x)
    center_y <- mean(y)

    # Matriz de covarianza
    cov_mat <- cov(cbind(x, y))

    # Radio basado en distribución chi-cuadrado con 2 grados de libertad
    # Similar a FactoMineR: sqrt(qchisq(level, df = 2))
    radius <- sqrt(stats::qchisq(level, df = 2))

    # Descomposición eigen para obtener ejes de la elipse
    eigen_decomp <- eigen(cov_mat)
    eigenvalues <- eigen_decomp$values
    eigenvectors <- eigen_decomp$vectors

    # Verificar que los eigenvalues sean positivos
    if (any(eigenvalues <= 0)) return(NULL)

    # Ángulos para parametrizar la elipse
    theta <- seq(0, 2 * pi, length.out = npoints + 1)

    # Semi-ejes de la elipse
    a <- radius * sqrt(eigenvalues[1])
    b <- radius * sqrt(eigenvalues[2])

    # Ángulo de rotación
    angle <- atan2(eigenvectors[2, 1], eigenvectors[1, 1])

    # Coordenadas de la elipse (parametrización)
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
# Función principal: PCA Highchart
# -----------------------------------------------------------------------------

#' PCA Plot Interactivo con Highcharts
#'
#' Genera un scatter plot de PCA con opción de mostrar elipses o convex hulls
#' por grupo, similar a factoextra::fviz_pca_ind.
#'
#' @param scores_df Data frame generado por build_pca_scores() con columnas:
#'   PC1, PC2, PC1_Perc, PC2_Perc, Subset, Condition, Replicate, SampleID
#' @param color_by Columna para colorear los puntos (default: "Condition")
#' @param group_order Vector con el orden de los grupos/condiciones (opcional)
#' @param palette Vector nombrado de colores o NULL para paleta automática
#' @param title Título del gráfico (opcional, usa Subset por defecto)
#' @param addEllipses Mostrar elipses/hulls alrededor de los grupos (default: TRUE)
#' @param ellipse_type Tipo de elipse: "convex" para convex hull o "confidence"
#'   para elipse de confianza basada en distribución normal (default: "convex")
#' @param ellipse_level Nivel de confianza para ellipse_type = "confidence"
#'   (default: 0.95). Valores típicos: 0.95, 0.90, 0.68
#' @param ellipse_fill_opacity Opacidad del relleno (0-1, default: 0.12)
#' @param ellipse_line_width Ancho de línea del contorno (default: 1)
#' @param ellipse_npoints Número de puntos para dibujar la elipse de confianza
#'   (default: 100). Solo aplica para ellipse_type = "confidence"
#' @param point_size Radio de los puntos (default: 5)
#' @param show_labels Mostrar etiquetas de los puntos (SampleID) (default: FALSE)
#' @param label_size Tamaño de fuente de las etiquetas en px (default: 10)
#'
#' @return Objeto highchart
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
  # 1) Validación de inputs y conversión a data.frame
  # ---------------------------------------------------------------------------

  # Convertir a data.frame para evitar problemas con tibbles
  scores_df <- as.data.frame(scores_df)

  required <- c("PC1", "PC2", "SampleID")
  missing <- setdiff(required, names(scores_df))
  if (length(missing) > 0) {
    stop("Columnas requeridas faltantes: ", paste(missing, collapse = ", "))
  }

  if (!is.null(color_by) && !(color_by %in% names(scores_df))) {
    stop("La columna '", color_by, "' no existe en scores_df.")
  }

  # ---------------------------------------------------------------------------
  # 2) Configurar etiquetas de ejes con % varianza
  # ---------------------------------------------------------------------------
  pc1p <- unique(scores_df$PC1_Perc)
  pc2p <- unique(scores_df$PC2_Perc)
  x_lab <- if (length(pc1p) == 1) paste0("PC1 (", pc1p, "%)") else "PC1"
  y_lab <- if (length(pc2p) == 1) paste0("PC2 (", pc2p, "%)") else "PC2"

  # ---------------------------------------------------------------------------
  # 3) Configurar niveles de grupos
  # ---------------------------------------------------------------------------
  if (!is.null(group_order)) {
    scores_df[[color_by]] <- factor(scores_df[[color_by]], levels = group_order)
  } else {
    scores_df[[color_by]] <- factor(scores_df[[color_by]])
  }
  lvls <- levels(scores_df[[color_by]])

  # ---------------------------------------------------------------------------
  # 4) Configurar paleta de colores
  # ---------------------------------------------------------------------------
  default_palette <- c(
    "#457B9D", "#E63946", "#2A9D8F", "#E9C46A",
    "#9B5DE5", "#F4A261", "#264653", "#00BBF9"
  )

  if (is.null(palette)) {
    pal <- grDevices::hcl.colors(length(lvls), "Dark 3")
    palette <- stats::setNames(pal, lvls)
  } else if (is.character(palette) && is.null(names(palette))) {
    if (length(palette) < length(lvls)) {
      palette <- rep(palette, length.out = length(lvls))
    }
    palette <- stats::setNames(palette[seq_along(lvls)], lvls)
  } else if (is.character(palette) && !is.null(names(palette))) {
    miss <- setdiff(lvls, names(palette))
    if (length(miss) > 0) {
      stop("Faltan colores para niveles: ", paste(miss, collapse = ", "))
    }
    palette <- palette[lvls]
  }

  # ---------------------------------------------------------------------------
  # 5) Título del gráfico
  # ---------------------------------------------------------------------------
  chart_title <- title %||% unique(scores_df$Subset)[1]

  # ---------------------------------------------------------------------------
  # 6) Construir highchart base
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
          lineWidth = 1,
          lineColor = "#FFFFFF",
          states = list(
            hover = list(
              radiusPlus = 2,
              lineWidthPlus = 1,
              lineColor = "#1D3557"
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
  # 7) Añadir series de scatter por grupo (con id para vincular elipses)
  # ---------------------------------------------------------------------------
  split_list <- split(scores_df, scores_df[[color_by]], drop = TRUE)

  for (g in lvls) {
    if (!(g %in% names(split_list))) next

    d <- as.data.frame(split_list[[g]])
    group_id <- paste0("scatter_", gsub("[^a-zA-Z0-9]", "_", as.character(g)))

    # Crear lista de puntos con valores escalares explícitos
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

    # Configurar dataLabels si show_labels = TRUE
    data_labels_config <- if (isTRUE(show_labels)) {
      list(
        enabled = TRUE,
        format = "{point.SampleID}",
        style = list(
          fontSize = paste0(label_size, "px"),
          fontWeight = "normal",
          color = "#1D3557",
          textOutline = "2px #FFFFFF"
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
        color = unname(palette[as.character(g)]),
        zIndex = 5,
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
  # 8) Añadir elipses/hulls vinculadas a las series de scatter
  # ---------------------------------------------------------------------------
  if (isTRUE(addEllipses)) {

    # Calcular coordenadas según el tipo de elipse
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
      rgba_color <- hex_to_rgba(base_color, ellipse_fill_opacity)

      pts <- lapply(seq_len(nrow(poly)), function(k) {
        list(x = poly$x[k], y = poly$y[k])
      })

      hc <- hc |>
        hc_add_series(
          data = pts,
          type = "area",
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
# Función wrapper: Generar lista de PCA plots
# -----------------------------------------------------------------------------

#' Generar Lista de PCA Plots para Múltiples Subsets
#'
#' Genera automáticamente PCA plots para "all", "any" y/o comparaciones específicas.
#' Similar a factoextra::fviz_pca_ind con opciones de elipses.
#'
#' @param pca_input Data frame en formato long (ver build_pca_scores para estructura)
#' @param modes Vector de modos a generar: "all", "any", y/o nombres de comparaciones
#'   (default: c("all", "any"))
#' @param alpha Umbral de significancia para proteínas DEPs (default: 0.05)
#' @param color_by Columna para colorear (default: "Condition")
#' @param group_order Orden de grupos/condiciones (opcional)
#' @param palette Vector nombrado de colores o NULL para automático
#' @param addEllipses Mostrar elipses/hulls alrededor de los grupos (default: TRUE)
#' @param ellipse_type Tipo de elipse: "convex" o "confidence" (default: "convex")
#' @param ellipse_level Nivel de confianza para ellipse_type = "confidence"
#'   (default: 0.95)
#' @param ellipse_fill_opacity Opacidad del relleno de elipses (0-1, default: 0.12)
#' @param ellipse_line_width Ancho de línea del contorno (default: 1)
#' @param ellipse_npoints Número de puntos para elipse de confianza (default: 100)
#' @param point_size Radio de los puntos (default: 5)
#' @param show_labels Mostrar etiquetas de los puntos (SampleID) (default: FALSE)
#' @param label_size Tamaño de fuente de las etiquetas en px (default: 10)
#' @param center Centrar datos antes de PCA (default: TRUE)
#' @param scale. Escalar datos antes de PCA (default: TRUE)
#' @param filter_samples_to_comparison Para comparaciones específicas, filtrar
#'   muestras solo a las condiciones involucradas (default: FALSE)
#'
#' @return Lista nombrada de objetos highchart
#'
#' @examples
#' \dontrun{
#' # Cargar datos
#' pca_input <- arrow::read_parquet("PCA_Input.parquet")
#'
#' # PCA con convex hull (default)
#' hc_pcas <- pca_highchart_list(
#'   pca_input   = pca_input,
#'   modes       = c("all", "any"),
#'   group_order = c("A", "B", "C", "D")
#' )
#'
#' # PCA con elipse de confianza 95%
#' hc_pcas <- pca_highchart_list(
#'   pca_input     = pca_input,
#'   modes         = c("all"),
#'   group_order   = c("A", "B", "C", "D"),
#'   ellipse_type  = "confidence",
#'   ellipse_level = 0.95
#' )
#'
#' # PCA con etiquetas visibles
#' hc_pcas <- pca_highchart_list(
#'   pca_input   = pca_input,
#'   modes       = c("all"),
#'   show_labels = TRUE,
#'   label_size  = 9
#' )
#'
#' # Visualizar
#' hc_pcas[["all"]]
#' }
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
  # 1) Validación de inputs
  # ---------------------------------------------------------------------------
  required_cols <- c("SampleID", "FeatureID", "Intensity")
  missing_cols <- setdiff(required_cols, names(pca_input))
  if (length(missing_cols) > 0) {
    stop("Columnas requeridas faltantes: ", paste(missing_cols, collapse = ", "))
  }

  if (!is.null(color_by) && !(color_by %in% names(pca_input))) {
    warning("'", color_by, "' no está en el data frame. Se ignorará.")
    color_by <- NULL
  }

  # ---------------------------------------------------------------------------
  # 2) Detectar comparaciones disponibles (columnas adjP_*)
  # ---------------------------------------------------------------------------
  adjp_cols <- grep("^adjP_", names(pca_input), value = TRUE)
  available_comparisons <- sub("^adjP_", "", adjp_cols)

  # ---------------------------------------------------------------------------
  # 3) Generar plots para cada modo
  # ---------------------------------------------------------------------------
  hc_list <- list()

  for (m in modes) {

    # Determinar el mode interno y la comparación (si aplica)
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
      # Es una comparación específica
      if (!(m %in% available_comparisons)) {
        warning("Comparación '", m, "' no encontrada. Se omite.")
        next
      }
      internal_mode <- "specific"
      comparison <- m
      subset_label <- paste0("DEPs (", m, ")")
      plot_title <- paste0("PCA (DEPs ", m, ")")
    }

    # Intentar construir scores (puede fallar si no hay suficientes proteínas)
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
      warning("Error generando PCA para '", m, "': ", e$message)
      return(NULL)
    })

    if (is.null(scores_df)) next

    # Generar plot
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
# EJEMPLOS DE USO
# =============================================================================

# --- Cargar datos ---
# pca_input <- arrow::read_parquet("PCA_Input.parquet")
# pca_input <- readr::read_tsv("PCA_Input.tsv")

# --- Ejemplo básico: PCA con convex hull (default) ---
# sc_all <- build_pca_scores(pca_input, mode = "all")
# p_all <- pca_highchart(sc_all, color_by = "Condition", group_order = c("A","B","C","D"))
# p_all

# --- PCA con elipse de confianza 95% ---
# sc_all <- build_pca_scores(pca_input, mode = "all")
# p_all <- pca_highchart(
#   sc_all,
#   color_by = "Condition",
#   group_order = c("A","B","C","D"),
#   ellipse_type = "confidence",
#   ellipse_level = 0.95
# )
# p_all

# --- PCA con elipse de confianza 68% (1 desviación estándar) ---
# p_68 <- pca_highchart(
#   sc_all,
#   color_by = "Condition",
#   ellipse_type = "confidence",
#   ellipse_level = 0.68,
#   ellipse_fill_opacity = 0.2
# )
# p_68

# --- PCA sin elipses (solo puntos) ---
# p_noellipse <- pca_highchart(
#   sc_all,
#   color_by = "Condition",
#   addEllipses = FALSE,
#   point_size = 6
# )
# p_noellipse

# --- Generar múltiples PCA plots con convex hull ---
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

# --- Generar PCA plots con elipse de confianza ---
# hc_pcas <- pca_highchart_list(
#   pca_input     = pca_input,
#   modes         = c("all", "any"),
#   group_order   = c("A", "B", "C", "D"),
#   ellipse_type  = "confidence",
#   ellipse_level = 0.95,
#   ellipse_fill_opacity = 0.15
# )
# hc_pcas[["all"]]

# --- Con paleta personalizada ---
# my_palette <- c(A = "#457B9D", B = "#E63946", C = "#2A9D8F", D = "#E9C46A")
# hc_pcas <- pca_highchart_list(
#   pca_input   = pca_input,
#   modes       = c("all", "any"),
#   group_order = c("A", "B", "C", "D"),
#   palette     = my_palette
# )

# --- Comparar diferentes niveles de confianza ---
# sc_all <- build_pca_scores(pca_input, mode = "all")
#
# # 68% (1 SD)
# p_68 <- pca_highchart(sc_all, ellipse_type = "confidence", ellipse_level = 0.68,
#                       title = "PCA - 68% CI")
# # 95% (2 SD aprox)
# p_95 <- pca_highchart(sc_all, ellipse_type = "confidence", ellipse_level = 0.95,
#                       title = "PCA - 95% CI")
# # 99%
# p_99 <- pca_highchart(sc_all, ellipse_type = "confidence", ellipse_level = 0.99,
#                       title = "PCA - 99% CI")

# --- PCA con etiquetas visibles (útil para exportar) ---
# p_labels <- pca_highchart(
#   sc_all,
#   color_by = "Condition",
#   group_order = c("A","B","C","D"),
#   show_labels = TRUE,
#   label_size = 9
# )
# p_labels

# --- PCA con etiquetas y elipse de confianza ---
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

# --- Usando pca_highchart_list con etiquetas ---
# hc_pcas <- pca_highchart_list(
#   pca_input   = pca_input,
#   modes       = c("all", "any"),
#   group_order = c("A", "B", "C", "D"),
#   show_labels = TRUE,
#   label_size  = 9
# )
# hc_pcas[["all"]]
