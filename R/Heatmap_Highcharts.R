# =============================================================================
# Heatmap con tidyHeatmap para Datos de Proteómica
# =============================================================================

library(tidyHeatmap)
library(dplyr)
library(tidyr)
library(tibble)


# -----------------------------------------------------------------------------
# Operador null-coalesce
# -----------------------------------------------------------------------------

`%||%` <- function(a, b) if (!is.null(a) && length(a) && !is.na(a[1])) a else b


# -----------------------------------------------------------------------------
# Función para normalizar color hex (eliminar canal alpha si existe)
# -----------------------------------------------------------------------------

normalize_hex <- function(hex) {
  hex <- gsub("^#", "", hex)
  if (nchar(hex) == 8) {
    hex <- substr(hex, 1, 6)
  }
  paste0("#", hex)
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
#' @param data Data frame en formato long con columnas FeatureID, sig_any, adjP_*
#' @param mode Modo de filtrado: "all", "any", o "target"
#' @param alpha Umbral de significancia para modo "target" (default: 0.05)
#' @param comparison Nombre de la comparación para modo "target"
#'
#' @return Vector de FeatureIDs que cumplen el criterio
get_feature_ids <- function(data,
                            mode = c("all", "any", "target"),
                            alpha = 0.05,
                            comparison = NULL) {

  mode <- match.arg(mode)
  data <- as.data.frame(data)

  # Obtener features únicos

  feat <- data[!duplicated(data$FeatureID), , drop = FALSE]

  if (mode == "all") {
    return(feat$FeatureID)
  }

  if (mode == "any") {
    if (!("sig_any" %in% names(feat))) {
      stop("La columna 'sig_any' es requerida para mode = 'any'")
    }
    return(feat$FeatureID[feat$sig_any == TRUE])
  }

  # mode == "target"
  if (is.null(comparison)) {
    stop("El argumento 'comparison' es requerido para mode = 'target'")
  }

  col <- adjp_col(comparison)
  if (!(col %in% names(feat))) {
    stop("No existe la columna: ", col)
  }

  feat$FeatureID[feat[[col]] <= alpha]
}


# -----------------------------------------------------------------------------
# Función para obtener paleta de colores para el heatmap
# -----------------------------------------------------------------------------

#' Obtener paleta de colores para valores del heatmap
#'
#' @param palette Especificación de paleta:
#'   - NULL: usa paleta por defecto (RdBu)
#'   - Vector de colores: usa esos colores directamente
#'   - String "paletteer::" (ej: "viridis::viridis"): usa paletteer
#'   - String "brewer:" (ej: "brewer:RdBu"): usa RColorBrewer
#' @param n Número de colores a generar
#' @param reverse Invertir la paleta (default: FALSE)
#'
#' @return Vector de colores o función colorRamp2
get_heatmap_palette <- function(palette = NULL,
                                n = 11,
                                reverse = FALSE) {

  # Paleta por defecto (divergente azul-blanco-rojo)
  if (is.null(palette)) {
    colors <- c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7",
                "#FDDBC7", "#F4A582", "#D6604D", "#B2182B")
    if (reverse) colors <- rev(colors)
    return(colors)
  }

  # Paleta de paletteer
  if (is.character(palette) && length(palette) == 1 && grepl("::", palette)) {
    if (!requireNamespace("paletteer", quietly = TRUE)) {
      stop("Para usar paletteer, instala con: install.packages('paletteer')")
    }

    # Intentar como paleta discreta primero
    colors <- tryCatch({
      raw_pal <- as.character(paletteer::paletteer_d(palette, n = n))
      sapply(raw_pal, normalize_hex, USE.NAMES = FALSE)
    }, error = function(e) {
      # Intentar como paleta continua
      tryCatch({
        raw_pal <- as.character(paletteer::paletteer_c(palette, n = n))
        sapply(raw_pal, normalize_hex, USE.NAMES = FALSE)
      }, error = function(e2) {
        stop("Error al cargar paleta '", palette, "': ", e2$message)
      })
    })

    if (reverse) colors <- rev(colors)
    return(colors)
  }

  # Paleta de RColorBrewer
  if (is.character(palette) && length(palette) == 1 && startsWith(palette, "brewer:")) {
    if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
      stop("Para usar RColorBrewer, instala con: install.packages('RColorBrewer')")
    }

    nm <- sub("^brewer:", "", palette)
    colors <- tryCatch({
      maxc <- RColorBrewer::brewer.pal.info[nm, "maxcolors"]
      RColorBrewer::brewer.pal(min(n, maxc), nm)
    }, error = function(e) {
      stop("Error al cargar paleta brewer '", nm, "': ", e$message)
    })

    if (reverse) colors <- rev(colors)
    return(colors)
  }

  # Vector de colores directo
  if (is.character(palette) && length(palette) > 1) {
    if (reverse) palette <- rev(palette)
    return(palette)
  }

  # Nombre de paleta simple (intentar RColorBrewer)
  if (is.character(palette) && length(palette) == 1) {
    if (requireNamespace("RColorBrewer", quietly = TRUE) &&
        palette %in% rownames(RColorBrewer::brewer.pal.info)) {
      maxc <- RColorBrewer::brewer.pal.info[palette, "maxcolors"]
      colors <- RColorBrewer::brewer.pal(min(n, maxc), palette)
      if (reverse) colors <- rev(colors)
      return(colors)
    }
  }

  # Por defecto
  colors <- c("#2166AC", "#4393C3", "#92C5DE", "#D1E5F0", "#F7F7F7",
              "#FDDBC7", "#F4A582", "#D6604D", "#B2182B")
  if (reverse) colors <- rev(colors)
  colors
}


#' Obtener paleta de colores para anotaciones categóricas
#'
#' @param levels Vector de niveles (categorías)
#' @param palette Especificación de paleta (mismo formato que get_heatmap_palette)
#'
#' @return Vector nombrado de colores
get_annotation_palette <- function(levels, palette = NULL) {

  n <- length(levels)

  # Paleta por defecto
  default_palette <- c(
    "#457B9D", "#E63946", "#2A9D8F", "#E9C46A",
    "#9B5DE5", "#F4A261", "#264653", "#00BBF9",
    "#FF6B6B", "#4ECDC4", "#45B7D1", "#96CEB4"
  )

  if (is.null(palette)) {
    pal <- grDevices::hcl.colors(n, "Dark 3")
    return(stats::setNames(pal, levels))
  }

  # Paleta de paletteer
  if (is.character(palette) && length(palette) == 1 && grepl("::", palette)) {
    if (!requireNamespace("paletteer", quietly = TRUE)) {
      stop("Para usar paletteer, instala con: install.packages('paletteer')")
    }
    colors <- tryCatch({
      raw_pal <- as.character(paletteer::paletteer_d(palette))
      sapply(raw_pal, normalize_hex, USE.NAMES = FALSE)
    }, error = function(e) {
      stop("Error al cargar paleta '", palette, "': ", e$message)
    })
    if (length(colors) < n) {
      colors <- rep(colors, length.out = n)
    }
    return(stats::setNames(colors[seq_len(n)], levels))
  }

  # Paleta de RColorBrewer
  if (is.character(palette) && length(palette) == 1 && startsWith(palette, "brewer:")) {
    if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
      stop("Para usar RColorBrewer, instala con: install.packages('RColorBrewer')")
    }
    nm <- sub("^brewer:", "", palette)
    colors <- tryCatch({
      maxc <- RColorBrewer::brewer.pal.info[nm, "maxcolors"]
      RColorBrewer::brewer.pal(max(3, min(n, maxc)), nm)
    }, error = function(e) {
      stop("Error al cargar paleta brewer '", nm, "': ", e$message)
    })
    if (length(colors) < n) {
      colors <- rep(colors, length.out = n)
    }
    return(stats::setNames(colors[seq_len(n)], levels))
  }

  # Vector nombrado de colores
  if (is.character(palette) && !is.null(names(palette))) {
    miss <- setdiff(levels, names(palette))
    if (length(miss) > 0) {
      stop("Faltan colores para niveles: ", paste(miss, collapse = ", "))
    }
    return(palette[levels])
  }

  # Vector de colores sin nombres
  if (is.character(palette) && length(palette) >= 1) {
    if (length(palette) < n) {
      palette <- rep(palette, length.out = n)
    }
    return(stats::setNames(palette[seq_len(n)], levels))
  }

  # Por defecto
  if (length(default_palette) < n) {
    default_palette <- rep(default_palette, length.out = n)
  }
  stats::setNames(default_palette[seq_len(n)], levels)
}


# -----------------------------------------------------------------------------
# Función principal: Preparar datos para heatmap
# -----------------------------------------------------------------------------

#' Preparar datos en formato largo para tidyHeatmap
#'
#' @param data Data frame en formato largo con columnas:
#'   - SampleID: Identificador de muestra
#'   - FeatureID: Identificador de proteína/feature
#'   - Intensity: Valor de intensidad (log2)
#'   - Condition: Condición experimental
#'   - Replicate: Número de réplica
#'   - sig_any: Lógico indicando significancia (para mode="any")
#'   - adjP_*: Columnas de p-valores ajustados (para mode="target")
#' @param mode Modo de filtrado: "all", "any", o "target"
#' @param alpha Umbral de significancia para modo "target" (default: 0.05)
#' @param comparison Nombre de la comparación para modo "target" (ej: "B-A")
#' @param scale_data Tipo de escalado: "none", "row", "column" (default: "row")
#' @param sample_order Orden de muestras: "clustering", "condition", o vector personalizado
#' @param condition_order Orden de condiciones cuando sample_order = "condition"
#'
#' @return Data frame en formato largo listo para tidyHeatmap
prepare_heatmap_data <- function(data,
                                 mode = c("all", "any", "target"),
                                 alpha = 0.05,
                                 comparison = NULL,
                                 scale_data = c("row", "none", "column"),
                                 sample_order = c("clustering", "condition", "custom"),
                                 condition_order = NULL) {

  mode <- match.arg(mode)
  scale_data <- match.arg(scale_data)

  if (is.character(sample_order) && length(sample_order) == 1) {
    sample_order <- match.arg(sample_order, c("clustering", "condition", "custom"))
  }

  # Validar columnas requeridas
  required_cols <- c("SampleID", "FeatureID", "Intensity", "Condition", "Replicate")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Columnas requeridas faltantes: ", paste(missing_cols, collapse = ", "))
  }

  data <- as.data.frame(data)

  # Obtener IDs de features según el modo
  ids <- get_feature_ids(data, mode = mode, alpha = alpha, comparison = comparison)
  if (length(ids) < 2) {
    stop("Subset '", mode, "' sin suficientes proteínas (mínimo 2). Encontradas: ", length(ids))
  }

  # Filtrar datos
  dt <- data[data$FeatureID %in% ids, , drop = FALSE]
  dt <- dt[is.finite(dt$Intensity) & !is.na(dt$Intensity), , drop = FALSE]

  # Crear data frame largo para heatmap
  hm_data <- dt %>%
    select(SampleID, FeatureID, Intensity, Condition, Replicate) %>%
    distinct()

  # Ordenar muestras según el criterio especificado
  if (is.character(sample_order) && sample_order == "condition") {
    # Ordenar por Condition y luego por Replicate
    if (!is.null(condition_order)) {
      hm_data <- hm_data %>%
        mutate(Condition = factor(Condition, levels = condition_order))
    }
    hm_data <- hm_data %>%
      arrange(Condition, Replicate) %>%
      mutate(SampleID = factor(SampleID, levels = unique(SampleID)))
  } else if (is.character(sample_order) && length(sample_order) > 1) {
    # Orden personalizado (vector de SampleIDs)
    if (!all(sample_order %in% unique(hm_data$SampleID))) {
      missing <- setdiff(sample_order, unique(hm_data$SampleID))
      warning("Muestras no encontradas en datos: ", paste(missing, collapse = ", "))
    }
    hm_data <- hm_data %>%
      mutate(SampleID = factor(SampleID, levels = sample_order))
  }
  # Si sample_order == "clustering", dejamos que tidyHeatmap haga el clustering

  # Aplicar escalado si se solicita
  if (scale_data != "none") {
    # Pivotar a matriz para escalar
    mat_wide <- hm_data %>%
      select(SampleID, FeatureID, Intensity) %>%
      pivot_wider(names_from = SampleID, values_from = Intensity) %>%
      as.data.frame()

    rownames_feat <- mat_wide$FeatureID
    mat <- as.matrix(mat_wide[, -1])
    rownames(mat) <- rownames_feat

    if (scale_data == "row") {
      mat <- t(scale(t(mat)))
    } else if (scale_data == "column") {
      mat <- scale(mat)
    }

    # Reemplazar NaN por NA
    mat[is.nan(mat)] <- NA

    # Volver a formato largo
    mat_df <- as.data.frame(mat)
    mat_df$FeatureID <- rownames(mat)

    scaled_long <- mat_df %>%
      pivot_longer(cols = -FeatureID, names_to = "SampleID", values_to = "Intensity")

    # Reintegrar metadata
    metadata <- hm_data %>%
      select(SampleID, Condition, Replicate) %>%
      distinct()

    hm_data <- scaled_long %>%
      left_join(metadata, by = "SampleID")

    # Restaurar orden de factores si aplica
    if (is.character(sample_order) && sample_order == "condition") {
      if (!is.null(condition_order)) {
        hm_data <- hm_data %>%
          mutate(Condition = factor(Condition, levels = condition_order))
      }
      hm_data <- hm_data %>%
        arrange(Condition, as.numeric(Replicate)) %>%
        mutate(SampleID = factor(SampleID, levels = unique(SampleID)))
    }
  }

  hm_data
}


# -----------------------------------------------------------------------------
# Función principal: Heatmap con tidyHeatmap
# -----------------------------------------------------------------------------

#' Crear Heatmap con tidyHeatmap para Proteómica
#'
#' @param data Data frame en formato largo (ver prepare_heatmap_data para estructura)
#' @param mode Modo de filtrado de proteínas: "all", "any", o "target"
#' @param alpha Umbral de significancia para proteínas DEPs (default: 0.05)
#' @param comparison Nombre de la comparación para modo "target" (ej: "B-A")
#' @param scale_data Tipo de escalado: "none", "row", "column" (default: "row")
#' @param sample_order Orden de muestras: "clustering", "condition", o vector personalizado
#' @param condition_order Orden de condiciones cuando sample_order = "condition"
#' @param cluster_rows Aplicar clustering a filas (default: FALSE)
#' @param cluster_columns Aplicar clustering a columnas (default: TRUE)
#' @param show_row_names Mostrar nombres de filas (default: TRUE si <= 50 proteínas)
#' @param show_column_names Mostrar nombres de columnas (default: TRUE)
#' @param palette_value Paleta para valores del heatmap (ver get_heatmap_palette)
#' @param palette_annotation Paleta para anotación de Condition
#' @param reverse_palette Invertir paleta de valores (default: FALSE)
#' @param row_title Título para las filas
#' @param column_title Título para las columnas
#' @param show_annotation Mostrar anotación de Condition (default: TRUE)
#' @param split_by_condition Dividir el heatmap por condición (default: FALSE)
#' @param row_names_size Tamaño de fuente de nombres de fila (default: 7)
#' @param column_names_size Tamaño de fuente de nombres de columna (default: 9)
#' @param column_names_rotation Rotación de nombres de columna en grados (default: 45)
#'
#' @return Objeto tidyHeatmap/ComplexHeatmap
#'
#' @examples
#' \dontrun{
#' # Cargar datos
#' hm_input <- arrow::read_parquet("PCA_Input.parquet")
#'
#' # Heatmap de todas las proteínas
#' hm_all <- proteomics_heatmap(
#'   data = hm_input,
#'   mode = "all",
#'   scale_data = "row",
#'   sample_order = "condition",
#'   condition_order = c("A", "B", "C", "D")
#' )
#'
#' # Heatmap con separación por condición
#' hm_split <- proteomics_heatmap(
#'   data = hm_input,
#'   mode = "any",
#'   scale_data = "row",
#'   split_by_condition = TRUE,
#'   condition_order = c("A", "B", "C", "D")
#' )
#'
#' # Heatmap de proteínas significativas (any)
#' hm_any <- proteomics_heatmap(
#'   data = hm_input,
#'   mode = "any",
#'   scale_data = "row"
#' )
#'
#' # Heatmap de proteínas significativas en comparación específica
#' hm_target <- proteomics_heatmap(
#'   data = hm_input,
#'   mode = "target",
#'   comparison = "B-A",
#'   scale_data = "row"
#' )
#' }
proteomics_heatmap <- function(data,
                               mode = c("all", "any", "target"),
                               alpha = 0.05,
                               comparison = NULL,
                               scale_data = c("row", "none", "column"),
                               sample_order = "clustering",
                               condition_order = NULL,
                               cluster_rows = FALSE,
                               cluster_columns = TRUE,
                               show_row_names = NULL,
                               show_column_names = TRUE,
                               palette_value = NULL,
                               palette_annotation = NULL,
                               reverse_palette = FALSE,
                               row_title = "Proteins",
                               column_title = "Samples",
                               show_annotation = TRUE,
                               split_by_condition = FALSE,
                               row_names_size = 7,
                               column_names_size = 9,
                               column_names_rotation = 45) {

  mode <- match.arg(mode)
  scale_data <- match.arg(scale_data)

  # ---------------------------------------------------------------------------
  # 1) Preparar datos
  # ---------------------------------------------------------------------------

  hm_data <- prepare_heatmap_data(
    data = data,
    mode = mode,
    alpha = alpha,
    comparison = comparison,
    scale_data = scale_data,
    sample_order = sample_order,
    condition_order = condition_order
  )

  # Determinar número de proteínas para auto-configuración
  n_proteins <- length(unique(hm_data$FeatureID))
  n_samples <- length(unique(hm_data$SampleID))

  # Auto-decidir si mostrar nombres de filas
  if (is.null(show_row_names)) {
    show_row_names <- n_proteins <= 50
  }

  # ---------------------------------------------------------------------------
  # 2) Obtener paletas de colores
  # ---------------------------------------------------------------------------

  # Paleta para valores
  value_colors <- get_heatmap_palette(
    palette = palette_value,
    n = 11,
    reverse = reverse_palette
  )

  # Crear función colorRamp2 para valores
  if (requireNamespace("circlize", quietly = TRUE)) {
    # Determinar rango de valores
    val_range <- range(hm_data$Intensity, na.rm = TRUE)

    if (scale_data != "none") {
      # Para datos escalados, usar rango simétrico
      abs_max <- max(abs(val_range), na.rm = TRUE)
      val_breaks <- seq(-abs_max, abs_max, length.out = length(value_colors))
    } else {
      val_breaks <- seq(val_range[1], val_range[2], length.out = length(value_colors))
    }

    palette_func <- circlize::colorRamp2(val_breaks, value_colors)
  } else {
    palette_func <- value_colors
  }

  # Paleta para anotación de Condition
  condition_levels <- unique(hm_data$Condition)
  if (!is.null(condition_order)) {
    condition_levels <- condition_order[condition_order %in% condition_levels]
  }
  annotation_colors <- get_annotation_palette(condition_levels, palette_annotation)

  # ---------------------------------------------------------------------------
  # 3) Configurar clustering
  # ---------------------------------------------------------------------------

  # Determinar si aplicar clustering a columnas
  cluster_cols_final <- cluster_columns
  if (is.character(sample_order) && sample_order != "clustering") {
    cluster_cols_final <- FALSE
  }

  # ---------------------------------------------------------------------------
  # 4) Crear heatmap con tidyHeatmap
  # ---------------------------------------------------------------------------

  # Convertir a tibble (requerido por tidyHeatmap)
  hm_data <- tibble::as_tibble(hm_data)

  # Crear heatmap base
  hm <- hm_data %>%
    heatmap(
      .row = FeatureID,
      .column = SampleID,
      .value = Intensity,
      scale = "none",  # Ya escalamos antes
      cluster_rows = cluster_rows,
      cluster_columns = cluster_cols_final,
      palette_value = palette_func,
      show_row_names = show_row_names,
      show_column_names = show_column_names,
      row_names_gp = grid::gpar(fontsize = row_names_size),
      column_names_gp = grid::gpar(fontsize = column_names_size),
      column_names_rot = column_names_rotation,
      row_title = row_title,
      column_title = column_title
    )

  # Añadir separación por condición si se solicita (usando annotation_group)
  if (split_by_condition) {
    hm <- hm %>%
      annotation_group(Condition)
  }

  # Añadir anotación de Condition como barra de color si se solicita
  if (show_annotation) {
    hm <- hm %>%
      annotation_tile(
        Condition,
        palette = annotation_colors
      )
  }

  hm
}


# -----------------------------------------------------------------------------
# Función wrapper: Generar lista de heatmaps
# -----------------------------------------------------------------------------

#' Generar Lista de Heatmaps para Múltiples Subsets
#'
#' Genera automáticamente heatmaps para "all", "any" y/o comparaciones específicas.
#'
#' @param data Data frame en formato long (ver prepare_heatmap_data para estructura)
#' @param modes Vector de modos a generar: "all", "any", y/o nombres de comparaciones
#'   (default: c("all", "any"))
#' @param alpha Umbral de significancia para proteínas DEPs (default: 0.05)
#' @param scale_data Tipo de escalado: "none", "row", "column" (default: "row")
#' @param sample_order Orden de muestras: "clustering", "condition", o vector personalizado
#' @param condition_order Orden de condiciones cuando sample_order = "condition"
#' @param cluster_rows Aplicar clustering a filas (default: FALSE)
#' @param cluster_columns Aplicar clustering a columnas (default: TRUE)
#' @param show_row_names Mostrar nombres de filas (default: auto)
#' @param show_column_names Mostrar nombres de columnas (default: TRUE)
#' @param palette_value Paleta para valores del heatmap
#' @param palette_annotation Paleta para anotación de Condition
#' @param reverse_palette Invertir paleta de valores (default: FALSE)
#' @param show_annotation Mostrar anotación de Condition (default: TRUE)
#' @param split_by_condition Dividir el heatmap por condición (default: FALSE)
#' @param row_names_size Tamaño de fuente de nombres de fila (default: 7)
#' @param column_names_size Tamaño de fuente de nombres de columna (default: 9)
#' @param column_names_rotation Rotación de nombres de columna (default: 45)
#'
#' @return Lista nombrada de objetos tidyHeatmap
#'
#' @examples
#' \dontrun{
#' # Cargar datos
#' hm_input <- arrow::read_parquet("PCA_Input.parquet")
#'
#' # Generar heatmaps para all y any
#' hm_list <- proteomics_heatmap_list(
#'   data = hm_input,
#'   modes = c("all", "any"),
#'   sample_order = "condition",
#'   condition_order = c("A", "B", "C", "D")
#' )
#'
#' # Visualizar
#' hm_list[["all"]]
#' hm_list[["any"]]
#'
#' # Generar heatmaps incluyendo comparaciones específicas
#' hm_list <- proteomics_heatmap_list(
#'   data = hm_input,
#'   modes = c("all", "any", "B-A", "C-A"),
#'   sample_order = "clustering"
#' )
#' hm_list[["B-A"]]
#'
#' # Con paleta personalizada de paletteer
#' hm_list <- proteomics_heatmap_list(
#'   data = hm_input,
#'   modes = c("all"),
#'   palette_value = "viridis::viridis"
#' )
#'
#' # Con paleta de RColorBrewer
#' hm_list <- proteomics_heatmap_list(
#'   data = hm_input,
#'   modes = c("all"),
#'   palette_value = "brewer:RdYlBu",
#'   palette_annotation = "brewer:Set1"
#' )
#' }
proteomics_heatmap_list <- function(data,
                                    modes = c("all", "any"),
                                    alpha = 0.05,
                                    scale_data = c("row", "none", "column"),
                                    sample_order = "clustering",
                                    condition_order = NULL,
                                    cluster_rows = FALSE,
                                    cluster_columns = TRUE,
                                    show_row_names = NULL,
                                    show_column_names = TRUE,
                                    palette_value = NULL,
                                    palette_annotation = NULL,
                                    reverse_palette = FALSE,
                                    show_annotation = TRUE,
                                    split_by_condition = FALSE,
                                    row_names_size = 7,
                                    column_names_size = 9,
                                    column_names_rotation = 45) {

  scale_data <- match.arg(scale_data)

  # ---------------------------------------------------------------------------
  # 1) Validación de inputs
  # ---------------------------------------------------------------------------

  required_cols <- c("SampleID", "FeatureID", "Intensity", "Condition", "Replicate")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Columnas requeridas faltantes: ", paste(missing_cols, collapse = ", "))
  }

  # ---------------------------------------------------------------------------
  # 2) Detectar comparaciones disponibles (columnas adjP_*)
  # ---------------------------------------------------------------------------

  adjp_cols <- grep("^adjP_", names(data), value = TRUE)
  available_comparisons <- sub("^adjP_", "", adjp_cols)

  # ---------------------------------------------------------------------------
  # 3) Generar heatmaps para cada modo
  # ---------------------------------------------------------------------------

  hm_list <- list()

  for (m in modes) {

    # Determinar el mode interno y la comparación (si aplica)
    if (m == "all") {
      internal_mode <- "all"
      comparison <- NULL
      row_title <- "All Proteins"
    } else if (m == "any") {
      internal_mode <- "any"
      comparison <- NULL
      row_title <- "DEPs (any comparison)"
    } else {
      # Es una comparación específica (modo "target")
      if (!(m %in% available_comparisons)) {
        warning("Comparación '", m, "' no encontrada. Se omite.")
        next
      }
      internal_mode <- "target"
      comparison <- m
      row_title <- paste0("DEPs (", m, ")")
    }

    # Intentar generar heatmap (puede fallar si no hay suficientes proteínas)
    hm <- tryCatch({
      proteomics_heatmap(
        data = data,
        mode = internal_mode,
        alpha = alpha,
        comparison = comparison,
        scale_data = scale_data,
        sample_order = sample_order,
        condition_order = condition_order,
        cluster_rows = cluster_rows,
        cluster_columns = cluster_columns,
        show_row_names = show_row_names,
        show_column_names = show_column_names,
        palette_value = palette_value,
        palette_annotation = palette_annotation,
        reverse_palette = reverse_palette,
        row_title = row_title,
        column_title = "Samples",
        show_annotation = show_annotation,
        split_by_condition = split_by_condition,
        row_names_size = row_names_size,
        column_names_size = column_names_size,
        column_names_rotation = column_names_rotation
      )
    }, error = function(e) {
      warning("Error generando heatmap para '", m, "': ", e$message)
      return(NULL)
    })

    if (!is.null(hm)) {
      hm_list[[m]] <- hm
    }
  }

  hm_list
}


# =============================================================================
# EJEMPLOS DE USO
# =============================================================================

# --- Cargar datos ---
# hm_input <- arrow::read_parquet("PCA_Input.parquet")
# hm_input <- readr::read_tsv("PCA_Input.tsv")

# --- Ejemplo básico: Heatmap de todas las proteínas ---
# hm_all <- proteomics_heatmap(
#   data = hm_input,
#   mode = "all",
#   scale_data = "row",
#   sample_order = "condition",
#   condition_order = c("A", "B", "C", "D")
# )
# hm_all

# --- Heatmap con clustering de columnas ---
# hm_cluster <- proteomics_heatmap(
#   data = hm_input,
#   mode = "all",
#   scale_data = "row",
#   sample_order = "clustering",
#   cluster_columns = TRUE,
#   cluster_rows = FALSE
# )
# hm_cluster

# --- Heatmap de proteínas significativas en cualquier comparación ---
# hm_any <- proteomics_heatmap(
#   data = hm_input,
#   mode = "any",
#   scale_data = "row",
#   sample_order = "condition",
#   condition_order = c("A", "B", "C", "D")
# )
# hm_any

# --- Heatmap de proteínas significativas en una comparación específica ---
# hm_target <- proteomics_heatmap(
#   data = hm_input,
#   mode = "target",
#   comparison = "B-A",
#   alpha = 0.05,
#   scale_data = "row"
# )
# hm_target

# --- Con paleta personalizada de paletteer ---
# hm_viridis <- proteomics_heatmap(
#   data = hm_input,
#   mode = "all",
#   scale_data = "row",
#   palette_value = "viridis::viridis",
#   reverse_palette = TRUE
# )
# hm_viridis

# --- Con paleta de RColorBrewer ---
# hm_brewer <- proteomics_heatmap(
#   data = hm_input,
#   mode = "all",
#   scale_data = "row",
#   palette_value = "brewer:RdYlBu",
#   palette_annotation = "brewer:Set1"
# )
# hm_brewer

# --- Generar múltiples heatmaps ---
# hm_list <- proteomics_heatmap_list(
#   data = hm_input,
#   modes = c("all", "any", "B-A", "C-A", "D-A"),
#   scale_data = "row",
#   sample_order = "condition",
#   condition_order = c("A", "B", "C", "D")
# )
# hm_list[["all"]]
# hm_list[["any"]]
# hm_list[["B-A"]]

# --- Orden personalizado de muestras ---
# custom_order <- c("A_1", "A_2", "B_1", "B_2", "C_1", "C_2", "D_1", "D_2")
# hm_custom <- proteomics_heatmap(
#   data = hm_input,
#   mode = "all",
#   scale_data = "row",
#   sample_order = custom_order,
#   cluster_columns = FALSE
# )
# hm_custom

# --- Paletas disponibles ---
# Divergentes (buenas para datos escalados):
#   - NULL (default RdBu-like)
#   - "brewer:RdBu", "brewer:RdYlBu", "brewer:PiYG", "brewer:BrBG"
#   - "viridis::plasma", "viridis::inferno"
#
# Secuenciales (buenas para datos no escalados):
#   - "viridis::viridis", "viridis::magma"
#   - "brewer:Blues", "brewer:Reds", "brewer:YlOrRd"
#
# Para anotaciones:
#   - "brewer:Set1", "brewer:Set2", "brewer:Dark2"
#   - "ggsci::category10_d3", "ggsci::nrc_npg"
