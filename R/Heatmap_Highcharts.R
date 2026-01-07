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
# Wrapper y método print para heatmaps con título personalizado
# -----------------------------------------------------------------------------

#' Crear wrapper para heatmap con título
#'
#' @param hm Objeto tidyHeatmap (InputHeatmap)
#' @param title Título del heatmap
#' @param title_size Tamaño de fuente del título
#' @param title_face Estilo de fuente del título
#'
#' @return Objeto proteomics_heatmap (lista S3)
wrap_heatmap_with_title <- function(hm, title, title_size = 14, title_face = "bold") {

  structure(
    list(
      heatmap = hm,
      title = title,
      title_size = title_size,
      title_face = title_face
    ),
    class = "proteomics_heatmap"
  )
}

#' @export
print.proteomics_heatmap <- function(x, ...) {
  hm <- x$heatmap

  # Convertir InputHeatmap a ComplexHeatmap usando el método de tidyHeatmap
  if (inherits(hm, "InputHeatmap")) {
    # Usar as.list para extraer los componentes y luego reconstruir
    # O simplemente convertir usando el método interno de tidyHeatmap
    ht <- tryCatch({
      # Intentar slot directo (versiones antiguas)
      methods::slot(hm, "ht")
    }, error = function(e) {
      tryCatch({
        # Intentar con input_heatmap (versiones más nuevas)
        methods::slot(hm, "input_heatmap")
      }, error = function(e2) {
        # Usar la conversión de tidyHeatmap a ComplexHeatmap
        tidyHeatmap::as_ComplexHeatmap(hm)
      })
    })
  } else {
    ht <- hm
  }

  # Dibujar con título
  ComplexHeatmap::draw(
    ht,
    column_title = x$title,
    column_title_gp = grid::gpar(fontsize = x$title_size, fontface = x$title_face),
    ...
  )

  invisible(x)
}


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
#' @param feature_ids Vector de FeatureIDs específicos (si se proporciona, ignora mode/alpha)
#' @param scale_data Tipo de escalado: "none", "row", "column" (default: "row")
#' @param sample_order Orden de muestras: "clustering", "condition", o vector personalizado
#' @param condition_order Orden de condiciones cuando sample_order = "condition"
#'
#' @return Data frame en formato largo listo para tidyHeatmap
prepare_heatmap_data <- function(data,
                                 mode = c("all", "any", "target"),
                                 alpha = 0.05,
                                 comparison = NULL,
                                 feature_ids = NULL,
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

  # Para mode = "target", filtrar muestras solo a las condiciones de la comparación
  if (mode == "target" && !is.null(comparison)) {
    # Extraer condiciones de la comparación (ej: "B-A" -> c("B", "A"))
    conds <- unique(trimws(strsplit(comparison, "[-|:]")[[1]]))
    if (length(conds) >= 2) {
      data <- data[data$Condition %in% conds, , drop = FALSE]
      # Actualizar condition_order para solo incluir las condiciones relevantes
      if (!is.null(condition_order)) {
        condition_order <- condition_order[condition_order %in% conds]
      } else {
        condition_order <- conds
      }
    }
  }

  # Obtener IDs de features según el modo o usar los proporcionados
  if (!is.null(feature_ids) && length(feature_ids) > 0) {
    # Usar IDs proporcionados directamente
    available_ids <- unique(data$FeatureID)
    ids <- feature_ids[feature_ids %in% available_ids]
    if (length(ids) < length(feature_ids)) {
      missing <- setdiff(feature_ids, available_ids)
      warning("FeatureIDs no encontrados en datos (ignorados): ",
              paste(head(missing, 5), collapse = ", "),
              if (length(missing) > 5) paste0(" ... y ", length(missing) - 5, " más"))
    }
  } else {
    # Usar filtrado por mode
    ids <- get_feature_ids(data, mode = mode, alpha = alpha, comparison = comparison)
  }

  if (length(ids) < 2) {
    stop("Subset sin suficientes proteínas (mínimo 2). Encontradas: ", length(ids))
  }

  # Filtrar datos
  dt <- data[data$FeatureID %in% ids, , drop = FALSE]
  dt <- dt[is.finite(dt$Intensity) & !is.na(dt$Intensity), , drop = FALSE]


  # Crear data frame largo para heatmap
  # Incluir adjP si mode = "target" para anotación de filas
  if (mode == "target" && !is.null(comparison)) {
    adjp_colname <- adjp_col(comparison)
    if (adjp_colname %in% names(dt)) {
      hm_data <- dt %>%
        select(SampleID, FeatureID, Intensity, Condition, Replicate, all_of(adjp_colname)) %>%
        rename(adjP = all_of(adjp_colname)) %>%
        distinct()
    } else {
      hm_data <- dt %>%
        select(SampleID, FeatureID, Intensity, Condition, Replicate) %>%
        distinct()
    }
  } else {
    hm_data <- dt %>%
      select(SampleID, FeatureID, Intensity, Condition, Replicate) %>%
      distinct()
  }

  # Ordenar muestras según el criterio especificado
  if (is.character(sample_order) && length(sample_order) == 1 && sample_order == "condition") {
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
    # Filtrar para incluir solo las muestras especificadas
    samples_in_data <- unique(hm_data$SampleID)
    valid_samples <- sample_order[sample_order %in% samples_in_data]

    if (length(valid_samples) == 0) {
      stop("Ninguna de las muestras especificadas en sample_order está en los datos")
    }

    if (!all(sample_order %in% samples_in_data)) {
      missing <- setdiff(sample_order, samples_in_data)
      warning("Muestras no encontradas en datos (ignoradas): ", paste(missing, collapse = ", "))
    }

    # Filtrar datos para incluir solo las muestras especificadas
    hm_data <- hm_data %>%
      filter(SampleID %in% valid_samples) %>%
      mutate(SampleID = factor(SampleID, levels = valid_samples))
  }
  # Si sample_order == "clustering", dejamos que tidyHeatmap haga el clustering

  # Aplicar escalado si se solicita
  if (scale_data != "none") {
    # Pivotar a matriz para escalar (usar mean para posibles duplicados)
    mat_wide <- hm_data %>%
      select(SampleID, FeatureID, Intensity) %>%
      pivot_wider(names_from = SampleID, values_from = Intensity, values_fn = mean) %>%
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

    # Reintegrar metadata (incluyendo adjP si existe)
    if ("adjP" %in% names(hm_data)) {
      metadata <- hm_data %>%
        select(SampleID, FeatureID, Condition, Replicate, adjP) %>%
        distinct()

      hm_data <- scaled_long %>%
        left_join(metadata, by = c("SampleID", "FeatureID"))
    } else {
      metadata <- hm_data %>%
        select(SampleID, Condition, Replicate) %>%
        distinct()

      hm_data <- scaled_long %>%
        left_join(metadata, by = "SampleID")
    }

    # Restaurar orden de factores si aplica
    if (is.character(sample_order) && length(sample_order) == 1 && sample_order == "condition") {
      if (!is.null(condition_order)) {
        hm_data <- hm_data %>%
          mutate(Condition = factor(Condition, levels = condition_order))
      }
      hm_data <- hm_data %>%
        arrange(Condition, as.numeric(Replicate)) %>%
        mutate(SampleID = factor(SampleID, levels = unique(SampleID)))
    } else if (is.character(sample_order) && length(sample_order) > 1) {
      # Restaurar orden personalizado
      hm_data <- hm_data %>%
        mutate(SampleID = factor(SampleID, levels = sample_order))
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
#' @param feature_ids Vector de FeatureIDs específicos a mostrar (default: NULL, usa filtrado por mode).
#'   Si se proporciona, solo se muestran estas proteínas, ignorando el filtrado por mode/alpha.
#' @param row_annotation Anotaciones para filas (FeatureIDs). Puede ser:
#'   - Ruta a archivo TSV con columna "FeatureID" y columnas categóricas adicionales
#'   - Data frame con la misma estructura
#'   - NULL: sin anotaciones de fila (default)
#' @param row_annotation_cols Vector de nombres de columnas a mostrar como anotaciones.
#'   Default: NULL (usa todas las columnas excepto FeatureID)
#' @param row_annotation_palette Lista nombrada de paletas para cada anotación.
#'   Ejemplo: list(Pathway = "brewer:Set1", Function = c("red", "blue", "green"))
#' @param row_annotation_size Ancho de las barras de anotación de filas. Puede ser:
#'   - Número: interpretado como centímetros (ej: 0.3 = 0.3cm)
#'   - Objeto unit: grid::unit(0.3, "cm")
#'   - NULL: usa el valor por defecto de tidyHeatmap
#' @param row_annotation_name_size Tamaño de fuente del nombre de las anotaciones de fila (default: 8)
#' @param row_order_by Nombre de columna de anotación para ordenar filas (default: NULL)
#' @param split_rows_by Nombre de columna de anotación para separar filas en grupos (default: NULL)
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
#' @param row_title_size Tamaño de fuente del título de filas (default: 10)
#' @param column_title_size Tamaño de fuente del título de columnas (default: 10)
#' @param show_row_title Mostrar título de filas (default: TRUE)
#' @param show_column_title Mostrar título de columnas (default: TRUE)
#' @param show_annotation Mostrar anotación de Condition (default: TRUE)
#' @param split_by_condition Dividir el heatmap por condición (default: FALSE)
#' @param show_adjp_annotation Mostrar anotación de adjP en filas para mode="target" (default: TRUE)
#' @param palette_adjp Paleta para anotación de adjP. Puede ser:
#'   - NULL: usa paleta por defecto (rojo-naranja-blanco)
#'   - Vector de 3 colores: c(color_0, color_medio, color_0.05)
#'   - String "brewer:nombre": usa paleta de RColorBrewer
#'   - String "paquete::paleta": usa paleta de paletteer
#' @param row_names_size Tamaño de fuente de nombres de fila (default: 7)
#' @param column_names_size Tamaño de fuente de nombres de columna (default: 9)
#' @param column_names_rotation Rotación de nombres de columna en grados (default: 45)
#' @param column_dend_height Altura del dendrograma de columnas. Puede ser:
#'   - Número: interpretado como milímetros (ej: 30 = 30mm)
#'   - Objeto unit: grid::unit(2, "cm")
#'   - NULL: usa el valor por defecto de ComplexHeatmap
#' @param row_dend_width Ancho del dendrograma de filas. Mismo formato que column_dend_height
#' @param show_heatmap_legend Mostrar leyenda del heatmap (default: TRUE)
#' @param show_annotation_legend Mostrar leyenda de las anotaciones (default: TRUE)
#' @param border_color Color del borde de las celdas del heatmap. Puede ser:
#'   - NULL o FALSE: sin borde (default)
#'   - TRUE: borde negro
#'   - String de color: color específico (ej: "black", "grey", "#CCCCCC")
#' @param heatmap_title Título principal del heatmap (default: NULL, sin título)
#' @param heatmap_title_size Tamaño de fuente del título principal (default: 14)
#' @param heatmap_title_face Estilo de fuente del título: "plain", "bold", "italic", "bold.italic" (default: "bold")
#' @param export_path Ruta para exportar los datos del heatmap a TSV (default: NULL, no exporta).
#'   El archivo incluirá FeatureID, valores de intensidad por muestra, y metadatos (adjP si aplica).
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
                               feature_ids = NULL,
                               row_annotation = NULL,
                               row_annotation_cols = NULL,
                               row_annotation_palette = NULL,
                               row_annotation_size = NULL,
                               row_annotation_name_size = 8,
                               row_order_by = NULL,
                               split_rows_by = NULL,
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
                               row_title_size = 10,
                               column_title_size = 10,
                               show_row_title = TRUE,
                               show_column_title = TRUE,
                               show_annotation = TRUE,
                               split_by_condition = FALSE,
                               show_adjp_annotation = TRUE,
                               palette_adjp = NULL,
                               row_names_size = 7,
                               column_names_size = 9,
                               column_names_rotation = 45,
                               column_dend_height = NULL,
                               row_dend_width = NULL,
                               show_heatmap_legend = TRUE,
                               show_annotation_legend = TRUE,
                               border_color = NULL,
                               heatmap_title = NULL,
                               heatmap_title_size = 14,
                               heatmap_title_face = c("bold", "plain", "italic", "bold.italic"),
                               export_path = NULL) {

  mode <- match.arg(mode)
  scale_data <- match.arg(scale_data)
  heatmap_title_face <- match.arg(heatmap_title_face)

  # Procesar border_color
  rect_gp <- NULL
  if (!is.null(border_color) && !isFALSE(border_color)) {
    if (isTRUE(border_color)) {
      rect_gp <- grid::gpar(col = "black")
    } else {
      rect_gp <- grid::gpar(col = border_color)
    }
  }

  # Procesar títulos (ocultar si show_*_title es FALSE)
  if (!show_row_title) {
    row_title <- NULL
  }
  if (!show_column_title) {
    column_title <- NULL
  }

  # Convertir tamaños de dendrogramas a unidades grid si son numéricos
  if (!is.null(column_dend_height) && is.numeric(column_dend_height)) {
    column_dend_height <- grid::unit(column_dend_height, "mm")
  }
  if (!is.null(row_dend_width) && is.numeric(row_dend_width)) {
    row_dend_width <- grid::unit(row_dend_width, "mm")
  }

  # ---------------------------------------------------------------------------
  # 1) Preparar datos
  # ---------------------------------------------------------------------------

  hm_data <- prepare_heatmap_data(
    data = data,
    mode = mode,
    alpha = alpha,
    comparison = comparison,
    feature_ids = feature_ids,
    scale_data = scale_data,
    sample_order = sample_order,
    condition_order = condition_order
  )

  # ---------------------------------------------------------------------------
  # 1b) Exportar datos a TSV si se solicita
  # ---------------------------------------------------------------------------

  if (!is.null(export_path) && nzchar(export_path)) {
    # Crear matriz wide con FeatureID como filas y SampleID como columnas
    export_wide <- hm_data %>%
      select(FeatureID, SampleID, Intensity) %>%
      tidyr::pivot_wider(
        names_from = SampleID,
        values_from = Intensity,
        values_fn = mean
      )

    # Añadir adjP si existe
    if ("adjP" %in% names(hm_data)) {
      adjp_data <- hm_data %>%
        select(FeatureID, adjP) %>%
        distinct()
      export_wide <- export_wide %>%
        left_join(adjp_data, by = "FeatureID")
    }

    # Exportar a TSV
    readr::write_tsv(export_wide, export_path)
    message("Datos exportados a: ", export_path)
  }

  # ---------------------------------------------------------------------------
  # 1c) Procesar anotaciones de fila
  # ---------------------------------------------------------------------------

  row_annot_data <- NULL
  row_annot_colors <- list()
  row_split_vector <- NULL

  if (!is.null(row_annotation)) {
    # Cargar anotaciones si es ruta a archivo
    if (is.character(row_annotation) && length(row_annotation) == 1 && file.exists(row_annotation)) {
      row_annot_data <- readr::read_tsv(row_annotation, show_col_types = FALSE)
    } else if (is.data.frame(row_annotation)) {
      row_annot_data <- as.data.frame(row_annotation)
    } else {
      warning("row_annotation debe ser ruta a archivo TSV o data.frame")
    }

    if (!is.null(row_annot_data)) {
      # Validar columna FeatureID
      if (!("FeatureID" %in% names(row_annot_data))) {
        stop("El archivo de anotaciones debe tener una columna 'FeatureID'")
      }

      # Determinar columnas a usar
      all_annot_cols <- setdiff(names(row_annot_data), "FeatureID")
      if (is.null(row_annotation_cols)) {
        row_annotation_cols <- all_annot_cols
      } else {
        missing_cols <- setdiff(row_annotation_cols, all_annot_cols)
        if (length(missing_cols) > 0) {
          warning("Columnas de anotación no encontradas: ", paste(missing_cols, collapse = ", "))
          row_annotation_cols <- intersect(row_annotation_cols, all_annot_cols)
        }
      }

      # Convertir valores vacíos al string "NA" (no NA_character_) para asignar color blanco
      for (col in row_annotation_cols) {
        values <- row_annot_data[[col]]
        if (is.character(values)) {
          # Convertir strings vacíos, espacios en blanco, y NA reales a string "NA"
          row_annot_data[[col]] <- ifelse(
            is.na(values) | trimws(values) == "",
            "NA",
            values
          )
        }
      }

      # Filtrar solo FeatureIDs presentes en los datos
      feature_ids_in_data <- unique(hm_data$FeatureID)
      row_annot_data <- row_annot_data[row_annot_data$FeatureID %in% feature_ids_in_data, , drop = FALSE]

      # Preparar paletas de colores para cada anotación (incluyendo "NA" con blanco)
      for (col in row_annotation_cols) {
        # Obtener niveles únicos
        all_levels <- unique(row_annot_data[[col]])
        has_na <- "NA" %in% all_levels

        # Separar niveles reales de "NA" y ordenar: primero reales, luego "NA"
        real_levels <- sort(all_levels[all_levels != "NA"])
        if (has_na) {
          ordered_levels <- c(real_levels, "NA")
        } else {
          ordered_levels <- real_levels
        }

        # Convertir la columna a factor con niveles ordenados
        row_annot_data[[col]] <- factor(row_annot_data[[col]], levels = ordered_levels)

        # Obtener colores para niveles reales
        if (!is.null(row_annotation_palette) && col %in% names(row_annotation_palette)) {
          real_colors <- get_annotation_palette(real_levels, row_annotation_palette[[col]])
        } else {
          real_colors <- get_annotation_palette(real_levels, NULL)
        }

        # Construir paleta final en el EXACTO orden de los niveles del factor
        # Esto asegura que tidyHeatmap asigne los colores correctamente
        final_palette <- character(length(ordered_levels))
        names(final_palette) <- ordered_levels
        for (lvl in ordered_levels) {
          if (lvl == "NA") {
            final_palette[lvl] <- "#FFFFFF"
          } else {
            final_palette[lvl] <- real_colors[lvl]
          }
        }
        row_annot_colors[[col]] <- final_palette
      }

      # Unir anotaciones con hm_data
      cols_to_join <- c("FeatureID", row_annotation_cols)
      hm_data <- hm_data %>%
        left_join(row_annot_data[, cols_to_join, drop = FALSE], by = "FeatureID")

      # Asegurar que los factores se mantienen después del join
      for (col in row_annotation_cols) {
        if (col %in% names(hm_data) && col %in% names(row_annot_data)) {
          hm_data[[col]] <- factor(hm_data[[col]], levels = levels(row_annot_data[[col]]))
        }
      }

      # Ordenar filas por anotación si se especifica
      if (!is.null(row_order_by) && row_order_by %in% row_annotation_cols) {
        # Crear un orden explícito basado en la columna de anotación
        # Primero ordenar row_annot_data por la columna especificada
        row_annot_sorted <- row_annot_data %>%
          arrange(!!rlang::sym(row_order_by))

        # Obtener el orden de FeatureIDs
        feature_order <- unique(row_annot_sorted$FeatureID)

        # Convertir FeatureID a factor con el orden correcto
        hm_data <- hm_data %>%
          mutate(FeatureID = factor(FeatureID, levels = feature_order)) %>%
          arrange(FeatureID)  # Ordenar explícitamente los datos
      }

      # Preparar row_split si se especifica (DESPUÉS de ordenar)
      if (!is.null(split_rows_by) && split_rows_by %in% row_annotation_cols) {
        # Obtener FeatureIDs únicos en el orden actual de hm_data
        if (!is.null(row_order_by)) {
          # Si hay orden, usar los niveles del factor
          ordered_features <- levels(hm_data$FeatureID)
        } else {
          ordered_features <- unique(as.character(hm_data$FeatureID))
        }

        # Crear mapping de FeatureID a valor de split
        feature_to_split <- row_annot_data %>%
          select(FeatureID, all_of(split_rows_by)) %>%
          distinct()
        feature_to_split <- stats::setNames(
          feature_to_split[[split_rows_by]],
          feature_to_split$FeatureID
        )

        # Crear vector de split en el orden correcto
        split_values <- feature_to_split[ordered_features]

        # Determinar orden de niveles del split (según aparición en datos ordenados)
        split_levels <- unique(split_values)
        split_levels <- split_levels[!is.na(split_levels)]

        row_split_vector <- factor(split_values, levels = split_levels)
      }
    }
  }

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
  if (is.character(sample_order) && (length(sample_order) > 1 || sample_order != "clustering")) {
    cluster_cols_final <- FALSE
  }

  # Determinar si aplicar clustering a filas
  cluster_rows_final <- cluster_rows
  if (!is.null(row_order_by) || !is.null(split_rows_by)) {
    cluster_rows_final <- FALSE
  }

  # ---------------------------------------------------------------------------
  # 4) Crear heatmap con tidyHeatmap
  # ---------------------------------------------------------------------------

  # Convertir a tibble (requerido por tidyHeatmap)
  hm_data <- tibble::as_tibble(hm_data)

  # Preparar column_split si se solicita (para separación visual entre condiciones)
  col_split_vector <- NULL
  if (split_by_condition) {
    # Obtener el orden de muestras únicas tal como aparecen en los datos
    sample_info <- hm_data %>%
      select(SampleID, Condition) %>%
      distinct()

    # Mantener el orden original de SampleID en los datos
    sample_order_vec <- unique(as.character(hm_data$SampleID))
    sample_info <- sample_info[match(sample_order_vec, as.character(sample_info$SampleID)), ]
    col_split_vector <- factor(sample_info$Condition, levels = unique(sample_info$Condition))
  }

  # Preparar argumentos opcionales
  extra_args <- list()
  if (!is.null(column_dend_height)) {
    extra_args$column_dend_height <- column_dend_height
  }
  if (!is.null(row_dend_width)) {
    extra_args$row_dend_width <- row_dend_width
  }
  if (!is.null(rect_gp)) {
    extra_args$rect_gp <- rect_gp
  }
  if (!is.null(row_split_vector)) {
    extra_args$row_split <- row_split_vector
  }

  # Crear heatmap base
  hm_args <- c(
    list(
      .data = hm_data,
      .row = rlang::sym("FeatureID"),
      .column = rlang::sym("SampleID"),
      .value = rlang::sym("Intensity"),
      scale = "none",
      cluster_rows = cluster_rows_final,
      cluster_columns = cluster_cols_final,
      palette_value = palette_func,
      show_row_names = show_row_names,
      show_column_names = show_column_names,
      row_names_gp = grid::gpar(fontsize = row_names_size),
      column_names_gp = grid::gpar(fontsize = column_names_size),
      column_names_rot = column_names_rotation,
      row_title = row_title,
      column_title = column_title,
      row_title_gp = grid::gpar(fontsize = row_title_size),
      column_title_gp = grid::gpar(fontsize = column_title_size),
      column_split = col_split_vector,
      show_heatmap_legend = show_heatmap_legend
    ),
    extra_args
  )

  hm <- do.call(tidyHeatmap::heatmap, hm_args)

  # Añadir anotación de Condition como barra de color si se solicita
  if (show_annotation) {
    hm <- hm %>%
      annotation_tile(
        Condition,
        palette = annotation_colors,
        show_legend = show_annotation_legend
      )
  }

  # Añadir anotaciones de fila personalizadas (suppressWarnings para manejar NA silenciosamente)
  if (!is.null(row_annot_data) && length(row_annotation_cols) > 0) {
    # Preparar size como unit si es numérico
    annot_size <- NULL
    if (!is.null(row_annotation_size)) {
      if (is.numeric(row_annotation_size)) {
        annot_size <- grid::unit(row_annotation_size, "cm")
      } else {
        annot_size <- row_annotation_size
      }
    }

    # Preparar annotation_name_gp
    annot_name_gp <- grid::gpar(fontsize = row_annotation_name_size)

    # Identificar las columnas válidas
    valid_cols <- row_annotation_cols[row_annotation_cols %in% names(hm_data)]
    n_cols <- length(valid_cols)

    for (i in seq_along(valid_cols)) {
      col <- valid_cols[i]
      # Solo aplicar size en la última anotación (evita warning de tidyHeatmap)
      use_size <- if (i == n_cols) annot_size else NULL

      hm <- hm %>%
        annotation_tile(
          !!rlang::sym(col),
          palette = row_annot_colors[[col]],
          size = use_size,
          annotation_name_gp = annot_name_gp,
          show_legend = show_annotation_legend
        )
    }
  }

  # Añadir anotación de adjP para mode = "target" (anotación de filas)
  if (show_adjp_annotation && mode == "target" && "adjP" %in% names(hm_data)) {
    # Obtener colores para la paleta de adjP
    adjp_colors <- c("#67001F", "#F4A582", "#F7F7F7")  # Default: rojo oscuro -> naranja -> blanco

    if (!is.null(palette_adjp)) {
      if (is.character(palette_adjp) && length(palette_adjp) >= 3) {
        # Vector de colores personalizado
        adjp_colors <- palette_adjp[1:3]
      } else if (is.character(palette_adjp) && length(palette_adjp) == 1) {
        # Paleta de paletteer o RColorBrewer
        if (grepl("::", palette_adjp)) {
          # Paletteer
          if (requireNamespace("paletteer", quietly = TRUE)) {
            raw_pal <- tryCatch({
              as.character(paletteer::paletteer_c(palette_adjp, n = 3))
            }, error = function(e) {
              tryCatch({
                as.character(paletteer::paletteer_d(palette_adjp, n = 3))
              }, error = function(e2) NULL)
            })
            if (!is.null(raw_pal) && length(raw_pal) >= 3) {
              adjp_colors <- sapply(raw_pal[1:3], normalize_hex, USE.NAMES = FALSE)
            }
          }
        } else if (startsWith(palette_adjp, "brewer:")) {
          # RColorBrewer
          if (requireNamespace("RColorBrewer", quietly = TRUE)) {
            nm <- sub("^brewer:", "", palette_adjp)
            raw_pal <- tryCatch({
              RColorBrewer::brewer.pal(3, nm)
            }, error = function(e) NULL)
            if (!is.null(raw_pal)) {
              adjp_colors <- raw_pal
            }
          }
        }
      }
    }

    # Crear paleta para p-valores (valores bajos = más significativos)
    adjp_palette <- circlize::colorRamp2(
      c(0, 0.01, 0.05),
      adjp_colors
    )

    hm <- hm %>%
      annotation_tile(
        adjP,
        palette = adjp_palette,
        show_legend = show_annotation_legend
      )
  }

  # ---------------------------------------------------------------------------
  # 5) Añadir título principal si se especifica
  # ---------------------------------------------------------------------------

  if (!is.null(heatmap_title) && nzchar(heatmap_title)) {
    hm <- wrap_heatmap_with_title(
      hm = hm,
      title = heatmap_title,
      title_size = heatmap_title_size,
      title_face = heatmap_title_face
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
#' @param feature_ids Vector de FeatureIDs específicos a mostrar (default: NULL)
#' @param row_annotation Anotaciones para filas (ruta TSV o data.frame)
#' @param row_annotation_cols Columnas a usar como anotaciones de fila
#' @param row_annotation_palette Lista de paletas para anotaciones de fila
#' @param row_annotation_size Ancho de las barras de anotación de filas (número en cm o unit)
#' @param row_annotation_name_size Tamaño de fuente del nombre de anotaciones de fila (default: 8)
#' @param row_order_by Columna de anotación para ordenar filas
#' @param split_rows_by Columna de anotación para separar filas en grupos
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
#' @param row_title_size Tamaño de fuente del título de filas (default: 10)
#' @param column_title_size Tamaño de fuente del título de columnas (default: 10)
#' @param show_row_title Mostrar título de filas (default: TRUE)
#' @param show_column_title Mostrar título de columnas (default: TRUE)
#' @param show_annotation Mostrar anotación de Condition (default: TRUE)
#' @param split_by_condition Dividir el heatmap por condición (default: FALSE)
#' @param show_adjp_annotation Mostrar anotación de adjP en filas para comparaciones (default: TRUE)
#' @param palette_adjp Paleta para anotación de adjP (ver proteomics_heatmap)
#' @param row_names_size Tamaño de fuente de nombres de fila (default: 7)
#' @param column_names_size Tamaño de fuente de nombres de columna (default: 9)
#' @param column_names_rotation Rotación de nombres de columna (default: 45)
#' @param column_dend_height Altura del dendrograma de columnas (número en mm o unit)
#' @param row_dend_width Ancho del dendrograma de filas (número en mm o unit)
#' @param show_heatmap_legend Mostrar leyenda del heatmap (default: TRUE)
#' @param show_annotation_legend Mostrar leyenda de las anotaciones (default: TRUE)
#' @param border_color Color del borde de las celdas (NULL, TRUE, o color)
#' @param heatmap_title Título principal del heatmap (default: NULL, sin título).
#'   Se puede usar "\{mode\}" como placeholder que será reemplazado por el nombre del modo
#' @param heatmap_title_size Tamaño de fuente del título principal (default: 14)
#' @param heatmap_title_face Estilo de fuente del título (default: "bold")
#' @param export_path Ruta base para exportar datos a TSV (default: NULL).
#'   Se añadirá el nombre del modo al archivo (ej: "export_all.tsv", "export_B-A.tsv")
#' @param export_modes Vector de modos a exportar (default: NULL, exporta todos).
#'   Solo aplica si export_path está definido. Ejemplo: c("all", "B-A")
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
                                    feature_ids = NULL,
                                    row_annotation = NULL,
                                    row_annotation_cols = NULL,
                                    row_annotation_palette = NULL,
                                    row_annotation_size = NULL,
                                    row_annotation_name_size = 8,
                                    row_order_by = NULL,
                                    split_rows_by = NULL,
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
                                    row_title_size = 10,
                                    column_title_size = 10,
                                    show_row_title = TRUE,
                                    show_column_title = TRUE,
                                    show_annotation = TRUE,
                                    split_by_condition = FALSE,
                                    show_adjp_annotation = TRUE,
                                    palette_adjp = NULL,
                                    row_names_size = 7,
                                    column_names_size = 9,
                                    column_names_rotation = 45,
                                    column_dend_height = NULL,
                                    row_dend_width = NULL,
                                    show_heatmap_legend = TRUE,
                                    show_annotation_legend = TRUE,
                                    border_color = NULL,
                                    heatmap_title = NULL,
                                    heatmap_title_size = 14,
                                    heatmap_title_face = "bold",
                                    export_path = NULL,
                                    export_modes = NULL) {

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

    # Procesar título (reemplazar {mode} si existe)
    current_title <- NULL
    if (!is.null(heatmap_title)) {
      current_title <- gsub("\\{mode\\}", m, heatmap_title)
    }

    # Procesar export_path (añadir modo al nombre del archivo)
    current_export_path <- NULL
    if (!is.null(export_path) && nzchar(export_path)) {
      # Verificar si este modo debe exportarse
      should_export <- is.null(export_modes) || m %in% export_modes
      if (should_export) {
        # Separar directorio, nombre y extensión
        dir_path <- dirname(export_path)
        base_name <- tools::file_path_sans_ext(basename(export_path))
        ext <- tools::file_ext(export_path)
        if (nzchar(ext)) ext <- paste0(".", ext) else ext <- ".tsv"
        current_export_path <- file.path(dir_path, paste0(base_name, "_", m, ext))
      }
    }

    # Intentar generar heatmap (puede fallar si no hay suficientes proteínas)
    hm <- tryCatch({
      proteomics_heatmap(
        data = data,
        mode = internal_mode,
        alpha = alpha,
        comparison = comparison,
        feature_ids = feature_ids,
        row_annotation = row_annotation,
        row_annotation_cols = row_annotation_cols,
        row_annotation_palette = row_annotation_palette,
        row_annotation_size = row_annotation_size,
        row_annotation_name_size = row_annotation_name_size,
        row_order_by = row_order_by,
        split_rows_by = split_rows_by,
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
        row_title_size = row_title_size,
        column_title_size = column_title_size,
        show_row_title = show_row_title,
        show_column_title = show_column_title,
        show_annotation = show_annotation,
        split_by_condition = split_by_condition,
        show_adjp_annotation = show_adjp_annotation,
        palette_adjp = palette_adjp,
        row_names_size = row_names_size,
        column_names_size = column_names_size,
        column_names_rotation = column_names_rotation,
        column_dend_height = column_dend_height,
        row_dend_width = row_dend_width,
        show_heatmap_legend = show_heatmap_legend,
        show_annotation_legend = show_annotation_legend,
        border_color = border_color,
        heatmap_title = current_title,
        heatmap_title_size = heatmap_title_size,
        heatmap_title_face = heatmap_title_face,
        export_path = current_export_path
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

# --- Heatmap con título personalizado ---
# hm_title <- proteomics_heatmap(
#   data = hm_input,
#   mode = "any",
#   scale_data = "row",
#   heatmap_title = "Differential Expression Analysis",
#   heatmap_title_size = 16,
#   heatmap_title_face = "bold"
# )
# hm_title

# --- Lista de heatmaps con títulos dinámicos ---
# hm_list <- proteomics_heatmap_list(
#   data = hm_input,
#   modes = c("all", "any", "B-A"),
#   heatmap_title = "Proteomics Heatmap: {mode}",
#   heatmap_title_size = 14
# )
# # Los títulos serán: "Proteomics Heatmap: all", "Proteomics Heatmap: any", "Proteomics Heatmap: B-A"

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
