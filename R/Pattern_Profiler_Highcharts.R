# =============================================================================
# Pattern Profiler (Clustering) con Highcharts para Datos de Proteómica
# =============================================================================
#
# Módulo para identificar patrones de expresión (clusters) en proteínas
# diferencialmente expresadas usando soft clustering (Mfuzz).
#
# Genera visualizaciones interactivas con Highcharts mostrando los perfiles
# de expresión por cluster a través de las condiciones experimentales.
#
# =============================================================================

library(highcharter)
library(dplyr)
library(tidyr)

# Dependencias de Bioconductor (requeridas para clustering)
if (!requireNamespace("Biobase", quietly = TRUE)) {
  stop("Instala 'Biobase' desde Bioconductor: BiocManager::install('Biobase')")
}
if (!requireNamespace("Mfuzz", quietly = TRUE)) {
  stop("Instala 'Mfuzz' desde Bioconductor: BiocManager::install('Mfuzz')")
}
library(Biobase)
library(Mfuzz)


# -----------------------------------------------------------------------------
# Operador null-coalesce
# -----------------------------------------------------------------------------

`%||%` <- function(a, b) if (!is.null(a) && length(a) && !is.na(a[1])) a else b


# -----------------------------------------------------------------------------
# Funciones auxiliares para colores
# -----------------------------------------------------------------------------

#' Normalizar color hex (eliminar canal alpha si existe)
#'
#' @param hex Color en formato hexadecimal
#' @return Color hex normalizado (#RRGGBB)
normalize_hex <- function(hex) {

hex <- gsub("^#", "", hex)
  if (nchar(hex) == 8) {
    hex <- substr(hex, 1, 6)
  }
  paste0("#", hex)
}


#' Convertir color hex a rgba
#'
#' @param hex Color en formato hexadecimal
#' @param alpha Opacidad (0-1)
#' @return String rgba
hex_to_rgba <- function(hex, alpha = 0.7) {
  hex <- normalize_hex(hex)
  hex <- gsub("^#", "", hex)
  r <- strtoi(substr(hex, 1, 2), base = 16)
  g <- strtoi(substr(hex, 3, 4), base = 16)
  b <- strtoi(substr(hex, 5, 6), base = 16)
  sprintf("rgba(%d, %d, %d, %.2f)", r, g, b, alpha)
}


#' Oscurecer un color hex
#'
#' @param hex Color en formato hexadecimal
#' @param factor Factor de oscurecimiento (0-1)
#' @return Color hex oscurecido
darken_hex <- function(hex, factor = 0.3) {
  hex <- normalize_hex(hex)
  hex <- gsub("^#", "", hex)
  r <- strtoi(substr(hex, 1, 2), base = 16)
  g <- strtoi(substr(hex, 3, 4), base = 16)
  b <- strtoi(substr(hex, 5, 6), base = 16)
  r <- max(0, round(r * (1 - factor)))
  g <- max(0, round(g * (1 - factor)))
  b <- max(0, round(b * (1 - factor)))
  sprintf("#%02X%02X%02X", r, g, b)
}


#' Interpolar color en gradiente basado en valor
#'
#' @param value Valor a mapear (típicamente z-score o membership)
#' @param low_color Color para valores bajos (hex o nombre)
#' @param mid_color Color para valores medios (hex o nombre)
#' @param high_color Color para valores altos (hex o nombre)
#' @param midpoint Punto medio del gradiente
#' @param limits Vector c(min, max) para los límites
#' @return Color hex interpolado
interpolate_color <- function(value, low_color = "#2166AC", mid_color = "#CCCCCC",
                               high_color = "#B2182B", midpoint = 0,
                               limits = c(-3, 3)) {

  # Convertir nombres de colores a hex si es necesario
  to_hex <- function(col) {
    if (grepl("^#", col)) return(col)
    rgb_vals <- col2rgb(col)
    sprintf("#%02X%02X%02X", rgb_vals[1], rgb_vals[2], rgb_vals[3])
  }

  low_color <- to_hex(low_color)
  mid_color <- to_hex(mid_color)
  high_color <- to_hex(high_color)

  # Función para parsear hex a RGB
  hex_to_rgb <- function(hex) {
    hex <- gsub("^#", "", hex)
    c(
      strtoi(substr(hex, 1, 2), base = 16),
      strtoi(substr(hex, 3, 4), base = 16),
      strtoi(substr(hex, 5, 6), base = 16)
    )
  }

  # Normalizar valor a los límites
  value <- max(limits[1], min(limits[2], value))

  # Evitar división por cero
  lower_range <- midpoint - limits[1]
  upper_range <- limits[2] - midpoint

  if (value <= midpoint) {
    t <- if (lower_range > 0) (value - limits[1]) / lower_range else 0.5
    col1 <- hex_to_rgb(low_color)
    col2 <- hex_to_rgb(mid_color)
  } else {
    t <- if (upper_range > 0) (value - midpoint) / upper_range else 0.5
    col1 <- hex_to_rgb(mid_color)
    col2 <- hex_to_rgb(high_color)
  }

  r <- round(col1[1] + t * (col2[1] - col1[1]))
  g <- round(col1[2] + t * (col2[2] - col1[2]))
  b <- round(col1[3] + t * (col2[3] - col1[3]))

  # Asegurar valores en rango válido
  r <- max(0, min(255, r))
  g <- max(0, min(255, g))
  b <- max(0, min(255, b))

  sprintf("#%02X%02X%02X", r, g, b)
}


# -----------------------------------------------------------------------------
# Funciones de configuración de paleta
# -----------------------------------------------------------------------------

#' Configurar paleta de colores para clusters
#'
#' @param n_clusters Número de clusters
#' @param palette Paleta de colores (NULL, "ggsci::palette", "brewer:Name", o vector)
#' @return Vector de colores
configure_cluster_palette <- function(n_clusters, palette = NULL) {

  default_palette <- c(
    "#457B9D", "#E63946", "#2A9D8F", "#E9C46A",
    "#9B5DE5", "#F4A261", "#264653", "#00BBF9",
    "#FF006E", "#8338EC", "#3A86FF", "#FB5607"
  )

  if (is.null(palette)) {
    pal <- default_palette
  } else if (is.character(palette) && length(palette) == 1 && grepl("::", palette)) {
    if (!requireNamespace("paletteer", quietly = TRUE)) {
      stop("Para usar paletteer, instala con: install.packages('paletteer')")
    }
    pal <- tryCatch({
      raw_pal <- as.character(paletteer::paletteer_d(palette))
      sapply(raw_pal, normalize_hex, USE.NAMES = FALSE)
    }, error = function(e) {
      warning("Error cargando paleta '", palette, "'. Usando default.")
      default_palette
    })
  } else if (is.character(palette) && length(palette) == 1 && startsWith(palette, "brewer:")) {
    if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
      stop("Para usar RColorBrewer, instala con: install.packages('RColorBrewer')")
    }
    nm <- sub("^brewer:", "", palette)
    pal <- tryCatch({
      maxc <- RColorBrewer::brewer.pal.info[nm, "maxcolors"]
      RColorBrewer::brewer.pal(maxc, nm)
    }, error = function(e) {
      warning("Error cargando paleta brewer '", nm, "'. Usando default.")
      default_palette
    })
  } else if (is.character(palette)) {
    pal <- palette
  } else {
    pal <- default_palette
  }

  if (length(pal) < n_clusters) {
    pal <- rep(pal, length.out = n_clusters)
  }

  pal[seq_len(n_clusters)]
}


# =============================================================================
# FUNCIONES DE PREPARACIÓN DE DATOS
# =============================================================================

#' Preparar datos para clustering desde archivo parquet/TSV
#'
#' Filtra proteínas significativas y agrega intensidades por condición.
#'
#' @param data Data frame en formato long con columnas:
#'   - SampleID: Identificador de muestra
#'   - FeatureID: Identificador de proteína
#'   - Intensity: Valor de intensidad (log2)
#'   - Condition: Condición experimental
#'   - sig_any: Lógico indicando significancia (opcional)
#'   - adjP_*: Columnas de p-valores ajustados (opcional)
#' @param filter_mode Modo de filtrado: "all", "any" (usa sig_any), o "specific"
#' @param alpha Umbral de significancia para modo "specific" (default: 0.05)
#' @param comparison Comparación específica para filtrar (ej: "B-A")
#' @param condition_order Orden de las condiciones (opcional)
#' @param aggregate Método de agregación por condición: "mean" o "median"
#'
#' @return Lista con:
#'   - matrix: Matriz de intensidades agregadas (proteínas x condiciones)
#'   - feature_info: Data frame con información de features
#'   - conditions: Vector de condiciones ordenadas
prepare_clustering_data <- function(data,
                                     filter_mode = c("any", "all", "specific"),
                                     alpha = 0.05,
                                     comparison = NULL,
                                     condition_order = NULL,
                                     aggregate = c("median", "mean")) {

  filter_mode <- match.arg(filter_mode)
  aggregate <- match.arg(aggregate)

  # ---------------------------------------------------------------------------
  # 1) Validación de columnas requeridas
  # ---------------------------------------------------------------------------
  required_cols <- c("SampleID", "FeatureID", "Intensity", "Condition")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Columnas requeridas faltantes: ", paste(missing_cols, collapse = ", "))
  }

  # Convertir a data.frame
  data <- as.data.frame(data)

  # ---------------------------------------------------------------------------
  # 2) Filtrar proteínas según el modo
  # ---------------------------------------------------------------------------
  # Obtener features únicos con su metadata
  feat_cols <- c("FeatureID", grep("^adjP_|^sig_", names(data), value = TRUE))
  feat_cols <- intersect(feat_cols, names(data))
  feat_info <- data[!duplicated(data$FeatureID), feat_cols, drop = FALSE]
  rownames(feat_info) <- feat_info$FeatureID

  if (filter_mode == "all") {
    selected_ids <- feat_info$FeatureID

  } else if (filter_mode == "any") {
    if (!("sig_any" %in% names(feat_info))) {
      stop("La columna 'sig_any' es requerida para filter_mode = 'any'")
    }
    selected_ids <- feat_info$FeatureID[feat_info$sig_any == TRUE]

  } else {
    # filter_mode == "specific"
    if (is.null(comparison)) {
      stop("El argumento 'comparison' es requerido para filter_mode = 'specific'")
    }
    col_name <- paste0("adjP_", comparison)
    if (!(col_name %in% names(feat_info))) {
      stop("No existe la columna: ", col_name)
    }
    selected_ids <- feat_info$FeatureID[
      !is.na(feat_info[[col_name]]) & feat_info[[col_name]] <= alpha
    ]
  }

  if (length(selected_ids) < 3) {
    stop("Se requieren al menos 3 proteínas para clustering. Encontradas: ",
         length(selected_ids))
  }

  # Filtrar datos
  data_filtered <- data[data$FeatureID %in% selected_ids, , drop = FALSE]
  feat_info <- feat_info[selected_ids, , drop = FALSE]

  # ---------------------------------------------------------------------------
  # 3) Determinar orden de condiciones
  # ---------------------------------------------------------------------------
  if (!is.null(condition_order)) {
    conditions <- condition_order[condition_order %in% unique(data_filtered$Condition)]
  } else {
    conditions <- sort(unique(data_filtered$Condition))
  }

  if (length(conditions) < 2) {
    stop("Se requieren al menos 2 condiciones para clustering.")
  }

  # ---------------------------------------------------------------------------
  # 4) Agregar intensidades por condición
  # ---------------------------------------------------------------------------
  agg_fun <- switch(aggregate,
    mean = function(x) mean(x, na.rm = TRUE),
    median = function(x) median(x, na.rm = TRUE)
  )

  # Pivotar a matriz wide y agregar
  data_wide <- data_filtered %>%
    filter(Condition %in% conditions) %>%
    group_by(FeatureID, Condition) %>%
    summarise(Intensity = agg_fun(Intensity), .groups = "drop") %>%
    pivot_wider(names_from = Condition, values_from = Intensity) %>%
    as.data.frame()

  rownames(data_wide) <- data_wide$FeatureID
  mat <- as.matrix(data_wide[, conditions, drop = FALSE])

  # Eliminar filas con NA
  complete_rows <- complete.cases(mat)
  if (sum(complete_rows) < 3) {
    stop("Menos de 3 proteínas con datos completos tras filtrado.")
  }
  mat <- mat[complete_rows, , drop = FALSE]
  feat_info <- feat_info[rownames(mat), , drop = FALSE]

  message(sprintf("Datos preparados: %d proteínas x %d condiciones",
                  nrow(mat), ncol(mat)))

  list(
    matrix = mat,
    feature_info = feat_info,
    conditions = conditions
  )
}


# =============================================================================
# FUNCIONES DE CLUSTERING (MFUZZ)
# =============================================================================

#' Crear ExpressionSet desde matriz de datos
#'
#' @param mat Matriz de intensidades (proteínas x condiciones)
#' @param feature_info Data frame con información de features (opcional)
#' @return Objeto ExpressionSet
create_expression_set <- function(mat, feature_info = NULL) {

  if (!requireNamespace("Biobase", quietly = TRUE)) {
    stop("Instala 'Biobase' desde Bioconductor: BiocManager::install('Biobase')")
  }

  # Asegurar que es matriz numérica
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  storage.mode(mat) <- "double"

  # Crear ExpressionSet
  eset <- Biobase::ExpressionSet(assayData = mat)

  # Añadir feature data si existe
 if (!is.null(feature_info) && nrow(feature_info) > 0) {
    fd <- feature_info[rownames(mat), , drop = FALSE]
    Biobase::fData(eset) <- as.data.frame(fd)
  }

  eset
}


#' Estandarizar ExpressionSet (z-score por fila)
#'
#' @param eset Objeto ExpressionSet
#' @return ExpressionSet estandarizado
standardize_eset <- function(eset) {

  if (!requireNamespace("Mfuzz", quietly = TRUE)) {
    stop("Instala 'Mfuzz' desde Bioconductor: BiocManager::install('Mfuzz')")
  }

  # Manejar NAs si existen
  if (any(is.na(Biobase::exprs(eset)))) {
    eset <- Mfuzz::filter.NA(eset, thres = 0.25)
    eset <- Mfuzz::fill.NA(eset, mode = "mean")
  }

  Mfuzz::standardise(eset)
}


# =============================================================================
# EVALUACIÓN AUTOMÁTICA DEL NÚMERO DE CLUSTERS
# =============================================================================

#' Calcular índice Xie-Beni para clustering fuzzy
#'
#' @param X Matriz de datos (n x p)
#' @param U Matriz de membership (n x c)
#' @param m Parámetro de fuzzificación
#' @return Valor del índice Xie-Beni (menor es mejor)
.xie_beni_index <- function(X, U, m) {
  n <- nrow(X)
  c_num <- ncol(U)

  Um <- U^m

  # Calcular centros de clusters
  V <- t(sapply(seq_len(c_num), function(k) {
    colSums(Um[, k] * X) / sum(Um[, k])
  }))

  # Distancia cuadrada de cada punto a cada centro
  dist2 <- sapply(seq_len(c_num), function(k) {
    rowSums((X - matrix(V[k, ], nrow = n, ncol = ncol(X), byrow = TRUE))^2)
  })

  # Numerador: suma ponderada de distancias
  num <- sum(Um * dist2)

  # Denominador: n * distancia mínima entre centros
  if (c_num > 1) {
    dcent <- as.matrix(dist(V))^2
    diag(dcent) <- Inf
    dmin <- min(dcent)
    if (!is.finite(dmin) || dmin <= 0) return(Inf)
  } else {
    return(Inf)
  }

  num / (n * dmin)
}


#' Calcular Fuzzy Partition Coefficient (FPC)
#'
#' @param U Matriz de membership (n x c)
#' @return Valor FPC (mayor es mejor, máximo = 1)
.fpc_index <- function(U) {
  sum(U^2) / nrow(U)
}


#' Calcular Average Maximum Membership (AMM)
#'
#' @param U Matriz de membership (n x c)
#' @return Valor AMM (mayor es mejor)
.amm_index <- function(U) {
  mean(apply(U, 1, max))
}


#' Evaluar rango de valores de c para Mfuzz
#'
#' Calcula múltiples métricas de calidad para diferentes números de clusters.
#'
#' @param eset_std ExpressionSet estandarizado
#' @param c_range Rango de valores de c a evaluar (default: 2:10)
#' @param m Parámetro de fuzzificación (NULL para estimar automáticamente)
#' @param seeds Número de semillas para estabilidad (default: 5)
#' @param verbose Mostrar progreso (default: TRUE)
#'
#' @return Data frame con métricas por valor de c
evaluate_cluster_range <- function(eset_std,
                                    c_range = 2:10,
                                    m = NULL,
                                    seeds = 5L,
                                    verbose = TRUE) {

  if (!requireNamespace("Mfuzz", quietly = TRUE)) {
    stop("Instala 'Mfuzz' desde Bioconductor")
  }

  X <- Biobase::exprs(eset_std)
  if (is.null(m)) m <- Mfuzz::mestimate(eset_std)

  # Calcular Dmin (distancia mínima entre centros - guía visual)
  dmin_values <- tryCatch({
    Mfuzz::Dmin(eset_std, m = m, crange = c_range)
  }, error = function(e) rep(NA_real_, length(c_range)))

  results <- lapply(seq_along(c_range), function(i) {
    c_val <- c_range[i]
    if (verbose) message("Evaluando c = ", c_val, "...")

    # Ejecutar clustering con múltiples semillas
    metrics <- lapply(seq_len(seeds), function(s) {
      set.seed(1000 * c_val + s)
      cl <- Mfuzz::mfuzz(eset_std, c = c_val, m = m)
      U <- cl$membership[match(rownames(X), rownames(cl$membership)), , drop = FALSE]

      list(
        XB = .xie_beni_index(X, U, m),
        FPC = .fpc_index(U),
        AMM = .amm_index(U)
      )
    })

    # Promediar métricas
    XB_vals <- sapply(metrics, `[[`, "XB")
    FPC_vals <- sapply(metrics, `[[`, "FPC")
    AMM_vals <- sapply(metrics, `[[`, "AMM")

    data.frame(
      c = c_val,
      Dmin = dmin_values[i],
      XB_mean = mean(XB_vals[is.finite(XB_vals)], na.rm = TRUE),
      XB_sd = sd(XB_vals[is.finite(XB_vals)], na.rm = TRUE),
      FPC_mean = mean(FPC_vals, na.rm = TRUE),
      AMM_mean = mean(AMM_vals, na.rm = TRUE),
      m = m,
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, results)
}


#' Seleccionar número óptimo de clusters automáticamente
#'
#' Usa una combinación de métricas (Xie-Beni, FPC, AMM) para determinar
#' el número óptimo de clusters.
#'
#' @param eset_std ExpressionSet estandarizado
#' @param c_range Rango de valores de c a evaluar (default: 2:10)
#' @param m Parámetro de fuzzificación (NULL para estimar)
#' @param method Método de selección: "xb" (Xie-Beni mínimo),
#'   "consensus" (combinación de métricas), "elbow" (codo de Dmin)
#' @param verbose Mostrar información (default: TRUE)
#'
#' @return Lista con:
#'   - optimal_c: Número óptimo de clusters
#'   - metrics: Data frame con todas las métricas
#'   - m: Parámetro m usado
select_optimal_clusters <- function(eset_std,
                                     c_range = 2:10,
                                     m = NULL,
                                     method = c("xb", "consensus", "elbow"),
                                     verbose = TRUE) {

  method <- match.arg(method)

  if (!requireNamespace("Mfuzz", quietly = TRUE)) {
    stop("Instala 'Mfuzz' desde Bioconductor")
  }

  if (is.null(m)) m <- Mfuzz::mestimate(eset_std)

  # Evaluar rango de clusters
  metrics <- evaluate_cluster_range(eset_std, c_range = c_range, m = m,
                                     seeds = 5, verbose = verbose)

  # Seleccionar c óptimo según método
  if (method == "xb") {
    # Mínimo Xie-Beni
    optimal_c <- metrics$c[which.min(metrics$XB_mean)]

  } else if (method == "elbow") {
    # Método del codo usando Dmin
    dmin <- metrics$Dmin
    if (all(is.na(dmin))) {
      warning("Dmin no disponible. Usando método XB.")
      optimal_c <- metrics$c[which.min(metrics$XB_mean)]
    } else {
      # Calcular segunda derivada para detectar codo
      d1 <- diff(dmin)
      d2 <- diff(d1)
      elbow_idx <- which.max(d2) + 1
      optimal_c <- metrics$c[min(elbow_idx, length(metrics$c))]
    }

  } else {
    # method == "consensus"
    # Combinar rankings de múltiples métricas
    rank_xb <- rank(metrics$XB_mean)  # Menor es mejor
    rank_fpc <- rank(-metrics$FPC_mean)  # Mayor es mejor
    rank_amm <- rank(-metrics$AMM_mean)  # Mayor es mejor

    # Promedio de rankings (menor promedio = mejor)
    avg_rank <- (rank_xb + rank_fpc + rank_amm) / 3
    optimal_c <- metrics$c[which.min(avg_rank)]
  }

  if (verbose) {
    message(sprintf("\nNúmero óptimo de clusters (método '%s'): %d", method, optimal_c))
    message(sprintf("Parámetro m estimado: %.2f", m))
  }

  list(
    optimal_c = optimal_c,
    metrics = metrics,
    m = m
  )
}


# =============================================================================
# EJECUTAR CLUSTERING
# =============================================================================

#' Ejecutar clustering Mfuzz
#'
#' @param eset_std ExpressionSet estandarizado
#' @param c Número de clusters
#' @param m Parámetro de fuzzificación (NULL para estimar)
#' @param min_membership Pertenencia mínima para asignar a cluster (default: 0.5)
#'
#' @return Lista con:
#'   - eset_std: ExpressionSet estandarizado
#'   - cl: Objeto de clustering Mfuzz
#'   - m: Parámetro m usado
#'   - membership_table: Tabla con asignaciones y memberships
run_mfuzz_clustering <- function(eset_std, c, m = NULL, min_membership = 0.5) {

  if (!requireNamespace("Mfuzz", quietly = TRUE)) {
    stop("Instala 'Mfuzz' desde Bioconductor")
  }

  if (is.null(m)) m <- Mfuzz::mestimate(eset_std)

  set.seed(123)
  cl <- Mfuzz::mfuzz(eset_std, c = c, m = m)

  # Construir tabla de memberships
  mem <- cl$membership
  ids <- rownames(mem)

  membership_table <- data.frame(
    FeatureID = ids,
    Cluster = as.integer(cl$cluster),
    MaxMembership = apply(mem, 1, max),
    stringsAsFactors = FALSE
  )

  # Añadir membership individual por cluster
  for (k in seq_len(c)) {
    membership_table[[paste0("Membership_C", k)]] <- mem[, k]
  }

  # Añadir z-scores por condición
  zscores <- as.data.frame(Biobase::exprs(eset_std), check.names = FALSE)
  zscores$FeatureID <- rownames(zscores)
  membership_table <- merge(membership_table, zscores, by = "FeatureID", sort = FALSE)

  # Marcar los que pasan el filtro de membership
  membership_table$PassFilter <- membership_table$MaxMembership >= min_membership

  message(sprintf("Clustering completado: %d clusters, %d proteínas",
                  c, nrow(membership_table)))
  message(sprintf("Proteínas con membership >= %.2f: %d (%.1f%%)",
                  min_membership,
                  sum(membership_table$PassFilter),
                  100 * mean(membership_table$PassFilter)))

  list(
    eset_std = eset_std,
    cl = cl,
    m = m,
    c = c,
    min_membership = min_membership,
    membership_table = membership_table
  )
}


# =============================================================================
# FUNCIÓN WRAPPER PRINCIPAL
# =============================================================================

#' Pattern Profiler: Clustering de proteínas significativas
#'
#' Función principal que ejecuta todo el pipeline de clustering:
#' preparación de datos, selección automática de clusters, y clustering.
#'
#' @param data Data frame en formato long (parquet/TSV)
#' @param filter_mode Modo de filtrado: "any", "all", o "specific"
#' @param alpha Umbral de significancia (default: 0.05)
#' @param comparison Comparación específica para filter_mode = "specific"
#' @param condition_order Orden de condiciones (opcional)
#' @param aggregate Método de agregación: "median" o "mean"
#' @param c_range Rango de clusters a evaluar (default: 2:10)
#' @param auto_select_c Seleccionar c automáticamente (default: TRUE)
#' @param c Número de clusters fijo (si auto_select_c = FALSE)
#' @param selection_method Método de selección: "xb", "consensus", "elbow"
#' @param min_membership Pertenencia mínima (default: 0.5)
#' @param verbose Mostrar información de progreso (default: TRUE)
#'
#' @return Lista con todos los resultados del clustering
#'
#' @examples
#' \dontrun{
#' # Cargar datos
#' data <- arrow::read_parquet("data/PCA_Input.parquet")
#'
#' # Ejecutar clustering con selección automática
#' result <- pattern_profiler(
#'   data = data,
#'   filter_mode = "any",
#'   condition_order = c("A", "B", "C", "D"),
#'   auto_select_c = TRUE
#' )
#'
#' # Ver número óptimo de clusters
#' result$optimal_c
#'
#' # Ver tabla de resultados
#' head(result$membership_table)
#' }
pattern_profiler <- function(data,
                              filter_mode = c("any", "all", "specific"),
                              alpha = 0.05,
                              comparison = NULL,
                              condition_order = NULL,
                              aggregate = c("median", "mean"),
                              c_range = 2:10,
                              auto_select_c = TRUE,
                              c = NULL,
                              selection_method = c("xb", "consensus", "elbow"),
                              min_membership = 0.5,
                              verbose = TRUE) {

  filter_mode <- match.arg(filter_mode)
  aggregate <- match.arg(aggregate)
  selection_method <- match.arg(selection_method)

  # -------------------------------------------------------------------------
  # 1) Preparar datos
  # -------------------------------------------------------------------------
  if (verbose) message("=== Preparando datos para clustering ===")

  prep <- prepare_clustering_data(
    data = data,
    filter_mode = filter_mode,
    alpha = alpha,
    comparison = comparison,
    condition_order = condition_order,
    aggregate = aggregate
  )

  # -------------------------------------------------------------------------
  # 2) Crear y estandarizar ExpressionSet
  # -------------------------------------------------------------------------
  if (verbose) message("\n=== Estandarizando datos (z-score) ===")

  eset <- create_expression_set(prep$matrix, prep$feature_info)
  eset_std <- standardize_eset(eset)

  # -------------------------------------------------------------------------
  # 3) Seleccionar número de clusters
  # -------------------------------------------------------------------------
  if (auto_select_c) {
    if (verbose) message("\n=== Evaluando número óptimo de clusters ===")

    selection <- select_optimal_clusters(
      eset_std = eset_std,
      c_range = c_range,
      method = selection_method,
      verbose = verbose
    )
    optimal_c <- selection$optimal_c
    m <- selection$m
    metrics <- selection$metrics
  } else {
    if (is.null(c)) {
      stop("Si auto_select_c = FALSE, debes especificar 'c'")
    }
    optimal_c <- c
    m <- Mfuzz::mestimate(eset_std)
    metrics <- NULL
  }

  # -------------------------------------------------------------------------
  # 4) Ejecutar clustering
  # -------------------------------------------------------------------------
  if (verbose) message("\n=== Ejecutando clustering Mfuzz ===")

  clustering <- run_mfuzz_clustering(
    eset_std = eset_std,
    c = optimal_c,
    m = m,
    min_membership = min_membership
  )

  # -------------------------------------------------------------------------
  # 5) Preparar resultado final
  # -------------------------------------------------------------------------
  result <- list(
    # Datos preparados
    conditions = prep$conditions,
    n_features = nrow(prep$matrix),

    # Parámetros de clustering
    optimal_c = optimal_c,
    m = m,
    min_membership = min_membership,
    selection_metrics = metrics,

    # Resultados de clustering
    eset_std = clustering$eset_std,
    cl = clustering$cl,
    membership_table = clustering$membership_table,

    # Metadata
    filter_mode = filter_mode,
    aggregate = aggregate
  )

  class(result) <- c("pattern_profiler_result", class(result))

  if (verbose) message("\n=== Clustering completado ===")

  result
}


# =============================================================================
# VISUALIZACIÓN CON HIGHCHARTS
# =============================================================================

#' Preparar datos de perfiles para un cluster específico
#'
#' @param result Resultado de pattern_profiler()
#' @param cluster Número de cluster
#' @param min_membership Filtro de membership mínima (default: usa el de result)
#' @return Lista con datos preparados para visualización
prepare_cluster_profile_data <- function(result, cluster, min_membership = NULL) {

  if (is.null(min_membership)) {
    min_membership <- result$min_membership
  }

  # Filtrar proteínas del cluster con membership suficiente
  mt <- result$membership_table
  conditions <- result$conditions

  cluster_data <- mt[mt$Cluster == cluster & mt$MaxMembership >= min_membership, , drop = FALSE]

  if (nrow(cluster_data) == 0) {
    return(NULL)
  }

  # Extraer z-scores (columnas de condiciones)
  zscore_cols <- intersect(conditions, names(cluster_data))
  zscores <- as.matrix(cluster_data[, zscore_cols, drop = FALSE])
  rownames(zscores) <- cluster_data$FeatureID

  # Calcular centroide (media o mediana)
  centroid <- colMeans(zscores, na.rm = TRUE)

  list(
    cluster = cluster,
    n_proteins = nrow(cluster_data),
    conditions = conditions,
    zscores = zscores,
    centroid = centroid,
    memberships = cluster_data$MaxMembership,
    feature_ids = cluster_data$FeatureID,
    min_membership = min_membership
  )
}


#' Gráfico de Perfil de Cluster con Highcharts
#'
#' Genera un gráfico interactivo mostrando los perfiles de expresión
#' de las proteínas en un cluster específico.
#'
#' @param result Resultado de pattern_profiler()
#' @param cluster Número de cluster a visualizar
#' @param min_membership Filtro de membership mínima
#' @param color_mode Modo de coloración: "cluster" (color del cluster),
#'   "zscore" (gradiente por z-score), o "membership" (gradiente por membership)
#' @param show_centroid Mostrar línea central (default: TRUE)
#' @param centroid_summary Método para centroide: "mean" o "median"
#' @param palette_gradient Vector de 3 colores para gradiente c(low, mid, high)
#'   Solo usado cuando color_mode = "zscore" o "membership"
#' @param cluster_color Color base para el cluster (opcional, se auto-asigna)
#' @param line_width Ancho de líneas de perfil (default: 1.5)
#' @param centroid_width Ancho de línea del centroide (default: 3)
#' @param opacity_range Rango de opacidad c(min, max) basado en membership
#' @param title Título personalizado (opcional)
#' @param height Altura del gráfico en píxeles
#'
#' @return Objeto highchart
cluster_profile_highchart <- function(result,
                                       cluster,
                                       min_membership = NULL,
                                       color_mode = c("cluster", "zscore", "membership"),
                                       show_centroid = TRUE,
                                       centroid_summary = c("mean", "median"),
                                       palette_gradient = c("#2166AC", "#DDDDDD", "#B2182B"),
                                       cluster_color = NULL,
                                       line_width = 1.5,
                                       centroid_width = 3,
                                       opacity_range = c(0.3, 0.8),
                                       title = NULL,
                                       height = NULL) {

  color_mode <- match.arg(color_mode)
  centroid_summary <- match.arg(centroid_summary)

  # ---------------------------------------------------------------------------
  # 1) Preparar datos del cluster
  # ---------------------------------------------------------------------------
  profile_data <- prepare_cluster_profile_data(result, cluster, min_membership)

  if (is.null(profile_data)) {
    warning(sprintf("Cluster %d: sin proteínas con membership >= %.2f",
                    cluster, min_membership %||% result$min_membership))
    return(NULL)
  }

  conditions <- profile_data$conditions
  zscores <- profile_data$zscores
  memberships <- profile_data$memberships
  feature_ids <- profile_data$feature_ids
  n_proteins <- profile_data$n_proteins

  # ---------------------------------------------------------------------------
  # 2) Calcular centroide
  # ---------------------------------------------------------------------------
  if (centroid_summary == "mean") {
    centroid <- colMeans(zscores, na.rm = TRUE)
  } else {
    centroid <- apply(zscores, 2, median, na.rm = TRUE)
  }
  # Convertir a vector sin nombres para evitar warnings de jsonlite
  centroid <- as.numeric(centroid)

  # ---------------------------------------------------------------------------
  # 3) Configurar color del cluster
  # ---------------------------------------------------------------------------
  if (is.null(cluster_color)) {
    cluster_palette <- configure_cluster_palette(result$optimal_c)
    cluster_color <- unname(cluster_palette[cluster])
  }

  # ---------------------------------------------------------------------------
  # 4) Título del gráfico
  # ---------------------------------------------------------------------------
  if (is.null(title)) {
    title <- sprintf("Cluster %d (n = %d, membership >= %.2f)",
                     cluster, n_proteins, profile_data$min_membership)
  }

  # ---------------------------------------------------------------------------
  # 5) Construir series de líneas de perfil
  # ---------------------------------------------------------------------------
  # Calcular rango de memberships para normalización de opacidad
  # Usamos el rango desde min_membership hasta 1 (no el rango observado)
  # para que la opacidad refleje membership absoluta
  mem_min_threshold <- profile_data$min_membership
  mem_range_for_opacity <- 1 - mem_min_threshold

  # Calcular rango de z-scores para modo zscore
  if (color_mode == "zscore") {
    all_zscores <- as.numeric(zscores)
    zscore_q <- quantile(all_zscores, c(0.05, 0.95), na.rm = TRUE)
    zscore_range <- c(
      max(-3, min(-0.5, zscore_q[1])),  # Al menos -0.5 para tener contraste
      min(3, max(0.5, zscore_q[2]))      # Al menos 0.5 para tener contraste
    )
  }

  profile_series <- lapply(seq_len(n_proteins), function(i) {
    zscore_row <- as.numeric(zscores[i, ])
    membership_val <- memberships[i]
    feature_id <- feature_ids[i]

    # Calcular opacidad basada en membership absoluta (desde threshold hasta 1)
    mem_norm <- (membership_val - mem_min_threshold) / max(0.01, mem_range_for_opacity)
    mem_norm <- max(0, min(1, mem_norm))  # Asegurar en [0,1]
    opacity <- opacity_range[1] + mem_norm * (opacity_range[2] - opacity_range[1])

    # Determinar color según modo
    if (color_mode == "cluster") {
      # Usar color del cluster directamente
      line_color <- cluster_color

    } else if (color_mode == "zscore") {
      # Color basado en el z-score del último punto (o el más extremo)
      # Esto muestra la tendencia final del perfil
      last_zscore <- zscore_row[length(zscore_row)]
      # Alternativamente: usar el z-score más extremo
      # extreme_idx <- which.max(abs(zscore_row))
      # extreme_zscore <- zscore_row[extreme_idx]
      line_color <- interpolate_color(
        last_zscore,
        low_color = palette_gradient[1],
        mid_color = palette_gradient[2],
        high_color = palette_gradient[3],
        midpoint = 0,
        limits = zscore_range
      )

    } else {
      # color_mode == "membership"
      line_color <- interpolate_color(
        membership_val,
        low_color = palette_gradient[1],
        mid_color = palette_gradient[2],
        high_color = palette_gradient[3],
        midpoint = (mem_min_threshold + 1) / 2,
        limits = c(mem_min_threshold, 1)
      )
    }

    # Construir puntos de la línea (como lista, no vector nombrado)
    points <- lapply(seq_along(conditions), function(j) {
      list(
        x = as.integer(j - 1),
        y = round(zscore_row[j], 4),
        condition = conditions[j]
      )
    })

    # IMPORTANTE: Para líneas en Highcharts, usar 'opacity' como propiedad de serie
    # además del color (rgba no siempre funciona para líneas)
    list(
      name = feature_id,
      type = "line",
      data = points,
      color = line_color,
      opacity = opacity,
      lineWidth = line_width,
      marker = list(enabled = FALSE),
      enableMouseTracking = TRUE,
      showInLegend = FALSE,
      states = list(
        hover = list(
          lineWidth = line_width + 1.5,
          opacity = 1,
          enabled = TRUE
        )
      )
    )
  })

  # ---------------------------------------------------------------------------
  # 6) Serie del centroide
  # ---------------------------------------------------------------------------
  centroid_series <- NULL
  if (show_centroid) {
    centroid_points <- lapply(seq_along(conditions), function(j) {
      list(
        x = as.integer(j - 1),
        y = round(centroid[j], 4),
        condition = conditions[j]
      )
    })

    centroid_series <- list(
      name = paste0("Centroid (", centroid_summary, ")"),
      type = "line",
      data = centroid_points,
      color = darken_hex(cluster_color, 0.2),
      lineWidth = centroid_width,
      marker = list(
        enabled = TRUE,
        symbol = "circle",
        radius = 5,
        fillColor = cluster_color,
        lineColor = darken_hex(cluster_color, 0.3),
        lineWidth = 2
      ),
      zIndex = 10,
      showInLegend = TRUE
    )
  }

  # ---------------------------------------------------------------------------
  # 7) Construir highchart
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
    hc_subtitle(
      text = sprintf("Color by %s | Opacity by membership",
                     ifelse(color_mode == "zscore", "z-score", "membership")),
      style = list(
        fontSize = "12px",
        color = "#6C757D"
      )
    ) |>
    hc_xAxis(
      categories = conditions,
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
      opacity = series$opacity,
      lineWidth = series$lineWidth,
      marker = series$marker,
      enableMouseTracking = series$enableMouseTracking,
      showInLegend = series$showInLegend,
      states = series$states
    )
  }

  # Añadir centroide si corresponde
  if (!is.null(centroid_series)) {
    hc <- hc |> hc_add_series(
      name = centroid_series$name,
      type = centroid_series$type,
      data = centroid_series$data,
      color = centroid_series$color,
      lineWidth = centroid_series$lineWidth,
      marker = centroid_series$marker,
      zIndex = centroid_series$zIndex,
      showInLegend = centroid_series$showInLegend
    )
  }

  hc
}


#' Lista de Gráficos de Perfil de Clusters con Highcharts
#'
#' Genera una lista de gráficos interactivos para todos los clusters.
#'
#' @param result Resultado de pattern_profiler()
#' @param clusters Vector de clusters a visualizar (NULL = todos)
#' @param min_membership Filtro de membership mínima
#' @param color_mode Modo de coloración: "cluster", "zscore" o "membership"
#' @param show_centroid Mostrar línea central (default: TRUE)
#' @param centroid_summary Método para centroide: "mean" o "median"
#' @param palette_gradient Vector de 3 colores para gradiente (solo zscore/membership)
#' @param palette Paleta para colores de clusters (NULL, "ggsci::", "brewer:")
#' @param line_width Ancho de líneas de perfil
#' @param centroid_width Ancho de línea del centroide
#' @param opacity_range Rango de opacidad basado en membership
#' @param height Altura de cada gráfico en píxeles
#'
#' @return Lista nombrada de objetos highchart
#'
#' @examples
#' \dontrun{
#' # Ejecutar clustering
#' result <- pattern_profiler(
#'   data = data,
#'   filter_mode = "any",
#'   condition_order = c("A", "B", "C", "D")
#' )
#'
#' # Generar gráficos para todos los clusters (color por cluster)
#' hc_profiles <- cluster_profile_highchart_list(result)
#'
#' # Visualizar cluster 1
#' hc_profiles[["Cluster_1"]]
#'
#' # Con gradiente por z-score
#' hc_profiles <- cluster_profile_highchart_list(
#'   result,
#'   color_mode = "zscore",
#'   palette_gradient = c("blue", "gray", "red")
#' )
#' }
cluster_profile_highchart_list <- function(result,
                                            clusters = NULL,
                                            min_membership = NULL,
                                            color_mode = c("cluster", "zscore", "membership"),
                                            show_centroid = TRUE,
                                            centroid_summary = c("mean", "median"),
                                            palette_gradient = c("#2166AC", "#DDDDDD", "#B2182B"),
                                            palette = NULL,
                                            line_width = 1.5,
                                            centroid_width = 3,
                                            opacity_range = c(0.3, 0.8),
                                            height = NULL) {

  color_mode <- match.arg(color_mode)
  centroid_summary <- match.arg(centroid_summary)

  # Determinar clusters a visualizar
  if (is.null(clusters)) {
    clusters <- seq_len(result$optimal_c)
  }

  # Configurar paleta de colores para clusters
  cluster_colors <- configure_cluster_palette(result$optimal_c, palette)

  # Generar gráficos
  hc_list <- lapply(clusters, function(k) {
    hc <- tryCatch({
      cluster_profile_highchart(
        result = result,
        cluster = k,
        min_membership = min_membership,
        color_mode = color_mode,
        show_centroid = show_centroid,
        centroid_summary = centroid_summary,
        palette_gradient = palette_gradient,
        cluster_color = cluster_colors[k],
        line_width = line_width,
        centroid_width = centroid_width,
        opacity_range = opacity_range,
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


#' Gráfico de Resumen de Todos los Clusters (Centroides)
#'
#' Genera un gráfico comparativo mostrando los centroides de todos los clusters.
#'
#' @param result Resultado de pattern_profiler()
#' @param clusters Vector de clusters a incluir (NULL = todos)
#' @param min_membership Filtro de membership mínima
#' @param centroid_summary Método para centroide: "mean" o "median"
#' @param palette Paleta de colores para clusters
#' @param line_width Ancho de líneas
#' @param show_markers Mostrar marcadores en los puntos (default: TRUE)
#' @param title Título personalizado
#' @param height Altura del gráfico
#'
#' @return Objeto highchart
cluster_centroids_highchart <- function(result,
                                         clusters = NULL,
                                         min_membership = NULL,
                                         centroid_summary = c("mean", "median"),
                                         palette = NULL,
                                         line_width = 2.5,
                                         show_markers = TRUE,
                                         title = NULL,
                                         height = NULL) {

  centroid_summary <- match.arg(centroid_summary)

  if (is.null(min_membership)) {
    min_membership <- result$min_membership
  }

  # Determinar clusters
  if (is.null(clusters)) {
    clusters <- seq_len(result$optimal_c)
  }

  # Configurar colores (sin nombres para evitar warnings jsonlite)
  cluster_colors <- unname(configure_cluster_palette(result$optimal_c, palette))

  # Convertir conditions a vector sin nombres
  conditions <- as.character(result$conditions)

  # Título
  if (is.null(title)) {
    title <- sprintf("Cluster Centroids (%s)", centroid_summary)
  }

  # ---------------------------------------------------------------------------
  # Construir highchart base
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
    hc_subtitle(
      text = sprintf("Membership >= %.2f", min_membership),
      style = list(
        fontSize = "12px",
        color = "#6C757D"
      )
    ) |>
    hc_xAxis(
      categories = as.list(conditions),
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
      gridLineColor = "#F1F3F4",
      gridLineDashStyle = "Dot",
      plotLines = list(
        list(value = 0, color = "#ADB5BD", width = 1, dashStyle = "Dash")
      )
    ) |>
    hc_legend(
      enabled = TRUE,
      layout = "horizontal",
      align = "center",
      verticalAlign = "bottom",
      itemStyle = list(fontSize = "12px", fontWeight = "normal", color = "#495057")
    ) |>
    hc_tooltip(
      useHTML = TRUE,
      shared = FALSE,
      headerFormat = "",
      pointFormat = paste0(
        "<div style='padding: 4px;'>",
        "<b style='color: {series.color};'>{series.name}</b><br/>",
        "<span style='color: #6C757D;'>Condition:</span> <b>{point.condition}</b><br/>",
        "<span style='color: #6C757D;'>z-score:</span> <b>{point.y:.3f}</b><br/>",
        "<span style='color: #6C757D;'>n proteins:</span> <b>{point.n}</b>",
        "</div>"
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
  # Añadir serie por cluster
  # ---------------------------------------------------------------------------
  for (k in clusters) {
    profile_data <- prepare_cluster_profile_data(result, k, min_membership)

    if (is.null(profile_data)) next

    # Calcular centroide (convertir a vector numérico sin nombres)
    if (centroid_summary == "mean") {
      centroid <- as.numeric(colMeans(profile_data$zscores, na.rm = TRUE))
    } else {
      centroid <- as.numeric(apply(profile_data$zscores, 2, median, na.rm = TRUE))
    }

    n_prots <- as.integer(profile_data$n_proteins)
    k_color <- cluster_colors[k]

    # Construir puntos con x explícito (formato que Highcharts entiende bien)
    points <- lapply(seq_along(conditions), function(j) {
      list(
        x = as.integer(j - 1),
        y = round(centroid[j], 4),
        condition = conditions[j],
        n = n_prots
      )
    })

    hc <- hc |> hc_add_series(
      name = sprintf("Cluster %d (n=%d)", k, n_prots),
      type = "line",
      data = points,
      color = k_color,
      lineWidth = line_width,
      marker = list(
        enabled = show_markers,
        symbol = "circle",
        radius = 5,
        fillColor = k_color,
        lineColor = darken_hex(k_color, 0.2),
        lineWidth = 1
      )
    )
  }

  hc
}


# =============================================================================
# EJEMPLOS DE USO
# =============================================================================

# --- Cargar datos ---
# data <- arrow::read_parquet("data/PCA_Input.parquet")
# data <- readr::read_tsv("data/PCA_Input.tsv")

# --- Ejecutar clustering con selección automática ---
# result <- pattern_profiler(
#   data = data,
#   filter_mode = "any",
#   condition_order = c("A", "B", "C", "D"),
#   c_range = 2:8,
#   auto_select_c = TRUE,
#   selection_method = "xb",
#   min_membership = 0.5
# )

# --- Ver resultados ---
# result$optimal_c
# result$selection_metrics
# head(result$membership_table)

# --- Generar gráficos de perfiles por cluster (color del cluster) ---
# hc_profiles <- cluster_profile_highchart_list(result)
# hc_profiles[["Cluster_1"]]
# hc_profiles[["Cluster_2"]]

# --- Gráfico de un cluster específico ---
# hc_c1 <- cluster_profile_highchart(
#   result,
#   cluster = 1,
#   color_mode = "cluster",  # "cluster", "zscore", o "membership"
#   min_membership = 0.5
# )
# hc_c1

# --- Gráfico con gradiente por z-score ---
# hc_zscore <- cluster_profile_highchart(
#   result,
#   cluster = 1,
#   color_mode = "zscore",
#   palette_gradient = c("blue", "lightgray", "red")
# )
# hc_zscore

# --- Gráfico de centroides comparativo ---
# hc_centroids <- cluster_centroids_highchart(result)
# hc_centroids

# --- Con paleta de colores personalizada para clusters ---
# hc_profiles <- cluster_profile_highchart_list(
#   result,
#   palette = "ggsci::nrc_npg"
# )

# --- Con gradiente por membership ---
# hc_profiles <- cluster_profile_highchart_list(
#   result,
#   color_mode = "membership",
#   palette_gradient = c("#440154", "#21918c", "#fde725"),  # viridis
#   opacity_range = c(0.4, 0.9)
# )
