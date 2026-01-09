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
#' @param low_color Color para valores bajos
#' @param mid_color Color para valores medios
#' @param high_color Color para valores altos
#' @param midpoint Punto medio del gradiente
#' @param limits Vector c(min, max) para los límites
#' @return Color hex interpolado
interpolate_color <- function(value, low_color = "#2166AC", mid_color = "#F7F7F7",
                               high_color = "#B2182B", midpoint = 0,
                               limits = c(-3, 3)) {
  # Normalizar valor a rango 0-1
  value <- max(limits[1], min(limits[2], value))

  if (value <= midpoint) {
    t <- (value - limits[1]) / (midpoint - limits[1])
    col1 <- col2rgb(low_color) / 255
    col2 <- col2rgb(mid_color) / 255
  } else {
    t <- (value - midpoint) / (limits[2] - midpoint)
    col1 <- col2rgb(mid_color) / 255
    col2 <- col2rgb(high_color) / 255
  }

  r <- round((col1[1] + t * (col2[1] - col1[1])) * 255)
  g <- round((col1[2] + t * (col2[2] - col1[2])) * 255)
  b <- round((col1[3] + t * (col2[3] - col1[3])) * 255)

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
# EJEMPLOS DE USO (comentados)
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

# --- Ejecutar con número fijo de clusters ---
# result <- pattern_profiler(
#   data = data,
#   filter_mode = "any",
#   condition_order = c("A", "B", "C", "D"),
#   auto_select_c = FALSE,
#   c = 4,
#   min_membership = 0.7
# )

# --- Ver resultados ---
# result$optimal_c
# result$selection_metrics
# head(result$membership_table)
