# =============================================================================
# Pattern Profiler Analysis: Mfuzz Clustering desde SummarizedExperiment
# =============================================================================
#
# Este script procesa un objeto SummarizedExperiment junto con resultados de
# expresión diferencial (DEPs_results) para realizar soft-clustering con Mfuzz.
#
# Input:
#   - se_proc: SummarizedExperiment con assays de intensidad
#   - DEPs_results: DataFrame con adj.P.Val por comparación
#
# Output:
#   - Pattern_Profiler_Input.parquet: Tabla en formato LONG
#
# Autor: Sergio Ciordia
# Licencia: MIT
# =============================================================================

# -----------------------------------------------------------------------------
# DEPENDENCIAS
# -----------------------------------------------------------------------------





# =============================================================================
# FUNCIONES DE EXTRACCIÓN DE DATOS
# =============================================================================

#' Extraer datos de SummarizedExperiment
#'
#' Obtiene la matriz de intensidades, IDs de features y metadatos de muestras
#' desde un objeto SummarizedExperiment.
#'
#' @param se_proc Objeto SummarizedExperiment
#' @param assay_name Nombre del assay a usar (NULL = primero disponible)
#'
#' @return Lista con:
#'   - intensity_matrix: Matriz features x samples
#'   - feature_ids: Vector de IDs de proteínas
#'   - sample_metadata: DataFrame con SampleID, Condition, etc.
extract_se_data <- function(se_proc, assay_name = NULL) {


  # Validar input

  if (!inherits(se_proc, "SummarizedExperiment")) {
    stop("'se_proc' debe ser un objeto SummarizedExperiment")
  }


  # Obtener nombre del assay

  available_assays <- SummarizedExperiment::assayNames(se_proc)
  if (length(available_assays) == 0) {
    stop("El SummarizedExperiment no tiene assays")
  }

  if (is.null(assay_name)) {
    assay_name <- available_assays[1]
    message(sprintf("Usando assay: '%s'", assay_name))
  } else if (!(assay_name %in% available_assays)) {
    stop(sprintf("Assay '%s' no encontrado. Disponibles: %s",
                 assay_name, paste(available_assays, collapse = ", ")))
  }

  # Extraer matriz de intensidades
  intensity_matrix <- SummarizedExperiment::assay(se_proc, assay_name)

  # Obtener IDs de features desde rowData

  row_data <- as.data.frame(SummarizedExperiment::rowData(se_proc))

  # Buscar columna de ID (varios nombres posibles)
  id_candidates <- c("FeatureID", "Protein.IDs", "ID", "protein_id", "ProteinID")
  id_col <- intersect(id_candidates, names(row_data))[1]

  if (!is.na(id_col)) {
    feature_ids <- as.character(row_data[[id_col]])
  } else {
    # Usar rownames si no hay columna de ID
    feature_ids <- rownames(intensity_matrix)
    if (is.null(feature_ids)) {
      feature_ids <- paste0("Feature_", seq_len(nrow(intensity_matrix)))
    }
  }
  rownames(intensity_matrix) <- feature_ids

  # Obtener metadatos de muestras desde colData
  sample_metadata <- as.data.frame(SummarizedExperiment::colData(se_proc))
  sample_metadata$SampleID <- rownames(sample_metadata)

  # Validar que existe columna Condition

  if (!("Condition" %in% names(sample_metadata))) {
    # Intentar encontrar alternativas
    cond_candidates <- c("condition", "Group", "group", "Treatment", "treatment")
    cond_col <- intersect(cond_candidates, names(sample_metadata))[1]

    if (!is.na(cond_col)) {
      sample_metadata$Condition <- sample_metadata[[cond_col]]
      message(sprintf("Usando '%s' como columna de condición", cond_col))
    } else {
      stop("No se encontró columna 'Condition' en colData")
    }
  }

  list(
    intensity_matrix = intensity_matrix,
    feature_ids = feature_ids,
    sample_metadata = sample_metadata,
    assay_name = assay_name
  )
}


#' Fusionar información de significancia desde DEPs_results
#'
#' Convierte DEPs_results (formato largo) a formato ancho con columnas
#' adjP_* por cada comparación, y calcula sig_any.
#'
#' @param feature_ids Vector de IDs de features
#' @param DEPs_results DataFrame con columnas Protein.IDs, adj.P.Val, Comparison, Assay
#' @param assay_name Nombre del assay para filtrar DEPs_results
#' @param alpha Umbral de significancia (default: 0.05)
#'
#' @return DataFrame con FeatureID, adjP_*, sig_any
merge_significance_info <- function(feature_ids, DEPs_results, assay_name, alpha = 0.05) {

  DEPs_results <- as.data.frame(DEPs_results)

  # Validar columnas requeridas
  required_cols <- c("Protein.IDs", "adj.P.Val", "Comparison", "Assay")
  missing_cols <- setdiff(required_cols, names(DEPs_results))
  if (length(missing_cols) > 0) {
    stop("Columnas faltantes en DEPs_results: ", paste(missing_cols, collapse = ", "))
  }

  # Filtrar por Assay
  available_assays <- unique(DEPs_results$Assay)
  if (!(assay_name %in% available_assays)) {
    stop(sprintf("Assay '%s' no encontrado en DEPs_results. Disponibles: %s",
                 assay_name, paste(available_assays, collapse = ", ")))
  }

  DEPs_results <- DEPs_results[DEPs_results$Assay == assay_name, ]
  message(sprintf("   - Filtrado DEPs_results por Assay = '%s' (%d filas)",
                  assay_name, nrow(DEPs_results)))

  # Pivotar a formato ancho: una fila por proteína, columnas adjP_* por comparación
  sig_wide <- DEPs_results %>%
    dplyr::select(Protein.IDs, Comparison, adj.P.Val) %>%
    dplyr::distinct() %>%
    tidyr::pivot_wider(
      names_from = Comparison,
      values_from = adj.P.Val,
      names_prefix = "adjP_"
    )

  # Crear DataFrame base con todos los feature_ids

  feature_info <- data.frame(
    FeatureID = feature_ids,
    stringsAsFactors = FALSE
  )

  # Fusionar con significancia
  feature_info <- merge(
    feature_info,
    sig_wide,
    by.x = "FeatureID",
    by.y = "Protein.IDs",
    all.x = TRUE
  )

  # Calcular sig_any (significativo en cualquier comparación)
  adjP_cols <- grep("^adjP_", names(feature_info), value = TRUE)

  if (length(adjP_cols) > 0) {
    feature_info$sig_any <- apply(
      feature_info[, adjP_cols, drop = FALSE], 1,
      function(x) any(x <= alpha, na.rm = TRUE)
    )
  } else {
    feature_info$sig_any <- FALSE
    warning("No se encontraron columnas adjP_* en DEPs_results")
  }

  # Reemplazar NA en sig_any con FALSE
  feature_info$sig_any[is.na(feature_info$sig_any)] <- FALSE

  feature_info
}


#' Filtrar features por significancia
#'
#' @param feature_info DataFrame con columnas de significancia
#' @param filter_mode Modo de selección:
#'   - "any": features significativos en AL MENOS una comparación (usa sig_any).
#'   - "all": TODOS los features, sin filtrar por significancia (no es "significativo
#'            en todas las comparaciones").
#'   - "specific": significativos en la comparación indicada por `comparison`.
#' @param comparison Comparación específica (para mode="specific")
#' @param alpha Umbral de significancia
#'
#' @return Vector de FeatureIDs seleccionados
filter_significant_features <- function(feature_info,
                                         filter_mode = c("any", "all", "specific"),
                                         comparison = NULL,
                                         alpha = 0.05) {

  filter_mode <- match.arg(filter_mode)

  if (filter_mode == "all") {
    return(feature_info$FeatureID)
  }

  if (filter_mode == "any") {
    if (!("sig_any" %in% names(feature_info))) {
      stop("Columna 'sig_any' no encontrada. Ejecutar merge_significance_info() primero.")
    }
    selected <- feature_info$FeatureID[feature_info$sig_any == TRUE]
    message(sprintf("Filtro 'any': %d de %d features significativos",
                    length(selected), nrow(feature_info)))
    return(selected)
  }

  # filter_mode == "specific"
  if (is.null(comparison)) {
    stop("El argumento 'comparison' es requerido para filter_mode = 'specific'")
  }

  adjP_col <- paste0("adjP_", comparison)
  if (!(adjP_col %in% names(feature_info))) {
    available <- grep("^adjP_", names(feature_info), value = TRUE)
    stop(sprintf("Comparación '%s' no encontrada. Disponibles: %s",
                 comparison, paste(gsub("^adjP_", "", available), collapse = ", ")))
  }

  selected <- feature_info$FeatureID[
    !is.na(feature_info[[adjP_col]]) & feature_info[[adjP_col]] <= alpha
  ]
  message(sprintf("Filtro específico '%s' (alpha=%.3f): %d features",
                  comparison, alpha, length(selected)))

  selected
}


#' Construir matriz de clustering agregada por condición
#'
#' @param intensity_matrix Matriz de intensidades (features x samples)
#' @param sample_metadata DataFrame con SampleID, Condition
#' @param selected_features Vector de features a incluir
#' @param condition_order Orden de condiciones (NULL = orden alfabético)
#' @param aggregate Método de agregación: "median" o "mean"
#'
#' @return Matriz de intensidades agregadas (features x conditions)
build_clustering_matrix <- function(intensity_matrix,
                                     sample_metadata,
                                     selected_features,
                                     condition_order = NULL,
                                     aggregate = c("median", "mean")) {

  aggregate <- match.arg(aggregate)
  agg_fun <- if (aggregate == "median") median else mean

  # Filtrar features seleccionados
  intensity_matrix <- intensity_matrix[selected_features, , drop = FALSE]

  # Obtener condiciones únicas
  conditions <- unique(sample_metadata$Condition)

  if (is.null(condition_order)) {
    condition_order <- sort(conditions)
    message(sprintf("Orden de condiciones: %s", paste(condition_order, collapse = " -> ")))
  } else {
    # Validar que todas las condiciones existen
    missing <- setdiff(condition_order, conditions)
    if (length(missing) > 0) {
      stop("Condiciones no encontradas: ", paste(missing, collapse = ", "))
    }
  }

  # Agregar por condición
  agg_matrix <- matrix(
    NA_real_,
    nrow = length(selected_features),
    ncol = length(condition_order),
    dimnames = list(selected_features, condition_order)
  )

  for (cond in condition_order) {
    samples_in_cond <- sample_metadata$SampleID[sample_metadata$Condition == cond]
    samples_in_cond <- intersect(samples_in_cond, colnames(intensity_matrix))

    if (length(samples_in_cond) == 0) {
      warning(sprintf("No hay muestras para condición '%s'", cond))
      next
    }

    agg_matrix[, cond] <- apply(
      intensity_matrix[, samples_in_cond, drop = FALSE], 1,
      function(x) agg_fun(x, na.rm = TRUE)
    )
  }

  agg_matrix
}


# =============================================================================
# FUNCIONES DE CLUSTERING (Mfuzz)
# =============================================================================

#' Crear ExpressionSet para Mfuzz
#'
#' @param mat Matriz de expresión (features x conditions)
#' @param feature_info DataFrame con info de features (opcional)
#'
#' @return Objeto ExpressionSet
create_expression_set <- function(mat, feature_info = NULL) {

  # Asegurar que es matriz
  mat <- as.matrix(mat)

  # Crear AnnotatedDataFrame para features
  if (!is.null(feature_info) && nrow(feature_info) == nrow(mat)) {
    rownames(feature_info) <- rownames(mat)
    fData <- Biobase::AnnotatedDataFrame(data = feature_info)
  } else {
    fData <- Biobase::AnnotatedDataFrame(data = data.frame(
      FeatureID = rownames(mat),
      row.names = rownames(mat)
    ))
  }

  # Crear AnnotatedDataFrame para condiciones
  pData <- Biobase::AnnotatedDataFrame(data = data.frame(
    Condition = colnames(mat),
    row.names = colnames(mat)
  ))

  # Crear ExpressionSet
  Biobase::ExpressionSet(
    assayData = mat,
    featureData = fData,
    phenoData = pData
  )
}


#' Estandarizar ExpressionSet (z-score por fila)
#'
#' @param eset Objeto ExpressionSet
#'
#' @return ExpressionSet estandarizado
standardize_eset <- function(eset) {

  # Filtrar filas con muchos NAs
  eset_filtered <- Mfuzz::filter.NA(eset, thres = 0.25)

  n_removed <- nrow(eset) - nrow(eset_filtered)
  if (n_removed > 0) {
    message(sprintf("Eliminadas %d features con >25%% NAs", n_removed))
  }

  # Imputar NAs restantes
  eset_filled <- Mfuzz::fill.NA(eset_filtered, mode = "knn")

  # Estandarizar por filas (z-score)
  eset_std <- Mfuzz::standardise(eset_filled)

  # Descartar features constantes: sd = 0 produce (x - media)/0 = NaN tras
  # standardise, y esas filas romperian o degenerarian mfuzz. Se detectan por
  # filas no finitas en la matriz estandarizada.
  X_std <- Biobase::exprs(eset_std)
  finite_rows <- apply(X_std, 1, function(r) all(is.finite(r)))
  n_const <- sum(!finite_rows)
  if (n_const > 0) {
    message(sprintf("Eliminadas %d features constantes (sd = 0) tras estandarizar",
                    n_const))
    eset_std <- eset_std[finite_rows, ]
  }

  eset_std
}


# -----------------------------------------------------------------------------
# Métricas de evaluación de clusters
# -----------------------------------------------------------------------------

#' Índice Xie-Beni
#'
#' Mide compacidad frente a separación entre clusters. Menor es mejor.
#'
#' @param X Matriz de datos (features x condiciones).
#' @param U Matriz de pertenencias (features x clusters).
#' @param centers Matriz de centroides (clusters x condiciones).
#' @param m Exponente de difuminado (fuzzifier).
#' @return Valor numérico del índice.
#' @keywords internal
#' @noRd
.xie_beni_index <- function(X, U, centers, m = 2) {

  n <- nrow(X)
  c <- nrow(centers)

  # Compactness: suma ponderada de distancias intra-cluster
  compactness <- 0
  for (i in seq_len(n)) {
    for (j in seq_len(c)) {
      dist_sq <- sum((X[i, ] - centers[j, ])^2)
      compactness <- compactness + (U[i, j]^m) * dist_sq
    }
  }

  # Separation: mínima distancia entre centroides
  min_sep <- Inf
  for (j1 in seq_len(c - 1)) {
    for (j2 in (j1 + 1):c) {
      sep <- sum((centers[j1, ] - centers[j2, ])^2)
      if (sep < min_sep) min_sep <- sep
    }
  }

  if (min_sep == 0) min_sep <- .Machine$double.eps

  xb <- compactness / (n * min_sep)
  xb
}


#' Fuzzy Partition Coefficient (FPC)
#'
#' Mide la nitidez de la partición difusa. Mayor es mejor.
#'
#' @param U Matriz de pertenencias (features x clusters).
#' @return Valor numérico del coeficiente.
#' @keywords internal
#' @noRd
.fpc_index <- function(U) {
  n <- nrow(U)
  sum(U^2) / n
}


#' Average Maximum Membership (AMM)
#'
#' Promedio de la pertenencia máxima de cada feature. Mayor es mejor.
#'
#' @param U Matriz de pertenencias (features x clusters).
#' @return Valor numérico del promedio.
#' @keywords internal
#' @noRd
.amm_index <- function(U) {
  mean(apply(U, 1, max))
}


#' Evaluar rango de números de clusters
#'
#' @param eset_std ExpressionSet estandarizado
#' @param c_range Vector de números de clusters a evaluar
#' @param m Parámetro de fuzziness
#' @param seeds Semillas para reproducibilidad
#' @param verbose Mostrar progreso
#'
#' @return DataFrame con métricas por número de clusters
evaluate_cluster_range <- function(eset_std,
                                    c_range = 2:10,
                                    m = NULL,
                                    seeds = c(42, 123, 456),
                                    verbose = TRUE) {

  # Estimar m si no se proporciona
  if (is.null(m)) {
    m <- Mfuzz::mestimate(eset_std)
    if (verbose) message(sprintf("Parámetro m estimado: %.3f", m))
  }

  X <- Biobase::exprs(eset_std)
  results <- list()

  # Se fija una semilla por réplica dentro del bucle; el RNG del usuario se
  # restaura al salir de la función.
  old_rng <- .rng_state()
  on.exit(.rng_restore(old_rng), add = TRUE)

  for (c in c_range) {
    if (verbose) message(sprintf("Evaluando c = %d...", c))

    # NA (no 0): una semilla que falla no debe contar como 0 en la media,
    # porque XB se minimiza y un 0 espurio sesgaria la seleccion de c.
    xb_vals <- rep(NA_real_, length(seeds))
    fpc_vals <- rep(NA_real_, length(seeds))
    amm_vals <- rep(NA_real_, length(seeds))
    dmin_vals <- rep(NA_real_, length(seeds))

    for (s in seq_along(seeds)) {
      set.seed(seeds[s])

      cl <- tryCatch({
        Mfuzz::mfuzz(eset_std, c = c, m = m)
      }, error = function(e) NULL)

      if (is.null(cl)) next

      U <- cl$membership
      centers <- cl$centers

      xb_vals[s] <- .xie_beni_index(X, U, centers, m)
      fpc_vals[s] <- .fpc_index(U)
      amm_vals[s] <- .amm_index(U)

      # Dmin: distancia mínima entre centroides
      dists <- as.matrix(dist(centers))
      diag(dists) <- Inf
      dmin_vals[s] <- min(dists)
    }

    results[[as.character(c)]] <- data.frame(
      c = c,
      XB = mean(xb_vals, na.rm = TRUE),
      XB_sd = sd(xb_vals, na.rm = TRUE),
      FPC = mean(fpc_vals, na.rm = TRUE),
      AMM = mean(amm_vals, na.rm = TRUE),
      Dmin = mean(dmin_vals, na.rm = TRUE)
    )
  }

  do.call(rbind, results)
}


#' Seleccionar número óptimo de clusters
#'
#' @param eset_std ExpressionSet estandarizado
#' @param c_range Rango de clusters a evaluar
#' @param m Parámetro de fuzziness
#' @param method Método: "xb", "consensus", "elbow"
#' @param verbose Mostrar progreso
#'
#' @return Lista con optimal_c, metrics, m
select_optimal_clusters <- function(eset_std,
                                     c_range = 2:10,
                                     m = NULL,
                                     method = c("xb", "consensus", "elbow"),
                                     verbose = TRUE) {

  method <- match.arg(method)

  if (is.null(m)) {
    m <- Mfuzz::mestimate(eset_std)
  }

  metrics <- evaluate_cluster_range(eset_std, c_range, m, verbose = verbose)

  if (method == "xb") {
    # Mínimo Xie-Beni
    optimal_c <- metrics$c[which.min(metrics$XB)]

  } else if (method == "elbow") {
    # Método del codo sobre Dmin
    dmin_diff <- -diff(metrics$Dmin)
    elbow_idx <- which.max(dmin_diff) + 1
    optimal_c <- metrics$c[min(elbow_idx, nrow(metrics))]

  } else {
    # Consensus: promedio de rankings
    metrics$rank_XB <- rank(metrics$XB)
    metrics$rank_FPC <- rank(-metrics$FPC)
    metrics$rank_AMM <- rank(-metrics$AMM)
    metrics$rank_Dmin <- rank(-metrics$Dmin)
    metrics$avg_rank <- (metrics$rank_XB + metrics$rank_FPC +
                          metrics$rank_AMM + metrics$rank_Dmin) / 4
    optimal_c <- metrics$c[which.min(metrics$avg_rank)]
  }

  if (verbose) {
    message(sprintf("\nMétodo '%s': c óptimo = %d", method, optimal_c))
  }

  list(
    optimal_c = optimal_c,
    metrics = metrics,
    m = m
  )
}


#' Ejecutar clustering Mfuzz
#'
#' @param eset_std ExpressionSet estandarizado
#' @param c Número de clusters
#' @param m Parámetro de fuzziness
#' @param seed Semilla para reproducibilidad
#'
#' @return Objeto de clustering Mfuzz
run_mfuzz_clustering <- function(eset_std, c, m, seed = 42) {

  # El RNG del usuario se restaura al salir (mfuzz depende de la semilla).
  old_rng <- .rng_state()
  on.exit(.rng_restore(old_rng), add = TRUE)

  set.seed(seed)
  cl <- Mfuzz::mfuzz(eset_std, c = c, m = m)

  cl
}


# =============================================================================
# CONSTRUCCIÓN DE SALIDA EN FORMATO LONG
# =============================================================================

#' Construir tabla de salida en formato LONG
#'
#' Crea una fila por cada combinación FeatureID-Cluster donde
#' membership >= min_membership. Esto permite que una proteína
#' aparezca en múltiples clusters (soft-clustering).
#'
#' @param cl Objeto de clustering Mfuzz
#' @param eset_std ExpressionSet estandarizado (con z-scores)
#' @param conditions Vector de nombres de condiciones
#' @param min_membership Umbral mínimo de membership
#'
#' @return DataFrame en formato long
build_long_output <- function(cl, eset_std, conditions, min_membership) {

  mem_matrix <- cl$membership  # n_proteins x n_clusters
  zscores <- Biobase::exprs(eset_std)  # n_proteins x n_conditions

  n_proteins <- nrow(mem_matrix)
  n_clusters <- ncol(mem_matrix)

  # Preallocar lista para eficiencia

  long_rows <- vector("list", n_proteins * n_clusters)
  row_idx <- 0

  for (i in seq_len(n_proteins)) {
    feature_id <- rownames(mem_matrix)[i]
    memberships <- mem_matrix[i, ]

    # Encontrar clusters donde membership >= umbral
    qualifying_clusters <- which(memberships >= min_membership)

    for (k in qualifying_clusters) {
      row_idx <- row_idx + 1

      # Crear fila base
      row_data <- list(
        FeatureID = feature_id,
        Cluster = as.integer(k),
        Membership = round(memberships[k], 6)
      )

      # Añadir z-scores por condición
      for (cond in conditions) {
        row_data[[cond]] <- round(zscores[i, cond], 6)
      }

      long_rows[[row_idx]] <- as.data.frame(row_data, stringsAsFactors = FALSE)
    }
  }

  # Combinar todas las filas
  result <- do.call(rbind, long_rows[seq_len(row_idx)])

  # Ordenar por Cluster y Membership descendente
  result <- result[order(result$Cluster, -result$Membership), ]
  rownames(result) <- NULL

  result
}


# =============================================================================
# FUNCIÓN PRINCIPAL: PATTERN PROFILER ANALYSIS
# =============================================================================

#' Pattern Profiler Analysis
#'
#' Pipeline completo de clustering desde SummarizedExperiment.
#' Genera archivo parquet en formato LONG para visualización.
#'
#' @param se_proc SummarizedExperiment con datos de intensidad
#' @param DEPs_results DataFrame con resultados de expresión diferencial (debe tener columna 'Assay')
#' @param assay_name Nombre del assay a usar (default: "LoessCyc"). Se usa para filtrar DEPs_results también.
#' @param filter_mode Modo de filtrado: "any" (signif. en alguna comparación),
#'   "all" (TODOS los features, sin filtrar por significancia), "specific"
#'   (signif. en la comparación de `comparison`)
#' @param alpha Umbral de significancia (default: 0.05)
#' @param comparison Comparación específica (para filter_mode="specific")
#' @param condition_order Orden de condiciones para los perfiles
#' @param aggregate Método de agregación: "median" o "mean"
#' @param c_range Rango de clusters a evaluar
#' @param auto_select_c Selección automática de clusters (default: TRUE)
#' @param c Número fijo de clusters (si auto_select_c=FALSE)
#' @param selection_method Método de selección: "xb", "consensus", "elbow"
#' @param min_membership Umbral mínimo de membership para incluir en salida
#' @param output_file Ruta del archivo parquet de salida. Por defecto `NULL`,
#'   que no escribe nada en disco; los datos en formato largo se devuelven
#'   igualmente en el elemento `long_output` del resultado.
#' @param verbose Mostrar mensajes de progreso
#'
#' @return Lista (invisible) con los resultados del clustering: `optimal_c`, `m`,
#'   `conditions`, recuentos de features, `selection_metrics`, `cluster_counts`,
#'   el objeto `cl` de Mfuzz, el `eset_std` estandarizado y `long_output`, el
#'   data.frame en formato largo que se escribe cuando se indica `output_file`.
#'
#' @examples
#' \dontrun{
#' result <- pattern_profiler_analysis(
#'   se_proc = se_proc,
#'   DEPs_results = DEPs_results,
#'   filter_mode = "any",
#'   condition_order = c("A", "B", "C", "D"),
#'   min_membership = 0.25,
#'   output_file = "data-raw/Pattern_Profiler_Input.parquet"
#' )
#' }
#' @export
pattern_profiler_analysis <- function(se_proc,
                                       DEPs_results,
                                       assay_name = "LoessCyc",
                                       filter_mode = c("any", "all", "specific"),
                                       alpha = 0.05,
                                       comparison = NULL,
                                       condition_order = NULL,
                                       aggregate = c("median", "mean"),
                                       c_range = 2:10,
                                       auto_select_c = TRUE,
                                       c = NULL,
                                       selection_method = c("xb", "consensus", "elbow"),
                                       min_membership = 0.25,
                                       output_file = NULL,
                                       verbose = TRUE) {

  filter_mode <- match.arg(filter_mode)
  aggregate <- match.arg(aggregate)
  selection_method <- match.arg(selection_method)

  # Mfuzz llama a exprs() y cmeans() sin cualificar, así que necesita Biobase y
  # e1071 adjuntados en la ruta de búsqueda. Se adjuntan aquí y se sueltan al
  # salir, de modo que la sesión del usuario queda como estaba.
  .pp_adjuntados <- .mfuzz_deps_attach()
  on.exit(.mfuzz_deps_detach(.pp_adjuntados), add = TRUE)

  if (verbose) message("=== Pattern Profiler Analysis ===\n")

  # -------------------------------------------------------------------------
  # 1) Extraer datos del SummarizedExperiment
  # -------------------------------------------------------------------------
  if (verbose) message("1. Extrayendo datos del SummarizedExperiment...")

  se_data <- extract_se_data(se_proc, assay_name)

  if (verbose) {
    message(sprintf("   - Features: %d", length(se_data$feature_ids)))
    message(sprintf("   - Muestras: %d", nrow(se_data$sample_metadata)))
    message(sprintf("   - Condiciones: %s",
                    paste(unique(se_data$sample_metadata$Condition), collapse = ", ")))
  }

  # -------------------------------------------------------------------------
  # 2) Fusionar información de significancia
  # -------------------------------------------------------------------------
  if (verbose) message("\n2. Fusionando información de significancia...")

  feature_info <- merge_significance_info(
    se_data$feature_ids,
    DEPs_results,
    assay_name = se_data$assay_name,
    alpha = alpha
  )

  n_sig <- sum(feature_info$sig_any, na.rm = TRUE)
  if (verbose) {
    message(sprintf("   - Features significativos (alpha=%.3f): %d", alpha, n_sig))
  }

  # -------------------------------------------------------------------------
  # 3) Filtrar features
  # -------------------------------------------------------------------------
  if (verbose) message("\n3. Filtrando features...")

  selected_features <- filter_significant_features(
    feature_info,
    filter_mode = filter_mode,
    comparison = comparison,
    alpha = alpha
  )

  if (length(selected_features) < 10) {
    stop("Muy pocos features seleccionados (< 10). Ajustar filtros.")
  }

  # -------------------------------------------------------------------------
  # 4) Construir matriz de clustering
  # -------------------------------------------------------------------------
  if (verbose) message("\n4. Construyendo matriz de clustering...")

  clustering_matrix <- build_clustering_matrix(
    se_data$intensity_matrix,
    se_data$sample_metadata,
    selected_features,
    condition_order,
    aggregate
  )

  # Guardar orden de condiciones
  conditions <- colnames(clustering_matrix)

  if (verbose) {
    message(sprintf("   - Matriz: %d features x %d condiciones",
                    nrow(clustering_matrix), ncol(clustering_matrix)))
  }

  # -------------------------------------------------------------------------
  # 5) Crear y estandarizar ExpressionSet
  # -------------------------------------------------------------------------
  if (verbose) message("\n5. Estandarizando datos (z-score)...")

  eset <- create_expression_set(clustering_matrix)
  eset_std <- standardize_eset(eset)

  n_final <- nrow(eset_std)
  if (verbose) {
    message(sprintf("   - Features finales (post-filtro NA): %d", n_final))
  }

  # -------------------------------------------------------------------------
  # 6) Seleccionar número de clusters
  # -------------------------------------------------------------------------
  if (auto_select_c) {
    if (verbose) message("\n6. Seleccionando número óptimo de clusters...")

    selection <- select_optimal_clusters(
      eset_std,
      c_range = c_range,
      method = selection_method,
      verbose = verbose
    )

    optimal_c <- selection$optimal_c
    m <- selection$m
    selection_metrics <- selection$metrics

  } else {
    if (is.null(c)) {
      stop("Debe especificar 'c' cuando auto_select_c = FALSE")
    }
    optimal_c <- c
    m <- Mfuzz::mestimate(eset_std)
    selection_metrics <- NULL

    if (verbose) {
      message(sprintf("\n6. Usando c = %d (fijo), m = %.3f", optimal_c, m))
    }
  }

  # -------------------------------------------------------------------------
  # 7) Ejecutar clustering
  # -------------------------------------------------------------------------
  if (verbose) message(sprintf("\n7. Ejecutando Mfuzz clustering (c=%d)...", optimal_c))

  cl <- run_mfuzz_clustering(eset_std, c = optimal_c, m = m)

  # Contar proteínas por cluster (asignación hard)
  hard_assignment <- apply(cl$membership, 1, which.max)
  cluster_counts <- table(hard_assignment)

  if (verbose) {
    message("   Distribución de clusters (asignación hard):")
    for (k in seq_len(optimal_c)) {
      count <- ifelse(as.character(k) %in% names(cluster_counts),
                      cluster_counts[[as.character(k)]], 0)
      message(sprintf("     Cluster %d: %d proteínas", k, count))
    }
  }

  # -------------------------------------------------------------------------
  # 8) Construir salida en formato LONG
  # -------------------------------------------------------------------------
  if (verbose) message(sprintf("\n8. Construyendo tabla LONG (membership >= %.2f)...",
                               min_membership))

  long_output <- build_long_output(cl, eset_std, conditions, min_membership)

  n_unique_features <- length(unique(long_output$FeatureID))
  n_rows <- nrow(long_output)
  multi_cluster <- n_rows - n_unique_features

  if (verbose) {
    message(sprintf("   - Filas totales: %d", n_rows))
    message(sprintf("   - Features únicos: %d", n_unique_features))
    message(sprintf("   - Features en múltiples clusters: %d", multi_cluster))
  }

  # -------------------------------------------------------------------------
  # 9) Guardar archivo parquet (solo si se ha indicado una ruta)
  # -------------------------------------------------------------------------
  # Con output_file = NULL no se escribe nada: la función no debe crear archivos
  # en el espacio de trabajo del usuario sin que este indique dónde. El
  # data.frame en formato largo se devuelve igualmente en el resultado.
  if (is.null(output_file)) {
    if (verbose) message("\n9. Sin output_file: no se guarda ningún archivo")
  } else {
    if (verbose) message(sprintf("\n9. Guardando archivo: %s", output_file))

    output_dir <- dirname(output_file)
    if (output_dir != "." && !dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }

    arrow::write_parquet(long_output, output_file)
  }

  if (verbose) message("\n=== Análisis completado ===")

  # -------------------------------------------------------------------------
  # Retornar resultados (invisiblemente)
  # -------------------------------------------------------------------------
  result <- list(
    optimal_c = optimal_c,
    m = m,
    conditions = conditions,
    n_features_input = length(selected_features),
    n_features_final = n_final,
    n_rows_output = n_rows,
    min_membership = min_membership,
    selection_metrics = selection_metrics,
    cluster_counts = as.data.frame(cluster_counts),
    output_file = output_file,
    long_output = long_output,
    cl = cl,
    eset_std = eset_std
  )

  class(result) <- c("pattern_profiler_result", "list")

  invisible(result)
}


# =============================================================================
# EJEMPLOS DE USO
# =============================================================================

# --- Uso típico ---
# Los objetos 'se_proc' (SummarizedExperiment) y 'DEPs_results' (dataframe)
# ya están cargados en el environment desde pasos previos del pipeline.
#
# source("R/Pattern_Profiler_Analysis.R")
#
# # Ejecutar análisis (usa assay 'LoessCyc' por defecto)
# result <- pattern_profiler_analysis(
#   se_proc = se_proc,
#   DEPs_results = DEPs_results,
#   assay_name = "LoessCyc",  # default, filtra también DEPs_results por esta columna
#   filter_mode = "any",
#   condition_order = c("A", "B", "C", "D"),
#   c_range = 2:8,
#   selection_method = "xb",
#   min_membership = 0.25,
#   output_file = "data-raw/Pattern_Profiler_Input.parquet"
# )
#
# # Ver resultados
# result$optimal_c
# result$selection_metrics
#
# --- Usar un assay diferente ---
# result <- pattern_profiler_analysis(
#   se_proc = se_proc,
#   DEPs_results = DEPs_results,
#   assay_name = "log2",  # Usar log2 en lugar de LoessCyc
#   filter_mode = "any",
#   condition_order = c("A", "B", "C", "D")
# )
