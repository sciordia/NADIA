# =============================================================================
# Procesamiento de Datos Proteómicos
# =============================================================================
#
# Funciones para procesar datos proteómicos de Spectronaut:
#   - Filtrado por presencia en grupos
#   - Normalización (Cyclic Loess)
#   - Imputación mixta (MAR + MNAR)
#   - Análisis diferencial (limma)
#
# Dependencias requeridas:
#   - SummarizedExperiment, S4Vectors
#   - limma
#   - readr (opcional, para exportación)
#
# Dependencias opcionales:
#   - rrcovNA (para imputación MAR con impSeqRob)
#
# Nota: La imputación MNAR (método "min") NO requiere dependencias externas.
#
# Autor: Sergio Ciordia
# Licencia: MIT
# =============================================================================

# --- Operador %||% (coalescencia nula) ---
`%||%` <- function(a, b) if (is.null(a)) b else a

# =============================================================================
# FUNCIONES LINKER (internas)
# =============================================================================

#' Prepara metadata desde el objeto spectronaut_data
#'
#' @param preprocessing Lista de clase spectronaut_data
#' @return Data frame con columnas: Column, Condition, Replicate
#' @keywords internal
.prepare_metadata <- function(preprocessing) {
  stopifnot(inherits(preprocessing, "spectronaut_data"))


  md <- preprocessing$metadata

  # Renombrar columnas al formato esperado
  result <- data.frame(
    Column = md$Coding,
    Condition = md$R.Condition,
    Replicate = md$R.Replicate,
    stringsAsFactors = FALSE
  )

  rownames(result) <- result$Column
  result
}

#' Prepara datos de proteínas desde el objeto spectronaut_data
#'
#' @param preprocessing Lista de clase spectronaut_data
#' @return Data frame con columnas: ProteinGroups, GeneNames, UniqPepts, y columnas de intensidad
#' @keywords internal
.prepare_protein_data <- function(preprocessing) {
  stopifnot(inherits(preprocessing, "spectronaut_data"))

  pq <- preprocessing$protein_quant

  # Identificar columnas de cantidad (PG.Quantity_*)
  quantity_cols <- grep("^PG\\.Quantity_", names(pq), value = TRUE)
  if (length(quantity_cols) == 0) {
    stop("No se encontraron columnas PG.Quantity_* en protein_quant")
  }

  # Identificar columnas de péptidos únicos (PG.NrOfStrippedSequencesUsedForQuantification_*)
  pept_cols <- grep("^PG\\.NrOfStrippedSequencesUsedForQuantification_", names(pq), value = TRUE)

  # Calcular UniqPepts como el máximo por fila
  if (length(pept_cols) > 0) {
    pept_mat <- as.matrix(pq[, pept_cols, drop = FALSE])
    UniqPepts <- apply(pept_mat, 1, function(x) max(x, na.rm = TRUE))
    UniqPepts[!is.finite(UniqPepts)] <- 0L
  } else {
    UniqPepts <- rep(NA_integer_, nrow(pq))
  }

  # Construir resultado
  result <- data.frame(
    ProteinGroups = pq$PG.ProteinGroups,
    GeneNames = pq$PG.Genes,
    UniqPepts = as.integer(UniqPepts),
    stringsAsFactors = FALSE
  )

  # Agregar columnas de intensidad con nombres simplificados (Condicion_Replica)
  intensity_data <- pq[, quantity_cols, drop = FALSE]
  names(intensity_data) <- gsub("^PG\\.Quantity_", "", names(intensity_data))

  result <- cbind(result, intensity_data)
  result
}

# =============================================================================
# FUNCIONES REPLICADAS DE proteoDA (internas)
# =============================================================================

#' Convierte valores cero a NA
#'
#' Reemplaza valores 0 en una matriz o data frame con NA.
#' Útil para datos de DIA-NN donde 0 indica no detección.
#'
#' @param data Matriz o data frame numérico
#' @return Objeto del mismo tipo con 0 convertidos a NA
#' @keywords internal
.zero_to_missing <- function(data) {
  if (is.data.frame(data)) {
    numeric_cols <- sapply(data, is.numeric)
    data[numeric_cols] <- lapply(data[numeric_cols], function(x) {
      x[x == 0] <- NA
      x
    })
  } else if (is.matrix(data)) {
    data[data == 0] <- NA
  }
  data
}

#' Filtra proteínas por presencia en grupos
#'
#' Mantiene proteínas que tienen suficientes valores no-NA en al menos
#' un número mínimo de grupos experimentales.
#'
#' @param data Matrix de intensidades (proteínas × muestras)
#' @param metadata Data frame con información de muestras
#' @param min_reps Mínimo de réplicas con valores no-NA por grupo.
#'   Si NULL, se usa la mitad del tamaño del grupo más pequeño.
#' @param min_groups Mínimo de grupos que deben cumplir min_reps (default: 1)
#' @param grouping_column Nombre de la columna de agrupación en metadata (default: "Condition")
#' @return Lista con:
#'   - data: Matriz filtrada
#'   - keep: Vector lógico de filas conservadas
#'   - summary: Resumen del filtrado
#' @keywords internal
.filter_proteins_by_group <- function(
    data,
    metadata,
    min_reps = NULL,
    min_groups = 1,
    grouping_column = "Condition"
) {
  # Validaciones
  stopifnot(is.matrix(data) || is.data.frame(data))
  data <- as.matrix(data)

  if (!grouping_column %in% names(metadata)) {
    stop("La columna '", grouping_column, "' no existe en metadata")
  }

  # Alinear muestras
  if (!is.null(rownames(metadata))) {
    common_samples <- intersect(colnames(data), rownames(metadata))
    if (length(common_samples) == 0) {
      stop("No hay muestras en común entre data y metadata")
    }
    data <- data[, common_samples, drop = FALSE]
    metadata <- metadata[common_samples, , drop = FALSE]
  }

  groups <- as.factor(metadata[[grouping_column]])
  group_levels <- levels(groups)

  # Calcular min_reps automáticamente si no se especifica
  if (is.null(min_reps)) {
    group_sizes <- table(groups)
    min_reps <- max(1, floor(min(group_sizes) / 2))
  }

  # Calcular número de valores no-NA por grupo para cada proteína
  n_present_per_group <- sapply(group_levels, function(g) {
    cols <- which(groups == g)
    if (length(cols) == 0) return(rep(0, nrow(data)))
    rowSums(!is.na(data[, cols, drop = FALSE]))
  })

  if (!is.matrix(n_present_per_group)) {
    n_present_per_group <- matrix(n_present_per_group, ncol = 1)
  }

  # Contar grupos que cumplen el criterio
  groups_ok <- rowSums(n_present_per_group >= min_reps)
  keep <- groups_ok >= min_groups

  list(
    data = data[keep, , drop = FALSE],
    keep = keep,
    n_present_per_group = n_present_per_group,
    summary = list(
      n_total = nrow(data),
      n_keep = sum(keep),
      n_drop = sum(!keep),
      min_reps = min_reps,
      min_groups = min_groups
    )
  )
}

# =============================================================================
# FUNCIONES DE IMPUTACIÓN (internas)
# =============================================================================

#' Máscara MNAR por condición con evidencia en otras condiciones
#'
#' Identifica celdas candidatas a MNAR: NA en una condición pero con
#' presencia suficiente en otras condiciones.
#'
#' @param x Matriz de intensidades (proteínas × muestras)
#' @param condition Vector de condiciones alineado con columnas
#' @param prop_na_in_condition Proporción mínima de NA en la condición (default: 1.0 = 100%)
#' @param prop_present_in_other_condition Proporción mínima de no-NA en otra condición (default: 0.0)
#' @param min_present_in_other_condition Número mínimo absoluto de no-NA en otra condición (default: 1)
#' @param require_n_other_conditions Número de otras condiciones que deben cumplir (default: 1)
#' @param drop_empty_levels Eliminar niveles vacíos del factor (default: TRUE)
#' @return Matriz lógica de mismas dimensiones indicando celdas MNAR
#' @keywords internal
.mnar_mask_by_condition <- function(
    x, condition,
    prop_na_in_condition = 1.0,
    prop_present_in_other_condition = 0.0,
    min_present_in_other_condition = 1,
    require_n_other_conditions = 1,
    drop_empty_levels = TRUE
) {
  stopifnot(is.matrix(x) || is.data.frame(x))
  x <- as.matrix(x)
  stopifnot(ncol(x) == length(condition))
  stopifnot(prop_na_in_condition >= 0 && prop_na_in_condition <= 1)
  stopifnot(prop_present_in_other_condition >= 0 && prop_present_in_other_condition <= 1)
  stopifnot(min_present_in_other_condition >= 0)
  stopifnot(require_n_other_conditions >= 1)

  condition <- as.factor(condition)
  if (drop_empty_levels) condition <- droplevels(condition)

  mnar <- matrix(FALSE, nrow = nrow(x), ncol = ncol(x), dimnames = dimnames(x))

  for (g in levels(condition)) {
    jg <- which(condition == g)
    jn <- which(condition != g)
    if (length(jg) == 0) next

    # % NA dentro de g
    frac_na_g <- rowMeans(is.na(x[, jg, drop = FALSE]))
    cond_na_ok <- frac_na_g >= prop_na_in_condition

    # Evidencia en otras condiciones (evaluada por condición)
    present_ok_n <- integer(nrow(x))
    for (h in setdiff(levels(condition), g)) {
      jh <- which(condition == h)
      if (length(jh) == 0) next
      frac_present_h <- rowMeans(!is.na(x[, jh, drop = FALSE]))
      count_present_h <- rowSums(!is.na(x[, jh, drop = FALSE]))
      ok_h <- (frac_present_h >= prop_present_in_other_condition) &
        (count_present_h >= min_present_in_other_condition)
      present_ok_n <- present_ok_n + as.integer(ok_h)
    }
    cond_present_ok <- present_ok_n >= require_n_other_conditions

    rows_mnar_g <- cond_na_ok & cond_present_ok
    if (any(rows_mnar_g)) mnar[rows_mnar_g, jg] <- TRUE
  }
  mnar
}

#' Imputación MNAR con valor mínimo
#'
#' Reemplaza todos los NA con el valor mínimo de la matriz.
#' Implementación equivalente a MsCoreUtils::impute_min() sin dependencias.
#'
#' @param x Matriz numérica
#' @return Matriz con NA reemplazados por el mínimo global
#' @keywords internal
.impute_min <- function(x) {
  val <- min(x, na.rm = TRUE)
  x[is.na(x)] <- val
  x
}

#' Imputación mixta: MAR/MCAR con impSeqRob, MNAR con "min"
#'
#' Realiza imputación en dos etapas:
#' 1. MAR/MCAR: impSeqRob (rrcovNA) - robusto y recomendado
#' 2. MNAR: valor mínimo global (sin dependencias externas)
#'
#' @param x Matriz de intensidades log2
#' @param condition Vector de condiciones alineado con columnas
#' @param prop_na_in_condition Proporción de NA para clasificar MNAR (default: 1.0)
#' @param prop_present_in_other_condition Proporción presente en otras condiciones (default: 0.0)
#' @param min_present_in_other_condition Mínimo de valores presentes (default: 1)
#' @param require_n_other_conditions Condiciones requeridas con presencia (default: 1)
#' @param mar_method Método MAR: "impSeqRob" o "none" (default: "impSeqRob")
#' @param impSeqRob_args Lista de argumentos para impSeqRob (default: list(alpha = 0.9))
#' @return Lista con:
#'   - x_imputed: Matriz imputada
#'   - mnar_mask: Máscara de celdas MNAR
#'   - mar_mask: Máscara de celdas MAR
#'   - summary: Estadísticas de imputación
#' @keywords internal
.impute_mixed <- function(
    x, condition,
    prop_na_in_condition = 1.0,
    prop_present_in_other_condition = 0.0,
    min_present_in_other_condition = 1,
    require_n_other_conditions = 1,
    mar_method = c("impSeqRob", "none"),
    impSeqRob_args = list(alpha = 0.9, norm_impute = FALSE, check_data = FALSE, verbose = TRUE)
) {
  x <- as.matrix(x)
  stopifnot(ncol(x) == length(condition))
  mar_method <- match.arg(mar_method)

  # Construir máscaras
  mnar_mask <- .mnar_mask_by_condition(
    x, condition,
    prop_na_in_condition = prop_na_in_condition,
    prop_present_in_other_condition = prop_present_in_other_condition,
    min_present_in_other_condition = min_present_in_other_condition,
    require_n_other_conditions = require_n_other_conditions
  )
  mar_mask <- is.na(x) & !mnar_mask

  # ---- Etapa 1: Imputación MAR/MCAR con impSeqRob ----
  x_stage1 <- x
  if (mar_method == "impSeqRob") {
    if (!requireNamespace("rrcovNA", quietly = TRUE)) {
      stop("Para mar_method = 'impSeqRob' necesitas el paquete 'rrcovNA'.")
    }
    # Validar argumentos
    allowed <- c("alpha", "norm_impute", "check_data", "verbose")
    bad <- setdiff(names(impSeqRob_args), allowed)
    if (length(bad)) {
      warning("Argumentos no soportados para impSeqRob(): ",
              paste(bad, collapse = ", "), ". Se ignoran.")
    }
    args_final <- modifyList(
      list(alpha = 0.9, norm_impute = FALSE, check_data = FALSE, verbose = TRUE),
      impSeqRob_args[names(impSeqRob_args) %in% allowed]
    )

    imp1 <- do.call(rrcovNA::impSeqRob, c(list(x = x), args_final))
    x_imp1 <- if (is.list(imp1) && !is.null(imp1$x)) imp1$x else as.matrix(imp1)
    x_stage1[mar_mask] <- x_imp1[mar_mask]
  }
  # else: "none" -> deja MAR/MCAR como NA

  # ---- Etapa 2: Imputación MNAR con valor mínimo ----
  # Implementación simple equivalente a MsCoreUtils::impute_min()
  # Sin dependencia de MSnbase/Biobase
  x_min_all <- .impute_min(x_stage1)

  x_final <- x_stage1
  x_final[mnar_mask] <- x_min_all[mnar_mask]

  # Resumen de NA
  na0 <- mean(is.na(x))
  na1 <- mean(is.na(x_stage1))
  naF <- mean(is.na(x_final))

  list(
    x_imputed = x_final,
    mnar_mask = mnar_mask,
    mar_mask = mar_mask,
    summary = list(
      na_rate_initial = na0,
      na_rate_after_mar = na1,
      na_rate_final = naF,
      pct_na_marked_mnar = ifelse(any(is.na(x)), mean(mnar_mask[is.na(x)]), NA_real_)
    )
  )
}

#' Prefiltrado de proteínas por reglas MNAR
#'
#' Mantiene una proteína si:
#' (A) Tiene MNAR en ≥1 condición, O
#' (B) Presenta señal suficiente en ≥ require_n_other_conditions condiciones
#'
#' @param x Matriz de intensidades log2
#' @param condition Vector de condiciones
#' @param prop_na_in_condition Proporción NA para MNAR (default: 0.51)
#' @param prop_present_in_other_condition Proporción presente requerida (default: 0.5)
#' @param min_present_in_other_condition Mínimo de valores presentes (default: 1)
#' @param require_n_other_conditions Condiciones requeridas (default: 1)
#' @return Lista con:
#'   - keep: Vector lógico de filas a conservar
#'   - summary: Resumen del filtrado
#' @keywords internal
.prefilter_by_rules <- function(
    x, condition,
    prop_na_in_condition = 0.51,
    prop_present_in_other_condition = 0.5,
    min_present_in_other_condition = 1,
    require_n_other_conditions = 1
) {
  cond <- as.factor(condition)

  # 1) Máscara MNAR
  mnar_mask <- .mnar_mask_by_condition(
    x, cond,
    prop_na_in_condition = prop_na_in_condition,
    prop_present_in_other_condition = prop_present_in_other_condition,
    min_present_in_other_condition = min_present_in_other_condition,
    require_n_other_conditions = require_n_other_conditions
  )
  row_has_mnar <- rowSums(mnar_mask) > 0

  # 2) Presencia por condición
  levs <- levels(cond)
  present_ok_mat <- sapply(levs, function(g) {
    jg <- which(cond == g)
    if (length(jg) == 0) return(rep(FALSE, nrow(x)))
    frac_present_g <- rowMeans(!is.na(x[, jg, drop = FALSE]))
    count_present_g <- rowSums(!is.na(x[, jg, drop = FALSE]))
    (frac_present_g >= prop_present_in_other_condition) &
      (count_present_g >= min_present_in_other_condition)
  })
  if (!is.matrix(present_ok_mat)) present_ok_mat <- cbind(present_ok_mat)

  n_conditions_with_presence <- rowSums(present_ok_mat)
  row_has_presence <- n_conditions_with_presence >= require_n_other_conditions

  # 3) Vector final
  keep <- row_has_mnar | row_has_presence

  list(
    keep = keep,
    row_has_mnar = row_has_mnar,
    n_conditions_with_presence = n_conditions_with_presence,
    present_ok_by_condition = `colnames<-`(present_ok_mat, levs),
    mnar_mask = mnar_mask,
    summary = list(
      n_total = nrow(x),
      n_keep = sum(keep),
      n_drop = sum(!keep),
      n_with_MNAR = sum(row_has_mnar),
      n_with_presence_rule = sum(row_has_presence)
    )
  )
}

#' Renombra rownames usando columna de IDs
#'
#' @param x_df Matriz o data frame con rownames como ProteinGroups
#' @param rd Data frame de rowData con columna de IDs
#' @param id_col Nombre de la columna de IDs (default: "IDs")
#' @return Data frame con rownames renombrados
#' @keywords internal
.rename_rownames_from_rd <- function(x_df, rd, id_col = "IDs") {
  stopifnot(all(rownames(x_df) %in% rownames(rd)))
  ids <- rd[rownames(x_df), id_col, drop = TRUE] |> as.character()
  if (anyNA(ids)) {
    pg_missing <- rownames(x_df)[is.na(ids)]
    stop(sprintf("IDs faltantes para %d ProteinGroups (ejemplos: %s)",
                 length(pg_missing), paste(head(pg_missing, 10), collapse = ", ")))
  }
  if (anyDuplicated(ids)) {
    warning(sprintf("IDs duplicados detectados: %d. Se aplicará make.unique().",
                    sum(duplicated(ids))))
    ids <- make.unique(ids)
  }
  x_out <- as.data.frame(x_df, check.names = FALSE)
  rownames(x_out) <- ids
  x_out
}

# =============================================================================
# FUNCIONES REPLICADAS DE PRONE (internas)
# =============================================================================

#' Crea un SummarizedExperiment desde datos proteómicos
#'
#' @param data Data frame con proteínas y valores de intensidad
#' @param metadata Data frame con información de muestras
#' @param protein_column Nombre de columna con IDs de proteína
#' @param gene_column Nombre de columna con nombres de genes
#' @param condition_column Nombre de columna de condición en metadata
#' @param label_column Nombre de columna de etiqueta de muestra en metadata
#' @return SummarizedExperiment con assays: raw, log2
#' @keywords internal
.load_proteomics_data <- function(
    data,
    metadata,
    protein_column = "ProteinGroups",
    gene_column = "GeneNames",
    condition_column = "Condition",
    label_column = "Column"
) {
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Se requiere el paquete 'SummarizedExperiment'")
  }

  # Validar columnas requeridas
  if (!protein_column %in% names(data)) {
    stop("Columna '", protein_column, "' no encontrada en data")
  }
  if (!label_column %in% names(metadata)) {
    stop("Columna '", label_column, "' no encontrada en metadata")
  }
  if (!condition_column %in% names(metadata)) {
    stop("Columna '", condition_column, "' no encontrada en metadata")
  }

  # Identificar columnas de intensidad
  annotation_cols <- c(protein_column, gene_column, "UniqPepts")
  annotation_cols <- intersect(annotation_cols, names(data))
  intensity_cols <- setdiff(names(data), annotation_cols)

  # Extraer matrices
  annotation <- data[, annotation_cols, drop = FALSE]
  intensity <- as.matrix(data[, intensity_cols, drop = FALSE])
  storage.mode(intensity) <- "double"

  # IDs únicos
  protein_ids <- as.character(annotation[[protein_column]])
  protein_ids <- make.unique(protein_ids)
  rownames(intensity) <- protein_ids
  rownames(annotation) <- protein_ids

  # Alinear metadata con columnas de intensidad
  rownames(metadata) <- metadata[[label_column]]
  common_samples <- intersect(colnames(intensity), rownames(metadata))

  if (length(common_samples) == 0) {
    stop("No hay muestras en común entre data y metadata")
  }

  intensity <- intensity[, common_samples, drop = FALSE]
  metadata <- metadata[common_samples, , drop = FALSE]

  # Crear rowData
  row_data <- S4Vectors::DataFrame(annotation)
  names(row_data)[names(row_data) == protein_column] <- "Protein.IDs"
  if (gene_column %in% names(row_data)) {
    names(row_data)[names(row_data) == gene_column] <- "Gene.Names"
  }
  row_data$IDs <- rownames(row_data)

  # Crear colData
  col_data <- S4Vectors::DataFrame(metadata)

  # Crear assays: raw y log2
  raw_assay <- intensity
  log2_assay <- log2(intensity)

  # Crear SummarizedExperiment
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(raw = raw_assay, log2 = log2_assay),
    rowData = row_data,
    colData = col_data,
    metadata = list(
      condition = condition_column,
      label = label_column
    )
  )

  se
}

#' Resumen de valores NA en un SummarizedExperiment
#'
#' @param se SummarizedExperiment
#' @param assay_name Nombre del assay a evaluar (default: "log2")
#' @return Data frame con estadísticas de NA por muestra
#' @keywords internal
.get_NA_overview <- function(se, assay_name = "log2") {
  stopifnot(inherits(se, "SummarizedExperiment"))

  if (!assay_name %in% SummarizedExperiment::assayNames(se)) {
    stop("Assay '", assay_name, "' no encontrado")
  }

  x <- SummarizedExperiment::assay(se, assay_name)

  # Por muestra
  total_vals <- nrow(x)
  na_counts <- colSums(is.na(x))
  na_pct <- na_counts / total_vals * 100

  sample_stats <- data.frame(
    Sample = colnames(x),
    Total.Values = total_vals,
    NA.Values = na_counts,
    NA.Percentage = round(na_pct, 2),
    stringsAsFactors = FALSE
  )

  # Global
  total_cells <- length(x)
  total_na <- sum(is.na(x))

  attr(sample_stats, "global") <- list(
    Total.Cells = total_cells,
    Total.NA = total_na,
    NA.Percentage = round(total_na / total_cells * 100, 2)
  )

  sample_stats
}

#' Genera comparaciones para análisis diferencial
#'
#' @param se SummarizedExperiment
#' @param condition_column Columna de condición (si NULL, usa metadata del SE)
#' @param control Condición control. Si NULL, genera todas las comparaciones pareadas
#' @return Factor con comparaciones en formato "Tratamiento-Control"
#' @keywords internal
.specify_comparisons <- function(
    se,
    condition_column = NULL,
    control = NULL
) {
  stopifnot(inherits(se, "SummarizedExperiment"))

  # Obtener columna de condición
  if (is.null(condition_column)) {
    condition_column <- S4Vectors::metadata(se)$condition %||% "Condition"
  }

  cd <- as.data.frame(SummarizedExperiment::colData(se))
  if (!condition_column %in% names(cd)) {
    stop("Columna '", condition_column, "' no encontrada en colData")
  }

  conditions <- unique(as.character(cd[[condition_column]]))

  if (!is.null(control)) {
    if (!control %in% conditions) {
      stop("Control '", control, "' no está en las condiciones: ",
           paste(conditions, collapse = ", "))
    }
    # Comparaciones vs control
    others <- setdiff(conditions, control)
    comparisons <- paste0(others, "-", control)
  } else {
    # Todas las comparaciones pareadas
    comparisons <- character()
    for (i in seq_along(conditions)) {
      for (j in seq_along(conditions)) {
        if (i < j) {
          comparisons <- c(comparisons, paste0(conditions[j], "-", conditions[i]))
        }
      }
    }
  }

  factor(comparisons)
}

#' Ejecuta análisis limma
#'
#' @param data Matriz de intensidades log2
#' @param condition_vector Vector de condiciones alineado con columnas
#' @param comparisons Vector de comparaciones
#' @param covariate Covariable opcional para el modelo
#' @return Objeto fit de limma
#' @keywords internal
.perform_limma <- function(data, condition_vector, comparisons, covariate = NULL) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Se requiere el paquete 'limma'")
  }

  condition <- factor(condition_vector)

  # Crear matriz de diseño
  if (is.null(covariate)) {
    design <- model.matrix(~ 0 + condition)
    colnames(design) <- levels(condition)
  } else {
    design <- model.matrix(~ 0 + condition + covariate)
    colnames(design)[seq_along(levels(condition))] <- levels(condition)
  }

  # Crear matriz de contrastes
  contrast_strings <- as.character(comparisons)
  contrast_matrix <- limma::makeContrasts(
    contrasts = contrast_strings,
    levels = design
  )

  # Ajustar modelo
  fit <- limma::lmFit(data, design)
  fit <- limma::contrasts.fit(fit, contrast_matrix)
  fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)

  fit
}

#' Extrae resultados de limma fit
#'
#' @param fit Objeto fit de limma
#' @param comparisons Vector de comparaciones
#' @param logFC_up Umbral superior de logFC para "Up"
#' @param logFC_down Umbral inferior de logFC para "Down"
#' @param alpha Umbral de significancia
#' @param p_adj Usar p-valor ajustado (TRUE) o no ajustado (FALSE)
#' @return Data frame con resultados
#' @keywords internal
.extract_limma_results <- function(
    fit,
    comparisons,
    logFC_up = 1,
    logFC_down = -1,
    alpha = 0.05,
    p_adj = TRUE
) {
  results_list <- lapply(seq_along(comparisons), function(i) {
    comp <- as.character(comparisons[i])
    tt <- limma::topTable(fit, coef = i, number = Inf, sort.by = "none")

    # Columnas estándar
    df <- data.frame(
      Protein.IDs = rownames(tt),
      logFC = tt$logFC,
      P.Value = tt$P.Value,
      adj.P.Val = tt$adj.P.Val,
      stringsAsFactors = FALSE
    )

    # Clasificar cambios
    p_col <- if (p_adj) "adj.P.Val" else "P.Value"
    df$Change <- "No Change"
    df$Change[df$logFC >= logFC_up & df[[p_col]] < alpha] <- "Up"
    df$Change[df$logFC <= logFC_down & df[[p_col]] < alpha] <- "Down"
    df$Change <- factor(df$Change, levels = c("Up", "Down", "No Change"))

    df$Comparison <- comp
    df
  })

  do.call(rbind, results_list)
}

#' Ejecuta análisis de expresión diferencial
#'
#' Función principal que coordina el análisis DE con limma.
#'
#' @param se SummarizedExperiment con datos procesados
#' @param comparisons Comparaciones a realizar (resultado de .specify_comparisons)
#' @param assay_name Nombre del assay a usar. Si NULL, usa el último disponible
#' @param condition_column Columna de condición. Si NULL, usa metadata del SE
#' @param logFC Aplicar filtro por logFC (default: TRUE)
#' @param logFC_up Umbral superior de logFC (default: 1)
#' @param logFC_down Umbral inferior de logFC (default: -1)
#' @param p_adj Usar p-valor ajustado (default: TRUE)
#' @param alpha Umbral de significancia (default: 0.05)
#' @return Data frame con resultados DE
#' @keywords internal
.run_DE <- function(
    se,
    comparisons,
    assay_name = NULL,
    condition_column = NULL,
    logFC = TRUE,
    logFC_up = 1,
    logFC_down = -1,
    p_adj = TRUE,
    alpha = 0.05
) {
  stopifnot(inherits(se, "SummarizedExperiment"))

  # Determinar assay
  if (is.null(assay_name)) {
    assay_names <- SummarizedExperiment::assayNames(se)
    assay_name <- assay_names[length(assay_names)]
  }

  if (!assay_name %in% SummarizedExperiment::assayNames(se)) {
    stop("Assay '", assay_name, "' no encontrado")
  }

  # Obtener datos
  x <- SummarizedExperiment::assay(se, assay_name)
  cd <- as.data.frame(SummarizedExperiment::colData(se))
  rd <- as.data.frame(SummarizedExperiment::rowData(se))

  # Columna de condición
  if (is.null(condition_column)) {
    condition_column <- S4Vectors::metadata(se)$condition %||% "Condition"
  }

  condition_vec <- cd[[condition_column]]

  # Ejecutar limma
  fit <- .perform_limma(x, condition_vec, comparisons, covariate = NULL)

  # Extraer resultados
  if (!logFC) {
    logFC_up <- 0
    logFC_down <- 0
  }

  results <- .extract_limma_results(
    fit, comparisons,
    logFC_up = logFC_up,
    logFC_down = logFC_down,
    alpha = alpha,
    p_adj = p_adj
  )

  # Agregar información de genes
  if ("Gene.Names" %in% names(rd)) {
    gene_map <- rd[, c("Protein.IDs", "Gene.Names"), drop = FALSE]
    gene_map <- unique(gene_map)
    names(gene_map) <- c("Protein.IDs", "Gene.Names")
    results <- merge(results, gene_map, by = "Protein.IDs", all.x = TRUE, sort = FALSE)
  }

  # Agregar IDs
  if ("IDs" %in% names(rd)) {
    id_map <- rd[, c("Protein.IDs", "IDs"), drop = FALSE]
    id_map <- unique(id_map)
    results <- merge(results, id_map, by = "Protein.IDs", all.x = TRUE, sort = FALSE)
  }

  # Agregar columna Assay
  results$Assay <- assay_name

  # Reordenar columnas
  col_order <- c("Protein.IDs", "Gene.Names", "IDs", "logFC", "P.Value",
                 "adj.P.Val", "Change", "Comparison", "Assay")
  col_order <- intersect(col_order, names(results))
  results <- results[, col_order]

  results
}

# =============================================================================
# FUNCIÓN PRINCIPAL
# =============================================================================

#' Procesa datos proteómicos de Spectronaut
#'
#' Pipeline completo de procesamiento: filtrado, normalización, imputación
#' y análisis de expresión diferencial.
#'
#' @param preprocessing Lista de clase spectronaut_data (resultado de preprocess_spectronaut)
#' @param export_dir Directorio de salida para archivos exportados (default: "./results")
#' @param min_reps_filter Mínimo de réplicas para filtrado. Si NULL, se calcula automáticamente
#' @param min_groups_filter Mínimo de grupos para filtrado (default: 1)
#' @param normalization_method Método de normalización (default: "cyclicloess")
#' @param prop_na_mnar Proporción de NA para clasificar como MNAR (default: 0.51)
#' @param prop_present_mar Proporción presente para MAR (default: 0.5)
#' @param min_present_mar Mínimo de valores presentes para MAR (default: 1)
#' @param require_n_conditions Número de condiciones requeridas con presencia (default: 1)
#' @param mar_method Método de imputación MAR (default: "impSeqRob")
#' @param comparisons Comparaciones para DE. Si NULL, genera todas las pareadas
#' @param control Condición control para comparaciones. Si NULL, compara todas
#' @param logFC_threshold Umbral logFC para significancia (default: 0)
#' @param alpha Umbral p-valor ajustado (default: 0.05)
#' @param export_normalized Exportar matriz normalizada (default: TRUE)
#' @param export_imputed Exportar matriz imputada (default: TRUE)
#' @param verbose Mostrar mensajes de progreso (default: TRUE)
#'
#' @return Lista de clase proteomics_result con:
#'   \itemize{
#'     \item se_proc: SummarizedExperiment procesado con todos los assays
#'     \item DEPs_results: Data frame con resultados de expresión diferencial
#'     \item comparisons: Comparaciones realizadas
#'     \item parameters: Parámetros utilizados
#'   }
#'
#' @examples
#' \dontrun{
#' # 1. Preprocesar datos de Spectronaut
#' preprocessing <- preprocess_spectronaut(
#'   file_path = "data/Spectronaut_Report.tsv",
#'   condition_order = c("A", "B", "C", "D")
#' )
#'
#' # 2. Procesar (normalización, imputación, DE)
#' result <- process_proteomics(
#'   preprocessing = preprocessing,
#'   export_dir = "./results",
#'   alpha = 0.05
#' )
#'
#' # 3. Acceder a resultados
#' result$se_proc        # SummarizedExperiment procesado
#' result$DEPs_results   # Resultados diferenciales
#'
#' # 4. Visualización con funciones existentes
#' # volcano_highchart_list(result$DEPs_results, ...)
#' }
#'
#' @export
process_proteomics <- function(
    preprocessing,
    export_dir = "./results",
    min_reps_filter = NULL,
    min_groups_filter = 1,
    normalization_method = "cyclicloess",
    prop_na_mnar = 0.51,
    prop_present_mar = 0.5,
    min_present_mar = 1,
    require_n_conditions = 1,
    mar_method = "impSeqRob",
    comparisons = NULL,
    control = NULL,
    logFC_threshold = 0,
    alpha = 0.05,
    export_normalized = TRUE,
    export_imputed = TRUE,
    verbose = TRUE
) {
  # =========================================================================
  # VALIDACIONES
  # =========================================================================

  if (!inherits(preprocessing, "spectronaut_data")) {
    stop("El argumento 'preprocessing' debe ser resultado de preprocess_spectronaut()")
  }

  # Crear directorio de salida
  if (!dir.exists(export_dir)) {
    dir.create(export_dir, recursive = TRUE)
  }

  # Verificar paquetes requeridos
  required_packages <- c("SummarizedExperiment", "limma")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Se requiere el paquete '", pkg, "'. Instálalo con BiocManager::install('", pkg, "')")
    }
  }

  # =========================================================================
  # 1. PREPARAR DATOS DESDE PREPROCESSING
  # =========================================================================

  if (verbose) cat("=== PREPARANDO DATOS ===\n")

  metadata <- .prepare_metadata(preprocessing)
  protein_data <- .prepare_protein_data(preprocessing)

  if (verbose) {
    cat("- Metadatos:", nrow(metadata), "muestras\n")
    cat("- Proteínas:", nrow(protein_data), "proteínas iniciales\n")
  }

  # =========================================================================
  # 2. CONVERTIR CEROS A NA
  # =========================================================================

  if (verbose) cat("\n=== CONVIRTIENDO CEROS A NA ===\n")

  # Identificar columnas de intensidad
  annotation_cols <- c("ProteinGroups", "GeneNames", "UniqPepts")
  intensity_cols <- setdiff(names(protein_data), annotation_cols)

  # Convertir ceros
  intensity_mat <- as.matrix(protein_data[, intensity_cols])
  intensity_mat <- .zero_to_missing(intensity_mat)
  rownames(intensity_mat) <- protein_data$ProteinGroups

  n_zeros <- sum(protein_data[, intensity_cols] == 0, na.rm = TRUE)
  if (verbose) cat("- Ceros convertidos a NA:", n_zeros, "\n")

  # =========================================================================
  # 3. FILTRAR PROTEÍNAS POR GRUPO
  # =========================================================================

  if (verbose) cat("\n=== FILTRANDO PROTEÍNAS POR PRESENCIA ===\n")

  filtered <- .filter_proteins_by_group(
    data = intensity_mat,
    metadata = metadata,
    min_reps = min_reps_filter,
    min_groups = min_groups_filter,
    grouping_column = "Condition"
  )

  if (verbose) {
    cat("- Proteínas antes:", filtered$summary$n_total, "\n")
    cat("- Proteínas después:", filtered$summary$n_keep, "\n")
    cat("- Proteínas eliminadas:", filtered$summary$n_drop, "\n")
    cat("- Min réplicas:", filtered$summary$min_reps, "\n")
    cat("- Min grupos:", filtered$summary$min_groups, "\n")
  }

  # Filtrar protein_data
  protein_data_filtered <- protein_data[filtered$keep, , drop = FALSE]

  # =========================================================================
  # 4. CREAR SUMMARIZEDEXPERIMENT
  # =========================================================================

  if (verbose) cat("\n=== CREANDO SUMMARIZEDEXPERIMENT ===\n")

  se <- .load_proteomics_data(
    data = protein_data_filtered,
    metadata = metadata,
    protein_column = "ProteinGroups",
    gene_column = "GeneNames",
    condition_column = "Condition",
    label_column = "Column"
  )

  if (verbose) {
    na_overview <- .get_NA_overview(se, "log2")
    global_na <- attr(na_overview, "global")
    cat("- NA global:", global_na$NA.Percentage, "%\n")
  }

  # =========================================================================
  # 5. NORMALIZACIÓN
  # =========================================================================

  if (verbose) cat("\n=== NORMALIZANDO (", normalization_method, ") ===\n")

  x_log2 <- SummarizedExperiment::assay(se, "log2")

  x_norm <- limma::normalizeBetweenArrays(x_log2, method = normalization_method)
  rownames(x_norm) <- rownames(x_log2)

  # Exportar matriz normalizada
  if (export_normalized) {
    norm_file <- file.path(export_dir, paste0("matrix_log2_", normalization_method, ".tsv"))
    if (requireNamespace("readr", quietly = TRUE)) {
      readr::write_tsv(
        data.frame(ProteinGroups = rownames(x_norm), x_norm, check.names = FALSE),
        norm_file
      )
    } else {
      write.table(
        data.frame(ProteinGroups = rownames(x_norm), x_norm, check.names = FALSE),
        norm_file, sep = "\t", quote = FALSE, row.names = FALSE
      )
    }
    if (verbose) cat("- Exportado:", basename(norm_file), "\n")
  }

  # =========================================================================
  # 6. PREFILTRADO POR REGLAS MNAR
  # =========================================================================

  if (verbose) cat("\n=== PREFILTRADO POR REGLAS MNAR ===\n")

  cd <- as.data.frame(SummarizedExperiment::colData(se))
  condition_vec <- as.factor(cd$Condition)

  pf <- .prefilter_by_rules(
    x = x_norm,
    condition = condition_vec,
    prop_na_in_condition = prop_na_mnar,
    prop_present_in_other_condition = prop_present_mar,
    min_present_in_other_condition = min_present_mar,
    require_n_other_conditions = require_n_conditions
  )
  keep_rows <- pf$keep

  if (verbose) {
    cat("- Proteínas conservadas:", pf$summary$n_keep, "\n")
    cat("- Proteínas eliminadas:", pf$summary$n_drop, "\n")
  }

  x_norm_prefilt <- x_norm[keep_rows, , drop = FALSE]

  # =========================================================================
  # 7. IMPUTACIÓN MIXTA
  # =========================================================================

  if (verbose) cat("\n=== IMPUTACIÓN MIXTA (MAR + MNAR) ===\n")

  res_impute <- .impute_mixed(
    x = x_norm_prefilt,
    condition = condition_vec,
    prop_na_in_condition = prop_na_mnar,
    prop_present_in_other_condition = prop_present_mar,
    min_present_in_other_condition = min_present_mar,
    require_n_other_conditions = require_n_conditions,
    mar_method = mar_method
  )
  x_imputed <- res_impute$x_imputed

  if (verbose) {
    cat("- NA inicial:", round(res_impute$summary$na_rate_initial * 100, 2), "%\n")
    cat("- NA después MAR:", round(res_impute$summary$na_rate_after_mar * 100, 2), "%\n")
    cat("- NA final:", round(res_impute$summary$na_rate_final * 100, 2), "%\n")
  }

  # Exportar matriz imputada
  if (export_imputed) {
    imp_file <- file.path(export_dir, paste0("matrix_log2_", normalization_method, "_imputed.tsv"))
    if (requireNamespace("readr", quietly = TRUE)) {
      readr::write_tsv(
        data.frame(ProteinGroups = rownames(x_imputed), x_imputed, check.names = FALSE),
        imp_file
      )
    } else {
      write.table(
        data.frame(ProteinGroups = rownames(x_imputed), x_imputed, check.names = FALSE),
        imp_file, sep = "\t", quote = FALSE, row.names = FALSE
      )
    }
    if (verbose) cat("- Exportado:", basename(imp_file), "\n")
  }

  # =========================================================================
  # 8. ACTUALIZAR SUMMARIZEDEXPERIMENT
  # =========================================================================

  if (verbose) cat("\n=== ACTUALIZANDO SUMMARIZEDEXPERIMENT ===\n")

  # Filtrar SE a las filas procesadas
  rd <- as.data.frame(SummarizedExperiment::rowData(se))

  # Renombrar rownames de x_imputed usando IDs
  if ("IDs" %in% names(rd)) {
    rd_subset <- rd[rownames(x_imputed), , drop = FALSE]
    x_imputed_ids <- .rename_rownames_from_rd(x_imputed, rd_subset, id_col = "IDs")
  } else {
    x_imputed_ids <- x_imputed
  }

  # Alinear SE con matriz imputada
  common_ids <- intersect(rownames(se), rownames(x_imputed_ids))
  if (length(common_ids) == 0) {
    # Intentar con Protein.IDs
    common_ids <- intersect(rd$Protein.IDs, rownames(x_imputed))
    if (length(common_ids) > 0) {
      se_subset <- se[rd$Protein.IDs %in% common_ids, ]
      mat <- x_imputed[common_ids, colnames(se_subset), drop = FALSE]
    } else {
      stop("No hay IDs en común entre SE y matriz imputada")
    }
  } else {
    se_subset <- se[common_ids, ]
    mat <- as.matrix(x_imputed_ids[common_ids, colnames(se_subset), drop = FALSE])
  }

  storage.mode(mat) <- "double"

  # Agregar assay normalizado+imputado
  assay_name <- paste0(
    tools::toTitleCase(gsub("cyclic", "Cyc", normalization_method))
  )
  SummarizedExperiment::assay(se_subset, assay_name) <- mat

  se_proc <- se_subset

  if (verbose) {
    cat("- Assays disponibles:", paste(SummarizedExperiment::assayNames(se_proc), collapse = ", "), "\n")
    cat("- Proteínas finales:", nrow(se_proc), "\n")
  }

  # =========================================================================
  # 9. ANÁLISIS DIFERENCIAL
  # =========================================================================

  if (verbose) cat("\n=== ANÁLISIS DIFERENCIAL (limma) ===\n")

  # Generar comparaciones si no se especifican
  if (is.null(comparisons)) {
    comparisons <- .specify_comparisons(se_proc, condition_column = "Condition", control = control)
  }

  if (verbose) cat("- Comparaciones:", paste(comparisons, collapse = ", "), "\n")

  # Ejecutar análisis DE
  DEPs_results <- .run_DE(
    se = se_proc,
    comparisons = comparisons,
    assay_name = assay_name,
    condition_column = "Condition",
    logFC = TRUE,
    logFC_up = logFC_threshold,
    logFC_down = -logFC_threshold,
    p_adj = TRUE,
    alpha = alpha
  )

  if (verbose) {
    n_sig <- sum(DEPs_results$Change != "No Change")
    cat("- Proteínas diferenciales (total):", n_sig, "\n")

    for (comp in unique(DEPs_results$Comparison)) {
      subset <- DEPs_results[DEPs_results$Comparison == comp, ]
      n_up <- sum(subset$Change == "Up")
      n_down <- sum(subset$Change == "Down")
      cat("  ", comp, ": Up=", n_up, ", Down=", n_down, "\n", sep = "")
    }
  }

  # =========================================================================
  # 10. EXPORTAR RESULTADOS DE
  # =========================================================================

  if (verbose) cat("\n=== EXPORTANDO RESULTADOS ===\n")

  de_file <- file.path(export_dir, "DEPs_results.tsv")
  if (requireNamespace("readr", quietly = TRUE)) {
    readr::write_tsv(DEPs_results, de_file)
  } else {
    write.table(DEPs_results, de_file, sep = "\t", quote = FALSE, row.names = FALSE)
  }

  if (verbose) cat("- Exportado:", basename(de_file), "\n")

  # =========================================================================
  # 11. RETORNAR RESULTADO
  # =========================================================================

  result <- list(
    se_proc = se_proc,
    DEPs_results = DEPs_results,
    comparisons = comparisons,
    parameters = list(
      min_reps_filter = filtered$summary$min_reps,
      min_groups_filter = min_groups_filter,
      normalization_method = normalization_method,
      prop_na_mnar = prop_na_mnar,
      prop_present_mar = prop_present_mar,
      mar_method = mar_method,
      logFC_threshold = logFC_threshold,
      alpha = alpha,
      export_dir = export_dir
    )
  )

  class(result) <- c("proteomics_result", "list")

  if (verbose) cat("\n=== PROCESAMIENTO COMPLETADO ===\n")

  result
}

# =============================================================================
# MÉTODO PRINT
# =============================================================================

#' Método print para proteomics_result
#'
#' @param x Objeto proteomics_result
#' @param ... Argumentos adicionales (ignorados)
#' @export
print.proteomics_result <- function(x, ...) {
  cat("=== Resultado de Procesamiento Proteómico ===\n\n")

  # Resumen del SE
  se <- x$se_proc
  cat("SummarizedExperiment:\n")
  cat("  - Proteínas:", nrow(se), "\n")
  cat("  - Muestras:", ncol(se), "\n")
  cat("  - Assays:", paste(SummarizedExperiment::assayNames(se), collapse = ", "), "\n")

  # Resumen de condiciones
  cd <- as.data.frame(SummarizedExperiment::colData(se))
  if ("Condition" %in% names(cd)) {
    cat("  - Condiciones:", paste(unique(cd$Condition), collapse = ", "), "\n")
  }

  cat("\nResultados Diferenciales:\n")
  cat("  - Total filas:", nrow(x$DEPs_results), "\n")
  cat("  - Comparaciones:", paste(unique(x$DEPs_results$Comparison), collapse = ", "), "\n")

  # Resumen por comparación
  for (comp in unique(x$DEPs_results$Comparison)) {
    subset <- x$DEPs_results[x$DEPs_results$Comparison == comp, ]
    n_up <- sum(subset$Change == "Up", na.rm = TRUE)
    n_down <- sum(subset$Change == "Down", na.rm = TRUE)
    cat("    ", comp, ": Up=", n_up, ", Down=", n_down, "\n", sep = "")
  }

  cat("\nParámetros:\n")
  cat("  - Normalización:", x$parameters$normalization_method, "\n")
  cat("  - Alpha:", x$parameters$alpha, "\n")
  cat("  - logFC threshold:", x$parameters$logFC_threshold, "\n")
  cat("  - Directorio salida:", x$parameters$export_dir, "\n")

  invisible(x)
}

# =============================================================================
# EJEMPLOS DE USO (comentados)
# =============================================================================

# # Flujo completo de procesamiento
# #
# # 1. Preprocesar datos de Spectronaut
# # preprocessing <- preprocess_spectronaut(
# #   file_path = "data/Spectronaut_Report.tsv",
# #   condition_order = c("A", "B", "C", "D")
# # )
# #
# # 2. Procesar (normalización, imputación, DE)
# # result <- process_proteomics(
# #   preprocessing = preprocessing,
# #   export_dir = "./results",
# #   alpha = 0.05,
# #   logFC_threshold = 0.5
# # )
# #
# # 3. Acceder a resultados
# # result$se_proc        # SummarizedExperiment procesado
# # result$DEPs_results   # Resultados diferenciales
# #
# # 4. Usar funciones de visualización existentes
# # source("./R/Volcano_Plot_Highcharts.R")
# # volcano_highchart_list(
# #   result$DEPs_results,
# #   logFC_col = "logFC",
# #   pval_col = "adj.P.Val",
# #   gene_col = "Gene.Names",
# #   comparison_col = "Comparison"
# # )
