# =============================================================================
# Preprocesamiento de Datos de Spectronaut
# =============================================================================
#
# Convierte reportes de Spectronaut a formato estructurado para análisis
# proteómico downstream.
#
# Copyright 2025 Sergio Ciordia
# Licensed under MIT
# =============================================================================

# --- Dependencias ---

# =============================================================================
# Funciones Auxiliares Internas
# =============================================================================

#' Obtiene función agregadora por nombre
#' @noRd
.get_aggregator <- function(name) {
  nm <- tolower(name %||% "")
  switch(nm,
    "max"    = function(x) max(x, na.rm = TRUE),
    "mean"   = function(x) mean(x, na.rm = TRUE),
    "median" = function(x) stats::median(x, na.rm = TRUE),
    "min"    = function(x) min(x, na.rm = TRUE),
    stop("Agregador no soportado: '", name, "'. Usa: max, mean, median, min.")
  )
}

#' Retorna el primer valor no-NA de un vector
#' @noRd
.first_non_na <- function(x) {

  y <- na.omit(x)
  if (length(y) == 0) NA else y[1]
}

#' Limpia campos con valores separados por punto y coma
#' @description Convierte campos tipo "16,8%;16.8%" a numérico y agrega por grupos
#' @noRd
.clean_semicolon_numeric <- function(df, value_col, group_cols, out_col = value_col,
                                     agg_fun = .get_aggregator("max")) {
  df %>%
    select(all_of(c(group_cols, value_col))) %>%
    mutate(!!value_col := as.character(.data[[value_col]])) %>%
    separate_rows(all_of(value_col), sep = ";") %>%
    mutate(
      .val_num = suppressWarnings(
        as.numeric(gsub("[^0-9.]+", "", gsub(",", ".", .data[[value_col]])))
      )
    ) %>%
    group_by(across(all_of(group_cols))) %>%
    summarise(
      !!out_col := if (all(is.na(.val_num))) NA_real_ else agg_fun(.val_num),
      .groups = "drop"
    )
}

#' Normaliza nombres de columnas de Spectronaut
#' @description Estandariza variantes de nombres de columnas a formato consistente
#' @noRd
.normalize_column_names <- function(df) {
  nm <- names(df)


  # PG.Cscore.RunWise
  cscore_variants <- c("PG.Cscore..Run.Wise.", "PG.Cscore (Run-Wise)", "PG.Cscore_RunWise")
  for (variant in cscore_variants) {
    if (variant %in% nm) {
      names(df)[nm == variant] <- "PG.Cscore.RunWise"
      nm <- names(df)
      break
    }
  }


  # Columnas globales/experiment-wide
  global_mappings <- list(
    "PG.NrOfPrecursorsIdentified..Experiment.wide." = "PG.NrOfPrecursorsIdentified.Global",
    "PG.NrOfStrippedSequencesIdentified..Experiment.wide." = "PG.NrOfStrippedSequencesIdentified.Global",
    "PG.Coverage..Global." = "PG.Coverage.Global"
  )

  for (old_name in names(global_mappings)) {
    if (old_name %in% nm) {
      names(df)[nm == old_name] <- global_mappings[[old_name]]
      nm <- names(df)
    }
  }

  df
}

#' Construye columna Coding y ordena por condición/replicado
#' @noRd
.make_coding <- function(df, cond_order) {
  df %>%
    mutate(
      R.Condition = factor(R.Condition, levels = cond_order, ordered = TRUE),
      R.Replicate = as.integer(R.Replicate),
      Coding = paste(R.Condition, R.Replicate, sep = "_")
    ) %>%
    arrange(R.Condition, R.Replicate)
}

#' Valida columnas requeridas en dataframe
#' @noRd
.validate_spectronaut_columns <- function(df, required_cols) {
  missing <- setdiff(required_cols, names(df))
  if (length(missing) > 0) {
    stop(
      "Columnas requeridas faltantes en el archivo:\n",
      "  - ", paste(missing, collapse = "\n  - "), "\n",
      "Verifica que el archivo sea un reporte válido de Spectronaut."
    )
  }
  invisible(TRUE)
}

#' Extrae información base de proteínas (estática, 1 fila por proteína)
#' @noRd
.extract_base_info <- function(df, mw_clean) {
  static_cols <- c("PG.ProteinGroups", "PG.ProteinDescriptions", "PG.Genes", "PG.MolecularWeight")

  df %>%
    select(all_of(static_cols)) %>%
    group_by(PG.ProteinGroups) %>%
    summarise(
      PG.ProteinDescriptions = .first_non_na(PG.ProteinDescriptions),
      PG.Genes = .first_non_na(PG.Genes),
      PG.MolecularWeight = .first_non_na(PG.MolecularWeight),
      .groups = "drop"
    ) %>%
    select(-PG.MolecularWeight) %>%
    left_join(mw_clean, by = "PG.ProteinGroups")
}

# =============================================================================
# Función Principal
# =============================================================================

#' Preprocesa reportes de Spectronaut
#'
#' @description
#' Convierte un reporte de Spectronaut (formato TSV largo) a tres tablas
#' estructuradas: metadata de runs, identificación de proteínas y cuantificación.
#'
#' @param file_path Ruta al archivo TSV de Spectronaut.
#' @param condition_order Vector de caracteres con el orden de las condiciones
#'   experimentales (ej: `c("Control", "Tratado")`).
#' @param export_dir Directorio para exportar archivos TSV. Si es `NULL` (default),
#'   no se exportan archivos.
#' @param agg_coverage_run Método de agregación para PG.Coverage por muestra.
#'   Opciones: "max" (default), "mean", "median", "min".
#' @param agg_coverage_global Método de agregación para PG.Coverage.Global.
#'   Opciones: "max" (default), "mean", "median", "min".
#' @param agg_mw Método de agregación para PG.MolecularWeight.
#'   Opciones: "max" (default), "mean", "median", "min".
#' @param agg_cscore_runwise Método de agregación para PG.Cscore.RunWise.
#'   Opciones: "mean" (default), "max", "median", "min".
#' @param timestamp_suffix Lógico. Si `TRUE` (default), añade timestamp a nombres
#'   de archivos exportados.
#' @param verbose Lógico. Si `TRUE` (default), muestra mensajes de progreso.
#'
#' @return Lista con clase `spectronaut_data` conteniendo:
#'   \describe{
#'     \item{metadata}{Data frame con información de runs (1 fila por muestra)}
#'     \item{protein_id}{Data frame con métricas de identificación por proteína}
#'     \item{protein_quant}{Data frame con métricas de cuantificación por proteína}
#'   }
#'
#' @examples
#' \dontrun{
#' # Uso básico
#' result <- preprocess_spectronaut(
#'   file_path = "data/Spectronaut_Report.tsv",
#'   condition_order = c("Control", "Treatment")
#' )
#'
#' # Acceder a componentes
#' head(result$metadata)
#' head(result$protein_quant)
#'
#' # Con exportación
#' result <- preprocess_spectronaut(
#'   file_path = "data/Spectronaut_Report.tsv",
#'   condition_order = c("A", "B", "C", "D"),
#'   export_dir = "./results",
#'   agg_coverage_run = "mean",
#'   verbose = TRUE
#' )
#' }
#'
#' @export
preprocess_spectronaut <- function(
    file_path,
    condition_order,
    export_dir = NULL,
    agg_coverage_run = c("max", "mean", "median", "min"),
    agg_coverage_global = c("max", "mean", "median", "min"),
    agg_mw = c("max", "mean", "median", "min"),
    agg_cscore_runwise = c("mean", "max", "median", "min"),
    timestamp_suffix = TRUE,
    verbose = TRUE
) {

  # --- Validación de argumentos ---
  agg_coverage_run <- match.arg(agg_coverage_run)
  agg_coverage_global <- match.arg(agg_coverage_global)
  agg_mw <- match.arg(agg_mw)
  agg_cscore_runwise <- match.arg(agg_cscore_runwise)

  if (!file.exists(file_path)) {
    stop("Archivo no encontrado: ", file_path)
  }

  if (length(condition_order) == 0 || !is.character(condition_order)) {
    stop("condition_order debe ser un vector de caracteres no vacío.")
  }

  # --- Lectura del archivo ---
  if (verbose) message("Leyendo archivo: ", basename(file_path))

  df <- read.delim(
    file_path,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = TRUE
  )

  # Normalizar nombres de columnas

  df <- .normalize_column_names(df)

  # Validar columnas requeridas
  required_cols <- c(
    "R.FileName", "R.Condition", "R.Replicate",
    "PG.ProteinGroups", "PG.Quantity"
  )
  .validate_spectronaut_columns(df, required_cols)

  # Crear columna Coding
  df <- .make_coding(df, condition_order)

  # Obtener niveles de Coding ordenados (calculado una sola vez)
  coding_levels <- df %>%
    distinct(R.Condition, R.Replicate, Coding) %>%
    arrange(R.Condition, R.Replicate) %>%
    pull(Coding)

  # --- Crear run_summary (metadata) ---
  if (verbose) message("Generando metadata de runs...")

  summary_cols <- c(
    "R.FileName", "R.Condition", "R.Replicate", "Coding",
    "R.PrecursorsIdentified", "R.StrippedSequencesIdentified",
    "R.ProteinGroupsIdentified"
  )

  run_summary <- df %>%
    select(any_of(summary_cols)) %>%
    distinct() %>%
    arrange(R.FileName) %>%
    group_by(R.FileName) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    arrange(R.Condition, R.Replicate) %>%
    as.data.frame()

  rownames(run_summary) <- run_summary$Coding

  # --- Preparar agregadores ---
  cov_fun_run <- .get_aggregator(agg_coverage_run)
  cov_fun_global <- .get_aggregator(agg_coverage_global)
  mw_fun <- .get_aggregator(agg_mw)
  cscore_agg <- .get_aggregator(agg_cscore_runwise)

  # --- MW limpio (compartido entre protein_ID y protein_QUANT) ---
  mw_clean <- .clean_semicolon_numeric(
    df, "PG.MolecularWeight",
    group_cols = "PG.ProteinGroups",
    out_col = "PG.MolecularWeight",
    agg_fun = mw_fun
  )

  # ==========================================================================
  # protein_ID
  # ==========================================================================
  if (verbose) message("Procesando protein_ID...")

  # Coverage run-wise limpio

  coverage_clean <- .clean_semicolon_numeric(
    df, "PG.Coverage",
    group_cols = c("PG.ProteinGroups", "Coding"),
    out_col = "PG.Coverage",
    agg_fun = cov_fun_run
  )

  # Integrar coverage limpio y tipar numéricos
  df2 <- df %>%
    select(-PG.Coverage) %>%
    left_join(coverage_clean, by = c("PG.ProteinGroups", "Coding")) %>%
    mutate(
      PG.NrOfPrecursorsIdentified = suppressWarnings(as.numeric(PG.NrOfPrecursorsIdentified)),
      PG.NrOfStrippedSequencesIdentified = suppressWarnings(as.numeric(PG.NrOfStrippedSequencesIdentified)),
      PG.Cscore.RunWise = suppressWarnings(as.numeric(PG.Cscore.RunWise))
    )

  # Información base
  base_info <- .extract_base_info(df2, mw_clean)

  # Métricas run-wise a formato ancho
  runwise_cols <- c(
    "PG.NrOfPrecursorsIdentified", "PG.NrOfStrippedSequencesIdentified",
    "PG.Coverage", "PG.Cscore.RunWise"
  )

  runwise_wide <- df2 %>%
    select(PG.ProteinGroups, Coding, all_of(runwise_cols)) %>%
    pivot_longer(cols = all_of(runwise_cols), names_to = "metric", values_to = "value") %>%
    group_by(PG.ProteinGroups, Coding, metric) %>%
    summarise(
      value = if (all(is.na(value))) NA_real_ else
        if (unique(metric) == "PG.Cscore.RunWise") cscore_agg(value) else
          if (unique(metric) == "PG.Coverage") cov_fun_run(value) else
            max(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      Coding = factor(Coding, levels = coding_levels, ordered = TRUE),
      metric_coding = paste0(metric, "_", Coding)
    ) %>%
    select(PG.ProteinGroups, metric_coding, value) %>%
    pivot_wider(names_from = metric_coding, values_from = value, names_repair = "unique")

  # Ensamblar protein_ID
  protein_ID <- base_info %>%
    left_join(runwise_wide, by = "PG.ProteinGroups")

  # Ordenar columnas
  static_cols <- c("PG.ProteinGroups", "PG.ProteinDescriptions", "PG.Genes", "PG.MolecularWeight")
  metric_order <- c(
    "PG.NrOfPrecursorsIdentified", "PG.NrOfStrippedSequencesIdentified",
    "PG.Coverage", "PG.Cscore.RunWise"
  )
  runwise_order <- unlist(lapply(metric_order, function(m) paste0(m, "_", coding_levels)))
  final_cols <- c(static_cols, intersect(runwise_order, names(protein_ID)))

  protein_ID <- protein_ID %>%
    select(all_of(final_cols)) %>%
    arrange(PG.ProteinGroups)

  # Validar unicidad
  if (n_distinct(protein_ID$PG.ProteinGroups) != nrow(protein_ID)) {
    stop(
      "Error de integridad: protein_ID contiene filas duplicadas por PG.ProteinGroups. ",
      "Revisa los datos de entrada."
    )
  }

  # ==========================================================================
  # protein_QUANT
  # ==========================================================================
  if (verbose) message("Procesando protein_QUANT...")

  # Coverage.Global limpio
  coverage_global_clean <- .clean_semicolon_numeric(
    df, "PG.Coverage.Global",
    group_cols = "PG.ProteinGroups",
    out_col = "PG.Coverage.Global",
    agg_fun = cov_fun_global
  )

  # Preparar df para QUANT
  df4 <- df %>%
    select(-any_of("PG.Coverage.Global")) %>%
    left_join(coverage_global_clean, by = "PG.ProteinGroups") %>%
    mutate(
      PG.NrOfPrecursorsIdentified.Global = suppressWarnings(as.numeric(PG.NrOfPrecursorsIdentified.Global)),
      PG.NrOfStrippedSequencesIdentified.Global = suppressWarnings(as.numeric(PG.NrOfStrippedSequencesIdentified.Global)),
      PG.NrOfPrecursorsUsedForQuantification = suppressWarnings(as.numeric(PG.NrOfPrecursorsUsedForQuantification)),
      PG.NrOfStrippedSequencesUsedForQuantification = suppressWarnings(as.numeric(PG.NrOfStrippedSequencesUsedForQuantification)),
      PG.Quantity = suppressWarnings(as.numeric(PG.Quantity)),
      PG.Cscore = suppressWarnings(as.numeric(PG.Cscore))
    )

  # Información base (reutiliza mw_clean)
  base_info2 <- .extract_base_info(df4, mw_clean)

  # Métricas globales
  global_metrics <- df4 %>%
    group_by(PG.ProteinGroups) %>%
    summarise(
      PG.NrOfPrecursorsIdentified.Global = .first_non_na(PG.NrOfPrecursorsIdentified.Global),
      PG.NrOfStrippedSequencesIdentified.Global = .first_non_na(PG.NrOfStrippedSequencesIdentified.Global),
      PG.Coverage.Global = .first_non_na(PG.Coverage.Global),
      PG.Cscore = .first_non_na(PG.Cscore),
      .groups = "drop"
    )

  # Métricas por muestra a formato ancho
  pivot_metrics <- c(
    "PG.NrOfPrecursorsUsedForQuantification",
    "PG.NrOfStrippedSequencesUsedForQuantification",
    "PG.Quantity"
  )

  runwise_wide2 <- df4 %>%
    select(PG.ProteinGroups, Coding, all_of(pivot_metrics)) %>%
    distinct() %>%
    pivot_longer(cols = all_of(pivot_metrics), names_to = "metric", values_to = "value") %>%
    group_by(PG.ProteinGroups, Coding, metric) %>%
    summarise(
      value = if (all(is.na(value))) NA_real_ else .first_non_na(value),
      .groups = "drop"
    ) %>%
    mutate(
      Coding = factor(Coding, levels = coding_levels, ordered = TRUE),
      metric_coding = paste0(metric, "_", Coding)
    ) %>%
    select(PG.ProteinGroups, metric_coding, value) %>%
    pivot_wider(names_from = metric_coding, values_from = value, names_repair = "unique")

  # Ensamblar protein_QUANT
  protein_QUANT <- base_info2 %>%
    left_join(global_metrics, by = "PG.ProteinGroups") %>%
    left_join(runwise_wide2, by = "PG.ProteinGroups")

  # Ordenar columnas
  static_cols2 <- c("PG.ProteinGroups", "PG.ProteinDescriptions", "PG.Genes", "PG.MolecularWeight")
  global_cols <- c(
    "PG.NrOfPrecursorsIdentified.Global",
    "PG.NrOfStrippedSequencesIdentified.Global",
    "PG.Coverage.Global",
    "PG.Cscore"
  )
  runwise_order2 <- unlist(lapply(pivot_metrics, function(m) paste0(m, "_", coding_levels)))
  final_cols2 <- c(static_cols2, global_cols, intersect(runwise_order2, names(protein_QUANT)))

  protein_QUANT <- protein_QUANT %>%
    select(all_of(final_cols2)) %>%
    arrange(PG.ProteinGroups)

  # Validar unicidad
  if (n_distinct(protein_QUANT$PG.ProteinGroups) != nrow(protein_QUANT)) {
    stop(
      "Error de integridad: protein_QUANT contiene filas duplicadas por PG.ProteinGroups. ",
      "Revisa los datos de entrada."
    )
  }

  # ==========================================================================
  # Exportación opcional
  # ==========================================================================
  if (!is.null(export_dir)) {
    if (verbose) message("Exportando archivos a: ", export_dir)

    if (!dir.exists(export_dir)) {
      dir.create(export_dir, recursive = TRUE)
    }

    suffix <- if (timestamp_suffix) {
      paste0("_", format(Sys.time(), "%Y%m%d_%H%M%S"))
    } else {
      ""
    }

    readr::write_tsv(
      run_summary,
      file.path(export_dir, paste0("Metadata", suffix, ".tsv")),
      na = ""
    )
    readr::write_tsv(
      protein_ID,
      file.path(export_dir, paste0("Protein_ID", suffix, ".tsv")),
      na = ""
    )
    readr::write_tsv(
      protein_QUANT,
      file.path(export_dir, paste0("Protein_QUANT", suffix, ".tsv")),
      na = ""
    )

    if (verbose) message("Archivos exportados exitosamente.")
  }

  # ==========================================================================
  # Construir resultado
  # ==========================================================================
  if (verbose) {
    message(
      "Procesamiento completado:\n",
      "  - Runs: ", nrow(run_summary), "\n",
      "  - Proteínas (ID): ", nrow(protein_ID), "\n",
      "  - Proteínas (QUANT): ", nrow(protein_QUANT)
    )
  }

  result <- list(
    metadata = as.data.frame(run_summary),
    protein_id = as.data.frame(protein_ID),
    protein_quant = as.data.frame(protein_QUANT)
  )

  class(result) <- c("spectronaut_data", "proteomics_data", "list")
  return(result)
}

# =============================================================================
# Métodos para clase spectronaut_data
# =============================================================================

#' @export
print.spectronaut_data <- function(x, ...) {
  cat("Datos de Spectronaut preprocesados\n")
  cat("----------------------------------\n")
  cat("Runs (metadata):", nrow(x$metadata), "\n")
  cat("Proteínas (ID):", nrow(x$protein_id), "\n")
  cat("Proteínas (QUANT):", nrow(x$protein_quant), "\n")
  cat("\nCondiciones:", paste(unique(x$metadata$R.Condition), collapse = ", "), "\n")
  invisible(x)
}
