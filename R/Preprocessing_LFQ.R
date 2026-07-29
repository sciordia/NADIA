# =============================================================================
# Preprocesamiento de Datos LFQ de Proteome Discoverer
# =============================================================================
#
# Convierte exports de proteínas LFQ de Proteome Discoverer a la misma estructura
# que `preprocess_spectronaut()` / `preprocess_tmt()` para alimentar el pipeline
# downstream sin cambios (Processing.R y módulos asociados).
#
# A diferencia de TMT (métricas globales por proteína), un experimento LFQ trae
# métricas POR MUESTRA en el mismo archivo, como el report de Spectronaut:
#   - `Abundance: <muestra>`                       (intensidad)
#   - `# PSMs (by Search Engine): <muestra>`       (PSMs por muestra)
#   - `# Peptides (by Search Engine): <muestra>`   (péptidos por muestra)
#   - `Score Mascot: <muestra>`                    (score por muestra)
#
# El mapeo muestra->condición se toma del archivo de anotación (`_Annot`), con
# columnas `Column`, `Condition` (obligatorias) y `Experiment` (opcional). Las
# intensidades se mapean POR NOMBRE de muestra (no por posición): el orden de las
# columnas `Abundance:` puede diferir del de las columnas de métricas.
#
# Copyright 2025 Sergio Ciordia
# Licensed under MIT
# =============================================================================

# --- Dependencias ---
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(tibble)
  library(readr)
})

# =============================================================================
# Funciones Auxiliares Internas
# =============================================================================

#' Extrae Gene Name de la columna Description (formato UniProt embebido)
#' @description Busca el patrón `GN=<gene>` habitual en los exports de PD.
#' @noRd
.parse_gene_from_description <- function(x) {
  m <- stringr::str_match(x, "GN=([^ ]+)")
  m[, 2]
}

#' Extrae la especie (OS=...) de la columna Description
#' @description Toma el texto entre `OS=` y el siguiente token `XX=` (p.ej. OX=).
#' @noRd
.parse_species_from_description <- function(x) {
  m <- stringr::str_match(x, "OS=(.+?)\\s+[A-Za-z]+=")
  m[, 2]
}

#' Lookup tolerante a variantes habituales de nombres de columna de PD
#' @noRd
.pd_resolve_col <- function(df, ...) {
  candidates <- unlist(list(...), use.names = FALSE)
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) NA_character_ else hit[1]
}

#' Extrae un vector numérico de `df[[col]]` o NA si la columna no existe
#' @noRd
.pd_numeric_or_na <- function(df, col) {
  if (is.na(col) || !col %in% names(df)) return(rep(NA_real_, nrow(df)))
  suppressWarnings(as.numeric(df[[col]]))
}

#' Extrae un vector character de `df[[col]]` o NA si la columna no existe
#' @noRd
.pd_char_or_na <- function(df, col) {
  if (is.na(col) || !col %in% names(df)) return(rep(NA_character_, nrow(df)))
  as.character(df[[col]])
}

#' Valida columnas mínimas requeridas para un export de Proteome Discoverer
#' @noRd
.validate_pd_columns <- function(df) {
  required <- c("Accession", "Description")
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop(
      "Columnas requeridas faltantes en el archivo de datos:\n  - ",
      paste(missing, collapse = "\n  - "),
      "\nVerifica que el archivo sea un export de proteínas de Proteome Discoverer."
    )
  }
  invisible(TRUE)
}

#' Valida el archivo de anotación (diseño experimental)
#' @noRd
.validate_lfq_annot <- function(annot) {
  required <- c("Column", "Condition")
  missing <- setdiff(required, names(annot))
  if (length(missing) > 0) {
    stop(
      "Columnas requeridas faltantes en el archivo de anotación (_Annot):\n  - ",
      paste(missing, collapse = "\n  - "),
      "\nEl _Annot debe tener al menos 'Column' y 'Condition'."
    )
  }
  if (anyDuplicated(annot$Column) > 0) {
    stop("El archivo de anotación tiene valores duplicados en 'Column'.")
  }
  invisible(TRUE)
}

#' Mapea columnas de una familia (por prefijo) a nombres de muestra
#' @description Hace `grep(prefix_regex, ...)`, quita el prefijo para obtener el
#'   nombre de muestra y devuelve un vector nombrado (nombre de muestra -> nombre
#'   de columna) restringido a las muestras del diseño. El mapeo es POR NOMBRE,
#'   por lo que es robusto a diferencias de orden entre familias de columnas.
#' @return Vector character nombrado por muestra (subset de `samples`).
#' @noRd
.lfq_sample_cols <- function(df, prefix_regex, samples) {
  cols <- grep(prefix_regex, names(df), value = TRUE)
  if (length(cols) == 0) return(setNames(character(0), character(0)))
  smp <- sub(prefix_regex, "", cols)
  keep <- smp %in% samples
  setNames(cols[keep], smp[keep])
}

# =============================================================================
# Función Principal
# =============================================================================

#' Preprocesa exports LFQ de Proteome Discoverer
#'
#' @description
#' Convierte un export de proteínas LFQ de Proteome Discoverer (formato TSV ancho)
#' a la misma estructura de salida que `preprocess_spectronaut()` /
#' `preprocess_tmt()`: tres data.frames (`metadata`, `protein_id`, `protein_quant`)
#' listos para el pipeline downstream (`process_proteomics()`).
#'
#' El diseño experimental (muestra -> condición) se toma del archivo de anotación
#' `annot_path` (columnas `Column`, `Condition`, y opcional `Experiment`). Las
#' intensidades (`Abundance: <muestra>`) y las métricas por muestra
#' (`# PSMs (by Search Engine)`, `# Peptides (by Search Engine)`, `Score Mascot`)
#' se mapean POR NOMBRE de muestra.
#'
#' @param file_path Ruta al TSV de datos exportado de Proteome Discoverer.
#' @param annot_path Ruta al TSV de anotación con columnas `Column`, `Condition`
#'   (obligatorias) y `Experiment` (opcional).
#' @param condition_order Vector de caracteres con el orden de las condiciones.
#'   Si es `NULL` (default), se deriva de `Condition` (orden de aparición). Si se
#'   provee, fija los niveles del factor y descarta muestras cuya condición no
#'   esté en la lista.
#' @param export_dir Directorio para exportar TSVs. `NULL` (default) = no exporta.
#' @param timestamp_suffix Lógico. Si `TRUE` (default), añade timestamp a los
#'   archivos exportados.
#' @param verbose Lógico. Si `TRUE` (default), muestra mensajes de progreso.
#'
#' @return Lista con clase `c("lfq_data", "proteomics_data", "list")` conteniendo
#'   `metadata`, `protein_id` y `protein_quant`.
#'
#' @examples
#' \dontrun{
#' result <- preprocess_lfq(
#'   file_path  = "data/20260710_AIturrate_2659_LFQ_QUANT_onlyRAW.tsv",
#'   annot_path = "data/20260710_AIturrate_2659_LFQ_QUANT_onlyRAW_Annot.tsv",
#'   export_dir = "./results"
#' )
#' print(result)
#' }
#'
#' @export
preprocess_lfq <- function(
    file_path,
    annot_path,
    condition_order = NULL,
    export_dir = NULL,
    timestamp_suffix = TRUE,
    verbose = TRUE
) {

  # --- Validación de argumentos ---
  if (!file.exists(file_path)) stop("Archivo de datos no encontrado: ", file_path)
  if (!file.exists(annot_path)) stop("Archivo de anotación no encontrado: ", annot_path)
  if (!is.null(condition_order) &&
      (length(condition_order) == 0 || !is.character(condition_order))) {
    stop("condition_order debe ser NULL o un vector de caracteres no vacío.")
  }

  # --- Lectura (check.names = FALSE para preservar headers de PD) ---
  if (verbose) message("Leyendo archivo de datos: ", basename(file_path))
  df <- read.delim(file_path, header = TRUE, sep = "\t",
                   stringsAsFactors = FALSE, check.names = FALSE)
  .validate_pd_columns(df)

  if (verbose) message("Leyendo archivo de anotación: ", basename(annot_path))
  annot <- read.delim(annot_path, header = TRUE, sep = "\t",
                      stringsAsFactors = FALSE, check.names = FALSE)
  .validate_lfq_annot(annot)
  annot$Column    <- as.character(annot$Column)
  annot$Condition <- as.character(annot$Condition)
  annot$Experiment <- if ("Experiment" %in% names(annot)) {
    as.character(annot$Experiment)
  } else {
    NA_character_
  }

  # --- Resolver orden de condiciones y filtrar diseño ---
  if (is.null(condition_order)) {
    condition_order <- unique(annot$Condition)
  } else {
    keep_ann <- annot$Condition %in% condition_order
    if (!any(keep_ann)) {
      stop("Ninguna condición del _Annot coincide con condition_order = c(",
           paste0("'", condition_order, "'", collapse = ", "), ").\n",
           "Condiciones en el _Annot: ", paste(unique(annot$Condition), collapse = ", "))
    }
    annot <- annot[keep_ann, , drop = FALSE]
  }

  # --- Tabla de muestras (ordenada por condición y orden del Annot) ---
  annot$R.Condition <- factor(annot$Condition, levels = condition_order,
                              ordered = TRUE)
  annot <- annot[order(annot$R.Condition, seq_len(nrow(annot))), , drop = FALSE]
  # Replicado = índice secuencial dentro de cada condición (robusto a nombres)
  annot$R.Replicate <- as.integer(
    stats::ave(seq_len(nrow(annot)), annot$R.Condition,
               FUN = function(i) seq_along(i))
  )
  coding_levels <- annot$Column   # Coding = nombre de muestra del _Annot

  # --- Verificar que cada muestra del diseño tiene su columna Abundance ---
  abund_map <- .lfq_sample_cols(df, "^Abundance:\\s*", coding_levels)
  missing_abund <- setdiff(coding_levels, names(abund_map))
  if (length(missing_abund) > 0) {
    stop("No se encontró columna 'Abundance:' para estas muestras del _Annot:\n  - ",
         paste(missing_abund, collapse = "\n  - "))
  }
  # Avisar de columnas Abundance del TSV que no están en el diseño
  all_abund <- sub("^Abundance:\\s*", "",
                   grep("^Abundance:\\s*", names(df), value = TRUE))
  extra_abund <- setdiff(all_abund, coding_levels)
  if (length(extra_abund) > 0 && verbose) {
    message("Aviso: columnas 'Abundance:' ignoradas (no están en el _Annot): ",
            paste(extra_abund, collapse = ", "))
  }

  # ==========================================================================
  # metadata
  # ==========================================================================
  if (verbose) message("Generando metadata de muestras...")
  run_summary <- data.frame(
    R.FileName   = unname(abund_map[coding_levels]),
    R.Condition  = annot$R.Condition,
    R.Replicate  = annot$R.Replicate,
    Coding       = coding_levels,
    R.Experiment = annot$Experiment,
    stringsAsFactors = FALSE
  )
  rownames(run_summary) <- run_summary$Coding

  # ==========================================================================
  # Información base de proteínas
  # ==========================================================================
  if (verbose) message("Procesando información base de proteínas...")

  col_mw     <- .pd_resolve_col(df, "MW [kDa]", "MW (kDa)", "MW")
  col_pi     <- .pd_resolve_col(df, "calc. pI", "calc pI", "Calculated pI")
  col_master <- .pd_resolve_col(df, "Master")
  col_pgids  <- .pd_resolve_col(df, "Protein Group IDs")
  col_fdr    <- .pd_resolve_col(df, "Protein FDR Confidence: Combined",
                                    "Protein FDR Confidence")
  col_qvalue <- .pd_resolve_col(df, "Exp. q-value: Combined", "Exp. q-value")
  col_pep    <- .pd_resolve_col(df, "Sum PEP Score")
  col_psms   <- .pd_resolve_col(df, "# PSMs", "Number of PSMs")
  col_pepts  <- .pd_resolve_col(df, "# Peptides", "Number of Peptides")
  col_uniq   <- .pd_resolve_col(df, "# Unique Peptides", "Number of Unique Peptides")
  col_cov    <- .pd_resolve_col(df, "Coverage [%]", "Coverage (%)", "Coverage")

  base_info <- data.frame(
    PG.ProteinGroups       = as.character(df$Accession),
    PG.ProteinDescriptions = as.character(df$Description),
    PG.Genes               = .parse_gene_from_description(df$Description),
    PG.MolecularWeight     = .pd_numeric_or_na(df, col_mw),
    stringsAsFactors = FALSE
  )

  lfq_extras <- data.frame(
    calc.pI                = .pd_numeric_or_na(df, col_pi),
    Master                 = .pd_char_or_na(df, col_master),
    Protein.Group.IDs      = .pd_char_or_na(df, col_pgids),
    Protein.FDR.Confidence = .pd_char_or_na(df, col_fdr),
    Exp.q.value            = .pd_numeric_or_na(df, col_qvalue),
    NrUniquePeptides       = .pd_numeric_or_na(df, col_uniq),
    Species                = .parse_species_from_description(df$Description),
    stringsAsFactors = FALSE
  )

  # Métricas globales
  g_psms <- .pd_numeric_or_na(df, col_psms)
  g_pept <- .pd_numeric_or_na(df, col_pepts)
  g_cov  <- .pd_numeric_or_na(df, col_cov)
  g_pep  <- .pd_numeric_or_na(df, col_pep)

  # --- Helper: matriz por muestra (mapeada por nombre) para una familia ---
  sample_matrix <- function(prefix_regex, out_prefix) {
    map <- .lfq_sample_cols(df, prefix_regex, coding_levels)
    mat <- matrix(NA_real_, nrow = nrow(df), ncol = length(coding_levels),
                  dimnames = list(NULL, paste0(out_prefix, "_", coding_levels)))
    for (i in seq_along(coding_levels)) {
      col <- map[[coding_levels[i]]]
      if (!is.null(col) && !is.na(col)) {
        mat[, i] <- suppressWarnings(as.numeric(df[[col]]))
      }
    }
    as.data.frame(mat, stringsAsFactors = FALSE, check.names = FALSE)
  }
  # Réplica de una métrica global por muestra (para shape wide)
  wide_global <- function(values, out_prefix) {
    mat <- matrix(rep(values, length(coding_levels)),
                  nrow = length(values), ncol = length(coding_levels))
    colnames(mat) <- paste0(out_prefix, "_", coding_levels)
    as.data.frame(mat, stringsAsFactors = FALSE, check.names = FALSE)
  }

  # Familias por muestra
  psms_by <- sample_matrix("^# PSMs \\(by Search Engine\\):\\s*",
                           "PG.NrOfPrecursorsUsedForQuantification")
  pept_by <- sample_matrix("^# Peptides \\(by Search Engine\\):\\s*",
                           "PG.NrOfStrippedSequencesUsedForQuantification")
  mascot_by <- sample_matrix("^Score Mascot:\\s*", "PG.Cscore.RunWise")

  abund_mat <- df[, unname(abund_map[coding_levels]), drop = FALSE]
  abund_mat[] <- lapply(abund_mat, function(v) suppressWarnings(as.numeric(v)))
  names(abund_mat) <- paste0("PG.Quantity_", coding_levels)

  # ==========================================================================
  # protein_ID  (métricas de identificación por muestra donde existan)
  # ==========================================================================
  if (verbose) message("Procesando protein_ID...")

  id_psms <- psms_by; names(id_psms) <- paste0("PG.NrOfPrecursorsIdentified_", coding_levels)
  id_pept <- pept_by; names(id_pept) <- paste0("PG.NrOfStrippedSequencesIdentified_", coding_levels)
  id_cov  <- wide_global(g_cov, "PG.Coverage")
  id_score <- mascot_by  # PG.Cscore.RunWise_<coding>

  protein_ID <- cbind(base_info, lfq_extras, id_psms, id_pept, id_cov, id_score)

  static_cols  <- names(base_info)
  extras_cols  <- names(lfq_extras)
  metric_order <- c("PG.NrOfPrecursorsIdentified",
                    "PG.NrOfStrippedSequencesIdentified",
                    "PG.Coverage", "PG.Cscore.RunWise")
  runwise_order <- unlist(lapply(metric_order,
                                 function(m) paste0(m, "_", coding_levels)))
  protein_ID <- protein_ID[, c(static_cols, extras_cols,
                               intersect(runwise_order, names(protein_ID))),
                           drop = FALSE]
  protein_ID <- protein_ID[order(protein_ID$PG.ProteinGroups), , drop = FALSE]
  rownames(protein_ID) <- NULL
  if (anyDuplicated(protein_ID$PG.ProteinGroups) > 0) {
    stop("Error de integridad: protein_ID contiene Accession duplicados.")
  }

  # ==========================================================================
  # protein_QUANT
  # ==========================================================================
  if (verbose) message("Procesando protein_QUANT...")

  global_metrics <- data.frame(
    PG.NrOfPrecursorsIdentified.Global        = g_psms,
    PG.NrOfStrippedSequencesIdentified.Global = g_pept,
    PG.Coverage.Global                        = g_cov,
    PG.Cscore                                 = g_pep,
    stringsAsFactors = FALSE
  )

  protein_QUANT <- cbind(base_info, global_metrics, psms_by, pept_by, abund_mat)

  static_cols2 <- names(base_info)
  global_cols  <- names(global_metrics)
  q_metric_order <- c("PG.NrOfPrecursorsUsedForQuantification",
                      "PG.NrOfStrippedSequencesUsedForQuantification",
                      "PG.Quantity")
  q_runwise <- unlist(lapply(q_metric_order,
                             function(m) paste0(m, "_", coding_levels)))
  protein_QUANT <- protein_QUANT[, c(static_cols2, global_cols,
                                     intersect(q_runwise, names(protein_QUANT))),
                                 drop = FALSE]
  protein_QUANT <- protein_QUANT[order(protein_QUANT$PG.ProteinGroups), , drop = FALSE]
  rownames(protein_QUANT) <- NULL
  if (anyDuplicated(protein_QUANT$PG.ProteinGroups) > 0) {
    stop("Error de integridad: protein_QUANT contiene Accession duplicados.")
  }

  # ==========================================================================
  # Exportación opcional
  # ==========================================================================
  if (!is.null(export_dir)) {
    if (verbose) message("Exportando archivos a: ", export_dir)
    if (!dir.exists(export_dir)) dir.create(export_dir, recursive = TRUE)
    suffix <- if (timestamp_suffix) paste0("_", format(Sys.time(), "%Y%m%d_%H%M%S")) else ""
    readr::write_tsv(run_summary,
                     file.path(export_dir, paste0("Metadata", suffix, ".tsv")), na = "")
    readr::write_tsv(protein_ID,
                     file.path(export_dir, paste0("Protein_ID", suffix, ".tsv")), na = "")
    readr::write_tsv(protein_QUANT,
                     file.path(export_dir, paste0("Protein_QUANT", suffix, ".tsv")), na = "")
    if (verbose) message("Archivos exportados exitosamente.")
  }

  # ==========================================================================
  # Resultado
  # ==========================================================================
  if (verbose) {
    message("Procesamiento completado:\n",
            "  - Muestras: ", nrow(run_summary), "\n",
            "  - Proteínas (ID): ", nrow(protein_ID), "\n",
            "  - Proteínas (QUANT): ", nrow(protein_QUANT))
  }

  result <- list(
    metadata      = as.data.frame(run_summary),
    protein_id    = as.data.frame(protein_ID),
    protein_quant = as.data.frame(protein_QUANT)
  )
  class(result) <- c("lfq_data", "proteomics_data", "list")
  return(result)
}

# =============================================================================
# Métodos para clase lfq_data
# =============================================================================

#' @export
print.lfq_data <- function(x, ...) {
  cat("Datos LFQ (Proteome Discoverer) preprocesados\n")
  cat("---------------------------------------------\n")
  cat("Muestras (metadata):", nrow(x$metadata), "\n")
  cat("Proteínas (ID):", nrow(x$protein_id), "\n")
  cat("Proteínas (QUANT):", nrow(x$protein_quant), "\n")
  cat("\nCondiciones:",
      paste(levels(x$metadata$R.Condition) %||%
              unique(x$metadata$R.Condition), collapse = ", "), "\n")
  invisible(x)
}

# =============================================================================
# Ejemplos de Uso (no ejecutar)
# =============================================================================
if (FALSE) {
  result <- preprocess_lfq(
    file_path  = "data/20260710_AIturrate_2659_LFQ_QUANT_onlyRAW.tsv",
    annot_path = "data/20260710_AIturrate_2659_LFQ_QUANT_onlyRAW_Annot.tsv",
    export_dir = "./results"
  )
  print(result)

  source("R/Processing.R")
  res <- process_proteomics(
    preprocessing = result,
    norm_method   = "cycloess",
    imp_method    = "combo",
    control       = "WT",
    export_dir    = "./results_lfq"
  )
}
