# =============================================================================
# Preprocesamiento de Datos TMT de Proteome Discoverer
# =============================================================================
#
# Convierte exports de proteínas de Proteome Discoverer (TMT/TMTpro) a la
# misma estructura que `preprocess_spectronaut()` para alimentar el pipeline
# downstream sin cambios (Processing.R y módulos asociados).
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

# --- Operador %||% ---
if (!exists("%||%", mode = "function")) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
}

# =============================================================================
# Funciones Auxiliares Internas
# =============================================================================

#' Extrae Gene Name de la columna Description (formato UniProt embebido)
#' @description Busca el patrón `GN=<gene>` (separado por espacios) habitual en
#'   los exports de Proteome Discoverer.
#' @noRd
.parse_gene_from_description <- function(x) {
  m <- stringr::str_match(x, "GN=([^ ]+)")
  m[, 2]
}

#' Detecta columnas `Abundance:` y deriva Coding / Condition / Replicate
#' @noRd
.parse_abundance_columns <- function(df) {
  abund_cols <- grep("^Abundance:\\s*", names(df), value = TRUE)
  if (length(abund_cols) == 0) {
    stop("No se encontraron columnas 'Abundance:' en el archivo. ",
         "Verifica que sea un export de Proteome Discoverer.")
  }

  coding <- sub("^Abundance:\\s*", "", abund_cols)
  m <- stringr::str_match(coding, "^(.+)_(\\d+)$")
  if (any(is.na(m[, 1]))) {
    bad <- coding[is.na(m[, 1])]
    stop(
      "No se pudieron parsear los sufijos de Abundance (esperado <Condicion>_<Replicado>):\n  - ",
      paste(bad, collapse = "\n  - ")
    )
  }

  data.frame(
    abundance_col = abund_cols,
    Coding        = coding,
    R.Condition   = m[, 2],
    R.Replicate   = as.integer(m[, 3]),
    stringsAsFactors = FALSE
  )
}

#' Valida columnas mínimas requeridas para un export de Proteome Discoverer
#' @noRd
.validate_pd_columns <- function(df) {
  required <- c("Accession", "Description")
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop(
      "Columnas requeridas faltantes en el archivo:\n  - ",
      paste(missing, collapse = "\n  - "),
      "\nVerifica que el archivo sea un export de proteínas de Proteome Discoverer."
    )
  }
  invisible(TRUE)
}

#' Lookup tolerante a variantes habituales de nombres de columna de PD
#' @description PD exporta columnas con caracteres especiales (espacios, #, %,
#'   corchetes, ":"). Esta función prueba varias variantes y devuelve la primera
#'   que exista, o `NA_character_` si ninguna coincide.
#' @noRd
.tmt_resolve_col <- function(df, ...) {
  candidates <- unlist(list(...), use.names = FALSE)
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) NA_character_ else hit[1]
}

#' Extrae un vector numérico de `df[[col]]` o NA si la columna no existe
#' @noRd
.tmt_numeric_or_na <- function(df, col) {
  if (is.na(col) || !col %in% names(df)) {
    return(rep(NA_real_, nrow(df)))
  }
  suppressWarnings(as.numeric(df[[col]]))
}

#' Extrae un vector character de `df[[col]]` o NA si la columna no existe
#' @noRd
.tmt_char_or_na <- function(df, col) {
  if (is.na(col) || !col %in% names(df)) {
    return(rep(NA_character_, nrow(df)))
  }
  as.character(df[[col]])
}

# =============================================================================
# Función Principal
# =============================================================================

#' Preprocesa exports TMT de Proteome Discoverer
#'
#' @description
#' Convierte un export de proteínas de Proteome Discoverer (formato TSV ancho,
#' con columnas `Abundance: <Condicion>_<Replicado>`) a la misma estructura de
#' salida que `preprocess_spectronaut()`: tres data.frames (`metadata`,
#' `protein_id`, `protein_quant`) listos para el pipeline downstream
#' (`process_proteomics()` y módulos asociados).
#'
#' Las métricas globales de identificación de PD (`# PSMs`, `# Peptides`,
#' `# Unique Peptides`, `Coverage [%]`) se replican en columnas por canal para
#' encajar con el contrato wide de Spectronaut.
#'
#' @param file_path Ruta al TSV exportado de Proteome Discoverer.
#' @param condition_order Vector de caracteres con el orden de las condiciones
#'   experimentales (ej: `c("A","B","C","D","IS")`). Solo se conservan los
#'   canales cuya condición esté en este vector — útil para excluir Internal
#'   Standards omitiendo `"IS"`.
#' @param export_dir Directorio para exportar archivos TSV. Si es `NULL`
#'   (default), no se exportan archivos.
#' @param timestamp_suffix Lógico. Si `TRUE` (default), añade timestamp a los
#'   nombres de los archivos exportados.
#' @param verbose Lógico. Si `TRUE` (default), muestra mensajes de progreso.
#'
#' @return Lista con clase `c("tmt_data", "proteomics_data", "list")`
#'   conteniendo:
#'   \describe{
#'     \item{metadata}{Data frame con un registro por canal/muestra}
#'     \item{protein_id}{Data frame con métricas de identificación por proteína
#'       (formato wide, métricas globales replicadas por canal)}
#'     \item{protein_quant}{Data frame con métricas de cuantificación por
#'       proteína (incluye `PG.Quantity_<Coding>`)}
#'   }
#'
#' @examples
#' \dontrun{
#' result <- preprocess_tmt(
#'   file_path = "data/20260527_Q25_TMTpro_TMT1y2_10Fr_Static_3engines_onlyRAW.tsv",
#'   condition_order = c("A", "B", "C", "D", "IS"),
#'   export_dir = "./results"
#' )
#' print(result)
#' }
#'
#' @export
preprocess_tmt <- function(
    file_path,
    condition_order,
    export_dir = NULL,
    timestamp_suffix = TRUE,
    verbose = TRUE
) {

  # --- Validación de argumentos ---
  if (!file.exists(file_path)) {
    stop("Archivo no encontrado: ", file_path)
  }

  if (length(condition_order) == 0 || !is.character(condition_order)) {
    stop("condition_order debe ser un vector de caracteres no vacío.")
  }

  # --- Lectura del archivo (check.names = FALSE para preservar headers PD) ---
  if (verbose) message("Leyendo archivo: ", basename(file_path))

  df <- read.delim(
    file_path,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  # --- Validar columnas mínimas ---
  .validate_pd_columns(df)

  # --- Parsear columnas Abundance ---
  abund <- .parse_abundance_columns(df)

  # Filtrar canales según condition_order (excluye p.ej. "IS" si no está)
  keep <- abund$R.Condition %in% condition_order
  if (!any(keep)) {
    stop(
      "Ningún canal Abundance coincide con condition_order = c(",
      paste0("'", condition_order, "'", collapse = ", "), ").\n",
      "Condiciones detectadas en el archivo: ",
      paste(unique(abund$R.Condition), collapse = ", ")
    )
  }
  abund <- abund[keep, , drop = FALSE]

  # Orden final por condition_order y replicate
  abund$R.Condition <- factor(abund$R.Condition,
                              levels = condition_order, ordered = TRUE)
  abund <- abund[order(abund$R.Condition, abund$R.Replicate), , drop = FALSE]
  coding_levels <- abund$Coding

  # --- Resolver columnas opcionales de PD (nombres tolerantes) ---
  col_mw       <- .tmt_resolve_col(df, "MW [kDa]", "MW (kDa)", "MW")
  col_pi       <- .tmt_resolve_col(df, "calc. pI", "calc pI", "Calculated pI")
  col_master   <- .tmt_resolve_col(df, "Master")
  col_pgids    <- .tmt_resolve_col(df, "Protein Group IDs")
  col_fdr      <- .tmt_resolve_col(df, "Protein FDR Confidence: Combined",
                                       "Protein FDR Confidence")
  col_qvalue   <- .tmt_resolve_col(df, "Exp. q-value: Combined",
                                       "Exp. q-value")
  col_pep      <- .tmt_resolve_col(df, "Sum PEP Score")
  col_mascot   <- .tmt_resolve_col(df, "Score Mascot: Mascot", "Score: Mascot")
  col_sequest  <- .tmt_resolve_col(df, "Score Sequest HT: Sequest HT",
                                       "Score: Sequest HT")
  col_psms     <- .tmt_resolve_col(df, "# PSMs", "Number of PSMs")
  col_pepts    <- .tmt_resolve_col(df, "# Peptides", "Number of Peptides")
  col_uniq     <- .tmt_resolve_col(df, "# Unique Peptides",
                                       "Number of Unique Peptides")
  col_cov      <- .tmt_resolve_col(df, "Coverage [%]", "Coverage (%)",
                                       "Coverage")

  # ==========================================================================
  # metadata
  # ==========================================================================
  if (verbose) message("Generando metadata de canales...")

  run_summary <- data.frame(
    R.FileName  = abund$abundance_col,
    R.Condition = abund$R.Condition,
    R.Replicate = abund$R.Replicate,
    Coding      = abund$Coding,
    stringsAsFactors = FALSE
  )
  rownames(run_summary) <- run_summary$Coding

  # ==========================================================================
  # Información base de proteínas (compartida entre protein_ID y protein_QUANT)
  # ==========================================================================
  if (verbose) message("Procesando información base de proteínas...")

  base_info <- data.frame(
    PG.ProteinGroups       = as.character(df$Accession),
    PG.ProteinDescriptions = as.character(df$Description),
    PG.Genes               = .parse_gene_from_description(df$Description),
    PG.MolecularWeight     = .tmt_numeric_or_na(df, col_mw),
    stringsAsFactors = FALSE
  )

  # Extras TMT-específicos (se mantienen tras MW para no romper downstream)
  tmt_extras <- data.frame(
    calc.pI                = .tmt_numeric_or_na(df, col_pi),
    Master                 = .tmt_char_or_na(df, col_master),
    Protein.Group.IDs      = .tmt_char_or_na(df, col_pgids),
    Protein.FDR.Confidence = .tmt_char_or_na(df, col_fdr),
    Exp.q.value            = .tmt_numeric_or_na(df, col_qvalue),
    Score.Mascot           = .tmt_numeric_or_na(df, col_mascot),
    Score.SequestHT        = .tmt_numeric_or_na(df, col_sequest),
    NrUniquePeptides       = .tmt_numeric_or_na(df, col_uniq),
    stringsAsFactors = FALSE
  )

  # Métricas globales que se replicarán por canal
  g_psms <- .tmt_numeric_or_na(df, col_psms)
  g_pept <- .tmt_numeric_or_na(df, col_pepts)
  g_cov  <- .tmt_numeric_or_na(df, col_cov)
  g_pep  <- .tmt_numeric_or_na(df, col_pep)
  g_uniq <- .tmt_numeric_or_na(df, col_uniq)

  # ==========================================================================
  # protein_ID  (métricas globales replicadas por canal para shape wide)
  # ==========================================================================
  if (verbose) message("Procesando protein_ID...")

  mk_wide <- function(values, metric_name) {
    mat <- matrix(rep(values, length(coding_levels)),
                  nrow = length(values), ncol = length(coding_levels))
    colnames(mat) <- paste0(metric_name, "_", coding_levels)
    as.data.frame(mat, stringsAsFactors = FALSE)
  }

  id_psms <- mk_wide(g_psms, "PG.NrOfPrecursorsIdentified")
  id_pept <- mk_wide(g_pept, "PG.NrOfStrippedSequencesIdentified")
  id_cov  <- mk_wide(g_cov,  "PG.Coverage")
  id_pep  <- mk_wide(g_pep,  "PG.Cscore.RunWise")

  protein_ID <- cbind(base_info, tmt_extras, id_psms, id_pept, id_cov, id_pep)

  # Ordenar columnas: estáticas + extras + métricas en orden estable
  static_cols  <- c("PG.ProteinGroups", "PG.ProteinDescriptions",
                    "PG.Genes", "PG.MolecularWeight")
  extras_cols  <- names(tmt_extras)
  metric_order <- c("PG.NrOfPrecursorsIdentified",
                    "PG.NrOfStrippedSequencesIdentified",
                    "PG.Coverage", "PG.Cscore.RunWise")
  runwise_order <- unlist(lapply(metric_order,
                                 function(m) paste0(m, "_", coding_levels)))
  final_cols <- c(static_cols, extras_cols,
                  intersect(runwise_order, names(protein_ID)))
  protein_ID <- protein_ID[, final_cols, drop = FALSE]
  protein_ID <- protein_ID[order(protein_ID$PG.ProteinGroups), , drop = FALSE]
  rownames(protein_ID) <- NULL

  # Validar unicidad
  if (anyDuplicated(protein_ID$PG.ProteinGroups) > 0) {
    stop("Error de integridad: protein_ID contiene Accession duplicados. ",
         "Revisa los datos de entrada.")
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

  # Métricas por canal: NrOfPrecursorsUsedForQuantification (PSMs replicado),
  # NrOfStrippedSequencesUsedForQuantification (Unique Peptides replicado),
  # Quantity (Abundance: <Coding>)
  q_psms <- mk_wide(g_psms, "PG.NrOfPrecursorsUsedForQuantification")
  q_uniq <- mk_wide(g_uniq, "PG.NrOfStrippedSequencesUsedForQuantification")

  abund_mat <- df[, abund$abundance_col, drop = FALSE]
  abund_mat[] <- lapply(abund_mat, function(v) suppressWarnings(as.numeric(v)))
  names(abund_mat) <- paste0("PG.Quantity_", coding_levels)

  protein_QUANT <- cbind(base_info, global_metrics, q_psms, q_uniq, abund_mat)

  static_cols2 <- c("PG.ProteinGroups", "PG.ProteinDescriptions",
                    "PG.Genes", "PG.MolecularWeight")
  global_cols  <- c("PG.NrOfPrecursorsIdentified.Global",
                    "PG.NrOfStrippedSequencesIdentified.Global",
                    "PG.Coverage.Global", "PG.Cscore")
  q_metric_order <- c("PG.NrOfPrecursorsUsedForQuantification",
                      "PG.NrOfStrippedSequencesUsedForQuantification",
                      "PG.Quantity")
  q_runwise <- unlist(lapply(q_metric_order,
                             function(m) paste0(m, "_", coding_levels)))
  final_cols2 <- c(static_cols2, global_cols,
                   intersect(q_runwise, names(protein_QUANT)))
  protein_QUANT <- protein_QUANT[, final_cols2, drop = FALSE]
  protein_QUANT <- protein_QUANT[order(protein_QUANT$PG.ProteinGroups), ,
                                 drop = FALSE]
  rownames(protein_QUANT) <- NULL

  if (anyDuplicated(protein_QUANT$PG.ProteinGroups) > 0) {
    stop("Error de integridad: protein_QUANT contiene Accession duplicados. ",
         "Revisa los datos de entrada.")
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
      "  - Canales: ", nrow(run_summary), "\n",
      "  - Proteínas (ID): ", nrow(protein_ID), "\n",
      "  - Proteínas (QUANT): ", nrow(protein_QUANT)
    )
  }

  result <- list(
    metadata      = as.data.frame(run_summary),
    protein_id    = as.data.frame(protein_ID),
    protein_quant = as.data.frame(protein_QUANT)
  )

  class(result) <- c("tmt_data", "proteomics_data", "list")
  return(result)
}

# =============================================================================
# Métodos para clase tmt_data
# =============================================================================

#' @export
print.tmt_data <- function(x, ...) {
  cat("Datos TMT (Proteome Discoverer) preprocesados\n")
  cat("---------------------------------------------\n")
  cat("Canales (metadata):", nrow(x$metadata), "\n")
  cat("Proteínas (ID):", nrow(x$protein_id), "\n")
  cat("Proteínas (QUANT):", nrow(x$protein_quant), "\n")
  cat("\nCondiciones:",
      paste(levels(x$metadata$R.Condition) %||%
              unique(x$metadata$R.Condition), collapse = ", "),
      "\n")
  invisible(x)
}

# =============================================================================
# Ejemplos de Uso (no ejecutar)
# =============================================================================
if (FALSE) {
  # Ejemplo 1: incluyendo canales IS como condición separada
  result <- preprocess_tmt(
    file_path = "data/20260527_Q25_TMTpro_TMT1y2_10Fr_Static_3engines_onlyRAW.tsv",
    condition_order = c("A", "B", "C", "D", "IS")
  )
  print(result)

  # Ejemplo 2: descartando canales IS y exportando
  result <- preprocess_tmt(
    file_path = "data/20260527_Q25_TMTpro_TMT1y2_10Fr_Static_3engines_onlyRAW.tsv",
    condition_order = c("A", "B", "C", "D"),
    export_dir = "./results",
    timestamp_suffix = TRUE
  )

  # Encadenar con el pipeline downstream sin cambios
  source("R/Processing.R")
  res <- process_proteomics(
    preprocessing = result,
    norm_method   = "cycloess",
    imp_method    = "combo",
    export_dir    = "./results_tmt"
  )
}
