# =============================================================================
# DIA-NN Data Preprocessing
# =============================================================================
#
# Converts a DIA-NN protein-group matrix (`report.pg_matrix.tsv`) into the same
# structure as `preprocess_spectronaut()`, so the downstream pipeline
# (Processing.R and the associated modules) consumes it unchanged.
#
# DIA-NN is a DIA search engine like Spectronaut, but its matrix export is wide
# -- one row per protein group, one column per run -- so the shape of this file
# follows Preprocessing_TMT.R rather than Preprocessing.R.
#
# Copyright 2025 Sergio Ciordia
# License: GPL-3
# =============================================================================

# --- Dependencies ---

# =============================================================================
# Internal Helper Functions
# =============================================================================

# .parse_abundance_columns() and .nadia_stamp() live in R/utils.R. The first is
# shared with preprocess_tmt(): both formats can declare the design in the
# column names, using the same `Abundance: <Condition>_<Replicate>` convention.

#' Annotation columns of a DIA-NN matrix export
#'
#' The union of the columns DIA-NN 1.8 and 2.x write before the run columns.
#' Anything outside this set is taken to be a run, which is what makes the
#' detection survive a version that adds or drops an annotation column.
#'
#' @noRd
.DIANN_ANNOT_COLS <- c(
  "Protein.Group", "Protein.Ids", "Protein.Names", "Genes",
  "First.Protein.Description", "N.Sequences", "N.Proteotypic.Sequences"
)

#' Extract a numeric vector from `df[[col]]`, or NA if the column is absent
#' @noRd
.diann_numeric_or_na <- function(df, col) {
  if (is.na(col) || !col %in% names(df)) {
    return(rep(NA_real_, nrow(df)))
  }
  suppressWarnings(as.numeric(df[[col]]))
}

#' Extract a character vector from `df[[col]]`, or NA if the column is absent
#' @noRd
.diann_char_or_na <- function(df, col) {
  if (is.na(col) || !col %in% names(df)) {
    return(rep(NA_character_, nrow(df)))
  }
  as.character(df[[col]])
}

#' Work out which columns hold intensities, and how they map to samples
#'
#' Two routes, and deliberately no third one that guesses. Either the columns
#' already carry the `Abundance: <Condition>_<Replicate>` convention, or the
#' caller names them in file order through `sample_names`. Inferring the design
#' from the raw-file name would work on one export and fail silently on the
#' next, which is the worst of both.
#'
#' @param df           Data frame read from the export.
#' @param sample_names Character vector in the order of the run columns, or
#'   `NULL`.
#' @return A data frame with `abundance_col` (the original header),
#'   `Coding`, `R.Condition` and `R.Replicate`.
#' @noRd
.diann_resolve_design <- function(df, sample_names = NULL) {
  run_cols <- setdiff(names(df), .DIANN_ANNOT_COLS)
  if (length(run_cols) == 0) {
    stop("No run columns found in the file: every column is a known DIA-NN ",
         "annotation column. Check that it is a protein-group matrix export ",
         "(report.pg_matrix.tsv).")
  }

  if (is.null(sample_names)) {
    if (!any(grepl("^Abundance:\\s*", names(df)))) {
      stop(
        "Could not work out which sample each run column belongs to.\n",
        "Either rename them to the 'Abundance: <Condition>_<Replicate>' ",
        "convention, or pass sample_names in this order:\n  - ",
        paste(run_cols, collapse = "\n  - ")
      )
    }
    return(.parse_abundance_columns(df))
  }

  if (!is.character(sample_names)) {
    stop("sample_names must be a character vector.")
  }
  if (length(sample_names) != length(run_cols)) {
    stop("sample_names has ", length(sample_names), " values but the file has ",
         length(run_cols), " run columns. They are, in order:\n  - ",
         paste(run_cols, collapse = "\n  - "))
  }
  if (anyDuplicated(sample_names) > 0) {
    stop("sample_names has duplicated values:\n  - ",
         paste(unique(sample_names[duplicated(sample_names)]),
               collapse = "\n  - "))
  }

  m <- stringr::str_match(sample_names, "^(.+)_(\\d+)$")
  if (any(is.na(m[, 1]))) {
    stop(
      "Could not parse these sample_names (expected <Condition>_<Replicate>):\n  - ",
      paste(sample_names[is.na(m[, 1])], collapse = "\n  - ")
    )
  }

  data.frame(
    abundance_col = run_cols,
    Coding        = sample_names,
    R.Condition   = m[, 2],
    R.Replicate   = as.integer(m[, 3]),
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# Main Function
# =============================================================================

#' Preprocess DIA-NN protein-group matrices
#'
#' @description
#' Converts a DIA-NN protein-group matrix (`report.pg_matrix.tsv`, wide TSV
#' format) into the same output structure as `preprocess_spectronaut()`: three
#' data.frames (`metadata`, `protein_id`, `protein_quant`) ready for the
#' downstream pipeline ([process_proteomics()] and the associated modules).
#'
#' DIA-NN names its intensity columns after the raw file, which says nothing
#' about the experimental design, so the design has to be declared. There are
#' two ways, and neither needs a separate annotation file: rename the columns to
#' `Abundance: <Condition>_<Replicate>` as Proteome Discoverer writes them, or
#' pass `sample_names` in the order the run columns appear in the file. The
#' design is never inferred from the raw-file name -- a guess that happened to
#' be right for one export would be wrong and silent for the next.
#'
#' A protein-group matrix carries only intensities. The identification metrics
#' the wide contract expects are filled from what DIA-NN does report
#' (`N.Sequences`, `N.Proteotypic.Sequences`), left as `NA` where it reports
#' nothing (the precursor counts), or omitted where nothing downstream needs
#' them (coverage and the Spectronaut C-score).
#'
#' @param file_path Path to the `report.pg_matrix.tsv` exported by DIA-NN.
#' @param condition_order Character vector with the order of the experimental
#'   conditions (e.g. `c("A", "B", "D")`). Only the runs whose condition is in
#'   this vector are kept, so it doubles as a filter.
#' @param sample_names Character vector of `<Condition>_<Replicate>` labels, one
#'   per run column, **in the order the columns appear in the file**. Use it to
#'   read the report exactly as DIA-NN wrote it. If `NULL` (default), the
#'   columns are expected to follow the `Abundance:` convention.
#' @param export_dir Directory to export the TSV files to. If `NULL` (default),
#'   no files are exported.
#' @param timestamp_suffix Logical. If `TRUE` (default), appends a timestamp to
#'   the names of the exported files.
#' @param verbose Logical. If `TRUE` (default), shows progress messages.
#'
#' @return A list with class `c("diann_data", "proteomics_data", "list")`
#'   containing:
#'   \describe{
#'     \item{metadata}{Data frame with one record per run, including the number
#'       of protein groups quantified in each}
#'     \item{protein_id}{Data frame with the identification metrics per protein
#'       (wide format, the global metrics replicated per run)}
#'     \item{protein_quant}{Data frame with the quantification metrics per
#'       protein (includes `PG.Quantity_<Coding>`)}
#'   }
#'
#'   The call that produced the object and a fingerprint of the file it was read
#'   from (path, size, modification time and MD5) travel with it as the
#'   attributes `nadia_call` and `nadia_source`. They are attributes rather than
#'   list elements so that the three-element structure above is unchanged;
#'   [write_nadia()] records them as the provenance of an analysis.
#'
#' @examples
#' # A trimmed DIA-NN protein-group matrix ships with the package. Its run
#' # columns have been renamed to the "Abundance: <Condition>_<Replicate>"
#' # convention, so no further argument is needed.
#' report <- system.file("extdata", "nadia_diann_report.tsv.gz",
#'                       package = "NADIA")
#'
#' diann <- preprocess_diann(report, condition_order = c("A", "B", "D"),
#'                           verbose = FALSE)
#' diann
#' table(diann$metadata$R.Condition)
#'
#' # Straight from DIA-NN, where the columns are still raw-file paths, the
#' # design is given instead by naming the runs in the order they appear:
#' # preprocess_diann("report.pg_matrix.tsv",
#' #                  condition_order = c("A", "B", "D"),
#' #                  sample_names = c("A_1", "A_2", "A_3", "A_4",
#' #                                   "B_1", "B_2", "B_3", "B_4",
#' #                                   "D_1", "D_2", "D_3", "D_4"))
#'
#' @export
preprocess_diann <- function(
    file_path,
    condition_order,
    sample_names = NULL,
    export_dir = NULL,
    timestamp_suffix = TRUE,
    verbose = TRUE
) {
  # ==========================================================================
  # Validation
  # ==========================================================================
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path)
  }
  if (missing(condition_order) || length(condition_order) == 0 ||
      !is.character(condition_order)) {
    stop("condition_order must be a non-empty character vector.")
  }

  if (verbose) message("Reading file: ", file_path)

  df <- utils::read.delim(file_path, header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE, check.names = FALSE)

  if (!"Protein.Group" %in% names(df)) {
    stop("Required column missing from the data file:\n  - Protein.Group\n",
         "Check that the file is a DIA-NN protein-group matrix ",
         "(report.pg_matrix.tsv).")
  }

  # ==========================================================================
  # The design
  # ==========================================================================
  runs <- .diann_resolve_design(df, sample_names)

  keep <- runs$R.Condition %in% condition_order
  if (!any(keep)) {
    stop(
      "No run matches condition_order = c(",
      paste0("'", condition_order, "'", collapse = ", "), ").\n",
      "Conditions detected in the file: ",
      paste(unique(runs$R.Condition), collapse = ", ")
    )
  }
  runs <- runs[keep, , drop = FALSE]

  runs$R.Condition <- factor(runs$R.Condition,
                             levels = condition_order, ordered = TRUE)
  runs <- runs[order(runs$R.Condition, runs$R.Replicate), , drop = FALSE]
  coding_levels <- runs$Coding

  # ==========================================================================
  # The intensities
  # ==========================================================================
  # DIA-NN writes an empty cell where a protein group was not quantified, so
  # the coercion is the step that turns absence into NA.
  quant <- df[, runs$abundance_col, drop = FALSE]
  quant[] <- lapply(quant, function(v) suppressWarnings(as.numeric(v)))
  names(quant) <- paste0("PG.Quantity_", coding_levels)

  # ==========================================================================
  # metadata
  # ==========================================================================
  if (verbose) message("Generating run metadata...")

  run_summary <- data.frame(
    R.FileName  = runs$abundance_col,
    R.Condition = runs$R.Condition,
    R.Replicate = runs$R.Replicate,
    Coding      = runs$Coding,
    # The one per-run count a protein-group matrix can give: how many groups
    # carry an intensity in that run. summary_list_widget() picks it up.
    R.ProteinGroupsIdentified = as.numeric(colSums(!is.na(quant))),
    stringsAsFactors = FALSE
  )
  rownames(run_summary) <- run_summary$Coding

  # ==========================================================================
  # Base protein information (shared between protein_ID and protein_QUANT)
  # ==========================================================================
  if (verbose) message("Processing the base protein information...")

  col_descr <- if ("First.Protein.Description" %in% names(df)) {
    "First.Protein.Description"
  } else {
    NA_character_
  }

  base_info <- data.frame(
    PG.ProteinGroups       = as.character(df$Protein.Group),
    PG.ProteinDescriptions = .diann_char_or_na(df, col_descr),
    PG.Genes               = .diann_char_or_na(df, "Genes"),
    # DIA-NN does not report a molecular weight. The column is kept, all NA, so
    # that the object is shape-compatible with the other three formats; every
    # consumer of it already guards against a missing value.
    PG.MolecularWeight     = rep(NA_real_, nrow(df)),
    stringsAsFactors = FALSE
  )

  # DIA-NN-specific extras (kept after MW so as not to break downstream code)
  diann_extras <- data.frame(
    PG.ProteinNames  = .diann_char_or_na(df, "Protein.Names"),
    Protein.Ids      = .diann_char_or_na(df, "Protein.Ids"),
    NrUniquePeptides = .diann_numeric_or_na(df, "N.Proteotypic.Sequences"),
    stringsAsFactors = FALSE
  )

  # Global metrics that will be replicated per run. `N.Sequences` counts every
  # peptide sequence and `N.Proteotypic.Sequences` only the unique ones, which
  # is the same distinction TMT draws between "# Peptides" and
  # "# Unique Peptides"; they map onto the same two families.
  g_pept <- .diann_numeric_or_na(df, "N.Sequences")
  g_uniq <- .diann_numeric_or_na(df, "N.Proteotypic.Sequences")
  g_prec <- rep(NA_real_, nrow(df))

  # ==========================================================================
  # protein_ID  (global metrics replicated per run to get the wide shape)
  # ==========================================================================
  if (verbose) message("Processing protein_ID...")

  mk_wide <- function(values, metric_name) {
    mat <- matrix(rep(values, length(coding_levels)),
                  nrow = length(values), ncol = length(coding_levels))
    colnames(mat) <- paste0(metric_name, "_", coding_levels)
    as.data.frame(mat, stringsAsFactors = FALSE)
  }

  # The precursor family is all NA -- DIA-NN does not report it -- but the
  # column has to exist: protein_list_reactable() and quant_list_widget()
  # detect their sample columns by grepping for it and abort without it.
  # PG.Coverage and PG.Cscore.RunWise are left out instead of being filled with
  # NA: nothing requires them, every column builder skips a family it does not
  # find, and an empty filter box that blanks the table is worse than a column
  # that is simply not there.
  id_prec <- mk_wide(g_prec, "PG.NrOfPrecursorsIdentified")
  id_pept <- mk_wide(g_pept, "PG.NrOfStrippedSequencesIdentified")

  protein_ID <- cbind(base_info, diann_extras, id_prec, id_pept)

  static_cols  <- c("PG.ProteinGroups", "PG.ProteinDescriptions",
                    "PG.Genes", "PG.MolecularWeight")
  extras_cols  <- names(diann_extras)
  metric_order <- c("PG.NrOfPrecursorsIdentified",
                    "PG.NrOfStrippedSequencesIdentified")
  runwise_order <- unlist(lapply(metric_order,
                                 function(m) paste0(m, "_", coding_levels)))
  final_cols <- c(static_cols, extras_cols,
                  intersect(runwise_order, names(protein_ID)))
  protein_ID <- protein_ID[, final_cols, drop = FALSE]
  protein_ID <- protein_ID[order(protein_ID$PG.ProteinGroups), , drop = FALSE]
  rownames(protein_ID) <- NULL

  if (anyDuplicated(protein_ID$PG.ProteinGroups) > 0) {
    stop("Integrity check failed: protein_ID contains duplicated ",
         "Protein.Group values. Check the input data.")
  }

  # ==========================================================================
  # protein_QUANT
  # ==========================================================================
  if (verbose) message("Processing protein_QUANT...")

  global_metrics <- data.frame(
    PG.NrOfPrecursorsIdentified.Global        = g_prec,
    PG.NrOfStrippedSequencesIdentified.Global = g_pept,
    stringsAsFactors = FALSE
  )

  q_prec <- mk_wide(g_prec, "PG.NrOfPrecursorsUsedForQuantification")
  q_uniq <- mk_wide(g_uniq, "PG.NrOfStrippedSequencesUsedForQuantification")

  protein_QUANT <- cbind(base_info, global_metrics, q_prec, q_uniq, quant)

  static_cols2 <- static_cols
  global_cols  <- names(global_metrics)
  q_metric_order <- c("PG.NrOfPrecursorsUsedForQuantification",
                      "PG.NrOfStrippedSequencesUsedForQuantification",
                      "PG.Quantity")
  q_runwise <- unlist(lapply(q_metric_order,
                             function(m) paste0(m, "_", coding_levels)))
  final_cols2 <- c(static_cols2, global_cols,
                   intersect(q_runwise, names(protein_QUANT)))
  protein_QUANT <- protein_QUANT[, final_cols2, drop = FALSE]
  # Same row order as protein_ID, and for a reason: write_nadia() rebuilds the
  # static columns of protein_quant from the protein table it took out of
  # protein_id, so a different order would not survive the round trip.
  protein_QUANT <- protein_QUANT[order(protein_QUANT$PG.ProteinGroups), ,
                                 drop = FALSE]
  rownames(protein_QUANT) <- NULL

  if (anyDuplicated(protein_QUANT$PG.ProteinGroups) > 0) {
    stop("Integrity check failed: protein_QUANT contains duplicated ",
         "Protein.Group values. Check the input data.")
  }

  # ==========================================================================
  # Optional export
  # ==========================================================================
  if (!is.null(export_dir)) {
    if (verbose) message("Exporting files to: ", export_dir)

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

    if (verbose) message("Files exported successfully.")
  }

  # ==========================================================================
  # Build the result
  # ==========================================================================
  if (verbose) {
    message(
      "Processing complete:\n",
      "  - Runs: ", nrow(run_summary), "\n",
      "  - Proteins (ID): ", nrow(protein_ID), "\n",
      "  - Proteins (QUANT): ", nrow(protein_QUANT)
    )
  }

  result <- list(
    metadata      = as.data.frame(run_summary),
    protein_id    = as.data.frame(protein_ID),
    protein_quant = as.data.frame(protein_QUANT)
  )

  class(result) <- c("diann_data", "proteomics_data", "list")
  result <- .nadia_stamp(result, match.call(), c(report = file_path))
  return(result)
}

# =============================================================================
# Methods for the diann_data class
# =============================================================================

#' Print a summary of a DIA-NN preprocessing result
#'
#' Reports the number of runs, protein groups and conditions, so that the shape
#' of the experiment can be checked at a glance before processing it.
#'
#' @param x A `diann_data` object, as returned by [preprocess_diann()].
#' @param ... Ignored, present for compatibility with the `print` generic.
#' @return `x`, invisibly. Called for the summary it prints.
#'
#' @examples
#' diann <- preprocess_diann(
#'   system.file("extdata", "nadia_diann_report.tsv.gz", package = "NADIA"),
#'   condition_order = c("A", "B", "D"),
#'   verbose = FALSE)
#' print(diann)
#'
#' @export
print.diann_data <- function(x, ...) {
  cat("Preprocessed DIA-NN data\n")
  cat("------------------------\n")
  cat("Runs (metadata):", nrow(x$metadata), "\n")
  cat("Proteins (ID):", nrow(x$protein_id), "\n")
  cat("Proteins (QUANT):", nrow(x$protein_quant), "\n")
  cat("\nConditions:",
      paste(levels(x$metadata$R.Condition) %||%
              unique(x$metadata$R.Condition), collapse = ", "),
      "\n")
  invisible(x)
}
