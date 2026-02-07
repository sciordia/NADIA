# =============================================================================
# Normalization Module
# =============================================================================
#
# Functions for proteomics data normalization:
#   - Zero-to-NA conversion
#   - Protein filtering by group presence
#   - SummarizedExperiment creation
#   - Cyclic Loess normalization (limma)
#
# Dependencies:
#   - SummarizedExperiment, S4Vectors
#   - limma
#
# Author: Sergio Ciordia
# License: MIT
# =============================================================================

# --- Null coalescing operator ---
if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
}

# =============================================================================
# INTERNAL FUNCTIONS
# =============================================================================

#' Convert zero values to NA
#'
#' Replaces 0 values in a matrix or data frame with NA.
#' Useful for DIA-NN data where 0 indicates non-detection.
#'
#' @param data Numeric matrix or data frame
#' @return Object of same type with 0 converted to NA
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

#' Filter proteins by group presence
#'
#' Keeps proteins that have enough non-NA values in at least
#' a minimum number of experimental groups.
#'
#' @param data Intensity matrix (proteins x samples)
#' @param metadata Data frame with sample information
#' @param min_reps Minimum replicates with non-NA values per group.
#'   If NULL, uses half of the smallest group size.
#' @param min_groups Minimum groups meeting min_reps (default: 1)
#' @param grouping_column Name of grouping column in metadata (default: "Condition")
#' @return List with:
#'   - data: Filtered matrix
#'   - keep: Logical vector of kept rows
#'   - summary: Filtering summary
#' @keywords internal
.filter_proteins_by_group <- function(
    data,
    metadata,
    min_reps = NULL,
    min_groups = 1,
    grouping_column = "Condition"
) {
  # Validations
  stopifnot(is.matrix(data) || is.data.frame(data))
  data <- as.matrix(data)

  if (!grouping_column %in% names(metadata)) {
    stop("La columna '", grouping_column, "' no existe en metadata")
  }

  # Align samples
  if (!is.null(rownames(metadata))) {
    common_samples <- intersect(colnames(data), rownames(metadata))
    if (length(common_samples) == 0) {
      stop("No hay muestras en comun entre data y metadata")
    }
    data <- data[, common_samples, drop = FALSE]
    metadata <- metadata[common_samples, , drop = FALSE]
  }

  groups <- as.factor(metadata[[grouping_column]])
  group_levels <- levels(groups)

  # Auto-compute min_reps if not specified
  if (is.null(min_reps)) {
    group_sizes <- table(groups)
    min_reps <- max(1, floor(min(group_sizes) / 2))
  }

  # Compute non-NA counts per group for each protein
  n_present_per_group <- sapply(group_levels, function(g) {
    cols <- which(groups == g)
    if (length(cols) == 0) return(rep(0, nrow(data)))
    rowSums(!is.na(data[, cols, drop = FALSE]))
  })

  if (!is.matrix(n_present_per_group)) {
    n_present_per_group <- matrix(n_present_per_group, ncol = 1)
  }

  # Count groups meeting criteria
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

#' Create SummarizedExperiment from proteomics data
#'
#' @param data Data frame with proteins and intensity values
#' @param metadata Data frame with sample information
#' @param protein_column Protein ID column name
#' @param gene_column Gene name column name
#' @param condition_column Condition column name in metadata
#' @param label_column Sample label column name in metadata
#' @return SummarizedExperiment with assays: raw, log2
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

  # Validate required columns
  if (!protein_column %in% names(data)) {
    stop("Columna '", protein_column, "' no encontrada en data")
  }
  if (!label_column %in% names(metadata)) {
    stop("Columna '", label_column, "' no encontrada en metadata")
  }
  if (!condition_column %in% names(metadata)) {
    stop("Columna '", condition_column, "' no encontrada en metadata")
  }

  # Identify intensity columns
  annotation_cols <- c(protein_column, gene_column, "UniqPepts")
  annotation_cols <- intersect(annotation_cols, names(data))
  intensity_cols <- setdiff(names(data), annotation_cols)

  # Extract matrices
  annotation <- data[, annotation_cols, drop = FALSE]
  intensity <- as.matrix(data[, intensity_cols, drop = FALSE])
  storage.mode(intensity) <- "double"

  # Unique IDs
  protein_ids <- as.character(annotation[[protein_column]])
  protein_ids <- make.unique(protein_ids)
  rownames(intensity) <- protein_ids
  rownames(annotation) <- protein_ids

  # Align metadata with intensity columns
  rownames(metadata) <- metadata[[label_column]]
  common_samples <- intersect(colnames(intensity), rownames(metadata))

  if (length(common_samples) == 0) {
    stop("No hay muestras en comun entre data y metadata")
  }

  intensity <- intensity[, common_samples, drop = FALSE]
  metadata <- metadata[common_samples, , drop = FALSE]

  # Create rowData
  row_data <- S4Vectors::DataFrame(annotation)
  names(row_data)[names(row_data) == protein_column] <- "Protein.IDs"
  if (gene_column %in% names(row_data)) {
    names(row_data)[names(row_data) == gene_column] <- "Gene.Names"
  }
  row_data$IDs <- rownames(row_data)

  # Create colData
  col_data <- S4Vectors::DataFrame(metadata)

  # Create assays: raw and log2
  raw_assay <- intensity
  log2_assay <- log2(intensity)

  # Create SummarizedExperiment
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

#' NA overview for a SummarizedExperiment
#'
#' @param se SummarizedExperiment
#' @param assay_name Assay name to evaluate (default: "log2")
#' @return Data frame with NA statistics per sample
#' @keywords internal
.get_NA_overview <- function(se, assay_name = "log2") {
  stopifnot(inherits(se, "SummarizedExperiment"))

  if (!assay_name %in% SummarizedExperiment::assayNames(se)) {
    stop("Assay '", assay_name, "' no encontrado")
  }

  x <- SummarizedExperiment::assay(se, assay_name)

  # Per sample
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

# =============================================================================
# MAIN FUNCTION
# =============================================================================

#' Normalize proteomics data
#'
#' Complete normalization pipeline: zero-to-NA conversion, protein filtering
#' by group presence, SummarizedExperiment creation, and Cyclic Loess
#' normalization.
#'
#' @param data Data frame with ProteinGroups, GeneNames, UniqPepts + intensity columns
#' @param metadata Data frame with Column, Condition, Replicate
#' @param min_reps Minimum replicates for filtering. NULL = auto: floor(min_group_size / 2)
#' @param min_groups Minimum groups meeting min_reps (default: 1)
#' @param cyclic_loess_method Cyclic Loess method: "fast" or "pairs" (default: "fast")
#' @param cyclic_loess_iterations Number of iterations for Cyclic Loess (default: 3)
#' @param cyclic_loess_span Span parameter for Cyclic Loess (default: 0.7)
#' @param verbose Print progress messages (default: TRUE)
#'
#' @return List with:
#'   \itemize{
#'     \item se: SummarizedExperiment with assays raw, log2, normalized
#'     \item filter_summary: Filtering summary
#'     \item na_overview: NA statistics
#'   }
#'
#' @examples
#' \dontrun{
#' norm_result <- normalize_proteomics(
#'   data = protein_data,
#'   metadata = metadata,
#'   cyclic_loess_method = "fast",
#'   cyclic_loess_iterations = 3
#' )
#' }
#'
#' @export
normalize_proteomics <- function(
    data,
    metadata,
    min_reps = NULL,
    min_groups = 1,
    cyclic_loess_method = c("fast", "pairs"),
    cyclic_loess_iterations = 3,
    cyclic_loess_span = 0.7,
    verbose = TRUE
) {
  # Match method argument
  cyclic_loess_method <- match.arg(cyclic_loess_method)

  # Validate required packages
  if (!requireNamespace("SummarizedExperiment", quietly = TRUE)) {
    stop("Se requiere el paquete 'SummarizedExperiment'. ",
         "Instalalo con BiocManager::install('SummarizedExperiment')")
  }
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Se requiere el paquete 'limma'. ",
         "Instalalo con BiocManager::install('limma')")
  }

  # =========================================================================
  # 1. ZERO TO NA CONVERSION
  # =========================================================================

  if (verbose) cat("\n=== CONVIRTIENDO CEROS A NA ===\n")

  annotation_cols <- c("ProteinGroups", "GeneNames", "UniqPepts")
  intensity_cols <- setdiff(names(data), annotation_cols)

  intensity_mat <- as.matrix(data[, intensity_cols])
  n_zeros <- sum(intensity_mat == 0, na.rm = TRUE)
  intensity_mat <- .zero_to_missing(intensity_mat)
  rownames(intensity_mat) <- data$ProteinGroups

  if (verbose) cat("- Ceros convertidos a NA:", n_zeros, "\n")

  # =========================================================================
  # 2. FILTER PROTEINS BY GROUP PRESENCE
  # =========================================================================

  if (verbose) cat("\n=== FILTRANDO PROTEINAS POR PRESENCIA ===\n")

  filtered <- .filter_proteins_by_group(
    data = intensity_mat,
    metadata = metadata,
    min_reps = min_reps,
    min_groups = min_groups,
    grouping_column = "Condition"
  )

  if (verbose) {
    cat("- Proteinas antes:", filtered$summary$n_total, "\n")
    cat("- Proteinas despues:", filtered$summary$n_keep, "\n")
    cat("- Proteinas eliminadas:", filtered$summary$n_drop, "\n")
    cat("- Min replicas:", filtered$summary$min_reps, "\n")
    cat("- Min grupos:", filtered$summary$min_groups, "\n")
  }

  # Filter protein_data
  protein_data_filtered <- data[filtered$keep, , drop = FALSE]

  # =========================================================================
  # 3. CREATE SUMMARIZEDEXPERIMENT
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

  na_overview <- NULL
  if (verbose) {
    na_overview <- .get_NA_overview(se, "log2")
    global_na <- attr(na_overview, "global")
    cat("- NA global:", global_na$NA.Percentage, "%\n")
  }

  # =========================================================================
  # 4. CYCLIC LOESS NORMALIZATION
  # =========================================================================

  if (verbose) cat("\n=== NORMALIZANDO (Cyclic Loess: method=",
                   cyclic_loess_method, ", iterations=",
                   cyclic_loess_iterations, ", span=",
                   cyclic_loess_span, ") ===\n", sep = "")

  x_log2 <- SummarizedExperiment::assay(se, "log2")

  x_norm <- limma::normalizeCyclicLoess(
    x_log2,
    method = cyclic_loess_method,
    iterations = cyclic_loess_iterations,
    span = cyclic_loess_span
  )
  rownames(x_norm) <- rownames(x_log2)

  # Add normalized assay to SE
  SummarizedExperiment::assay(se, "normalized") <- x_norm

  if (verbose) cat("- Assays disponibles:",
                   paste(SummarizedExperiment::assayNames(se), collapse = ", "), "\n")

  # Compute NA overview if not done yet
  if (is.null(na_overview)) {
    na_overview <- .get_NA_overview(se, "log2")
  }

  # =========================================================================
  # RETURN
  # =========================================================================

  list(
    se = se,
    filter_summary = filtered$summary,
    na_overview = na_overview
  )
}
