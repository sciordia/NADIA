# =============================================================================
# Reading a .nadia file back into objects the rest of the package understands
# =============================================================================
#
# `read_nadia()` loads the file; `nadia_preprocessing()` and `nadia_result()`
# turn it back into the objects the pipeline produced, so that every plotting
# and table function in NADIA works from the file with no adapter in between.
#
# The tables themselves come out of the SQL views stored in the file, not from
# joins written here. This code only restores what SQL cannot carry: R classes,
# factor levels and row names. Keeping the reconstruction in the views means a
# client that is not R -- a browser running DuckDB-WASM, say -- gets the same
# tables from the same definitions.
#
# Author: Sergio Ciordia
# License: GPL-3
# =============================================================================


#' Restore the R classes a SQL table cannot carry
#'
#' DuckDB has no ordered factor and cannot distinguish an integer column from a
#' double one holding whole numbers. `nadia_schema` recorded both, so the
#' round trip can be exact rather than approximate.
#'
#' @param df     Data frame as read from the database.
#' @param schema The `nadia_schema` table.
#' @param table  Which canonical table this is.
#' @return `df` with its original classes and column order.
#' @keywords internal
#' @noRd
.db_restore_types <- function(df, schema, table) {
  s <- schema[schema$table_name == table, , drop = FALSE]
  if (nrow(s) == 0) return(df)
  s <- s[order(s$position), , drop = FALSE]

  for (i in seq_len(nrow(s))) {
    cl <- s$column_name[i]
    if (!cl %in% names(df)) next
    x <- df[[cl]]

    df[[cl]] <- switch(
      s$r_class[i],
      factor    = ,
      ordered   = factor(as.character(x),
                         levels = jsonlite::fromJSON(s$lvls[i]),
                         ordered = isTRUE(s$ordered[i])),
      integer   = as.integer(x),
      numeric   = as.numeric(x),
      character = as.character(x),
      logical   = as.logical(x),
      x
    )
  }
  df[, intersect(s$column_name, names(df)), drop = FALSE]
}

#' Read one of the stored views, with its classes restored
#'
#' @param con    A DBI connection.
#' @param view   View name.
#' @param schema The `nadia_schema` table.
#' @param table  Canonical table the view reproduces.
#' @return A data frame, or `NULL` if the view is not in the file.
#' @keywords internal
#' @noRd
.db_read_view <- function(con, view, schema, table) {
  ok <- DBI::dbGetQuery(con, paste0(
    "SELECT count(*) AS n FROM duckdb_views() WHERE NOT internal AND view_name = ",
    .db_lit(view)))$n > 0
  if (!ok) return(NULL)

  df <- DBI::dbGetQuery(con, paste0("SELECT * FROM ", view))
  df <- .db_restore_types(df, schema, table)
  rownames(df) <- NULL
  df
}

#' Rebuild a matrix from the long intensities table
#'
#' @param long    Data frame with `protein_id`, `sample`, `value`.
#' @param rows    Protein order.
#' @param cols    Sample order.
#' @return A numeric matrix.
#' @keywords internal
#' @noRd
.db_long_to_matrix <- function(long, rows, cols) {
  mat <- matrix(NA_real_, nrow = length(rows), ncol = length(cols),
                dimnames = list(rows, cols))
  idx <- cbind(match(long$protein_id, rows), match(long$sample, cols))
  keep <- !is.na(idx[, 1]) & !is.na(idx[, 2])
  mat[idx[keep, , drop = FALSE]] <- long$value[keep]
  mat
}


# =============================================================================
# Opening a file
# =============================================================================

#' Open a connection to a `.nadia` file
#'
#' For querying the file directly with SQL or dplyr. Most users want
#' [read_nadia()] instead; this is the way in when the file is larger than
#' memory or when a specific query beats loading everything.
#'
#' The caller owns the connection and must close it with
#' `DBI::dbDisconnect(con, shutdown = TRUE)`.
#'
#' @param file      Path to a `.nadia` file.
#' @param read_only Open without write access. `TRUE` by default, which is also
#'   what allows several processes to read the same file at once.
#'
#' @return A DBI connection.
#'
#' @seealso [read_nadia()], [nadia_tables()]
#'
#' @examples
#' if (requireNamespace("duckdb", quietly = TRUE)) {
#'     data(nadia_dia)
#'     f <- file.path(tempdir(), "query.nadia")
#'     write_nadia(f, nadia_dia,
#'                 process_proteomics(nadia_dia, verbose = FALSE),
#'                 verbose = FALSE)
#'
#'     con <- nadia_connect(f)
#'     DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM proteins")
#'     DBI::dbDisconnect(con, shutdown = TRUE)
#'
#'     unlink(f)
#' }
#' @export
nadia_connect <- function(file, read_only = TRUE) {
  .db_require()
  if (!file.exists(file)) {
    stop("'", file, "' does not exist.", call. = FALSE)
  }
  con <- .db_open(file, read_only = read_only)
  tryCatch(.db_check_version(con), error = function(e) {
    DBI::dbDisconnect(con, shutdown = TRUE)
    stop(e)
  })
  con
}

#' List the contents of a `.nadia` file without loading it
#'
#' @param file Path to a `.nadia` file.
#'
#' @return A data frame with `name`, `type` (`table` or `view`) and `rows`.
#'
#' @seealso [read_nadia()]
#'
#' @examples
#' if (requireNamespace("duckdb", quietly = TRUE)) {
#'     data(nadia_dia)
#'     f <- file.path(tempdir(), "tables.nadia")
#'     write_nadia(f, nadia_dia, verbose = FALSE)
#'
#'     nadia_tables(f)
#'
#'     unlink(f)
#' }
#' @export
nadia_tables <- function(file) {
  con <- nadia_connect(file)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  # `internal` filters out the system catalogue, which duckdb_views() lists
  # alongside ours -- some fifty pg_* and information_schema entries.
  tabs <- DBI::dbGetQuery(con,
    "SELECT table_name AS name, estimated_size AS rows FROM duckdb_tables()
     WHERE NOT internal ORDER BY table_name")
  tabs$type <- "table"
  views <- DBI::dbGetQuery(con,
    "SELECT view_name AS name FROM duckdb_views()
     WHERE NOT internal ORDER BY view_name")
  out <- tabs[, c("name", "type", "rows")]

  if (nrow(views) > 0) {
    vrows <- vapply(views$name, function(v) {
      as.numeric(DBI::dbGetQuery(con,
        paste0("SELECT count(*) AS n FROM ", v))$n)
    }, numeric(1))
    out <- rbind(out, data.frame(name = views$name, type = "view",
                                 rows = vrows, stringsAsFactors = FALSE))
  }
  out <- out[order(out$type, out$name), ]
  rownames(out) <- NULL
  out
}

#' Read a `.nadia` file
#'
#' Loads the whole analysis: the preprocessing tables, the processing results,
#' the intensity matrices, the parameters and the provenance. The canonical
#' tables come out of the SQL views stored in the file, with their R classes
#' and column order restored.
#'
#' To get objects the rest of the package can use directly, pass the result to
#' [nadia_preprocessing()] or [nadia_result()].
#'
#' @param file Path to a `.nadia` file.
#'
#' @return An object of class `nadia_db`: a list with `file`, `meta`,
#'   `parameters`, `calls`, `packages`, `source_files`, the canonical tables
#'   (`metadata`, `protein_id`, `protein_quant` and, when present,
#'   `DEPs_results`, `BoxPlot_Input`, `PCA_Input`, `pattern_profiler`), the
#'   `matrices` and the `imputed_values` mask.
#'
#' @seealso [write_nadia()], [nadia_result()], [nadia_preprocessing()]
#'
#' @examples
#' if (requireNamespace("duckdb", quietly = TRUE)) {
#'     data(nadia_dia)
#'     res <- process_proteomics(nadia_dia, verbose = FALSE)
#'     f <- file.path(tempdir(), "read.nadia")
#'     write_nadia(f, nadia_dia, res, verbose = FALSE)
#'
#'     db <- read_nadia(f)
#'     db
#'     identical(db$protein_quant, nadia_dia$protein_quant)
#'
#'     unlink(f)
#' }
#' @export
read_nadia <- function(file) {
  con <- nadia_connect(file)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  meta   <- .db_meta(con)
  schema <- DBI::dbGetQuery(con, "SELECT * FROM nadia_schema")

  db <- list(
    file         = normalizePath(file),
    meta         = meta,
    parameters   = DBI::dbGetQuery(con, "SELECT * FROM parameters"),
    calls        = DBI::dbGetQuery(con, "SELECT * FROM nadia_calls"),
    packages     = DBI::dbGetQuery(con, "SELECT * FROM nadia_packages"),
    source_files = DBI::dbGetQuery(con, "SELECT * FROM nadia_source_files"),
    schema       = schema
  )

  db$metadata      <- .db_read_view(con, "v_metadata",      schema, "metadata")
  db$protein_id    <- .db_read_view(con, "v_protein_id",    schema, "protein_id")
  db$protein_quant <- .db_read_view(con, "v_protein_quant", schema, "protein_quant")

  # `metadata` is indexed by sample code everywhere in the package.
  if (!is.null(db$metadata)) {
    rownames(db$metadata) <- as.character(db$metadata$Coding)
  }

  db$DEPs_results  <- .db_read_view(con, "v_deps_results",  schema, "DEPs_results")
  db$BoxPlot_Input <- .db_read_view(con, "v_boxplot_input", schema, "BoxPlot_Input")
  db$PCA_Input     <- .db_read_view(con, "v_pca_input",     schema, "PCA_Input")
  db$pattern_profiler <- .db_read_view(con, "v_pattern_profiler", schema,
                                       "pattern_profiler")

  # The raw dimension table too, not just the `metadata` view: covariates
  # passed through `covariate_df` reach the analysis via colData and appear
  # here, but not in the preprocessing metadata the view reproduces.
  db$samples <- DBI::dbGetQuery(con, "SELECT * FROM samples ORDER BY position")

  # The matrices, in the order the dimension tables define.
  smp   <- db$samples$sample
  prot  <- DBI::dbGetQuery(con,
    "SELECT protein_id, in_analysis FROM proteins ORDER BY position")

  db$matrices <- list()
  assays <- DBI::dbGetQuery(con, "SELECT DISTINCT assay AS a FROM intensities")$a
  for (a in assays) {
    long <- DBI::dbGetQuery(con, paste0(
      "SELECT protein_id, sample, value FROM intensities WHERE assay = ",
      .db_lit(a)))
    rows <- if (identical(a, "raw")) {
      prot$protein_id
    } else {
      prot$protein_id[!is.na(prot$in_analysis) & prot$in_analysis]
    }
    db$matrices[[a]] <- .db_long_to_matrix(long, rows, smp)
  }

  # `log2` is never stored, being exactly log2 of the raw assay with zeros read
  # as missing. It is rebuilt here rather than in nadia_result() because when
  # `norm_method = "log2"` there is no separate normalized assay and the
  # imputation started from this one: the delta below has nothing to sit on
  # otherwise.
  in_analysis <- prot$protein_id[!is.na(prot$in_analysis) & prot$in_analysis]
  if (!is.null(db$matrices[["raw"]]) && length(in_analysis) > 0) {
    r <- db$matrices[["raw"]][in_analysis, , drop = FALSE]
    r[r == 0] <- NA_real_
    l <- log2(r)
    l[is.infinite(l)] <- NA
    db$matrices[["log2"]] <- l
  }

  # The imputed assay: the normalized one with the filled cells substituted in.
  impn  <- meta[["imputed_assay"]]
  store <- meta[["imputation_storage"]]
  if (!is.na(impn) && identical(store, "delta") &&
      DBI::dbExistsTable(con, "imputed_values")) {
    db$imputed_values <- DBI::dbGetQuery(con, "SELECT * FROM imputed_values")
    from <- if ("BERT" %in% names(db$matrices)) "BERT" else meta[["norm_assay"]]
    mat  <- db$matrices[[from]]
    if (!is.null(mat) && nrow(db$imputed_values) > 0) {
      idx <- cbind(match(db$imputed_values$protein_id, rownames(mat)),
                   match(db$imputed_values$sample,     colnames(mat)))
      mat[idx] <- db$imputed_values$value
    }
    db$matrices[[impn]] <- mat
  }

  class(db) <- c("nadia_db", "list")
  db
}

#' Print a summary of a `.nadia` file
#'
#' @param x   A `nadia_db` object.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @export
print.nadia_db <- function(x, ...) {
  m <- x$meta
  cat("NADIA analysis file\n")
  cat("  File:        ", basename(x$file), "\n", sep = "")
  cat("  Written:     ", m[["created_utc"]], " UTC by NADIA ",
      m[["nadia_version"]], "\n", sep = "")
  cat("  Schema:      version ", m[["schema_version"]],
      " (DuckDB >= ", sub("^v", "", m[["storage_compat"]]), ")\n", sep = "")
  cat("  Experiment:  ", toupper(m[["experiment_type"]]), ", ",
      m[["n_samples"]], " samples, ", m[["n_proteins"]], " proteins\n", sep = "")

  if (identical(m[["has_processing"]], "TRUE")) {
    cat("  Processing:  ", m[["norm_assay"]], " -> ", m[["imputed_assay"]],
        " (", m[["imputation_storage"]], ")\n", sep = "")
    if (!is.null(x$DEPs_results)) {
      cat("  Comparisons: ",
          paste(unique(x$DEPs_results$Comparison), collapse = ", "), "\n", sep = "")
    }
  } else {
    cat("  Processing:  not stored\n")
  }
  if (identical(m[["has_pattern_profiler"]], "TRUE")) {
    cat("  Clusters:    ", length(unique(x$pattern_profiler$Cluster)), "\n", sep = "")
  }
  cat("\nUse nadia_result() and nadia_preprocessing() to get plottable objects.\n")
  invisible(x)
}


# =============================================================================
# Turning the file back into objects
# =============================================================================

#' Rebuild the preprocessing object from a `.nadia` file
#'
#' Returns an object indistinguishable from what `preprocess_spectronaut()`,
#' `preprocess_tmt()` or `preprocess_lfq()` produced, so it can be fed straight
#' to [process_proteomics()] or to the table functions.
#'
#' @param db A `nadia_db` from [read_nadia()].
#'
#' @return A `proteomics_data` object with `metadata`, `protein_id` and
#'   `protein_quant`.
#'
#' @seealso [read_nadia()], [nadia_result()]
#'
#' @examples
#' if (requireNamespace("duckdb", quietly = TRUE)) {
#'     data(nadia_dia)
#'     f <- file.path(tempdir(), "pre.nadia")
#'     write_nadia(f, nadia_dia, verbose = FALSE)
#'
#'     pre <- nadia_preprocessing(read_nadia(f))
#'     identical(pre$protein_id, nadia_dia$protein_id)
#'
#'     unlink(f)
#' }
#' @export
nadia_preprocessing <- function(db) {
  stopifnot(inherits(db, "nadia_db"))

  out <- list(metadata      = db$metadata,
              protein_id    = db$protein_id,
              protein_quant = db$protein_quant)

  type <- switch(db$meta[["experiment_type"]],
                 dia = "spectronaut_data",
                 tmt = "tmt_data",
                 lfq = "lfq_data",
                 "spectronaut_data")
  class(out) <- c(type, "proteomics_data", "list")

  cl <- db$calls$call[db$calls$step == "preprocess"]
  if (length(cl) > 0) attr(out, "nadia_call") <- str2lang(cl[1])
  if (nrow(db$source_files) > 0) attr(out, "nadia_source") <- db$source_files

  out
}

#' Rebuild the processing result from a `.nadia` file
#'
#' Returns a `proteomics_result` with everything [process_proteomics()] returns,
#' the `SummarizedExperiment` included, so that every plotting function in the
#' package works from the file directly:
#'
#' ```
#' db  <- read_nadia("experiment.nadia")
#' res <- nadia_result(db)
#' volcano_highchart_list(res$DEPs_results)
#' pca_highchart_list(res$PCA_Input)
#' ```
#'
#' @param db A `nadia_db` from [read_nadia()].
#'
#' @return A `proteomics_result` object, or an error if the file holds only the
#'   preprocessing step.
#'
#' @seealso [read_nadia()], [nadia_preprocessing()]
#'
#' @examples
#' if (requireNamespace("duckdb", quietly = TRUE)) {
#'     data(nadia_dia)
#'     res <- process_proteomics(nadia_dia, verbose = FALSE)
#'     f <- file.path(tempdir(), "result.nadia")
#'     write_nadia(f, nadia_dia, res, verbose = FALSE)
#'
#'     back <- nadia_result(read_nadia(f))
#'     identical(back$DEPs_results, res$DEPs_results)
#'     class(volcano_highchart_list(back$DEPs_results)[[1]])
#'
#'     unlink(f)
#' }
#' @export
nadia_result <- function(db) {
  stopifnot(inherits(db, "nadia_db"))
  if (!identical(db$meta[["has_processing"]], "TRUE")) {
    stop("This file holds only the preprocessing step; there is no ",
         "processing result to rebuild. Run process_proteomics() on ",
         "nadia_preprocessing(db).", call. = FALSE)
  }

  norm <- db$meta[["norm_assay"]]
  impn <- db$meta[["imputed_assay"]]

  # --- The SummarizedExperiment ---------------------------------------------
  ref     <- db$matrices[[if (!is.na(impn)) impn else norm]]
  keep    <- rownames(ref)
  samples <- colnames(ref)

  assays <- list()
  raw <- db$matrices[["raw"]][keep, samples, drop = FALSE]
  # A zero in the report means "not quantified", as .zero_to_missing() has it.
  raw[raw == 0] <- NA_real_
  assays$raw <- raw

  log2_assay <- log2(raw)
  log2_assay[is.infinite(log2_assay)] <- NA
  assays$log2 <- log2_assay

  # The pipeline adds assays in a fixed order and `assayNames()` reflects it, so
  # the rebuilt object has to follow the same one rather than whatever order the
  # rows came back in.
  rest <- intersect(c(norm, "BERT", impn),
                    setdiff(names(db$matrices), c("raw", "log2")))
  rest <- c(rest, setdiff(names(db$matrices), c("raw", "log2", rest)))
  for (a in rest) {
    assays[[a]] <- db$matrices[[a]][keep, samples, drop = FALSE]
  }

  md <- db$metadata[match(samples, as.character(db$metadata$Coding)), ,
                    drop = FALSE]
  col_data <- data.frame(Column = as.character(md$Coding),
                         Condition = md$R.Condition,
                         Replicate = md$R.Replicate,
                         stringsAsFactors = FALSE)

  # Covariates come from the `samples` table: they were never part of the
  # preprocessing metadata, so the `metadata` view does not carry them.
  covar <- setdiff(names(db$samples),
                   c("sample", "position", "condition_level",
                     names(db$metadata)))
  if (length(covar) > 0) {
    si <- match(samples, db$samples$sample)
    for (nm in covar) col_data[[nm]] <- db$samples[[nm]][si]
  }
  rownames(col_data) <- col_data$Column

  # UniqPepts is the per-protein maximum of the peptide counts used for
  # quantification, exactly as .prepare_protein_data() computes it.
  pq   <- db$protein_quant
  pept <- grep("^PG\\.NrOfStrippedSequencesUsedForQuantification_", names(pq),
               value = TRUE)
  uniq <- if (length(pept) > 0) {
    v <- suppressWarnings(apply(as.matrix(pq[, pept, drop = FALSE]), 1,
                                max, na.rm = TRUE))
    v[!is.finite(v)] <- 0L
    as.integer(v)
  } else {
    rep(NA_integer_, nrow(pq))
  }
  idx <- match(keep, as.character(pq$PG.ProteinGroups))
  row_data <- S4Vectors::DataFrame(
    Protein.IDs = as.character(pq$PG.ProteinGroups)[idx],
    Gene.Names  = as.character(pq$PG.Genes)[idx],
    UniqPepts   = uniq[idx],
    IDs         = keep,
    row.names   = keep)

  se <- SummarizedExperiment::SummarizedExperiment(
    assays   = assays,
    rowData  = row_data,
    colData  = S4Vectors::DataFrame(col_data),
    metadata = list(condition = "Condition", label = "Column"))

  # --- Parameters ------------------------------------------------------------
  pars <- .db_deserialize_step(db$parameters, "process")
  cl   <- db$calls$call[db$calls$step == "process"]
  if (length(cl) > 0) pars$call <- str2lang(cl[1])

  comparisons <- unique(as.character(db$DEPs_results$Comparison))

  out <- list(se_proc       = se,
              DEPs_results  = db$DEPs_results,
              BoxPlot_Input = db$BoxPlot_Input,
              PCA_Input     = db$PCA_Input,
              comparisons   = comparisons,
              parameters    = pars)
  class(out) <- c("proteomics_result", "list")
  out
}

#' Rebuild the Pattern Profiler table from a `.nadia` file
#'
#' @param db A `nadia_db` from [read_nadia()].
#'
#' @return The `long_output` data frame, ready for
#'   [cluster_profile_highchart_list()], or `NULL` if the file holds no
#'   clustering.
#'
#' @seealso [read_nadia()], [nadia_add_pattern_profiler()]
#'
#' @examples
#' if (requireNamespace("duckdb", quietly = TRUE)) {
#'     data(nadia_dia)
#'     f <- file.path(tempdir(), "nopp.nadia")
#'     write_nadia(f, nadia_dia, verbose = FALSE)
#'
#'     is.null(nadia_pattern_profiler(read_nadia(f)))
#'
#'     unlink(f)
#' }
#' @export
nadia_pattern_profiler <- function(db) {
  stopifnot(inherits(db, "nadia_db"))
  db$pattern_profiler
}

#' Turn the stored parameters of one step back into a list
#'
#' @param pars A `parameters` table.
#' @param step Which step to extract.
#' @return A named list.
#' @keywords internal
#' @noRd
.db_deserialize_step <- function(pars, step) {
  p <- pars[pars$step == step, , drop = FALSE]
  if (nrow(p) == 0) return(list())

  out <- lapply(seq_len(nrow(p)), function(i) {
    v <- p$value[i]
    switch(
      p$type[i],
      null      = NULL,
      json      = tryCatch(jsonlite::fromJSON(v), error = function(e) v),
      call      = tryCatch(str2lang(v), error = function(e) v),
      factor    = tryCatch({
                    f <- jsonlite::fromJSON(v)
                    factor(f$values, levels = f$levels, ordered = isTRUE(f$ordered))
                  }, error = function(e) v),
      logical   = as.logical(v),
      integer   = as.integer(v),
      numeric   = as.numeric(v),
      character = v,
      if (startsWith(p$type[i], "json:")) {
        tryCatch(jsonlite::fromJSON(v), error = function(e) v)
      } else {
        v
      }
    )
  })
  names(out) <- p$name
  out[!vapply(out, is.null, logical(1))]
}


# =============================================================================
# Getting out again
# =============================================================================

#' Export the tables of a `.nadia` file to Parquet
#'
#' A way out to an open, stable format. DuckDB reads its own files back a long
#' way, but Parquet is the safer bet for something that has to be readable in
#' ten years by software nobody has written yet, and it is what a repository
#' will expect alongside a publication.
#'
#' @param db      A `nadia_db` from [read_nadia()], or the path to a file.
#' @param dir     Directory to write to; created if it does not exist.
#' @param views   Export the reconstructed canonical tables (`TRUE`, the
#'   default) or the underlying storage tables (`FALSE`).
#' @param verbose Print progress.
#'
#' @return The paths written, invisibly.
#'
#' @seealso [read_nadia()], [write_nadia()]
#'
#' @examples
#' if (requireNamespace("duckdb", quietly = TRUE) &&
#'     requireNamespace("arrow", quietly = TRUE)) {
#'     data(nadia_dia)
#'     f <- file.path(tempdir(), "export.nadia")
#'     write_nadia(f, nadia_dia, verbose = FALSE)
#'
#'     out <- file.path(tempdir(), "nadia_parquet_export")
#'     nadia_export_parquet(f, out, verbose = FALSE)
#'     list.files(out)
#'
#'     unlink(c(f, out), recursive = TRUE)
#' }
#' @export
nadia_export_parquet <- function(db, dir, views = TRUE, verbose = TRUE) {
  .db_require()
  file <- if (inherits(db, "nadia_db")) db$file else db
  if (!file.exists(file)) stop("'", file, "' does not exist.", call. = FALSE)

  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)

  con <- nadia_connect(file)
  on.exit(DBI::dbDisconnect(con, shutdown = TRUE), add = TRUE)

  names_ <- if (views) {
    c(DBI::dbGetQuery(con,
        "SELECT view_name AS n FROM duckdb_views() WHERE NOT internal")$n,
      c("nadia_meta", "parameters", "nadia_calls", "nadia_packages",
        "nadia_source_files"))
  } else {
    DBI::dbGetQuery(con,
      "SELECT table_name AS n FROM duckdb_tables() WHERE NOT internal")$n
  }

  # DuckDB writes the Parquet itself, so `arrow` is not needed to export.
  paths <- vapply(names_, function(nm) {
    p <- file.path(dir, paste0(nm, ".parquet"))
    DBI::dbExecute(con, paste0("COPY (SELECT * FROM ", nm, ") TO ", .db_lit(p),
                               " (FORMAT PARQUET)"))
    if (verbose) message("- ", basename(p))
    p
  }, character(1))

  invisible(unname(paths))
}
