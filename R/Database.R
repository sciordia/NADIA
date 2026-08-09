# =============================================================================
# The .nadia format: one DuckDB file holding a complete analysis
# =============================================================================
#
# An analysis used to leave behind a dozen loose files -- Metadata_*.tsv,
# Protein_ID_*.tsv, Protein_QUANT_*.tsv, matrix_log2_*.tsv, the three
# visualization inputs and the Pattern Profiler parquet -- which had to be kept
# together by hand and which, taken apart, said nothing about how they were
# produced. A `.nadia` file is a DuckDB database holding all of it, plus the
# parameters and the provenance of the run.
#
# Two principles shape the schema:
#
#   1. Every fact is stored once. Nothing that can be derived is stored: the
#      log2 assay is exactly log2 of the raw one, PCA_Input is the imputed assay
#      joined to the DE table, and the raw assay is already the quantification
#      column of protein_quant. What is genuinely irreducible -- the normalized
#      assay and the imputed cells -- is what the file keeps.
#
#   2. The reconstruction lives in the file, as SQL views, not in this R code.
#      That way there is one implementation instead of one per client, and the
#      file is self-describing: anything that speaks DuckDB sees the canonical
#      tables without knowing anything about NADIA.
#
# Nothing inside a `.nadia` depends on R. No serialized objects, no columns only
# R can interpret. Parameters go in as text and JSON so that any client -- a
# browser running DuckDB-WASM, for instance -- can read the whole thing.
#
# Author: Sergio Ciordia
# License: GPL-3
# =============================================================================


# The version of the schema, not of the package. It goes into `nadia_meta` and
# the reader refuses a file it does not know how to read, rather than silently
# misinterpreting it. Bump it whenever a table or a view changes shape.
.NADIA_SCHEMA_VERSION <- 1L

# Written files declare compatibility with this DuckDB release, so that any
# client from 0.10.2 onwards can open them. It is the current default, fixed
# here explicitly so that a change of default does not quietly strand files
# already in circulation.
.NADIA_STORAGE_COMPAT <- "v0.10.2"

# The per-protein columns that carry no sample dimension and are shared by
# protein_id and protein_quant.
.NADIA_PROTEIN_KEY <- "PG.ProteinGroups"

# Tables whose column order, classes and factor levels are recorded in
# `nadia_schema`. That is type metadata, not data: without it a round trip
# cannot restore an ordered factor or tell an integer from a double.
.NADIA_CANONICAL <- c("metadata", "protein_id", "protein_quant",
                      "DEPs_results", "BoxPlot_Input", "PCA_Input",
                      "pattern_profiler")


# =============================================================================
# Connection handling
# =============================================================================

#' Fail with an actionable message when the database packages are missing
#'
#' @return `TRUE` invisibly.
#' @keywords internal
#' @noRd
.db_require <- function() {
  missing <- c("duckdb", "DBI")[
    !vapply(c("duckdb", "DBI"), requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop("The .nadia format needs ", paste(missing, collapse = " and "),
         ". Install with: install.packages(c(\"duckdb\", \"DBI\"))",
         call. = FALSE)
  }
  invisible(TRUE)
}

#' Open a connection to a .nadia file
#'
#' @param file      Path to the file.
#' @param read_only Open without write access.
#' @return A DBI connection.
#' @keywords internal
#' @noRd
.db_open <- function(file, read_only = FALSE) {
  .db_require()

  # `shared_home = FALSE` keeps DuckDB from creating ~/.duckdb to cache
  # extensions and secrets. NADIA loads no extensions, so the directory would be
  # an empty side effect in the user's home -- and a package should not leave
  # things there. The argument is only passed when the installed duckdb has it.
  args <- list()
  if ("shared_home" %in% names(formals(duckdb::duckdb))) {
    args$shared_home <- FALSE
  }
  drv <- do.call(duckdb::duckdb, args)

  con <- DBI::dbConnect(drv, dbdir = file, read_only = read_only)
  if (!read_only) {
    DBI::dbExecute(con, paste0("SET storage_compatibility_version = '",
                               .NADIA_STORAGE_COMPAT, "'"))
  }
  con
}

#' Close a connection, checkpointing first when it was open for writing
#'
#' @param con       A DBI connection.
#' @param read_only Whether it was opened read-only.
#' @return `NULL` invisibly.
#' @keywords internal
#' @noRd
.db_close <- function(con, read_only = FALSE) {
  if (!read_only) {
    try(DBI::dbExecute(con, "CHECKPOINT"), silent = TRUE)
  }
  DBI::dbDisconnect(con, shutdown = TRUE)
  invisible(NULL)
}

#' Quote an identifier for DuckDB
#'
#' Column names such as `PG.Cscore.RunWise_A_1` must be quoted or the dots read
#' as schema separators.
#'
#' @param x Character vector of identifiers.
#' @return The quoted identifiers.
#' @keywords internal
#' @noRd
.db_id <- function(x) paste0("\"", gsub("\"", "\"\"", x), "\"")

#' Quote a string literal for DuckDB
#'
#' @param x Character vector.
#' @return The quoted literals.
#' @keywords internal
#' @noRd
.db_lit <- function(x) paste0("'", gsub("'", "''", x), "'")


# =============================================================================
# Serializing parameters and schemas
# =============================================================================

#' Turn an arbitrary argument into a text value plus a type tag
#'
#' Parameters range from a single number to a nested list (`method_args`) or an
#' unevaluated call. Everything becomes text, with a tag saying how to read it
#' back, so that a client which is not R can still display it.
#'
#' @param x The value.
#' @return A list with `value` and `type`.
#' @keywords internal
#' @noRd
.db_serialize <- function(x) {
  if (is.null(x)) {
    return(list(value = NA_character_, type = "null"))
  }
  if (is.call(x) || is.name(x) || is.function(x)) {
    return(list(value = paste(deparse(x), collapse = " "), type = "call"))
  }
  if (is.factor(x)) {
    return(list(
      value = as.character(jsonlite::toJSON(
        list(values = as.character(x), levels = levels(x),
             ordered = is.ordered(x)), auto_unbox = TRUE)),
      type = "factor"))
  }
  if (length(x) == 1L && is.atomic(x)) {
    return(list(value = as.character(x), type = class(x)[1]))
  }
  # Everything else -- vectors, lists, data frames -- goes in as JSON. One
  # encoding instead of a home-made separator, and one any client can read.
  list(value = as.character(jsonlite::toJSON(x, auto_unbox = TRUE,
                                             null = "null", na = "null",
                                             digits = NA)),
       type = if (is.atomic(x)) paste0("json:", class(x)[1]) else "json")
}

#' Describe a data frame: column order, classes and factor levels
#'
#' This is what makes a round trip exact. DuckDB has no notion of an ordered
#' factor, and it cannot tell whether a column of whole numbers was an integer
#' or a double; recording the classes costs a few hundred rows and removes the
#' guesswork from the reader.
#'
#' @param df    The data frame.
#' @param table Name to record it under.
#' @return A data frame for the `nadia_schema` table.
#' @keywords internal
#' @noRd
.db_describe <- function(df, table) {
  if (is.null(df) || ncol(df) == 0) return(NULL)
  data.frame(
    table_name  = table,
    position    = seq_along(df),
    column_name = names(df),
    r_class     = vapply(df, function(x) class(x)[1], character(1)),
    lvls        = vapply(df, function(x) {
      if (is.factor(x)) {
        as.character(jsonlite::toJSON(levels(x)))
      } else {
        NA_character_
      }
    }, character(1)),
    ordered = vapply(df, function(x) is.ordered(x), logical(1)),
    stringsAsFactors = FALSE, row.names = NULL
  )
}


# =============================================================================
# Building the tables
# =============================================================================

#' Split column names into the metric and the sample they belong to
#'
#' Per-sample columns are named `<metric>_<Coding>`. Matching against the known
#' sample codes rather than against a regular expression is what makes this work
#' with sample names that themselves contain underscores.
#'
#' @param cols    Column names.
#' @param samples Sample codes.
#' @return A data frame with `column`, `metric` and `sample` for the per-sample
#'   columns only.
#' @keywords internal
#' @noRd
.db_split_sample_cols <- function(cols, samples) {
  out <- lapply(cols, function(cl) {
    hit <- samples[endsWith(cl, paste0("_", samples))]
    if (length(hit) == 0) return(NULL)
    # With nested sample names ("A_1" and "X_A_1") the longest match is right.
    hit <- hit[which.max(nchar(hit))]
    data.frame(column = cl,
               metric = substr(cl, 1, nchar(cl) - nchar(hit) - 1L),
               sample = hit, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, out)
  if (is.null(out)) {
    data.frame(column = character(0), metric = character(0),
               sample = character(0), stringsAsFactors = FALSE)
  } else {
    out
  }
}

#' Build the `samples` table
#'
#' Merges the preprocessing metadata with the `colData` of the processed object,
#' which are two views of the same twelve rows. `position` preserves the column
#' order of the original matrices and `condition_level` the levels of the
#' condition factor, neither of which survives a plain SQL table.
#'
#' @param preprocessing Preprocessing object.
#' @param result        Result of `process_proteomics()`, or `NULL`.
#' @return A data frame.
#' @keywords internal
#' @noRd
.db_samples <- function(preprocessing, result = NULL) {
  md <- preprocessing$metadata
  out <- data.frame(sample = as.character(md$Coding),
                    position = seq_len(nrow(md)),
                    stringsAsFactors = FALSE)

  for (nm in names(md)) {
    x <- md[[nm]]
    out[[nm]] <- if (is.factor(x)) as.character(x) else x
  }
  cond <- md$R.Condition
  out$condition_level <- if (is.factor(cond)) as.integer(cond) else NA_integer_

  # Covariates reach the analysis through colData and exist nowhere else.
  if (!is.null(result) && !is.null(result$se_proc)) {
    cd <- as.data.frame(SummarizedExperiment::colData(result$se_proc))
    extra <- setdiff(names(cd), c("Column", "Condition", "Replicate", names(out)))
    if (length(extra) > 0) {
      idx <- match(out$sample, rownames(cd))
      for (nm in extra) {
        x <- cd[[nm]]
        out[[nm]] <- if (is.factor(x)) as.character(x[idx]) else x[idx]
      }
    }
  }
  out
}

#' Build the `proteins` table
#'
#' Everything that does not depend on the sample, taken from both `protein_id`
#' and `protein_quant` (whose shared columns are identical) and deduplicated.
#'
#' @param preprocessing Preprocessing object.
#' @param result        Result of `process_proteomics()`, or `NULL`.
#' @return A data frame.
#' @keywords internal
#' @noRd
.db_proteins <- function(preprocessing, result = NULL) {
  id  <- preprocessing$protein_id
  qt  <- preprocessing$protein_quant
  smp <- as.character(preprocessing$metadata$Coding)

  glob_id <- setdiff(names(id), .db_split_sample_cols(names(id), smp)$column)
  glob_qt <- setdiff(names(qt), .db_split_sample_cols(names(qt), smp)$column)

  out <- id[, glob_id, drop = FALSE]
  add <- setdiff(glob_qt, glob_id)
  if (length(add) > 0) {
    idx <- match(out[[.NADIA_PROTEIN_KEY]], qt[[.NADIA_PROTEIN_KEY]])
    out <- cbind(out, qt[idx, add, drop = FALSE])
  }

  names(out)[names(out) == .NADIA_PROTEIN_KEY] <- "protein_id"
  out$position <- seq_len(nrow(out))

  # Which proteins survived the group-presence filter of the normalization.
  out$in_analysis <- if (!is.null(result) && !is.null(result$se_proc)) {
    out$protein_id %in% rownames(result$se_proc)
  } else {
    NA
  }
  rownames(out) <- NULL
  out
}

#' Build the `protein_sample_metrics` table
#'
#' Long format. Absorbs every per-sample column of `protein_id` and
#' `protein_quant` except the quantity itself, which is an intensity and belongs
#' in `intensities`.
#'
#' @param preprocessing Preprocessing object.
#' @return A data frame with `protein_id`, `sample`, `metric`, `value`.
#' @keywords internal
#' @noRd
.db_protein_sample_metrics <- function(preprocessing) {
  smp <- as.character(preprocessing$metadata$Coding)
  res <- list()

  for (src in c("protein_id", "protein_quant")) {
    df  <- preprocessing[[src]]
    map <- .db_split_sample_cols(names(df), smp)
    map <- map[map$metric != "PG.Quantity", , drop = FALSE]
    if (nrow(map) == 0) next

    keys <- as.character(df[[.NADIA_PROTEIN_KEY]])
    for (i in seq_len(nrow(map))) {
      res[[length(res) + 1L]] <- data.frame(
        protein_id = keys,
        sample     = map$sample[i],
        metric     = map$metric[i],
        value      = as.numeric(df[[map$column[i]]]),
        stringsAsFactors = FALSE)
    }
  }
  if (length(res) == 0) {
    return(data.frame(protein_id = character(0), sample = character(0),
                      metric = character(0), value = numeric(0),
                      stringsAsFactors = FALSE))
  }
  out <- do.call(rbind, res)
  rownames(out) <- NULL
  out
}

#' Turn a matrix into long form
#'
#' @param mat   A matrix with dimnames.
#' @param assay Name to tag the values with.
#' @return A data frame with `protein_id`, `sample`, `assay`, `value`.
#' @keywords internal
#' @noRd
.db_matrix_long <- function(mat, assay) {
  data.frame(
    protein_id = rep(rownames(mat), times = ncol(mat)),
    sample     = rep(colnames(mat), each  = nrow(mat)),
    assay      = assay,
    value      = as.numeric(mat),
    stringsAsFactors = FALSE)
}

#' Build the `intensities` table
#'
#' Only the assays that cannot be derived. `raw` comes from `protein_quant`
#' rather than from the SummarizedExperiment, because that is the primary source
#' and it still holds the proteins the filter later dropped; `log2` is omitted
#' entirely, being exactly `log2(raw)`.
#'
#' @param preprocessing Preprocessing object.
#' @param result        Result of `process_proteomics()`, or `NULL`.
#' @param assays        Names of the extra assays to store.
#' @return A data frame.
#' @keywords internal
#' @noRd
.db_intensities <- function(preprocessing, result = NULL, assays = character(0)) {
  qt  <- preprocessing$protein_quant
  smp <- as.character(preprocessing$metadata$Coding)
  map <- .db_split_sample_cols(names(qt), smp)
  map <- map[map$metric == "PG.Quantity", , drop = FALSE]

  raw <- as.matrix(qt[, map$column, drop = FALSE])
  colnames(raw) <- map$sample
  rownames(raw) <- as.character(qt[[.NADIA_PROTEIN_KEY]])
  raw <- raw[, smp[smp %in% colnames(raw)], drop = FALSE]

  out <- list(.db_matrix_long(raw, "raw"))
  for (a in assays) {
    out[[length(out) + 1L]] <-
      .db_matrix_long(SummarizedExperiment::assay(result$se_proc, a), a)
  }
  res <- do.call(rbind, out)
  rownames(res) <- NULL
  res
}

#' Build the `imputed_values` table
#'
#' Only the cells that were filled, with the branch that filled them. This is
#' both the imputed data and the MAR/MNAR mask: on a dataset with 8 % missing
#' values it stores 8 % of the cells instead of a second full matrix.
#'
#' The branch comes from the masks that `process_proteomics()` now keeps in the
#' object's metadata. Note that `mnar_mask` also covers observed cells -- it
#' marks whole conditions -- so it has to be intersected with the cells that
#' were actually missing.
#'
#' @param result    Result of `process_proteomics()`.
#' @param norm_name Name of the assay imputation started from.
#' @param imp_name  Name of the imputed assay.
#' @return A data frame with `protein_id`, `sample`, `value`, `branch`.
#' @keywords internal
#' @noRd
.db_imputed_values <- function(result, norm_name, imp_name) {
  se  <- result$se_proc
  nrm <- SummarizedExperiment::assay(se, norm_name)
  imp <- SummarizedExperiment::assay(se, imp_name)

  filled <- is.na(nrm) & !is.na(imp)
  if (!any(filled)) {
    return(data.frame(protein_id = character(0), sample = character(0),
                      value = numeric(0), branch = character(0),
                      stringsAsFactors = FALSE))
  }

  md     <- S4Vectors::metadata(se)
  branch <- matrix("filled", nrow(nrm), ncol(nrm), dimnames = dimnames(nrm))
  if (!is.null(md$mar_mask)) {
    branch[md$mar_mask & filled] <- "MAR"
  }
  if (!is.null(md$mnar_mask)) {
    branch[md$mnar_mask & filled] <- "MNAR"
  }

  idx <- which(filled, arr.ind = TRUE)
  data.frame(
    protein_id = rownames(nrm)[idx[, 1]],
    sample     = colnames(nrm)[idx[, 2]],
    value      = imp[filled],
    branch     = branch[filled],
    stringsAsFactors = FALSE)
}

#' Build the `de_results` table
#'
#' Drops the two columns that are stored elsewhere: `Gene.Names` lives in
#' `proteins` and `Assay` is constant across the table and goes into
#' `nadia_meta`.
#'
#' @param result Result of `process_proteomics()`.
#' @return A data frame.
#' @keywords internal
#' @noRd
.db_de_results <- function(result) {
  de  <- result$DEPs_results
  out <- de[, setdiff(names(de), c("Gene.Names", "Assay")), drop = FALSE]
  names(out)[names(out) == "Protein.IDs"] <- "protein_id"
  out$Change <- as.character(out$Change)
  rownames(out) <- NULL
  out
}

#' Build the two Pattern Profiler tables
#'
#' The z-score profile depends only on the protein, so in `long_output` it is
#' repeated once per cluster the protein belongs to. Splitting membership from
#' profile removes that duplication.
#'
#' @param pp Result of `pattern_profiler_analysis()`.
#' @return A list with `pp_membership` and `pp_profile`.
#' @keywords internal
#' @noRd
.db_pattern_profiler <- function(pp) {
  lo    <- pp$long_output
  conds <- setdiff(names(lo), c("FeatureID", "Cluster", "Membership"))

  mem <- data.frame(protein_id = as.character(lo$FeatureID),
                    Cluster    = as.integer(lo$Cluster),
                    Membership = as.numeric(lo$Membership),
                    position   = seq_len(nrow(lo)),
                    stringsAsFactors = FALSE)

  uniq <- lo[!duplicated(lo$FeatureID), , drop = FALSE]
  prof <- do.call(rbind, lapply(seq_along(conds), function(i) {
    data.frame(protein_id = as.character(uniq$FeatureID),
               Condition  = conds[i],
               position   = i,
               ZScore     = as.numeric(uniq[[conds[i]]]),
               stringsAsFactors = FALSE)
  }))
  rownames(prof) <- NULL
  list(pp_membership = mem, pp_profile = prof)
}

#' Collect the parameters of every step into one long table
#'
#' @param preprocessing Preprocessing object.
#' @param result        Result of `process_proteomics()`, or `NULL`.
#' @param pp            Result of `pattern_profiler_analysis()`, or `NULL`.
#' @return A data frame with `step`, `name`, `value`, `type`.
#' @keywords internal
#' @noRd
.db_parameters <- function(preprocessing, result = NULL, pp = NULL) {

  # Returns rows rather than assigning into an enclosing list, which keeps
  # superassignment out of the package altogether.
  as_rows <- function(step, lst) {
    if (length(lst) == 0) return(NULL)
    do.call(rbind, lapply(names(lst), function(nm) {
      s <- .db_serialize(lst[[nm]])
      data.frame(step = step, name = nm, value = s$value, type = s$type,
                 stringsAsFactors = FALSE)
    }))
  }

  # An argument taken from a call: evaluated where it can be, kept as the
  # expression itself where it cannot.
  from_call <- function(cl, drop = character(0)) {
    if (is.null(cl)) return(list())
    args <- as.list(cl)[-1]
    args <- args[!names(args) %in% drop]
    lapply(args, function(a) {
      if (is.call(a) || is.name(a)) a else tryCatch(eval(a), error = function(e) a)
    })
  }

  rows <- list()

  # The preprocessors record no parameters, so the call is the only source.
  rows$pre <- as_rows("preprocess",
                      from_call(attr(preprocessing, "nadia_call")))

  if (!is.null(result)) {
    rows$proc <- as_rows("process", result$parameters)
    rows$comp <- as_rows("process",
                         list(comparisons = as.character(result$comparisons)))
  }

  if (!is.null(pp)) {
    keep <- c("optimal_c", "m", "conditions", "min_membership",
              "n_features_input", "n_features_final", "n_rows_output")
    rows$pp <- as_rows("pattern_profiler", pp[intersect(keep, names(pp))])
    # `se_proc` and `DEPs_results` are whole objects; the name they were passed
    # under is informative, their contents are already in the file.
    rows$pp_call <- as_rows("pattern_profiler",
                            from_call(pp$call, c("se_proc", "DEPs_results")))
  }

  rows <- rows[!vapply(rows, is.null, logical(1))]

  if (length(rows) == 0) {
    return(data.frame(step = character(0), name = character(0),
                      value = character(0), type = character(0),
                      stringsAsFactors = FALSE))
  }
  out <- do.call(rbind, rows)
  out[!duplicated(out[, c("step", "name")]), , drop = FALSE]
}

#' Which optional packages actually took part in the run
#'
#' Not a `sessionInfo()`: only the packages whose version could change the
#' numbers in the file.
#'
#' @param result Result of `process_proteomics()`, or `NULL`.
#' @return A data frame with `package` and `version`.
#' @keywords internal
#' @noRd
.db_packages <- function(result = NULL) {
  pkgs <- c("NADIA", "duckdb", "DBI")

  if (!is.null(result)) {
    p <- result$parameters
    pkgs <- c(pkgs, if (identical(p$de_method, "limpa")) "limpa" else "limma")
    if (isTRUE(p$batch_correct)) pkgs <- c(pkgs, "BERT")

    by_method <- c(bpca = "pcaMethods", knn = "impute", mice = "mice",
                   missForest = "missForest", Impseq = "rrcovNA",
                   Impseqrob = "rrcovNA", QRILC = "imputeLCMD",
                   MLE = "norm", MinProb = "imputeLCMD", limpa = "limpa")
    used <- unlist(p[c("imp_method", "mar_method", "mnar_method")])
    pkgs <- c(pkgs, unname(by_method[intersect(used, names(by_method))]))
    if (identical(p$norm_method, "vsn")) pkgs <- c(pkgs, "vsn")
    if (identical(p$norm_method, "Rlr")) pkgs <- c(pkgs, "MASS")
  }

  pkgs <- unique(pkgs)
  ver  <- vapply(pkgs, function(p) {
    tryCatch(as.character(utils::packageVersion(p)),
             error = function(e) NA_character_)
  }, character(1))
  ok <- !is.na(ver)
  data.frame(package = pkgs[ok], version = ver[ok],
             stringsAsFactors = FALSE, row.names = NULL)
}


# =============================================================================
# Writing
# =============================================================================

#' Write a complete analysis to a single `.nadia` file
#'
#' Stores the preprocessing tables, the processing results and, optionally, the
#' Pattern Profiler output in one DuckDB database, together with the parameters
#' and the provenance of the run. Nothing is duplicated: everything that can be
#' derived is left out and rebuilt on reading, through SQL views stored in the
#' file itself.
#'
#' Use [read_nadia()] to read it back, and [nadia_result()] to get an object the
#' rest of the package can plot directly.
#'
#' @section What goes in:
#' Twelve tables. `samples` and `proteins` are the dimensions;
#' `protein_sample_metrics`, `intensities`, `imputed_values` and `de_results`
#' the facts; `pp_membership` and `pp_profile` the optional clustering; and
#' `nadia_meta`, `nadia_packages`, `nadia_source_files`, `nadia_calls`,
#' `parameters` and `nadia_schema` the provenance.
#'
#' Three things are deliberately *not* stored, because they are exactly
#' derivable: the `log2` assay (`log2` of the raw one), the `PCA_Input` table
#' (the imputed assay joined to the DE results) and the observed cells of the
#' imputed assay (identical to the normalized ones). The last of these is
#' checked rather than assumed, since a model-based method such as `limpa` can
#' revise observed values; when it does, the full assay is stored instead.
#'
#' @section Where not to put it:
#' Not inside a synchronised folder such as iCloud Drive, Dropbox or OneDrive
#' while you are working. DuckDB keeps the file open with locks and partial
#' writes, and a synchronising client can corrupt it. Write it locally and copy
#' it afterwards.
#'
#' @param file             Path to write to. The `.nadia` extension is a
#'   convention; the file is a DuckDB database whatever it is called.
#' @param preprocessing    Object from `preprocess_spectronaut()`,
#'   `preprocess_tmt()` or `preprocess_lfq()`. Required: it is the base of the
#'   file.
#' @param result           Object from [process_proteomics()], or `NULL` to
#'   store only the preprocessing.
#' @param pattern_profiler Object from [pattern_profiler_analysis()], or `NULL`.
#'   Can also be added later with [nadia_add_pattern_profiler()].
#' @param overwrite        Replace `file` if it already exists.
#' @param verbose          Print progress.
#'
#' @return The path, invisibly.
#'
#' @seealso [read_nadia()], [nadia_result()], [nadia_add_pattern_profiler()]
#'
#' @examples
#' if (requireNamespace("duckdb", quietly = TRUE)) {
#'     data(nadia_dia)
#'     res <- process_proteomics(nadia_dia, verbose = FALSE)
#'
#'     f <- file.path(tempdir(), "example.nadia")
#'     write_nadia(f, nadia_dia, res, verbose = FALSE)
#'
#'     db <- read_nadia(f)
#'     db
#'
#'     unlink(f)
#' }
#' @export
write_nadia <- function(file,
                        preprocessing,
                        result = NULL,
                        pattern_profiler = NULL,
                        overwrite = FALSE,
                        verbose = TRUE) {
  .db_require()

  if (!inherits(preprocessing, "proteomics_data")) {
    stop("`preprocessing` must be the object returned by preprocess_spectronaut(), ",
         "preprocess_tmt() or preprocess_lfq().", call. = FALSE)
  }
  if (!is.null(result) && !inherits(result, "proteomics_result")) {
    stop("`result` must be the object returned by process_proteomics().",
         call. = FALSE)
  }
  if (!is.null(pattern_profiler) &&
      !inherits(pattern_profiler, "pattern_profiler_result")) {
    stop("`pattern_profiler` must be the object returned by ",
         "pattern_profiler_analysis().", call. = FALSE)
  }

  if (file.exists(file)) {
    if (!overwrite) {
      stop("'", file, "' already exists. Pass overwrite = TRUE to replace it.",
           call. = FALSE)
    }
    unlink(c(file, paste0(file, ".wal")), force = TRUE)
  }
  dir <- dirname(file)
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)

  con <- .db_open(file)
  on.exit(.db_close(con), add = TRUE)

  if (verbose) message("=== WRITING ", basename(file), " ===")

  # ---------------------------------------------------------------------------
  # Which assays exist, and how the imputed one is stored
  # ---------------------------------------------------------------------------
  norm_name <- imp_name <- NA_character_
  extra     <- character(0)
  imp_store <- "none"

  if (!is.null(result)) {
    all_assays <- SummarizedExperiment::assayNames(result$se_proc)
    norm_name  <- result$parameters$norm_method
    if (!norm_name %in% all_assays) norm_name <- "log2"

    imp_name <- setdiff(all_assays, c("raw", "log2", norm_name, "BERT"))
    imp_name <- if (length(imp_name) > 0) imp_name[1] else NA_character_

    # The assay imputation actually started from: batch correction, when used,
    # replaces the normalized assay as the input.
    imp_from <- if ("BERT" %in% all_assays) "BERT" else norm_name

    extra <- setdiff(intersect(c(norm_name, "BERT"), all_assays), "log2")

    if (!is.na(imp_name)) {
      nrm <- SummarizedExperiment::assay(result$se_proc, imp_from)
      imp <- SummarizedExperiment::assay(result$se_proc, imp_name)
      obs <- !is.na(nrm)
      # Checked, not assumed: a model-based method may revise observed values,
      # in which case only the whole assay is faithful.
      imp_store <- if (isTRUE(all.equal(nrm[obs], imp[obs]))) "delta" else "full"
      if (identical(imp_store, "full")) extra <- c(extra, imp_name)
    }
  }

  # ---------------------------------------------------------------------------
  # The tables
  # ---------------------------------------------------------------------------
  tables <- list(
    samples               = .db_samples(preprocessing, result),
    proteins              = .db_proteins(preprocessing, result),
    protein_sample_metrics = .db_protein_sample_metrics(preprocessing),
    intensities           = .db_intensities(preprocessing, result, extra),
    parameters            = .db_parameters(preprocessing, result, pattern_profiler),
    nadia_packages        = .db_packages(result)
  )

  if (!is.null(result)) {
    tables$de_results <- .db_de_results(result)
    if (!is.na(imp_name) && identical(imp_store, "delta")) {
      imp_from <- if ("BERT" %in% SummarizedExperiment::assayNames(result$se_proc)) {
        "BERT"
      } else {
        norm_name
      }
      tables$imputed_values <- .db_imputed_values(result, imp_from, imp_name)
    }
  }

  if (!is.null(pattern_profiler)) {
    tables <- c(tables, .db_pattern_profiler(pattern_profiler))
  }

  # Provenance -----------------------------------------------------------------
  src <- attr(preprocessing, "nadia_source")
  tables$nadia_source_files <- if (is.null(src)) {
    data.frame(role = character(0), path = character(0),
               size_bytes = numeric(0), mtime = character(0),
               md5 = character(0), stringsAsFactors = FALSE)
  } else {
    src
  }

  calls <- list()
  cl <- attr(preprocessing, "nadia_call")
  if (!is.null(cl)) calls$preprocess <- paste(deparse(cl), collapse = " ")
  if (!is.null(result$parameters$call)) {
    calls$process <- paste(deparse(result$parameters$call), collapse = " ")
  }
  if (!is.null(pattern_profiler$call)) {
    calls$pattern_profiler <- paste(deparse(pattern_profiler$call), collapse = " ")
  }
  # Built column by column rather than from the list, so that a run with no
  # recorded calls still produces a two-column table with no rows instead of a
  # data frame with no columns at all, which DuckDB rejects.
  tables$nadia_calls <- data.frame(
    step = as.character(names(calls) %||% character(0)),
    call = as.character(unlist(calls, use.names = FALSE) %||% character(0)),
    stringsAsFactors = FALSE, row.names = NULL)

  # The column order, classes and factor levels of the canonical tables.
  schema <- do.call(rbind, list(
    .db_describe(preprocessing$metadata,      "metadata"),
    .db_describe(preprocessing$protein_id,    "protein_id"),
    .db_describe(preprocessing$protein_quant, "protein_quant"),
    .db_describe(result$DEPs_results,         "DEPs_results"),
    .db_describe(result$BoxPlot_Input,        "BoxPlot_Input"),
    .db_describe(result$PCA_Input,            "PCA_Input"),
    .db_describe(pattern_profiler$long_output, "pattern_profiler")))
  tables$nadia_schema <- if (is.null(schema)) {
    data.frame(table_name = character(0), position = integer(0),
               column_name = character(0), r_class = character(0),
               lvls = character(0), ordered = logical(0),
               stringsAsFactors = FALSE)
  } else {
    schema
  }

  type <- sub("_data$", "", class(preprocessing)[1])
  meta <- c(
    schema_version     = as.character(.NADIA_SCHEMA_VERSION),
    nadia_version      = as.character(utils::packageVersion("NADIA")),
    created_utc        = format(Sys.time(), tz = "UTC", "%Y-%m-%d %H:%M:%S"),
    r_version          = paste(R.version$major, R.version$minor, sep = "."),
    platform           = R.version$platform,
    duckdb_version     = as.character(utils::packageVersion("duckdb")),
    storage_compat     = .NADIA_STORAGE_COMPAT,
    experiment_type    = if (identical(type, "spectronaut")) "dia" else type,
    n_samples          = as.character(nrow(tables$samples)),
    n_proteins         = as.character(nrow(tables$proteins)),
    norm_assay         = norm_name,
    imputed_assay      = imp_name,
    imputation_storage = imp_store,
    has_processing     = as.character(!is.null(result)),
    has_pattern_profiler = as.character(!is.null(pattern_profiler)),
    de_assay           = if (!is.null(result)) result$DEPs_results$Assay[1] else NA_character_,
    alpha              = if (!is.null(result)) as.character(result$parameters$alpha) else NA_character_
  )
  tables$nadia_meta <- data.frame(key = names(meta), value = unname(meta),
                                  stringsAsFactors = FALSE, row.names = NULL)

  # ---------------------------------------------------------------------------
  # Write
  # ---------------------------------------------------------------------------
  for (nm in names(tables)) {
    DBI::dbWriteTable(con, nm, as.data.frame(tables[[nm]]), overwrite = TRUE)
    if (verbose) {
      message(sprintf("- %-24s %8d rows", nm, nrow(tables[[nm]])))
    }
  }

  .db_create_views(con)
  if (verbose) message("- views created")

  invisible(file)
}


#' Add a Pattern Profiler run to an existing `.nadia` file
#'
#' The clustering is an optional extra step, so it is written separately: this
#' adds it to a file that already holds the analysis, and refuses to create one
#' from scratch.
#'
#' @param file             Path to an existing `.nadia` file.
#' @param pattern_profiler Object from [pattern_profiler_analysis()].
#' @param overwrite        Replace a clustering already in the file.
#' @param verbose          Print progress.
#'
#' @return The path, invisibly.
#'
#' @seealso [write_nadia()]
#'
#' @examples
#' if (requireNamespace("duckdb", quietly = TRUE)) {
#'     data(nadia_dia)
#'     res <- process_proteomics(nadia_dia, verbose = FALSE)
#'     f <- file.path(tempdir(), "pp.nadia")
#'     write_nadia(f, nadia_dia, res, verbose = FALSE)
#'
#'     # pp <- pattern_profiler_analysis(res$se_proc, res$DEPs_results)
#'     # nadia_add_pattern_profiler(f, pp)
#'
#'     unlink(f)
#' }
#' @export
nadia_add_pattern_profiler <- function(file, pattern_profiler,
                                       overwrite = FALSE, verbose = TRUE) {
  .db_require()

  if (!file.exists(file)) {
    stop("'", file, "' does not exist. A Pattern Profiler run is added to a ",
         "file that already holds the analysis; create it with write_nadia() ",
         "first.", call. = FALSE)
  }
  if (!inherits(pattern_profiler, "pattern_profiler_result")) {
    stop("`pattern_profiler` must be the object returned by ",
         "pattern_profiler_analysis().", call. = FALSE)
  }

  con <- .db_open(file)
  on.exit(.db_close(con), add = TRUE)
  .db_check_version(con)

  if (DBI::dbExistsTable(con, "pp_membership") && !overwrite) {
    stop("'", file, "' already holds a Pattern Profiler run. Pass ",
         "overwrite = TRUE to replace it.", call. = FALSE)
  }

  tabs <- .db_pattern_profiler(pattern_profiler)
  for (nm in names(tabs)) {
    DBI::dbWriteTable(con, nm, tabs[[nm]], overwrite = TRUE)
    if (verbose) message(sprintf("- %-24s %8d rows", nm, nrow(tabs[[nm]])))
  }

  sch <- .db_describe(pattern_profiler$long_output, "pattern_profiler")
  DBI::dbExecute(con, "DELETE FROM nadia_schema WHERE table_name = 'pattern_profiler'")
  DBI::dbAppendTable(con, "nadia_schema", sch)

  pars <- .db_parameters(list(), NULL, pattern_profiler)
  if (nrow(pars) > 0) {
    DBI::dbExecute(con, "DELETE FROM parameters WHERE step = 'pattern_profiler'")
    DBI::dbAppendTable(con, "parameters", pars)
  }
  if (!is.null(pattern_profiler$call)) {
    DBI::dbExecute(con, "DELETE FROM nadia_calls WHERE step = 'pattern_profiler'")
    DBI::dbAppendTable(con, "nadia_calls", data.frame(
      step = "pattern_profiler",
      call = paste(deparse(pattern_profiler$call), collapse = " "),
      stringsAsFactors = FALSE))
  }
  DBI::dbExecute(con,
    "UPDATE nadia_meta SET value = 'TRUE' WHERE key = 'has_pattern_profiler'")

  .db_create_views(con)
  invisible(file)
}


# =============================================================================
# The views: the reconstruction, stored inside the file
# =============================================================================
#
# Every canonical table is rebuilt by a view rather than by R code. That buys
# two things. There is one implementation instead of one per client, so an
# application reading the file with DuckDB-WASM in a browser gets exactly the
# table `pca_highchart_list()` consumes without reimplementing a single pivot.
# And the file explains itself: open it with any DuckDB client and the tables
# are there, whether or not you have ever heard of NADIA.
#
# The SQL is generated from `nadia_schema`, which holds the column order of the
# originals, so the views reproduce that order rather than an arbitrary one.

#' Refuse a file written by a newer schema
#'
#' @param con A DBI connection.
#' @return The schema version, invisibly.
#' @keywords internal
#' @noRd
.db_check_version <- function(con) {
  if (!DBI::dbExistsTable(con, "nadia_meta")) {
    stop("This is not a NADIA file: it has no 'nadia_meta' table.", call. = FALSE)
  }
  v <- DBI::dbGetQuery(con, "SELECT value FROM nadia_meta WHERE key = 'schema_version'")
  v <- suppressWarnings(as.integer(v$value[1]))
  if (is.na(v)) stop("This file declares no schema version.", call. = FALSE)
  if (v > .NADIA_SCHEMA_VERSION) {
    stop("This file uses .nadia schema version ", v, ", and this version of ",
         "NADIA reads up to ", .NADIA_SCHEMA_VERSION,
         ". Upgrade NADIA to open it.", call. = FALSE)
  }
  invisible(v)
}

#' Read the key/value metadata table
#'
#' @param con A DBI connection.
#' @return A named character vector.
#' @keywords internal
#' @noRd
.db_meta <- function(con) {
  m <- DBI::dbGetQuery(con, "SELECT key, value FROM nadia_meta")
  stats::setNames(m$value, m$key)
}

#' Column names of a canonical table, in their original order
#'
#' @param schema The `nadia_schema` table.
#' @param table  Table name.
#' @return A character vector, empty if the table was not stored.
#' @keywords internal
#' @noRd
.db_cols_of <- function(schema, table) {
  s <- schema[schema$table_name == table, , drop = FALSE]
  if (nrow(s) == 0) return(character(0))
  s$column_name[order(s$position)]
}

#' A `MAX(CASE WHEN ...)` pivot expression
#'
#' Written out by hand rather than with DuckDB's `PIVOT`, because the column
#' order has to match the original table exactly and `PIVOT` decides it itself.
#'
#' @param metric Metric to select, or `NULL` for a single-metric source.
#' @param sample Sample to select.
#' @param alias  Name of the resulting column.
#' @return A SQL fragment.
#' @keywords internal
#' @noRd
.db_pivot_expr <- function(metric, sample, alias) {
  cond <- paste0("sample = ", .db_lit(sample))
  if (!is.null(metric)) {
    cond <- paste0(cond, " AND metric = ", .db_lit(metric))
  }
  paste0("MAX(CASE WHEN ", cond, " THEN value END) AS ", .db_id(alias))
}

#' Create every view the file can support
#'
#' Views are created only for the tables whose sources are present, so a file
#' holding just the preprocessing gets `v_metadata`, `v_protein_id` and
#' `v_protein_quant` and nothing else.
#'
#' @param con A DBI connection open for writing.
#' @return `NULL` invisibly.
#' @keywords internal
#' @noRd
.db_create_views <- function(con) {
  meta   <- .db_meta(con)
  schema <- DBI::dbGetQuery(con, "SELECT * FROM nadia_schema")
  smp    <- DBI::dbGetQuery(con, "SELECT * FROM samples ORDER BY position")
  samples <- smp$sample

  drop_view <- function(v) DBI::dbExecute(con, paste0("DROP VIEW IF EXISTS ", v))
  make_view <- function(v, sql) {
    drop_view(v)
    DBI::dbExecute(con, paste0("CREATE VIEW ", v, " AS ", sql))
  }

  # --- v_metadata ------------------------------------------------------------
  md_cols <- .db_cols_of(schema, "metadata")
  if (length(md_cols) > 0) {
    make_view("v_metadata", paste0(
      "SELECT ", paste(.db_id(md_cols), collapse = ", "),
      " FROM samples ORDER BY position"))
  }

  # --- v_protein_id / v_protein_quant ----------------------------------------
  prot_cols <- DBI::dbGetQuery(con, "SELECT * FROM proteins LIMIT 0")

  build_protein_view <- function(view, table, quantity) {
    cols <- .db_cols_of(schema, table)
    if (length(cols) == 0) return(invisible(NULL))

    metrics  <- .db_split_sample_cols(cols, samples)
    metrics  <- metrics[metrics$metric != "PG.Quantity", , drop = FALSE]
    quantity_cols <- .db_split_sample_cols(cols, samples)
    quantity_cols <- quantity_cols[quantity_cols$metric == "PG.Quantity", ,
                                   drop = FALSE]

    sel <- vapply(cols, function(cl) {
      if (identical(cl, .NADIA_PROTEIN_KEY)) {
        paste0("p.protein_id AS ", .db_id(cl))
      } else if (cl %in% metrics$column) {
        paste0("m.", .db_id(cl))
      } else if (cl %in% quantity_cols$column) {
        paste0("q.", .db_id(cl))
      } else {
        paste0("p.", .db_id(cl))
      }
    }, character(1))

    sql <- paste0("SELECT ", paste(sel, collapse = ", "), " FROM proteins p")

    if (nrow(metrics) > 0) {
      piv <- vapply(seq_len(nrow(metrics)), function(i) {
        .db_pivot_expr(metrics$metric[i], metrics$sample[i], metrics$column[i])
      }, character(1))
      sql <- paste0(sql,
        " LEFT JOIN (SELECT protein_id, ", paste(piv, collapse = ", "),
        " FROM protein_sample_metrics GROUP BY protein_id) m",
        " ON m.protein_id = p.protein_id")
    }
    if (quantity && nrow(quantity_cols) > 0) {
      piv <- vapply(seq_len(nrow(quantity_cols)), function(i) {
        .db_pivot_expr(NULL, quantity_cols$sample[i], quantity_cols$column[i])
      }, character(1))
      sql <- paste0(sql,
        " LEFT JOIN (SELECT protein_id, ", paste(piv, collapse = ", "),
        " FROM intensities WHERE assay = 'raw' GROUP BY protein_id) q",
        " ON q.protein_id = p.protein_id")
    }
    make_view(view, paste0(sql, " ORDER BY p.position"))
  }

  build_protein_view("v_protein_id",    "protein_id",    quantity = FALSE)
  build_protein_view("v_protein_quant", "protein_quant", quantity = TRUE)

  # --- Everything below needs the processing step ----------------------------
  if (!identical(meta[["has_processing"]], "TRUE")) return(invisible(NULL))

  norm  <- meta[["norm_assay"]]
  impn  <- meta[["imputed_assay"]]
  store <- meta[["imputation_storage"]]
  alpha <- suppressWarnings(as.numeric(meta[["alpha"]]))
  if (is.na(alpha)) alpha <- 0.05

  # The assay imputation started from; batch correction replaces the normalized
  # one as its input.
  has_bert <- nrow(DBI::dbGetQuery(con,
    "SELECT 1 FROM intensities WHERE assay = 'BERT' LIMIT 1")) > 0
  imp_from <- if (has_bert) "BERT" else norm

  # `raw` keeps the zeros of the original report; every derived assay treats a
  # zero as missing, exactly as .zero_to_missing() does in R.
  sql_log2 <- "log2(NULLIF(i.value, 0))"

  # Where an assay comes from. `log2` is never stored, so it is derived here --
  # which matters when `norm_method = "log2"`, since then imputation started
  # from it and there is no separate normalized assay to build on.
  assay_src <- function(a) {
    if (identical(a, "log2")) {
      paste0("SELECT i.protein_id, i.sample, ", sql_log2, " AS value",
             " FROM intensities i",
             " JOIN proteins p ON p.protein_id = i.protein_id AND p.in_analysis",
             " WHERE i.assay = 'raw'")
    } else {
      paste0("SELECT protein_id, sample, value FROM intensities WHERE assay = ",
             .db_lit(a))
    }
  }

  # The imputed assay: the normalized one with the filled cells substituted in,
  # unless the whole assay had to be stored.
  imputed_src <- if (identical(store, "delta")) {
    paste0(
      "SELECT n.protein_id, n.sample, COALESCE(v.value, n.value) AS value",
      " FROM (", assay_src(imp_from), ") n",
      " LEFT JOIN imputed_values v",
      " ON v.protein_id = n.protein_id AND v.sample = n.sample")
  } else {
    assay_src(impn)
  }

  # --- v_deps_results --------------------------------------------------------
  de_cols <- .db_cols_of(schema, "DEPs_results")
  if (length(de_cols) > 0) {
    de_assay <- meta[["de_assay"]]
    sel <- vapply(de_cols, function(cl) {
      if (identical(cl, "Protein.IDs")) {
        "d.protein_id AS \"Protein.IDs\""
      } else if (identical(cl, "Gene.Names")) {
        "p.\"PG.Genes\" AS \"Gene.Names\""
      } else if (identical(cl, "Assay")) {
        paste0(.db_lit(de_assay), " AS \"Assay\"")
      } else {
        paste0("d.", .db_id(cl))
      }
    }, character(1))
    make_view("v_deps_results", paste0(
      "SELECT ", paste(sel, collapse = ", "),
      " FROM de_results d LEFT JOIN proteins p ON p.protein_id = d.protein_id",
      " ORDER BY d.rowid"))
  }

  # --- v_boxplot_input -------------------------------------------------------
  # Two assays stacked: log2, recomputed from raw, and the imputed one. The
  # ordering reproduces .se_to_long(): assay, then sample, then protein.
  if (!is.na(impn)) {
    make_view("v_boxplot_input", paste0(
      "WITH imp AS (", imputed_src, ") ",
      "SELECT s.sample AS \"Column\", 'log2' AS \"Assay\", ", sql_log2,
      " AS \"Intensity\", s.\"R.Condition\" AS \"Condition\",",
      " s.\"R.Replicate\" AS \"Replicate\", i.protein_id AS \"Protein.IDs\",",
      " 1 AS ao, s.position AS sp, p.position AS pp",
      " FROM intensities i",
      " JOIN proteins p ON p.protein_id = i.protein_id AND p.in_analysis",
      " JOIN samples s ON s.sample = i.sample",
      " WHERE i.assay = 'raw'",
      " UNION ALL ",
      "SELECT s.sample, ", .db_lit(impn), ", i.value, s.\"R.Condition\",",
      " s.\"R.Replicate\", i.protein_id, 2, s.position, p.position",
      " FROM imp i",
      " JOIN proteins p ON p.protein_id = i.protein_id",
      " JOIN samples s ON s.sample = i.sample",
      " ORDER BY ao, sp, pp"))

    # --- v_pca_input ---------------------------------------------------------
    comps <- DBI::dbGetQuery(con,
      "SELECT DISTINCT \"Comparison\" AS c FROM de_results")$c
    adjp <- vapply(comps, function(cp) paste0(
      "MAX(CASE WHEN \"Comparison\" = ", .db_lit(cp), " THEN \"adj.P.Val\" END) AS ",
      .db_id(paste0("adjP_", cp))), character(1))
    sig <- paste(vapply(comps, function(cp) paste0(
      "COALESCE(a.", .db_id(paste0("adjP_", cp)), " < ", alpha, ", FALSE)"),
      character(1)), collapse = " OR ")

    make_view("v_pca_input", paste0(
      "WITH imp AS (", imputed_src, "), ",
      "adj AS (SELECT protein_id, ", paste(adjp, collapse = ", "),
      " FROM de_results GROUP BY protein_id) ",
      "SELECT s.sample AS \"SampleID\", i.protein_id AS \"FeatureID\",",
      " i.value AS \"Intensity\", s.\"R.Condition\" AS \"Condition\",",
      " s.\"R.Replicate\" AS \"Replicate\", ",
      paste(paste0("a.", .db_id(paste0("adjP_", comps))), collapse = ", "),
      ", ", if (length(comps) > 0) sig else "FALSE", " AS \"sig_any\"",
      " FROM imp i",
      " JOIN proteins p ON p.protein_id = i.protein_id",
      " JOIN samples s ON s.sample = i.sample",
      " LEFT JOIN adj a ON a.protein_id = i.protein_id",
      " ORDER BY s.position, p.position"))

    # --- v_matrix_norm / v_matrix_imputed ------------------------------------
    piv_named <- function(src) {
      piv <- vapply(samples, function(s) .db_pivot_expr(NULL, s, s), character(1))
      paste0("SELECT p.protein_id AS \"ProteinGroups\", ",
             paste(paste0("m.", .db_id(samples)), collapse = ", "),
             " FROM proteins p JOIN (SELECT protein_id, ",
             paste(piv, collapse = ", "), " FROM (", src,
             ") GROUP BY protein_id) m ON m.protein_id = p.protein_id",
             " ORDER BY p.position")
    }
    make_view("v_matrix_norm", piv_named(assay_src(imp_from)))
    make_view("v_matrix_imputed", piv_named(imputed_src))
  }

  # --- v_pattern_profiler ----------------------------------------------------
  if (DBI::dbExistsTable(con, "pp_membership")) {
    conds <- DBI::dbGetQuery(con,
      "SELECT \"Condition\" AS c FROM pp_profile GROUP BY \"Condition\", position ORDER BY position")$c
    piv <- vapply(conds, function(cd) paste0(
      "MAX(CASE WHEN \"Condition\" = ", .db_lit(cd), " THEN \"ZScore\" END) AS ",
      .db_id(cd)), character(1))
    make_view("v_pattern_profiler", paste0(
      "SELECT m.protein_id AS \"FeatureID\", m.\"Cluster\", m.\"Membership\", ",
      paste(paste0("z.", .db_id(conds)), collapse = ", "),
      " FROM pp_membership m LEFT JOIN (SELECT protein_id, ",
      paste(piv, collapse = ", "),
      " FROM pp_profile GROUP BY protein_id) z ON z.protein_id = m.protein_id",
      " ORDER BY m.position"))
  }

  invisible(NULL)
}
