
# Load the libraries


# =============================================================================
# Internal helpers
# =============================================================================

# --- Shared numeric filter: hidden input + JS with AND/OR operators ---
.numeric_filter_hidden <- reactable::JS("function() { return null; }")

.numeric_filter_method <- reactable::JS("function(rows, columnId, filterValue) {
  if (!filterValue) return rows;
  var orGroups = filterValue.split('|').map(function(s) { return s.trim(); }).filter(Boolean);
  return rows.filter(function(row) {
    var v = row.values[columnId];
    if (v == null || isNaN(v)) return false;
    return orGroups.some(function(group) {
      var conditions = group.split(',').map(function(s) { return s.trim(); }).filter(Boolean);
      return conditions.every(function(cond) {
        var match = cond.match(/^(>=|<=|!=|<>|>|<|=)?\\s*(.+)$/);
        if (!match) return true;
        var op = match[1] || '=';
        var val = parseFloat(match[2]);
        if (isNaN(val)) return true;
        switch(op) {
          case '>=': return v >= val;
          case '<=': return v <= val;
          case '>':  return v > val;
          case '<':  return v < val;
          case '!=': case '<>': return v !== val;
          case '=':  return v === val;
          default:   return true;
        }
      });
    });
  });
}")


#' Load and validate differential expression data
#' @param input Data frame or path to a TSV/Parquet file
#' @return Validated data frame
#' @noRd
.rl_load_data <- function(input) {
  if (is.character(input) && length(input) == 1) {
    if (!file.exists(input)) {
      stop("File not found: ", input)
    }
    ext <- tolower(tools::file_ext(input))
    if (ext == "tsv" || ext == "txt") {
      if (requireNamespace("readr", quietly = TRUE)) {
        input <- readr::read_tsv(input, show_col_types = FALSE)
      } else {
        input <- read.delim(input, sep = "\t", stringsAsFactors = FALSE)
      }
    } else if (ext == "parquet") {
      if (!requireNamespace("arrow", quietly = TRUE)) {
        stop("The 'arrow' package is required to read Parquet files")
      }
      input <- arrow::read_parquet(input)
    } else if (ext == "csv") {
      if (requireNamespace("readr", quietly = TRUE)) {
        input <- readr::read_csv(input, show_col_types = FALSE)
      } else {
        input <- read.csv(input, stringsAsFactors = FALSE)
      }
    } else {
      stop("Unsupported format: ", ext, ". Use TSV, CSV or Parquet.")
    }
  }

  df <- as.data.frame(input)

  # Validate the required columns
  required <- c("Protein.IDs", "Gene.Names", "logFC", "P.Value", "adj.P.Val",
                 "Change", "Comparison")
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop("Missing required columns: ", paste(missing, collapse = ", "))
  }

  # Type coercion
  df$logFC     <- as.numeric(df$logFC)
  df$P.Value   <- as.numeric(df$P.Value)
  df$adj.P.Val <- as.numeric(df$adj.P.Val)

  df
}


#' Join with protein_quant to add Description and Quant_Pepts
#' @param df Data frame of DE results
#' @param protein_quant Data frame (preprocessing$protein_quant) or path to a TSV.
#'   NULL to skip the join.
#' @return Data frame with the Description and Quant_Pepts columns added
#' @noRd
.rl_join_protein_info <- function(df, protein_quant) {
  if (is.null(protein_quant)) return(df)

  # Load it when a path is given
  if (is.character(protein_quant) && length(protein_quant) == 1) {
    pq <- .rl_load_data(protein_quant)
  } else {
    pq <- as.data.frame(protein_quant)
  }

  # Check the key column
  if (!"PG.ProteinGroups" %in% names(pq)) {
    stop("protein_quant must contain the 'PG.ProteinGroups' column")
  }

  # Extract Description
  desc_col <- if ("PG.ProteinDescriptions" %in% names(pq)) pq$PG.ProteinDescriptions else NA_character_
  info <- data.frame(
    Protein.IDs = pq$PG.ProteinGroups,
    Description = desc_col,
    stringsAsFactors = FALSE
  )

  # Compute the maximum Quant_Pepts
  pept_cols <- grep("^PG\\.NrOfStrippedSequencesUsedForQuantification_", names(pq), value = TRUE)
  if (length(pept_cols) > 0) {
    pept_mat <- as.matrix(pq[, pept_cols, drop = FALSE])
    storage.mode(pept_mat) <- "numeric"
    info$Quant_Pepts <- apply(pept_mat, 1, function(x) {
      vals <- x[!is.na(x)]
      if (length(vals) == 0) 0L else as.integer(max(vals))
    })
  } else {
    info$Quant_Pepts <- NA_integer_
  }

  # Join by Protein.IDs (match preserves the order and avoids duplication)
  idx <- match(df$Protein.IDs, info$Protein.IDs)
  df$Description  <- info$Description[idx]
  df$Quant_Pepts  <- info$Quant_Pepts[idx]

  df
}


#' Teal/green reactable theme
#' @return A reactableTheme object
#' @noRd
.rl_theme <- function() {
  reactableTheme(
    cellPadding = "8px 12px",
    highlightColor = "rgba(2,144,82,0.1)",
    stripedColor = "rgba(180, 220, 210, 0.2)",
    rowSelectedStyle = list(
      backgroundColor = "rgba(2,144,82,0.6)",
      color = "#000000",
      boxShadow = "inset 2px 0 0 0 #ffa62d"
    )
  )
}


#' Dark blue reactable theme (Protein_ID widget)
#' @return A reactableTheme object
#' @noRd
.pl_theme <- function() {
  reactableTheme(
    cellPadding = "8px 12px",
    highlightColor = "rgba(31, 78, 120, 0.10)",
    stripedColor = "rgba(180, 200, 230, 0.20)",
    rowSelectedStyle = list(
      backgroundColor = "rgba(31, 78, 120, 0.6)",
      color = "#ffffff",
      boxShadow = "inset 2px 0 0 0 #ffa62d"
    )
  )
}


#' Dark red reactable theme (Protein_QUANT widget)
#' @return A reactableTheme object
#' @noRd
.ql_theme <- function() {
  reactableTheme(
    cellPadding = "8px 12px",
    highlightColor = "rgba(102, 5, 5, 0.10)",
    stripedColor = "rgba(230, 180, 180, 0.12)",
    rowSelectedStyle = list(
      backgroundColor = "rgba(102, 5, 5, 0.6)",
      color = "#ffffff",
      boxShadow = "inset 2px 0 0 0 #ffa62d"
    )
  )
}


#' Embedded CSS for badges, bars, detail rows, filters and export
#' @return A tags$style object
#' @noRd
.rl_css <- function() {
  ruta <- system.file("css", "results_list.css", package = "NADIA")
  if (!nzchar(ruta) || !file.exists(ruta)) {
    warning("inst/css/results_list.css was not found; the tables will be shown ",
            "without styling.")
    return(htmltools::tags$style(htmltools::HTML("")))
  }
  htmltools::tags$style(
    htmltools::HTML(paste(readLines(ruta, warn = FALSE), collapse = "\n")))
}


#' JavaScript libraries backing the "Export to Excel" button
#'
#' ExcelJS builds the workbook and PapaParse parses the TSV that
#' `Reactable.getDataCSV()` produces. Both are shipped inside the package rather
#' than pulled from a CDN, so the button works without an internet connection and
#' `htmlwidgets::saveWidget(selfcontained = TRUE)` embeds them in the HTML.
#'
#' PapaParse is kept, small as it is, because `getDataCSV()` quotes any field
#' containing the separator, a quote or a newline; splitting the TSV by hand would
#' corrupt those rows.
#'
#' Both libraries are under the MIT licence; their terms travel with them in
#' `inst/js/*/LICENSE`.
#'
#' @return A list of two htmlDependency objects.
#' @noRd
.rl_export_deps <- function() {
  list(
    htmltools::htmlDependency(
      name       = "exceljs",
      version    = "4.4.0",
      src        = c(file = system.file("js", "exceljs", package = "NADIA")),
      script     = "exceljs.min.js",
      all_files  = FALSE),
    htmltools::htmlDependency(
      name       = "papaparse",
      version    = "5.4.1",
      src        = c(file = system.file("js", "papaparse", package = "NADIA")),
      script     = "papaparse.min.js",
      all_files  = FALSE)
  )
}


#' Interface language strings
#' @return A reactableLang object
#' @noRd
.rl_lang <- function() {
  reactableLang(
    searchPlaceholder = "Search...",
    pagePrevious      = "Previous",
    pageNext          = "Next",
    noData            = "No data to display",
    pageSizeOptions   = "Show {rows}",
    pageInfo          = "{rowStart}\u2013{rowEnd} of {rows} Proteins"
  )
}


#' Detail renderer for expandable rows (JS renderer)
#' @param has_assay Logical, whether the data frame has an Assay column
#' @param has_description Logical, whether the data frame has a Description column
#' @return A JS object for the reactable details parameter
#' @noRd
.rl_detail_row <- function(has_assay = TRUE, has_description = FALSE) {
  assay_block <- if (has_assay) {
    "
    var assay = row['Assay'] || '';
    if (assay) {
      html += '<div class=\"detail-row\"><span class=\"detail-label\">Method:</span> ' + assay + '</div>';
    }
    "
  } else {
    ""
  }

  desc_block <- if (has_description) {
    "
    var desc = row['Description'] || '';
    if (desc) {
      html += '<div class=\"detail-row\"><span class=\"detail-label\">Description:</span> ' + desc + '</div>';
    }
    "
  } else {
    ""
  }

  JS(sprintf("function(rowInfo) {
    var row = rowInfo.row;
    var html = '<div class=\"rl-detail\">';

    // Protein.IDs with UniProt links
    var pids = (row['Protein.IDs'] || '').split(';').map(function(s) { return s.trim(); }).filter(Boolean);
    var links = pids.map(function(pid) {
      return '<a href=\"https://www.uniprot.org/uniprot/' + pid + '\" target=\"_blank\">' + pid + '</a>';
    }).join(' \\u00b7 ');
    html += '<div class=\"detail-row\"><span class=\"detail-label\">Proteins:</span> ' + links + '</div>';

    // Full Gene.Names
    var genes = row['Gene.Names'] || '';
    html += '<div class=\"detail-row\"><span class=\"detail-label\">Genes:</span> ' + genes + '</div>';

    // Description
    %s

    // P-value and FDR at full precision
    var pval = row['P.Value'];
    var fdr = row['adj.P.Val'];
    html += '<div class=\"detail-row\"><span class=\"detail-label\">P-value:</span> ' + (pval != null ? pval.toExponential(4) : '') + '</div>';
    html += '<div class=\"detail-row\"><span class=\"detail-label\">FDR:</span> ' + (fdr != null ? fdr.toExponential(4) : '') + '</div>';

    // Assay (when present)
    %s

    html += '</div>';
    return React.createElement('div', { dangerouslySetInnerHTML: { __html: html } });
  }", desc_block, assay_block))
}


#' Build the shared column definitions
#' @param max_abs_lfc Maximum absolute logFC, used to scale the bars
#' @param alpha Significance threshold
#' @param has_assay Logical, whether an Assay column is present
#' @param single_assay Logical, whether there is a single assay
#' @param show_missing Logical, whether to show the Missing% columns
#' @param has_description Logical, whether a Description column is present
#' @param has_quant_pepts Logical, whether a Quant_Pepts column is present
#' @return List of colDef
#' @noRd
.rl_build_columns <- function(max_abs_lfc, alpha, has_assay, single_assay,
                               show_missing, has_description, has_quant_pepts) {

  # JS renderer for Missing% in rating style with a coloured circle
  .missing_cell_js <- JS("function(cellInfo) {
    var pct = cellInfo.value;
    if (pct == null || isNaN(pct)) return '';
    var rounded = Math.round(pct);
    var color;
    if (pct === 0) color = '#02905A';
    else if (pct <= 12.5) color = '#f5c842';
    else if (pct <= 25) color = '#e8a735';
    else if (pct <= 37.5) color = '#e07b3c';
    else color = '#d94545';
    return '\\u25cf ' + rounded;
  }")

  .missing_style_js <- JS("function(rowInfo, column) {
    var pct = rowInfo.row[column.id];
    if (pct == null || isNaN(pct)) return {};
    var color;
    if (pct === 0) color = '#02905A';
    else if (pct <= 12.5) color = '#f5c842';
    else if (pct <= 25) color = '#e8a735';
    else if (pct <= 37.5) color = '#e07b3c';
    else color = '#d94545';
    return { color: color, fontWeight: 600 };
  }")

  cols <- list(
    # --- Comparison (1st) ---
    Comparison = colDef(
      name = "Comparison",
      width = 120,
      align = "center"
    ),

    # --- Protein Groups (2nd) ---
    Protein.IDs = colDef(
      name = "Protein Groups",
      width = 160,
      html = TRUE,
      cell = JS("function(cellInfo) {
        var val = cellInfo.value || '';
        var pids = val.split(';').map(function(s) { return s.trim(); }).filter(Boolean);
        var first = pids[0] || val;
        if (pids.length > 1) {
          return '<strong>' + first + '</strong> <span class=\"protein-count\">+' + (pids.length - 1) + '</span>';
        }
        return '<strong>' + first + '</strong>';
      }")
    )
  )

  # --- Description (3rd, conditional) ---
  if (has_description) {
    cols$Description <- colDef(
      name = "Description",
      minWidth = 250
    )
  }

  # --- Gene Names (4th) ---
  cols$Gene.Names <- colDef(
    name = "Gene Names",
    minWidth = 140,
    html = TRUE,
    cell = JS("function(cellInfo) {
      var val = cellInfo.value || '';
      var genes = val.split(';').map(function(s) { return s.trim(); }).filter(Boolean);
      var first = genes[0] || val;
      if (genes.length > 1) {
        return first + ' <span class=\"protein-count\">+' + (genes.length - 1) + '</span>';
      }
      return first;
    }"),
    style = list(alignItems = "center")
  )

  # --- Quant_Pepts (5th, conditional) ---
  if (has_quant_pepts) {
    cols$Quant_Pepts <- colDef(
      name = "Quant Pepts",
      width = 110,
      align = "center",
      filterable = TRUE,
      filterInput = .numeric_filter_hidden,
      filterMethod = .numeric_filter_method
    )
  }

  # --- Change (badge) ---
  cols$Change <- colDef(
    name = "Change",
    width = 130,
    align = "center",
    html = TRUE,
    cell = JS("function(cellInfo) {
      var val = cellInfo.value;
      var cls = 'tag status-grey';
      if (val === 'Up') cls = 'tag status-green';
      else if (val === 'Down') cls = 'tag status-red';
      return '<span class=\"' + cls + '\">' + val + '</span>';
    }")
  )

  # --- logFC (value + bar coloured by Change) ---
  cols$logFC <- colDef(
    name = "log\u2082 FC",
    width = 160,
    align = "center",
    html = TRUE,
    filterable = TRUE,
    filterInput = .numeric_filter_hidden,
    filterMethod = .numeric_filter_method,
    cell = JS(sprintf("function(cellInfo) {
      var val = cellInfo.value;
      var change = cellInfo.row['Change'];
      var maxLfc = %s;
      var pct = Math.abs(val) / maxLfc * 100;
      var barClass = 'lfc-bar nochange';
      var valueClass = 'lfc-value nochange';
      if (change === 'Up') {
        barClass = 'lfc-bar up';
        valueClass = 'lfc-value up';
      } else if (change === 'Down') {
        barClass = 'lfc-bar down';
        valueClass = 'lfc-value down';
      }
      var formatted = (val > 0 ? '+' : '') + val.toFixed(3);
      return '<div class=\"lfc-bar-container\">' +
        '<span class=\"' + valueClass + '\">' + formatted + '</span>' +
        '<div class=\"lfc-bar-wrapper\">' +
        '<div class=\"' + barClass + '\" style=\"width:' + pct.toFixed(1) + '%%\"></div>' +
        '</div></div>';
    }", max_abs_lfc))
  )

  # --- P-value ---
  cols$P.Value <- colDef(
    name = "P-value",
    width = 110,
    align = "right",
    html = TRUE,
    cell = JS("function(cellInfo) {
      var val = cellInfo.value;
      if (val == null || isNaN(val)) return '';
      return val.toExponential(2);
    }")
  )

  # --- FDR ---
  cols$adj.P.Val <- colDef(
    name = "FDR",
    width = 120,
    align = "right",
    html = TRUE,
    filterable = TRUE,
    filterInput = .numeric_filter_hidden,
    filterMethod = .numeric_filter_method,
    cell = JS(sprintf("function(cellInfo) {
      var val = cellInfo.value;
      if (val == null || isNaN(val)) return '';
      var formatted = val.toExponential(2);
      if (val < %s) {
        return '<strong style=\"color: #0E6655;\">' + formatted + '</strong>';
      }
      return formatted;
    }", alpha))
  )

  # --- Assay ---
  if (has_assay) {
    cols$Assay <- colDef(name = "Method", width = 130, show = !single_assay)
  }

  # --- Missing% in rating style with a coloured circle ---
  if (show_missing) {
    cols$MissingGlobal <- colDef(
      name = "% Missing", width = 110, align = "center", class = "missing-cell",
      header = function(value) htmltools::tags$span(title = "Global percentage of NAs in the comparison", value),
      filterable = TRUE, filterInput = .numeric_filter_hidden, filterMethod = .numeric_filter_method,
      cell = .missing_cell_js, style = .missing_style_js
    )
    cols$MissingPCT1 <- colDef(
      name = "% Group 1", width = 110, align = "center", class = "missing-cell",
      header = function(value) htmltools::tags$span(title = "Percentage of NAs in the numerator", value),
      filterable = TRUE, filterInput = .numeric_filter_hidden, filterMethod = .numeric_filter_method,
      cell = .missing_cell_js, style = .missing_style_js
    )
    cols$MissingPCT2 <- colDef(
      name = "% Group 2", width = 110, align = "center", class = "missing-cell",
      header = function(value) htmltools::tags$span(title = "Percentage of NAs in the denominator", value),
      filterable = TRUE, filterInput = .numeric_filter_hidden, filterMethod = .numeric_filter_method,
      cell = .missing_cell_js, style = .missing_style_js
    )
  }

  cols
}


# =============================================================================
# Main function
# =============================================================================

#' Interactive Reactable Table for Differential Expression Results
#'
#' Builds a professionally formatted reactable table to display protein
#' differential expression results. Works both in Shiny (it returns a reactable
#' object) and standalone via \code{results_list_widget()}.
#'
#' @param data Data frame or path to a TSV/CSV/Parquet file with DE results.
#'   Required columns: Protein.IDs, Gene.Names, logFC, P.Value, adj.P.Val,
#'   Change, Comparison. Optional: Assay, MissingGlobal, MissingPCT1, MissingPCT2
#' @param protein_quant Data frame (preprocessing$protein_quant) or path to a
#'   Protein_QUANT_*.tsv file. If not NULL, the Description and Quant_Pepts
#'   columns are added.
#' @param comparisons Vector of comparisons to include (NULL = all of them)
#' @param ain Vector of assays to keep (NULL = all of them)
#' @param alpha Significance threshold used to highlight the FDR (default: 0.05)
#' @param lfc_thr log2 fold-change threshold (default: 0, reserved for future use)
#' @param page_size Rows per page (default: 15)
#' @param height Table height in pixels (default: 720)
#' @param show_missing Show the missingness percentage columns (default: TRUE)
#' @param element_id Element ID for the Reactable JS API (default: NULL)
#' @param selection Selection type: "single", "multiple", or NULL (default: NULL)
#' @param searchable Enable the internal search box (default: TRUE)
#'
#' @return A reactable object
#'
#' @examples
#' data(nadia_dia)
#' res <- process_proteomics(nadia_dia, verbose = FALSE)
#'
#' # protein_quant adds the Description and Quant_Pepts columns
#' tbl <- results_list_reactable(
#'   res$DEPs_results,
#'   protein_quant = nadia_dia$protein_quant,
#'   comparisons   = "B-A"
#' )
#' class(tbl)
#'
#' @export
results_list_reactable <- function(
    data,
    protein_quant = NULL,
    comparisons = NULL,
    ain = NULL,
    alpha = 0.05,
    lfc_thr = 0,
    page_size = 15,
    height = 720,
    show_missing = TRUE,
    element_id = NULL,
    selection = NULL,
    searchable = TRUE
) {

  # --- Load and validate ---
  df <- .rl_load_data(data)

  # --- Join with protein_quant ---
  df <- .rl_join_protein_info(df, protein_quant)

  # --- Filtering ---
  if (!is.null(ain) && "Assay" %in% names(df)) {
    df <- df[df$Assay %in% ain, , drop = FALSE]
  }
  if (!is.null(comparisons)) {
    df <- df[df$Comparison %in% comparisons, , drop = FALSE]
  }
  if (nrow(df) == 0) {
    stop("No data left after applying the comparison/assay filters")
  }

  if (nrow(df) > 15000 && is.null(comparisons)) {
    message("Note: ", format(nrow(df), big.mark = ","),
            " rows. Consider filtering by 'comparisons' for better performance.")
  }

  # --- Detect the optional columns ---
  has_missing <- all(c("MissingGlobal", "MissingPCT1", "MissingPCT2") %in% names(df))
  show_missing <- show_missing && has_missing
  has_assay <- "Assay" %in% names(df)
  single_assay <- has_assay && length(unique(df$Assay)) == 1
  has_description <- "Description" %in% names(df)
  has_quant_pepts <- "Quant_Pepts" %in% names(df)

  max_abs_lfc <- max(abs(df$logFC), na.rm = TRUE)
  if (max_abs_lfc == 0) max_abs_lfc <- 1

  cols <- .rl_build_columns(max_abs_lfc, alpha, has_assay, single_assay,
                             show_missing, has_description, has_quant_pepts)

  # --- Reorder the data frame columns (reactable uses this visual order) ---
  desired_order <- c("Comparison", "Protein.IDs", "Description", "Gene.Names",
                     "Quant_Pepts", "Change", "logFC", "P.Value", "adj.P.Val",
                     "Assay", "MissingGlobal", "MissingPCT1", "MissingPCT2")
  desired_order <- intersect(desired_order, names(df))
  df <- df[, c(desired_order, setdiff(names(df), desired_order)), drop = FALSE]

  reactable(
    df,
    elementId     = element_id,
    defaultSorted = "adj.P.Val",
    defaultPageSize    = page_size,
    showPageSizeOptions = TRUE,
    pageSizeOptions = c(15, 30, 50, 100),
    resizable   = TRUE,
    selection   = selection,
    onClick     = if (!is.null(selection)) "select" else NULL,
    defaultColDef = colDef(
      align = "left",
      headerStyle = list(
        background  = "rgba(14, 102, 85, 0.9)",
        color       = "#ffffff",
        height      = "45px",
        display     = "flex",
        alignItems  = "center"
      ),
      class = JS("function(rowInfo, column, state) {
        for (let i = 0; i < state.sorted.length; i++) {
          if (state.sorted[i].id === column.id) {
            return 'sorted'
          }
        }
      }"),
      style = list(height = "48px", display = "flex", alignItems = "center")
    ),
    columns    = cols,
    wrap       = FALSE,
    class      = "rl-table",
    rowStyle   = if (!is.null(selection)) list(cursor = "pointer") else NULL,
    highlight  = TRUE,
    searchable = searchable,
    height     = height,
    striped    = TRUE,
    theme      = .rl_theme(),
    language   = .rl_lang(),
    details    = .rl_detail_row(has_assay, has_description)
  )
}


# =============================================================================
# Wrapper standalone / Quarto
# =============================================================================

#' Full Widget with Filters, Search, Export and CSS
#'
#' Wraps \code{results_list_reactable()} with interactive crosstalk filters
#' (Comparison, Change, Method), a search box, an export-to-Excel button and
#' embedded CSS. The result is browsable: printing it at the console opens it
#' automatically in the RStudio Viewer or in the browser.
#'
#' @inheritParams results_list_reactable
#' @param element_id Element ID (default: "deps_table")
#'
#' @return A browsable htmltools object
#'
#' @examples
#' if (requireNamespace("crosstalk", quietly = TRUE)) {
#'   data(nadia_dia)
#'   res <- process_proteomics(nadia_dia, verbose = FALSE)
#'
#'   # Printing the widget opens it in the Viewer / browser
#'   w <- results_list_widget(res$DEPs_results,
#'                            protein_quant = nadia_dia$protein_quant)
#'   print(class(w))
#' }
#'
#' @export
results_list_widget <- function(
    data,
    protein_quant = NULL,
    comparisons = NULL,
    ain = NULL,
    alpha = 0.05,
    lfc_thr = 0,
    page_size = 15,
    height = 720,
    show_missing = TRUE,
    element_id = "deps_table",
    selection = NULL,
    searchable = TRUE
) {

  if (!requireNamespace("crosstalk", quietly = TRUE)) {
    stop("The 'crosstalk' package is required for the interactive filters. ",
         "Install it with install.packages('crosstalk')")
  }

  # --- Load, join and pre-filter ---
  df <- .rl_load_data(data)
  df <- .rl_join_protein_info(df, protein_quant)

  if (!is.null(ain) && "Assay" %in% names(df)) {
    df <- df[df$Assay %in% ain, , drop = FALSE]
  }
  if (!is.null(comparisons)) {
    df <- df[df$Comparison %in% comparisons, , drop = FALSE]
  }
  if (nrow(df) == 0) {
    stop("No data left after applying the comparison/assay filters")
  }

  # --- Detect the optional columns ---
  has_missing <- all(c("MissingGlobal", "MissingPCT1", "MissingPCT2") %in% names(df))
  show_missing_cols <- show_missing && has_missing
  has_assay <- "Assay" %in% names(df)
  single_assay <- has_assay && length(unique(df$Assay)) == 1
  has_description <- "Description" %in% names(df)
  has_quant_pepts <- "Quant_Pepts" %in% names(df)

  max_abs_lfc <- max(abs(df$logFC), na.rm = TRUE)
  if (max_abs_lfc == 0) max_abs_lfc <- 1

  # --- Reorder the data frame columns (reactable uses this visual order) ---
  desired_order <- c("Comparison", "Protein.IDs", "Description", "Gene.Names",
                     "Quant_Pepts", "Change", "logFC", "P.Value", "adj.P.Val",
                     "Assay", "MissingGlobal", "MissingPCT1", "MissingPCT2")
  desired_order <- intersect(desired_order, names(df))
  df <- df[, c(desired_order, setdiff(names(df), desired_order)), drop = FALSE]

  # --- SharedData for crosstalk ---
  shared_data <- crosstalk::SharedData$new(df)

  # --- CSS ---
  css <- .rl_css()

  # --- ExcelJS and PapaParse, shipped with the package (see .rl_export_deps) ---
  export_deps <- .rl_export_deps()

  # --- JavaScript: toggle filters, clear filters, export to Excel ---
  js_code <- tags$script(HTML(sprintf("
    var rlFiltersVisible = false;

    function rlToggleFilters() {
      var container = document.querySelector('.rl-filters-container');
      var btn = document.querySelector('.rl-btn-toggle-filters');
      if (rlFiltersVisible) {
        container.style.maxHeight = '0';
        container.style.opacity = '0';
        container.style.marginBottom = '0';
        container.style.padding = '0';
        container.style.borderWidth = '0';
        container.style.boxShadow = 'none';
        container.style.overflow = 'hidden';
        btn.classList.add('filters-hidden');
      } else {
        container.style.maxHeight = '500px';
        container.style.opacity = '1';
        container.style.marginBottom = '1rem';
        container.style.padding = '1.5rem';
        container.style.borderWidth = '1px';
        container.style.boxShadow = '0 2px 4px rgba(0, 0, 0, 0.05)';
        container.style.overflow = 'visible';
        btn.classList.remove('filters-hidden');
      }
      rlFiltersVisible = !rlFiltersVisible;
    }

    function rlApplyNumericFilter(columnId, value) {
      Reactable.setFilter('%s', columnId, value || undefined);
    }

    function rlClearFilters() {
      var selects = document.querySelectorAll('.rl-filter-item .selectized');
      selects.forEach(function(sel) {
        if (sel.selectize) sel.selectize.clear();
      });
      var numInputs = document.querySelectorAll('.rl-numeric-input');
      numInputs.forEach(function(inp) {
        inp.value = '';
        var col = inp.getAttribute('data-column');
        if (col) Reactable.setFilter('%s', col, undefined);
      });
      var searchInput = document.querySelector('.rl-search-input');
      if (searchInput) {
        searchInput.value = '';
        Reactable.setSearch('%s', '');
      }
    }

    async function rlExportExcel() {
      try {
        var tsv = Reactable.getDataCSV('%s', { sep: '\\t' });
        var parseResult = Papa.parse(tsv, { header: true, delimiter: '\\t', skipEmptyLines: true });
        var rows = parseResult.data;
        var headers = parseResult.meta.fields;

        var headerMap = {
          'Gene.Names': 'Gene Names',
          'Comparison': 'Comparison',
          'Change': 'Change',
          'logFC': 'log2 FC',
          'adj.P.Val': 'FDR',
          'P.Value': 'P-value',
          'Protein.IDs': 'Protein Groups',
          'Description': 'Description',
          'Quant_Pepts': 'Quant Pepts',
          'Assay': 'Method',
          'MissingGlobal': '%%Missing',
          'MissingPCT1': '%%Group1',
          'MissingPCT2': '%%Group2'
        };

        var wb = new ExcelJS.Workbook();
        var ws = wb.addWorksheet('DE Results');

        var visibleHeaders = headers;

        ws.columns = visibleHeaders.map(function(h) {
          return {
            header: headerMap[h] || h,
            key: h,
            width: (h === 'Protein.IDs' || h === 'Gene.Names') ? 25 :
                   (h === 'Description') ? 40 : 15
          };
        });

        var headerRow = ws.getRow(1);
        headerRow.font = { bold: true, color: { argb: 'FFFFFFFF' }, size: 11 };
        headerRow.fill = { type: 'pattern', pattern: 'solid', fgColor: { argb: 'FF0E6655' } };
        headerRow.alignment = { vertical: 'middle', horizontal: 'center' };
        headerRow.height = 30;

        rows.forEach(function(row) {
          var rowData = {};
          visibleHeaders.forEach(function(h) { rowData[h] = row[h]; });
          var addedRow = ws.addRow(rowData);

          ['logFC', 'P.Value', 'adj.P.Val', 'MissingGlobal', 'MissingPCT1', 'MissingPCT2', 'Quant_Pepts'].forEach(function(col) {
            if (visibleHeaders.indexOf(col) === -1) return;
            var cell = addedRow.getCell(col);
            if (cell && cell.value) {
              var num = parseFloat(cell.value);
              if (!isNaN(num)) cell.value = num;
            }
          });
        });

        ws.eachRow(function(row) {
          row.eachCell(function(cell) {
            cell.border = {
              top: { style: 'thin' }, left: { style: 'thin' },
              bottom: { style: 'thin' }, right: { style: 'thin' }
            };
          });
        });

        var buffer = await wb.xlsx.writeBuffer();
        var blob = new Blob([buffer], { type: 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet' });
        var url = window.URL.createObjectURL(blob);
        var a = document.createElement('a');
        a.href = url;
        a.download = 'DEPs_Results_' + new Date().toISOString().split('T')[0] + '.xlsx';
        a.click();
        window.URL.revokeObjectURL(url);
      } catch(e) {
        console.error('Error while exporting:', e);
      }
    }
  ", element_id, element_id, element_id, element_id)))

  # --- Search bar + action buttons ---
  search_actions <- div(class = "rl-search-actions",
    tags$input(
      type = "search",
      placeholder = "Search...",
      class = "rl-search-input",
      oninput = sprintf("Reactable.setSearch('%s', this.value)", element_id)
    ),
    div(class = "rl-action-buttons",
      tags$button(
        class = "rl-btn-action rl-btn-toggle-filters filters-hidden",
        onclick = "rlToggleFilters()",
        title = "Show/Hide filters",
        HTML('<svg xmlns="http://www.w3.org/2000/svg" width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><polygon points="22 3 2 3 10 12.46 10 19 14 21 14 12.46 22 3"></polygon></svg>')
      ),
      tags$button(
        class = "rl-btn-action",
        onclick = "rlClearFilters()",
        title = "Clear filters",
        HTML('<svg xmlns="http://www.w3.org/2000/svg" width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M3 6h18"></path><path d="M8 6V4c0-1 1-2 2-2h4c1 0 2 1 2 2v2"></path><path d="M19 6l-1 14c0 1-1 2-2 2H8c-1 0-2-1-2-2L5 6"></path><line x1="1" y1="1" x2="23" y2="23" stroke="#E63946" stroke-width="2"></line></svg>')
      ),
      tags$button(
        class = "rl-btn-action",
        onclick = "rlExportExcel()",
        title = "Export to Excel",
        HTML('<svg xmlns="http://www.w3.org/2000/svg" width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4"></path><polyline points="7 10 12 15 17 10"></polyline><line x1="12" y1="15" x2="12" y2="3"></line></svg>')
      )
    )
  )

  # --- Crosstalk + numeric filter panel ---
  filter_items <- list(
    div(class = "rl-filter-item",
      tags$label(class = "rl-filter-label", "Comparison"),
      crosstalk::filter_select(
        id = "rl_filter_comparison", label = NULL,
        sharedData = shared_data, group = ~Comparison, multiple = TRUE
      )
    ),
    div(class = "rl-filter-item",
      tags$label(class = "rl-filter-label", "Change"),
      crosstalk::filter_select(
        id = "rl_filter_change", label = NULL,
        sharedData = shared_data, group = ~Change, multiple = TRUE
      )
    )
  )

  if (has_assay && !single_assay) {
    filter_items <- c(filter_items, list(
      div(class = "rl-filter-item",
        tags$label(class = "rl-filter-label", "Method"),
        crosstalk::filter_select(
          id = "rl_filter_assay", label = NULL,
          sharedData = shared_data, group = ~Assay, multiple = TRUE
        )
      )
    ))
  }

  # --- Numeric filters (same row, same styles) ---
  .make_numeric_filter <- function(label, column_id, placeholder) {
    div(class = "rl-filter-item",
      tags$label(class = "rl-filter-label", label),
      tags$input(
        type = "text",
        class = "rl-numeric-input",
        `data-column` = column_id,
        placeholder = placeholder,
        oninput = sprintf("rlApplyNumericFilter('%s', this.value)", column_id)
      ),
      tags$span(class = "rl-filter-hint", "AND: >=2, <=5 \u00b7 OR: <=-1 | >=1")
    )
  }

  filter_items <- c(filter_items, list(
    .make_numeric_filter("log\u2082 FC", "logFC", "\u2264 -1, \u2265 2 ..."),
    .make_numeric_filter("FDR", "adj.P.Val", "\u2264 0.05 ...")
  ))

  if (has_quant_pepts) {
    filter_items <- c(filter_items, list(
      .make_numeric_filter("Quant Pepts", "Quant_Pepts", "\u2265 3 ...")
    ))
  }

  if (show_missing_cols) {
    filter_items <- c(filter_items, list(
      .make_numeric_filter("% Missing", "MissingGlobal", "\u2264 25 ..."),
      .make_numeric_filter("% Group 1", "MissingPCT1", "\u2264 50 ..."),
      .make_numeric_filter("% Group 2", "MissingPCT2", "\u2264 50 ...")
    ))
  }

  filters_panel <- div(class = "rl-filters-container",
    div(class = "rl-filters-row", filter_items)
  )

  # --- Columns ---
  cols <- .rl_build_columns(max_abs_lfc, alpha, has_assay, single_assay,
                             show_missing_cols, has_description, has_quant_pepts)

  # --- Reactable table with SharedData ---
  tbl <- reactable(
    shared_data,
    elementId     = element_id,
    defaultSorted = "adj.P.Val",
    defaultPageSize    = page_size,
    showPageSizeOptions = TRUE,
    pageSizeOptions = c(15, 30, 50, 100),
    resizable   = TRUE,
    selection   = selection,
    onClick     = if (!is.null(selection)) "select" else NULL,
    defaultColDef = colDef(
      align = "left",
      headerStyle = list(
        background  = "rgba(14, 102, 85, 0.9)",
        color       = "#ffffff",
        height      = "45px",
        display     = "flex",
        alignItems  = "center"
      ),
      class = JS("function(rowInfo, column, state) {
        for (let i = 0; i < state.sorted.length; i++) {
          if (state.sorted[i].id === column.id) {
            return 'sorted'
          }
        }
      }"),
      style = list(height = "48px", display = "flex", alignItems = "center")
    ),
    columns    = cols,
    wrap       = FALSE,
    class      = "rl-table",
    rowStyle   = if (!is.null(selection)) list(cursor = "pointer") else NULL,
    highlight  = TRUE,
    searchable = searchable,
    height     = height,
    striped    = TRUE,
    theme      = .rl_theme(),
    language   = .rl_lang(),
    details    = .rl_detail_row(has_assay, has_description)
  )

  browsable(tagList(
    css,
    export_deps,
    js_code,
    search_actions,
    filters_panel,
    tbl
  ))
}


# =============================================================================
# =============================================================================
# Helpers Protein_ID (pl-*)
# =============================================================================
# =============================================================================

#' Load and validate a Protein_ID file
#' @param input Data frame or path to a TSV/CSV/Parquet file in Protein_ID format
#' @return Validated data frame
#' @noRd
.pl_load_data <- function(input) {
  if (is.character(input) && length(input) == 1) {
    if (!file.exists(input)) stop("File not found: ", input)
    ext <- tolower(tools::file_ext(input))
    if (ext %in% c("tsv", "txt")) {
      if (requireNamespace("readr", quietly = TRUE)) {
        input <- readr::read_tsv(input, show_col_types = FALSE)
      } else {
        input <- utils::read.delim(input, sep = "\t", stringsAsFactors = FALSE)
      }
    } else if (ext == "csv") {
      if (requireNamespace("readr", quietly = TRUE)) {
        input <- readr::read_csv(input, show_col_types = FALSE)
      } else {
        input <- utils::read.csv(input, stringsAsFactors = FALSE)
      }
    } else if (ext == "parquet") {
      if (!requireNamespace("arrow", quietly = TRUE)) {
        stop("The 'arrow' package is required to read Parquet files")
      }
      input <- arrow::read_parquet(input)
    } else {
      stop("Unsupported format: ", ext, ". Use TSV, CSV or Parquet.")
    }
  }

  df <- as.data.frame(input)

  if (!"PG.ProteinGroups" %in% names(df)) {
    stop("Required column 'PG.ProteinGroups' not found. ",
         "The file must be a Protein_ID exported by preprocess_spectronaut().")
  }
  if (!any(grepl("^PG\\.NrOfPrecursorsIdentified_", names(df)))) {
    stop("No sample columns were detected (PG.NrOfPrecursorsIdentified_*). ",
         "Check that the file really is a Protein_ID.")
  }
  df
}


#' Parse Protein_ID column names into (metric, condition, replicate)
#' @param colnames_vec Vector of data frame column names
#' @return Sorted data frame: column, metric, condition, replicate, coding
#' @noRd
.pl_parse_samples <- function(colnames_vec) {
  metrics <- c(
    "PG.NrOfPrecursorsIdentified",
    "PG.NrOfStrippedSequencesIdentified",
    "PG.Coverage",
    "PG.Cscore.RunWise"
  )

  rows <- list()
  for (m in metrics) {
    pat <- paste0("^", gsub("\\.", "\\\\.", m), "_(.+?)_(\\d+)$")
    hits <- regmatches(colnames_vec, regexec(pat, colnames_vec))
    for (i in seq_along(hits)) {
      h <- hits[[i]]
      if (length(h) == 3) {
        rows[[length(rows) + 1]] <- data.frame(
          column = colnames_vec[i],
          metric = m,
          condition = h[2],
          replicate = as.integer(h[3]),
          coding = paste(h[2], h[3], sep = "_"),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(rows) == 0) {
    return(data.frame(column = character(0), metric = character(0),
                      condition = character(0), replicate = integer(0),
                      coding = character(0), stringsAsFactors = FALSE))
  }

  df <- do.call(rbind, rows)
  df$metric <- factor(df$metric, levels = metrics)
  df <- df[order(df$metric, df$condition, df$replicate), ]
  df$metric <- as.character(df$metric)
  rownames(df) <- NULL
  df
}


#' Colour palette per condition
#' @param conditions Vector of condition codes (e.g., c("A","B","C","D"))
#' @return Named vector of hex colours
#' @noRd
.pl_condition_palette <- function(conditions) {
  defaults <- c(
    A = "#4F81BD", B = "#9BBB59", C = "#F79646", D = "#8064A2",
    E = "#4BACC6", F = "#C0504D", G = "#9F8A76", H = "#646464",
    I = "#FFC000", J = "#E377C2", K = "#2E7D32", L = "#6A1B9A"
  )
  conditions <- unique(conditions)
  out <- setNames(rep(NA_character_, length(conditions)), conditions)
  for (cc in conditions) {
    if (cc %in% names(defaults)) out[cc] <- defaults[[cc]]
  }
  missing_idx <- which(is.na(out))
  if (length(missing_idx) > 0) {
    if (requireNamespace("scales", quietly = TRUE)) {
      out[missing_idx] <- scales::hue_pal()(length(missing_idx))
    } else {
      out[missing_idx] <- rep("#6c757d", length(missing_idx))
    }
  }
  out
}


#' Build the colDefs for the Protein_ID table
#' @param sample_map Output of .pl_parse_samples()
#' @param max_per_col Named numeric vector with the per-column maximum (for the data bars)
#' @return List of colDef
#' @noRd
.pl_build_columns <- function(sample_map, max_per_col) {
  cols <- list()

  cols$PG.ProteinGroups <- colDef(
    name = "Protein Groups",
    minWidth = 160,
    sticky = "left",
    html = TRUE,
    cell = JS("function(cellInfo) {
      var val = cellInfo.value || '';
      var pids = val.split(';').map(function(s){return s.trim();}).filter(Boolean);
      var first = pids[0] || val;
      if (pids.length > 1) {
        return '<strong>' + first + '</strong> <span class=\"protein-count\">+' + (pids.length - 1) + '</span>';
      }
      return '<strong>' + first + '</strong>';
    }")
  )

  cols$PG.ProteinDescriptions <- colDef(
    name = "Descriptions",
    minWidth = 220,
    sticky = "left"
  )

  cols$PG.Genes <- colDef(
    name = "Gene Names",
    minWidth = 120,
    sticky = "left",
    html = TRUE,
    cell = JS("function(cellInfo) {
      var val = cellInfo.value || '';
      var genes = val.split(';').map(function(s){return s.trim();}).filter(Boolean);
      var first = genes[0] || val;
      if (genes.length > 1) {
        return first + ' <span class=\"protein-count\">+' + (genes.length - 1) + '</span>';
      }
      return first;
    }")
  )

  cols$PG.MolecularWeight <- colDef(
    name = "MW [kDa]",
    width = 95,
    align = "right",
    cell = JS("function(cellInfo) {
      var val = cellInfo.value;
      if (val == null || isNaN(val)) return '';
      return (val / 1000).toFixed(2);
    }")
  )

  # First column of each metric -> anchor for the border and the aggregated filter
  first_cols_by_metric <- vapply(
    split(sample_map$column, sample_map$metric),
    function(x) x[1], character(1)
  )

  for (i in seq_len(nrow(sample_map))) {
    col_id    <- sample_map$column[i]
    cond      <- sample_map$condition[i]
    replicate <- sample_map$replicate[i]
    metric    <- sample_map$metric[i]
    max_val   <- max_per_col[[col_id]]
    if (!is.finite(max_val) || max_val <= 0) max_val <- 1

    fmt_js <- switch(
      metric,
      "PG.Coverage"        = "val.toFixed(1)",
      "PG.Cscore.RunWise"  = "val.toFixed(2)",
      "Math.round(val)"
    )

    cell_js <- JS(sprintf("function(cellInfo) {
      var val = cellInfo.value;
      if (val == null || isNaN(val)) return '<span class=\"pl-plain-value pl-bar-empty\">\u2013</span>';
      return '<span class=\"pl-plain-value\">' + (%s) + '</span>';
    }", fmt_js))

    is_anchor <- col_id %in% first_cols_by_metric
    cell_class <- if (is_anchor) "pl-sample-cell pl-group-boundary" else "pl-sample-cell"
    hdr_class  <- if (is_anchor) paste0("pl-hdr-", cond, " pl-group-boundary") else paste0("pl-hdr-", cond)
    fm <- if (is_anchor) {
      .pl_agg_filter_method(sample_map$column[sample_map$metric == metric])
    } else {
      .numeric_filter_method
    }

    cols[[col_id]] <- colDef(
      name = paste(cond, replicate, sep = "_"),
      width = 72,
      align = "center",
      html = TRUE,
      class = cell_class,
      headerClass = hdr_class,
      filterable = TRUE,
      filterInput = .numeric_filter_hidden,
      filterMethod = fm,
      cell = cell_js
    )
  }

  cols
}


#' Build the columnGroups for the Protein_ID table
#' @param sample_map Output of .pl_parse_samples()
#' @param static_cols Vector of the static column names that are present
#' @return List of colGroup
#' @noRd
.pl_build_column_groups <- function(sample_map, static_cols) {
  groups <- list()

  if (length(static_cols) > 0) {
    groups[[length(groups) + 1]] <- colGroup(
      name = "PROTEIN ANNOTATION",
      columns = static_cols,
      sticky = "left"
    )
  }

  metric_labels <- list(
    "PG.NrOfPrecursorsIdentified"        = "# Uniq. PSMs Identified",
    "PG.NrOfStrippedSequencesIdentified" = "# Uniq. Pepts Identified",
    "PG.Coverage"                        = "Coverage [%]",
    "PG.Cscore.RunWise"                  = "Spectronaut Cscore"
  )

  for (m in names(metric_labels)) {
    cols <- sample_map$column[sample_map$metric == m]
    if (length(cols) == 0) next
    groups[[length(groups) + 1]] <- colGroup(
      name = metric_labels[[m]],
      columns = cols
    )
  }

  groups
}


#' Expandable detail row for Protein_ID
#' @return A JS object
#' @noRd
.pl_detail_row <- function() {
  JS("function(rowInfo) {
    var row = rowInfo.row;
    var html = '<div class=\"rl-detail\">';

    var pids = (row['PG.ProteinGroups'] || '').split(';').map(function(s){return s.trim();}).filter(Boolean);
    var links = pids.map(function(pid){
      return '<a href=\"https://www.uniprot.org/uniprot/' + pid + '\" target=\"_blank\">' + pid + '</a>';
    }).join(' \\u00b7 ');
    html += '<div class=\"detail-row\"><span class=\"detail-label\">Proteins:</span> ' + links + '</div>';

    html += '<div class=\"detail-row\"><span class=\"detail-label\">Genes:</span> ' + (row['PG.Genes'] || '') + '</div>';
    html += '<div class=\"detail-row\"><span class=\"detail-label\">Description:</span> ' + (row['PG.ProteinDescriptions'] || '') + '</div>';

    var mw = row['PG.MolecularWeight'];
    if (mw != null && !isNaN(mw)) {
      html += '<div class=\"detail-row\"><span class=\"detail-label\">MW:</span> ' + (mw / 1000).toFixed(2) + ' kDa</div>';
    }

    html += '</div>';
    return React.createElement('div', { dangerouslySetInnerHTML: { __html: html } });
  }")
}


#' Aggregated OR numeric filter across several columns
#' @description Builds a filterMethod that evaluates the numeric expression
#'   (>=, <=, >, <, =, !=) against ANY of the given columns: a single column
#'   satisfying it is enough to keep the row. It supports the compound AND
#'   (comma) and OR (pipe) operators, just like \code{.numeric_filter_method}.
#' @param col_ids Vector of column names over which the OR is applied.
#' @return A JS object to be used as filterMethod in colDef.
#' @noRd
.pl_agg_filter_method <- function(col_ids) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("'jsonlite' is required for the aggregated filters.")
  }
  cols_json <- jsonlite::toJSON(col_ids)
  JS(sprintf("function(rows, columnId, filterValue) {
    var cols = %s;
    if (!filterValue) return rows;
    var orGroups = filterValue.split('|').map(function(s){return s.trim();}).filter(Boolean);
    return rows.filter(function(row) {
      return orGroups.some(function(group) {
        var conditions = group.split(',').map(function(s){return s.trim();}).filter(Boolean);
        return conditions.every(function(cond) {
          var m = cond.match(/^(>=|<=|!=|<>|>|<|=)?\\s*(.+)$/);
          if (!m) return true;
          var op = m[1] || '=';
          var val = parseFloat(m[2]);
          if (isNaN(val)) return true;
          return cols.some(function(cid) {
            var v = row.values[cid];
            if (v == null || isNaN(v)) return false;
            switch(op) {
              case '>=': return v >= val;
              case '<=': return v <= val;
              case '>':  return v > val;
              case '<':  return v < val;
              case '!=': case '<>': return v !== val;
              case '=':  return v === val;
              default:   return true;
            }
          });
        });
      });
    });
  }", cols_json))
}


# =============================================================================
# Helpers Protein_QUANT (quant_list_widget)
# =============================================================================

#' Load and validate a Protein_QUANT file
#' @param input Data frame or path to a TSV/CSV/Parquet file in Protein_QUANT format
#' @return Validated data frame
#' @noRd
.ql_load_quant <- function(input) {
  if (is.character(input) && length(input) == 1) {
    if (!file.exists(input)) stop("Protein_QUANT file not found: ", input)
    ext <- tolower(tools::file_ext(input))
    if (ext %in% c("tsv", "txt")) {
      if (requireNamespace("readr", quietly = TRUE)) {
        input <- readr::read_tsv(input, show_col_types = FALSE)
      } else {
        input <- utils::read.delim(input, sep = "\t", stringsAsFactors = FALSE)
      }
    } else if (ext == "csv") {
      if (requireNamespace("readr", quietly = TRUE)) {
        input <- readr::read_csv(input, show_col_types = FALSE)
      } else {
        input <- utils::read.csv(input, stringsAsFactors = FALSE)
      }
    } else if (ext == "parquet") {
      if (!requireNamespace("arrow", quietly = TRUE)) {
        stop("The 'arrow' package is required to read Parquet files")
      }
      input <- arrow::read_parquet(input)
    } else {
      stop("Unsupported format: ", ext, ". Use TSV, CSV or Parquet.")
    }
  }

  df <- as.data.frame(input)

  if (!"PG.ProteinGroups" %in% names(df)) {
    stop("Required column 'PG.ProteinGroups' not found. ",
         "The file must be a Protein_QUANT exported by preprocess_spectronaut().")
  }
  if (!any(grepl("^PG\\.NrOfPrecursorsUsedForQuantification_", names(df)))) {
    stop("No sample columns were detected (PG.NrOfPrecursorsUsedForQuantification_*). ",
         "Check that the file really is a Protein_QUANT.")
  }
  df
}


#' Load the normalized/imputed log2 matrix
#' @param input Data frame or path to a TSV/CSV/Parquet file with a ProteinGroups column + samples
#' @return Validated data frame
#' @noRd
.ql_load_matrix <- function(input) {
  if (is.character(input) && length(input) == 1) {
    if (!file.exists(input)) stop("Matrix file not found: ", input)
    ext <- tolower(tools::file_ext(input))
    if (ext %in% c("tsv", "txt")) {
      if (requireNamespace("readr", quietly = TRUE)) {
        input <- readr::read_tsv(input, show_col_types = FALSE)
      } else {
        input <- utils::read.delim(input, sep = "\t", stringsAsFactors = FALSE)
      }
    } else if (ext == "csv") {
      if (requireNamespace("readr", quietly = TRUE)) {
        input <- readr::read_csv(input, show_col_types = FALSE)
      } else {
        input <- utils::read.csv(input, stringsAsFactors = FALSE)
      }
    } else if (ext == "parquet") {
      if (!requireNamespace("arrow", quietly = TRUE)) {
        stop("The 'arrow' package is required to read Parquet files")
      }
      input <- arrow::read_parquet(input)
    } else {
      stop("Unsupported format: ", ext, ". Use TSV, CSV or Parquet.")
    }
  }

  df <- as.data.frame(input)

  if (!"ProteinGroups" %in% names(df)) {
    stop("Required column 'ProteinGroups' not found in the matrix. ",
         "Check that it is a valid matrix_log2_<...>.tsv.")
  }

  sample_cols <- setdiff(names(df), "ProteinGroups")
  ok <- grepl("^[A-Za-z0-9]+_\\d+$", sample_cols)
  if (!any(ok)) {
    stop("The matrix has no sample columns following the <cond>_<rep> pattern.")
  }
  df
}


#' Left-join the log2 matrix onto Protein_QUANT by Protein Groups
#' @param quant_df Protein_QUANT data frame
#' @param matrix_df log2 matrix data frame
#' @return quant_df with 16 extra `NormImp.Log2_<cond>_<rep>` columns
#' @noRd
.ql_join_matrix <- function(quant_df, matrix_df) {
  sample_cols <- setdiff(names(matrix_df), "ProteinGroups")
  ok <- grepl("^[A-Za-z0-9]+_\\d+$", sample_cols)
  sample_cols <- sample_cols[ok]

  idx <- match(quant_df$PG.ProteinGroups, matrix_df$ProteinGroups)
  for (sc in sample_cols) {
    new_name <- paste0("NormImp.Log2_", sc)
    quant_df[[new_name]] <- as.numeric(matrix_df[[sc]][idx])
  }
  quant_df
}


#' Parse Protein_QUANT column names into (metric, condition, replicate)
#' @param colnames_vec Vector of column names
#' @return Sorted data frame: column, metric, condition, replicate, coding
#' @noRd
.ql_parse_samples <- function(colnames_vec) {
  metrics <- c(
    "PG.NrOfPrecursorsUsedForQuantification",
    "PG.NrOfStrippedSequencesUsedForQuantification",
    "PG.Quantity",
    "NormImp.Log2"
  )

  rows <- list()
  for (m in metrics) {
    pat <- paste0("^", gsub("\\.", "\\\\.", m), "_(.+?)_(\\d+)$")
    hits <- regmatches(colnames_vec, regexec(pat, colnames_vec))
    for (i in seq_along(hits)) {
      h <- hits[[i]]
      if (length(h) == 3) {
        rows[[length(rows) + 1]] <- data.frame(
          column = colnames_vec[i],
          metric = m,
          condition = h[2],
          replicate = as.integer(h[3]),
          coding = paste(h[2], h[3], sep = "_"),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(rows) == 0) {
    return(data.frame(column = character(0), metric = character(0),
                      condition = character(0), replicate = integer(0),
                      coding = character(0), stringsAsFactors = FALSE))
  }

  df <- do.call(rbind, rows)
  df$metric <- factor(df$metric, levels = metrics)
  df <- df[order(df$metric, df$condition, df$replicate), ]
  df$metric <- as.character(df$metric)
  rownames(df) <- NULL
  df
}


#' Build the colDefs for the Protein_QUANT table
#' @param sample_map Output of .ql_parse_samples()
#' @param df Final data frame (used to detect which global columns are present)
#' @return List of colDef indexed by column name
#' @noRd
.ql_build_columns <- function(sample_map, df) {
  cols <- list()

  cols$PG.ProteinGroups <- colDef(
    name = "Protein Groups",
    minWidth = 160,
    sticky = "left",
    html = TRUE,
    cell = JS("function(cellInfo) {
      var val = cellInfo.value || '';
      var pids = val.split(';').map(function(s){return s.trim();}).filter(Boolean);
      var first = pids[0] || val;
      if (pids.length > 1) {
        return '<strong>' + first + '</strong> <span class=\"protein-count\">+' + (pids.length - 1) + '</span>';
      }
      return '<strong>' + first + '</strong>';
    }")
  )

  cols$PG.ProteinDescriptions <- colDef(
    name = "Descriptions",
    minWidth = 220,
    sticky = "left"
  )

  cols$PG.Genes <- colDef(
    name = "Gene Names",
    minWidth = 120,
    sticky = "left",
    html = TRUE,
    cell = JS("function(cellInfo) {
      var val = cellInfo.value || '';
      var genes = val.split(';').map(function(s){return s.trim();}).filter(Boolean);
      var first = genes[0] || val;
      if (genes.length > 1) {
        return first + ' <span class=\"protein-count\">+' + (genes.length - 1) + '</span>';
      }
      return first;
    }")
  )

  cols$PG.MolecularWeight <- colDef(
    name = "MW [kDa]",
    width = 95,
    align = "right",
    cell = JS("function(cellInfo) {
      var val = cellInfo.value;
      if (val == null || isNaN(val)) return '';
      return (val / 1000).toFixed(2);
    }")
  )

  # Global columns (single-col, plain values with an individual filter)
  global_specs <- list(
    list(id = "PG.NrOfPrecursorsIdentified.Global",        name = "# PSMs",       fmt = "Math.round(val)", width = 95),
    list(id = "PG.NrOfStrippedSequencesIdentified.Global", name = "# Uniq Pepts", fmt = "Math.round(val)", width = 115),
    list(id = "PG.Coverage.Global",                        name = "Coverage [%]", fmt = "val.toFixed(1)",  width = 115),
    list(id = "PG.Cscore",                                 name = "Cscore",       fmt = "val.toFixed(2)",  width = 95)
  )
  for (gs in global_specs) {
    if (!gs$id %in% names(df)) next
    cell_js <- JS(sprintf("function(cellInfo) {
      var val = cellInfo.value;
      if (val == null || isNaN(val)) return '<span class=\"pl-plain-value pl-bar-empty\">-</span>';
      return '<span class=\"pl-plain-value\">' + (%s) + '</span>';
    }", gs$fmt))
    cols[[gs$id]] <- colDef(
      name = gs$name,
      width = gs$width,
      align = "center",
      html = TRUE,
      filterable = TRUE,
      filterInput = .numeric_filter_hidden,
      filterMethod = .numeric_filter_method,
      cell = cell_js
    )
  }

  # Samples: 4 metrics x 16 conditions/reps
  first_cols_by_metric <- vapply(
    split(sample_map$column, sample_map$metric),
    function(x) x[1], character(1)
  )

  for (i in seq_len(nrow(sample_map))) {
    col_id    <- sample_map$column[i]
    cond      <- sample_map$condition[i]
    replicate <- sample_map$replicate[i]
    metric    <- sample_map$metric[i]

    fmt_js <- switch(
      metric,
      "PG.Quantity"  = "val.toExponential(2)",
      "NormImp.Log2" = "val.toFixed(2)",
      "Math.round(val)"
    )

    cell_js <- JS(sprintf("function(cellInfo) {
      var val = cellInfo.value;
      if (val == null || isNaN(val)) return '<span class=\"pl-plain-value pl-bar-empty\">-</span>';
      return '<span class=\"pl-plain-value\">' + (%s) + '</span>';
    }", fmt_js))

    is_anchor <- col_id %in% first_cols_by_metric
    cell_class <- if (is_anchor) "pl-sample-cell pl-group-boundary" else "pl-sample-cell"
    hdr_class  <- if (is_anchor) paste0("pl-hdr-", cond, " pl-group-boundary") else paste0("pl-hdr-", cond)
    fm <- if (is_anchor) {
      .pl_agg_filter_method(sample_map$column[sample_map$metric == metric])
    } else {
      .numeric_filter_method
    }

    col_width <- if (metric == "PG.Quantity") 88 else 72

    cols[[col_id]] <- colDef(
      name = paste(cond, replicate, sep = "_"),
      width = col_width,
      align = "center",
      html = TRUE,
      class = cell_class,
      headerClass = hdr_class,
      filterable = TRUE,
      filterInput = .numeric_filter_hidden,
      filterMethod = fm,
      cell = cell_js
    )
  }

  cols
}


#' Build the columnGroups for the Protein_QUANT table
#' @param sample_map Output of .ql_parse_samples()
#' @param static_cols Vector of the static columns that are present
#' @param global_cols Vector of the global columns that are present
#' @return List of colGroup
#' @noRd
.ql_build_column_groups <- function(sample_map, static_cols, global_cols) {
  groups <- list()

  if (length(static_cols) > 0) {
    groups[[length(groups) + 1]] <- colGroup(
      name = "PROTEIN ANNOTATION",
      columns = static_cols,
      sticky = "left"
    )
  }

  if (length(global_cols) > 0) {
    groups[[length(groups) + 1]] <- colGroup(
      name = "GLOBAL IDENTIFICATION DATA",
      columns = global_cols
    )
  }

  metric_labels <- list(
    "PG.NrOfPrecursorsUsedForQuantification"        = "# Uniq. PSMs Quantified",
    "PG.NrOfStrippedSequencesUsedForQuantification" = "# Uniq. Pepts Quantified",
    "PG.Quantity"                                   = "Raw Abundances",
    "NormImp.Log2"                                  = "Normalized and Imputed Abundances (Log2)"
  )
  for (m in names(metric_labels)) {
    cc <- sample_map$column[sample_map$metric == m]
    if (length(cc) == 0) next
    groups[[length(groups) + 1]] <- colGroup(
      name = metric_labels[[m]],
      columns = cc
    )
  }

  groups
}


# =============================================================================
# Helpers Summary-List (Metadata) widget
# =============================================================================

#' Load and validate the sample metadata
#' @param input Data frame or path to a TSV/CSV file
#' @return Validated data frame with the display columns
#' @noRd
.sl_load_metadata <- function(input) {
  if (is.character(input) && length(input) == 1) {
    if (!file.exists(input)) stop("File not found: ", input)
    ext <- tolower(tools::file_ext(input))
    if (ext %in% c("tsv", "txt")) {
      if (requireNamespace("readr", quietly = TRUE)) {
        df <- as.data.frame(readr::read_tsv(input, show_col_types = FALSE))
      } else {
        df <- read.delim(input, stringsAsFactors = FALSE, check.names = FALSE)
      }
    } else if (ext == "csv") {
      if (requireNamespace("readr", quietly = TRUE)) {
        df <- as.data.frame(readr::read_csv(input, show_col_types = FALSE))
      } else {
        df <- read.csv(input, stringsAsFactors = FALSE, check.names = FALSE)
      }
    } else {
      stop("Unsupported extension: ", ext)
    }
  } else if (is.data.frame(input)) {
    df <- as.data.frame(input)
  } else {
    stop("input must be a data.frame or a file path")
  }

  # Mandatory: present in every source (Spectronaut, LFQ, TMT).
  required <- c("R.FileName", "R.Condition", "R.Replicate", "Coding")
  missing_cols <- setdiff(required, names(df))
  if (length(missing_cols) > 0) {
    stop("Missing columns in the metadata: ", paste(missing_cols, collapse = ", "))
  }

  out <- data.frame(
    FileName         = as.character(df$R.FileName),
    Condition        = as.character(df$R.Condition),
    Replicate        = as.integer(df$R.Replicate),
    Coding           = as.character(df$Coding),
    check.names      = FALSE,
    stringsAsFactors = FALSE
  )

  # Optional count columns: they only exist in Spectronaut metadata. They are
  # absent for LFQ/TMT, and the widget degrades to showing just the base columns.
  count_map <- list(
    `# Unique PSMs`     = "R.PrecursorsIdentified",
    `# Unique Peptides` = "R.StrippedSequencesIdentified",
    `# Protein Groups`  = "R.ProteinGroupsIdentified"
  )
  for (disp in names(count_map)) {
    src <- count_map[[disp]]
    if (src %in% names(df)) out[[disp]] <- as.numeric(df[[src]])
  }

  out <- out[order(out$Condition, out$Replicate), , drop = FALSE]
  rownames(out) <- NULL
  out
}


#' Light tint of a hex colour (used to tint rows)
#' @param hex Base hex colour (e.g. "#4F81BD")
#' @param alpha Colour opacity over a white background (0-1). Default 0.12.
#' @return Tinted hex string "#RRGGBB"
#' @noRd
.sl_tint_color <- function(hex, alpha = 0.12) {
  if (is.null(hex) || is.na(hex) || !nzchar(hex)) return("#FFFFFF")
  h <- gsub("^#", "", hex)
  if (nchar(h) != 6) return("#FFFFFF")
  r <- strtoi(substr(h, 1, 2), 16L)
  g <- strtoi(substr(h, 3, 4), 16L)
  b <- strtoi(substr(h, 5, 6), 16L)
  tr <- round(r * alpha + 255 * (1 - alpha))
  tg <- round(g * alpha + 255 * (1 - alpha))
  tb <- round(b * alpha + 255 * (1 - alpha))
  sprintf("#%02X%02X%02X", tr, tg, tb)
}


#' Black reactable theme (Summary widget)
#' @return A reactableTheme object
#' @noRd
.sl_theme <- function() {
  reactableTheme(
    cellPadding = "8px 12px",
    highlightColor = "rgba(0, 0, 0, 0.06)",
    stripedColor = "rgba(0, 0, 0, 0.0)",
    rowSelectedStyle = list(
      backgroundColor = "rgba(0, 0, 0, 0.55)",
      color = "#ffffff",
      boxShadow = "inset 2px 0 0 0 #ffa62d"
    )
  )
}


#' Build the colDefs for the Summary table
#' @param df Metadata data frame (output of .sl_load_metadata())
#' @param palette Named vector: condition -> strong hex colour
#' @return List of colDef
#' @noRd
.sl_build_columns <- function(df, palette) {
  palette_json <- jsonlite::toJSON(as.list(palette), auto_unbox = TRUE)

  chip_render_condition <- reactable::JS(sprintf("
    function(cellInfo) {
      var pal = %s;
      var cond = String(cellInfo.value);
      var color = pal[cond] || '#6c757d';
      return '<span class=\"sl-chip\" style=\"background-color:' + color + '\">' +
             cond + '</span>';
    }
  ", palette_json))

  chip_render_coding <- reactable::JS(sprintf("
    function(cellInfo) {
      var pal = %s;
      var coding = String(cellInfo.value);
      var cond = (coding.split('_')[0]) || coding;
      var color = pal[cond] || '#6c757d';
      return '<span class=\"sl-chip\" style=\"background-color:' + color + '\">' +
             coding + '</span>';
    }
  ", palette_json))

  number_render <- reactable::JS("
    function(cellInfo) {
      var v = cellInfo.value;
      if (v == null || isNaN(v)) return '';
      return Math.round(v).toLocaleString('en-US');
    }
  ")

  cols <- list()

  cols$FileName <- colDef(
    name     = "FileName",
    sticky   = "left",
    minWidth = 360,
    align    = "left",
    vAlign   = "center",
    headerClass = "sl-filename-header",
    headerStyle = list(
      background     = "#1a1a1a",
      color          = "#ffffff",
      height         = "45px",
      display        = "flex",
      alignItems     = "center"
    ),
    style    = list(fontSize   = "14.5px",
                    fontWeight = "500",
                    display    = "flex",
                    alignItems = "center"),
    filterable = TRUE
  )

  cols$Condition <- colDef(
    name     = "Condition",
    minWidth = 105,
    align    = "center",
    html     = TRUE,
    cell     = chip_render_condition,
    filterable = FALSE
  )

  cols$Replicate <- colDef(
    name     = "Replicate",
    minWidth = 90,
    align    = "center",
    cell     = reactable::JS("function(c) { return (c.value == null) ? '' : String(c.value); }"),
    filterable = FALSE
  )

  cols$Coding <- colDef(
    name     = "Coding",
    minWidth = 105,
    align    = "center",
    html     = TRUE,
    cell     = chip_render_coding,
    filterable = FALSE
  )

  # Only the count columns that are present (absent from LFQ/TMT metadata)
  for (nm in intersect(c("# Unique PSMs", "# Unique Peptides", "# Protein Groups"),
                       names(df))) {
    cols[[nm]] <- colDef(
      name        = nm,
      minWidth    = 145,
      align       = "center",
      cell        = number_render,
      filterable  = TRUE,
      filterInput = .numeric_filter_hidden,
      filterMethod = .numeric_filter_method
    )
  }

  cols
}


#' rowStyle JS function that tints the rows by condition
#' @return A JS object for reactable::rowStyle
#' @noRd
.sl_build_row_style <- function() {
  reactable::JS("
    function(rowInfo) {
      if (!rowInfo || !rowInfo.values) return null;
      var cond = rowInfo.values['Condition'];
      var pal = (typeof SL_TINT_PALETTE !== 'undefined') ? SL_TINT_PALETTE : null;
      if (!pal) return null;
      var t = pal[cond];
      return t ? { backgroundColor: t } : null;
    }
  ")
}


# =============================================================================
# Main function -- Protein_ID
# =============================================================================

#' Interactive Reactable Table for Protein_ID Data (post-Spectronaut)
#'
#' Builds a reactable table with 2-level grouped headers, static columns pinned
#' (sticky) to the left and data bars coloured by condition inside each numeric
#' cell. It reproduces the usual Excel layout and improves on it.
#'
#' @param data Data frame or path to a Protein_ID TSV/CSV/Parquet file (output of
#'   preprocess_spectronaut() -> protein_id).
#' @param metadata Optional: metadata data frame (run_summary) with a Coding
#'   column that fixes the order of the sample columns.
#' @param page_size Rows per page (default 15).
#' @param height Table height in px (default 720).
#' @param element_id Reactable element ID (default NULL).
#' @param selection "single" | "multiple" | NULL.
#' @param searchable Enable the search box (default TRUE).
#'
#' @return A reactable object.
#'
#' @examples
#' data(nadia_dia)
#'
#' # metadata$Coding fixes the order of the sample columns
#' tbl <- protein_list_reactable(nadia_dia$protein_id,
#'                               metadata = nadia_dia$metadata)
#' class(tbl)
#'
#' @export
protein_list_reactable <- function(
    data,
    metadata = NULL,
    page_size = 15,
    height = 720,
    element_id = NULL,
    selection = NULL,
    searchable = TRUE
) {
  df <- .pl_load_data(data)
  sample_map <- .pl_parse_samples(names(df))

  if (nrow(sample_map) == 0) {
    stop("No sample columns were found in the data frame. ",
         "Expected format: PG.<metric>_<condition>_<replicate>")
  }

  if (!is.null(metadata) && is.data.frame(metadata) && "Coding" %in% names(metadata)) {
    coding_order <- as.character(metadata$Coding)
    sample_map$.idx <- match(sample_map$coding, coding_order)
    metric_levels <- c("PG.NrOfPrecursorsIdentified",
                       "PG.NrOfStrippedSequencesIdentified",
                       "PG.Coverage",
                       "PG.Cscore.RunWise")
    sample_map <- sample_map[order(factor(sample_map$metric, levels = metric_levels),
                                   sample_map$.idx), ]
    sample_map$.idx <- NULL
    rownames(sample_map) <- NULL
  }

  max_per_col <- vapply(sample_map$column, function(cc) {
    suppressWarnings(max(as.numeric(df[[cc]]), na.rm = TRUE))
  }, numeric(1))
  max_per_col[!is.finite(max_per_col)] <- 1
  names(max_per_col) <- sample_map$column

  static_cols <- intersect(
    c("PG.ProteinGroups", "PG.ProteinDescriptions", "PG.Genes", "PG.MolecularWeight"),
    names(df)
  )
  ordered_cols <- c(static_cols, sample_map$column)
  df <- df[, ordered_cols, drop = FALSE]

  cols <- .pl_build_columns(sample_map, max_per_col)
  groups <- .pl_build_column_groups(sample_map, static_cols)

  reactable(
    df,
    elementId     = element_id,
    defaultSorted = list(PG.ProteinGroups = "asc"),
    defaultPageSize    = page_size,
    showPageSizeOptions = TRUE,
    pageSizeOptions = c(15, 30, 50, 100),
    resizable   = TRUE,
    selection   = selection,
    onClick     = if (!is.null(selection)) "select" else NULL,
    defaultColDef = colDef(
      align = "left",
      headerStyle = list(
        background  = "rgba(31, 78, 120, 0.9)",
        color       = "#ffffff",
        height      = "45px",
        display     = "flex",
        alignItems  = "center",
        justifyContent = "center"
      ),
      style = list(height = "48px", display = "flex", alignItems = "center")
    ),
    columns      = cols,
    columnGroups = groups,
    wrap         = FALSE,
    class        = "rl-table pl-table",
    rowStyle     = if (!is.null(selection)) list(cursor = "pointer") else NULL,
    highlight    = TRUE,
    searchable   = searchable,
    height       = height,
    striped      = TRUE,
    theme        = .pl_theme(),
    language     = .rl_lang(),
    details      = .pl_detail_row()
  )
}


# =============================================================================
# Standalone / Quarto wrapper for Protein_ID
# =============================================================================

#' Full Widget for the Protein_ID Table with Filters, Search and Excel Export
#'
#' Wraps \code{protein_list_reactable()} with:
#'  - A search bar and buttons (toggle filters, clear, export to Excel).
#'  - Per-condition chips that hide/show the 16 associated columns.
#'  - A numeric filter on MW in Da with AND/OR operators.
#'  - Excel export with 2-level headers (merged groups + sub-labels) and a
#'    per-condition tint, reproducing the manual layout.
#'
#' @inheritParams protein_list_reactable
#' @param element_id Element ID (default "protein_id_table").
#'
#' @return A browsable htmltools object.
#'
#' @examples
#' if (requireNamespace("jsonlite", quietly = TRUE)) {
#'   data(nadia_dia)
#'
#'   w <- protein_list_widget(nadia_dia$protein_id,
#'                            metadata = nadia_dia$metadata)
#'   print(class(w))   # print(w) itself opens it in the Viewer / browser
#' }
#'
#' @export
protein_list_widget <- function(
    data,
    metadata = NULL,
    page_size = 15,
    height = 720,
    element_id = "protein_id_table",
    selection = NULL,
    searchable = TRUE
) {

  df <- .pl_load_data(data)
  sample_map <- .pl_parse_samples(names(df))

  if (nrow(sample_map) == 0) {
    stop("No sample columns were found in the data frame.")
  }

  if (!is.null(metadata) && is.data.frame(metadata) && "Coding" %in% names(metadata)) {
    coding_order <- as.character(metadata$Coding)
    sample_map$.idx <- match(sample_map$coding, coding_order)
    metric_levels <- c("PG.NrOfPrecursorsIdentified",
                       "PG.NrOfStrippedSequencesIdentified",
                       "PG.Coverage",
                       "PG.Cscore.RunWise")
    sample_map <- sample_map[order(factor(sample_map$metric, levels = metric_levels),
                                   sample_map$.idx), ]
    sample_map$.idx <- NULL
    rownames(sample_map) <- NULL
  }

  max_per_col <- vapply(sample_map$column, function(cc) {
    suppressWarnings(max(as.numeric(df[[cc]]), na.rm = TRUE))
  }, numeric(1))
  max_per_col[!is.finite(max_per_col)] <- 1
  names(max_per_col) <- sample_map$column

  static_cols <- intersect(
    c("PG.ProteinGroups", "PG.ProteinDescriptions", "PG.Genes", "PG.MolecularWeight"),
    names(df)
  )
  ordered_cols <- c(static_cols, sample_map$column)
  df <- df[, ordered_cols, drop = FALSE]

  cols   <- .pl_build_columns(sample_map, max_per_col)
  groups <- .pl_build_column_groups(sample_map, static_cols)

  conditions <- unique(sample_map$condition)
  palette    <- .pl_condition_palette(conditions)

  # --- JSON structures for the client side ---
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("The 'jsonlite' package is required by protein_list_widget().")
  }

  cond_struct <- setNames(
    lapply(conditions, function(cc) sample_map$column[sample_map$condition == cc]),
    conditions
  )
  cond_struct_json <- jsonlite::toJSON(cond_struct, auto_unbox = FALSE)

  static_headers_map <- c(
    PG.ProteinGroups       = "Protein Groups",
    PG.ProteinDescriptions = "Descriptions",
    PG.Genes               = "Gene Names",
    PG.MolecularWeight     = "MW [kDa]"
  )
  group_struct <- list()
  if (length(static_cols) > 0) {
    group_struct[[length(group_struct) + 1]] <- list(
      name    = jsonlite::unbox("PROTEIN ANNOTATION"),
      columns = static_cols,
      headers = unname(static_headers_map[static_cols])
    )
  }
  metric_labels <- c(
    "PG.NrOfPrecursorsIdentified"        = "# Uniq. PSMs Identified",
    "PG.NrOfStrippedSequencesIdentified" = "# Uniq. Pepts Identified",
    "PG.Coverage"                        = "Coverage [%]",
    "PG.Cscore.RunWise"                  = "Spectronaut Cscore"
  )
  for (m in names(metric_labels)) {
    cc <- sample_map$column[sample_map$metric == m]
    if (length(cc) == 0) next
    hh <- sample_map$coding[sample_map$metric == m]
    group_struct[[length(group_struct) + 1]] <- list(
      name    = jsonlite::unbox(metric_labels[[m]]),
      columns = cc,
      headers = hh
    )
  }
  group_struct_json <- jsonlite::toJSON(group_struct)

  palette_json <- jsonlite::toJSON(as.list(palette), auto_unbox = TRUE)

  # --- CSS + CDN scripts ---
  css <- .rl_css()
  # --- ExcelJS and PapaParse, shipped with the package (see .rl_export_deps) ---
  export_deps <- .rl_export_deps()

  # --- JavaScript (toggle filters, condition chips, clear, export) ---
  js_code <- tags$script(HTML(sprintf("
    var plFiltersVisible = false;
    var plHiddenConditions = {};
    var plCondStruct = %s;
    var plGroupStruct = %s;
    var plPalette = %s;

    function plToggleFilters() {
      var container = document.querySelector('.pl-filters-container');
      var btn = document.querySelector('.pl-btn-toggle-filters');
      if (plFiltersVisible) {
        container.style.maxHeight = '0';
        container.style.opacity = '0';
        container.style.marginBottom = '0';
        container.style.padding = '0';
        container.style.borderWidth = '0';
        container.style.boxShadow = 'none';
        container.style.overflow = 'hidden';
        btn.classList.add('filters-hidden');
      } else {
        container.style.maxHeight = '500px';
        container.style.opacity = '1';
        container.style.marginBottom = '1rem';
        container.style.padding = '1.5rem';
        container.style.borderWidth = '1px';
        container.style.boxShadow = '0 2px 4px rgba(0, 0, 0, 0.05)';
        container.style.overflow = 'visible';
        btn.classList.remove('filters-hidden');
      }
      plFiltersVisible = !plFiltersVisible;
    }

    function plApplyNumericFilter(columnId, value) {
      Reactable.setFilter('%s', columnId, value || undefined);
    }

    function plToggleCondition(cond) {
      plHiddenConditions[cond] = !plHiddenConditions[cond];
      var chip = document.querySelector('.pl-cond-chip[data-cond=\"' + cond + '\"]');
      if (chip) chip.classList.toggle('off', plHiddenConditions[cond]);
      var hiddenCols = [];
      Object.keys(plHiddenConditions).forEach(function(c) {
        if (plHiddenConditions[c] && plCondStruct[c]) {
          hiddenCols = hiddenCols.concat(plCondStruct[c]);
        }
      });
      Reactable.setHiddenColumns('%s', hiddenCols);
    }

    function plClearFilters() {
      var numInputs = document.querySelectorAll('.pl-filters-container .rl-numeric-input');
      numInputs.forEach(function(inp) {
        inp.value = '';
        var col = inp.getAttribute('data-column');
        if (col) Reactable.setFilter('%s', col, undefined);
      });
      Object.keys(plHiddenConditions).forEach(function(c) {
        plHiddenConditions[c] = false;
        var chip = document.querySelector('.pl-cond-chip[data-cond=\"' + c + '\"]');
        if (chip) chip.classList.remove('off');
      });
      Reactable.setHiddenColumns('%s', []);
      var searchInput = document.querySelector('.pl-search-input');
      if (searchInput) {
        searchInput.value = '';
        Reactable.setSearch('%s', '');
      }
    }

    function plHex2Argb(hex) {
      var h = (hex || '').replace('#', '');
      if (h.length !== 6) return 'FF808080';
      return 'FF' + h.toUpperCase();
    }

    async function plExportExcel() {
      try {
        var tsv = Reactable.getDataCSV('%s', { sep: '\\t' });
        var parseResult = Papa.parse(tsv, { header: true, delimiter: '\\t', skipEmptyLines: true });
        var rows = parseResult.data;

        var wb = new ExcelJS.Workbook();
        var ws = wb.addWorksheet('Protein-List_ID');

        var flatCols = [];
        var flatHeaders = [];
        var colCondMap = {};
        plGroupStruct.forEach(function(g) {
          g.columns.forEach(function(c, i) {
            flatCols.push(c);
            flatHeaders.push(g.headers[i]);
            var mm = (g.headers[i] || '').match(/^([A-Za-z0-9]+)_[0-9]+$/);
            if (mm) colCondMap[c] = mm[1];
          });
        });

        flatCols.forEach(function(c, i) {
          var width;
          if (c === 'PG.ProteinGroups' || c === 'PG.Genes') width = 22;
          else if (c === 'PG.ProteinDescriptions') width = 42;
          else if (c === 'PG.MolecularWeight') width = 12;
          else width = 8;
          ws.getColumn(i + 1).width = width;
        });

        // Row 1: group headers (merged)
        var colOffset = 1;
        plGroupStruct.forEach(function(g) {
          if (!g.columns || g.columns.length === 0) return;
          var startCol = colOffset;
          var endCol = colOffset + g.columns.length - 1;
          var cell = ws.getCell(1, startCol);
          cell.value = g.name;
          if (endCol > startCol) ws.mergeCells(1, startCol, 1, endCol);
          cell.font = { bold: true, color: { argb: 'FFFFFFFF' }, size: 12 };
          cell.fill = { type: 'pattern', pattern: 'solid', fgColor: { argb: 'FF1F5078' } };
          cell.alignment = { vertical: 'middle', horizontal: 'center' };
          colOffset = endCol + 1;
        });
        ws.getRow(1).height = 28;

        // Row 2: per-column headers tinted by condition
        flatCols.forEach(function(c, i) {
          var cell = ws.getCell(2, i + 1);
          cell.value = flatHeaders[i];
          var cond = colCondMap[c];
          var fillColor = cond && plPalette[cond] ? plHex2Argb(plPalette[cond]) : 'FF595959';
          cell.fill = { type: 'pattern', pattern: 'solid', fgColor: { argb: fillColor } };
          cell.font = { bold: true, color: { argb: 'FFFFFFFF' }, size: 11 };
          cell.alignment = { vertical: 'middle', horizontal: 'center' };
        });
        ws.getRow(2).height = 24;

        // Data rows starting at 3
        var nextRow = 3;
        rows.forEach(function(row) {
          var excelRow = ws.getRow(nextRow);
          flatCols.forEach(function(c, i) {
            var cell = excelRow.getCell(i + 1);
            var v = row[c];
            if (v != null && v !== '') {
              var num = parseFloat(v);
              cell.value = !isNaN(num) ? num : v;
            }
          });
          excelRow.commit();
          nextRow++;
        });

        // Thin borders across used range
        for (var r = 1; r < nextRow; r++) {
          var row = ws.getRow(r);
          for (var c = 1; c <= flatCols.length; c++) {
            row.getCell(c).border = {
              top: { style: 'thin' }, left: { style: 'thin' },
              bottom: { style: 'thin' }, right: { style: 'thin' }
            };
          }
        }

        // Freeze header rows + first 4 cols
        ws.views = [{ state: 'frozen', xSplit: Math.min(4, flatCols.length), ySplit: 2 }];

        var buffer = await wb.xlsx.writeBuffer();
        var blob = new Blob([buffer], { type: 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet' });
        var url = window.URL.createObjectURL(blob);
        var a = document.createElement('a');
        a.href = url;
        a.download = 'Protein_List_ID_' + new Date().toISOString().split('T')[0] + '.xlsx';
        a.click();
        window.URL.revokeObjectURL(url);
      } catch(e) {
        console.error('Error while exporting:', e);
      }
    }
  ", cond_struct_json, group_struct_json, palette_json,
      element_id, element_id, element_id, element_id, element_id, element_id)))

  # --- Search bar + buttons ---
  search_actions <- div(class = "rl-search-actions",
    tags$input(
      type = "search",
      placeholder = "Search proteins, genes, descriptions...",
      class = "rl-search-input pl-search-input",
      oninput = sprintf("Reactable.setSearch('%s', this.value)", element_id)
    ),
    div(class = "rl-action-buttons",
      tags$button(
        class = "rl-btn-action pl-btn-toggle-filters filters-hidden",
        onclick = "plToggleFilters()",
        title = "Show/Hide filters",
        HTML('<svg xmlns="http://www.w3.org/2000/svg" width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><polygon points="22 3 2 3 10 12.46 10 19 14 21 14 12.46 22 3"></polygon></svg>')
      ),
      tags$button(
        class = "rl-btn-action",
        onclick = "plClearFilters()",
        title = "Clear filters",
        HTML('<svg xmlns="http://www.w3.org/2000/svg" width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M3 6h18"></path><path d="M8 6V4c0-1 1-2 2-2h4c1 0 2 1 2 2v2"></path><path d="M19 6l-1 14c0 1-1 2-2 2H8c-1 0-2-1-2-2L5 6"></path><line x1="1" y1="1" x2="23" y2="23" stroke="#E63946" stroke-width="2"></line></svg>')
      ),
      tags$button(
        class = "rl-btn-action",
        onclick = "plExportExcel()",
        title = "Export to Excel",
        HTML('<svg xmlns="http://www.w3.org/2000/svg" width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4"></path><polyline points="7 10 12 15 17 10"></polyline><line x1="12" y1="15" x2="12" y2="3"></line></svg>')
      )
    )
  )

  # --- Filter panel ---
  cond_chips <- lapply(conditions, function(cc) {
    tags$span(
      class = "pl-cond-chip",
      `data-cond` = cc,
      style = sprintf("background-color:%s;", palette[[cc]]),
      onclick = sprintf("plToggleCondition('%s')", cc),
      cc
    )
  })

  # --- Aggregated filter anchors (first column of each metric) ---
  .first_of <- function(mm) {
    cc <- sample_map$column[sample_map$metric == mm]
    if (length(cc) == 0) NA_character_ else cc[1]
  }
  metric_filters <- list(
    list(label = "# Uniq. PSMs Identified",  anchor = .first_of("PG.NrOfPrecursorsIdentified"),        placeholder = "\u2265 2 ..."),
    list(label = "# Uniq. Pepts Identified", anchor = .first_of("PG.NrOfStrippedSequencesIdentified"), placeholder = "\u2265 2 ..."),
    list(label = "Coverage [%]",             anchor = .first_of("PG.Coverage"),                        placeholder = "\u2265 10 ..."),
    list(label = "Spectronaut Cscore",       anchor = .first_of("PG.Cscore.RunWise"),                  placeholder = "\u2265 2 ...")
  )
  metric_filter_items <- lapply(Filter(function(f) !is.na(f$anchor), metric_filters), function(f) {
    div(class = "rl-filter-item",
      tags$label(class = "rl-filter-label", f$label),
      tags$input(
        type = "text",
        class = "rl-numeric-input",
        `data-column` = f$anchor,
        placeholder = f$placeholder,
        oninput = sprintf("plApplyNumericFilter('%s', this.value)", f$anchor)
      ),
      tags$span(class = "rl-filter-hint", "OR across the samples in the block")
    )
  })

  filters_panel <- div(class = "rl-filters-container pl-filters-container",
    div(class = "rl-filters-row",
      div(class = "rl-filter-item", style = "flex: 2 1 300px;",
        tags$label(class = "rl-filter-label", "Conditions (click to hide/show 16 cols)"),
        div(class = "pl-cond-chips", cond_chips)
      ),
      metric_filter_items
    )
  )

  # --- Table ---
  tbl <- reactable(
    df,
    elementId     = element_id,
    defaultSorted = list(PG.ProteinGroups = "asc"),
    defaultPageSize    = page_size,
    showPageSizeOptions = TRUE,
    pageSizeOptions = c(15, 30, 50, 100),
    resizable   = TRUE,
    selection   = selection,
    onClick     = if (!is.null(selection)) "select" else NULL,
    defaultColDef = colDef(
      align = "left",
      headerStyle = list(
        background  = "rgba(31, 78, 120, 0.9)",
        color       = "#ffffff",
        height      = "45px",
        display     = "flex",
        alignItems  = "center",
        justifyContent = "center"
      ),
      style = list(height = "48px", display = "flex", alignItems = "center")
    ),
    columns      = cols,
    columnGroups = groups,
    wrap         = FALSE,
    class        = "rl-table pl-table",
    rowStyle     = if (!is.null(selection)) list(cursor = "pointer") else NULL,
    highlight    = TRUE,
    searchable   = searchable,
    height       = height,
    striped      = TRUE,
    theme        = .pl_theme(),
    language     = .rl_lang(),
    details      = .pl_detail_row()
  )

  browsable(tagList(
    css,
    export_deps,
    js_code,
    search_actions,
    filters_panel,
    tbl
  ))
}


# =============================================================================
# Main function -- Protein_QUANT (quant_list_widget)
# =============================================================================

#' Full Widget for the Protein_QUANT Table with the log2 Matrix Appended
#'
#' Reactable table that reproduces the layout of the Protein-List_QUANT Excel
#' sheet, with 2-level headers, sticky columns, a blue palette, aggregated
#' filters and Excel export. It appends to Protein_QUANT the 16 columns of the
#' normalized/imputed log2 matrix, matching by Protein Groups.
#'
#' @param data Data frame or path to a Protein_QUANT TSV/CSV/Parquet file.
#' @param matrix_data Data frame or path to a TSV/CSV/Parquet file with the
#'   normalized/imputed log2 matrix (a ProteinGroups column + <cond>_<rep> samples).
#' @param metadata Optional: metadata data frame with a Coding column that fixes
#'   the sample order.
#' @param page_size Rows per page (default 15).
#' @param height Table height in px (default 720).
#' @param element_id Reactable element ID.
#' @param selection "single" | "multiple" | NULL.
#' @param searchable Enable the search box (default TRUE).
#'
#' @return A browsable object.
#'
#' @examples
#' if (requireNamespace("jsonlite", quietly = TRUE)) {
#'   data(nadia_dia)
#'   res <- process_proteomics(nadia_dia, verbose = FALSE)
#'
#'   # The log2 matrix: a ProteinGroups column plus one column per sample
#'   mat <- SummarizedExperiment::assay(res$se_proc, "Impseqrob_min")
#'   mat_df <- data.frame(ProteinGroups = rownames(mat), mat,
#'                        check.names = FALSE)
#'
#'   w <- quant_list_widget(nadia_dia$protein_quant, matrix_data = mat_df,
#'                          metadata = nadia_dia$metadata)
#'   print(class(w))
#' }
#'
#' @export
quant_list_widget <- function(
    data        = "data-raw/Protein_QUANT_20260423_142504.tsv",
    matrix_data = "data-raw/matrix_log2_cycloess_Impseq_min.tsv",
    metadata    = NULL,
    page_size   = 15,
    height      = 720,
    element_id  = "quant_id_table",
    selection   = NULL,
    searchable  = TRUE
) {

  df_quant  <- .ql_load_quant(data)
  df_matrix <- .ql_load_matrix(matrix_data)
  df <- .ql_join_matrix(df_quant, df_matrix)

  sample_map <- .ql_parse_samples(names(df))

  if (nrow(sample_map) == 0) {
    stop("No sample columns were found in the combined data frame.")
  }

  if (!is.null(metadata) && is.data.frame(metadata) && "Coding" %in% names(metadata)) {
    coding_order <- as.character(metadata$Coding)
    sample_map$.idx <- match(sample_map$coding, coding_order)
    metric_levels <- c("PG.NrOfPrecursorsUsedForQuantification",
                       "PG.NrOfStrippedSequencesUsedForQuantification",
                       "PG.Quantity",
                       "NormImp.Log2")
    sample_map <- sample_map[order(factor(sample_map$metric, levels = metric_levels),
                                   sample_map$.idx), ]
    sample_map$.idx <- NULL
    rownames(sample_map) <- NULL
  }

  static_cols <- intersect(
    c("PG.ProteinGroups", "PG.ProteinDescriptions", "PG.Genes", "PG.MolecularWeight"),
    names(df)
  )
  global_cols <- intersect(
    c("PG.NrOfPrecursorsIdentified.Global",
      "PG.NrOfStrippedSequencesIdentified.Global",
      "PG.Coverage.Global",
      "PG.Cscore"),
    names(df)
  )

  ordered_cols <- c(static_cols, global_cols, sample_map$column)
  df <- df[, ordered_cols, drop = FALSE]

  cols   <- .ql_build_columns(sample_map, df)
  groups <- .ql_build_column_groups(sample_map, static_cols, global_cols)

  conditions <- unique(sample_map$condition)
  palette    <- .pl_condition_palette(conditions)

  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("The 'jsonlite' package is required by quant_list_widget().")
  }

  cond_struct <- setNames(
    lapply(conditions, function(cc) sample_map$column[sample_map$condition == cc]),
    conditions
  )
  cond_struct_json <- jsonlite::toJSON(cond_struct, auto_unbox = FALSE)

  static_headers_map <- c(
    PG.ProteinGroups       = "Protein Groups",
    PG.ProteinDescriptions = "Descriptions",
    PG.Genes               = "Gene Names",
    PG.MolecularWeight     = "MW [kDa]"
  )
  global_headers_map <- c(
    PG.NrOfPrecursorsIdentified.Global        = "# PSMs",
    PG.NrOfStrippedSequencesIdentified.Global = "# Uniq Pepts",
    PG.Coverage.Global                        = "Coverage [%]",
    PG.Cscore                                 = "Cscore"
  )

  group_struct <- list()
  if (length(static_cols) > 0) {
    group_struct[[length(group_struct) + 1]] <- list(
      name    = jsonlite::unbox("PROTEIN ANNOTATION"),
      columns = static_cols,
      headers = unname(static_headers_map[static_cols])
    )
  }
  if (length(global_cols) > 0) {
    group_struct[[length(group_struct) + 1]] <- list(
      name    = jsonlite::unbox("GLOBAL IDENTIFICATION DATA"),
      columns = global_cols,
      headers = unname(global_headers_map[global_cols])
    )
  }
  metric_labels <- c(
    "PG.NrOfPrecursorsUsedForQuantification"        = "# Uniq. PSMs Quantified",
    "PG.NrOfStrippedSequencesUsedForQuantification" = "# Uniq. Pepts Quantified",
    "PG.Quantity"                                   = "Raw Abundances",
    "NormImp.Log2"                                  = "Normalized and Imputed Abundances (Log2)"
  )
  for (m in names(metric_labels)) {
    cc <- sample_map$column[sample_map$metric == m]
    if (length(cc) == 0) next
    hh <- sample_map$coding[sample_map$metric == m]
    group_struct[[length(group_struct) + 1]] <- list(
      name    = jsonlite::unbox(metric_labels[[m]]),
      columns = cc,
      headers = hh
    )
  }
  group_struct_json <- jsonlite::toJSON(group_struct)

  palette_json <- jsonlite::toJSON(as.list(palette), auto_unbox = TRUE)

  css <- .rl_css()
  # --- ExcelJS and PapaParse, shipped with the package (see .rl_export_deps) ---
  export_deps <- .rl_export_deps()

  js_code <- tags$script(HTML(sprintf("
    var qlFiltersVisible = false;
    var qlHiddenConditions = {};
    var qlCondStruct = %s;
    var qlGroupStruct = %s;
    var qlPalette = %s;

    function qlToggleFilters() {
      var container = document.querySelector('.ql-filters-container');
      var btn = document.querySelector('.ql-btn-toggle-filters');
      if (qlFiltersVisible) {
        container.style.maxHeight = '0';
        container.style.opacity = '0';
        container.style.marginBottom = '0';
        container.style.padding = '0';
        container.style.borderWidth = '0';
        container.style.boxShadow = 'none';
        container.style.overflow = 'hidden';
        btn.classList.add('filters-hidden');
      } else {
        container.style.maxHeight = '500px';
        container.style.opacity = '1';
        container.style.marginBottom = '1rem';
        container.style.padding = '1.5rem';
        container.style.borderWidth = '1px';
        container.style.boxShadow = '0 2px 4px rgba(0, 0, 0, 0.05)';
        container.style.overflow = 'visible';
        btn.classList.remove('filters-hidden');
      }
      qlFiltersVisible = !qlFiltersVisible;
    }

    function qlApplyNumericFilter(columnId, value) {
      Reactable.setFilter('%s', columnId, value || undefined);
    }

    function qlToggleCondition(cond) {
      qlHiddenConditions[cond] = !qlHiddenConditions[cond];
      var chip = document.querySelector('.ql-cond-chip[data-cond=\"' + cond + '\"]');
      if (chip) chip.classList.toggle('off', qlHiddenConditions[cond]);
      var hiddenCols = [];
      Object.keys(qlHiddenConditions).forEach(function(c) {
        if (qlHiddenConditions[c] && qlCondStruct[c]) {
          hiddenCols = hiddenCols.concat(qlCondStruct[c]);
        }
      });
      Reactable.setHiddenColumns('%s', hiddenCols);
    }

    function qlClearFilters() {
      var numInputs = document.querySelectorAll('.ql-filters-container .rl-numeric-input');
      numInputs.forEach(function(inp) {
        inp.value = '';
        var col = inp.getAttribute('data-column');
        if (col) Reactable.setFilter('%s', col, undefined);
      });
      Object.keys(qlHiddenConditions).forEach(function(c) {
        qlHiddenConditions[c] = false;
        var chip = document.querySelector('.ql-cond-chip[data-cond=\"' + c + '\"]');
        if (chip) chip.classList.remove('off');
      });
      Reactable.setHiddenColumns('%s', []);
      var searchInput = document.querySelector('.ql-search-input');
      if (searchInput) {
        searchInput.value = '';
        Reactable.setSearch('%s', '');
      }
    }

    function qlHex2Argb(hex) {
      var h = (hex || '').replace('#', '');
      if (h.length !== 6) return 'FF808080';
      return 'FF' + h.toUpperCase();
    }

    async function qlExportExcel() {
      try {
        var tsv = Reactable.getDataCSV('%s', { sep: '\\t' });
        var parseResult = Papa.parse(tsv, { header: true, delimiter: '\\t', skipEmptyLines: true });
        var rows = parseResult.data;

        var wb = new ExcelJS.Workbook();
        var ws = wb.addWorksheet('Protein-List_QUANT');

        var flatCols = [];
        var flatHeaders = [];
        var colCondMap = {};
        qlGroupStruct.forEach(function(g) {
          g.columns.forEach(function(c, i) {
            flatCols.push(c);
            flatHeaders.push(g.headers[i]);
            var mm = (g.headers[i] || '').match(/^([A-Za-z0-9]+)_[0-9]+$/);
            if (mm) colCondMap[c] = mm[1];
          });
        });

        flatCols.forEach(function(c, i) {
          var width;
          if (c === 'PG.ProteinGroups' || c === 'PG.Genes') width = 22;
          else if (c === 'PG.ProteinDescriptions') width = 42;
          else if (c === 'PG.MolecularWeight') width = 12;
          else if (c.indexOf('.Global') !== -1 || c === 'PG.Cscore') width = 12;
          else if (c.indexOf('PG.Quantity_') === 0) width = 14;
          else width = 9;
          ws.getColumn(i + 1).width = width;
        });

        // Row 1: group headers (merged) - blue
        var colOffset = 1;
        qlGroupStruct.forEach(function(g) {
          if (!g.columns || g.columns.length === 0) return;
          var startCol = colOffset;
          var endCol = colOffset + g.columns.length - 1;
          var cell = ws.getCell(1, startCol);
          cell.value = g.name;
          if (endCol > startCol) ws.mergeCells(1, startCol, 1, endCol);
          cell.font = { bold: true, color: { argb: 'FFFFFFFF' }, size: 12 };
          cell.fill = { type: 'pattern', pattern: 'solid', fgColor: { argb: 'FF660505' } };
          cell.alignment = { vertical: 'middle', horizontal: 'center' };
          colOffset = endCol + 1;
        });
        ws.getRow(1).height = 28;

        // Row 2: per-column headers; condition tint for samples, plain for static + global
        flatCols.forEach(function(c, i) {
          var cell = ws.getCell(2, i + 1);
          cell.value = flatHeaders[i];
          var cond = colCondMap[c];
          var fillColor = cond && qlPalette[cond] ? qlHex2Argb(qlPalette[cond]) : 'FF595959';
          cell.fill = { type: 'pattern', pattern: 'solid', fgColor: { argb: fillColor } };
          cell.font = { bold: true, color: { argb: 'FFFFFFFF' }, size: 11 };
          cell.alignment = { vertical: 'middle', horizontal: 'center' };
        });
        ws.getRow(2).height = 24;

        // Data rows starting at 3
        var nextRow = 3;
        rows.forEach(function(row) {
          var excelRow = ws.getRow(nextRow);
          flatCols.forEach(function(c, i) {
            var cell = excelRow.getCell(i + 1);
            var v = row[c];
            if (v != null && v !== '') {
              var num = parseFloat(v);
              cell.value = !isNaN(num) ? num : v;
            }
          });
          excelRow.commit();
          nextRow++;
        });

        // Thin borders
        for (var r = 1; r < nextRow; r++) {
          var row = ws.getRow(r);
          for (var c = 1; c <= flatCols.length; c++) {
            row.getCell(c).border = {
              top: { style: 'thin' }, left: { style: 'thin' },
              bottom: { style: 'thin' }, right: { style: 'thin' }
            };
          }
        }

        ws.views = [{ state: 'frozen', xSplit: Math.min(4, flatCols.length), ySplit: 2 }];

        var buffer = await wb.xlsx.writeBuffer();
        var blob = new Blob([buffer], { type: 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet' });
        var url = window.URL.createObjectURL(blob);
        var a = document.createElement('a');
        a.href = url;
        a.download = 'Protein_List_QUANT_' + new Date().toISOString().split('T')[0] + '.xlsx';
        a.click();
        window.URL.revokeObjectURL(url);
      } catch(e) {
        console.error('Error while exporting:', e);
      }
    }
  ", cond_struct_json, group_struct_json, palette_json,
      element_id, element_id, element_id, element_id, element_id, element_id)))

  search_actions <- div(class = "rl-search-actions",
    tags$input(
      type = "search",
      placeholder = "Search proteins, genes, descriptions...",
      class = "rl-search-input ql-search-input",
      oninput = sprintf("Reactable.setSearch('%s', this.value)", element_id)
    ),
    div(class = "rl-action-buttons",
      tags$button(
        class = "rl-btn-action ql-btn-toggle-filters filters-hidden",
        onclick = "qlToggleFilters()",
        title = "Show/Hide filters",
        HTML('<svg xmlns="http://www.w3.org/2000/svg" width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><polygon points="22 3 2 3 10 12.46 10 19 14 21 14 12.46 22 3"></polygon></svg>')
      ),
      tags$button(
        class = "rl-btn-action",
        onclick = "qlClearFilters()",
        title = "Clear filters",
        HTML('<svg xmlns="http://www.w3.org/2000/svg" width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M3 6h18"></path><path d="M8 6V4c0-1 1-2 2-2h4c1 0 2 1 2 2v2"></path><path d="M19 6l-1 14c0 1-1 2-2 2H8c-1 0-2-1-2-2L5 6"></path><line x1="1" y1="1" x2="23" y2="23" stroke="#E63946" stroke-width="2"></line></svg>')
      ),
      tags$button(
        class = "rl-btn-action",
        onclick = "qlExportExcel()",
        title = "Export to Excel",
        HTML('<svg xmlns="http://www.w3.org/2000/svg" width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4"></path><polyline points="7 10 12 15 17 10"></polyline><line x1="12" y1="15" x2="12" y2="3"></line></svg>')
      )
    )
  )

  cond_chips <- lapply(conditions, function(cc) {
    tags$span(
      class = "pl-cond-chip ql-cond-chip",
      `data-cond` = cc,
      style = sprintf("background-color:%s;", palette[[cc]]),
      onclick = sprintf("qlToggleCondition('%s')", cc),
      cc
    )
  })

  .first_of <- function(mm) {
    cc <- sample_map$column[sample_map$metric == mm]
    if (length(cc) == 0) NA_character_ else cc[1]
  }
  metric_filters <- list(
    list(label = "# Uniq. PSMs Quantified",  anchor = .first_of("PG.NrOfPrecursorsUsedForQuantification"),        placeholder = ">= 2 ..."),
    list(label = "# Uniq. Pepts Quantified", anchor = .first_of("PG.NrOfStrippedSequencesUsedForQuantification"), placeholder = ">= 2 ..."),
    list(label = "Raw Abundance",            anchor = .first_of("PG.Quantity"),                                   placeholder = ">= 1e5 ..."),
    list(label = "Norm. Imp. Log2",          anchor = .first_of("NormImp.Log2"),                                  placeholder = ">= 10 ...")
  )
  metric_filter_items <- lapply(Filter(function(f) !is.na(f$anchor), metric_filters), function(f) {
    div(class = "rl-filter-item",
      tags$label(class = "rl-filter-label", f$label),
      tags$input(
        type = "text",
        class = "rl-numeric-input",
        `data-column` = f$anchor,
        placeholder = f$placeholder,
        oninput = sprintf("qlApplyNumericFilter('%s', this.value)", f$anchor)
      ),
      tags$span(class = "rl-filter-hint", "OR across the samples in the block")
    )
  })

  filters_panel <- div(class = "rl-filters-container ql-filters-container pl-filters-container",
    div(class = "rl-filters-row",
      div(class = "rl-filter-item", style = "flex: 2 1 300px;",
        tags$label(class = "rl-filter-label", "Conditions (click to hide/show 16 cols)"),
        div(class = "pl-cond-chips", cond_chips)
      ),
      metric_filter_items
    )
  )

  tbl <- reactable(
    df,
    elementId     = element_id,
    defaultSorted = list(PG.ProteinGroups = "asc"),
    defaultPageSize    = page_size,
    showPageSizeOptions = TRUE,
    pageSizeOptions = c(15, 30, 50, 100),
    resizable   = TRUE,
    selection   = selection,
    onClick     = if (!is.null(selection)) "select" else NULL,
    defaultColDef = colDef(
      align = "left",
      headerStyle = list(
        background  = "rgba(102, 5, 5, 0.9)",
        color       = "#ffffff",
        height      = "45px",
        display     = "flex",
        alignItems  = "center",
        justifyContent = "center"
      ),
      style = list(height = "48px", display = "flex", alignItems = "center")
    ),
    columns      = cols,
    columnGroups = groups,
    wrap         = FALSE,
    class        = "rl-table pl-table ql-table",
    rowStyle     = if (!is.null(selection)) list(cursor = "pointer") else NULL,
    highlight    = TRUE,
    searchable   = searchable,
    height       = height,
    striped      = TRUE,
    theme        = .ql_theme(),
    language     = .rl_lang(),
    details      = .pl_detail_row()
  )

  browsable(tagList(
    css,
    export_deps,
    js_code,
    search_actions,
    filters_panel,
    tbl
  ))
}


# =============================================================================
# Main function -- Summary (summary_list_widget)
# =============================================================================

#' Interactive Reactable Table for the Sample Metadata
#'
#' Builds a reactable table that reproduces the layout of the metadata Excel
#' sheet (FileName, Condition, Replicate, Coding, # Unique PSMs, # Unique
#' Peptides, # Protein Groups). Black header, rows tinted by condition in light
#' shades and Condition/Coding cells with a strong colour chip (the same scheme
#' used by protein_list_widget() / quant_list_widget()).
#'
#' @param data Data frame or path to a metadata TSV/CSV file. It must have the
#'   columns R.FileName, R.Condition, R.Replicate, Coding, R.PrecursorsIdentified,
#'   R.StrippedSequencesIdentified, R.ProteinGroupsIdentified.
#' @param page_size Page size (default 16, i.e. all the rows).
#' @param height Height in px of the table container (default 540).
#' @param element_id Widget id in the DOM (default "summary_table").
#' @param selection Selection type ("multiple" or NULL).
#' @param searchable Enable the global search box (default FALSE).
#'
#' @return An htmltools object (browsable tagList) ready to be rendered.
#'
#' @examples
#' if (requireNamespace("jsonlite", quietly = TRUE)) {
#'   data(nadia_dia)
#'
#'   # LFQ/TMT metadata lack the identification counts; the widget then shows
#'   # only FileName/Condition/Replicate/Coding
#'   w <- summary_list_widget(nadia_dia$metadata)
#'   print(class(w))
#' }
#'
#' @export
summary_list_widget <- function(
    data       = "data-raw/Metadata_20260423_142504.tsv",
    page_size  = 16,
    height     = 540,
    element_id = "summary_table",
    selection  = NULL,
    searchable = FALSE
) {

  df <- .sl_load_metadata(data)
  conditions   <- sort(unique(df$Condition))
  palette      <- .pl_condition_palette(conditions)
  tint_palette <- setNames(
    vapply(palette, .sl_tint_color, character(1), 0.12),
    names(palette)
  )

  cols   <- .sl_build_columns(df, palette)
  rstyle <- .sl_build_row_style()

  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("The 'jsonlite' package is required by summary_list_widget().")
  }
  palette_json      <- jsonlite::toJSON(as.list(palette), auto_unbox = TRUE)
  tint_palette_json <- jsonlite::toJSON(as.list(tint_palette), auto_unbox = TRUE)
  conditions_json   <- jsonlite::toJSON(as.list(conditions), auto_unbox = FALSE)

  # Columns to export (only the ones present; LFQ/TMT do not carry the counts)
  sl_width_map <- c(FileName = 45, Condition = 14, Replicate = 12, Coding = 14,
                    `# Unique PSMs` = 18, `# Unique Peptides` = 18,
                    `# Protein Groups` = 18)
  sl_export_cols    <- names(df)
  sl_export_widths  <- unname(sl_width_map[sl_export_cols])
  sl_export_widths[is.na(sl_export_widths)] <- 16
  export_cols_json   <- jsonlite::toJSON(sl_export_cols)
  export_widths_json <- jsonlite::toJSON(sl_export_widths)

  # --- CSS + CDN scripts ---
  css <- .rl_css()
  # --- ExcelJS and PapaParse, shipped with the package (see .rl_export_deps) ---
  export_deps <- .rl_export_deps()

  # --- Embedded JavaScript (slToggleFilters / slToggleCondition / slClearFilters / slExportExcel) ---
  js_code <- tags$script(HTML(sprintf("
    var slFiltersVisible = false;
    var slHiddenConditions = {};
    var SL_PALETTE = %s;
    var SL_TINT_PALETTE = %s;
    var SL_CONDITIONS = %s;
    var SL_EXPORT_COLS = %s;
    var SL_EXPORT_WIDTHS = %s;

    function slToggleFilters() {
      var container = document.querySelector('.sl-filters-container');
      var btn = document.querySelector('.sl-btn-toggle-filters');
      if (!container) return;
      if (slFiltersVisible) {
        container.style.maxHeight = '0';
        container.style.opacity = '0';
        container.style.marginBottom = '0';
        container.style.padding = '0';
        container.style.borderWidth = '0';
        container.style.boxShadow = 'none';
        container.style.overflow = 'hidden';
        if (btn) btn.classList.add('filters-hidden');
      } else {
        container.style.maxHeight = '500px';
        container.style.opacity = '1';
        container.style.marginBottom = '1rem';
        container.style.padding = '1.5rem';
        container.style.borderWidth = '1px';
        container.style.boxShadow = '0 2px 4px rgba(0, 0, 0, 0.05)';
        container.style.overflow = 'visible';
        if (btn) btn.classList.remove('filters-hidden');
      }
      slFiltersVisible = !slFiltersVisible;
    }

    function slApplyNumericFilter(columnId, value) {
      Reactable.setFilter('%s', columnId, value || undefined);
    }

    function slApplyConditionFilter() {
      var hidden = [];
      Object.keys(slHiddenConditions).forEach(function(c) {
        if (slHiddenConditions[c]) hidden.push(c);
      });
      Reactable.setFilter('%s', 'Condition', hidden.length > 0 ? hidden.join(',') : undefined);
    }

    function slToggleCondition(cond) {
      slHiddenConditions[cond] = !slHiddenConditions[cond];
      var chip = document.querySelector('.sl-cond-chip[data-cond=\"' + cond + '\"]');
      if (chip) chip.classList.toggle('off', slHiddenConditions[cond]);
      slApplyConditionFilter();
    }

    function slClearFilters() {
      var numInputs = document.querySelectorAll('.sl-filters-container .rl-numeric-input');
      numInputs.forEach(function(inp) {
        inp.value = '';
        var col = inp.getAttribute('data-column');
        if (col) Reactable.setFilter('%s', col, undefined);
      });
      Object.keys(slHiddenConditions).forEach(function(c) {
        slHiddenConditions[c] = false;
        var chip = document.querySelector('.sl-cond-chip[data-cond=\"' + c + '\"]');
        if (chip) chip.classList.remove('off');
      });
      Reactable.setFilter('%s', 'Condition', undefined);
      var searchInput = document.querySelector('.sl-search-input');
      if (searchInput) {
        searchInput.value = '';
        Reactable.setSearch('%s', '');
      }
    }

    function slHex2Argb(hex) {
      var h = (hex || '').replace('#', '');
      if (h.length !== 6) return 'FF808080';
      return 'FF' + h.toUpperCase();
    }

    async function slExportExcel() {
      try {
        var tsv = Reactable.getDataCSV('%s', { sep: '\\t' });
        var parseResult = Papa.parse(tsv, { header: true, delimiter: '\\t', skipEmptyLines: true });
        var rows = parseResult.data;

        var wb = new ExcelJS.Workbook();
        var ws = wb.addWorksheet('Summary');

        var flatCols = SL_EXPORT_COLS;
        var colWidths = SL_EXPORT_WIDTHS;
        flatCols.forEach(function(c, i) { ws.getColumn(i + 1).width = colWidths[i]; });

        // Row 1: header (black / white bold)
        flatCols.forEach(function(c, i) {
          var cell = ws.getCell(1, i + 1);
          cell.value = c;
          cell.fill = { type: 'pattern', pattern: 'solid', fgColor: { argb: 'FF1A1A1A' } };
          cell.font = { bold: true, color: { argb: 'FFFFFFFF' }, size: 11 };
          cell.alignment = { vertical: 'middle', horizontal: 'center' };
        });
        ws.getRow(1).height = 28;

        // Data rows: per-row tint + strong chip on Condition + Coding
        var nextRow = 2;
        rows.forEach(function(row) {
          var cond = String(row['Condition'] || '');
          var tintArgb = SL_TINT_PALETTE[cond] ? slHex2Argb(SL_TINT_PALETTE[cond]) : 'FFFFFFFF';
          var strongArgb = SL_PALETTE[cond] ? slHex2Argb(SL_PALETTE[cond]) : 'FF808080';

          var excelRow = ws.getRow(nextRow);
          flatCols.forEach(function(c, i) {
            var cell = excelRow.getCell(i + 1);
            var v = row[c];
            if (v != null && v !== '') {
              if (i >= 4) { // numeric columns (5/6/7)
                var num = parseFloat(v);
                cell.value = !isNaN(num) ? num : v;
                cell.numFmt = '#,##0';
              } else {
                cell.value = v;
              }
            }
            // tinted background across the whole row
            cell.fill = { type: 'pattern', pattern: 'solid', fgColor: { argb: tintArgb } };
            cell.alignment = { vertical: 'middle',
                               horizontal: (i === 0 ? 'left' : (i >= 4 ? 'right' : 'center')) };
            // override on Condition + Coding (i = 1 or 3)
            if (i === 1 || i === 3) {
              cell.fill = { type: 'pattern', pattern: 'solid', fgColor: { argb: strongArgb } };
              cell.font = { bold: true, color: { argb: 'FFFFFFFF' } };
              cell.alignment = { vertical: 'middle', horizontal: 'center' };
            }
          });
          excelRow.commit();
          nextRow++;
        });

        // Borders
        for (var r = 1; r < nextRow; r++) {
          var rr = ws.getRow(r);
          for (var c = 1; c <= flatCols.length; c++) {
            rr.getCell(c).border = {
              top:    { style: 'thin', color: { argb: 'FFBFBFBF' } },
              left:   { style: 'thin', color: { argb: 'FFBFBFBF' } },
              bottom: { style: 'thin', color: { argb: 'FFBFBFBF' } },
              right:  { style: 'thin', color: { argb: 'FFBFBFBF' } }
            };
          }
        }
        ws.views = [{ state: 'frozen', xSplit: 1, ySplit: 1 }];

        var buffer = await wb.xlsx.writeBuffer();
        var blob = new Blob([buffer], { type: 'application/vnd.openxmlformats-officedocument.spreadsheetml.sheet' });
        var url = window.URL.createObjectURL(blob);
        var a = document.createElement('a');
        a.href = url;
        a.download = 'Summary_' + new Date().toISOString().split('T')[0] + '.xlsx';
        a.click();
        window.URL.revokeObjectURL(url);
      } catch(e) {
        console.error('Error while exporting:', e);
      }
    }
  ", palette_json, tint_palette_json, conditions_json,
      export_cols_json, export_widths_json,
      element_id, element_id, element_id, element_id, element_id, element_id)))

  # --- Search bar + buttons ---
  search_input <- if (isTRUE(searchable)) {
    tags$input(
      type = "search",
      placeholder = "Search files, codings, conditions...",
      class = "rl-search-input sl-search-input",
      oninput = sprintf("Reactable.setSearch('%s', this.value)", element_id)
    )
  } else NULL

  search_actions <- div(class = "rl-search-actions",
    search_input,
    div(class = "rl-action-buttons",
      tags$button(
        class = "rl-btn-action sl-btn-toggle-filters filters-hidden",
        onclick = "slToggleFilters()",
        title = "Show/Hide filters",
        HTML('<svg xmlns="http://www.w3.org/2000/svg" width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><polygon points="22 3 2 3 10 12.46 10 19 14 21 14 12.46 22 3"></polygon></svg>')
      ),
      tags$button(
        class = "rl-btn-action",
        onclick = "slClearFilters()",
        title = "Clear filters",
        HTML('<svg xmlns="http://www.w3.org/2000/svg" width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M3 6h18"></path><path d="M8 6V4c0-1 1-2 2-2h4c1 0 2 1 2 2v2"></path><path d="M19 6l-1 14c0 1-1 2-2 2H8c-1 0-2-1-2-2L5 6"></path><line x1="1" y1="1" x2="23" y2="23" stroke="#E63946" stroke-width="2"></line></svg>')
      ),
      tags$button(
        class = "rl-btn-action",
        onclick = "slExportExcel()",
        title = "Export to Excel",
        HTML('<svg xmlns="http://www.w3.org/2000/svg" width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4"></path><polyline points="7 10 12 15 17 10"></polyline><line x1="12" y1="15" x2="12" y2="3"></line></svg>')
      )
    )
  )

  # --- Filter panel: condition chips + 3 numeric inputs ---
  cond_chips <- lapply(conditions, function(cc) {
    tags$span(
      class = "pl-cond-chip sl-cond-chip",
      `data-cond` = cc,
      style = sprintf("background-color:%s;", palette[[cc]]),
      onclick = sprintf("slToggleCondition('%s')", cc),
      cc
    )
  })

  numeric_filters <- list(
    list(label = "# Unique PSMs",     col = "# Unique PSMs",     placeholder = ">= 120000 ..."),
    list(label = "# Unique Peptides", col = "# Unique Peptides", placeholder = ">= 90000 ..."),
    list(label = "# Protein Groups",  col = "# Protein Groups",  placeholder = ">= 9000 ...")
  )
  # Only keep filters for the count columns present (absent for LFQ/TMT)
  numeric_filters <- Filter(function(f) f$col %in% names(df), numeric_filters)
  numeric_filter_items <- lapply(numeric_filters, function(f) {
    div(class = "rl-filter-item",
      tags$label(class = "rl-filter-label", f$label),
      tags$input(
        type = "text",
        class = "rl-numeric-input",
        `data-column` = f$col,
        placeholder = f$placeholder,
        oninput = sprintf("slApplyNumericFilter('%s', this.value)", f$col)
      ),
      tags$span(class = "rl-filter-hint", "Operators: >=, <=, >, <, =, !=")
    )
  })

  filters_panel <- div(class = "rl-filters-container sl-filters-container",
    div(class = "rl-filters-row",
      div(class = "rl-filter-item", style = "flex: 2 1 300px;",
        tags$label(class = "rl-filter-label", "Conditions (click to show/hide rows)"),
        div(class = "pl-cond-chips sl-cond-chips", cond_chips)
      ),
      numeric_filter_items
    )
  )

  # --- Condition filter (hidden): drops rows whose condition is in the list ---
  cols$Condition$filterable   <- TRUE
  cols$Condition$filterInput  <- .numeric_filter_hidden
  cols$Condition$filterMethod <- reactable::JS("
    function(rows, columnId, filterValue) {
      if (!filterValue) return rows;
      var exclude = String(filterValue).split(',').map(function(s){return s.trim();}).filter(Boolean);
      if (exclude.length === 0) return rows;
      return rows.filter(function(row) {
        return exclude.indexOf(String(row.values[columnId])) === -1;
      });
    }
  ")

  # --- Table ---
  tbl <- reactable(
    df,
    elementId           = element_id,
    defaultPageSize     = page_size,
    showPageSizeOptions = FALSE,
    pagination          = nrow(df) > page_size,
    resizable           = TRUE,
    selection           = selection,
    onClick             = if (!is.null(selection)) "select" else NULL,
    defaultColDef = colDef(
      align = "left",
      headerStyle = list(
        background     = "#1a1a1a",
        color          = "#ffffff",
        height         = "45px",
        display        = "flex",
        alignItems     = "center",
        justifyContent = "center"
      ),
      style = list(height = "44px", display = "flex", alignItems = "center")
    ),
    columns      = cols,
    columnGroups = NULL,
    wrap         = FALSE,
    class        = "rl-table sl-table",
    rowStyle     = rstyle,
    highlight    = TRUE,
    searchable   = searchable,
    height       = height,
    striped      = FALSE,
    theme        = .sl_theme(),
    language     = .rl_lang()
  )

  browsable(tagList(
    css,
    export_deps,
    js_code,
    search_actions,
    filters_panel,
    tbl
  ))
}
