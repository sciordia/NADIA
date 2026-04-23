
# Cargar las librerias
library(reactable)
library(htmltools)

# Operador null-coalesce
if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
}


# =============================================================================
# Helpers internos
# =============================================================================

# --- Filtro numerico compartido: input oculto + JS con operadores AND/OR ---
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


#' Cargar y validar datos de expresion diferencial
#' @param input Data frame o ruta a archivo TSV/Parquet
#' @return Data frame validado
#' @noRd
.rl_load_data <- function(input) {
  if (is.character(input) && length(input) == 1) {
    if (!file.exists(input)) {
      stop("Archivo no encontrado: ", input)
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
        stop("El paquete 'arrow' es necesario para leer archivos Parquet")
      }
      input <- arrow::read_parquet(input)
    } else if (ext == "csv") {
      if (requireNamespace("readr", quietly = TRUE)) {
        input <- readr::read_csv(input, show_col_types = FALSE)
      } else {
        input <- read.csv(input, stringsAsFactors = FALSE)
      }
    } else {
      stop("Formato no soportado: ", ext, ". Usa TSV, CSV o Parquet.")
    }
  }

  df <- as.data.frame(input)

  # Validar columnas requeridas
  required <- c("Protein.IDs", "Gene.Names", "logFC", "P.Value", "adj.P.Val",
                 "Change", "Comparison")
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop("Columnas requeridas faltantes: ", paste(missing, collapse = ", "))
  }

  # Coercion de tipos
  df$logFC     <- as.numeric(df$logFC)
  df$P.Value   <- as.numeric(df$P.Value)
  df$adj.P.Val <- as.numeric(df$adj.P.Val)

  df
}


#' Join con protein_quant para anadir Description y Quant_Pepts
#' @param df Data frame de resultados DE
#' @param protein_quant Data frame (preprocessing$protein_quant) o ruta a TSV.
#'   NULL para no hacer join.
#' @return Data frame con columnas Description y Quant_Pepts anadidas
#' @noRd
.rl_join_protein_info <- function(df, protein_quant) {
  if (is.null(protein_quant)) return(df)

  # Cargar si es ruta
  if (is.character(protein_quant) && length(protein_quant) == 1) {
    pq <- .rl_load_data(protein_quant)
  } else {
    pq <- as.data.frame(protein_quant)
  }

  # Verificar columna clave
  if (!"PG.ProteinGroups" %in% names(pq)) {
    stop("protein_quant debe contener la columna 'PG.ProteinGroups'")
  }

  # Extraer Description
  desc_col <- if ("PG.ProteinDescriptions" %in% names(pq)) pq$PG.ProteinDescriptions else NA_character_
  info <- data.frame(
    Protein.IDs = pq$PG.ProteinGroups,
    Description = desc_col,
    stringsAsFactors = FALSE
  )

  # Calcular max Quant_Pepts
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

  # Join por Protein.IDs (match para preservar orden y evitar duplicacion)
  idx <- match(df$Protein.IDs, info$Protein.IDs)
  df$Description  <- info$Description[idx]
  df$Quant_Pepts  <- info$Quant_Pepts[idx]

  df
}


#' Tema reactable estilo teal/green
#' @return Objeto reactableTheme
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


#' CSS embebido para badges, barras, detalle, filtros y export
#' @return Objeto tags$style
#' @noRd
.rl_css <- function() {
  tags$style(HTML("
    /* Etiquetas de estado */
    .tag {
      display: inline-block;
      padding: 0.1rem 0.5rem;
      border-radius: 10px;
      font-weight: bold;
      font-size: 14px;
      background-color: #eee;
      line-height: 1.5;
    }
    .status-green {
      color: hsl(121, 33%, 23%);
      background-color: hsl(121, 33%, 86%);
    }
    .status-red {
      color: hsl(0, 75%, 32%);
      background-color: hsl(0, 75%, 92%);
    }
    .status-grey {
      color: hsl(0, 0%, 35%);
      background-color: hsl(0, 0%, 90%);
    }

    /* Columnas ordenadas */
    .sorted {
      background: rgba(14, 102, 85, 0.1);
    }

    /* Tabla principal */
    .rl-table {
      width: 100% !important;
      overflow-x: auto;
      font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif;
      font-size: 14px;
      margin-top: 0.25rem;
      border: 1px solid hsl(213, 33%, 93%);
      border-radius: 4px;
      box-shadow: 0 4px 8px 0 rgba(0, 0, 0, 0.1);
    }

    /* Barra logFC */
    .lfc-bar-container {
      display: flex;
      align-items: center;
      gap: 6px;
      width: 100%;
    }
    .lfc-value {
      min-width: 60px;
      text-align: right;
      font-variant-numeric: tabular-nums;
      font-weight: 700;
    }
    .lfc-value.up       { color: #02905A; }
    .lfc-value.down     { color: #E63946; }
    .lfc-value.nochange { color: #6C757D; }
    .lfc-bar-wrapper {
      flex: 1;
      height: 14px;
      background: #f0f0f0;
      border-radius: 2px;
      overflow: hidden;
      position: relative;
    }
    .lfc-bar {
      height: 100%;
      border-radius: 2px;
      min-width: 2px;
    }
    .lfc-bar.up { background-color: #02905A; }
    .lfc-bar.down { background-color: #E63946; }
    .lfc-bar.nochange { background-color: #ADB5BD; }

    /* Panel de detalle expandido */
    .rl-detail {
      padding: 12px 20px;
      margin: 8px 40px;
      background: #f8f9fa;
      border-left: 4px solid rgba(204, 153, 0, 0.9);
      border-radius: 4px;
      font-size: 0.9rem;
      line-height: 1.6;
    }
    .rl-detail .detail-label {
      font-weight: 600;
      color: #0E6655;
      margin-right: 6px;
    }
    .rl-detail .detail-row {
      margin-bottom: 4px;
    }
    .rl-detail a {
      color: #0E6655;
      text-decoration: none;
    }
    .rl-detail a:hover {
      text-decoration: underline;
    }

    /* Protein count badge */
    .protein-count {
      display: inline-block;
      padding: 0 5px;
      border-radius: 8px;
      font-size: 11px;
      font-weight: 600;
      color: #666;
      background-color: #e9ecef;
      margin-left: 4px;
    }

    /* Ocultar busqueda por defecto y fila de filtros inline de reactable */
    .rt-search {
      display: none !important;
    }
    .rl-table .rt-thead.-filters,
    .rl-table .rt-tr.-filters,
    .rl-table .rt-thead .rt-tr-filters,
    .rl-table [class*='filterRow'],
    .rl-table .rt-thead .rt-th.-filter {
      display: none !important;
      height: 0 !important;
      overflow: hidden !important;
      padding: 0 !important;
      margin: 0 !important;
      border: none !important;
    }

    /* Contenedor de busqueda + botones de accion */
    .rl-search-actions {
      display: flex;
      gap: 0.75rem;
      margin-bottom: 0.5rem;
      align-items: center;
      width: 100%;
    }
    .rl-search-input {
      flex: 1;
      font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif;
      font-weight: 400;
      font-size: 14px;
      padding: 8px 12px;
      border: 2px solid #e0e0e0;
      border-radius: 8px;
      transition: all 0.3s ease;
      background-color: #ffffff;
      box-shadow: 0 2px 4px rgba(0, 0, 0, 0.05);
      height: 44px;
      box-sizing: border-box;
    }
    .rl-search-input:focus {
      border-color: #0E6655;
      box-shadow: 0 0 8px rgba(14, 102, 85, 0.3);
      outline: none;
      background-color: #fafafa;
    }
    .rl-search-input:hover {
      border-color: #0E6655;
      box-shadow: 0 2px 8px rgba(14, 102, 85, 0.1);
    }

    /* Botones de accion (filtros, export, limpiar) */
    .rl-action-buttons {
      display: flex;
      gap: 0.5rem;
      flex-shrink: 0;
    }
    .rl-btn-action {
      background-color: #ffffff;
      border: 2px solid #e0e0e0;
      border-radius: 8px;
      width: 44px;
      height: 44px;
      display: flex;
      align-items: center;
      justify-content: center;
      cursor: pointer;
      transition: all 0.3s ease;
      color: #495057;
      position: relative;
    }
    .rl-btn-action:hover {
      border-color: #0E6655;
      color: #0E6655;
      background-color: rgba(14, 102, 85, 0.05);
      transform: translateY(-1px);
      box-shadow: 0 2px 4px rgba(0, 0, 0, 0.1);
    }
    .rl-btn-action.filters-hidden {
      opacity: 0.5;
    }

    /* Tooltip en botones */
    .rl-btn-action:hover::before {
      content: attr(title);
      position: absolute;
      bottom: calc(100% + 8px);
      left: 50%;
      transform: translateX(-50%);
      padding: 6px 12px;
      background-color: rgba(33, 37, 41, 0.9);
      color: white;
      font-size: 12px;
      border-radius: 4px;
      white-space: nowrap;
      pointer-events: none;
      z-index: 1000;
    }

    /* Contenedor de filtros (colapsable) */
    .rl-filters-container {
      background: #f8f9fa;
      border: 0 solid #dee2e6;
      border-radius: 8px;
      padding: 0;
      margin-bottom: 0;
      box-shadow: none;
      max-height: 0;
      opacity: 0;
      overflow: hidden;
      transition: max-height 0.4s ease, opacity 0.3s ease,
                  margin-bottom 0.3s ease, padding 0.3s ease,
                  border-width 0.3s ease, box-shadow 0.3s ease;
    }
    .rl-filters-row {
      display: flex;
      flex-wrap: wrap;
      gap: 0.75rem;
    }
    .rl-filters-row > .rl-filter-item {
      flex: 1 1 0;
      min-width: 100px;
    }
    .rl-filter-item {
      display: flex;
      flex-direction: column;
      gap: 0.15rem;
      font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif;
      font-size: 12.5px;
    }
    .rl-filter-label {
      font-weight: normal;
      color: #495057;
      font-size: 12.5px;
      margin-bottom: 0.1rem;
    }

    /* Selectize dropdowns dentro de filtros */
    .rl-filter-item .selectize-input {
      border: 2px solid #e0e0e0;
      border-radius: 6px;
      padding: 5px 8px;
      font-size: 12.5px;
      min-height: 32px;
      transition: all 0.3s ease;
    }
    .rl-filter-item .selectize-input:hover {
      border-color: #0E6655;
    }
    .rl-filter-item .selectize-input.focus {
      border-color: #0E6655;
      box-shadow: 0 0 0 0.2rem rgba(14, 102, 85, 0.25);
    }
    .rl-filter-item .selectize-dropdown {
      border: 2px solid #0E6655;
      border-top: none;
      border-radius: 0 0 6px 6px;
      box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
      font-size: 12.5px;
    }
    .rl-filter-item .selectize-dropdown .active {
      background-color: rgba(14, 102, 85, 0.1);
      color: #0E6655;
    }
    .selectize-control.multi .selectize-input > div {
      background: #0E6655;
      color: white;
      border-radius: 4px;
      padding: 2px 8px;
      margin-right: 5px;
    }
    .selectize-control.multi .selectize-input > div .remove {
      color: rgba(255, 255, 255, 0.7);
    }
    .selectize-control.multi .selectize-input > div .remove:hover {
      color: #ffa62d;
    }

    /* Radio buttons estilo teal */
    .rl-table .rt-select-input[type='radio'] {
      background-color: #ffffff;
      border: 1px solid #0E6655;
      appearance: none;
      width: 16px;
      height: 16px;
      border-radius: 50%;
    }
    .rl-table .rt-select-input[type='radio']:checked {
      background-color: #0E6655;
      border: 1px solid #0E6655;
      box-shadow: inset 0 0 0 2px #fff;
    }

    /* Inputs numericos en panel de filtros */
    .rl-numeric-input {
      font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif;
      font-size: 12.5px;
      width: 100%;
      padding: 5px 8px;
      border: 2px solid #e0e0e0;
      border-radius: 6px;
      min-height: 32px;
      transition: all 0.3s ease;
      box-sizing: border-box;
      outline: none;
    }
    .rl-numeric-input:hover {
      border-color: #0E6655;
    }
    .rl-numeric-input:focus {
      border-color: #0E6655;
      box-shadow: 0 0 0 0.2rem rgba(14, 102, 85, 0.25);
    }
    .rl-filter-hint {
      font-size: 11px;
      color: #999;
      margin-top: 2px;
    }

    /* Centrado vertical en celdas Missing% */
    .rl-table .rt-td.missing-cell {
      display: flex !important;
      align-items: center;
      justify-content: center;
    }

    /* Pagination focus */
    .rl-table .rt-page-size-select:focus {
      border-color: #0E6655;
      box-shadow: 0 0 5px #0E6655;
      outline: none;
    }

    /* ==========================================================
       Protein List (pl-*) — tabla Protein_ID (post-preprocessing)
       ========================================================== */

    /* Data bar + valor para celdas compactas de 64 muestras */
    .pl-bar-wrapper {
      position: relative;
      width: 100%;
      height: 24px;
      background: #f5f5f7;
      border-radius: 3px;
      overflow: hidden;
      display: flex;
      align-items: center;
      justify-content: flex-end;
    }
    .pl-bar {
      position: absolute;
      left: 0; top: 0; bottom: 0;
      border-radius: 3px;
      opacity: 0.45;
    }
    .pl-bar-value {
      position: relative;
      z-index: 1;
      padding: 0 6px;
      font-variant-numeric: tabular-nums;
      font-size: 14px;
      font-weight: 600;
      color: #212529;
    }
    .pl-plain-value {
      font-variant-numeric: tabular-nums;
      font-size: 14px;
      font-weight: 600;
      color: #212529;
    }
    .pl-bar-empty {
      color: #bbb;
      font-weight: 400;
    }

    /* Colores por condicion (default paleta Office) */
    .pl-cond-A { background-color: #4F81BD; }
    .pl-cond-B { background-color: #9BBB59; }
    .pl-cond-C { background-color: #F79646; }
    .pl-cond-D { background-color: #8064A2; }
    .pl-cond-E { background-color: #4BACC6; }
    .pl-cond-F { background-color: #C0504D; }
    .pl-cond-G { background-color: #9F8A76; }
    .pl-cond-H { background-color: #646464; }

    /* Header styling: primera fila (grupos) y segunda (muestras) */
    .pl-table .rt-thead {
      font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif;
    }
    .pl-table .rt-tr-groups .rt-th {
      background: rgba(14, 102, 85, 0.9);
      color: #ffffff;
      font-weight: 600;
      border-right: 1px solid rgba(255,255,255,0.15);
    }
    .pl-table .rt-tr-header .rt-th {
      font-size: 14px;
    }

    /* Separador vertical entre bloques metricos */
    .pl-table .rt-td.pl-group-boundary,
    .pl-table .rt-th.pl-group-boundary {
      border-left: 2px solid rgba(33, 37, 41, 0.35) !important;
    }
    .pl-hdr-A { background: rgba(79, 129, 189, 0.85) !important; color: #ffffff !important; }
    .pl-hdr-B { background: rgba(155, 187, 89, 0.85) !important; color: #ffffff !important; }
    .pl-hdr-C { background: rgba(247, 150, 70, 0.85) !important; color: #ffffff !important; }
    .pl-hdr-D { background: rgba(128, 100, 162, 0.85) !important; color: #ffffff !important; }
    .pl-hdr-E { background: rgba(75, 172, 198, 0.85) !important; color: #ffffff !important; }
    .pl-hdr-F { background: rgba(192, 80, 77, 0.85)  !important; color: #ffffff !important; }
    .pl-hdr-G { background: rgba(159, 138, 118, 0.85) !important; color: #ffffff !important; }
    .pl-hdr-H { background: rgba(100, 100, 100, 0.85) !important; color: #ffffff !important; }

    /* Sticky: sombra lateral para marcar separacion */
    .pl-table .rt-td-sticky,
    .pl-table .rt-th-sticky {
      background-color: #ffffff !important;
      box-shadow: 2px 0 6px -3px rgba(0, 0, 0, 0.2);
    }
    .pl-table .rt-tr-striped .rt-td-sticky {
      background-color: #fafbfc !important;
    }
    .pl-table .rt-tr:hover .rt-td-sticky {
      background-color: rgba(2, 144, 82, 0.08) !important;
    }

    /* Celdas de muestra compactas (menos padding) */
    .pl-table .rt-td.pl-sample-cell {
      padding: 4px 4px !important;
    }

    /* Contenedor de chips para ocultar condiciones */
    .pl-cond-chips {
      display: flex;
      gap: 0.35rem;
      flex-wrap: wrap;
    }
    .pl-cond-chip {
      cursor: pointer;
      padding: 4px 12px;
      border-radius: 999px;
      font-size: 12.5px;
      font-weight: 600;
      color: #ffffff;
      border: 2px solid transparent;
      user-select: none;
      transition: all 0.15s ease;
    }
    .pl-cond-chip.off {
      opacity: 0.35;
      filter: grayscale(0.4);
    }
    .pl-cond-chip:hover {
      transform: translateY(-1px);
      box-shadow: 0 2px 4px rgba(0,0,0,0.15);
    }
  "))
}


#' Traduccion al espanol (interfaz)
#' @return Objeto reactableLang
#' @noRd
.rl_lang <- function() {
  reactableLang(
    searchPlaceholder = "Buscar...",
    pagePrevious      = "Anterior",
    pageNext          = "Siguiente",
    noData            = "No data to display",
    pageSizeOptions   = "Show {rows}",
    pageInfo          = "{rowStart}\u2013{rowEnd} of {rows} Proteins"
  )
}


#' Funcion de detalle para filas expandibles (JS renderer)
#' @param has_assay Logico, si el data frame tiene columna Assay
#' @param has_description Logico, si el data frame tiene columna Description
#' @return Objeto JS para el parametro details de reactable
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

    // Protein.IDs con enlaces UniProt
    var pids = (row['Protein.IDs'] || '').split(';').map(function(s) { return s.trim(); }).filter(Boolean);
    var links = pids.map(function(pid) {
      return '<a href=\"https://www.uniprot.org/uniprot/' + pid + '\" target=\"_blank\">' + pid + '</a>';
    }).join(' \\u00b7 ');
    html += '<div class=\"detail-row\"><span class=\"detail-label\">Proteins:</span> ' + links + '</div>';

    // Gene.Names completo
    var genes = row['Gene.Names'] || '';
    html += '<div class=\"detail-row\"><span class=\"detail-label\">Genes:</span> ' + genes + '</div>';

    // Description
    %s

    // P-valor y FDR con precision completa
    var pval = row['P.Value'];
    var fdr = row['adj.P.Val'];
    html += '<div class=\"detail-row\"><span class=\"detail-label\">P-value:</span> ' + (pval != null ? pval.toExponential(4) : '') + '</div>';
    html += '<div class=\"detail-row\"><span class=\"detail-label\">FDR:</span> ' + (fdr != null ? fdr.toExponential(4) : '') + '</div>';

    // Assay (si existe)
    %s

    html += '</div>';
    return React.createElement('div', { dangerouslySetInnerHTML: { __html: html } });
  }", desc_block, assay_block))
}


#' Construir definiciones de columnas compartidas
#' @param max_abs_lfc Valor maximo absoluto de logFC para escalar barras
#' @param alpha Umbral de significancia
#' @param has_assay Logico, si hay columna Assay
#' @param single_assay Logico, si solo hay un assay
#' @param show_missing Logico, si mostrar columnas Missing%
#' @param has_description Logico, si hay columna Description
#' @param has_quant_pepts Logico, si hay columna Quant_Pepts
#' @return Lista de colDef
#' @noRd
.rl_build_columns <- function(max_abs_lfc, alpha, has_assay, single_assay,
                               show_missing, has_description, has_quant_pepts) {

  # JS renderer para Missing% estilo rating con circulo de color
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
    # --- Comparison (1ro) ---
    Comparison = colDef(
      name = "Comparison",
      width = 120,
      align = "center"
    ),

    # --- Protein Groups (2do) ---
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

  # --- Description (3ro, condicional) ---
  if (has_description) {
    cols$Description <- colDef(
      name = "Description",
      minWidth = 250
    )
  }

  # --- Gene Names (4to) ---
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

  # --- Quant_Pepts (5to, condicional) ---
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

  # --- logFC (valor + barra coloreada por Change) ---
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

  # --- Missing% estilo rating con circulo de color ---
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
# Funcion principal
# =============================================================================

#' Tabla Reactable Interactiva para Resultados de Expresion Diferencial
#'
#' Genera una tabla reactable con formato profesional para visualizar
#' resultados de expresion diferencial de proteinas. Compatible con Shiny
#' (devuelve objeto reactable) y con uso standalone via \code{results_list_widget()}.
#'
#' @param data Data frame o ruta a archivo TSV/CSV/Parquet con resultados DE.
#'   Columnas requeridas: Protein.IDs, Gene.Names, logFC, P.Value, adj.P.Val,
#'   Change, Comparison. Opcionales: Assay, MissingGlobal, MissingPCT1, MissingPCT2
#' @param protein_quant Data frame (preprocessing$protein_quant) o ruta a archivo
#'   Protein_QUANT_*.tsv. Si no es NULL, anade columnas Description y Quant_Pepts.
#' @param comparisons Vector de comparaciones a incluir (NULL = todas)
#' @param ain Vector de assays a filtrar (NULL = todos)
#' @param alpha Umbral de significancia para resaltar FDR (default: 0.05)
#' @param lfc_thr Umbral de log2 fold-change (default: 0, reservado para uso futuro)
#' @param page_size Filas por pagina (default: 15)
#' @param height Altura de la tabla en pixeles (default: 720)
#' @param show_missing Mostrar columnas de porcentaje de ausencia (default: TRUE)
#' @param element_id ID del elemento para Reactable JS API (default: NULL)
#' @param selection Tipo de seleccion: "single", "multiple", o NULL (default: NULL)
#' @param searchable Habilitar busqueda interna (default: TRUE)
#'
#' @return Objeto reactable
#'
#' @examples
#' # Desde archivo
#' tbl <- results_list_reactable("results/VolcanoPlot_Input_cycloess_Impseq_min.tsv")
#'
#' # Con protein_quant para Description y Quant_Pepts
#' tbl <- results_list_reactable(de_res, protein_quant = preprocessing$protein_quant)
#'
#' # En Shiny
#' # output$tabla <- renderReactable({
#' #   results_list_reactable(data(), protein_quant = pq, element_id = "tabla")
#' # })
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

  # --- Carga y validacion ---
  df <- .rl_load_data(data)

  # --- Join con protein_quant ---
  df <- .rl_join_protein_info(df, protein_quant)

  # --- Filtrado ---
  if (!is.null(ain) && "Assay" %in% names(df)) {
    df <- df[df$Assay %in% ain, , drop = FALSE]
  }
  if (!is.null(comparisons)) {
    df <- df[df$Comparison %in% comparisons, , drop = FALSE]
  }
  if (nrow(df) == 0) {
    stop("No hay datos tras aplicar los filtros de comparaciones/assays")
  }

  if (nrow(df) > 15000 && is.null(comparisons)) {
    message("Nota: ", format(nrow(df), big.mark = "."),
            " filas. Considera filtrar por 'comparisons' para mejor rendimiento.")
  }

  # --- Detectar columnas opcionales ---
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

  # --- Reordenar columnas del data frame (reactable usa este orden visual) ---
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

#' Widget Completo con Filtros, Busqueda, Export y CSS
#'
#' Envuelve \code{results_list_reactable()} con filtros interactivos crosstalk
#' (Comparison, Change, Method), campo de busqueda, boton de exportar a Excel
#' y CSS embebido. El resultado es browsable: al imprimirlo en consola se abre
#' automaticamente en el Viewer de RStudio o en el navegador.
#'
#' @inheritParams results_list_reactable
#' @param element_id ID del elemento (default: "deps_table")
#'
#' @return Objeto htmltools browsable
#'
#' @examples
#' # Uso standalone
#' results_list_widget("results/VolcanoPlot_Input_cycloess_Impseq_min.tsv")
#'
#' # Con protein_quant
#' results_list_widget(de_res, protein_quant = preprocessing$protein_quant)
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
    stop("El paquete 'crosstalk' es necesario para filtros interactivos. ",
         "Inst\u00e1lalo con install.packages('crosstalk')")
  }

  # --- Carga, join y filtrado previo ---
  df <- .rl_load_data(data)
  df <- .rl_join_protein_info(df, protein_quant)

  if (!is.null(ain) && "Assay" %in% names(df)) {
    df <- df[df$Assay %in% ain, , drop = FALSE]
  }
  if (!is.null(comparisons)) {
    df <- df[df$Comparison %in% comparisons, , drop = FALSE]
  }
  if (nrow(df) == 0) {
    stop("No hay datos tras aplicar los filtros de comparaciones/assays")
  }

  # --- Detectar columnas opcionales ---
  has_missing <- all(c("MissingGlobal", "MissingPCT1", "MissingPCT2") %in% names(df))
  show_missing_cols <- show_missing && has_missing
  has_assay <- "Assay" %in% names(df)
  single_assay <- has_assay && length(unique(df$Assay)) == 1
  has_description <- "Description" %in% names(df)
  has_quant_pepts <- "Quant_Pepts" %in% names(df)

  max_abs_lfc <- max(abs(df$logFC), na.rm = TRUE)
  if (max_abs_lfc == 0) max_abs_lfc <- 1

  # --- Reordenar columnas del data frame (reactable usa este orden visual) ---
  desired_order <- c("Comparison", "Protein.IDs", "Description", "Gene.Names",
                     "Quant_Pepts", "Change", "logFC", "P.Value", "adj.P.Val",
                     "Assay", "MissingGlobal", "MissingPCT1", "MissingPCT2")
  desired_order <- intersect(desired_order, names(df))
  df <- df[, c(desired_order, setdiff(names(df), desired_order)), drop = FALSE]

  # --- SharedData para crosstalk ---
  shared_data <- crosstalk::SharedData$new(df)

  # --- CSS ---
  css <- .rl_css()

  # --- CDN scripts para ExcelJS y PapaParse ---
  cdn_scripts <- tagList(
    tags$script(src = "https://cdnjs.cloudflare.com/ajax/libs/exceljs/4.4.0/exceljs.min.js"),
    tags$script(src = "https://cdnjs.cloudflare.com/ajax/libs/PapaParse/5.4.1/papaparse.min.js")
  )

  # --- JavaScript: toggle filtros, limpiar filtros, exportar Excel ---
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
        console.error('Error exportando:', e);
      }
    }
  ", element_id, element_id, element_id, element_id)))

  # --- Barra de busqueda + botones de accion ---
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

  # --- Panel de filtros crosstalk + numericos ---
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

  # --- Filtros numericos (misma fila, mismos estilos) ---
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

  # --- Columnas ---
  cols <- .rl_build_columns(max_abs_lfc, alpha, has_assay, single_assay,
                             show_missing_cols, has_description, has_quant_pepts)

  # --- Tabla reactable con SharedData ---
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
    cdn_scripts,
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

#' Cargar y validar archivo Protein_ID
#' @param input Data frame o ruta a TSV/CSV/Parquet con formato Protein_ID
#' @return Data frame validado
#' @noRd
.pl_load_data <- function(input) {
  if (is.character(input) && length(input) == 1) {
    if (!file.exists(input)) stop("Archivo no encontrado: ", input)
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
        stop("El paquete 'arrow' es necesario para leer archivos Parquet")
      }
      input <- arrow::read_parquet(input)
    } else {
      stop("Formato no soportado: ", ext, ". Usa TSV, CSV o Parquet.")
    }
  }

  df <- as.data.frame(input)

  if (!"PG.ProteinGroups" %in% names(df)) {
    stop("Columna requerida 'PG.ProteinGroups' no encontrada. ",
         "El archivo debe ser un Protein_ID exportado por preprocess_spectronaut().")
  }
  if (!any(grepl("^PG\\.NrOfPrecursorsIdentified_", names(df)))) {
    stop("No se detectaron columnas de muestra (PG.NrOfPrecursorsIdentified_*). ",
         "Verifica que el archivo sea un Protein_ID válido.")
  }
  df
}


#' Parsear nombres de columnas Protein_ID en (metric, condition, replicate)
#' @param colnames_vec Vector de nombres de columnas del data frame
#' @return Data frame ordenado: column, metric, condition, replicate, coding
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


#' Paleta de colores por condicion
#' @param conditions Vector de codigos de condicion (e.g., c("A","B","C","D"))
#' @return Vector nombrado hex
#' @noRd
.pl_condition_palette <- function(conditions) {
  defaults <- c(
    A = "#4F81BD", B = "#9BBB59", C = "#F79646", D = "#8064A2",
    E = "#4BACC6", F = "#C0504D", G = "#9F8A76", H = "#646464"
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


#' Construir colDefs para la tabla Protein_ID
#' @param sample_map Salida de .pl_parse_samples()
#' @param max_per_col Vector numerico nombrado con max por columna (para data bars)
#' @return Lista de colDef
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

  # Primera columna de cada metrica -> ancla para borde y filtro agregado
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
      "PG.Cscore.RunWise"  = "val.toFixed(3)",
      "Math.round(val)"
    )

    if (metric == "PG.Coverage") {
      cell_js <- JS(sprintf("function(cellInfo) {
        var val = cellInfo.value;
        if (val == null || isNaN(val)) {
          return '<div class=\"pl-bar-wrapper\"><span class=\"pl-bar-value pl-bar-empty\">–</span></div>';
        }
        var pct = Math.min(100, Math.max(0, val / %s * 100));
        var formatted = %s;
        return '<div class=\"pl-bar-wrapper\">' +
          '<div class=\"pl-bar pl-cond-%s\" style=\"width:' + pct.toFixed(1) + '%%\"></div>' +
          '<span class=\"pl-bar-value\">' + formatted + '</span>' +
          '</div>';
      }", max_val, fmt_js, cond))
    } else {
      cell_js <- JS(sprintf("function(cellInfo) {
        var val = cellInfo.value;
        if (val == null || isNaN(val)) return '<span class=\"pl-plain-value pl-bar-empty\">–</span>';
        return '<span class=\"pl-plain-value\">' + (%s) + '</span>';
      }", fmt_js))
    }

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


#' Construir columnGroups para la tabla Protein_ID
#' @param sample_map Salida de .pl_parse_samples()
#' @param static_cols Vector de nombres de columnas estaticas presentes
#' @return Lista de colGroup
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


#' Detail row expandible para Protein_ID
#' @return Objeto JS
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


#' Filtro numerico agregado OR sobre varias columnas
#' @description Genera un filterMethod que evalua la expresion numerica (>=, <=, >, <,
#'   =, !=) contra CUALQUIERA de las columnas indicadas: basta con que una cumpla
#'   para conservar la fila. Soporta operadores compuestos AND (coma) y OR (barra),
#'   igual que \code{.numeric_filter_method}.
#' @param col_ids Vector de nombres de columnas sobre las que aplicar OR.
#' @return Objeto JS para usar como filterMethod en colDef.
#' @noRd
.pl_agg_filter_method <- function(col_ids) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("'jsonlite' es necesario para filtros agregados.")
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
# Funcion principal Protein_ID
# =============================================================================

#' Tabla Reactable Interactiva para datos Protein_ID (post-Spectronaut)
#'
#' Genera una tabla reactable con cabeceras agrupadas de 2 niveles, columnas
#' estaticas fijas (sticky) a la izquierda y data bars coloreadas por condicion
#' dentro de cada celda numerica. Reproduce el formato Excel habitual y lo mejora.
#'
#' @param data Data frame o ruta a TSV/CSV/Parquet de Protein_ID (salida de
#'   preprocess_spectronaut() -> protein_id).
#' @param metadata Opcional: data frame de metadata (run_summary) con columna
#'   Coding para fijar el orden de las columnas de muestra.
#' @param page_size Filas por pagina (default 15).
#' @param height Altura de la tabla en px (default 720).
#' @param element_id ID del elemento Reactable (default NULL).
#' @param selection "single" | "multiple" | NULL.
#' @param searchable Habilitar busqueda (default TRUE).
#'
#' @return Objeto reactable.
#'
#' @examples
#' \dontrun{
#' res <- preprocess_spectronaut(
#'   file_path = "data/Curso_Q24_DIA_Spectronaut_v20_Report.tsv",
#'   condition_order = c("A","B","C","D")
#' )
#' protein_list_reactable(res$protein_id, metadata = res$metadata)
#' }
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
    stop("No se encontraron columnas de muestra en el data frame. ",
         "Formato esperado: PG.<metric>_<condition>_<replicate>")
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
        background  = "rgba(14, 102, 85, 0.9)",
        color       = "#ffffff",
        height      = "40px",
        display     = "flex",
        alignItems  = "center",
        justifyContent = "center"
      ),
      style = list(height = "40px", display = "flex", alignItems = "center")
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
    theme        = .rl_theme(),
    language     = .rl_lang(),
    details      = .pl_detail_row()
  )
}


# =============================================================================
# Wrapper standalone / Quarto para Protein_ID
# =============================================================================

#' Widget Completo para tabla Protein_ID con filtros, busqueda y export Excel
#'
#' Envuelve \code{protein_list_reactable()} con:
#'  - Barra de busqueda y botones (toggle filtros, limpiar, export a Excel).
#'  - Chips por condicion que ocultan/muestran los 16 columnas asociadas.
#'  - Filtro numerico sobre MW [Da] con operadores AND/OR.
#'  - Export a Excel con cabeceras de 2 niveles (grupos mergeados + sub-labels)
#'    y tinte por condicion, reproduciendo el formato manual.
#'
#' @inheritParams protein_list_reactable
#' @param element_id ID del elemento (default "protein_id_table").
#'
#' @return Objeto htmltools browsable.
#'
#' @examples
#' \dontrun{
#' protein_list_widget("data/Protein_ID_20260423_142504.tsv")
#' }
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
    stop("No se encontraron columnas de muestra en el data frame.")
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

  # --- Estructuras JSON para el cliente ---
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("El paquete 'jsonlite' es necesario para protein_list_widget().")
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

  # --- CSS + scripts CDN ---
  css <- .rl_css()
  cdn_scripts <- tagList(
    tags$script(src = "https://cdnjs.cloudflare.com/ajax/libs/exceljs/4.4.0/exceljs.min.js"),
    tags$script(src = "https://cdnjs.cloudflare.com/ajax/libs/PapaParse/5.4.1/papaparse.min.js")
  )

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
          cell.fill = { type: 'pattern', pattern: 'solid', fgColor: { argb: 'FF0E6655' } };
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
        console.error('Error exportando:', e);
      }
    }
  ", cond_struct_json, group_struct_json, palette_json,
      element_id, element_id, element_id, element_id, element_id, element_id)))

  # --- Barra de busqueda + botones ---
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

  # --- Panel de filtros ---
  cond_chips <- lapply(conditions, function(cc) {
    tags$span(
      class = "pl-cond-chip",
      `data-cond` = cc,
      style = sprintf("background-color:%s;", palette[[cc]]),
      onclick = sprintf("plToggleCondition('%s')", cc),
      cc
    )
  })

  # --- Anclas de filtro agregado (primera columna de cada metrica) ---
  .first_of <- function(mm) {
    cc <- sample_map$column[sample_map$metric == mm]
    if (length(cc) == 0) NA_character_ else cc[1]
  }
  metric_filters <- list(
    list(label = "# Uniq. PSMs Identified",  anchor = .first_of("PG.NrOfPrecursorsIdentified"),        placeholder = "≥ 2 ..."),
    list(label = "# Uniq. Pepts Identified", anchor = .first_of("PG.NrOfStrippedSequencesIdentified"), placeholder = "≥ 2 ..."),
    list(label = "Coverage [%]",             anchor = .first_of("PG.Coverage"),                        placeholder = "≥ 10 ..."),
    list(label = "Spectronaut Cscore",       anchor = .first_of("PG.Cscore.RunWise"),                  placeholder = "≥ 2 ...")
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
      tags$span(class = "rl-filter-hint", "OR entre las muestras del bloque")
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

  # --- Tabla ---
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
        background  = "rgba(14, 102, 85, 0.9)",
        color       = "#ffffff",
        height      = "40px",
        display     = "flex",
        alignItems  = "center",
        justifyContent = "center"
      ),
      style = list(height = "40px", display = "flex", alignItems = "center")
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
    theme        = .rl_theme(),
    language     = .rl_lang(),
    details      = .pl_detail_row()
  )

  browsable(tagList(
    css,
    cdn_scripts,
    js_code,
    search_actions,
    filters_panel,
    tbl
  ))
}


# =============================================================================
# Ejemplos de Uso (no ejecutar)
# =============================================================================
if (FALSE) {
  source("R/Preprocessing.R")
  source("R/Results_List_reactable.R")

  # Desde TSV directo
  protein_list_widget("data/Protein_ID_20260423_142504.tsv")

  # Desde el resultado de preprocess_spectronaut()
  res <- preprocess_spectronaut(
    file_path = "data/Curso_Q24_DIA_Spectronaut_v20_Report.tsv",
    condition_order = c("A", "B", "C", "D")
  )
  protein_list_widget(res$protein_id, metadata = res$metadata)

  # Version sin filtros (Shiny puede envolverla)
  protein_list_reactable(res$protein_id, element_id = "tabla_id")
}
