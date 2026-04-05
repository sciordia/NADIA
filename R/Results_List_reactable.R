
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
    }
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
    .lfc-bar.positive { background-color: #E63946; }
    .lfc-bar.negative { background-color: #457B9D; }

    /* Barra porcentaje ausencia */
    .missing-bar-container {
      display: flex;
      align-items: center;
      gap: 4px;
      width: 100%;
    }
    .missing-value {
      min-width: 35px;
      text-align: right;
      font-variant-numeric: tabular-nums;
    }
    .missing-bar-wrapper {
      flex: 1;
      height: 12px;
      background: #f0f0f0;
      border-radius: 2px;
      overflow: hidden;
    }
    .missing-bar {
      height: 100%;
      border-radius: 2px;
    }

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

    /* Ocultar busqueda por defecto de reactable */
    .rt-search {
      display: none !important;
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
      display: grid;
      grid-template-columns: repeat(3, 1fr);
      gap: 1.5rem;
    }
    @media (max-width: 1200px) {
      .rl-filters-row { grid-template-columns: repeat(2, 1fr); }
    }
    @media (max-width: 768px) {
      .rl-filters-row { grid-template-columns: 1fr; }
    }
    .rl-filter-item {
      display: flex;
      flex-direction: column;
      gap: 0.25rem;
    }
    .rl-filter-label {
      font-weight: 600;
      color: #495057;
      font-size: 0.9rem;
      margin-bottom: 0.25rem;
    }

    /* Selectize dropdowns dentro de filtros */
    .rl-filter-item .selectize-input {
      border: 2px solid #e0e0e0;
      border-radius: 6px;
      padding: 8px 12px;
      font-size: 0.9rem;
      min-height: 38px;
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

    /* Pagination focus */
    .rl-table .rt-page-size-select:focus {
      border-color: #0E6655;
      box-shadow: 0 0 5px #0E6655;
      outline: none;
    }
  "))
}


#' Traduccion al espanol
#' @return Objeto reactableLang
#' @noRd
.rl_lang <- function() {
  reactableLang(
    searchPlaceholder = "Buscar...",
    pagePrevious      = "Anterior",
    pageNext          = "Siguiente",
    noData            = "No hay datos para mostrar",
    pageSizeOptions   = "Mostrar {rows}",
    pageInfo          = "{rowStart}\u2013{rowEnd} de {rows} Prote\u00ednas"
  )
}


#' Funcion de detalle para filas expandibles (JS renderer)
#' @param has_assay Logico, si el data frame tiene columna Assay
#' @return Objeto JS para el parametro details de reactable
#' @noRd
.rl_detail_row <- function(has_assay = TRUE) {
  assay_block <- if (has_assay) {
    "
    var assay = row['Assay'] || '';
    if (assay) {
      html += '<div class=\"detail-row\"><span class=\"detail-label\">M\\u00e9todo:</span> ' + assay + '</div>';
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
    html += '<div class=\"detail-row\"><span class=\"detail-label\">Prote\\u00ednas:</span> ' + links + '</div>';

    // Gene.Names completo
    var genes = row['Gene.Names'] || '';
    html += '<div class=\"detail-row\"><span class=\"detail-label\">Genes:</span> ' + genes + '</div>';

    // P-valor y FDR con precision completa
    var pval = row['P.Value'];
    var fdr = row['adj.P.Val'];
    html += '<div class=\"detail-row\"><span class=\"detail-label\">P-valor:</span> ' + (pval != null ? pval.toExponential(4) : '') + '</div>';
    html += '<div class=\"detail-row\"><span class=\"detail-label\">FDR:</span> ' + (fdr != null ? fdr.toExponential(4) : '') + '</div>';

    // Assay (si existe)
    %s

    html += '</div>';
    return React.createElement('div', { dangerouslySetInnerHTML: { __html: html } });
  }", assay_block))
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
#' # Con filtro de comparacion
#' tbl <- results_list_reactable(de_res, comparisons = c("B-A", "C-A"))
#'
#' # En Shiny
#' # output$tabla <- renderReactable({
#' #   results_list_reactable(data(), comparisons = input$comp, element_id = "tabla")
#' # })
results_list_reactable <- function(
    data,
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

  # Aviso si hay muchas filas sin filtro
  if (nrow(df) > 15000 && is.null(comparisons)) {
    message("Nota: ", format(nrow(df), big.mark = "."),
            " filas. Considera filtrar por 'comparisons' para mejor rendimiento.")
  }

  # --- Detectar columnas opcionales ---
  has_missing <- all(c("MissingGlobal", "MissingPCT1", "MissingPCT2") %in% names(df))
  show_missing <- show_missing && has_missing

  has_assay <- "Assay" %in% names(df)
  single_assay <- has_assay && length(unique(df$Assay)) == 1

  # --- Pre-calcular max abs logFC para escalar barras ---
  max_abs_lfc <- max(abs(df$logFC), na.rm = TRUE)
  if (max_abs_lfc == 0) max_abs_lfc <- 1

  # --- Definicion de columnas (JS renderers para rendimiento con +60K filas) ---
  cols <- list(

    # Gen (primer gen, negrita + badge count)
    Gene.Names = colDef(
      name = "Gen",
      minWidth = 140,
      html = TRUE,
      cell = JS("function(cellInfo) {
        var val = cellInfo.value || '';
        var genes = val.split(';').map(function(s) { return s.trim(); }).filter(Boolean);
        var first = genes[0] || val;
        if (genes.length > 1) {
          return '<strong>' + first + '</strong> <span class=\"protein-count\">+' + (genes.length - 1) + '</span>';
        }
        return '<strong>' + first + '</strong>';
      }"),
      style = list(alignItems = "center")
    ),

    # Comparacion
    Comparison = colDef(
      name = "Comparaci\u00f3n",
      width = 120,
      align = "center"
    ),

    # Cambio (badge)
    Change = colDef(
      name = "Cambio",
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
    ),

    # logFC (valor + barra)
    logFC = colDef(
      name = "log\u2082 FC",
      width = 160,
      align = "center",
      html = TRUE,
      cell = JS(sprintf("function(cellInfo) {
        var val = cellInfo.value;
        var maxLfc = %s;
        var pct = Math.abs(val) / maxLfc * 100;
        var barClass = val >= 0 ? 'lfc-bar positive' : 'lfc-bar negative';
        var formatted = val.toFixed(3);
        return '<div class=\"lfc-bar-container\">' +
          '<span class=\"lfc-value\">' + formatted + '</span>' +
          '<div class=\"lfc-bar-wrapper\">' +
          '<div class=\"' + barClass + '\" style=\"width:' + pct.toFixed(1) + '%%\"></div>' +
          '</div></div>';
      }", max_abs_lfc))
    ),

    # FDR (notacion cientifica, negrita si significativo)
    adj.P.Val = colDef(
      name = "FDR",
      width = 120,
      align = "right",
      html = TRUE,
      cell = JS(sprintf("function(cellInfo) {
        var val = cellInfo.value;
        if (val == null || isNaN(val)) return '';
        var formatted = val.toExponential(2);
        if (val < %s) {
          return '<strong style=\"color: #0E6655;\">' + formatted + '</strong>';
        }
        return formatted;
      }", alpha))
    ),

    # P-valor (oculto por defecto, visible en detalle)
    P.Value = colDef(
      name = "P-valor",
      width = 110,
      align = "right",
      show = FALSE,
      cell = JS("function(cellInfo) {
        var val = cellInfo.value;
        if (val == null || isNaN(val)) return '';
        return val.toExponential(2);
      }")
    ),

    # Protein.IDs (truncado + count badge)
    Protein.IDs = colDef(
      name = "Prote\u00ednas",
      minWidth = 160,
      html = TRUE,
      cell = JS("function(cellInfo) {
        var val = cellInfo.value || '';
        var pids = val.split(';').map(function(s) { return s.trim(); }).filter(Boolean);
        var first = pids[0] || val;
        if (pids.length > 1) {
          return first + ' <span class=\"protein-count\">+' + (pids.length - 1) + '</span>';
        }
        return first;
      }")
    )
  )

  # Assay: ocultar si hay un solo valor
  if (has_assay) {
    cols$Assay <- colDef(
      name = "M\u00e9todo",
      width = 130,
      show = !single_assay
    )
  }

  # Columnas de Missing% (JS renderer)
  if (show_missing) {
    .missing_col_js <- JS("function(cellInfo) {
      var pct = cellInfo.value;
      if (pct == null || isNaN(pct)) return '';
      var color = '#02905A';
      if (pct > 25) color = '#E63946';
      else if (pct > 0) color = '#ffa62d';
      return '<div class=\"missing-bar-container\">' +
        '<span class=\"missing-value\">' + pct + '%</span>' +
        '<div class=\"missing-bar-wrapper\">' +
        '<div class=\"missing-bar\" style=\"width:' + pct + '%; background-color:' + color + ';\"></div>' +
        '</div></div>';
    }")

    cols$MissingGlobal <- colDef(
      name = "% Ausencia",
      width = 120,
      align = "center",
      html = TRUE,
      cell = .missing_col_js
    )
    cols$MissingPCT1 <- colDef(
      name = "% Grupo 1",
      width = 110,
      align = "center",
      html = TRUE,
      cell = .missing_col_js
    )
    cols$MissingPCT2 <- colDef(
      name = "% Grupo 2",
      width = 110,
      align = "center",
      html = TRUE,
      cell = .missing_col_js
    )
  }

  # --- Construir reactable ---
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
      style = list(
        height     = "48px",
        display    = "flex",
        alignItems = "center"
      )
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
    details    = .rl_detail_row(has_assay)
  )
}


# =============================================================================
# Wrapper standalone / Quarto
# =============================================================================

#' Widget Completo con Filtros, Busqueda, Export y CSS
#'
#' Envuelve \code{results_list_reactable()} con filtros interactivos crosstalk
#' (Comparacion, Cambio, Metodo), campo de busqueda, boton de exportar a Excel
#' y CSS embebido. Ideal para documentos Quarto o uso interactivo en RStudio.
#' El resultado es browsable: al imprimirlo en consola se abre automaticamente
#' en el Viewer de RStudio o en el navegador.
#'
#' @inheritParams results_list_reactable
#' @param element_id ID del elemento (default: "deps_table")
#'
#' @return Objeto htmltools browsable (se muestra automaticamente en RStudio Viewer)
#'
#' @examples
#' # Uso standalone (se abre en RStudio Viewer con filtros)
#' results_list_widget("results/VolcanoPlot_Input_cycloess_Impseq_min.tsv")
#'
#' # Con filtro previo de comparacion (ademas de filtros interactivos)
#' results_list_widget(de_res, comparisons = "B-A")
results_list_widget <- function(
    data,
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

  # --- Carga y filtrado previo (misma logica que results_list_reactable) ---
  df <- .rl_load_data(data)

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

  # --- Pre-calcular max abs logFC ---
  max_abs_lfc <- max(abs(df$logFC), na.rm = TRUE)
  if (max_abs_lfc == 0) max_abs_lfc <- 1

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

    function rlClearFilters() {
      // Limpiar selectize (crosstalk)
      var selects = document.querySelectorAll('.rl-filter-item .selectized');
      selects.forEach(function(sel) {
        if (sel.selectize) sel.selectize.clear();
      });
      // Limpiar busqueda
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
          'Gene.Names': 'Gen',
          'Comparison': 'Comparaci\\u00f3n',
          'Change': 'Cambio',
          'logFC': 'log2 FC',
          'adj.P.Val': 'FDR',
          'P.Value': 'P-valor',
          'Protein.IDs': 'Prote\\u00ednas',
          'Assay': 'M\\u00e9todo',
          'MissingGlobal': '%%Ausencia',
          'MissingPCT1': '%%Grupo1',
          'MissingPCT2': '%%Grupo2'
        };

        var wb = new ExcelJS.Workbook();
        var ws = wb.addWorksheet('Resultados DE');

        // Columnas visibles (excluir las ocultas)
        var visibleHeaders = headers.filter(function(h) { return h !== 'P.Value'; });

        ws.columns = visibleHeaders.map(function(h) {
          return {
            header: headerMap[h] || h,
            key: h,
            width: (h === 'Protein.IDs' || h === 'Gene.Names') ? 30 : 15
          };
        });

        // Estilo header
        var headerRow = ws.getRow(1);
        headerRow.font = { bold: true, color: { argb: 'FFFFFFFF' }, size: 11 };
        headerRow.fill = { type: 'pattern', pattern: 'solid', fgColor: { argb: 'FF0E6655' } };
        headerRow.alignment = { vertical: 'middle', horizontal: 'center' };
        headerRow.height = 30;

        // Datos
        rows.forEach(function(row) {
          var rowData = {};
          visibleHeaders.forEach(function(h) { rowData[h] = row[h]; });
          var addedRow = ws.addRow(rowData);

          // Convertir numericos (solo columnas visibles)
          ['logFC', 'adj.P.Val', 'MissingGlobal', 'MissingPCT1', 'MissingPCT2'].forEach(function(col) {
            if (visibleHeaders.indexOf(col) === -1) return;
            var cell = addedRow.getCell(col);
            if (cell && cell.value) {
              var num = parseFloat(cell.value);
              if (!isNaN(num)) cell.value = num;
            }
          });
        });

        // Bordes
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
  ", element_id, element_id)))

  # --- Barra de busqueda + botones de accion ---
  search_actions <- div(class = "rl-search-actions",
    # Busqueda
    tags$input(
      type = "search",
      placeholder = "Buscar...",
      class = "rl-search-input",
      oninput = sprintf("Reactable.setSearch('%s', this.value)", element_id)
    ),
    # Botones
    div(class = "rl-action-buttons",
      # Toggle filtros
      tags$button(
        class = "rl-btn-action rl-btn-toggle-filters filters-hidden",
        onclick = "rlToggleFilters()",
        title = "Mostrar/Ocultar filtros",
        HTML('<svg xmlns="http://www.w3.org/2000/svg" width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><polygon points="22 3 2 3 10 12.46 10 19 14 21 14 12.46 22 3"></polygon></svg>')
      ),
      # Limpiar filtros
      tags$button(
        class = "rl-btn-action",
        onclick = "rlClearFilters()",
        title = "Limpiar filtros",
        HTML('<svg xmlns="http://www.w3.org/2000/svg" width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M3 6h18"></path><path d="M8 6V4c0-1 1-2 2-2h4c1 0 2 1 2 2v2"></path><path d="M19 6l-1 14c0 1-1 2-2 2H8c-1 0-2-1-2-2L5 6"></path><line x1="1" y1="1" x2="23" y2="23" stroke="#E63946" stroke-width="2"></line></svg>')
      ),
      # Exportar Excel
      tags$button(
        class = "rl-btn-action",
        onclick = "rlExportExcel()",
        title = "Exportar a Excel",
        HTML('<svg xmlns="http://www.w3.org/2000/svg" width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4"></path><polyline points="7 10 12 15 17 10"></polyline><line x1="12" y1="15" x2="12" y2="3"></line></svg>')
      )
    )
  )

  # --- Panel de filtros crosstalk ---
  filter_items <- list(
    div(class = "rl-filter-item",
      tags$label(class = "rl-filter-label", "Comparaci\u00f3n"),
      crosstalk::filter_select(
        id = "rl_filter_comparison",
        label = NULL,
        sharedData = shared_data,
        group = ~Comparison,
        multiple = TRUE
      )
    ),
    div(class = "rl-filter-item",
      tags$label(class = "rl-filter-label", "Cambio"),
      crosstalk::filter_select(
        id = "rl_filter_change",
        label = NULL,
        sharedData = shared_data,
        group = ~Change,
        multiple = TRUE
      )
    )
  )

  # Filtro Assay solo si hay mas de un valor
  if (has_assay && !single_assay) {
    filter_items <- c(filter_items, list(
      div(class = "rl-filter-item",
        tags$label(class = "rl-filter-label", "M\u00e9todo"),
        crosstalk::filter_select(
          id = "rl_filter_assay",
          label = NULL,
          sharedData = shared_data,
          group = ~Assay,
          multiple = TRUE
        )
      )
    ))
  }

  filters_panel <- div(class = "rl-filters-container",
    div(class = "rl-filters-row", filter_items)
  )

  # --- Construir columnas (replica la logica de results_list_reactable) ---
  cols <- list(
    Gene.Names = colDef(
      name = "Gen", minWidth = 140, html = TRUE,
      cell = JS("function(cellInfo) {
        var val = cellInfo.value || '';
        var genes = val.split(';').map(function(s) { return s.trim(); }).filter(Boolean);
        var first = genes[0] || val;
        if (genes.length > 1) {
          return '<strong>' + first + '</strong> <span class=\"protein-count\">+' + (genes.length - 1) + '</span>';
        }
        return '<strong>' + first + '</strong>';
      }"),
      style = list(alignItems = "center")
    ),
    Comparison = colDef(name = "Comparaci\u00f3n", width = 120, align = "center"),
    Change = colDef(
      name = "Cambio", width = 130, align = "center", html = TRUE,
      cell = JS("function(cellInfo) {
        var val = cellInfo.value;
        var cls = 'tag status-grey';
        if (val === 'Up') cls = 'tag status-green';
        else if (val === 'Down') cls = 'tag status-red';
        return '<span class=\"' + cls + '\">' + val + '</span>';
      }")
    ),
    logFC = colDef(
      name = "log\u2082 FC", width = 160, align = "center", html = TRUE,
      cell = JS(sprintf("function(cellInfo) {
        var val = cellInfo.value;
        var maxLfc = %s;
        var pct = Math.abs(val) / maxLfc * 100;
        var barClass = val >= 0 ? 'lfc-bar positive' : 'lfc-bar negative';
        var formatted = val.toFixed(3);
        return '<div class=\"lfc-bar-container\">' +
          '<span class=\"lfc-value\">' + formatted + '</span>' +
          '<div class=\"lfc-bar-wrapper\">' +
          '<div class=\"' + barClass + '\" style=\"width:' + pct.toFixed(1) + '%%\"></div>' +
          '</div></div>';
      }", max_abs_lfc))
    ),
    adj.P.Val = colDef(
      name = "FDR", width = 120, align = "right", html = TRUE,
      cell = JS(sprintf("function(cellInfo) {
        var val = cellInfo.value;
        if (val == null || isNaN(val)) return '';
        var formatted = val.toExponential(2);
        if (val < %s) {
          return '<strong style=\"color: #0E6655;\">' + formatted + '</strong>';
        }
        return formatted;
      }", alpha))
    ),
    P.Value = colDef(
      name = "P-valor", width = 110, align = "right", show = FALSE,
      cell = JS("function(cellInfo) {
        var val = cellInfo.value;
        if (val == null || isNaN(val)) return '';
        return val.toExponential(2);
      }")
    ),
    Protein.IDs = colDef(
      name = "Prote\u00ednas", minWidth = 160, html = TRUE,
      cell = JS("function(cellInfo) {
        var val = cellInfo.value || '';
        var pids = val.split(';').map(function(s) { return s.trim(); }).filter(Boolean);
        var first = pids[0] || val;
        if (pids.length > 1) {
          return first + ' <span class=\"protein-count\">+' + (pids.length - 1) + '</span>';
        }
        return first;
      }")
    )
  )

  if (has_assay) {
    cols$Assay <- colDef(name = "M\u00e9todo", width = 130, show = !single_assay)
  }

  if (show_missing_cols) {
    .missing_js <- JS("function(cellInfo) {
      var pct = cellInfo.value;
      if (pct == null || isNaN(pct)) return '';
      var color = '#02905A';
      if (pct > 25) color = '#E63946';
      else if (pct > 0) color = '#ffa62d';
      return '<div class=\"missing-bar-container\">' +
        '<span class=\"missing-value\">' + pct + '%</span>' +
        '<div class=\"missing-bar-wrapper\">' +
        '<div class=\"missing-bar\" style=\"width:' + pct + '%; background-color:' + color + ';\"></div>' +
        '</div></div>';
    }")
    cols$MissingGlobal <- colDef(name = "% Ausencia", width = 120, align = "center", html = TRUE, cell = .missing_js)
    cols$MissingPCT1   <- colDef(name = "% Grupo 1", width = 110, align = "center", html = TRUE, cell = .missing_js)
    cols$MissingPCT2   <- colDef(name = "% Grupo 2", width = 110, align = "center", html = TRUE, cell = .missing_js)
  }

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
    details    = .rl_detail_row(has_assay)
  )

  # --- Ensamblar widget ---
  browsable(tagList(
    css,
    cdn_scripts,
    js_code,
    search_actions,
    filters_panel,
    tbl
  ))
}
