
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


#' CSS embebido para badges, barras y detalle
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
      margin-top: 1rem;
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

    /* Contenedor de busqueda */
    .rl-search-container {
      font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif;
      font-weight: 400;
      font-size: 14px;
      border: 2px solid #e0e0e0;
      border-radius: 8px;
      transition: all 0.3s ease;
      background-color: #ffffff;
      box-shadow: 0 2px 4px rgba(0, 0, 0, 0.05);
    }
    .rl-search-container:focus {
      border-color: #0E6655;
      box-shadow: 0 0 8px rgba(14, 102, 85, 0.3);
      outline: none;
      background-color: #fafafa;
    }
    .rl-search-container:hover {
      border-color: #0E6655;
      box-shadow: 0 2px 8px rgba(14, 102, 85, 0.1);
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


#' Funcion de detalle para filas expandibles
#' @param df Data frame con los datos (el mismo pasado a reactable)
#' @return Funcion para el parametro details de reactable
#' @noRd
.rl_detail_row <- function(df) {
  function(index) {
    row <- df[index, ]

    # Protein.IDs con enlaces UniProt
    pids <- trimws(unlist(strsplit(as.character(row$Protein.IDs), ";")))
    protein_links <- lapply(pids, function(pid) {
      tags$span(
        tags$a(
          href = paste0("https://www.uniprot.org/uniprot/", pid),
          target = "_blank",
          pid
        ),
        " "
      )
    })

    # Gene.Names completo
    genes <- as.character(row$Gene.Names)

    # Construir panel de detalle
    div(class = "rl-detail",
      div(class = "detail-row",
        span(class = "detail-label", "Prote\u00ednas:"),
        tagList(protein_links)
      ),
      div(class = "detail-row",
        span(class = "detail-label", "Genes:"),
        genes
      ),
      div(class = "detail-row",
        span(class = "detail-label", "P-valor:"),
        formatC(row$P.Value, format = "e", digits = 4)
      ),
      div(class = "detail-row",
        span(class = "detail-label", "FDR:"),
        formatC(row$adj.P.Val, format = "e", digits = 4)
      ),
      if ("Assay" %in% names(row)) {
        div(class = "detail-row",
          span(class = "detail-label", "M\u00e9todo:"),
          as.character(row$Assay)
        )
      }
    )
  }
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
    details    = .rl_detail_row(df)
  )
}


# =============================================================================
# Wrapper standalone / Quarto
# =============================================================================

#' Widget Completo con Busqueda y CSS para Uso Standalone/Quarto
#'
#' Envuelve \code{results_list_reactable()} con CSS embebido y campo de busqueda
#' externo. Ideal para documentos Quarto o uso interactivo en RStudio.
#'
#' @inheritParams results_list_reactable
#' @param element_id ID del elemento (default: "deps_table")
#'
#' @return Objeto htmltools tagList
#'
#' @examples
#' # Uso standalone
#' widget <- results_list_widget("results/VolcanoPlot_Input_cycloess_Impseq_min.tsv")
#'
#' # Con filtro
#' widget <- results_list_widget(de_res, comparisons = "B-A")
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

  # CSS embebido
  css <- .rl_css()

  # Campo de busqueda externo
  search_input <- div(
    style = "margin-bottom: 0.75rem",
    tags$input(
      type = "search",
      placeholder = "Buscar...",
      class = "rl-search-container",
      style = "padding: 8px 12px; width: 100%; margin-bottom: 10px; border-radius: 6px;",
      oninput = sprintf("Reactable.setSearch('%s', this.value)", element_id)
    )
  )

  # Tabla reactable
  tbl <- results_list_reactable(
    data         = data,
    comparisons  = comparisons,
    ain          = ain,
    alpha        = alpha,
    lfc_thr      = lfc_thr,
    page_size    = page_size,
    height       = height,
    show_missing = show_missing,
    element_id   = element_id,
    selection    = selection,
    searchable   = searchable
  )

  tagList(css, search_input, tbl)
}
