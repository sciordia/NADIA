# =============================================================================
# Utilidades internas compartidas
# =============================================================================
#
# Definiciones comunes a todos los módulos de NADIA. Antes cada módulo redefinía
# `%||%` por su cuenta (18 copias con 3 semánticas distintas), de modo que la
# última definición cargada ganaba y cambiaba el comportamiento de los módulos
# ya cargados. Ahora hay una sola definición canónica.
#
# Los módulos que necesitan los helpers de RNG cargan este archivo mediante el
# bloque `if (!exists(".rng_state", mode = "function"))` de su cabecera. Al
# convertir el proyecto en paquete ese bloque desaparece: todos los archivos de
# R/ comparten un mismo namespace.
#
# Author: Sergio Ciordia
# License: MIT
# =============================================================================


#' Operador de coalescencia nula
#'
#' Devuelve `a` salvo que sea `NULL`, en cuyo caso devuelve `b`. Misma semántica
#' que `base::\%||\%` (disponible desde R 4.4) y que `rlang::\%||\%`.
#'
#' Se define aquí en lugar de importarlo de `base` para que el paquete siga
#' funcionando bajo el `Depends: R (>= 4.4)` declarado sin depender de que el
#' operador esté exportado en esa versión concreta. La definición es idéntica.
#'
#' @param a Valor a comprobar.
#' @param b Valor alternativo si `a` es `NULL`.
#' @return `a` si no es `NULL`; en caso contrario `b`.
#' @keywords internal
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a


#' Capturar el estado del generador de números aleatorios
#'
#' Las funciones que llaman a `set.seed()` alteran el RNG de la sesión del
#' usuario, de modo que el código que se ejecute después deja de ser
#' reproducible. Estos dos helpers permiten dejar el RNG como estaba:
#'
#' ```
#' old_rng <- .rng_state()
#' on.exit(.rng_restore(old_rng), add = TRUE)
#' set.seed(seed)
#' ```
#'
#' @return El valor de `.Random.seed`, o `NULL` si el RNG aún no se ha usado en
#'   esta sesión.
#' @keywords internal
#' @noRd
.rng_state <- function() {
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  } else {
    NULL
  }
}

#' Restaurar el estado del generador de números aleatorios
#'
#' @param state Valor devuelto previamente por `.rng_state()`.
#' @return `NULL`, de forma invisible. Se llama por su efecto secundario.
#' @keywords internal
#' @noRd
.rng_restore <- function(state) {
  if (is.null(state)) {
    # El RNG no se había inicializado antes de la llamada: se elimina la
    # semilla creada por set.seed() para no dejar rastro.
    if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  } else {
    assign(".Random.seed", state, envir = globalenv())
  }
  invisible(NULL)
}


# =============================================================================
# Helpers de color
# =============================================================================
#
# Estaban duplicados en Boxplot_Highcharts_Final.R, PCA_Highcharts_Final.R,
# Pattern_Profiler_Highcharts.R y Heatmap_tidyHeatmap.R (8 definiciones de 3
# funciones). Al vivir todos los módulos en el entorno global, la copia activa
# era la del último módulo cargado; en un paquete ganaría la del último archivo
# por orden alfabético. Se unifican aquí.
#
# Verificado antes de unificar: para entradas hexadecimales las tres copias de
# `hex_to_rgba` y las dos de `darken_hex` daban salidas idénticas a igual
# alpha/factor. Solo diferían en el valor por defecto, y los 7 sitios de llamada
# lo pasan siempre de forma explícita. Por eso aquí NO se declara default: una
# llamada sin el argumento debe fallar de forma visible en vez de tomar un valor
# arbitrario.


#' Normalizar un color hexadecimal
#'
#' Quita el `#` inicial, descarta el canal alpha si el color viene en formato
#' `RRGGBBAA` y devuelve `#RRGGBB` en mayúsculas.
#'
#' @param hex Color en formato hexadecimal, con o sin `#`.
#' @return Cadena `#RRGGBB` en mayúsculas.
#' @keywords internal
#' @noRd
.normalize_hex <- function(hex) {
  hex <- gsub("^#", "", hex)
  if (nchar(hex) == 8) {
    hex <- substr(hex, 1, 6)
  }
  paste0("#", toupper(hex))
}


#' Convertir un color hexadecimal a cadena `rgba()`
#'
#' @param hex Color hexadecimal.
#' @param alpha Opacidad entre 0 y 1. Sin valor por defecto a propósito: todos
#'   los sitios de llamada lo especifican.
#' @return Cadena `"rgba(r, g, b, a)"` apta para Highcharts.
#' @keywords internal
#' @noRd
.hex_to_rgba <- function(hex, alpha) {
  hex <- .normalize_hex(hex)
  rgb_vals <- grDevices::col2rgb(hex)
  sprintf("rgba(%d, %d, %d, %.2f)",
          rgb_vals[1], rgb_vals[2], rgb_vals[3], alpha)
}


#' Oscurecer un color hexadecimal
#'
#' @param hex Color hexadecimal.
#' @param factor Fracción de oscurecimiento entre 0 y 1. Sin valor por defecto a
#'   propósito: todos los sitios de llamada lo especifican.
#' @return Color `#RRGGBB` oscurecido.
#' @keywords internal
#' @noRd
.darken_hex <- function(hex, factor) {
  hex <- .normalize_hex(hex)
  rgb_vals <- grDevices::col2rgb(hex)
  rgb_dark <- pmax(0, rgb_vals * (1 - factor))
  sprintf("#%02X%02X%02X",
          round(rgb_dark[1]), round(rgb_dark[2]), round(rgb_dark[3]))
}


# =============================================================================
# Helpers de filtrado de features
# =============================================================================

#' Construir el nombre de la columna de p-valor ajustado de una comparación
#'
#' @param comparison Nombre de la comparación, p.ej. `"B-A"`.
#' @return Nombre de columna, p.ej. `"adjP_B-A"`.
#' @keywords internal
#' @noRd
.adjp_col <- function(comparison) {
  paste0("adjP_", comparison)
}


#' Obtener los IDs de features según el modo de filtrado
#'
#' Unifica las dos copias que había en `Heatmap_tidyHeatmap.R` y
#' `PCA_Highcharts_Final.R`. Diferían en dos cosas:
#'
#' * El nombre del tercer modo: `"target"` en la del heatmap y `"specific"` en la
#'   del PCA. Aquí se aceptan **ambos** como sinónimos, porque las funciones
#'   públicas de cada módulo propagan su propio vocabulario y cualquiera de las
#'   dos habría roto a la otra al fusionarse en un único namespace.
#' * El filtrado de `mode = "any"`: la del PCA usaba `sig_any == TRUE`, que cuela
#'   `FeatureID` `NA` cuando `sig_any` tiene `NA`. Se conserva el `which()` de la
#'   del heatmap, que no lo hace.
#'
#' @param data Data frame en formato long con columnas `FeatureID`, `sig_any` y
#'   `adjP_*`.
#' @param mode `"all"` (todos), `"any"` (significativo en alguna comparación) o
#'   `"target"`/`"specific"` (significativo en `comparison`).
#' @param alpha Umbral de significancia para el modo dirigido.
#' @param comparison Comparación a usar en el modo dirigido.
#' @return Vector de `FeatureID` que cumplen el criterio.
#' @keywords internal
#' @noRd
.get_feature_ids <- function(data,
                             mode = c("all", "any", "target", "specific"),
                             alpha = 0.05,
                             comparison = NULL) {

  mode <- match.arg(mode)
  data <- as.data.frame(data)

  feat <- data[!duplicated(data$FeatureID), , drop = FALSE]

  if (mode == "all") {
    return(feat$FeatureID)
  }

  if (mode == "any") {
    if (!("sig_any" %in% names(feat))) {
      stop("La columna 'sig_any' es requerida para mode = 'any'")
    }
    # which() evita colar FeatureID NA cuando sig_any tiene NA (a diferencia de
    # feat$sig_any == TRUE, que devolvería filas NA).
    return(feat$FeatureID[which(feat$sig_any)])
  }

  # modo dirigido: "target" y "specific" son sinónimos
  if (is.null(comparison)) {
    stop("El argumento 'comparison' es requerido para mode = '", mode, "'")
  }

  col <- .adjp_col(comparison)
  if (!(col %in% names(feat))) {
    stop("No existe la columna: ", col)
  }

  feat$FeatureID[which(feat[[col]] <= alpha)]
}


# =============================================================================
# Dependencias de Mfuzz que exigen estar adjuntadas
# =============================================================================
#
# `Mfuzz` declara `Depends: Biobase, e1071` y llama a `exprs()` y `cmeans()` sin
# cualificar, así que sus funciones solo resuelven esos nombres si ambos paquetes
# están en la ruta de búsqueda. `requireNamespace()` no basta: carga el namespace
# pero no lo adjunta. Antes lo conseguía el `library(Mfuzz)` de la cabecera del
# módulo, que arrastraba sus Depends; un paquete no puede hacer eso.
#
# La solución es adjuntarlos solo mientras se ejecuta el Pattern Profiler y
# dejar la ruta de búsqueda como estaba, para no alterar la sesión del usuario.

#' Adjuntar las dependencias que Mfuzz necesita en la ruta de búsqueda
#'
#' @return Vector con los paquetes que esta llamada ha adjuntado (posiblemente
#'   vacío si ya lo estaban). Debe pasarse a `.mfuzz_deps_detach()`.
#' @keywords internal
#' @noRd
.mfuzz_deps_attach <- function() {
  necesarios <- c("Mfuzz", "Biobase", "e1071")
  faltan <- necesarios[!vapply(necesarios, requireNamespace, logical(1),
                               quietly = TRUE)]
  if (length(faltan) > 0) {
    stop("El Pattern Profiler necesita ", paste(faltan, collapse = ", "),
         ". Instálalos con BiocManager::install(c(",
         paste(sprintf('"%s"', faltan), collapse = ", "), ")).",
         call. = FALSE)
  }
  # Solo Biobase y e1071 hacen falta adjuntados; Mfuzz se usa con Mfuzz::
  adjuntables <- c("Biobase", "e1071")
  ya_estaban <- paste0("package:", adjuntables) %in% search()
  for (p in adjuntables[!ya_estaban]) attachNamespace(asNamespace(p))
  adjuntables[!ya_estaban]
}

#' Deshacer lo que hizo `.mfuzz_deps_attach()`
#'
#' @param pkgs Vector devuelto por `.mfuzz_deps_attach()`.
#' @return `NULL`, de forma invisible. Se llama por su efecto secundario.
#' @keywords internal
#' @noRd
.mfuzz_deps_detach <- function(pkgs) {
  for (p in pkgs) {
    nombre <- paste0("package:", p)
    if (nombre %in% search()) {
      detach(nombre, character.only = TRUE, unload = FALSE)
    }
  }
  invisible(NULL)
}


# =============================================================================
# Helpers de Proteome Discoverer
# =============================================================================
#
# Estaban duplicados entre Preprocessing_LFQ.R y Preprocessing_TMT.R. Las dos
# copias de `.parse_gene_from_description` eran idénticas; las de
# `.validate_pd_columns` solo diferían en el texto del mensaje de error.

#' Extraer el Gene Name de la columna Description (formato UniProt embebido)
#'
#' Busca el patrón `GN=<gene>` habitual en los exports de Proteome Discoverer.
#'
#' @param x Vector de descripciones.
#' @return Vector de nombres de gen, `NA` donde no haya coincidencia.
#' @keywords internal
#' @noRd
.parse_gene_from_description <- function(x) {
  m <- stringr::str_match(x, "GN=([^ ]+)")
  m[, 2]
}

#' Validar las columnas mínimas de un export de Proteome Discoverer
#'
#' @param df Data frame leído del export.
#' @return `TRUE` de forma invisible; aborta con un error si falta alguna
#'   columna requerida.
#' @keywords internal
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
