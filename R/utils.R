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
#' La definición solo se aplica en R < 4.4; a partir de esa versión el operador
#' está en `base` y se usa directamente.
#'
#' @param a Valor a comprobar.
#' @param b Valor alternativo si `a` es `NULL`.
#' @return `a` si no es `NULL`; en caso contrario `b`.
#' @keywords internal
#' @noRd
if (!exists("%||%", mode = "function")) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
}


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
