# =============================================================================
# NADIA — Instalación de dependencias
# =============================================================================
#
# Sustituye a renv. Instala los paquetes que NADIA necesita desde CRAN y
# Bioconductor, saltándose los que ya están presentes.
#
# Uso:
#   source("install_dependencies.R")            # duras + opcionales (recomendado)
#   install_nadia_deps(optional = FALSE)        # solo las imprescindibles
#   install_nadia_deps(dry_run = TRUE)          # solo informa de lo que falta
#
# Las dependencias "duras" las cargan los módulos al hacer source(); sin ellas
# el módulo no arranca. Las "opcionales" están protegidas con requireNamespace()
# y solo hacen falta si se usa el método concreto que las invoca (por ejemplo
# `mice` únicamente si se elige imp_method = "mice").
# =============================================================================

# --- Dependencias imprescindibles --------------------------------------------

# Se corresponden con el campo Imports: de DESCRIPTION. `grid`, `methods`,
# `stats`, `tools` y `utils` no aparecen porque vienen con R.
.nadia_cran_hard <- c(
  "dplyr", "tidyr", "tidyselect", "tibble", "stringr", "readr", "rlang",
  "magrittr", "ggplot2", "highcharter", "htmlwidgets", "paletteer",
  "RColorBrewer", "reactable", "htmltools", "jsonlite"
)

.nadia_bioc_hard <- c(
  "SummarizedExperiment", "S4Vectors"
)

# --- Dependencias opcionales (guardadas con requireNamespace) -----------------

# Se corresponden con el campo Suggests: de DESCRIPTION.
.nadia_cran_soft <- c(
  "arrow",       # lectura y escritura de Parquet
  "scales",      # escalas de color en las tablas reactable
  "crosstalk",   # enlazado entre widgets
  "pROC",        # AUC y pAUC en benchmarking
  "vegan",       # métrica PSS de imputación
  "cluster",     # métricas de agrupamiento en normalización
  "MASS",        # normalización Rlr
  "mice",        # imputación mice
  "missForest",  # imputación missForest
  "norm",        # imputación MLE
  "rrcovNA",     # imputación Impseq / Impseqrob
  "imputeLCMD",  # imputación QRILC / MinProb
  "lme4",        # descomposición de varianza PVCA
  "circlize",    # paletas continuas en heatmaps
  "tidyHeatmap", # heatmaps estáticos
  "e1071",       # c-means; Mfuzz lo necesita adjuntado
  "ggsci",       # paletas que sirve paletteer
  "viridis",     # paletas que sirve paletteer
  "testthat",    # tests
  "knitr",       # viñeta
  "rmarkdown"    # viñeta
)

.nadia_bioc_soft <- c(
  "limma",       # expresión diferencial y normalización cyclic loess / quantile
  "limpa",       # imputación probabilística + DE (de_method = "limpa")
  "impute",      # imputación knn
  "pcaMethods",  # imputación bpca
  "vsn",         # normalización vsn
  "BERT",        # corrección de lote (HarmonizR)
  "Biobase",     # ExpressionSet; Mfuzz lo necesita adjuntado
  "Mfuzz",       # soft-clustering del Pattern Profiler
  "ComplexHeatmap", # motor de los heatmaps
  "BiocStyle",   # estilo de la viñeta
  "MSstats"      # solo para los scripts 7_MSstats_Workflow*.R
)

# --- Motor -------------------------------------------------------------------

.nadia_missing <- function(pkgs) {
  pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
}

# En sesiones no interactivas (Rscript) el mirror de CRAN no está fijado y
# install.packages() aborta con "trying to use CRAN without setting a mirror".
.nadia_set_repo <- function() {
  repo <- getOption("repos")
  if (is.null(repo[["CRAN"]]) || is.na(repo[["CRAN"]]) ||
      identical(unname(repo[["CRAN"]]), "@CRAN@")) {
    repo[["CRAN"]] <- "https://cloud.r-project.org"
    options(repos = repo)
    message("  [previo] mirror de CRAN fijado en ", repo[["CRAN"]])
  }
  invisible(NULL)
}

.nadia_install <- function(pkgs, source, dry_run) {
  missing <- .nadia_missing(pkgs)

  if (length(missing) == 0L) {
    message("  [", source, "] todo presente (", length(pkgs), " paquetes)")
    return(character(0))
  }

  message("  [", source, "] faltan ", length(missing), ": ",
          paste(missing, collapse = ", "))
  if (dry_run) return(missing)

  if (identical(source, "CRAN")) {
    install.packages(missing)
  } else {
    BiocManager::install(missing, update = FALSE, ask = FALSE)
  }

  # Lo que siga faltando tras el intento de instalación
  .nadia_missing(missing)
}

#' Instala las dependencias de NADIA
#'
#' @param optional Si TRUE (por defecto) instala también los paquetes de los
#'   métodos opcionales. Si FALSE, solo lo imprescindible para cargar los módulos.
#' @param dry_run Si TRUE, informa de lo que falta sin instalar nada.
#' @return De forma invisible, el vector de paquetes que siguen sin instalarse.
install_nadia_deps <- function(optional = TRUE, dry_run = FALSE) {

  message("NADIA — comprobando dependencias (R ", getRversion(), ")")
  if (dry_run) message("MODO SIMULACIÓN: no se instalará nada\n")

  if (!dry_run) .nadia_set_repo()

  # BiocManager es el requisito previo para todo lo de Bioconductor
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    if (dry_run) {
      message("  [previo] falta BiocManager")
    } else {
      message("  [previo] instalando BiocManager")
      install.packages("BiocManager")
    }
  }

  failed <- character(0)

  message("\nImprescindibles:")
  failed <- c(failed, .nadia_install(.nadia_cran_hard, "CRAN", dry_run))
  failed <- c(failed, .nadia_install(.nadia_bioc_hard, "Bioconductor", dry_run))

  if (optional) {
    message("\nOpcionales:")
    failed <- c(failed, .nadia_install(.nadia_cran_soft, "CRAN", dry_run))
    failed <- c(failed, .nadia_install(.nadia_bioc_soft, "Bioconductor", dry_run))
  } else {
    message("\nOpcionales: omitidas (optional = FALSE)")
  }

  if (length(failed) > 0L && !dry_run) {
    warning("No se pudieron instalar: ", paste(failed, collapse = ", "),
            call. = FALSE)
  } else if (!dry_run) {
    message("\nListo. Carga los módulos con source(\"R/<Modulo>.R\").")
  }

  invisible(failed)
}

# Al hacer source() de este archivo se lanza la instalación completa.
# Comenta la línea siguiente si prefieres llamar a la función a mano.
install_nadia_deps()
