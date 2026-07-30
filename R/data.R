# =============================================================================
# Documentación de los datos que acompañan al paquete
# =============================================================================

#' Experimento DIA de ejemplo, ya preprocesado
#'
#' Salida de [preprocess_spectronaut()] sobre un report de Spectronaut recortado,
#' lista para alimentar [process_proteomics()] y el resto del pipeline sin tener
#' que volver a leer el fichero de partida.
#'
#' @format Objeto S3 de clase `c("spectronaut_data", "proteomics_data", "list")`
#'   con tres elementos:
#'   \describe{
#'     \item{`metadata`}{Data frame de 12 filas y 7 columnas, una por inyección:
#'       `R.FileName`, `R.Condition`, `R.Replicate`, `Coding` y los tres
#'       recuentos de identificaciones que reporta Spectronaut
#'       (`R.PrecursorsIdentified`, `R.StrippedSequencesIdentified`,
#'       `R.ProteinGroupsIdentified`).}
#'     \item{`protein_id`}{Data frame de 2.000 filas por 56 columnas con la
#'       información de identificación de cada grupo de proteínas.}
#'     \item{`protein_quant`}{Data frame de 2.000 filas por 44 columnas: ocho de
#'       anotación más, por cada una de las 12 muestras, el número de precursores
#'       y de secuencias usados para cuantificar y la intensidad
#'       (`PG.Quantity_<condición>_<réplica>`).}
#'   }
#'
#' @details
#' Procede de un experimento de proteómica cuantitativa adquirido en un Orbitrap
#' Exploris en el Servicio de Proteómica del Centro Nacional de Biotecnología
#' (CNB-CSIC) y procesado con Spectronaut v20. El experimento original tiene 16
#' inyecciones del mismo digerido repartidas en cuatro condiciones (A-D) de
#' cuatro réplicas, y 10.437 grupos de proteínas.
#'
#' Para que los ejemplos del paquete se ejecuten con rapidez se conservan tres
#' condiciones (A, B y D) y las 2.000 proteínas con más valores cuantificados. Se
#' mantienen tres y no dos condiciones porque el soft-clustering de
#' [pattern_profiler_analysis()] agrupa perfiles a lo largo de las condiciones, y
#' con solo dos puntos por perfil no hay patrón que agrupar.
#'
#' La receta completa está en `system.file("scripts", "make_extdata.R", package =
#' "NADIA")`. El report recortado del que sale este objeto se distribuye también
#' como fichero, en `system.file("extdata", "nadia_dia_report.tsv.gz", package =
#' "NADIA")`, para poder ejecutar [preprocess_spectronaut()] de principio a fin.
#'
#' @source Servicio de Proteómica, Centro Nacional de Biotecnología (CNB-CSIC).
#'
#' @examples
#' data(nadia_dia)
#' nadia_dia
#'
#' # Las tres condiciones y sus cuatro réplicas
#' table(nadia_dia$metadata$R.Condition)
#'
#' # Proteínas y columnas de intensidad
#' dim(nadia_dia$protein_quant)
#' grep("^PG.Quantity_", colnames(nadia_dia$protein_quant), value = TRUE)
"nadia_dia"
