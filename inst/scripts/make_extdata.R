# =============================================================================
# Procedencia y preparación de inst/extdata/ y data/nadia_dia.rda
# =============================================================================
#
# Bioconductor exige documentar cómo se obtuvieron y prepararon los datos que
# acompañan al paquete. Este script es la receta completa y reproducible.
#
# Los ficheros de partida no se distribuyen con el paquete: son reports
# completos de entre 0,6 y 29 MB que viven en `data-raw/` del repositorio de
# desarrollo (https://github.com/sciordia/NADIA) y quedan fuera del tarball vía
# `.Rbuildignore`. Ejecuta este script desde la raíz del repositorio.
#
# Origen de los datos
# -------------------
# Todos proceden de experimentos de proteómica cuantitativa adquiridos en el
# Servicio de Proteómica del Centro Nacional de Biotecnología (CNB-CSIC).
#
# * DIA (Spectronaut): 16 inyecciones de un mismo digerido, repartidas en cuatro
#   condiciones (A-D) de cuatro réplicas cada una, adquiridas en un Orbitrap
#   Exploris y procesadas con Spectronaut v20. Report a nivel de grupo de
#   proteínas y run (formato largo).
# * DIA spike-in: el mismo report anotado con la especie de cada grupo de
#   proteínas, para el benchmarking basado en organismos añadidos.
# * TMT (Proteome Discoverer): experimento TMTpro de dos mezclas y diez
#   fracciones, buscado con tres motores.
# * LFQ (Proteome Discoverer): experimento sin marcaje, con su fichero de
#   anotación muestra -> condición.
#
# Recorte aplicado
# ----------------
# El objetivo es que `R CMD check` se mantenga muy por debajo de los 10 minutos
# con más de cien ejemplos ejecutándose, no ahorrar espacio: el report DIA
# completo ya cabría en el límite de tamaño. Se conservan tres de las cuatro
# condiciones (A, B y D) porque con solo dos el soft-clustering del Pattern
# Profiler queda degenerado —agrupa perfiles a lo largo de las condiciones, y con
# dos puntos por perfil no hay patrón que agrupar—, y las 2.000 proteínas con más
# valores cuantificados de las 10.437 originales.
# =============================================================================

CONDICIONES <- c("A", "B", "D")
N_PROTEINAS <- 2000L

RAW      <- "data-raw"
EXTDATA  <- "inst/extdata"
dir.create(EXTDATA, recursive = TRUE, showWarnings = FALSE)
dir.create("data",  recursive = TRUE, showWarnings = FALSE)

gz <- function(df, destino) {
  con <- gzfile(destino, "wt", compression = 9)
  on.exit(close(con), add = TRUE)
  utils::write.table(df, con, sep = "\t", row.names = FALSE, quote = FALSE,
                     na = "")
  message(sprintf("  %-32s %6.2f MB", basename(destino),
                  file.size(destino) / 1024^2))
}

# --- 1. Report DIA de Spectronaut -------------------------------------------
message("1. nadia_dia_report.tsv.gz")
dia <- utils::read.delim(file.path(RAW, "Curso_Q24_DIA_Spectronaut_v20_Report.tsv"),
                         sep = "\t", header = TRUE, check.names = FALSE,
                         stringsAsFactors = FALSE)
dia <- dia[dia$R.Condition %in% CONDICIONES, , drop = FALSE]

# Las 2.000 proteínas con más valores cuantificados en esas condiciones
cuantificadas <- table(dia$PG.ProteinGroups[!is.na(dia$PG.Quantity) &
                                              dia$PG.Quantity > 0])
elegidas <- names(sort(cuantificadas, decreasing = TRUE))[seq_len(N_PROTEINAS)]
dia <- dia[dia$PG.ProteinGroups %in% elegidas, , drop = FALSE]
gz(dia, file.path(EXTDATA, "nadia_dia_report.tsv.gz"))

# --- 2. Report DIA con especies (spike-in, para benchmarking) ---------------
message("2. nadia_dia_spikein.tsv.gz")
spk <- utils::read.delim(
  file.path(RAW, "Curso_Q24_DIA_Spectronaut_v20_Report_Species.tsv"),
  sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
gz(spk, file.path(EXTDATA, "nadia_dia_spikein.tsv.gz"))

# --- 3. Report TMT de Proteome Discoverer -----------------------------------
message("3. nadia_tmt_report.tsv.gz")
tmt <- utils::read.delim(
  file.path(RAW, "20260527_Q25_TMTpro_TMT1y2_10Fr_Static_3engines_onlyRAW.tsv"),
  sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
gz(tmt, file.path(EXTDATA, "nadia_tmt_report.tsv.gz"))

# --- 4. Report LFQ de Proteome Discoverer y su anotación --------------------
message("4. nadia_lfq_report.tsv.gz + nadia_lfq_annotation.tsv")
lfq <- utils::read.delim(
  file.path(RAW, "20260710_AIturrate_2659_LFQ_QUANT_onlyRAW.tsv"),
  sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
gz(lfq, file.path(EXTDATA, "nadia_lfq_report.tsv.gz"))

file.copy(file.path(RAW, "20260710_AIturrate_2659_LFQ_QUANT_onlyRAW_Annot.tsv"),
          file.path(EXTDATA, "nadia_lfq_annotation.tsv"), overwrite = TRUE)

# --- 5. Objeto preprocesado data/nadia_dia.rda ------------------------------
# Se genera desde el report ya recortado de inst/extdata, de modo que el objeto
# y el fichero de ejemplo describan exactamente el mismo experimento.
message("5. data/nadia_dia.rda")
nadia_dia <- NADIA::preprocess_spectronaut(
  file_path       = file.path(EXTDATA, "nadia_dia_report.tsv.gz"),
  condition_order = CONDICIONES,
  verbose         = FALSE)

save(nadia_dia, file = "data/nadia_dia.rda",
     compress = "xz", compression_level = 9)
message(sprintf("  %-32s %6.2f MB", "nadia_dia.rda",
                file.size("data/nadia_dia.rda") / 1024^2))
message(sprintf("  protein_quant: %d x %d | %d muestras | %d condiciones",
                nrow(nadia_dia$protein_quant), ncol(nadia_dia$protein_quant),
                nrow(nadia_dia$metadata),
                length(unique(nadia_dia$metadata$R.Condition))))
