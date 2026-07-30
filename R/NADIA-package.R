# =============================================================================
# Documentación del paquete e importaciones
# =============================================================================
#
# Todas las directivas `@importFrom` del paquete viven aquí. El código ya estaba
# mayoritariamente cualificado con `pkg::` (976 llamadas, entre ellas las 528 de
# ggplot2 y las 102 de SummarizedExperiment, que están al 100 %), así que solo
# hace falta importar los nombres que se usan sin cualificar. Concentrarlos en un
# único bloque evita reescribir 450 sitios de llamada.
#
# Los paquetes que solo hacen falta para un método concreto están en `Suggests:`
# y se comprueban en el punto de uso con `requireNamespace()`; sus llamadas van
# siempre cualificadas con `pkg::`, nunca importadas aquí.
#
# Author: Sergio Ciordia
# License: MIT
# =============================================================================

#' NADIA: análisis de proteómica DIA consciente de los valores ausentes
#'
#' Pipeline completo para el análisis de expresión diferencial de proteínas, con
#' especial atención al tratamiento de los valores ausentes característicos de la
#' adquisición independiente de datos. Cubre el recorrido íntegro desde el report
#' crudo hasta la figura interactiva.
#'
#' El paquete se organiza en cuatro bloques:
#'
#' \describe{
#'   \item{Preprocesado}{[preprocess_spectronaut()] para Spectronaut/DIA-NN,
#'     [preprocess_tmt()] y [preprocess_lfq()] para Proteome Discoverer. Los tres
#'     devuelven el mismo objeto S3 `proteomics_data`, de modo que el resto del
#'     pipeline los consume sin cambios.}
#'   \item{Procesado}{[process_proteomics()] coordina normalización (13 métodos),
#'     corrección de lote opcional, imputación (19 métodos, incluidos los híbridos
#'     MAR/MNAR) y expresión diferencial con limma o limpa.}
#'   \item{Métricas y benchmarking}{[normalization_metrics()],
#'     [imputation_metrics()], [benchmarking_proteomics()] y
#'     [benchmarking_multiple()] evalúan y ordenan las combinaciones de métodos.}
#'   \item{Visualización}{gráficos interactivos con Highcharts
#'     ([boxplot_highchart_list()], [volcano_highchart_list()],
#'     [pca_highchart_list()]), estáticos con ggplot2 y ComplexHeatmap
#'     ([proteomics_heatmap()]), y tablas con reactable
#'     ([results_list_reactable()]).}
#' }
#'
#' @name NADIA-package
#' @aliases NADIA
#' @keywords internal
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data :=
#'
#' @importFrom dplyr across arrange bind_cols bind_rows case_when distinct filter group_by if_else inner_join left_join mutate n n_distinct pull rename row_number select slice slice_head summarise ungroup
#' @importFrom tidyr pivot_longer pivot_wider separate_rows unite
#' @importFrom tidyselect all_of any_of
#'
#' @importFrom stats aggregate approx as.dist as.formula cmdscale coef complete.cases cor cov dist ecdf IQR lm mad median model.matrix na.omit p.adjust prcomp predict quantile rnorm runif sd setNames var weights
#' @importFrom utils combn head modifyList read.csv read.delim tail write.table
#' @importFrom methods new
#' @importFrom grDevices col2rgb colorRampPalette
#'
#' @importFrom highcharter highchart hc_add_series hc_chart hc_colors hc_exporting hc_legend hc_plotOptions hc_subtitle hc_title hc_tooltip hc_xAxis hc_yAxis hcaes
#' @importFrom htmltools browsable div span tagList tags HTML
#' @importFrom htmlwidgets JS
#' @importFrom reactable colDef colGroup reactable reactableLang reactableTheme
"_PACKAGE"

# Nota: `JS` lo reexportan tanto highcharter como htmlwidgets. Se importa solo de
# htmlwidgets para no provocar un conflicto de nombres en el NAMESPACE.


# =============================================================================
# Nombres de columna usados con evaluación no estándar
# =============================================================================
#
# Los verbos de dplyr/tidyr y las estéticas de ggplot2 reciben nombres de columna
# sin entrecomillar, así que `R CMD check` los toma por variables globales no
# definidas y emite una NOTE por cada uno (327 en total, 78 nombres distintos).
# Declararlos aquí es la solución habitual; la alternativa —prefijar los 327
# sitios con `.data$`— no cambia el comportamiento y sí el riesgo.
utils::globalVariables(c(
  ".val_num", "ACC_OI", "AUC", "Assay", "Category", "Change", "Cluster",
  "Coding", "Column", "Comparison", "Condition", "Correlation", "Count",
  "FeatureID", "Gene.Names", "Intensity", "Label", "MDS1", "MDS1_VarPct",
  "MDS2", "Membership", "Method", "Metric", "NRMSE", "Normalization",
  "PC1", "PC1_VarPct", "PC2", "PCV", "PEV", "PG.Coverage",
  "PG.Coverage.Global", "PG.Cscore", "PG.Cscore.RunWise", "PG.Genes",
  "PG.MolecularWeight", "PG.NrOfPrecursorsIdentified",
  "PG.NrOfPrecursorsIdentified.Global",
  "PG.NrOfPrecursorsUsedForQuantification",
  "PG.NrOfStrippedSequencesIdentified",
  "PG.NrOfStrippedSequencesIdentified.Global",
  "PG.NrOfStrippedSequencesUsedForQuantification",
  "PG.ProteinDescriptions", "PG.ProteinGroups", "PG.Quantity", "PMAD",
  "PSS", "Percentage", "Protein", "Protein.IDs", "R.Condition",
  "R.FileName", "R.Replicate", "Rank", "Rank_Final", "Replicate", "SOR",
  "Sample", "SampleID", "Species", "Value", "Value_scaled", "adj.P.Val",
  "adjP", "category", "fill_hex", "label", "label_text", "logFC", "metric",
  "metric_coding", "minusLog10P", "n_clusters", "pct_label", "pval_fmt",
  "rank_final", "text_color", "value"
))
