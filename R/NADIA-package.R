# =============================================================================
# Package documentation and imports
# =============================================================================
#
# Every `@importFrom` directive in the package lives here. The code was already
# largely qualified with `pkg::` (976 calls, among them the 528 ggplot2 and the
# 102 SummarizedExperiment ones, which are qualified 100 %), so only the names
# used unqualified need to be imported. Concentrating them in a single block
# avoids rewriting 450 call sites.
#
# Packages needed only for one specific method are listed in `Suggests:` and
# checked at the point of use with `requireNamespace()`; their calls are always
# qualified with `pkg::`, never imported here.
#
# Author: Sergio Ciordia
# License: GPL-3
# =============================================================================

#' NADIA: Missing Value-Aware DIA Proteomics Analysis
#'
#' A complete pipeline for differential protein expression analysis, with
#' particular attention to the missing values that characterise
#' data-independent acquisition. It covers the whole journey from the raw report
#' to the interactive figure.
#'
#' The package is organised in four blocks:
#'
#' \describe{
#'   \item{Preprocessing}{[preprocess_spectronaut()] for Spectronaut/DIA-NN,
#'     [preprocess_tmt()] and [preprocess_lfq()] for Proteome Discoverer. All
#'     three return the same `proteomics_data` S3 object, so the rest of the
#'     pipeline consumes them unchanged.}
#'   \item{Processing}{[process_proteomics()] coordinates normalization (13
#'     methods), optional batch correction, imputation (19 methods, including
#'     the MAR/MNAR hybrids) and differential expression with limma or limpa.}
#'   \item{Metrics and benchmarking}{[normalization_metrics()],
#'     [imputation_metrics()], [benchmarking_proteomics()] and
#'     [benchmarking_multiple()] evaluate and rank method combinations.}
#'   \item{Visualization}{interactive plots with Highcharts
#'     ([boxplot_highchart_list()], [volcano_highchart_list()],
#'     [pca_highchart_list()]), static ones with ggplot2 and ComplexHeatmap
#'     ([proteomics_heatmap()]), and tables with reactable
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

# Note: `JS` is re-exported by both highcharter and htmlwidgets. It is imported
# only from htmlwidgets so as not to cause a name clash in the NAMESPACE.


# =============================================================================
# Column names used with non-standard evaluation
# =============================================================================
#
# The dplyr/tidyr verbs and the ggplot2 aesthetics receive unquoted column
# names, so `R CMD check` takes them for undefined global variables and emits a
# NOTE for each one (327 in total, 78 distinct names). Declaring them here is
# the usual solution; the alternative -- prefixing the 327 sites with `.data$`
# -- does not change the behaviour but does add risk.
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
