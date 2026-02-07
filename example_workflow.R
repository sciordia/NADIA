# =============================================================================
# Example Workflow: Complete Proteomics Analysis Pipeline
# =============================================================================
#
# This script demonstrates the full analysis flow from raw Spectronaut data
# to final visualizations, integrating all modules in this project.
#
# Author: Sergio Ciordia
# License: MIT
# =============================================================================

# ===== 1. PREPROCESSING =====
# Parse and structure raw Spectronaut report data

source("R/Preprocessing.R")

preprocessing <- preprocess_spectronaut(
  file_path = "data/Spectronaut_Report.tsv",
  condition_order = c("A", "B", "C", "D"),
  export_dir = "./results"
)

# ===== 2. PROCESSING (coordinator) =====
# Loads Normalization.R, Imputation.R, DEAnalysis.R automatically

source("R/Processing.R")

result <- process_proteomics(
  preprocessing = preprocessing,
  export_dir = "./results",
  cyclic_loess_method = "fast",
  cyclic_loess_iterations = 3,
  cyclic_loess_span = 0.7,
  control = "A",
  alpha = 0.05,
  export_format = "both"
)

# Access results
print(result)
result$se_proc        # Processed SummarizedExperiment
result$DEPs_results   # Differential expression results

# ===== 2b. STANDALONE MODULE USAGE (alternative) =====
# Each module can be sourced and used independently:
#
# source("R/Normalization.R")
# norm_result <- normalize_proteomics(
#   data = protein_data,
#   metadata = metadata,
#   cyclic_loess_method = "pairs",
#   cyclic_loess_iterations = 5
# )
#
# source("R/Imputation.R")
# imp_result <- impute_proteomics(
#   se = norm_result$se,
#   normalized_assay_name = "normalized",
#   imputed_assay_name = "Cycloess"
# )
#
# source("R/DEAnalysis.R")
# de_result <- de_analysis_proteomics(
#   se = imp_result$se,
#   assay_name = "Cycloess",
#   control = "A",
#   alpha = 0.05
# )

# ===== 3. BOXPLOT =====
# Interactive boxplots of intensity distributions

source("R/Boxplot_Highcharts_Final.R")

boxplot_data <- arrow::read_parquet("./results/BoxPlot_Input.parquet")

hc_boxplots <- boxplot_highchart_list(
  data = boxplot_data,
  intensity_col = "Intensity",
  group_col = "Column",
  assay_col = "Assay",
  condition_col = "Condition"
)

# Display first boxplot
hc_boxplots[[1]]

# ===== 4. VOLCANO PLOT =====
# Interactive volcano plots for differential expression

source("R/Volcano_Plot_Highcharts_Final.R")

hc_volcanos <- volcano_highchart_list(
  de_res = result$DEPs_results,
  logFC_col = "logFC",
  pval_col = "adj.P.Val",
  gene_col = "Gene.Names",
  comparison_col = "Comparison"
)

# Display first volcano plot
hc_volcanos[[1]]

# ===== 5. PCA =====
# Interactive PCA plots

source("R/PCA_Highcharts_Final.R")

pca_input <- arrow::read_parquet("./results/PCA_Input.parquet")

hc_pcas <- pca_highchart_list(
  pca_input = pca_input,
  sample_col = "SampleID",
  feature_col = "FeatureID",
  intensity_col = "Intensity",
  condition_col = "Condition"
)

# Display first PCA plot
hc_pcas[[1]]

# ===== 6. HEATMAP =====
# Static heatmaps with clustering

source("R/Heatmap_tidyHeatmap.R")

hm <- proteomics_heatmap(
  data = pca_input,
  mode = "any"
)

# ===== 7. PATTERN PROFILER =====
# Pattern profiling analysis and visualization

source("R/Pattern_Profiler_Analysis.R")

pp_result <- pattern_profiler_analysis(
  se = result$se_proc,
  DEPs_results = result$DEPs_results
)

source("R/Pattern_Profiler_Highcharts.R")

hc_patterns <- pattern_profiler_highchart_list(
  data = pp_result
)

# Display first pattern plot
hc_patterns[[1]]
