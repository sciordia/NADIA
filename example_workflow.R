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

library(NADIA)

preprocessing <- preprocess_spectronaut(
  file_path = "data-raw/Curso_Q24_DIA_Spectronaut_v20_Report.tsv",
  condition_order = c("A", "B", "C", "D"),
  export_dir = "./results"
)

# ----- Alternativa: datos TMT (Proteome Discoverer) -----
# El resto del pipeline (Processing, Boxplot, Volcano, PCA, etc.) funciona
# sin cambios al recibir el output de preprocess_tmt().
#
# source("R/Preprocessing_TMT.R")
#
# preprocessing <- preprocess_tmt(
#   file_path = "data-raw/20260527_Q25_TMTpro_TMT1y2_10Fr_Static_3engines_onlyRAW.tsv",
#   condition_order = c("A", "B", "C", "D", "IS"),  # omite "IS" para descartar Internal Standards
#   export_dir = "./results"
# )

# ----- Alternativa: datos LFQ (Proteome Discoverer) -----
# The experimental design (sample -> condition) comes from the _Annot file
# (columns Column, Condition, Experiment). Per-sample intensities and metrics
# muestra se mapean POR NOMBRE. El resto del pipeline funciona sin cambios.
#
# source("R/Preprocessing_LFQ.R")
#
# preprocessing <- preprocess_lfq(
#   file_path  = "data-raw/20260710_AIturrate_2659_LFQ_QUANT_onlyRAW.tsv",
#   annot_path = "data-raw/20260710_AIturrate_2659_LFQ_QUANT_onlyRAW_Annot.tsv",
#   condition_order = c("WT", "MUT"),  # opcional; si NULL se deriva del _Annot
#   export_dir = "./results"
# )

# ===== 2. PROCESSING (coordinator) =====
# Loads Normalization.R, Imputation.R, DEAnalysis.R automatically


result <- process_proteomics(
  preprocessing = preprocessing,
  norm_method = "quantile",
  export_dir = "./results",
  min_reps_filter = 3,
  cyclic_loess_method = "fast",
  cyclic_loess_iterations = 3,
  cyclic_loess_span = 0.7,
  control = NULL,
  alpha = 0.05,
  export_format = "both"
)

# Access results
print(result)
result$se_proc        # Processed SummarizedExperiment
result$DEPs_results   # Differential expression results

# ===== 2b. STANDALONE MODULE USAGE (alternative) =====
# Each module can be sourced and used independently:

# source("R/Normalization.R")
# 
# # Extraer los dos objetos que necesita normalize_proteomics
# metadata     <- .prepare_metadata(preprocessing)
# protein_data <- .prepare_protein_data(preprocessing)
# 
# #Run Normalization
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
#   normalized_assay_name = "cycloess",
#   imputed_assay_name = "ImpSeqRob_Min"
# )
#
# source("R/DEAnalysis.R")
# de_result <- de_analysis_proteomics(
#   se = imp_result$se,
#   assay_name = "ImpSeqRob_Min",
#   control = "A",
#   alpha = 0.05
# )

# ===== 3. BOXPLOT =====
# Interactive boxplots of intensity distributions


boxplot_data <- arrow::read_parquet("./results/BoxPlot_Input.parquet")

hc_boxplots <- boxplot_highchart_list(
  data        = boxplot_data,
  assays      = c("log2", "ImpSeqRob_Min"),
  color_by    = "Condition",
  palette  = "ggsci::category10_d3",
  group_order = c("A", "B", "C", "D"),
  title = "Boxplot: {assay}",
  box_width = 20,
  horizontal = FALSE
)

hc_boxplots[["log2"]]
hc_boxplots[["ImpSeqRob_Min"]]


# ===== 4. VOLCANO PLOT =====
# Interactive volcano plots for differential expression


# --- 4.1 Volcano Plots ---                                                                                          
volcano_plots <- volcano_highchart_list(                                                                             
  de_res = result$DEPs_results,                                                                                                           
  lfc_thr = 0,                                                                                             
  alpha = 0.05,
  point_size = 3,
  p_col = "adj.P.Val",
  show_top_genes = 10,
  title = "Volcano Plot: {comparison}"                                                                                    
)                                                                                                                    

# Mostrar un volcano plot                                                                                            
volcano_plots[["B-A"]]   


# ===== 5. PCA =====
# Interactive PCA plots


pca_input <- arrow::read_parquet("./results/PCA_Input.parquet")

# --- Generar PCA plots con elipse de confianza ---
hc_pcas <- pca_highchart_list(
  pca_input     = pca_input,
  modes         = c("all", "any", "B-A", "C-A", "D-A", "C-B", "D-B", "D-C"),
  group_order   = c("A", "B", "C", "D"),
  ellipse_type  = "confidence",
  ellipse_level = 0.95,
  ellipse_fill_opacity = 0.15,
  point_size = 5,
  filter_samples_to_comparison = TRUE,
  show_labels = TRUE,
  label_size = 12
)
hc_pcas[["all"]]
hc_pcas[["any"]]
hc_pcas[["B-A"]]
hc_pcas[["C-A"]]
hc_pcas[["D-A"]]
hc_pcas[["C-B"]]
hc_pcas[["D-B"]]
hc_pcas[["D-C"]]


# ===== 6. HEATMAP =====
# Static heatmaps with clustering


hm <- proteomics_heatmap(
  data = pca_input,
  mode = "any"
)

# ===== 7. PATTERN PROFILER =====
# Pattern profiling analysis and visualization

# =============================================================================
# STEP 1: CLUSTERING ANALYSIS (Mfuzz soft clustering)
# =============================================================================


conditions = c("A", "B", "C", "D")

# Requiere: se_proc (SummarizedExperiment) y DEPs_results (dataframe DE)
pp_result <- pattern_profiler_analysis(
  se_proc         = result$se_proc,
  DEPs_results    = result$DEPs_results,
  assay_name      = "ImpSeqRob_Min",       # assay del SE a usar
  filter_mode     = "any",                 # "any", "all", o "specific"
  alpha           = 0.05,                  # umbral de significancia
  condition_order = conditions,            # orden del eje X
  aggregate       = "median",              # aggregation per condition
  c_range         = 2:8,                   # rango de clusters a evaluar
  auto_select_c   = TRUE,                  # automatic choice of c
  selection_method = "xb",                 # "xb", "consensus", o "elbow"
  min_membership  = 0.25,                  # membership threshold
  output_file     = "data-raw/Pattern_Profiler_Input.parquet"
)

# Inspeccionar resultados
pp_result$optimal_c          # optimal number of clusters
pp_result$selection_metrics  # evaluation metrics (XB, FPC, AMM, Dmin)
pp_result$cluster_counts     # proteins per cluster


# =============================================================================
# STEP 2: INTERACTIVE VISUALIZATION (Highcharts)
# =============================================================================


# Leer datos desde el parquet generado
data <- read_pattern_profiler_data("data-raw/Pattern_Profiler_Input.parquet")

# Statistical summary
summary <- summarize_pattern_profiler(data)
print(summary$cluster_summary)

# --- Option A: plot for ONE specific cluster ---
hc_c1 <- cluster_profile_highchart(data,
                                   cluster = 3,
                                   conditions = conditions,
                                   line_opacity = 0.6,
                                   line_width = 1.2,
                                   centroid_width = 4)
hc_c1


# --- Option B: list of plots for ALL clusters ---
hc_profiles <- cluster_profile_highchart_list(
  data       = data,
  conditions = conditions,
  palette    = "ggsci::nrc_npg"  # paleta personalizada (opcional)
)
hc_profiles[["Cluster_1"]]
hc_profiles[["Cluster_2"]]

# --- Option C: comparative plot of CENTROIDS ---
hc_centroids <- cluster_centroids_highchart(
  data       = data,
  conditions = conditions
)
hc_centroids


hc_patterns <- cluster_profile_highchart(
  data = pp_result,
  conditions = c("A", "B", "C", "D"),
  cluster = 1,
  line_opacity = 0.6,
  line_width = 1.2,
  centroid_width = 4
)

hc_patterns


# Display first pattern plot
hc_patterns[[1]]
