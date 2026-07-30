# =============================================================================
# EJEMPLO COMPLETO: Pattern Profiler con filter_mode = "any"
# =============================================================================

# -----------------------------------------------------------------------------
# PARTE 1: ANÁLISIS (Pattern_Profiler_Analysis.R)
# -----------------------------------------------------------------------------

# Cargar el script de análisis
library(NADIA)

# Cargar datos (se asume que ya existen en el entorno)
se_proc <- readRDS("data-raw/se_proc.rds")           # SummarizedExperiment
DEPs_results <- readr::read_tsv("data-raw/DEPs_results.tsv") # Resultados DEPs

# Definir el orden de condiciones (debe coincidir con tu diseño experimental)
condition_order <- c("A", "B", "C", "D")

# Ejecutar análisis con filter_mode = "any"
# (incluye proteínas significativas en AL MENOS UNA comparación)
result <- pattern_profiler_analysis(
  se_proc = se_proc,
  DEPs_results = DEPs_results,
  assay_name = "LoessCyc",
  filter_mode = "any",
  alpha = 0.05,
  condition_order = condition_order,
  aggregate = "median",
  c_range = 2:8,
  auto_select_c = TRUE,
  selection_method = "xb",
  min_membership = 0.7,
  output_file = "data-raw/Pattern_Profiler_Input.parquet",
  verbose = TRUE
)

# Ver métricas de selección de clusters
print(result$selection_metrics)

# -----------------------------------------------------------------------------
# PARTE 2: VISUALIZACIÓN (Pattern_Profiler_Highcharts.R)
# -----------------------------------------------------------------------------

# Cargar el script de visualización

# Leer los datos del parquet generado
data <- read_pattern_profiler_data("data-raw/Pattern_Profiler_Input.parquet")

# Ver resumen de los datos
summary <- summarize_pattern_profiler(data)
print(summary$cluster_summary)
print(paste("Features en múltiples clusters:", summary$n_multi_cluster_features))

# Definir condiciones (mismo orden que en análisis)
conditions <- c("A", "B", "C", "D")

# --- Gráfico de perfil de un cluster específico ---
hc_cluster1 <- cluster_profile_highchart(
  data = data,
  cluster = 1,
  conditions = conditions,
  show_centroid = TRUE,
  centroid_summary = "median"
)
hc_cluster1

# --- Lista de gráficos para todos los clusters ---
hc_profiles <- cluster_profile_highchart_list(
  data = data,
  conditions = conditions,
  min_membership = NULL,  # Usar todos (ya filtrados en análisis)
  show_centroid = TRUE,
  centroid_summary = "mean",
  palette = NULL  # Paleta por defecto
)

# Mostrar cluster 1
hc_profiles[["Cluster_1"]]

# Mostrar cluster 2
hc_profiles[["Cluster_2"]]

# --- Gráfico comparativo de centroides ---
hc_centroids <- cluster_centroids_highchart(
  data = data,
  conditions = conditions,
  centroid_summary = "median",
  line_width = 2.5,
  show_markers = TRUE
)
hc_centroids

# -----------------------------------------------------------------------------
# OPCIONES ADICIONALES
# -----------------------------------------------------------------------------

# Con paleta de colores personalizada (ggsci)
hc_profiles_color <- cluster_profile_highchart_list(
  data = data,
  conditions = conditions,
  palette = "ggsci::nrc_npg"
)

hc_profiles_color

# Con filtro adicional de membership más estricto
hc_profiles_strict <- cluster_profile_highchart_list(
  data = data,
  conditions = conditions,
  min_membership = 0.5  # Solo proteínas con membership >= 0.5
)

hc_profiles_strict

# Seleccionar clusters específicos
hc_selected <- cluster_profile_highchart_list(
  data = data,
  conditions = conditions,
  clusters = c(1, 2, 3)  # Solo clusters 1, 3 y 5
)

hc_selected
