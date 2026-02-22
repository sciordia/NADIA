
# =========================================================
# Código de ejemplo: Volcano-Plot
# =========================================================

# DEPs_results <- read_tsv("DEPs_results.tsv")
DEPs_results <- arrow::read_parquet("./results/processing/VolcanoPlot_Input.parquet")

# Cargar las funciones para el Volcano-Plot
source("./R/Volcano_Plot_Highcharts_Final.R")

# --- Ejemplo básico ---
hc_volcanos <- volcano_highchart_list(
  de_res      = DEPs_results,
  ain         = "Cycloess",
  comparisons = c("D-A"),
  alpha       = 0.05,
  p_col       = "adj.P.Val",
  point_size  = 3
)

hc_volcanos[["D-A"]]

# --- Combinando todas las opciones ---
hc_volcanos <- volcano_highchart_list(
  de_res          = DEPs_results,
  ain             = "Cycloess",
  comparisons     = c("B-A", "C-A", "D-A"),
  lfc_thr         = 0,
  alpha           = 0.05,
  point_size      = 3,
  show_top_genes  = 5,
  highlight_genes = c("EGFR", "plaP", "SEC6"),
  title           = "Análisis Diferencial",
  palette         = "ggsci::nrc_npg"
)

hc_volcanos[["B-A"]]



# =========================================================
# Código de ejemplo: Box-Plot
# =========================================================

se_proc <- arrow::read_parquet("./results/processing/BoxPlot_Input.parquet")

# Cargar las funciones para el Box-Plot
source("./R/Boxplot_Highcharts_Final.R")

hc_boxplots <- boxplot_highchart_list(
  data        = se_proc,
  assays      = c("log2", "Cycloess"),
  color_by    = "Condition",
  palette  = "ggsci::category10_d3",
  group_order = c("A", "B", "C", "D"),
  box_width = 20
)

hc_boxplots[["log2"]]
hc_boxplots[["Cycloess"]]


# =========================================================
# Código de ejemplo: PCA-Plot
# =========================================================

pca_input <- arrow::read_parquet("./results/processing/PCA_Input.parquet")

# Cargar las funciones para el PCA-Plot
source("./R/PCA_Highcharts_Final.R")

# --- Ejemplo básico: PCA con convex hull (default) ---
sc_all <- build_pca_scores(pca_input, mode = "all")
p_all <- pca_highchart(sc_all, color_by = "Condition", group_order = c("A","B","C","D"))
p_all

# --- PCA con elipse de confianza 95% ---
sc_all <- build_pca_scores(pca_input, mode = "all")
p_all <- pca_highchart(
  sc_all,
  color_by = "Condition",
  group_order = c("A","B","C","D"),
  ellipse_type = "confidence",
  ellipse_level = 0.95
)
p_all

# --- PCA con elipse de confianza 68% (1 desviación estándar) ---
p_68 <- pca_highchart(
  sc_all,
  color_by = "Condition",
  ellipse_type = "confidence",
  ellipse_level = 0.68,
  ellipse_fill_opacity = 0.2
)
p_68

# --- PCA sin elipses (solo puntos) ---
p_noellipse <- pca_highchart(
  sc_all,
  color_by = "Condition",
  addEllipses = FALSE,
  point_size = 6
)
p_noellipse

# --- Generar múltiples PCA plots con convex hull ---
hc_pcas <- pca_highchart_list(
  pca_input   = pca_input,
  modes       = c("all", "any", "B-A", "C-A", "D-A"),
  alpha       = 0.05,
  group_order = c("A", "B", "C", "D"),
  ellipse_type = "convex"
)
hc_pcas[["all"]]
hc_pcas[["any"]]
hc_pcas[["B-A"]]

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

# --- Con paleta personalizada I---
my_palette <- c(A = "#457B9D", B = "#E63946", C = "#2A9D8F", D = "#E9C46A")
hc_pcas <- pca_highchart_list(
  pca_input   = pca_input,
  modes       = c("all", "any"),
  group_order = c("A", "B", "C", "D"),
  palette     = my_palette
)

hc_pcas[["all"]]
hc_pcas[["any"]]

# --- Con paleta personalizada II: paletteer ---
hc_pcas <- pca_highchart_list(
  pca_input   = pca_input,
  modes       = c("all", "any"),
  group_order = c("A", "B", "C", "D"),
  palette     = "ggsci::category10_d3"
)

hc_pcas[["all"]]
hc_pcas[["any"]]

# --- Con paleta personalizada III: R2ColorBrewer ---
hc_pcas <- pca_highchart_list(
  pca_input   = pca_input,
  modes       = c("all", "any"),
  group_order = c("A", "B", "C", "D"),
  palette     = "brewer:Set1"
)

hc_pcas[["all"]]
hc_pcas[["any"]]

# --- Comparar diferentes niveles de confianza ---
sc_all <- build_pca_scores(pca_input, mode = "all")

# 68% (1 SD)
p_68 <- pca_highchart(sc_all, ellipse_type = "confidence", ellipse_level = 0.68,
                      title = "PCA - 68% CI")
# 95% (2 SD aprox)
p_95 <- pca_highchart(sc_all, ellipse_type = "confidence", ellipse_level = 0.95,
                      title = "PCA - 95% CI")
# 99%
p_99 <- pca_highchart(sc_all, ellipse_type = "confidence", ellipse_level = 0.99,
                      title = "PCA - 99% CI")

p_68
p_95
p_99


# =========================================================
# Código de ejemplo: Heatmap
# =========================================================

# --- Cargar datos ---
hm_input <- arrow::read_parquet("./results/processing/PCA_Input.parquet")

# Cargar las funciones para el Heatmap-Plot
source("./R/Heatmap_Highcharts.R")


# --- Ejemplo básico: Heatmap de todas las proteínas ---
hm_all <- proteomics_heatmap(
  data = hm_input,
  mode = "all",
  scale_data = "row",
  sample_order = "condition",
  condition_order = c("A", "B", "C", "D")
)
hm_all

# --- Heatmap con clustering de columnas ---
hm_cluster <- proteomics_heatmap(
  data = hm_input,
  mode = "all",
  scale_data = "row",
  sample_order = "clustering",
  cluster_columns = TRUE,
  cluster_rows = FALSE
)
hm_cluster

# --- Heatmap de proteínas significativas en cualquier comparación ---
hm_any <- proteomics_heatmap(
  data = hm_input,
  mode = "any",
  scale_data = "row",
  sample_order = "condition",
  condition_order = c("A", "B", "C", "D")
)
hm_any

# --- Heatmap de proteínas significativas en cualquier comparación (por clustering) ---
hm_any <- proteomics_heatmap(
  data = hm_input,
  mode = "any",
  scale_data = "none",
  sample_order = "clustering",
  condition_order = c("A", "B", "C", "D")
)
hm_any

# --- Heatmap de proteínas significativas en una comparación específica ---
hm_target <- proteomics_heatmap(
  data = hm_input,
  mode = "target",
  comparison = "B-A",
  alpha = 0.05,
  scale_data = "row"
)
hm_target

# --- Con paleta personalizada de paletteer ---
hm_viridis <- proteomics_heatmap(
  data = hm_input,
  mode = "all",
  scale_data = "row",
  palette_value = "viridis::viridis",
  reverse_palette = TRUE
)
hm_viridis

# --- Con paleta de RColorBrewer ---
hm_brewer <- proteomics_heatmap(
  data = hm_input,
  mode = "all",
  scale_data = "row",
  palette_value = "brewer:RdYlBu",
  palette_annotation = "brewer:Set1",
  show_annotation = TRUE
)
hm_brewer

# --- Generar múltiples heatmaps ---
hm_list <- proteomics_heatmap_list(
  data = hm_input,
  modes = c("all", "any", "B-A", "C-A", "D-A"),
  scale_data = "row",
  sample_order = "condition",
  split_by_condition = TRUE,
  condition_order = c("A", "B", "C", "D"),
)
hm_list[["all"]]
hm_list[["any"]]
hm_list[["B-A"]]

# --- Orden personalizado de muestras ---
custom_order <- c("A_1", "A_2", "B_1", "B_2", "C_1", "C_2", "D_1", "D_2")
hm_custom <- proteomics_heatmap(
  data = hm_input,
  mode = "all",
  scale_data = "row",
  sample_order = custom_order,
  palette_annotation = c(A = "#E63946", B = "#457B9D", C = "#2A9D8F", D = "#E9C46A"),
  cluster_columns = FALSE,
  split_by_condition = TRUE
)
hm_custom

# Paleta personalizada con 3 colores
hm <- proteomics_heatmap(
  data = hm_input,
  mode = "target",
  comparison = "D-A",
  cluster_rows = TRUE,
  split_by_condition = TRUE,
  show_adjp_annotation = TRUE,
  palette_adjp = c("#006837", "#addd8e", "#FFFFFF")
)

hm

# Título simple
hm <- proteomics_heatmap(
  data = hm_input,
  mode = "any",
  heatmap_title = "Differential Expression Analysis",
  heatmap_title_size = 16,
  heatmap_title_face = "bold",
  cluster_columns = FALSE
)

hm

# En proteomics_heatmap_list con placeholder {mode}
hm_list <- proteomics_heatmap_list(
  data = hm_input,
  modes = c("all", "any", "B-A"),
  heatmap_title = "Proteomics: {mode}",
  cluster_rows = TRUE
)

hm_list[["all"]]
hm_list[["any"]]
hm_list[["B-A"]]

# Con número (interpretado como milímetros)
hm <- proteomics_heatmap(
  data = hm_input,
  mode = "any",
  cluster_columns = TRUE,
  column_dend_height = 30,  # 30mm
  row_dend_width = 20,       # 20mm (si cluster_rows = TRUE)
  cluster_rows = TRUE
)

hm

# Con unidades grid explícitas
hm <- proteomics_heatmap(
  data = hm_input,
  mode = "any",
  column_dend_height = grid::unit(2, "cm"),
  row_dend_width = grid::unit(1.5, "cm")
)

hm

# Ocultar leyenda del heatmap
hm <- proteomics_heatmap(
  data = hm_input,
  mode = "any",
  show_heatmap_legend = FALSE,
  show_annotation = TRUE
)

hm

# Mostrar leyenda (default)
hm <- proteomics_heatmap(
  data = hm_input,
  mode = "any",
  show_heatmap_legend = TRUE
)

hm

# Ocultar todas las leyendas
hm <- proteomics_heatmap(
  data = hm_input,
  mode = "any",
  show_heatmap_legend = FALSE,
  show_annotation_legend = FALSE
)

hm

# Solo leyenda del heatmap, sin leyenda de anotaciones
hm <- proteomics_heatmap(
  data = hm_input,
  mode = "any",
  show_heatmap_legend = TRUE,
  show_annotation_legend = TRUE
)

hm

# Borde negro en cada celda del heatmap (si hay muchas proteínas se ve todo negro)
hm <- proteomics_heatmap(
  data = hm_input,
  mode = "any",
  border_color = TRUE
)

hm

# Borde gris claro en cada celda del heatmap (si hay muchas proteínas se ve todo gris)
hm <- proteomics_heatmap(
  data = hm_input,
  mode = "any",
  border_color = "grey80"
)

hm

# Sin borde (default)
hm <- proteomics_heatmap(
  data = hm_input,
  mode = "any",
  border_color = NULL
)

hm

# Ajustar tamaño de títulos
hm <- proteomics_heatmap(
  data = hm_input,
  mode = "any",
  row_title = "Proteins",
  column_title = "Samples",
  row_title_size = 12,
  column_title_size = 14
)

hm

# Ocultar ambos títulos
hm <- proteomics_heatmap(
  data = hm_input,
  mode = "any",
  show_row_title = FALSE,
  show_column_title = FALSE
)

hm

# Solo mostrar título de columnas
hm <- proteomics_heatmap(
  data = hm_input,
  mode = "any",
  show_row_title = FALSE,
  show_column_title = TRUE
)

hm

# Mostrar solo proteínas específicas en modo target
my_proteins <- c("P07278", "Q9Y680", "P08814;P20962")

hm <- proteomics_heatmap(
  data = hm_input,
  mode = "target",
  comparison = "B-A",
  feature_ids = my_proteins,  # Solo estas proteínas
  scale_data = "row",
  split_by_condition = TRUE
)

hm

# También funciona con otros modos (ignora el filtrado por mode/alpha)
hm <- proteomics_heatmap(
  data = hm_input,
  mode = "all",
  feature_ids = my_proteins,
  sample_order = "condition"
  
)

hm

# Comportamiento:
# - Si feature_ids es NULL (default): usa el filtrado normal por mode/alpha
# - Si feature_ids contiene IDs: muestra solo esas proteínas, ignorando mode/alpha
# - Si algunos IDs no existen en los datos, muestra un warning y los ignora

# Exportar datos del heatmap a TSV
hm <- proteomics_heatmap(
  data = hm_input,
  mode = "target",
  comparison = "B-A",
  feature_ids = c("P07278", "Q9Y680", "P08814;P20962"),
  export_path = "./export/heatmap_data.tsv"
)
# Mensaje: "Datos exportados a: heatmap_data.tsv"

# Con proteomics_heatmap_list (genera archivos por modo)
hm_list <- proteomics_heatmap_list(
  data = hm_input,
  modes = c("all", "any", "B-A"),
  export_path = "./export/export.tsv"
)
# Genera: export_all.tsv, export_any.tsv, export_B-A.tsv

# Exportar solo ciertos modos
hm_list <- proteomics_heatmap_list(
  data = hm_input,
  modes = c("all", "any", "B-A", "C-A"),
  export_path = "./export/export.tsv",
  export_modes = c("B-A", "C-A")  # Solo exporta estos dos
)
# Genera: export_B-A.tsv, export_C-A.tsv

# Exportar solo uno
hm_list <- proteomics_heatmap_list(
  data = hm_input,
  modes = c("all", "any", "B-A"),
  export_path = "./export/export.tsv",
  export_modes = "B-A"  # Solo exporta B-A
)
# Genera: export_B-A.tsv

# Exportar todos (comportamiento por defecto)
hm_list <- proteomics_heatmap_list(
  data = hm_input,
  modes = c("all", "any"),
  export_path = "./export/export.tsv",
  export_modes = NULL  # Exporta all y any
)

## Test File Annotation Group ##

# 1. Exportar FeatureIDs del heatmap:
  hm <- proteomics_heatmap(
    data = hm_input,
    mode = "target",
    comparison = "B-A",
    export_path = "./export/mis_proteinas.tsv"
  )

# 2. Editar el TSV añadiendo columnas categóricas:
#  FeatureID    Pathway       Function      Cluster
# P12345       Glycolysis    Enzyme        A
# Q67890       TCA_Cycle     Transporter   B
# P11111       Glycolysis    Receptor      A
# ...

# 3. Usar las anotaciones en el heatmap:
  hm <- proteomics_heatmap(
    data = hm_input,
    mode = "target",
    comparison = "B-A",
    
    # Cargar anotaciones
    row_annotation = "./export/mis_proteinas_mod.tsv",
    
    # Especificar qué columnas mostrar (opcional, default: todas)
    row_annotation_cols = c("Specie", "Process", "Function"),
    
    # Paletas personalizadas por anotación (opcional)
    row_annotation_palette = list(
      Specie = "brewer:Set1",
      Process = "brewer:Set1",
      Function = c("blue", "red", "green")
    ),
    
    # Ordenar filas por una anotación
    row_order_by = "Specie",
    
    # Separar filas por grupos (como column_split pero para filas)
    split_rows_by = "Specie"
  )

hm

#Ejemplo de ordenamiento por anotaciones múltiples + adjP:
hm <- proteomics_heatmap(
  data = hm_input,
  mode = "target",
  comparison = "B-A",
  row_annotation = "./export/mis_proteinas_mod.tsv",
  row_annotation_cols = c("Specie", "Process", "Function"),
  row_annotation_size = 0.8,
  row_annotation_name_size = 7,
  row_order_by = c("Specie", "Process", "Function", "adjP"),
  split_rows_by = "Specie"
)

hm

# Exportar a PNG (3000x2400px, 300 dpi)
proteomics_heatmap(
  data = hm_input,
  mode = "target",
  comparison = "B-A",
  export_file = "./export/heatmap.png",
  plot_width = 10,    # pulgadas
  plot_height = 8,    # pulgadas
  export_dpi = 300    # resolución
)

# Exportar a SVG (vectorial)
proteomics_heatmap(
  data = hm_input,
  mode = "target",
  comparison = "B-A",
  export_file = "./export/heatmap.svg",
  plot_width = 10,    # pulgadas
  plot_height = 8,    # pulgadas
  export_dpi = 300    # resolución
)

# Exportar a PDF (vectorial)
proteomics_heatmap(
  data = hm_input,
  mode = "target",
  comparison = "B-A",
  export_file = "./export/heatmap.pdf",
  plot_width = 10,    # pulgadas
  plot_height = 8,    # pulgadas
  export_dpi = 300    # resolución
)

# proteomics_heatmap_list():
# Exporta múltiples heatmaps con nombres automáticos
# (heatmap_all.png, heatmap_any.png, heatmap_B-A.png, etc.)
  
  proteomics_heatmap_list(
    data = hm_input,
    modes = c("all", "any", "B-A"),
    export_file = "export/heatmap.png",
    plot_width = 10,    # pulgadas
    plot_height = 8,    # pulgadas
    export_dpi = 300    # resolución
  )

# Parámetros:
# - export_file: Ruta del archivo (formato por extensión: .png, .svg, .pdf)
# - plot_width: Ancho en pulgadas (default: 10)
# - plot_height: Alto en pulgadas (default: 8)
# - export_dpi: Resolución para PNG (default: 300)


# =========================================================
# Código de ejemplo: Pattern Profileer
# =========================================================

  # =============================================================================
  # TEST: Pattern Profiler - Todas las funcionalidades
  # =============================================================================
  
  # Cargar el script
  source("R/Pattern_Profiler_Highcharts.R")
  
  # -----------------------------------------------------------------------------
  # 1) Cargar datos
  # -----------------------------------------------------------------------------
  data <- arrow::read_parquet("data/PCA_Input.parquet")
  
  # Ver estructura de los datos
  cat("Dimensiones:", nrow(data), "x", ncol(data), "\n")
  cat("Columnas:", paste(names(data), collapse = ", "), "\n")
  cat("Condiciones únicas:", paste(unique(data$Condition), collapse = ", "), "\n")
  
  # -----------------------------------------------------------------------------
  # 2) Ejecutar clustering con selección automática de clusters
  # -----------------------------------------------------------------------------
  result <- pattern_profiler_analysis(
    data = data,
    filter_mode = "any",
    condition_order = c("A", "B", "C", "D"),
    c_range = 2:8,
    auto_select_c = TRUE,
    selection_method = "xb",
    min_membership = 0.4,
    verbose = TRUE
  )
  
  # -----------------------------------------------------------------------------
  # 3) Ver resultados del clustering
  # -----------------------------------------------------------------------------
  cat("\n=== RESULTADOS ===\n")
  cat("Número óptimo de clusters:", result$optimal_c, "\n")
  cat("Parámetro m:", round(result$m, 3), "\n")
  cat("Total proteínas:", result$n_features, "\n")
  
  # Métricas de selección
  cat("\nMétricas de evaluación:\n")
  print(result$selection_metrics)
  
  # Tabla de membership (primeras filas)
  cat("\nTabla de membership (primeras 10 filas):\n")
  print(head(result$membership_table, 10))
  
  # Resumen por cluster
  cat("\nProteínas por cluster:\n")
  print(table(result$membership_table$Cluster))
  
  # -----------------------------------------------------------------------------
  # 4) Gráfico de un cluster específico
  # -----------------------------------------------------------------------------
  # Cluster 1 con configuración por defecto
  hc_c1 <- cluster_profile_highchart(
    result,
    cluster = 1
  )
  hc_c1
  
  # Cluster 2 con más opacidad
  hc_c2 <- cluster_profile_highchart(
    result,
    cluster = 2,
    line_opacity = 0.6,
    line_width = 1.2,
    centroid_width = 4
  )
  hc_c2
  
  # -----------------------------------------------------------------------------
  # 5) Lista de gráficos para todos los clusters
  # -----------------------------------------------------------------------------
  # Con configuración por defecto
  hc_profiles <- cluster_profile_highchart_list(result)
  
  # Ver nombres de los gráficos generados
  cat("\nGráficos generados:", paste(names(hc_profiles), collapse = ", "), "\n")
  
  # Visualizar cada uno
  hc_profiles[["Cluster_1"]]
  hc_profiles[["Cluster_2"]]
  hc_profiles[["Cluster_3"]]
  
  # Con paleta personalizada
  hc_profiles_custom <- cluster_profile_highchart_list(
    result,
    palette = "ggsci::nrc_npg",
    line_opacity = 0.5,
    line_width = 1
  )
  hc_profiles_custom[["Cluster_1"]]
  
  # -----------------------------------------------------------------------------
  # 6) Gráfico comparativo de centroides
  # -----------------------------------------------------------------------------
  # Todos los clusters
  hc_centroids <- cluster_centroids_highchart(result)
  hc_centroids
  
  # Con mediana en lugar de media
  hc_centroids_med <- cluster_centroids_highchart(
    result,
    centroid_summary = "median",
    title = "Centroides de Clusters (mediana)"
  )
  hc_centroids_med
  
  # Solo algunos clusters
  hc_centroids_subset <- cluster_centroids_highchart(
    result,
    clusters = c(1, 2),
    line_width = 3
  )
  hc_centroids_subset
  
  # -----------------------------------------------------------------------------
  # 7) Exportar resultados a TSV
  # -----------------------------------------------------------------------------
  # Con z-scores
  export_clustering_results(
    result,
    file_path = "./export/clustering_results_full.tsv",
    include_zscores = TRUE
  )
  
  # Sin z-scores (solo membership)
  export_clustering_results(
    result,
    file_path = "./export/clustering_results_simple.tsv",
    include_zscores = FALSE
  )
  
  # -----------------------------------------------------------------------------
  # 8) Pruebas adicionales
  # -----------------------------------------------------------------------------
  
  # --- Clustering con número fijo de clusters ---
  result_fixed <- pattern_profiler(
    data = data,
    filter_mode = "any",
    condition_order = c("A", "B", "C", "D"),
    auto_select_c = FALSE,
    c = 4,
    min_membership = 0.4
  )
  cat("\nClustering con c=4 fijo:", result_fixed$optimal_c, "clusters\n")
  
  # --- Filtrar por comparación específica ---
  result_specific <- pattern_profiler(
    data = data,
    filter_mode = "specific",
    comparison = "B-A",
    alpha = 0.01,
    condition_order = c("A", "B", "C", "D"),
    auto_select_c = TRUE
  )
  
  cat("\nClustering solo sig. en B-A:", result_specific$n_features, "proteínas\n")
  
  # --- Diferentes métodos de selección ---
  result_consensus <- pattern_profiler(
    data = data,
    filter_mode = "any",
    condition_order = c("A", "B", "C", "D"),
    selection_method = "consensus",
    verbose = FALSE
  )
  cat("\nMétodo consensus: c =", result_consensus$optimal_c, "\n")
  
  result_elbow <- pattern_profiler(
    data = data,
    filter_mode = "any",
    condition_order = c("A", "B", "C", "D"),
    selection_method = "elbow",
    verbose = FALSE
  )
  cat("Método elbow: c =", result_elbow$optimal_c, "\n")
  
  # -----------------------------------------------------------------------------
  # 9) Visualización en grid (requiere htmltools)
  # -----------------------------------------------------------------------------
  if (requireNamespace("htmltools", quietly = TRUE)) {
    library(htmltools)
    
    # Crear grid de todos los clusters
    grid <- browsable(
      tagList(
        tags$div(
          style = "display: grid; grid-template-columns: repeat(2, 1fr); gap: 10px;",
          lapply(hc_profiles, function(hc) {
            tags$div(hc)
          })
        )
      )
    )
    
    # Guardar como HTML
    htmltools::save_html(grid, "./export/cluster_profiles_grid.html")
    cat("\nGrid guardado en: export/cluster_profiles_grid.html\n")
  }
  
  cat("\n=== TESTS COMPLETADOS ===\n")
  
  # Este código prueba:
  #   1. Carga de datos desde TSV
  # 2. Clustering automático con método Xie-Beni
  # 3. Inspección de resultados (métricas, tabla de membership)
  # 4. Gráficos individuales con diferentes configuraciones
  # 5. Lista de gráficos para todos los clusters
  # 6. Gráfico de centroides comparativo
  # 7. Exportación a TSV (con/sin z-scores)
  # 8. Casos adicionales: clusters fijos, filtros específicos, diferentes métodos
  # 9. Grid HTML opcional para ver todos los clusters juntos