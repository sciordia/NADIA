
# =========================================================
# Código de ejemplo: Heatmap
# =========================================================

# --- Cargar datos ---
hm_input <- arrow::read_parquet("./results/processing/PCA_Input.parquet")

# Cargar las funciones para el Heatmap-Plot
source("./R/Heatmap_tidyHeatmap.R")


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

  
  ## Test 2 File Annotation Group: ANY ##
  
  # 1. Exportar FeatureIDs del heatmap:
  hm <- proteomics_heatmap(
    data = hm_input,
    mode = "any",
    export_path = "./export/mis_proteinas_any.tsv"
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
    mode = "any",
    
    # Cargar anotaciones
    row_annotation = "./export/clustering_results_any.tsv",
    
    # Especificar qué columnas mostrar (opcional, default: todas)
    row_annotation_cols = c("Cluster", "Membership"),
    
    # Paletas personalizadas por anotación (opcional)
    row_annotation_palette = list(
      Cluster = "brewer:Set1",
      Membership = c("blue", "red", "green")
    ),
    
    # Ordenar filas por una anotación
    row_order_by = "Cluster",
    
    # Separar filas por grupos (como column_split pero para filas)
    split_rows_by = "Cluster"
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
