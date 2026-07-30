
# =========================================================
# Código de ejemplo: Preprocessing
# =========================================================

# Cargar las funciones para el Heatmap-Plot
library(NADIA)

# Ejemplo 1: Uso básico
result <- preprocess_spectronaut(
  file_path = "data-raw/Curso_Q24_DIA_Spectronaut_v20_Report.tsv",
  condition_order = c("A", "B", "C", "D")
)

# Verificar estructura
print(result)
names(result)
head(result$metadata)
head(result$protein_id)
head(result$protein_quant)

# Ejemplo 2: Con exportación y agregadores personalizados
result <- preprocess_spectronaut(
  file_path = "data-raw/Curso_Q24_DIA_Spectronaut_v20_Report.tsv",
  condition_order = c("A", "B", "C", "D"),
  export_dir = "./results",
  timestamp_suffix = TRUE,
  verbose = TRUE
)

# Ejemplo 3: Sin timestamp en nombres de archivo
result <- preprocess_spectronaut(
  file_path = "data-raw/Spectronaut_Report.tsv",
  condition_order = c("Control", "Treatment"),
  export_dir = "./output",
  timestamp_suffix = FALSE
)
