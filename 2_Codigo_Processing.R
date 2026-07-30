
# =========================================================
# Código de ejemplo: Processing
# =========================================================

# Cargar las funciones para el Heatmap-Plot
library(NADIA)

# =============================================================================                                      
# EJEMPLO COMPLETO: Pipeline de Procesamiento Proteómico                                                             
# =============================================================================                                      
#                                                                                                                    
# Este ejemplo muestra el flujo completo desde un reporte de Spectronaut                                             
# hasta los resultados de expresión diferencial.                                                                     
#                                                                                                                    
# Prerequisitos:                                                                                                     
#   - Archivo TSV de Spectronaut                                                                                     
#   - Paquetes: SummarizedExperiment, limma, MSnbase, Biobase, rrcovNA                                               
# =============================================================================                                      

# --- 1. Cargar scripts ---                                                                                          

# =============================================================================                                      
# PASO 1: PREPROCESAMIENTO DE DATOS DE SPECTRONAUT                                                                   
# =============================================================================                                      

# Definir la ruta al archivo de Spectronaut                                                                          
file_path <- "data-raw/Curso_Q24_DIA_Spectronaut_v20_Report.tsv"                                                                         

# Definir el orden de las condiciones experimentales                                                                 
# IMPORTANTE: El orden determina las comparaciones (primera condición suele ser control)                             
condition_order <- c("A", "B", "C", "D")                                         

# Ejecutar preprocesamiento                                                                                          
preprocessing <- preprocess_spectronaut(                                                                             
  file_path = file_path,                                                                                             
  condition_order = condition_order,                                                                                 
  export_dir = "./results/preprocessing",  # Exportar archivos intermedios                                           
  agg_coverage_run = "max",                 # Agregación de coverage por run                                         
  agg_mw = "max",                           # Agregación de peso molecular                                           
  timestamp_suffix = TRUE,                  # Añadir timestamp a archivos                                            
  
  verbose = TRUE                            # Mostrar progreso                                                         
)                                                                                                                    

# --- Explorar el objeto preprocessing ---                                                                           
print(preprocessing)                                                                                                 

# Estructura del objeto spectronaut_data:                                                                            
preprocessing$metadata       # Información de runs/muestras                                                        
preprocessing$protein_id     # Métricas de identificación                                                          
preprocessing$protein_quant  # Cuantificación de proteínas                                                         

# Ver metadata de muestras                                                                                           
head(preprocessing$metadata)                                                                                         
#   R.FileName       R.Condition R.Replicate     Coding                                                              
# 1 sample_01.raw    Control     1           Control_1                                                               
# 2 sample_02.raw    Control     2           Control_2                                                               
# ...                                                                                                                

# Ver estructura de cuantificación                                                                                   
str(preprocessing$protein_quant)                                                                                     
# Columnas incluyen: PG.ProteinGroups, PG.Genes, PG.Quantity_Control_1, etc.                                         

# =============================================================================                                      
# PASO 2: PROCESAMIENTO (Normalización, Imputación, Análisis DE)                                                     
# =============================================================================                                      

# --- 2A. Uso básico con parámetros por defecto ---                                                                  
processing <- process_proteomics(                                                                                        
  preprocessing = preprocessing,                                                                                     
  export_dir = "./results/processing",                                                                               
  verbose = TRUE                                                                                                     
)                                                                                                                    

# --- 2B. Uso personalizado con todos los parámetros ---                                                             
processing <- process_proteomics(                                                                                        
  preprocessing = preprocessing,                                                                                     
  
  # Directorio de salida                                                                                             
  export_dir = "./results/processing",
  export_format = "parquet", # Formato de salida
  
  # --- Filtrado de proteínas ---                                                                                    
  min_reps_filter = 3,       # Mínimo 3 réplicas con datos por grupo                                                 
  min_groups_filter = 1,     # Al menos 1 grupo debe cumplir el criterio                                             
  
  # --- Normalización ---                                                                                            
  normalization_method = "cyclicloess",  # Opciones: "cyclicloess", "quantile", "scale"                              
  
  # --- Imputación MNAR ---                                                                                          
  prop_na_mnar = 0.51,       # >51% NA en condición = candidato a MNAR                                               
  prop_present_mar = 0.5,    # ≥50% presente en otras condiciones                                                    
  min_present_mar = 1,       # Mínimo 1 valor presente                                                               
  require_n_conditions = 1,  # Requerido en al menos 1 condición                                                     
  mar_method = "impSeqRob",  # Método MAR: "impSeqRob", "knn", "none"                                                
  
  # --- Análisis Diferencial ---                                                                                     
  comparisons = NULL,        # NULL = todas las comparaciones pareadas                                               
  control = NULL,            # Condición control (comparar todo vs Control)                                          
  logFC_threshold = 0,       # Umbral de logFC para significancia                                                    
  alpha = 0.05,              # Umbral de p-valor ajustado                                                            
  
  # --- Exportación ---                                                                                              
  export_normalized = TRUE,  # Exportar matriz normalizada                                                           
  export_imputed = TRUE,     # Exportar matriz imputada                                                              
  
  # --- Mensajes ---                                                                                                 
  verbose = TRUE                                                                                                     
)                                                                                                                    

# =============================================================================                                      
# PASO 3: EXPLORAR RESULTADOS                                                                                        
# =============================================================================                                      

# --- Ver resumen del resultado ---                                                                                  
print(processing)                                                                                                        
# === Resultado de Procesamiento Proteómico ===                                                                      
#                                                                                                                    
# SummarizedExperiment:                                                                                              
#   - Proteínas: 2847                                                                                                
#   - Muestras: 16                                                                                                   
#   - Assays: raw, log2, Cycloess                                                                                    
#   - Condiciones: Control, Treatment_A, Treatment_B, Treatment_C                                                    
#                                                                                                                    
# Resultados Diferenciales:                                                                                          
#   - Total filas: 8541                                                                                              
#   - Comparaciones: Treatment_A-Control, Treatment_B-Control, Treatment_C-Control                                   
#     Treatment_A-Control: Up=234, Down=187                                                                          
#     Treatment_B-Control: Up=156, Down=203                                                                          
#     Treatment_C-Control: Up=89, Down=67                                                                            

# --- Acceder a componentes ---                                                                                      

# 3.1 SummarizedExperiment procesado                                                                                 
se <- result$se_proc                                                                                                 

# Ver assays disponibles                                                                                             
SummarizedExperiment::assayNames(se)                                                                                 
# [1] "raw" "log2" "Cycloess"                                                                                        

# Extraer matriz normalizada+imputada                                                                                
mat_final <- SummarizedExperiment::assay(se, "Cycloess")                                                             
head(mat_final[, 1:4])                                                                                               

# Ver metadata de proteínas                                                                                          
row_data <- as.data.frame(SummarizedExperiment::rowData(se))                                                         
head(row_data)                                                                                                       
#   Protein.IDs Gene.Names UniqPepts IDs                                                                             
# 1 P12345      GENE1      15        P12345                                                                          
# 2 Q67890      GENE2      8         Q67890                                                                          

# Ver metadata de muestras                                                                                           
col_data <- as.data.frame(SummarizedExperiment::colData(se))                                                         
print(col_data)                                                                                                      

# 3.2 Resultados de Expresión Diferencial                                                                            
DEPs <- result$DEPs_results                                                                                          
head(DEPs)                                                                                                           
#   Protein.IDs Gene.Names IDs    logFC   P.Value adj.P.Val   Change          Comparison  Assay                      
# 1 P12345      GENE1      P12345  2.34   1.2e-05 3.4e-04   Up       Treatment_A-Control Cycloess                    
# 2 Q67890      GENE2      Q67890 -1.89   2.1e-04 8.9e-03   Down     Treatment_A-Control Cycloess                    

# 3.3 Filtrar proteínas significativas                                                                               
DEPs_sig <- DEPs[DEPs$Change != "No Change", ]                                                                       
nrow(DEPs_sig)                                                                                                       

# Por comparación específica                                                                                         
DEPs_BvsA <- DEPs[DEPs$Comparison == "B-A", ]                                                        
DEPs_BvsA_up <- DEPs_BvsA[DEPs_BvsA$Change == "Up", ]                                                                
DEPs_BvsA_down <- DEPs_BvsA[DEPs_BvsA$Change == "Down", ]                                                            

cat("B vs A:\n")                                                                                     
cat("  - Up-regulated:", nrow(DEPs_BvsA_up), "\n")                                                                   
cat("  - Down-regulated:", nrow(DEPs_BvsA_down), "\n")                                                               

# 3.4 Parámetros utilizados                                                                                          
print(result$parameters)                                                                                             
# $min_reps_filter                                                                                                   
# [1] 3                                                                                                              
# $normalization_method                                                                                              
# [1] "cyclicloess"                                                                                                  
# $alpha                                                                                                             
# [1] 0.05                                                                                                           
# ...                                                                                                                

# =============================================================================                                      
# PASO 4: VISUALIZACIÓN (usando funciones del proyecto)                                                              
# =============================================================================                                      

# Cargar funciones de visualización                                                                                  

# --- 4.1 Volcano Plots ---                                                                                          
volcano_plots <- volcano_highchart_list(                                                                             
  DEPs,                                                                                                              
  logFC_col = "logFC",                                                                                               
  pval_col = "adj.P.Val",                                                                                            
  gene_col = "Gene.Names",                                                                                           
  comparison_col = "Comparison",                                                                                     
  logFC_threshold = 0,                                                                                             
  pval_threshold = 0.05,                                                                                             
  title_prefix = "Volcano Plot: "                                                                                    
)                                                                                                                    

# Mostrar un volcano plot                                                                                            
volcano_plots[["Treatment_A-Control"]]                                                                               

# --- 4.2 Heatmap de proteínas significativas ---                                                                    
# Obtener IDs de proteínas significativas                                                                            
sig_proteins <- unique(DEPs_sig$Protein.IDs)                                                                         

# Extraer matriz para heatmap                                                                                        
mat_sig <- mat_final[rownames(mat_final) %in% sig_proteins, ]                                                        

# Crear heatmap                                                                                                      
heatmap_obj <- proteomics_heatmap(                                                                                   
  mat_sig,                                                                                                           
  metadata = col_data,                                                                                               
  cluster_rows = TRUE,                                                                                               
  cluster_columns = FALSE,                                                                                           
  show_row_names = FALSE,                                                                                            
  annotation_col = "Condition"                                                                                       
)                                                                                                                    
print(heatmap_obj)                                                                                                   

# =============================================================================                                      
# PASO 5: EXPORTACIÓN ADICIONAL                                                                                      
# =============================================================================                                      

# --- 5.1 Exportar tabla completa de DEPs ---                                                                        
readr::write_tsv(DEPs, "./results/DEPs_all_comparisons.tsv")                                                         

# --- 5.2 Exportar solo proteínas significativas ---                                                                 
readr::write_tsv(DEPs_sig, "./results/DEPs_significant.tsv")                                                         

# --- 5.3 Exportar por comparación ---                                                                               
comparisons <- unique(DEPs$Comparison)                                                                               
for (comp in comparisons) {                                                                                          
  comp_data <- DEPs[DEPs$Comparison == comp, ]                                                                       
  file_name <- paste0("./results/DEPs_", gsub("-", "_vs_", comp), ".tsv")                                            
  readr::write_tsv(comp_data, file_name)                                                                             
}                                                                                                                    

# --- 5.4 Guardar objeto completo para uso posterior ---                                                             
saveRDS(result, "./results/proteomics_result.rds")                                                                   

# Para cargar después:                                                                                               
# result <- readRDS("./results/proteomics_result.rds")                                                               

# =============================================================================                                      
# CASOS DE USO ESPECIALES                                                                                            
# =============================================================================                                      

# --- Caso A: Solo comparaciones vs control ---                                                                      
result_vs_control <- process_proteomics(                                                                             
  preprocessing = preprocessing,                                                                                     
  control = "Control",  # Todas las comparaciones serán vs Control                                                   
  verbose = TRUE                                                                                                     
)                                                                                                                    

# --- Caso B: Comparaciones específicas ---                                                                          
result_specific <- process_proteomics(                                                                               
  preprocessing = preprocessing,                                                                                     
  comparisons = c("Treatment_A-Control", "Treatment_B-Treatment_A"),                                                 
  verbose = TRUE                                                                                                     
)                                                                                                                    

# --- Caso C: Sin imputación (solo normalización) ---                                                                
result_no_impute <- process_proteomics(                                                                              
  preprocessing = preprocessing,                                                                                     
  mar_method = "none",  # No imputar MAR                                                                             
  verbose = TRUE                                                                                                     
)                                                                                                                    

# --- Caso D: Umbral estricto de logFC ---                                                                           
result_strict <- process_proteomics(                                                                                 
  preprocessing = preprocessing,                                                                                     
  logFC_threshold = 1.0,  # Solo |logFC| >= 1                                                                        
  alpha = 0.01,           # p-valor más estricto                                                                     
  verbose = TRUE                                                                                                     
)                                                                                                                    

# =============================================================================                                      
# ARCHIVOS GENERADOS                                                                                                 
# =============================================================================                                      
#                                                                                                                    
# ./results/processing/                                                                                              
#   ├── matrix_log2_cyclicloess.tsv       # Matriz normalizada                                                       
#   ├── matrix_log2_cyclicloess_imputed.tsv  # Matriz imputada                                                       
#   └── DEPs_results.tsv                  # Resultados DE                                                            
#                                                                                                                    
# =============================================================================      