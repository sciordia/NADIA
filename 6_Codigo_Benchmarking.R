
# =============================================================================                                      
# Ejemplo: Benchmarking con datos spike-in 3 especies
# =============================================================================                                      

library(NADIA)

library(dplyr)

# --- 1. Cargar datos DE ---
de_all <- read.delim("data-raw/DE_ALL_with_LoessCyc_bySample_Candidates_FIXED.tsv",
                     stringsAsFactors = FALSE)

# --- 2. Preparar columnas para el módulo ---
de_res <- de_all %>%
  transmute(
    Protein.IDs = ProteinGroups,
    Gene.Names  = Gene.Names,
    logFC       = `AVG.Log2.Ratio`,
    P.Value     = P.Value,
    adj.P.Val   = Qvalue,
    Change      = Change,
    Comparison  = Comparison,
    Species     = gsub("^_", "", ProteinNames)   # "_ECOLI" -> "ECOLI"
  )

# OPCIONAL1: Quitar los espacios si los hay en la columna 'Comparison'
de_res$Comparison <- gsub("\\s+", "", de_res$Comparison)

# OPCIONAL2: Cambiar las "/" por "-"
de_res$Comparison <- gsub("/", "-", de_res$Comparison)

# Verificar
table(de_res$Species)
table(de_res$Comparison)

# --- 3. Definir ground truth (spike-in design) ---
# Ajustar estos valores a tu diseño experimental real
expected <- tribble(
  ~Comparison, ~Species, ~expected_logFC,
  "B-A", "YEAST", -0.58,
  "B-A", "ECOLI", 1,
  "C-A", "YEAST", -1.60,
  "C-A", "ECOLI", 1.58,
  "D-A", "YEAST", -3.3,
  "D-A", "ECOLI", 2,
  "C-B", "YEAST", -1,
  "C-B", "ECOLI", 0.58,
  "D-B", "YEAST", -3.3,
  "D-B", "ECOLI", 1,
  "D-C", "YEAST", -3.3,
  "D-C", "ECOLI", 0.41
)

# HUMAN no aparece en expected_values -> truth = 0 (no debe cambiar)

# --- 4. Ejecutar benchmarking completo ---
benchmarking <- benchmarking_proteomics(
  de_res          = de_res,
  expected_values = expected,
  alpha           = 0.05,
  lfc_thr         = 0,
  p_col           = "adj.P.Val",
  comparisons     = c("B-A", "C-A", "D-A", "C-B", "D-B", "D-C"),
  output_dir      = "results/benchmark",
  verbose         = TRUE
)

# --- 5. Ver resultados ---
# Tabla de métricas clásicas
benchmarking$metrics_table

# Métricas OpDEA (pAUC, nMCC, G-mean)
benchmarking$opdea_metrics

# Dispersión por especie
benchmarking$dispersion_metrics

# Matrices de Confusión
benchmarking$confusion_by_species
benchmarking$confusion_overall

# --- 6. Visualizaciones ggplot2 ---
benchmarking$gg_heatmap
benchmarking$gg_confusion_by_species
benchmarking$gg_confusion_overall
benchmarking$gg_auc_bars
benchmarking$gg_metrics_bars
benchmarking$gg_signif_bars

# --- 7. Volcano interactivo (Highcharter) ---
benchmarking$hc_volcano_list[["B-A"]]
benchmarking$hc_volcano_list[["C-A"]]
benchmarking$hc_volcano_list[["D-A"]]

# --- 8. Uso individual de funciones ---
# Solo métricas clásicas (sin visualizaciones)
metrics <- compute_benchmark_metrics(
  de_res = de_res,
  ev     = expected,
  alpha  = 0.05,
  comparisons = c("B-A", "C-A")
)
metrics

# Solo métricas OpDEA (pAUC, nMCC, G-mean)
opdea <- compute_opdea_metrics(
  de_res = de_res,
  ev     = expected,
  alpha  = 0.05,
  comparisons = c("B-A", "C-A")
)
opdea

# Solo dispersión
disp <- compute_dispersion_metrics(
  de_res = de_res,
  ev     = expected,
  alpha  = 0.05
)
disp

# Heatmap personalizado
benchmark_heatmap_gg(
  metrics,
  metrics_to_show = c("AUC", "Sensitivity", "Specificity", "F1"),
  title = "Spike-in Benchmark"
)

# Curvas ROC y zoom (0-10% FPR)
# Automático con benchmarking_proteomics():                                                                          
benchmarking$gg_roc
benchmarking$gg_roc_zoom

# O standalone:
benchmark_roc_gg(benchmarking$classified_df, p_col = "adj.P.Val")
benchmark_roc_gg(benchmarking$classified_df, p_col = "adj.P.Val", zoom = TRUE)


# Notas importantes:
#   
# 1. expected_logFC: Los valores del ejemplo son aproximados para un diseño spike-in típico. Ajusta los ratios a tu
# diseño experimental real (ratios de mezcla ECOLI/YEAST/HUMAN).
# 2. Species: Se extrae de la columna ProteinNames quitando el guión bajo. HUMAN actúa como "negativo" (no aparece en
#                                                                                                       expected_values, así que su truth = 0).
# 3. comparisons: Solo incluí las vs control (A). Si quieres benchmarkear C / B o D / C, necesitas añadir sus
# expected_logFC correspondientes a la tabla expected.
# 4. Exportaciones: Con output_dir = "data-raw/benchmark" se generan 4 TSV + 5 PNG automáticamente.
# 5. AUC/pAUC: Requiere install.packages("pROC"). Si no está instalado, AUC y pAUC serán NA pero el resto funciona.
# 6. OpDEA metrics: nMCC y G_mean se calculan siempre (no requieren dependencias extra).
#    pAUC usa corrección McClish (partial.auc.correct = TRUE en pROC). Rango normalizado: 0.5-1.