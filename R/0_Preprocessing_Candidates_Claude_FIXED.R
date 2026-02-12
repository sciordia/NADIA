################################################################################
# ANÁLISIS PROTEÓMICO COMPLETO CON MÉTRICAS DE CLASIFICACIÓN - VERSIÓN FINAL
# Interpretación correcta: Positivo = Cambio Significativo, Negativo = Sin Cambio
# ECOLI/YEAST esperan ser POSITIVOS (TP), HUMAN espera ser NEGATIVO (TN)
################################################################################

# ══════════════════════════════════════════════════════════════════════════════
# 1) PAQUETES
# ══════════════════════════════════════════════════════════════════════════════
packages <- c(
  "dplyr","stringr","tidyr","rlang","ggplot2","tibble",
  "pROC","readr","scales","gridExtra","RColorBrewer","paletteer"
)

to_install <- packages[!(packages %in% installed.packages()[,"Package"])]
if (length(to_install)) install.packages(to_install, dependencies = TRUE)

library(dplyr)
library(stringr)
library(tidyr)
library(rlang)
library(ggplot2)
library(tibble)
library(pROC)
library(readr)
library(scales)
library(gridExtra)
library(RColorBrewer)

# SETTINGS
## Orden de las comparaciones mostradas en los gráficos
orden_comparaciones <- c("B/A", "C/A", "D/A", "C/B", "D/B", "D/C")
fc_threshold  <- log2(1)          # umbral mínimo de efecto
alpha <- 0.05                      # umbral de significancia
infile <- "./Q24_DIA/outputs_Exploris_Spectronaut_v20_5comp/DE_ALL_with_LoessCyc_bySample_Candidates_FIXED.tsv"
outdir <- "./Q24_DIA/outputs_Exploris_Spectronaut_v20_5comp/"

# ══════════════════════════════════════════════════════════════════════════════
# 2) FUNCIONES AUXILIARES
# ══════════════════════════════════════════════════════════════════════════════

# Especie desde ProteinNames (vectorizada, soporta múltiples entradas)
collapse_species <- function(x) {
  has_human <- grepl("_HUMAN\\b", x, ignore.case = TRUE)
  has_bovin <- grepl("_BOVIN\\b", x, ignore.case = TRUE)
  has_ecoli <- grepl("_ECOLI\\b", x, ignore.case = TRUE)
  has_yeast <- grepl("_YEAST\\b", x, ignore.case = TRUE)
  
  # HUMAN gana a BOVIN si coinciden
  has_bovin <- has_bovin & !has_human
  
  n_species <- has_human + has_bovin + has_ecoli + has_yeast
  
  out <- ifelse(n_species == 0, NA_character_,
                ifelse(n_species > 1, "MIX",
                       ifelse(has_human, "HUMAN",
                              ifelse(has_ecoli, "ECOLI",
                                     ifelse(has_yeast, "YEAST", "BOVIN")))))
  out
}

# Normaliza comparaciones a formato "A/B"
normalize_comparison <- function(x) {
  x %>%
    str_trim() %>%
    str_replace_all("(?i)\\bvs\\b", "/") %>%
    str_replace_all("[-_\\\\]", "/") %>%
    str_replace_all("\\s*/\\s*", "/") %>%
    str_replace_all("\\s+", " ") %>%
    toupper()
}

# Resumen por comparación (conteos UP/DOWN por especie)
summarize_proteins_by_comparison <- function(df, alpha = 0.05) {
  df %>%
    group_by(Comparison) %>%
    summarise(
      n_nonsignif = sum(is.na(Qvalue) | Qvalue > alpha, na.rm = TRUE),
      
      n_signif_up_ecoli = sum(!is.na(Qvalue) & Qvalue <= alpha & !is.na(`AVG Log2 Ratio`) &
                                `AVG Log2 Ratio` > 0 & Species == "ECOLI", na.rm = TRUE),
      n_signif_up_human = sum(!is.na(Qvalue) & Qvalue <= alpha & !is.na(`AVG Log2 Ratio`) &
                                `AVG Log2 Ratio` > 0 & Species == "HUMAN", na.rm = TRUE),
      n_signif_up_yeast = sum(!is.na(Qvalue) & Qvalue <= alpha & !is.na(`AVG Log2 Ratio`) &
                                `AVG Log2 Ratio` > 0 & Species == "YEAST", na.rm = TRUE),
      
      n_signif_down_ecoli = sum(!is.na(Qvalue) & Qvalue <= alpha & !is.na(`AVG Log2 Ratio`) &
                                  `AVG Log2 Ratio` < 0 & Species == "ECOLI", na.rm = TRUE),
      n_signif_down_human = sum(!is.na(Qvalue) & Qvalue <= alpha & !is.na(`AVG Log2 Ratio`) &
                                  `AVG Log2 Ratio` < 0 & Species == "HUMAN", na.rm = TRUE),
      n_signif_down_yeast = sum(!is.na(Qvalue) & Qvalue <= alpha & !is.na(`AVG Log2 Ratio`) &
                                  `AVG Log2 Ratio` < 0 & Species == "YEAST", na.rm = TRUE),
      
      n_signif_total = sum(!is.na(Qvalue) & Qvalue <= alpha, na.rm = TRUE),
      .groups = "drop"
    )
}

# SD/CV con trimming
trimmed_sd_cv <- function(x, trim = 0.1) {
  x <- sort(x); n <- length(x)
  if (n == 0) return(c(sd_trim = NA_real_, cv_trim = NA_real_))
  lo <- floor(n * trim) + 1
  hi <- max(lo, ceiling(n * (1 - trim)))
  x_trim <- x[lo:hi]
  sd_t  <- sd(x_trim, na.rm = TRUE)
  med_t <- median(x_trim, na.rm = TRUE)
  cv_t  <- sd_t / max(abs(med_t), 1e-8)
  c(sd_trim = sd_t, cv_trim = cv_t)
}

# ══════════════════════════════════════════════════════════════════════════════
# 3) CARGA Y PREPARACIÓN
# ══════════════════════════════════════════════════════════════════════════════

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("ANÁLISIS PROTEÓMICO - VERSIÓN FINAL\n")
cat("Interpretación: Positivo = Cambio Significativo\n")
cat("               Negativo = Sin Cambio\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

df <- read.delim(infile, sep = "\t", check.names = FALSE)
cat("Dimensiones del dataset:", nrow(df), "filas x", ncol(df), "columnas\n\n")

# Aplicar especie y filtros
df <- df %>%
  mutate(Species = collapse_species(ProteinNames)) %>%
  filter(!is.na(Species) & Species != "MIX") %>%
  mutate(Species = ifelse(Species == "BOVIN", "HUMAN", Species)) %>%
  filter(abs(`AVG Log2 Ratio`) >= fc_threshold)

cat("Distribución de especies después del filtrado:\n")
print(table(df$Species)); cat("\n")

# Detección robusta del nombre de columna de comparación
cmp_col <- grep("(?i)comparison.*group", names(df), value = TRUE)
if (length(cmp_col) == 0) cmp_col <- grep("(?i)comparison", names(df), value = TRUE)
if (length(cmp_col) == 0) stop("❌ No se encontró ninguna columna de comparación en el dataset")
if (length(cmp_col) > 1) {
  message("⚠️ Varias columnas coinciden: ", paste(cmp_col, collapse = ", "))
  message("   Se usará la primera: ", cmp_col[1])
}
cmp_col <- cmp_col[1]

# Renombrado correcto usando nombre dinámico
df <- df %>%
  dplyr::rename(Comparison_raw = !!rlang::sym(cmp_col)) %>%
  dplyr::mutate(Comparison = normalize_comparison(Comparison_raw))

cat("✅ Columna de comparación renombrada desde: ", cmp_col, " → 'Comparison_raw'\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# 4) DEFINICIÓN DE VERDAD TERRENO (Ground Truth) - VERSIÓN CORREGIDA
# ══════════════════════════════════════════════════════════════════════════════
# ECOLI: Esperamos POSITIVOS (cambio UP)
# YEAST: Esperamos POSITIVOS (cambio DOWN)  
# HUMAN: Esperamos NEGATIVOS (sin cambio)

cat("📊 DEFINICIÓN DE GROUND TRUTH:\n")
cat("   ECOLI → Esperamos POSITIVOS (cambios UP significativos)\n")
cat("   YEAST → Esperamos POSITIVOS (cambios DOWN significativos)\n")
cat("   HUMAN → Esperamos NEGATIVOS (sin cambios significativos)\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# 5) DATASET CON VERDAD Y PREDICCIONES - INTERPRETACIÓN CORRECTA
# ══════════════════════════════════════════════════════════════════════════════

df_truth <- df %>%
  mutate(
    # TRUTH: ¿Esperamos un cambio significativo?
    # 1 = Esperamos cambio (positivo), 0 = NO esperamos cambio (negativo)
    truth = case_when(
      Species == "ECOLI" ~ 1,  # Esperamos cambio (positivo)
      Species == "YEAST" ~ 1,  # Esperamos cambio (positivo)
      Species == "HUMAN" ~ 0,  # NO esperamos cambio (negativo)
      TRUE ~ NA_real_
    ),
    
    # PREDICTED: ¿Detectamos un cambio significativo?
    # 1 = Cambio detectado, 0 = Sin cambio detectado
    predicted = case_when(
      # ECOLI: cambio significativo si Q≤0.05 Y dirección correcta (UP)
      Species == "ECOLI" & Qvalue <= alpha & `AVG Log2 Ratio` > 0 ~ 1,
      Species == "ECOLI" & (Qvalue > alpha | `AVG Log2 Ratio` <= 0) ~ 0,
      
      # YEAST: cambio significativo si Q≤0.05 Y dirección correcta (DOWN)
      Species == "YEAST" & Qvalue <= alpha & `AVG Log2 Ratio` < 0 ~ 1,
      Species == "YEAST" & (Qvalue > alpha | `AVG Log2 Ratio` >= 0) ~ 0,
      
      # HUMAN: cambio significativo si Q≤0.05 (cualquier dirección)
      Species == "HUMAN" & Qvalue <= alpha ~ 1,
      Species == "HUMAN" & Qvalue > alpha ~ 0,
      
      TRUE ~ NA_real_
    ),
    
    # Score para ROC (mayor score = más significativo)
    score = -log10(pmax(Qvalue, 1e-300)),
    
    # Clasificación para debugging
    classification = case_when(
      predicted == 1 & truth == 1 ~ "TP",
      predicted == 1 & truth == 0 ~ "FP", 
      predicted == 0 & truth == 0 ~ "TN",
      predicted == 0 & truth == 1 ~ "FN",
      TRUE ~ NA_character_
    )
  )

cat("✅ df_truth creado con interpretación correcta\n")
cat("   Positivo = Cambio Significativo Esperado\n")
cat("   Negativo = Sin Cambio Esperado\n\n")

# Verificación de la lógica
cat("📋 VERIFICACIÓN DE CLASIFICACIONES:\n")
verification_table <- df_truth %>%
  group_by(Species, classification) %>%
  summarise(
    n = n(),
    median_qvalue = median(Qvalue, na.rm = TRUE),
    median_fc = median(`AVG Log2 Ratio`, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(Species, classification)

print(verification_table)
cat("\n")

# ══════════════════════════════════════════════════════════════════════════════
# 6) MÉTRICAS POR COMPARACIÓN
# ══════════════════════════════════════════════════════════════════════════════

metrics_table <- df_truth %>%
  filter(!is.na(truth) & !is.na(predicted)) %>%
  group_by(Comparison) %>%
  summarise(
    TP = sum(classification == "TP", na.rm = TRUE),
    FP = sum(classification == "FP", na.rm = TRUE),
    TN = sum(classification == "TN", na.rm = TRUE),
    FN = sum(classification == "FN", na.rm = TRUE),
    
    # Métricas derivadas
    Sensitivity = ifelse((TP + FN) > 0, TP / (TP + FN), NA_real_),  # TPR
    Specificity = ifelse((TN + FP) > 0, TN / (TN + FP), NA_real_),  # TNR
    Precision   = ifelse((TP + FP) > 0, TP / (TP + FP), NA_real_),  # PPV
    NPV         = ifelse((TN + FN) > 0, TN / (TN + FN), NA_real_),
    Accuracy    = ifelse((TP + FP + TN + FN) > 0, 
                         (TP + TN) / (TP + FP + TN + FN), NA_real_),
    F1_Score    = ifelse(Precision > 0 & Sensitivity > 0,
                         2 * (Precision * Sensitivity) / (Precision + Sensitivity),
                         NA_real_),
    MCC = {
      tp <- as.numeric(TP); fp <- as.numeric(FP)
      tn <- as.numeric(TN); fn <- as.numeric(FN)
      den <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
      ifelse(den > 0, (tp * tn - fp * fn) / den, NA_real_)
    },
    Total_Samples = TP + FP + TN + FN,
    .groups = "drop"
  )

# Cálculo de AUC
auc_table <- df_truth %>%
  filter(!is.na(truth) & !is.na(score)) %>%
  group_by(Comparison) %>%
  summarise(
    AUC = {
      d <- data.frame(truth = truth, score = score)
      d <- d[complete.cases(d), ]
      if (length(unique(d$truth)) < 2 || nrow(d) < 10) NA_real_ else
        tryCatch(as.numeric(pROC::roc(d$truth, d$score, quiet = TRUE)$auc),
                 error = function(e) NA_real_)
    },
    .groups = "drop"
  )

summary_table <- auc_table %>% 
  left_join(metrics_table, by = "Comparison") %>% 
  arrange(desc(AUC))

cat("✅ Métricas calculadas para ", nrow(summary_table), " comparaciones\n\n", sep = "")

# Interpretación de métricas
cat("📊 INTERPRETACIÓN DE MÉTRICAS:\n")
cat("   Sensibilidad: Capacidad de detectar cambios verdaderos (ECOLI↑, YEAST↓)\n")
cat("   Especificidad: Capacidad de identificar no-cambios (HUMAN)\n")
cat("   Precisión: De los cambios detectados, cuántos son verdaderos\n")
cat("   FP: Principalmente HUMAN incorrectamente clasificados como cambiados\n\n")

# Mostrar resumen
cat("📈 RESUMEN DE MÉTRICAS PRINCIPALES:\n")
summary_table %>%
  select(Comparison, AUC, Sensitivity, Specificity, Precision, Accuracy) %>%
  mutate(across(where(is.numeric), ~round(., 3))) %>%
  print()
cat("\n")

# ══════════════════════════════════════════════════════════════════════════════
# 7) SALIDAS: DIRECTORIO Y TABLAS
# ══════════════════════════════════════════════════════════════════════════════

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

write_tsv(summary_table, file.path(outdir, "metricas_completas_por_comparacion.tsv"))
cat("📁 Guardado: ", file.path(outdir, "metricas_completas_por_comparacion.tsv"), "\n\n", sep = "")

# Tabla de debugging por especie
debug_table <- df_truth %>%
  filter(!is.na(truth) & !is.na(predicted)) %>%
  group_by(Comparison, Species) %>%
  summarise(
    n_total = n(),
    n_truth_positive = sum(truth == 1),
    n_truth_negative = sum(truth == 0),
    n_predicted_positive = sum(predicted == 1),
    n_predicted_negative = sum(predicted == 0),
    n_TP = sum(classification == "TP"),
    n_FP = sum(classification == "FP"),
    n_TN = sum(classification == "TN"),
    n_FN = sum(classification == "FN"),
    accuracy = (n_TP + n_TN) / n_total,
    .groups = "drop"
  )

write_tsv(debug_table, file.path(outdir, "debug_confusion_by_species.tsv"))
cat("📁 Debug table guardada: ", file.path(outdir, "debug_confusion_by_species.tsv"), "\n\n", sep = "")

# ══════════════════════════════════════════════════════════════════════════════
# 8) HEATMAP DE MÉTRICAS
# ══════════════════════════════════════════════════════════════════════════════

metrics_for_heatmap <- summary_table %>%
  select(Comparison, AUC, Sensitivity, Specificity, Precision, F1_Score, Accuracy) %>%
  pivot_longer(-Comparison, names_to = "Metric", values_to = "Value") %>%
  mutate(
    Metric = factor(Metric, levels = c("AUC","Sensitivity","Specificity",
                                       "Precision","F1_Score","Accuracy")),
    Comparison = factor(Comparison, levels = rev(orden_comparaciones))
  )

p_heatmap <- ggplot(metrics_for_heatmap, aes(x = Metric, y = Comparison, fill = Value)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(is.na(Value), "NA", sprintf("%.3f", Value))), 
            color = "black", size = 5) +
  scale_fill_gradient2(low = "#d73027", mid = "#fee090", high = "#1a9850",
                       midpoint = 0.5, limits = c(0,1), na.value = "grey90", 
                       name = "Valor") +
  labs(title = "Heatmap de Métricas de Clasificación por Comparación",
       subtitle = "Positivo=Cambio Esperado (ECOLI/YEAST), Negativo=Sin Cambio (HUMAN)",
       x = "Métrica", y = "Comparación") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.title = element_text(face = "bold", size = 14),
        panel.grid = element_blank(),
        legend.position = "right")

ggsave(file.path(outdir, "heatmap_metricas.png"), p_heatmap, width = 10, height = 6, dpi = 300)
print(p_heatmap)

# ══════════════════════════════════════════════════════════════════════════════
# 9) CURVAS ROC
# ══════════════════════════════════════════════════════════════════════════════

roc_list <- lapply(split(df_truth, df_truth$Comparison), function(d) {
  d <- d[complete.cases(d$truth, d$score), ]
  if (length(unique(d$truth)) < 2 || nrow(d) < 10) return(NULL)
  tryCatch(pROC::roc(response = d$truth, predictor = d$score, quiet = TRUE),
           error = function(e) NULL)
})

roc_list <- roc_list[!vapply(roc_list, is.null, logical(1))]

if (length(roc_list) > 0) {
  aucs <- vapply(roc_list, function(r) as.numeric(pROC::auc(r)), numeric(1))
  labels <- paste0(names(roc_list), " (AUC=", sprintf("%.3f", aucs), ")")
  
  p_roc <- pROC::ggroc(roc_list, legacy.axes = TRUE, size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
    scale_color_brewer(palette = "Set1", labels = labels) +
    labs(title = "Curvas ROC por Comparación",
         subtitle = "Sensibilidad vs (1-Especificidad). Diagonal = clasificador aleatorio",
         x = "False Positive Rate (1 - Specificity)", 
         y = "True Positive Rate (Sensitivity)",
         color = "Comparación") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 14),
          legend.position = "right",
          legend.text = element_text(size = 9),
          panel.grid.minor = element_blank()) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1))
  
  ggsave(file.path(outdir, "curvas_roc_completas.png"), p_roc, width = 10, height = 7, dpi = 300)
  print(p_roc)
  
  # ROC con zoom
  p_roc_zoom <- p_roc +
    coord_cartesian(xlim = c(0, 0.1), ylim = c(0, 1)) +
    labs(title = "Curvas ROC por Comparación - ZOOM",
         subtitle = "Vista ampliada de la región de baja tasa de falsos positivos")
  
  ggsave(file.path(outdir, "curvas_roc_zoom.png"), p_roc_zoom, width = 10, height = 7, dpi = 300)
  print(p_roc_zoom)
}

# ══════════════════════════════════════════════════════════════════════════════
# 10) BARRAS DE AUC
# ══════════════════════════════════════════════════════════════════════════════

p_auc_bars <- ggplot(summary_table, aes(x = reorder(Comparison, AUC), y = AUC, fill = AUC)) +
  geom_col() +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "red", linewidth = 0.5) +
  geom_text(aes(label = ifelse(is.na(AUC), "NA", sprintf("%.3f", AUC))),
            hjust = -0.1, size = 3) +
  scale_fill_gradient2(low = "#d73027", mid = "#fee090", high = "#1a9850",
                       midpoint = 0.75, limits = c(0,1)) +
  coord_flip() +
  ylim(0, 1.05) +
  labs(title = "AUC por Comparación",
       subtitle = "Área bajo la curva ROC. Línea roja = clasificador aleatorio (0.5)",
       x = "Comparación", y = "AUC", fill = "AUC") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14),
        legend.position = "right")

ggsave(file.path(outdir, "barras_auc.png"), p_auc_bars, width = 8, height = 6, dpi = 300)
print(p_auc_bars)

# ══════════════════════════════════════════════════════════════════════════════
# 11) COMPARACIÓN DE MÚLTIPLES MÉTRICAS
# ══════════════════════════════════════════════════════════════════════════════

metrics_comparison <- summary_table %>%
  select(Comparison, Sensitivity, Specificity, Precision, F1_Score, Accuracy) %>%
  pivot_longer(-Comparison, names_to = "Metric", values_to = "Value") %>%
  mutate(Metric = factor(Metric, levels = c("Sensitivity","Specificity",
                                            "Precision","F1_Score","Accuracy")))

p_metrics_comparison <- ggplot(metrics_comparison, aes(x = Comparison, y = Value, fill = Metric)) +
  geom_col(position = "dodge", color = "black", linewidth = 0.2) +
  scale_fill_brewer(palette = "Set2") +
  coord_flip() + ylim(0,1) +
  labs(title = "Comparación de Métricas de Clasificación",
       subtitle = "Sens=Detección cambios (ECOLI/YEAST), Spec=Detección no-cambios (HUMAN)",
       x = "Comparación", y = "Valor de la Métrica", fill = "Métrica") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 14),
        legend.position = "bottom",
        panel.grid.major.y = element_blank())

ggsave(file.path(outdir, "comparacion_metricas.png"), p_metrics_comparison, width = 10, height = 7, dpi = 300)
print(p_metrics_comparison)

# ══════════════════════════════════════════════════════════════════════════════
# 12) MATRIZ DE CONFUSIÓN VISUAL
# ══════════════════════════════════════════════════════════════════════════════

confusion_matrices <- summary_table %>%
  select(Comparison, TP, FP, TN, FN) %>%
  mutate(
    Total = TP + FP + TN + FN,
    TP_pct = TP / pmax(Total, 1) * 100,
    FP_pct = FP / pmax(Total, 1) * 100,
    TN_pct = TN / pmax(Total, 1) * 100,
    FN_pct = FN / pmax(Total, 1) * 100
  ) %>%
  select(Comparison, TP_pct, FP_pct, TN_pct, FN_pct) %>%
  pivot_longer(-Comparison, names_to = "Category", values_to = "Percentage") %>%
  mutate(
    Category = factor(Category,
                      levels = c("TP_pct","FP_pct","FN_pct","TN_pct"),
                      labels = c("True Positive\n(Cambios Detectados)",
                                 "False Positive\n(HUMAN mal clasificado)",
                                 "False Negative\n(Cambios perdidos)",
                                 "True Negative\n(HUMAN correcto)")),
    Comparison = factor(Comparison, levels = rev(orden_comparaciones))
  )

p_confusion <- ggplot(confusion_matrices, aes(x = Category, y = Comparison, fill = Percentage)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.1f%%", Percentage)), color = "black", size = 4) +
  scale_fill_gradient2(low = "white", high = "#3182bd", limits = c(0,100)) +
  labs(title = "Distribución de Clasificaciones - Matriz de Confusión",
       subtitle = "TP=ECOLI/YEAST detectados, TN=HUMAN sin cambio, FP=HUMAN con cambio falso",
       x = "Categoría de Clasificación", y = "Comparación", fill = "Porcentaje") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        plot.title = element_text(face = "bold", size = 14),
        panel.grid = element_blank())

ggsave(file.path(outdir, "matriz_confusion.png"), p_confusion, width = 10, height = 6, dpi = 300)
print(p_confusion)

# ══════════════════════════════════════════════════════════════════════════════
# 13) ANÁLISIS DE PROTEÍNAS SIGNIFICATIVAS
# ══════════════════════════════════════════════════════════════════════════════

resumen <- summarize_proteins_by_comparison(df, alpha)

resumen_ext <- resumen %>%
  mutate(
    total_proteins      = n_nonsignif + n_signif_total,
    pct_signif          = 100 * n_signif_total / pmax(total_proteins, 1),
    n_signif_up_total   = n_signif_up_ecoli + n_signif_up_human + n_signif_up_yeast,
    n_signif_down_total = n_signif_down_ecoli + n_signif_down_human + n_signif_down_yeast,
    up_down_ratio       = ifelse(n_signif_down_total == 0, NA_real_, 
                                 n_signif_up_total / n_signif_down_total)
  )

print(resumen_ext)
write_tsv(resumen_ext, file.path(outdir, "resumen_proteinas_significativas.tsv"))
cat("\n📁 Guardado: ", file.path(outdir, "resumen_proteinas_significativas.tsv"), "\n\n", sep = "")

# Gráfico de proteínas significativas
plot_df <- resumen %>%
  transmute(
    Comparison,
    up_ECOLI   = n_signif_up_ecoli,
    up_HUMAN   = n_signif_up_human,
    up_YEAST   = n_signif_up_yeast,
    down_ECOLI = n_signif_down_ecoli,
    down_HUMAN = n_signif_down_human,
    down_YEAST = n_signif_down_yeast
  ) %>%
  pivot_longer(
    cols = -Comparison,
    names_to = "dir_species",
    values_to = "n"
  ) %>%
  mutate(
    Direction = ifelse(grepl("^up_", dir_species), "UP", "DOWN"),
    Species   = toupper(sub("^(up_|down_)", "", dir_species)),
    Species   = factor(Species, levels = c("HUMAN", "ECOLI", "YEAST"))
  )

plot_df <- plot_df %>%
  group_by(Comparison, Direction) %>%
  mutate(
    total_dir = sum(n, na.rm = TRUE),
    pct = ifelse(total_dir > 0, 100 * n / total_dir, 0)
  ) %>%
  ungroup() %>%
  mutate(Comparison = factor(Comparison, levels = rev(orden_comparaciones)))

# Gráfico de números absolutos
p_signif_n <- ggplot(plot_df, aes(x = Comparison, y = n, fill = Species)) +
  geom_col() +
  facet_wrap(~ Direction, nrow = 2, scales = "free_y") +
  scale_fill_manual(values = c("HUMAN" = "#E41A1C", "ECOLI" = "#377EB8", "YEAST" = "#4DAF4A")) +
  labs(
    title = "Proteínas Significativas por Comparación y Especie",
    subtitle = "HUMAN (rojo) = FP esperados, ECOLI/YEAST (azul/verde) = TP esperados",
    x = "Comparación",
    y = "Número de proteínas significativas",
    fill = "Especie"
  ) +
  coord_flip() +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "bold", size = 11)
  )

ggsave(file.path(outdir, "proteinas_significativas_barras_totales.png"),
       p_signif_n, width = 10, height = 8, dpi = 300)
print(p_signif_n)

# ══════════════════════════════════════════════════════════════════════════════
# 14) ANÁLISIS DE DISPERSIÓN POR ESPECIE
# ══════════════════════════════════════════════════════════════════════════════

# Filtro para proteínas significativas con dirección esperada
df_sig_dir <- df_truth %>%
  filter(
    Species %in% c("ECOLI", "YEAST"),
    !is.na(Qvalue), Qvalue <= alpha,
    !is.na(`AVG Log2 Ratio`),
    (Species == "ECOLI" & `AVG Log2 Ratio` > 0) |
      (Species == "YEAST" & `AVG Log2 Ratio` < 0)
  )

# Métricas de dispersión
metrics_by_cmp_species <- df_sig_dir %>%
  group_by(Comparison, Species) %>%
  summarise(
    n              = n(),
    median_log2    = median(`AVG Log2 Ratio`, na.rm = TRUE),
    sd_log2        = sd(`AVG Log2 Ratio`, na.rm = TRUE),
    cv_log2        = sd_log2 / pmax(abs(median_log2), 1e-8),
    mad_log2       = mad(`AVG Log2 Ratio`, center = median_log2, constant = 1, na.rm = TRUE),
    rcv_mad        = mad_log2 / pmax(abs(median_log2), 1e-8),
    iqr_log2       = IQR(`AVG Log2 Ratio`, na.rm = TRUE),
    sd_cv_trim     = list(trimmed_sd_cv(`AVG Log2 Ratio`, trim = 0.1)),
    .groups = "drop"
  ) %>%
  mutate(
    # Extraer valores de la lista de manera segura
    sd_trim_log2 = sapply(sd_cv_trim, function(x) {
      if(is.null(x) || length(x) == 0) return(NA_real_)
      x[["sd_trim"]]
    }),
    cv_trim_log2 = sapply(sd_cv_trim, function(x) {
      if(is.null(x) || length(x) == 0) return(NA_real_)
      x[["cv_trim"]]
    })
  ) %>%
  select(-sd_cv_trim)

# Ranking de homogeneidad
homogeneity_ranking <- metrics_by_cmp_species %>%
  group_by(Comparison) %>%
  summarise(
    n_total       = sum(n, na.rm = TRUE),
    mean_cv       = mean(cv_log2, na.rm = TRUE),
    mean_rcv_mad  = mean(rcv_mad, na.rm = TRUE),
    mean_cv_trim  = mean(cv_trim_log2, na.rm = TRUE),
    mean_iqr      = mean(iqr_log2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(mean_cv)

write_tsv(metrics_by_cmp_species, file.path(outdir, "dispersion_cmp_species_signif_dir.tsv"))
write_tsv(homogeneity_ranking, file.path(outdir, "dispersion_ranking_comparisons.tsv"))

cat("📁 Guardado análisis de dispersión:\n")
cat("   - ", file.path(outdir, "dispersion_cmp_species_signif_dir.tsv"), "\n")
cat("   - ", file.path(outdir, "dispersion_ranking_comparisons.tsv"), "\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# 15) TABLA RESUMEN FINAL
# ══════════════════════════════════════════════════════════════════════════════

summary_final <- summary_table %>%
  mutate(
    across(c(AUC, Sensitivity, Specificity, Precision, NPV, Accuracy, F1_Score, MCC), 
           ~ round(., 3)),
    Performance = case_when(
      AUC >= 0.9 ~ "Excelente",
      AUC >= 0.8 ~ "Muy Bueno",
      AUC >= 0.7 ~ "Bueno",
      AUC >= 0.6 ~ "Aceptable",
      TRUE       ~ "Pobre"
    ),
    # Añadir interpretación
    Interpretation = paste0(
      "Sens=", round(Sensitivity, 2), " (cambios detectados), ",
      "Spec=", round(Specificity, 2), " (no-cambios correctos)"
    )
  ) %>%
  select(Comparison, AUC, Performance, Sensitivity, Specificity, Precision, 
         F1_Score, Accuracy, MCC, Interpretation, everything())

cat("\n════════════════════════════════════════════════════════════════\n")
cat("TABLA RESUMEN FINAL\n")
cat("════════════════════════════════════════════════════════════════\n")
print(summary_final %>% select(-Interpretation))

write.csv(summary_final, file.path(outdir, "metricas_completas_resumen_final.csv"), 
          row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# 16) FUNCIÓN DE VALIDACIÓN FINAL
# ══════════════════════════════════════════════════════════════════════════════

validate_results <- function() {
  cat("\n\n════════════════════════════════════════════════════════════════\n")
  cat("VALIDACIÓN DE RESULTADOS FINALES\n")
  cat("════════════════════════════════════════════════════════════════\n\n")
  
  # Verificación por especie
  cat("📊 Verificación de Clasificaciones por Especie:\n")
  cat("─────────────────────────────────────────────\n")
  
  validation_summary <- df_truth %>%
    filter(!is.na(truth) & !is.na(predicted)) %>%
    group_by(Species) %>%
    summarise(
      n_total = n(),
      
      # Expected
      expected_label = case_when(
        Species %in% c("ECOLI", "YEAST") ~ "POSITIVO (cambio)",
        Species == "HUMAN" ~ "NEGATIVO (sin cambio)",
        TRUE ~ "Unknown"
      ),
      n_expected_positive = sum(truth == 1),
      n_expected_negative = sum(truth == 0),
      
      # Predicted
      n_predicted_positive = sum(predicted == 1),
      n_predicted_negative = sum(predicted == 0),
      
      # Classification
      n_TP = sum(classification == "TP"),
      n_FP = sum(classification == "FP"),
      n_TN = sum(classification == "TN"),
      n_FN = sum(classification == "FN"),
      
      # Metrics
      accuracy = (n_TP + n_TN) / n_total,
      
      # Interpretation
      main_classification = case_when(
        Species == "ECOLI" ~ paste0("TP=", n_TP, ", FN=", n_FN),
        Species == "YEAST" ~ paste0("TP=", n_TP, ", FN=", n_FN),
        Species == "HUMAN" ~ paste0("TN=", n_TN, ", FP=", n_FP),
        TRUE ~ "Unknown"
      ),
      
      .groups = "drop"
    )
  
  print(validation_summary %>% select(Species, expected_label, n_total, main_classification, accuracy))
  
  cat("\n✅ CONFIRMACIÓN DE INTERPRETACIÓN:\n")
  cat("   • ECOLI: Esperamos POSITIVOS → TP cuando detectados, FN cuando perdidos\n")
  cat("   • YEAST: Esperamos POSITIVOS → TP cuando detectados, FN cuando perdidos\n")
  cat("   • HUMAN: Esperamos NEGATIVOS → TN cuando sin cambio, FP cuando cambio falso\n")
  
  # Ejemplos específicos
  cat("\n📋 Ejemplos de Clasificación (5 de cada tipo):\n")
  cat("─────────────────────────────────────────────\n")
  
  examples <- df_truth %>%
    filter(!is.na(classification)) %>%
    group_by(Species, classification) %>%
    slice_head(n = 2) %>%
    select(Species, Qvalue, `AVG Log2 Ratio`, truth, predicted, classification) %>%
    mutate(
      Qvalue = round(Qvalue, 4),
      `AVG Log2 Ratio` = round(`AVG Log2 Ratio`, 2)
    )
  
  print(examples, n = Inf)
  
  cat("\n════════════════════════════════════════════════════════════════\n")
}

# Ejecutar validación
validate_results()

# ══════════════════════════════════════════════════════════════════════════════
# FIN DEL ANÁLISIS
# ══════════════════════════════════════════════════════════════════════════════

cat("\n✅ ANÁLISIS COMPLETADO - VERSIÓN FINAL CORREGIDA\n")
cat("\n📋 RESUMEN DE INTERPRETACIÓN:\n")
cat("   • Positivo = Cambio Significativo Esperado (ECOLI↑, YEAST↓)\n")
cat("   • Negativo = Sin Cambio Esperado (HUMAN)\n")
cat("   • Sensibilidad = Detección de cambios verdaderos\n")
cat("   • Especificidad = Identificación correcta de no-cambios\n")
cat("   • FP = Principalmente HUMAN con cambios falsos\n")
cat("   • TN = HUMAN correctamente sin cambios\n")

cat("\nArchivos generados en ", outdir, ":\n", sep = "")
list.files(outdir, pattern = "\\.(tsv|csv|png)$") %>% 
  paste("  -", .) %>% 
  cat(sep = "\n")

cat("\n\n═══════════════════════════════════════════════════════════════\n")
cat("FIN DEL ANÁLISIS - INTERPRETACIÓN CORRECTA\n")
cat("═══════════════════════════════════════════════════════════════\n")

