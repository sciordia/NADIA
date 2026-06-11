# =============================================================================
# Benchmarking de métodos de NORMALIZACIÓN con BERT + covariables
# Dataset: muestra controlada TMT (spike-in con ground truth conocido)
#
# Los módulos R/Benchmarking_Single.R y R/Benchmarking_Multiple.R son agnósticos
# al origen de los datos: solo consumen el `de_res`. Por tanto BERT y las
# covariables se aplican ÍNTEGRAMENTE dentro de process_proteomics() (sección 2b
# de R/Processing.R) y las llamadas de benchmarking no cambian.
# =============================================================================

# ===== 1. PREPROCESSING TMT =====
source("R/Preprocessing_TMT.R")

preprocessing <- preprocess_tmt(
  file_path       = "<RUTA_AL_REPORT_TMT>",          # <-- AJUSTAR
  condition_order = c("A", "B", "C", "D"),           # <-- AJUSTAR al diseño real
  export_dir      = "./results/Q25_TMTpro_TMT1y2_PD/benchmark/"
)

library(dplyr)

# ===== 2. covariate_df =====
# Debe contener:
#   Column    -> coincide con preprocessing$metadata$Coding (clave del merge)
#   TMT       -> asignación de batch para BERT (>= 2 valores únicos)
#   Injection -> covariable de diseño para el modelo limma (covariate_column)
#   Patient   -> réplica biológica para duplicateCorrelation (bio_replicate_column)
# NO incluir TMT como covariate_column (ya se elimina con BERT -> doble corrección).
covariate_df <- read.delim("<RUTA_COVARIATES_TMT>", stringsAsFactors = FALSE)  # <-- AJUSTAR
# Ejemplo de saneamiento (descomentar/ajustar a tus nombres de columna):
# covariate_df <- covariate_df |>
#   transmute(Column = Coding, TMT = TMT, Injection = Injection, Patient = Patient)

# ===== 3. species_df (Protein.IDs + Species) =====
species_df <- read.delim("<RUTA_SPECIES_TMT>", stringsAsFactors = FALSE) |>  # <-- AJUSTAR
  transmute(
    Protein.IDs = PG.ProteinGroups,
    Species     = PG.OrganismId
  )

# ===== 4. Ground truth (diseño spike-in del dataset TMT) =====
# AJUSTAR Comparison/Species/expected_logFC al diseño real del TMT controlado.
expected <- tibble::tribble(
  ~Comparison, ~Species, ~expected_logFC,
  "B-A", "YEAST", -0.58,
  "B-A", "ECOLI",  1.00,
  "C-A", "YEAST", -1.60,
  "C-A", "ECOLI",  1.58,
  "D-A", "YEAST", -3.30,
  "D-A", "ECOLI",  2.00,
  "C-B", "YEAST", -1.00,
  "C-B", "ECOLI",  0.58,
  "D-B", "YEAST", -3.30,
  "D-B", "ECOLI",  1.00,
  "D-C", "YEAST", -3.30,
  "D-C", "ECOLI",  0.41
)

# =============================================================================
# 5. LOOP DE COMBOS: normalización + BERT + covariables
# =============================================================================
source("R/Processing.R")
source("R/Benchmarking_Single.R")
source("R/Benchmarking_Multiple.R")

`%||%` <- function(a, b) if (is.null(a)) b else a

# Métodos de normalización a comparar (imputación y batch fijos para aislar el efecto de la norma)
combos <- list(
  list(norm = "cycloess",        imp = "combo", mar = "Impseq", mnar = "min"),
  list(norm = "Rlr",             imp = "combo", mar = "Impseq", mnar = "min"),
  list(norm = "vsn",             imp = "combo", mar = "Impseq", mnar = "min"),
  list(norm = "eqmedians",       imp = "combo", mar = "Impseq", mnar = "min"),
  list(norm = "quantile",        imp = "combo", mar = "Impseq", mnar = "min"),
  list(norm = "quantile.robust", imp = "combo", mar = "Impseq", mnar = "min"),
  list(norm = "GlobalMean",      imp = "combo", mar = "Impseq", mnar = "min"),
  list(norm = "GlobalMedian",    imp = "combo", mar = "Impseq", mnar = "min"),
  list(norm = "MAD",             imp = "combo", mar = "Impseq", mnar = "min"),
  list(norm = "medianNorm",      imp = "combo", mar = "Impseq", mnar = "min"),
  list(norm = "meanNorm",        imp = "combo", mar = "Impseq", mnar = "min")
)

# Helper para nombrar el Assay (sufijo _BERT para distinguir de corridas sin batch correction)
make_assay_name <- function(combo) {
  base <- if (combo$imp == "none") {
    paste0(combo$norm, "_none")
  } else if (combo$imp == "combo") {
    paste(combo$norm, combo$mar, combo$mnar, sep = "_")
  } else if (combo$imp == "softHybrid") {
    paste(combo$norm, "softHybrid", combo$mar, combo$mnar, sep = "_")
  } else if (combo$imp == "limpa") {
    paste0(combo$norm, "_limpa")
  } else {
    paste(combo$norm, combo$imp, sep = "_")
  }
  paste0(base, "_BERT")
}

opdea_list      <- list()
confusion_list  <- list()
classified_list <- list()
bench_met_list  <- list()

for (combo in combos) {
  assay_name <- make_assay_name(combo)
  message("\n================ ", assay_name, " ================")

  combo_ok <- tryCatch({
  result <- process_proteomics(
    preprocessing           = preprocessing,
    min_reps_filter         = 3,
    norm_method             = combo$norm,
    cyclic_loess_method     = "fast",
    cyclic_loess_iterations = 3,
    cyclic_loess_span       = 0.7,
    # --- BERT: idéntico en TODOS los combos (corre DESPUÉS de cada normalización) ---
    batch_correct           = TRUE,
    batch_column            = "TMT",
    batch_algorithm         = "ComBat",
    batch_ComBat_mode       = 1,
    batch_covariates        = NULL,        # pocas muestras por combinación
    # --- imputación ---
    imp_method              = combo$imp,
    mar_method              = combo$mar,
    mnar_method             = combo$mnar,
    # --- modelo DE con covariables ---
    de_method               = combo$de_method %||% "limma",
    control                 = NULL,
    alpha                   = 0.05,
    covariate_df            = covariate_df,
    covariate_column        = "Injection",  # covariable de diseño (limma)
    bio_replicate_column    = "Patient",    # duplicateCorrelation
    export_format           = "both"
  )

  bench <- benchmarking_proteomics(
    de_res          = result$DEPs_results,
    species_df      = species_df,
    expected_values = expected,
    alpha           = 0.05,
    output_dir      = paste0("results/Q25_TMTpro/benchmark_", assay_name)
  )

  opdea <- bench$opdea_metrics;        opdea$Assay <- assay_name
  conf  <- bench$confusion_overall;    conf$Assay  <- assay_name
  cls   <- bench$classified_df;        cls$Assay   <- assay_name
  met   <- bench$metrics_table;        met$Assay   <- assay_name

  opdea_list[[length(opdea_list) + 1]]           <- opdea
  confusion_list[[length(confusion_list) + 1]]   <- conf
  classified_list[[length(classified_list) + 1]] <- cls
  bench_met_list[[length(bench_met_list) + 1]]   <- met
  TRUE
  }, error = function(e) {
    message("  [SALTADO] ", assay_name, " falló: ", conditionMessage(e))
    FALSE
  })
  if (!isTRUE(combo_ok)) next
}

opdea_all      <- do.call(rbind, opdea_list)
confusion_all  <- do.call(rbind, confusion_list)
classified_all <- do.call(rbind, classified_list)
bench_met_all  <- do.call(rbind, bench_met_list)

# =============================================================================
# 6. RANKING MULTI-MÉTODO (sin cambios respecto al flujo estándar)
# =============================================================================
bm_result <- benchmarking_multiple(
  opdea_combined         = opdea_all,
  confusion_combined     = confusion_all,
  classified_combined    = classified_all,
  bench_metrics_combined = bench_met_all,
  output_dir             = "results/Q25_TMTpro/bm_multiple_norm_BERT_ALL",
  verbose                = TRUE
)

# --- Resultados clave ---
bm_result$mean_ranking
bm_result$gg_ranking_heatmap_mean
bm_result$gg_metrics_comparison
bm_result$extended_ranking
bm_result$gg_extended_ranking_bars
bm_result$gg_extended_ranking_heatmap

# ROC por comparación (cordura del ground truth)
# bm_result$gg_roc_zoom_by_comp[["B-A"]]
