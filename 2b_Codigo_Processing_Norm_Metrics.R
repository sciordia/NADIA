
# =========================================================
# Código de ejemplo: Normalization Metrics + Plots
# =========================================================

# Cargar las funciones para el Heatmap-Plot
source("./R/Normalization_Metrics.R")

# =============================================================================                                      
# PASO 1: PREPROCESAMIENTO DE DATOS DE SPECTRONAUT                                                                   
# =============================================================================                                      

# Importar matrices normalizadas                                                                                          
se_nm <- import_norm_matrices(
  tsv_dir       = "./results/tsv_norm/",
  metadata_path = "./results/Metadata.tsv",
  condition_col = "Condition",
  sample_col    = "Column"
)                                                                                                                    

# --- Explorar el objeto preprocessing ---                                                                           
plots_nm <- normalization_metrics(se_nm)

# Visualizar los plots
plots_nm$boxplot
plots_nm$density
plots_nm$rle
plots_nm$pca
plots_nm$pcv
plots_nm$pmad
plots_nm$pev
plots_nm$correlation
plots_nm$mds
plots_nm$dendrogram
plots_nm$ma
plots_nm$meansd
plots_nm$cv_intensity


nm_plot_density(se_nm, assay_names = "cycloess")

nm_plot_rle(se_nm, assay_names = c("cycloess", "Quantile"))