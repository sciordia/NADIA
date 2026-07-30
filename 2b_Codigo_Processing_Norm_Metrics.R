
# =========================================================
# Código de ejemplo: Normalization Metrics + Plots
# =========================================================

# Cargar las funciones para las métricas de normalización
library(NADIA)

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
plots_nm$scatter
plots_nm$pca
plots_nm$pcv
plots_nm$pmad
plots_nm$pev
plots_nm$correlation
plots_nm$mds
plots_nm$qq
plots_nm$metrics


nm_plot_pev(se_nm, diff = TRUE, baseline = "log2Norm")
nm_plot_pmad(se_nm, diff = TRUE, baseline = "log2Norm")
nm_plot_pcv(se_nm, diff = TRUE, baseline = "log2Norm")

nm_plot_scatter(se_nm, assay_names = c("log2Norm",
                                   "eqmedians",
                                   "quantile.robust",
                                   "cycloess",
                                   "quantile",
                                   "vsn"))

nm_plot_pca(se_nm, assay_names = c("log2Norm",
                                   "Rlr",
                                   "GlobalMean",
                                   "cycloess",
                                   "GlobalMedian",
                                   "vsn"))

nm_plot_boxplot(se_nm, assay_names = c("log2Norm",
                                   "eqmedians",
                                   "quantile.robust",
                                   "cycloess",
                                   "quantile",
                                   "vsn"))

nm_plot_density(se_nm, assay_names = c("log2Norm",
                                       "eqmedians",
                                       "quantile.robust",
                                       "cycloess",
                                       "quantile",
                                       "vsn"))

nm_plot_pmad(se_nm, assay_names = c("log2Norm",
                                       "eqmedians",
                                       "quantile.robust",
                                       "cycloess",
                                       "quantile",
                                       "vsn"))

nm_plot_pev(se_nm, assay_names = c("log2Norm",
                                    "eqmedians",
                                    "quantile.robust",
                                    "cycloess",
                                    "quantile",
                                    "vsn"))
nm_plot_correlation(se_nm, assay_names = c("log2Norm",
                                   "eqmedians",
                                   "quantile.robust",
                                   "cycloess",
                                   "quantile",
                                   "vsn"))
nm_plot_qq(se_nm, assay_names = c("log2Norm",
                                   "eqmedians",
                                   "quantile.robust",
                                   "cycloess",
                                   "quantile",
                                   "vsn"))
nm_plot_pcv(se_nm, assay_names = c("log2Norm",
                                            "eqmedians",
                                            "quantile.robust",
                                            "cycloess",
                                            "quantile",
                                            "vsn"))
nm_plot_mds(se_nm, assay_names = c("log2Norm",
                                            "eqmedians",
                                            "quantile.robust",
                                            "cycloess",
                                            "quantile",
                                            "vsn"))
