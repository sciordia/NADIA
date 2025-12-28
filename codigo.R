
# =========================================================
# Código de ejemplo: Volcano-Plot
# =========================================================

# DEPs_results <- read_tsv("DEPs_results.tsv")
DEPs_results <- arrow::read_parquet("./data/VolcanoPlot_Input.parquet")

# Cargar las funciones para el Volcano-Plot
source("./R/Volcano_Plot_Highcharts_Final.R")

# --- Ejemplo básico ---
hc_volcanos <- volcano_highchart_list(
  de_res      = DEPs_results,
  ain         = "LoessCyc",
  comparisons = c("B-A"),
  alpha       = 0.05,
  p_col       = "adj.P.Val",
  point_size  = 3
)

hc_volcanos[["B-A"]]

# --- Combinando todas las opciones ---
hc_volcanos <- volcano_highchart_list(
  de_res          = DEPs_results,
  ain             = "LoessCyc",
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

se_proc <- arrow::read_parquet("./data/BoxPlot_Input.parquet")

# Cargar las funciones para el Box-Plot
source("./R/Boxplot_Highcharts_Final.R")

hc_boxplots <- boxplot_highchart_list(
  data        = se_proc,
  assays      = c("log2", "LoessCyc"),
  color_by    = "Condition",
  palette  = "ggsci::category10_d3",
  group_order = c("A", "B", "C", "D"),
  box_width = 20
)

hc_boxplots[["log2"]]
hc_boxplots[["LoessCyc"]]


