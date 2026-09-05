# Benchmark Proteomics Differential Expression Results

Evaluates the reliability of DE results by comparing them against known
ground truth (multi-species spike-in experiments). Computes
classification metrics, dispersion statistics, and generates
visualizations.

## Usage

``` r
benchmarking_proteomics(
  de_res,
  expected_values,
  alpha = 0.05,
  lfc_thr = 0,
  p_col = "adj.P.Val",
  comparisons = NULL,
  assay = NULL,
  output_dir = NULL,
  verbose = TRUE,
  species_colors = NULL,
  species_df = NULL
)
```

## Arguments

- de_res:

  Data frame with DE results. Required columns: Protein.IDs, logFC,
  P.Value, adj.P.Val, Change, Comparison. Optional: Gene.Names, Assay,
  Species (if not using species_df).

- expected_values:

  Data frame with expected values. Required columns: Comparison,
  Species, expected_logFC.

- alpha:

  Significance threshold (default: 0.05)

- lfc_thr:

  Log fold-change threshold (default: 0)

- p_col:

  P-value column to use (default: "adj.P.Val")

- comparisons:

  Comparisons to include (NULL = all)

- assay:

  Assay to filter (NULL = all)

- output_dir:

  Directory for exported files (NULL = no export)

- verbose:

  Print progress messages (default: TRUE)

- species_colors:

  Named vector of colors per species (optional)

- species_df:

  Data frame with Protein.IDs and Species columns (optional). If NULL,
  de_res must already contain a Species column.

## Value

List with:

- metrics_table: Classification metrics per comparison

- confusion_by_species: Confusion matrix per Comparison x Species

- confusion_overall: Confusion matrix aggregated per Comparison

- dispersion_metrics: Dispersion stats per Comparison x Species

- classified_df: Full classified data frame

- gg_heatmap: ggplot2 performance heatmap

- gg_confusion_by_species: ggplot2 confusion matrix heatmap (by species)

- gg_confusion_overall: ggplot2 confusion matrix heatmap (aggregated)

- gg_auc_bars: ggplot2 AUC bar chart

- gg_metrics_bars: ggplot2 grouped metrics bar chart

- gg_signif_bars: ggplot2 significant proteins stacked bars (absolute)

- hc_volcano_list: Named list of Highcharter volcano plots

- parameters: List of parameters used

## Examples

``` r
data(nadia_dia)
de <- process_proteomics(nadia_dia, verbose = FALSE)$DEPs_results
sp <- utils::read.delim(system.file("extdata", "nadia_dia_spikein.tsv.gz",
                                    package = "NADIA"))
species_df <- data.frame(Protein.IDs = sp$PG.ProteinGroups,
                         Species = sp$PG.OrganismId)
ev <- data.frame(Comparison = rep(c("B-A", "D-A"), each = 2),
                 Species = c("ECOLI", "YEAST"),
                 expected_logFC = c(1, -0.58, 2, -3.3))
bench <- benchmarking_proteomics(de, ev, species_df = species_df,
                                 verbose = FALSE)
bench$metrics_table
#>   Comparison  TP FP   TN  FN Sensitivity Specificity Precision    NPV Accuracy
#> 1        B-A 469 27 1227 274      0.6312      0.9785    0.9456 0.8175   0.8493
#> 2        D-A 693 60 1194  50      0.9327      0.9522    0.9203 0.9598   0.9449
#>       F1    MCC    AUC Performance
#> 1 0.7571 0.6821 0.8323        Good
#> 2 0.9265 0.8825 0.9779   Excellent
bench$gg_heatmap
```
