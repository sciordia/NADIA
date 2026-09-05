# PVCA Analysis – Orchestrator

Runs the full PVCA workflow: compute variance components, generate plot,
and optionally export results.

## Usage

``` r
pvca_analysis(
  se,
  assay_name = NULL,
  technical_factors = character(0),
  biological_factors = character(0),
  pca_covariates = NULL,
  pca_threshold = 0.6,
  variance_threshold = 0.01,
  na_action = "complete",
  fill_value = -1,
  colors = NULL,
  verbose = TRUE,
  output_dir = NULL,
  export_plots = TRUE,
  export_tables = TRUE,
  plot_width = 10,
  plot_height = 6,
  plot_dpi = 150
)
```

## Arguments

- se:

  SummarizedExperiment object.

- assay_name:

  Character. Assay to analyze. NULL = second assay (normalized,
  pre-imputation).

- technical_factors:

  Character vector. Factors representing technical variation (e.g.,
  "Injection", "Digestion").

- biological_factors:

  Character vector. Factors representing biological variation (e.g.,
  "Condition", "Gender", "Age").

- pca_covariates:

  Character vector or NULL. If provided, generates PCA plots colored by
  each covariate (grid + individual plots). These are included in the
  returned list and exported as PNGs. Default NULL (skip).

- pca_threshold:

  Numeric (0-1). Cumulative variance threshold for PCA (default 0.6).

- variance_threshold:

  Numeric (0-1). Minimum weight for individual display (default 0.01).

- na_action:

  Character: "complete", "fill", or "none" (default "complete").

- fill_value:

  Numeric. NA replacement when na_action = "fill" (default -1).

- colors:

  Named character vector for category colors (default NULL = built-in
  palette).

- verbose:

  Logical (default TRUE).

- output_dir:

  Character. Directory for exports. NULL = no export.

- export_plots:

  Logical (default TRUE).

- export_tables:

  Logical (default TRUE).

- plot_width:

  Numeric in inches (default 10).

- plot_height:

  Numeric in inches (default 6).

- plot_dpi:

  Numeric (default 150).

## Value

Named list:

- variance_components:

  data.frame with label, weights, category

- plot:

  ggplot2 object

- n_proteins:

  Number of proteins used in the analysis

- assay_name:

  Assay analyzed

- parameters:

  List of parameters used

## Examples

``` r
data(nadia_dia)

# nadia_dia has no batch information, so the example builds a plausible one:
# two digestion batches crossed with the three conditions.
cov <- data.frame(Column = nadia_dia$metadata$Coding,
                  Batch = rep(c("b1", "b2"),
                              length.out = nrow(nadia_dia$metadata)),
                  stringsAsFactors = FALSE)
res <- process_proteomics(nadia_dia, covariate_df = cov, verbose = FALSE)

# Variance decomposition: how much of the signal each factor explains
pv <- pvca_analysis(res$se_proc, assay_name = "Impseqrob_min",
                    technical_factors = "Batch",
                    biological_factors = "Condition",
                    verbose = FALSE)
#> boundary (singular) fit: see help('isSingular')
pv$pvca_table
#> NULL
```
