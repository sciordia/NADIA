# Prepare Long-Format Data for tidyHeatmap

Prepare Long-Format Data for tidyHeatmap

## Usage

``` r
prepare_heatmap_data(
  data,
  mode = c("all", "any", "target"),
  alpha = 0.05,
  comparison = NULL,
  feature_ids = NULL,
  scale_data = c("row", "none", "column"),
  sample_order = c("clustering", "condition", "custom"),
  condition_order = NULL
)
```

## Arguments

- data:

  Data frame in long format with columns:

  - SampleID: sample identifier

  - FeatureID: protein/feature identifier

  - Intensity: intensity value (log2)

  - Condition: experimental condition

  - Replicate: replicate number

  - sig_any: logical flagging significance (for mode="any")

  - adjP\_\*: adjusted p-value columns (for mode="target")

- mode:

  Filtering mode: "all", "any", or "target"

- alpha:

  Significance threshold for mode "target" (default: 0.05)

- comparison:

  Name of the comparison for mode "target" (e.g. "B-A")

- feature_ids:

  Vector of specific FeatureIDs (when supplied, mode/alpha are ignored)

- scale_data:

  Scaling type: "none", "row", "column" (default: "row")

- sample_order:

  Sample order: "clustering", "condition", or a custom vector

- condition_order:

  Condition order when sample_order = "condition"

## Value

Data frame in long format, ready for tidyHeatmap
