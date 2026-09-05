# Filter features by significance

Filter features by significance

## Usage

``` r
filter_significant_features(
  feature_info,
  filter_mode = c("any", "all", "specific"),
  comparison = NULL,
  alpha = 0.05,
  verbose = TRUE
)
```

## Arguments

- feature_info:

  DataFrame with significance columns

- filter_mode:

  Selection mode:

  - "any": features significant in AT LEAST one comparison (uses
    sig_any).

  - "all": ALL features, with no significance filtering (it does NOT
    mean "significant in every comparison").

  - "specific": significant in the comparison given by `comparison`.

- comparison:

  Specific comparison (for mode="specific")

- alpha:

  Significance threshold

- verbose:

  Print progress messages (default: TRUE)

## Value

Vector of the selected FeatureIDs
