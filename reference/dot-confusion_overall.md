# Compute overall confusion matrix per comparison

Aggregates TP/FP/TN/FN across all species for each comparison.

## Usage

``` r
.confusion_overall(confusion_by_species_df)
```

## Arguments

- confusion_by_species_df:

  Data frame from .confusion_by_species()

## Value

Data frame with one row per Comparison
