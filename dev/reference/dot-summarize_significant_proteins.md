# Summarize significant proteins per comparison and species

Generates a summary table with counts of significant and non-significant
proteins broken down by species and direction (up/down), plus derived
percentages and ratios.

## Usage

``` r
.summarize_significant_proteins(classified_df)
```

## Arguments

- classified_df:

  Data frame from .classify_all_comparisons() with columns: Comparison,
  Species, logFC, truth, predicted, classification, is_significant and
  direction_error

## Value

Data frame with per-comparison summary

## Details

Significance is read from the `is_significant` column rather than
recomputed, so this table can never disagree with the classification it
summarizes. It used to apply only `p <= alpha`, ignoring `lfc_thr`
entirely, which made the counts here larger than the ones behind the
metrics whenever a fold-change threshold was in force.
