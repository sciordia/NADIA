# Compute partial AUC with McClish correction

Uses pROC::roc() with partial.auc and partial.auc.correct = TRUE (native
McClish correction). Normalized range: 0.5-1.

## Usage

``` r
.compute_pauc(classified_df, p_col, max_fpr = 0.1)
```

## Arguments

- classified_df:

  Classified data frame for one comparison

- p_col:

  P-value column name

- max_fpr:

  Maximum false positive rate threshold (default: 0.1)

## Value

Corrected partial AUC value or NA if pROC not available
