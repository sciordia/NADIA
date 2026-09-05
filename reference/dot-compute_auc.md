# Compute AUC using pROC (optional)

Compute AUC using pROC (optional)

## Usage

``` r
.compute_auc(classified_df, p_col)
```

## Arguments

- classified_df:

  Classified data frame for one comparison

- p_col:

  P-value column name

## Value

AUC value or NA if pROC not available
