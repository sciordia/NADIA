# Compute AUC using pROC

Compute AUC using pROC

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

AUC value, or NA if the comparison lacks both classes or has fewer than
three proteins
