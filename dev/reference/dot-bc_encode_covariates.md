# Encode protected covariates as a numeric design for BERT

BERT performs no encoding of its own:
[`BERT::BERT()`](https://rdrr.io/pkg/BERT/man/BERT.html) collects the
`Cov_*` columns with
`mod <- data.frame(data[, grepl("Cov", names(data))])` and hands them
straight to `sva::ComBat(mod = )` and
`limma::removeBatchEffect(design = )`. Its documentation therefore
requires integer columns. A character or factor column reaches ComBat as
text, which is coerced to double and becomes NA ("NA/NaN/Inf in foreign
function call"), and reaches limma as a non-numeric design ("design must
be a numeric matrix").

## Usage

``` r
.bc_encode_covariates(df)
```

## Arguments

- df:

  data.frame of covariate columns (rows = samples).

## Value

Numeric data.frame with the same rows and one or more columns.

## Details

Numeric columns are passed through unchanged, so a genuinely continuous
covariate keeps its meaning. Categorical columns become k-1
treatment-contrast indicators, without an intercept: ComBat builds its
own batch design, which already spans the grand mean, and adding an
intercept here would make the combined design rank deficient. For a
two-level covariate this reproduces the single 0/1 column used in BERT's
own examples.
