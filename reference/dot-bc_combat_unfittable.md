# Identify features ComBat cannot fit

Flags rows that, in any batch, have fewer than 2 finite observations or
zero within-batch variance. ComBat's parametric empirical-Bayes
estimation divides by the within-batch variance, so such features yield
NaN and abort the whole correction (`while (change > conv)` receives
NA). These features cannot be batch-corrected anyway (constant within a
batch) and should pass through unadjusted.

## Usage

``` r
.bc_combat_unfittable(mat, batch_vec)
```

## Arguments

- mat:

  Numeric matrix (features x samples).

- batch_vec:

  Batch assignment per column.

## Value

Logical vector (length nrow(mat)); TRUE = unfittable by ComBat.
