# Introduce artificial NAs into a complete-case matrix

Takes a matrix with no NAs (ground truth) and introduces NAs at random
or mimicking the NA pattern from a reference matrix.

## Usage

``` r
.im_introduce_na(
  mat,
  na_prop = 0.2,
  seed = 42L,
  pattern = "random",
  ref_mat = NULL
)
```

## Arguments

- mat:

  Numeric matrix (proteins x samples), no NAs (ground truth)

- na_prop:

  Numeric (0-1). Proportion of values to set as NA. Default 0.20. Only
  used when `pattern = "random"`; ignored when `pattern = "from_data"`.

- seed:

  Integer. Random seed. Default 42.

- pattern:

  Character: `"random"` for uniform random NAs, or `"from_data"` to
  replicate NAguideR's simulation strategy using the actual NA structure
  from `ref_mat`.

- ref_mat:

  Numeric matrix. Reference matrix with real NAs, used when
  `pattern = "from_data"`. Ignored otherwise.

## Value

Named list:

- mat_with_na:

  Matrix with artificial NAs introduced

- na_mask:

  Logical matrix (TRUE = artificially set to NA)

- true_mat:

  Original complete matrix (unchanged)

- summary:

  List with n_total, n_na, na_rate, rows_affected, row_ratio (from_data
  only)

## Details

When `pattern = "from_data"`, replicates the NAguideR strategy (Wang et
al., DOI:10.1093/nar/gkz903): the proportion of affected rows and the
per-column NA distribution are derived entirely from `ref_mat`, so
`na_prop` is ignored.
