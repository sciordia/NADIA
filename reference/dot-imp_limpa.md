# limpa: probabilistic imputation with detection probability curve (limpa)

Uses dpc/dpcCN to estimate a detection probability curve, then
dpcQuantByRow to impute protein-level abundances with standard errors.
The returned EList is attached as an attribute for downstream dpcDE.

## Usage

``` r
.imp_limpa(x, args = list())
```

## Arguments

- args:

  list with optional `use_dpcCN` (default FALSE), `model` (`"on"` or
  `"cn"`; defaults to `"cn"` when `use_dpcCN` is TRUE and `"on"`
  otherwise), `dpc.slope` (default 0.8), `dpc.start` (default NULL),
  `iterations` (default 2), `subset` (default 2000), `robust` (default
  TRUE), `chunk` (default 1000), `verbose` (default FALSE)
