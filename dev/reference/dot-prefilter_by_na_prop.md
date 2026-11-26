# Pre-filter proteins by NA proportion

Removes proteins with more than `max_na_prop` fraction of NAs. With
`max_na_prop = NULL` (the default) nothing is removed, so every
imputation method sees the same set of proteins: `combo` never applied
this filter, and having the single methods apply it made the same
argument mean different things depending on `imp_method`.

## Usage

``` r
.prefilter_by_na_prop(x, max_na_prop = NULL)
```

## Arguments

- x:

  Numeric matrix

- max_na_prop:

  Maximum NA proportion, or `NULL` to disable the filter (default).

## Value

List with keep (logical vector), summary
