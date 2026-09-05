# Rebuild the Pattern Profiler table from a `.nadia` file

Rebuild the Pattern Profiler table from a `.nadia` file

## Usage

``` r
nadia_pattern_profiler(db)
```

## Arguments

- db:

  A `nadia_db` from
  [`read_nadia()`](https://sciordia.github.io/NADIA/reference/read_nadia.md).

## Value

The `long_output` data frame, ready for
[`cluster_profile_highchart_list()`](https://sciordia.github.io/NADIA/reference/cluster_profile_highchart_list.md),
or `NULL` if the file holds no clustering.

## See also

[`read_nadia()`](https://sciordia.github.io/NADIA/reference/read_nadia.md),
[`nadia_add_pattern_profiler()`](https://sciordia.github.io/NADIA/reference/nadia_add_pattern_profiler.md)

## Examples

``` r
if (requireNamespace("duckdb", quietly = TRUE)) {
    data(nadia_dia)
    f <- file.path(tempdir(), "nopp.nadia")
    write_nadia(f, nadia_dia, verbose = FALSE)

    is.null(nadia_pattern_profiler(read_nadia(f)))

    unlink(f)
}
```
