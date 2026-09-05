# List the contents of a `.nadia` file without loading it

List the contents of a `.nadia` file without loading it

## Usage

``` r
nadia_tables(file)
```

## Arguments

- file:

  Path to a `.nadia` file.

## Value

A data frame with `name`, `type` (`table` or `view`) and `rows`.

## See also

[`read_nadia()`](https://sciordia.github.io/NADIA/reference/read_nadia.md)

## Examples

``` r
if (requireNamespace("duckdb", quietly = TRUE)) {
    data(nadia_dia)
    f <- file.path(tempdir(), "tables.nadia")
    write_nadia(f, nadia_dia, verbose = FALSE)

    nadia_tables(f)

    unlink(f)
}
```
