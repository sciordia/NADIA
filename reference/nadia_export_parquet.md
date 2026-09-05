# Export the tables of a `.nadia` file to Parquet

A way out to an open, stable format. DuckDB reads its own files back a
long way, but Parquet is the safer bet for something that has to be
readable in ten years by software nobody has written yet, and it is what
a repository will expect alongside a publication.

## Usage

``` r
nadia_export_parquet(db, dir, views = TRUE, verbose = TRUE)
```

## Arguments

- db:

  A `nadia_db` from
  [`read_nadia()`](https://sciordia.github.io/NADIA/reference/read_nadia.md),
  or the path to a file.

- dir:

  Directory to write to; created if it does not exist.

- views:

  Export the reconstructed canonical tables (`TRUE`, the default) or the
  underlying storage tables (`FALSE`).

- verbose:

  Print progress.

## Value

The paths written, invisibly.

## See also

[`read_nadia()`](https://sciordia.github.io/NADIA/reference/read_nadia.md),
[`write_nadia()`](https://sciordia.github.io/NADIA/reference/write_nadia.md)

## Examples

``` r
if (requireNamespace("duckdb", quietly = TRUE) &&
    requireNamespace("arrow", quietly = TRUE)) {
    data(nadia_dia)
    f <- file.path(tempdir(), "export.nadia")
    write_nadia(f, nadia_dia, verbose = FALSE)

    out <- file.path(tempdir(), "nadia_parquet_export")
    nadia_export_parquet(f, out, verbose = FALSE)
    list.files(out)

    unlink(c(f, out), recursive = TRUE)
}
```
