# Open a connection to a `.nadia` file

For querying the file directly with SQL or dplyr. Most users want
[`read_nadia()`](https://sciordia.github.io/NADIA/reference/read_nadia.md)
instead; this is the way in when the file is larger than memory or when
a specific query beats loading everything.

## Usage

``` r
nadia_connect(file, read_only = TRUE)
```

## Arguments

- file:

  Path to a `.nadia` file.

- read_only:

  Open without write access. `TRUE` by default, which is also what
  allows several processes to read the same file at once.

## Value

A DBI connection.

## Details

The caller owns the connection and must close it with
`DBI::dbDisconnect(con, shutdown = TRUE)`.

## See also

[`read_nadia()`](https://sciordia.github.io/NADIA/reference/read_nadia.md),
[`nadia_tables()`](https://sciordia.github.io/NADIA/reference/nadia_tables.md)

## Examples

``` r
if (requireNamespace("duckdb", quietly = TRUE)) {
    data(nadia_dia)
    f <- file.path(tempdir(), "query.nadia")
    write_nadia(f, nadia_dia,
                process_proteomics(nadia_dia, verbose = FALSE),
                verbose = FALSE)

    con <- nadia_connect(f)
    DBI::dbGetQuery(con, "SELECT COUNT(*) AS n FROM proteins")
    DBI::dbDisconnect(con, shutdown = TRUE)

    unlink(f)
}
```
