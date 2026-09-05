# Rebuild the preprocessing object from a `.nadia` file

Returns an object indistinguishable from what the `preprocess_*()`
function that wrote the file produced – Spectronaut, DIA-NN, TMT or
label-free – so it can be fed straight to
[`process_proteomics()`](https://sciordia.github.io/NADIA/reference/process_proteomics.md)
or to the table functions.

## Usage

``` r
nadia_preprocessing(db)
```

## Arguments

- db:

  A `nadia_db` from
  [`read_nadia()`](https://sciordia.github.io/NADIA/reference/read_nadia.md).

## Value

A `proteomics_data` object with `metadata`, `protein_id` and
`protein_quant`.

## See also

[`read_nadia()`](https://sciordia.github.io/NADIA/reference/read_nadia.md),
[`nadia_result()`](https://sciordia.github.io/NADIA/reference/nadia_result.md)

## Examples

``` r
if (requireNamespace("duckdb", quietly = TRUE)) {
    data(nadia_dia)
    f <- file.path(tempdir(), "pre.nadia")
    write_nadia(f, nadia_dia, verbose = FALSE)

    pre <- nadia_preprocessing(read_nadia(f))
    identical(pre$protein_id, nadia_dia$protein_id)

    unlink(f)
}
```
