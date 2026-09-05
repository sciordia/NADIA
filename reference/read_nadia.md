# Read a `.nadia` file

Loads the whole analysis: the preprocessing tables, the processing
results, the intensity matrices, the parameters and the provenance. The
canonical tables come out of the SQL views stored in the file, with
their R classes and column order restored.

## Usage

``` r
read_nadia(file)
```

## Arguments

- file:

  Path to a `.nadia` file.

## Value

An object of class `nadia_db`: a list with `file`, `meta`, `parameters`,
`calls`, `packages`, `source_files`, the canonical tables (`metadata`,
`protein_id`, `protein_quant` and, when present, `DEPs_results`,
`BoxPlot_Input`, `PCA_Input`, `pattern_profiler`), the `matrices` and
the `imputed_values` mask.

## Details

To get objects the rest of the package can use directly, pass the result
to
[`nadia_preprocessing()`](https://sciordia.github.io/NADIA/reference/nadia_preprocessing.md)
or
[`nadia_result()`](https://sciordia.github.io/NADIA/reference/nadia_result.md).

## See also

[`write_nadia()`](https://sciordia.github.io/NADIA/reference/write_nadia.md),
[`nadia_result()`](https://sciordia.github.io/NADIA/reference/nadia_result.md),
[`nadia_preprocessing()`](https://sciordia.github.io/NADIA/reference/nadia_preprocessing.md)

## Examples

``` r
if (requireNamespace("duckdb", quietly = TRUE)) {
    data(nadia_dia)
    res <- process_proteomics(nadia_dia, verbose = FALSE)
    f <- file.path(tempdir(), "read.nadia")
    write_nadia(f, nadia_dia, res, verbose = FALSE)

    db <- read_nadia(f)
    db
    identical(db$protein_quant, nadia_dia$protein_quant)

    unlink(f)
}
```
