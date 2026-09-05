# Write a complete analysis to a single `.nadia` file

Stores the preprocessing tables, the processing results and, optionally,
the Pattern Profiler output in one DuckDB database, together with the
parameters and the provenance of the run. Nothing is duplicated:
everything that can be derived is left out and rebuilt on reading,
through SQL views stored in the file itself.

## Usage

``` r
write_nadia(
  file,
  preprocessing,
  result = NULL,
  pattern_profiler = NULL,
  overwrite = FALSE,
  verbose = TRUE
)
```

## Arguments

- file:

  Path to write to. The `.nadia` extension is a convention; the file is
  a DuckDB database whatever it is called.

- preprocessing:

  A `proteomics_data` object, from any of the `preprocess_*()`
  functions. Required: it is the base of the file.

- result:

  Object from
  [`process_proteomics()`](https://sciordia.github.io/NADIA/reference/process_proteomics.md),
  or `NULL` to store only the preprocessing.

- pattern_profiler:

  Object from
  [`pattern_profiler_analysis()`](https://sciordia.github.io/NADIA/reference/pattern_profiler_analysis.md),
  or `NULL`. Can also be added later with
  [`nadia_add_pattern_profiler()`](https://sciordia.github.io/NADIA/reference/nadia_add_pattern_profiler.md).

- overwrite:

  Replace `file` if it already exists.

- verbose:

  Print progress.

## Value

The path, invisibly.

## Details

Use
[`read_nadia()`](https://sciordia.github.io/NADIA/reference/read_nadia.md)
to read it back, and
[`nadia_result()`](https://sciordia.github.io/NADIA/reference/nadia_result.md)
to get an object the rest of the package can plot directly.

## What goes in

Twelve tables. `samples` and `proteins` are the dimensions;
`protein_sample_metrics`, `intensities`, `imputed_values` and
`de_results` the facts; `pp_membership` and `pp_profile` the optional
clustering; and `nadia_meta`, `nadia_packages`, `nadia_source_files`,
`nadia_calls`, `parameters` and `nadia_schema` the provenance.

Three things are deliberately *not* stored, because they are exactly
derivable: the `log2` assay (`log2` of the raw one), the `PCA_Input`
table (the imputed assay joined to the DE results) and the observed
cells of the imputed assay (identical to the normalized ones). The last
of these is checked rather than assumed, since a model-based method such
as `limpa` can revise observed values; when it does, the full assay is
stored instead.

## Where not to put it

Not inside a synchronised folder such as iCloud Drive, Dropbox or
OneDrive while you are working. DuckDB keeps the file open with locks
and partial writes, and a synchronising client can corrupt it. Write it
locally and copy it afterwards.

## See also

[`read_nadia()`](https://sciordia.github.io/NADIA/reference/read_nadia.md),
[`nadia_result()`](https://sciordia.github.io/NADIA/reference/nadia_result.md),
[`nadia_add_pattern_profiler()`](https://sciordia.github.io/NADIA/reference/nadia_add_pattern_profiler.md)

## Examples

``` r
if (requireNamespace("duckdb", quietly = TRUE)) {
    data(nadia_dia)
    res <- process_proteomics(nadia_dia, verbose = FALSE)

    f <- file.path(tempdir(), "example.nadia")
    write_nadia(f, nadia_dia, res, verbose = FALSE)

    db <- read_nadia(f)
    db

    unlink(f)
}
```
