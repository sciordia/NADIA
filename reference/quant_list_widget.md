# Full Widget for the Protein_QUANT Table with the log2 Matrix Appended

Reactable table that reproduces the layout of the Protein-List_QUANT
Excel sheet, with 2-level headers, sticky columns, a blue palette,
aggregated filters and Excel export. It appends to Protein_QUANT the 16
columns of the normalized/imputed log2 matrix, matching by Protein
Groups.

## Usage

``` r
quant_list_widget(
  data = NULL,
  matrix_data = NULL,
  metadata = NULL,
  page_size = 15,
  height = 720,
  element_id = "quant_id_table",
  selection = NULL,
  searchable = TRUE
)
```

## Arguments

- data:

  Data frame or path to a Protein_QUANT TSV/CSV/Parquet file. Required;
  there is no default.

- matrix_data:

  Data frame or path to a TSV/CSV/Parquet file with the
  normalized/imputed log2 matrix (a ProteinGroups column + \_ samples).

- metadata:

  Optional: metadata data frame with a Coding column that fixes the
  sample order.

- page_size:

  Rows per page (default 15).

- height:

  Table height in px (default 720).

- element_id:

  Reactable element ID.

- selection:

  "single" \| "multiple" \| NULL.

- searchable:

  Enable the search box (default TRUE).

## Value

A browsable object.

## Examples

``` r
if (requireNamespace("jsonlite", quietly = TRUE)) {
  data(nadia_dia)
  res <- process_proteomics(nadia_dia, verbose = FALSE)

  # The log2 matrix: a ProteinGroups column plus one column per sample
  mat <- SummarizedExperiment::assay(res$se_proc, "Impseqrob_min")
  mat_df <- data.frame(ProteinGroups = rownames(mat), mat,
                       check.names = FALSE)

  w <- quant_list_widget(nadia_dia$protein_quant, matrix_data = mat_df,
                         metadata = nadia_dia$metadata)
  print(class(w))
}
#> [1] "shiny.tag.list" "list"          
```
