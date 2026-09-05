# Full Widget for the Protein_ID Table with Filters, Search and Excel Export

Wraps
[`protein_list_reactable()`](https://sciordia.github.io/NADIA/reference/protein_list_reactable.md)
with:

- A search bar and buttons (toggle filters, clear, export to Excel).

- Per-condition chips that hide/show the 16 associated columns.

- A numeric filter on MW in Da with AND/OR operators.

- Excel export with 2-level headers (merged groups + sub-labels) and a
  per-condition tint, reproducing the manual layout.

## Usage

``` r
protein_list_widget(
  data,
  metadata = NULL,
  page_size = 15,
  height = 720,
  element_id = "protein_id_table",
  selection = NULL,
  searchable = TRUE
)
```

## Arguments

- data:

  Data frame or path to a Protein_ID TSV/CSV/Parquet file (output of any
  preprocess\_\*() function -\> protein_id).

- metadata:

  Optional: metadata data frame (run_summary) with a Coding column that
  fixes the order of the sample columns.

- page_size:

  Rows per page (default 15).

- height:

  Table height in px (default 720).

- element_id:

  Element ID (default "protein_id_table").

- selection:

  "single" \| "multiple" \| NULL.

- searchable:

  Enable the search box (default TRUE).

## Value

A browsable htmltools object.

## Examples

``` r
if (requireNamespace("jsonlite", quietly = TRUE)) {
  data(nadia_dia)

  w <- protein_list_widget(nadia_dia$protein_id,
                           metadata = nadia_dia$metadata)
  print(class(w))   # print(w) itself opens it in the Viewer / browser
}
#> [1] "shiny.tag.list" "list"          
```
