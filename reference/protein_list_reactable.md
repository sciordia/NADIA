# Interactive Reactable Table for Protein_ID Data (post-Spectronaut)

Builds a reactable table with 2-level grouped headers, static columns
pinned (sticky) to the left and data bars coloured by condition inside
each numeric cell. It reproduces the usual Excel layout and improves on
it.

## Usage

``` r
protein_list_reactable(
  data,
  metadata = NULL,
  page_size = 15,
  height = 720,
  element_id = NULL,
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

  Reactable element ID (default NULL).

- selection:

  "single" \| "multiple" \| NULL.

- searchable:

  Enable the search box (default TRUE).

## Value

A reactable object.

## Examples

``` r
data(nadia_dia)

# metadata$Coding fixes the order of the sample columns
tbl <- protein_list_reactable(nadia_dia$protein_id,
                              metadata = nadia_dia$metadata)
class(tbl)
#> [1] "reactable"  "htmlwidget"
```
