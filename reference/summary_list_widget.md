# Interactive Reactable Table for the Sample Metadata

Builds a reactable table that reproduces the layout of the metadata
Excel sheet (FileName, Condition, Replicate, Coding, \# Unique PSMs, \#
Unique Peptides, \# Protein Groups). Black header, rows tinted by
condition in light shades and Condition/Coding cells with a strong
colour chip (the same scheme used by protein_list_widget() /
quant_list_widget()).

## Usage

``` r
summary_list_widget(
  data = NULL,
  page_size = 16,
  height = 540,
  element_id = "summary_table",
  selection = NULL,
  searchable = FALSE
)
```

## Arguments

- data:

  Data frame or path to a metadata TSV/CSV file. Required; there is no
  default. It must have the columns R.FileName, R.Condition,
  R.Replicate, Coding, R.PrecursorsIdentified,
  R.StrippedSequencesIdentified, R.ProteinGroupsIdentified.

- page_size:

  Page size (default 16, i.e. all the rows).

- height:

  Height in px of the table container (default 540).

- element_id:

  Widget id in the DOM (default "summary_table").

- selection:

  Selection type ("multiple" or NULL).

- searchable:

  Enable the global search box (default FALSE).

## Value

An htmltools object (browsable tagList) ready to be rendered.

## Examples

``` r
if (requireNamespace("jsonlite", quietly = TRUE)) {
  data(nadia_dia)

  # LFQ/TMT metadata lack the identification counts; the widget then shows
  # only FileName/Condition/Replicate/Coding
  w <- summary_list_widget(nadia_dia$metadata)
  print(class(w))
}
#> [1] "shiny.tag.list" "list"          
```
