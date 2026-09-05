# Interactive Reactable Table for Differential Expression Results

Builds a professionally formatted reactable table to display protein
differential expression results. Works both in Shiny (it returns a
reactable object) and standalone via
[`results_list_widget()`](https://sciordia.github.io/NADIA/reference/results_list_widget.md).

## Usage

``` r
results_list_reactable(
  data,
  protein_quant = NULL,
  comparisons = NULL,
  ain = NULL,
  alpha = 0.05,
  lfc_thr = 0,
  page_size = 15,
  height = 720,
  show_missing = TRUE,
  element_id = NULL,
  selection = NULL,
  searchable = TRUE
)
```

## Arguments

- data:

  Data frame or path to a TSV/CSV/Parquet file with DE results. Required
  columns: Protein.IDs, Gene.Names, logFC, P.Value, adj.P.Val, Change,
  Comparison. Optional: Assay, and the four missingness columns
  MissGlobal, MissComp, MissCND1 and MissCND2, which are shown only if
  all four are present

- protein_quant:

  Data frame (preprocessing\$protein_quant) or path to a
  Protein_QUANT\_\*.tsv file. If not NULL, the Description and
  Quant_Pepts columns are added.

- comparisons:

  Vector of comparisons to include (NULL = all of them)

- ain:

  Vector of assays to keep (NULL = all of them)

- alpha:

  Significance threshold used to highlight the FDR (default: 0.05)

- lfc_thr:

  log2 fold-change threshold (default: 0, reserved for future use)

- page_size:

  Rows per page (default: 15)

- height:

  Table height in pixels (default: 720)

- show_missing:

  Show the missingness percentage columns (default: TRUE)

- element_id:

  Element ID for the Reactable JS API (default: NULL)

- selection:

  Selection type: "single", "multiple", or NULL (default: NULL)

- searchable:

  Enable the internal search box (default: TRUE)

## Value

A reactable object

## Examples

``` r
data(nadia_dia)
res <- process_proteomics(nadia_dia, verbose = FALSE)

# protein_quant adds the Description and Quant_Pepts columns
tbl <- results_list_reactable(
  res$DEPs_results,
  protein_quant = nadia_dia$protein_quant,
  comparisons   = "B-A"
)
class(tbl)
#> [1] "reactable"  "htmlwidget"
```
