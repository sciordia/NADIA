# Confusion Matrix Heatmap (ggplot2)

Heatmap of Comparisons x (TP%, FP%, FN%, TN%) with white-to-blue
gradient.

## Usage

``` r
benchmark_confusion_gg(
  confusion_df,
  title = "Confusion Matrix by Comparison x Species",
  text_size = 3.5,
  axis_text_size = 10
)
```

## Arguments

- confusion_df:

  Data frame from .confusion_by_species()

- title:

  Plot title

- text_size:

  Size of cell text labels

- axis_text_size:

  Size of axis text

## Value

ggplot2 object

## Examples

``` r
data(nadia_dia)
de <- process_proteomics(nadia_dia, verbose = FALSE)$DEPs_results
sp <- utils::read.delim(system.file("extdata", "nadia_dia_spikein.tsv.gz",
                                    package = "NADIA"))
de$Species <- sp$PG.OrganismId[match(de$Protein.IDs, sp$PG.ProteinGroups)]
ev <- data.frame(Comparison = rep(c("B-A", "D-A"), each = 2),
                 Species = c("ECOLI", "YEAST"),
                 expected_logFC = c(1, -0.58, 2, -3.3))
bench <- benchmarking_proteomics(de, ev, verbose = FALSE)
benchmark_confusion_gg(bench$confusion_by_species)
```
