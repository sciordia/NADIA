# NADIA: Missing Value-Aware DIA Proteomics Analysis

A complete pipeline for differential protein expression analysis, with
particular attention to the missing values that characterise
data-independent acquisition. It covers the whole journey from the raw
report to the interactive figure.

## Details

The package is organised in four blocks:

- Preprocessing:

  [`preprocess_spectronaut()`](https://sciordia.github.io/NADIA/reference/preprocess_spectronaut.md)
  for Spectronaut,
  [`preprocess_diann()`](https://sciordia.github.io/NADIA/reference/preprocess_diann.md)
  for DIA-NN protein-group matrices,
  [`preprocess_tmt()`](https://sciordia.github.io/NADIA/reference/preprocess_tmt.md)
  and
  [`preprocess_lfq()`](https://sciordia.github.io/NADIA/reference/preprocess_lfq.md)
  for Proteome Discoverer. All four return the same `proteomics_data` S3
  object, so the rest of the pipeline consumes them unchanged.

- Processing:

  [`process_proteomics()`](https://sciordia.github.io/NADIA/reference/process_proteomics.md)
  coordinates normalization (13 methods), optional batch correction,
  imputation (20 methods, including the MAR/MNAR hybrids) and
  differential expression with limma or limpa.

- Metrics and benchmarking:

  [`normalization_metrics()`](https://sciordia.github.io/NADIA/reference/normalization_metrics.md),
  [`imputation_metrics()`](https://sciordia.github.io/NADIA/reference/imputation_metrics.md),
  [`benchmarking_proteomics()`](https://sciordia.github.io/NADIA/reference/benchmarking_proteomics.md)
  and
  [`benchmarking_multiple()`](https://sciordia.github.io/NADIA/reference/benchmarking_multiple.md)
  evaluate and rank method combinations.

- Visualization:

  interactive plots with Highcharts
  ([`boxplot_highchart_list()`](https://sciordia.github.io/NADIA/reference/boxplot_highchart_list.md),
  [`volcano_highchart_list()`](https://sciordia.github.io/NADIA/reference/volcano_highchart_list.md),
  [`pca_highchart_list()`](https://sciordia.github.io/NADIA/reference/pca_highchart_list.md)),
  static ones with ggplot2 and ComplexHeatmap
  ([`proteomics_heatmap()`](https://sciordia.github.io/NADIA/reference/proteomics_heatmap.md)),
  and tables with reactable
  ([`results_list_reactable()`](https://sciordia.github.io/NADIA/reference/results_list_reactable.md)).

## See also

Useful links:

- <https://github.com/sciordia/NADIA>

- <https://sciordia.github.io/NADIA/>

- Report bugs at <https://github.com/sciordia/NADIA/issues>

## Author

**Maintainer**: Sergio Ciordia <sciordia@gmail.com>
([ORCID](https://orcid.org/0000-0002-5726-853X))

Authors:

- Sergio Ciordia <sciordia@gmail.com>
  ([ORCID](https://orcid.org/0000-0002-5726-853X))
