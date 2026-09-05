# Package index

## Import and preprocessing

Four readers, one output contract. Each converts a quantification report
into the same `proteomics_data` object, so everything downstream is
independent of the search engine.

- [`preprocess_spectronaut()`](https://sciordia.github.io/NADIA/reference/preprocess_spectronaut.md)
  : Preprocess Spectronaut reports
- [`preprocess_diann()`](https://sciordia.github.io/NADIA/reference/preprocess_diann.md)
  : Preprocess DIA-NN protein-group matrices
- [`preprocess_tmt()`](https://sciordia.github.io/NADIA/reference/preprocess_tmt.md)
  : Preprocess Proteome Discoverer TMT exports
- [`preprocess_lfq()`](https://sciordia.github.io/NADIA/reference/preprocess_lfq.md)
  : Preprocess Proteome Discoverer LFQ exports

## The analysis pipeline

[`process_proteomics()`](https://sciordia.github.io/NADIA/reference/process_proteomics.md)
runs normalisation, optional batch correction, imputation and
differential abundance in order. The individual stages are exported too,
for workflows that need them separately.

- [`process_proteomics()`](https://sciordia.github.io/NADIA/reference/process_proteomics.md)
  : Process proteomics data
- [`normalize_proteomics()`](https://sciordia.github.io/NADIA/reference/normalize_proteomics.md)
  : Normalize proteomics data
- [`batch_correct_proteomics()`](https://sciordia.github.io/NADIA/reference/batch_correct_proteomics.md)
  : Batch Correction with BERT
- [`impute_proteomics()`](https://sciordia.github.io/NADIA/reference/impute_proteomics.md)
  : Impute proteomics data
- [`de_analysis_proteomics()`](https://sciordia.github.io/NADIA/reference/de_analysis_proteomics.md)
  : Perform differential expression analysis on proteomics data

## Batch effects

Deciding whether a batch effect matters before correcting it: variance
decomposition and PCA against the covariates.

- [`pvca_analysis()`](https://sciordia.github.io/NADIA/reference/pvca_analysis.md)
  : PVCA Analysis – Orchestrator
- [`pvca_compute()`](https://sciordia.github.io/NADIA/reference/pvca_compute.md)
  : Compute PVCA variance components
- [`pvca_plot()`](https://sciordia.github.io/NADIA/reference/pvca_plot.md)
  : Plot PVCA variance components
- [`pca_covariates_plot()`](https://sciordia.github.io/NADIA/reference/pca_covariates_plot.md)
  : PCA Plots Colored by Covariates

## Choosing a normalisation method

Ranks methods on properties of the data themselves, with no ground truth
required. Pass `methods` to compare the twelve methods; without it the
assays already in the object are ranked instead.

- [`normalization_metrics()`](https://sciordia.github.io/NADIA/reference/normalization_metrics.md)
  : Generate quality metric plots for normalization comparison
- [`nm_compute_metrics()`](https://sciordia.github.io/NADIA/reference/nm_compute_metrics.md)
  : Compute quantitative group-separation metrics per normalization
  method
- [`nm_prepare_se()`](https://sciordia.github.io/NADIA/reference/nm_prepare_se.md)
  : Prepare a SummarizedExperiment from a proteomics_data object
- [`nm_run_normalizations()`](https://sciordia.github.io/NADIA/reference/nm_run_normalizations.md)
  : Run multiple normalization methods from a baseline assay
- [`import_norm_matrices()`](https://sciordia.github.io/NADIA/reference/import_norm_matrices.md)
  : Import normalized matrices into a SummarizedExperiment
- [`nm_rank_final()`](https://sciordia.github.io/NADIA/reference/nm_rank_final.md)
  : Compute a combined final ranking across PCV, PMAD, PEV, Correlation
  and group-separation (PC1 F-ratio)
- [`nm_rank_pcv()`](https://sciordia.github.io/NADIA/reference/nm_rank_pcv.md)
  : Rank normalization methods by median PCV (ascending – lower is
  better)
- [`nm_rank_pmad()`](https://sciordia.github.io/NADIA/reference/nm_rank_pmad.md)
  : Rank normalization methods by median PMAD (ascending – lower is
  better)
- [`nm_rank_pev()`](https://sciordia.github.io/NADIA/reference/nm_rank_pev.md)
  : Rank normalization methods by median PEV (ascending – lower is
  better)
- [`nm_rank_cor()`](https://sciordia.github.io/NADIA/reference/nm_rank_cor.md)
  : Rank normalization methods by median intragroup correlation
  (descending – higher is better)
- [`nm_rank_pc1()`](https://sciordia.github.io/NADIA/reference/nm_rank_pc1.md)
  : Rank normalization methods by PC1 variance explained
- [`nm_rank_mds1()`](https://sciordia.github.io/NADIA/reference/nm_rank_mds1.md)
  : Rank normalization methods by MDS1 variance explained
- [`nm_plot_final_ranking()`](https://sciordia.github.io/NADIA/reference/nm_plot_final_ranking.md)
  : Horizontal bar chart of combined final ranking
- [`nm_plot_metrics()`](https://sciordia.github.io/NADIA/reference/nm_plot_metrics.md)
  : Bar chart of group-separation metrics per normalization method
- [`nm_plot_boxplot()`](https://sciordia.github.io/NADIA/reference/nm_plot_boxplot.md)
  : Intensity boxplot per method
- [`nm_plot_density()`](https://sciordia.github.io/NADIA/reference/nm_plot_density.md)
  : Density plot per method
- [`nm_plot_correlation()`](https://sciordia.github.io/NADIA/reference/nm_plot_correlation.md)
  : Intra-group correlation boxplot per method (PRONE-style)
- [`nm_plot_pcv()`](https://sciordia.github.io/NADIA/reference/nm_plot_pcv.md)
  : PCV boxplot per method (PRONE-style)
- [`nm_plot_pmad()`](https://sciordia.github.io/NADIA/reference/nm_plot_pmad.md)
  : PMAD boxplot per method (PRONE-style)
- [`nm_plot_pev()`](https://sciordia.github.io/NADIA/reference/nm_plot_pev.md)
  : PEV boxplot per method (PRONE-style)
- [`nm_plot_pca()`](https://sciordia.github.io/NADIA/reference/nm_plot_pca.md)
  : PCA scatter plot per method
- [`nm_plot_mds()`](https://sciordia.github.io/NADIA/reference/nm_plot_mds.md)
  : MDS 2D scatter per method
- [`nm_plot_pc1_ranking()`](https://sciordia.github.io/NADIA/reference/nm_plot_pc1_ranking.md)
  : Horizontal bar chart of PC1 variance ranking
- [`nm_plot_mds1_ranking()`](https://sciordia.github.io/NADIA/reference/nm_plot_mds1_ranking.md)
  : Horizontal bar chart of MDS1 variance ranking
- [`nm_plot_qq()`](https://sciordia.github.io/NADIA/reference/nm_plot_qq.md)
  : Q-Q plot per method (NormalyzerDE style)
- [`nm_plot_scatter()`](https://sciordia.github.io/NADIA/reference/nm_plot_scatter.md)
  : Sample-vs-sample scatter plot per method (NormalyzerDE style)

## Choosing an imputation method

Simulated missingness on complete rows. The masking hides observed
values, so the ranking measures MAR recovery: use it for the MAR stage,
not to choose an MNAR one.

- [`imputation_metrics()`](https://sciordia.github.io/NADIA/reference/imputation_metrics.md)
  : Generate imputation quality metric plots and ranking
- [`im_compute_metrics()`](https://sciordia.github.io/NADIA/reference/im_compute_metrics.md)
  : Compute imputation quality metrics (NRMSE, SOR, PSS, ACC_OI)
- [`im_prepare_se()`](https://sciordia.github.io/NADIA/reference/im_prepare_se.md)
  : Build a SummarizedExperiment for imputation benchmarking
- [`import_imp_matrices()`](https://sciordia.github.io/NADIA/reference/import_imp_matrices.md)
  : Import imputed matrices into a SummarizedExperiment
- [`im_plot_ranking()`](https://sciordia.github.io/NADIA/reference/im_plot_ranking.md)
  : Ranking heatmap of imputation methods
- [`im_plot_metrics()`](https://sciordia.github.io/NADIA/reference/im_plot_metrics.md)
  : Faceted bar chart of all 4 imputation quality metrics
- [`im_plot_nrmse()`](https://sciordia.github.io/NADIA/reference/im_plot_nrmse.md)
  : NRMSE bar plot per imputation method
- [`im_plot_sor()`](https://sciordia.github.io/NADIA/reference/im_plot_sor.md)
  : SOR bar plot per imputation method
- [`im_plot_pss()`](https://sciordia.github.io/NADIA/reference/im_plot_pss.md)
  : PSS bar plot per imputation method
- [`im_plot_acc_oi()`](https://sciordia.github.io/NADIA/reference/im_plot_acc_oi.md)
  : ACC_OI bar plot per imputation method

## Benchmarking one pipeline against a spike-in

With a known ground truth, how well a pipeline recovers the expected
changes. A positive counts only when it is significant and in the
expected direction.

- [`benchmarking_proteomics()`](https://sciordia.github.io/NADIA/reference/benchmarking_proteomics.md)
  : Benchmark Proteomics Differential Expression Results
- [`compute_benchmark_metrics()`](https://sciordia.github.io/NADIA/reference/compute_benchmark_metrics.md)
  : Compute benchmark metrics for all comparisons
- [`compute_opdea_metrics()`](https://sciordia.github.io/NADIA/reference/compute_opdea_metrics.md)
  : Compute OpDEA metrics for all comparisons
- [`compute_dispersion_metrics()`](https://sciordia.github.io/NADIA/reference/compute_dispersion_metrics.md)
  : Compute dispersion metrics by Comparison x Species
- [`benchmark_heatmap_gg()`](https://sciordia.github.io/NADIA/reference/benchmark_heatmap_gg.md)
  : Performance Heatmap (ggplot2)
- [`benchmark_confusion_gg()`](https://sciordia.github.io/NADIA/reference/benchmark_confusion_gg.md)
  : Confusion Matrix Heatmap (ggplot2)
- [`benchmark_confusion_overall_gg()`](https://sciordia.github.io/NADIA/reference/benchmark_confusion_overall_gg.md)
  : Overall Confusion Matrix Heatmap (ggplot2)
- [`benchmark_metrics_bars_gg()`](https://sciordia.github.io/NADIA/reference/benchmark_metrics_bars_gg.md)
  : Grouped Metrics Bar Chart (ggplot2)
- [`benchmark_auc_bars_gg()`](https://sciordia.github.io/NADIA/reference/benchmark_auc_bars_gg.md)
  : AUC Bar Chart (ggplot2)
- [`benchmark_roc_gg()`](https://sciordia.github.io/NADIA/reference/benchmark_roc_gg.md)
  : ROC Curves by Comparison (ggplot2)
- [`benchmark_signif_bars_gg()`](https://sciordia.github.io/NADIA/reference/benchmark_signif_bars_gg.md)
  : Significant Proteins Stacked Bars (ggplot2)
- [`benchmark_volcano_hc()`](https://sciordia.github.io/NADIA/reference/benchmark_volcano_hc.md)
  : Benchmark Volcano Plot (Highcharter)
- [`benchmark_volcano_hc_list()`](https://sciordia.github.io/NADIA/reference/benchmark_volcano_hc_list.md)
  : List of Benchmark Volcano Plots (Highcharter)

## Ranking several pipelines

Comparing normalisation and imputation combinations against each other,
with the OpDEA ranking and an extended eleven-metric one.

- [`benchmarking_multiple()`](https://sciordia.github.io/NADIA/reference/benchmarking_multiple.md)
  : Multiple-method OpDEA Benchmarking
- [`bm_compute_ranking()`](https://sciordia.github.io/NADIA/reference/bm_compute_ranking.md)
  : Compute OpDEA ranking across multiple methods
- [`bm_compute_extended_ranking()`](https://sciordia.github.io/NADIA/reference/bm_compute_extended_ranking.md)
  : Compute extended ranking across 11 classification + OpDEA metrics
- [`bm_compute_ranking_by_comparison()`](https://sciordia.github.io/NADIA/reference/bm_compute_ranking_by_comparison.md)
  : Compute OpDEA ranking separately for each comparison
- [`import_opdea_results()`](https://sciordia.github.io/NADIA/reference/import_opdea_results.md)
  : Import OpDEA metrics from multiple benchmark result folders
- [`import_benchmark_metrics()`](https://sciordia.github.io/NADIA/reference/import_benchmark_metrics.md)
  : Import benchmark metrics from multiple benchmark result folders
- [`import_classified_results()`](https://sciordia.github.io/NADIA/reference/import_classified_results.md)
  : Import classified results from multiple benchmark result folders
- [`import_confusion_results()`](https://sciordia.github.io/NADIA/reference/import_confusion_results.md)
  : Import confusion matrices from multiple benchmark result folders
- [`bm_plot_ranking_bars()`](https://sciordia.github.io/NADIA/reference/bm_plot_ranking_bars.md)
  : Bar chart of final ranks
- [`bm_plot_ranking_heatmap()`](https://sciordia.github.io/NADIA/reference/bm_plot_ranking_heatmap.md)
  : Heatmap of metric ranks by method
- [`bm_plot_ranking_bars_by_comp()`](https://sciordia.github.io/NADIA/reference/bm_plot_ranking_bars_by_comp.md)
  : Bar chart of final ranks for a single comparison
- [`bm_plot_ranking_heatmap_by_comp()`](https://sciordia.github.io/NADIA/reference/bm_plot_ranking_heatmap_by_comp.md)
  : Heatmap of metric ranks for a single comparison
- [`bm_plot_extended_ranking_bars()`](https://sciordia.github.io/NADIA/reference/bm_plot_extended_ranking_bars.md)
  : Horizontal bar chart of extended combined ranking
- [`bm_plot_extended_ranking_heatmap()`](https://sciordia.github.io/NADIA/reference/bm_plot_extended_ranking_heatmap.md)
  : Heatmap of extended ranking (11 metrics)
- [`bm_plot_metrics_comparison()`](https://sciordia.github.io/NADIA/reference/bm_plot_metrics_comparison.md)
  : Boxplot comparison of metric distributions across methods
- [`bm_plot_metrics_heatmap()`](https://sciordia.github.io/NADIA/reference/bm_plot_metrics_heatmap.md)
  : Heatmap of aggregated metric values (not ranks)
- [`bm_plot_confusion_stacked()`](https://sciordia.github.io/NADIA/reference/bm_plot_confusion_stacked.md)
  : Stacked bar chart of TP/FP/FN/TN (OpDEA Figure 5 style)
- [`bm_plot_roc()`](https://sciordia.github.io/NADIA/reference/bm_plot_roc.md)
  : ROC Curves by Method for a Single Comparison

## Interactive figures

Highcharts volcano plots, boxplots and PCA. Each returns a named list,
one entry per comparison or mode.

- [`volcano_highchart_list()`](https://sciordia.github.io/NADIA/reference/volcano_highchart_list.md)
  : Interactive Volcano Plot with Highcharter
- [`boxplot_highchart_list()`](https://sciordia.github.io/NADIA/reference/boxplot_highchart_list.md)
  : Interactive Highcharts Boxplot for Proteomics
- [`pca_highchart()`](https://sciordia.github.io/NADIA/reference/pca_highchart.md)
  : Interactive PCA Plot with Highcharts
- [`pca_highchart_list()`](https://sciordia.github.io/NADIA/reference/pca_highchart_list.md)
  : Build a List of PCA Plots for Multiple Subsets
- [`build_pca_scores()`](https://sciordia.github.io/NADIA/reference/build_pca_scores.md)
  : Build a PCA scores data frame from data in long format

## Static figures

Heatmaps through tidyHeatmap and ComplexHeatmap.

- [`proteomics_heatmap()`](https://sciordia.github.io/NADIA/reference/proteomics_heatmap.md)
  : Create a tidyHeatmap Heatmap for Proteomics
- [`proteomics_heatmap_list()`](https://sciordia.github.io/NADIA/reference/proteomics_heatmap_list.md)
  : Build a List of Heatmaps for Several Subsets

## Expression patterns

Fuzzy clustering of profiles across conditions, and the join that turns
it into the input for a functional analysis.

- [`pattern_profiler_analysis()`](https://sciordia.github.io/NADIA/reference/pattern_profiler_analysis.md)
  : Pattern Profiler Analysis
- [`summarize_pattern_profiler()`](https://sciordia.github.io/NADIA/reference/summarize_pattern_profiler.md)
  : Summary of the Pattern Profiler data
- [`read_pattern_profiler_data()`](https://sciordia.github.io/NADIA/reference/read_pattern_profiler_data.md)
  : Read Pattern Profiler data from a parquet file
- [`cluster_profile_highchart()`](https://sciordia.github.io/NADIA/reference/cluster_profile_highchart.md)
  : Cluster profile plot with Highcharts
- [`cluster_profile_highchart_list()`](https://sciordia.github.io/NADIA/reference/cluster_profile_highchart_list.md)
  : List of cluster profile plots with Highcharts
- [`cluster_centroids_highchart()`](https://sciordia.github.io/NADIA/reference/cluster_centroids_highchart.md)
  : Centroid plot for all the clusters
- [`deps_with_clusters()`](https://sciordia.github.io/NADIA/reference/deps_with_clusters.md)
  : Differential abundance with the Pattern Profiler cluster of each
  protein

## Interactive tables

Reactable tables for reviewing results. `*_reactable()` returns the
table alone; `*_widget()` returns a complete page with a toolbar and
Excel export.

- [`results_list_reactable()`](https://sciordia.github.io/NADIA/reference/results_list_reactable.md)
  : Interactive Reactable Table for Differential Expression Results
- [`results_list_widget()`](https://sciordia.github.io/NADIA/reference/results_list_widget.md)
  : Full Widget with Filters, Search, Export and CSS
- [`protein_list_reactable()`](https://sciordia.github.io/NADIA/reference/protein_list_reactable.md)
  : Interactive Reactable Table for Protein_ID Data (post-Spectronaut)
- [`protein_list_widget()`](https://sciordia.github.io/NADIA/reference/protein_list_widget.md)
  : Full Widget for the Protein_ID Table with Filters, Search and Excel
  Export
- [`quant_list_widget()`](https://sciordia.github.io/NADIA/reference/quant_list_widget.md)
  : Full Widget for the Protein_QUANT Table with the log2 Matrix
  Appended
- [`summary_list_widget()`](https://sciordia.github.io/NADIA/reference/summary_list_widget.md)
  : Interactive Reactable Table for the Sample Metadata

## The .nadia archive

One DuckDB file holding a whole analysis: data, parameters and
provenance, with the reconstruction stored as SQL views rather than
duplicated rows.

- [`write_nadia()`](https://sciordia.github.io/NADIA/reference/write_nadia.md)
  :

  Write a complete analysis to a single `.nadia` file

- [`read_nadia()`](https://sciordia.github.io/NADIA/reference/read_nadia.md)
  :

  Read a `.nadia` file

- [`nadia_result()`](https://sciordia.github.io/NADIA/reference/nadia_result.md)
  :

  Rebuild the processing result from a `.nadia` file

- [`nadia_preprocessing()`](https://sciordia.github.io/NADIA/reference/nadia_preprocessing.md)
  :

  Rebuild the preprocessing object from a `.nadia` file

- [`nadia_tables()`](https://sciordia.github.io/NADIA/reference/nadia_tables.md)
  :

  List the contents of a `.nadia` file without loading it

- [`nadia_connect()`](https://sciordia.github.io/NADIA/reference/nadia_connect.md)
  :

  Open a connection to a `.nadia` file

- [`nadia_export_parquet()`](https://sciordia.github.io/NADIA/reference/nadia_export_parquet.md)
  :

  Export the tables of a `.nadia` file to Parquet

- [`nadia_add_pattern_profiler()`](https://sciordia.github.io/NADIA/reference/nadia_add_pattern_profiler.md)
  :

  Add a Pattern Profiler run to an existing `.nadia` file

- [`nadia_pattern_profiler()`](https://sciordia.github.io/NADIA/reference/nadia_pattern_profiler.md)
  :

  Rebuild the Pattern Profiler table from a `.nadia` file

- [`nadia_deps_with_clusters()`](https://sciordia.github.io/NADIA/reference/nadia_deps_with_clusters.md)
  :

  Differential abundance with clusters, from a `.nadia` file

## Running NADIA from an assistant

The agent skill shipped with the package.

- [`nadia_skill_path()`](https://sciordia.github.io/NADIA/reference/nadia_skill_path.md)
  : Locate the NADIA agent skill

## Example data

A preprocessed DIA dataset: three conditions, four replicates, 2,000
protein groups, 8.1 % missing. The trimmed reports for all four input
formats live in `inst/extdata/`.

- [`nadia_dia`](https://sciordia.github.io/NADIA/reference/nadia_dia.md)
  : Example DIA experiment, already preprocessed
