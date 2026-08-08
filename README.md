# NADIA

<div align="justify">

*Differential abundance analysis in quantitative proteomics, with explicit
treatment of missing values.*

NADIA is an R package that provides a complete workflow for the differential
abundance analysis of proteins. It processes quantitative data from DIA, TMT and
label-free (DDA) experiments, with particular attention to the assessment and
treatment of the missing values that routinely appear in DIA proteomics data.

The name is drawn from the package title, *Missing Value-Aware Differential
Abundance Analysis of DIA Proteomics Data*, and it reads on two levels: NA is how
R represents a missing value, and DIA is data-independent acquisition, the
setting in which the package was first developed.

## Why missing values matter

DIA proteomics datasets routinely contain an appreciable proportion of missing
values. These values do not always represent random measurement failures. A
protein may fall below the limit of detection in one condition and be
quantifiable in another, producing precisely the pattern that a differential
abundance analysis sets out to identify.

For this reason, removing every protein that contains missing values may discard
relevant biological signal. Equally, imputing them without considering the likely
mechanism of missingness may attenuate real differences or introduce artificial
ones.

This relationship can be observed in the example dataset included with NADIA. The
table below crosses the percentage of missing values of each protein in a
comparison with its call in the differential abundance analysis:

| Missing values | Up | Down | No Change |
|---|---|---|---|
| 0 % | 795 | 209 | 3761 |
| 1–25 % | 90 | 57 | 239 |
| 26–50 % | 70 | 619 | 54 |
| > 50 % | 30 | 52 | 15 |

Most proteins with no missing values show no significant change. Among proteins
with 26–50 % missing values, by contrast, there is a marked asymmetry towards
negative changes: 619 are called Down, against 70 called Up.

This pattern is consistent with a substantial abundance-dependent component of
missingness: when a protein falls below the limit of detection in one condition,
the missing values concentrate precisely in the group with the lower abundance.
NADIA is designed to retain, examine and explicitly handle this information,
rather than to remove it or impute it without assessing the consequences.

## What NADIA covers

From the quantification report to reproducible results, ready to explore or
export:

1. *Import and preprocessing* — Reads quantification reports from Spectronaut
   and DIA-NN, as well as TMT and label-free (DDA) experiments processed with
   Proteome Discoverer. Every format is converted into a common
   `proteomics_data` object, so that the rest of the workflow is independent of
   the source software.
2. *Processing and differential abundance* — Normalises the data through 13
   methods, allows batch effects to be diagnosed and corrected, and imputes
   missing values with 19 individual methods or with hybrid strategies that treat
   values assumed to be MAR and MNAR separately. Differential abundance is
   analysed with `limma` or `limpa`.
3. *Method assessment and selection* — Compares normalisation and imputation
   methods using metrics computed on the data themselves. Where a known reference
   is available, such as a *spike-in* experiment, NADIA evaluates sensitivity,
   specificity and the recovery of the expected changes, and ranks normalisation
   and imputation combinations through the OpDEA approach.
4. *Result analysis and visualisation* — Produces interactive figures with
   Highcharts, including *volcano plots*, boxplots, PCA and protein cluster
   profiles. It also produces static figures with ggplot2 and ComplexHeatmap, and
   identifies proteins with similar profiles through fuzzy clustering.
5. *Tables and export* — Presents the results in interactive tables and exports
   the processed matrices, the differential abundance results and the data used
   by the visualisations in formats suited to archiving or further analysis.

## Method choice is not cosmetic

Normalisation and imputation are not merely technical steps. They can
substantially change the outcome of a differential abundance analysis.

The table below compares four normalisation and imputation combinations applied
to the same dataset, using the same statistical model and the same significance
and fold-change thresholds:

| Pipeline | Up | Down | No Change |
|---|---|---|---|
| `cycloess` + `combo` | 985 | 937 | 4069 |
| `cycloess` + `knn` | 872 | 179 | 4940 |
| `quantile` + `combo` | 902 | 2825 | 2264 |
| `log2` + `min` | 761 | 778 | 4452 |

*Note*: the counts correspond to protein–comparison results; the same protein may
appear in more than one comparison.

The number of results called Down ranges from 179 to 2,825, a difference of
roughly sixteen-fold. The total number of results called differentially abundant
also changes considerably: from around 1,050 with `cycloess` + `knn` to more than
3,700 with `quantile` + `combo`.

This variability does not by itself indicate which pipeline is correct. It shows
that the choice of normalisation and imputation method should be assessed on the
data, rather than adopted by convention alone. It also shows why a reproducible
analysis must state which methods and parameters produced the final table:
without that information, two analyses of the same quantification matrix can lead
to very different conclusions.

`vignette("choosing-methods")` shows how to compare the alternatives when no
known reference is available. When the experiment contains a spike-in with
expected changes, `vignette("benchmarking")` evaluates them against that
reference.

## Installation

NADIA requires R 4.6 or later. While the package is under development and not yet
available on Bioconductor, it can be installed from its GitHub repository:

```r
# install.packages("remotes")
remotes::install_github("sciordia/NADIA", build_vignettes = TRUE)
library(NADIA)
```

The mandatory dependencies declared in the `Imports:` field of `DESCRIPTION` are
installed automatically. The dependencies declared in `Suggests:` are only needed
for particular optional functions or methods. `mice`, for instance, is required
when that imputation method is selected, `pROC` computes the AUC and pAUC
benchmarking metrics, and `Mfuzz` is used by the Pattern Profiler.

When an optional dependency is missing, NADIA names the required package and
explains how to install it. All the package dependencies can be checked or
installed through the helper functions included in the repository:

```r
source("install_dependencies.R")         # from a clone of the repository
install_nadia_deps(dry_run = TRUE)       # only report what is missing
install_nadia_deps(optional = FALSE)     # only the essentials
```

Once NADIA is available on Bioconductor, the recommended installation will be:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("NADIA")
```

### A note on Highcharts licensing

NADIA uses the R package `highcharter`, distributed under the MIT licence, to
generate some of its interactive visualisations. highcharter acts as an interface
to the Highcharts JavaScript library, which is distributed under its own
licensing terms.

NADIA's licence neither grants nor implies a licence to use Highcharts. Depending
on the context in which the visualisations are used — personal, educational,
academic, institutional, governmental or commercial — a specific Highcharts
licence may be required. Users should consult the current
[Highsoft terms](https://www.highcharts.com/license) and ensure that their use
complies with the applicable licence.

NADIA's statistical analysis does not depend on producing interactive figures.
The static visualisations produced with ggplot2 and ComplexHeatmap can be used
instead, and the complete workflow can be run without generating a single
Highcharts figure.

## A first analysis

NADIA includes `nadia_dia`, a preprocessed DIA dataset of three conditions, four
replicates per condition and 2,000 protein groups. Roughly 8 % of its intensities
are missing values.

The example below normalises the data with `cycloess`, applies the hybrid `combo`
imputation strategy and runs the differential abundance analysis with `limma`:

```r
library(NADIA)

# Preprocessed example dataset
data(nadia_dia)

res <- process_proteomics(
    nadia_dia,
    norm_method = "cycloess",
    imp_method  = "combo",
    mar_method  = "Impseqrob",   # default value in combo imp_method
    mnar_method = "min",         # default value in combo imp_method
    de_method   = "limma"
)

head(res$DEPs_results)

# Or starting from an exported Spectronaut quantification report
prep <- preprocess_spectronaut(
  system.file("extdata", "nadia_dia_report.tsv.gz", package = "NADIA"),
  condition_order = c("A", "B", "D"))

res_from_report <- process_proteomics(
    prep,
    norm_method = "cycloess",
    imp_method  = "combo",
    de_method   = "limma"
)
```

In this call, `combo` activates a two-stage hybrid imputation strategy. If
neither `mar_method` nor `mnar_method` is specified, `process_proteomics()` uses
its defaults: `Impseqrob` for the values classified as MAR and `min` for those
classified as MNAR. The resulting imputed matrix is stored in the assay
`Impseqrob_min`. Both methods can be replaced through the `mar_method` and
`mnar_method` arguments.

`process_proteomics()` coordinates the main stages of the analysis and returns a
list containing the processed object, the differential abundance results and the
tables used by the visualisation modules. `res$DEPs_results` holds one row per
protein group and comparison analysed.

The function writes no files unless a directory is given through `export_dir`.

Alongside `logFC`, `P.Value`, `adj.P.Val` and the call recorded in `Change`, the
results table retains information about the values that were missing before
imputation:

- `MissGlobal` — the percentage of missing values of that protein group across
  every sample in the experiment, including the conditions that take no part in
  the comparison.
- `MissComp` — the percentage across the replicates of the two conditions being
  compared, and only those. With more than two conditions this is not the same
  figure as the one above.
- `MissCND1` and `MissCND2` — the percentage of missing values in each of the two
  conditions separately.

These columns do not on their own determine whether a change is biologically
real, but they show how far its estimate rests on observed intensities and how
far it depends on imputed values. A protein entirely absent in one condition and
present in the other, for example, represents a presence–absence pattern or
quantification below the limit of detection. In that case both the `logFC` and
its significance must be interpreted in the light of the imputation method used.

The second part of the example shows how to build the same kind of object from a
quantification report exported by Spectronaut. The resulting `prep` object can be
passed directly to `process_proteomics()` to continue with the same workflow.

## Documentation and suggested route

NADIA includes nine vignettes documenting the complete workflow and its main
modules. As a first approach, start with `vignette("NADIA")`, which walks through
an analysis from the quantification data to the differential abundance results
and their visualisations.

The remaining vignettes can be read independently, according to the question at
hand:

| Vignette | Question it answers |
|---|---|
| `NADIA` | How do I run a complete analysis with NADIA? |
| `input-formats` | How do I import DIA, TMT or label-free (DDA) data? |
| `missing-values` | What do the missing values mean and how should I impute them? |
| `choosing-methods` | How do I compare normalisation and imputation methods when I do not know the true answer? |
| `benchmarking` | How do I evaluate pipelines when I have a *spike-in* experiment with expected changes? |
| `batch-correction` | How do I diagnose, correct and verify a possible batch effect? |
| `pattern-profiler` | How do I identify proteins with similar profiles across conditions? |
| `visualization` | How do I customise the interactive and static figures? |
| `results-and-export` | How do I export, archive and share the analysis results? |

Every vignette contains executable code and generates its results from the
datasets shipped with the package, so the examples can be reproduced and adapted
to other experiments.

The installed vignettes can be browsed with `browseVignettes("NADIA")`.

## Package architecture

NADIA is organised modularly, following the main stages of the workflow. Each
stage has its own functions, which can be used independently in a custom
workflow, while `process_proteomics()` coordinates them to run a complete
analysis.

### Import and preprocessing

The preprocessing functions convert the different input formats into a common
structure:

- `preprocess_spectronaut()` handles Spectronaut and DIA-NN quantification
  reports in long format.
- `preprocess_tmt()` handles TMT reports exported by Proteome Discoverer.
- `preprocess_lfq()` handles label-free reports exported by Proteome Discoverer,
  together with their experimental annotation.

All three return a `proteomics_data` object with the same basic structure: sample
metadata, protein annotation and the quantification matrix. The later stages
therefore apply in the same way regardless of the source software or
quantification design.

### Main processing

`process_proteomics()` coordinates the central stages of the analysis:

1. normalisation
2. optional batch-effect correction
3. imputation of missing values
4. differential abundance analysis

From normalisation onwards the data live in a `SummarizedExperiment`, the
container most Bioconductor packages use for omics data: assays of identical
shape with protein groups in rows and samples in columns, `rowData` with one row
of annotation per protein group, and `colData` with one row of metadata per
sample, all indexed together so that subsetting the object subsets the annotation
with it. NADIA keeps one assay per pipeline stage, so the input matrix, its
logarithm, the normalised matrix and the imputed matrix all sit in the same
object. The result also returns the differential abundance table and the inputs
used by the visualisation modules.

The main stages can be run individually through `normalize_proteomics()`,
`batch_correct_proteomics()`, `impute_proteomics()` and
`de_analysis_proteomics()`.

### Assessment and benchmarking

The assessment modules are separate from the main processing:

- `normalization_metrics()` compares normalisation methods using observable
  properties of the data.
- `imputation_metrics()` evaluates imputation methods by simulating missing
  values.
- `benchmarking_proteomics()` compares one pipeline against the expected changes
  of a *spike-in* experiment.
- `benchmarking_multiple()` ranks multiple normalisation and imputation
  combinations.

### Visualisation and presentation of results

The visualisation functions mainly accept tables and data frames, so they can be
used both with the results of `process_proteomics()` and with results produced
elsewhere.

NADIA provides interactive volcano, boxplot, PCA and cluster-profile figures,
static figures through ggplot2 and ComplexHeatmap, and interactive tables built
on reactable.

## Example data and provenance

The example data shipped with NADIA come from quantitative proteomics experiments
acquired at the Proteomics Facility of the Centro Nacional de Biotecnología
(CNB-CSIC).

The package includes:

- `nadia_dia`, a preprocessed `proteomics_data` object that allows the main
  workflow to be run directly;
- trimmed DIA, TMT and label-free (DDA) quantification reports, stored in
  `inst/extdata/`, which demonstrate the import and preprocessing functions from
  their respective input formats;
- the additional information needed for the *spike-in* benchmarking examples.

The distributed reports contain a random subset of 2,000 protein groups drawn
from the full experiments. The selection was made at random, rather than keeping
only the best-covered proteins, so that the examples retain a realistic structure
of missing values. Selecting only the best-quantified proteins would have removed
much of the missingness NADIA exists to examine and handle.

The full experiments are unpublished and are part of neither the package nor the
repository. The script `inst/scripts/make_extdata.R` documents how the trimmed
reports and the `nadia_dia` object were derived from the original data, so the
operations used to prepare the example sets are on record even though the full
experiments are not distributed.

These data are provided solely to demonstrate, test and reproduce the behaviour
of the package. They should not be treated as reference datasets from which to
draw biological conclusions about the original experiments.

## License

NADIA is free and open-source software distributed under the terms of the GNU
General Public License, version 3 or any later version. The code may be used,
modified and redistributed in accordance with the conditions of that licence. See
the `LICENSE` file and the `License:` field of `DESCRIPTION` for further details.

Copyright © 2025–2026 Sergio Ciordia.

</div>
