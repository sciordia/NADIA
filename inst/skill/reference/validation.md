# Validation: justifying the pipeline you chose

Three tools. Which one applies depends on whether the experiment has a known
ground truth.

| Situation | Tool |
|---|---|
| No ground truth | `normalization_metrics()`, `imputation_metrics()` |
| Spike-in with expected changes, one pipeline | `benchmarking_proteomics()` |
| Spike-in, several pipelines to rank | `benchmarking_multiple()` |

## Without ground truth

### Normalisation

```r
methods <- c("log2Norm", "GlobalMedian", "GlobalMean", "eqmedians", "vsn",
             "medianNorm", "meanNorm", "quantile", "Rlr", "MAD", "cycloess",
             "quantile.robust")

nm <- normalization_metrics(
  se,                          # a SummarizedExperiment; res$se_proc works
  methods       = methods,     # WITHOUT THIS IT COMPARES SOMETHING ELSE. See below.
  condition_col = "Condition",
  base_assay    = "log2",
  output_dir    = "results/norm_metrics",
  seed          = 42L,
  verbose       = FALSE
)

nm$final_rank      # rank 1 = best
nm$final_ranking   # the plot
```

**`methods` is what makes this a comparison of normalisation methods.** It is
`NULL` by default, and then the function ranks *the assays already in the
object* instead — on a processed result that means `raw`, `log2`, the
normalised assay, `BERT` if you corrected batches, and the imputed one. That
ranking is meaningless for choosing a method: those are stages of one pipeline,
not competing alternatives, and the batch-corrected and imputed assays win it
by construction because they have the least within-group variability. Pass
`methods` and the 12 are generated and compared properly.

With `methods` supplied it compares them on:

- **PCV / PMAD / PEV** — within-group variability, per protein. Lower is better.
- **Intragroup correlation** — higher is better.
- **Rank_Sep** — group separation along PC1, as an F-ratio.

`nm_rank_final` is the mean of the PCV/PMAD/PEV/correlation ranks plus the
separation rank. Report the ranking, not just the winner: methods often tie
within noise, and then the choice can be made on other grounds.

### Imputation

```r
im <- imputation_metrics(
  se,
  assay_name = "cycloess",     # the NORMALISED assay, pre-imputation
  na_prop    = 0.2,            # fraction of observed cells masked
  seed       = 42L,
  output_dir = "results/imp_metrics",
  verbose    = FALSE
)

im$metrics_table       # the metrics, with Rank_Mean as the summary column
```

Ground-truth simulation: it takes the complete rows, hides a known fraction of
values, re-imputes, and compares. Metrics:

- **NRMSE** — recovery error. Lower is better.
- **SOR** — penalises features left un-imputed.
- **PSS** — structure preservation (needs `vegan`).
- **ACC_OI** — agreement on the masked cells only.
- **Rank_Mean** — the summary column of `im$metrics_table`.

16 methods by default.

**Read this ranking as a choice of MAR method, and nothing else.** The masking
hides values that were *observed*, so the artificial gaps are missing at random
by construction. MNAR methods are built for a kind of missingness the
simulation never creates, and they are punished accordingly. On the DIA-NN
example the table looks like this:

| Method | NRMSE | Rank_Mean |
|---|---|---|
| Impseqrob | 0.077 | 2.75 |
| knn | 0.101 | 3.00 |
| … | | |
| min | 3.514 | 14.00 |
| halfmin | 3.896 | 15.00 |
| zero | 8.941 | 16.00 |

`min` placing 14th of 16 does **not** mean `min` is a bad MNAR stage. It means
this benchmark cannot evaluate one. Reading the table as a straight ranking
leads to dropping exactly the MNAR stage that the missingness structure of the
same dataset calls for — the two halves of the evidence contradict each other
unless you know why.

So use it like this:

- **MAR stage of `combo`** — take the winner among `Impseqrob`, `Impseq`,
  `knn`, `bpca`, `missForest`, `mice`, `MLE`. That is what the table measures.
- **MNAR stage** — choose it from the missingness structure (see
  `reference/method-choice.md`) or, when a spike-in exists, from the benchmark
  below. Never from this table.
- **A single method for everything** — only defensible when the missingness
  shows no intensity dependence. Otherwise use a hybrid.

`ACC_OI` comes back `NA` for `min`, `halfmin` and `zero` because they fill
every masked cell with the same constant, leaving nothing to correlate. That is
a further sign the exercise does not fit them, not a defect.

## With a spike-in

`benchmarking_proteomics()` scores one pipeline against known changes.

```r
bm <- benchmarking_proteomics(
  de_res          = res$DEPs_results,
  expected_values = expected,     # the expected direction per species/comparison
  species_df      = species_map,  # Protein.IDs + Species; or a Species column in de_res
  alpha           = 0.05,
  lfc_thr         = 0,
  output_dir      = "results/benchmark",
  verbose         = FALSE
)

bm$metrics_table     # Sensitivity, Specificity, Precision, F1, AUC, MCC, ...
bm$opdea_metrics     # nMCC, G_mean, pAUC at FPR 0.01 / 0.05 / 0.10
```

The classification convention matters when you interpret the output:

- `truth` is the **immutable biological identity**: spike-in species = 1,
  background = 0. It is never contaminated by the prediction, so the AUC is
  honest.
- A positive counts as a true positive only if it is significant **and** in the
  expected direction. Significant with the wrong sign is a **false negative** —
  a detection failure — not a false positive.

**Full AUC and pAUC can disagree, and the low-FPR region is usually the one you
care about.** A method can win on the whole ROC curve and lose in the first
1 % of false positives, which is the regime an experiment actually operates in.
Report both.

### Ranking several pipelines

Run each combination into its own `output_dir`, then:

```r
bmm <- benchmarking_multiple(
  results_dir = "results",       # reads benchmark_opdea_metrics.tsv from subfolders
  output_dir  = "results/ranking",
  verbose     = FALSE
)

bmm$mean_ranking       # OpDEA ranking, 5 metrics
bmm$extended_ranking   # 11 metrics, when the per-comparison tables are present
```

The OpDEA ranking (Peng et al., *Nature Communications* 2024) averages the
ranks of nMCC, G_mean and the three pAUCs. An extended 11-metric ranking is
also produced when the per-comparison metrics are available.

`pROC` is needed for AUC/pAUC and `vegan` for PSS; both are in `Suggests`.

## What to say afterwards

State the criterion, not only the winner:

> Normalisation and imputation were selected with `normalization_metrics()`
> and `imputation_metrics()` on this dataset; `vsn` ranked first of the 12
> normalisation methods and `Impseqrob` first of 16 imputation methods
> (`Rank_Mean`). No spike-in was available, so recovery of true changes could
> not be measured directly.

Do not assume the default wins. On the TMT example shipped with the package,
`vsn` ranks first on all five criteria and `cycloess`, the default, ranks
sixth of twelve.

That last sentence is not a caveat to be dropped. Without ground truth these
metrics measure internal consistency, not correctness.
