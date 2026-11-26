# Compute missing-value percentages per protein

Four percentages, all measured on the normalized assay, that is, before
imputation, and all going from the widest scope to the narrowest:

- MissGlobal: % NAs over every sample in the experiment, including the
  conditions that take no part in the comparison. It does not depend on
  the comparison, so a given protein carries the same value in all of
  them.

- MissComp: % NAs over the replicates of the two conditions being
  compared only. With more than two conditions this is not the same as
  MissGlobal.

- MissCND1: % NAs over the samples of Cond1 (numerator, before "-")

- MissCND2: % NAs over the samples of Cond2 (denominator, after "-")

## Usage

``` r
.compute_missing_pct(se, DEPs_results, norm_assay_name)
```

## Arguments

- se:

  SummarizedExperiment with normalized assay

- DEPs_results:

  Data frame with DE results (must have Protein.IDs, Comparison)

- norm_assay_name:

  Name of the normalized assay

## Value

DEPs_results with MissGlobal, MissComp, MissCND1, MissCND2 columns added
