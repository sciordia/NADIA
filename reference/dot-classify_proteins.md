# Classify proteins for a single comparison

Assigns truth labels (immutable, biological identity) and predicted
labels.

- Species in expected_values -\> truth = 1 (spike-in / expected change)

- Species NOT in expected_values -\> truth = 0 (background)

- Positive species: predicted = 1 only if significant AND direction
  matches expected_logFC sign (TP); significant with wrong sign stays
  predicted = 0 (FN, a detection failure – NOT relabelled as a false
  positive).

- Negative species: predicted = 1 if significant (any direction -\> FP)
  `truth` is never modified by the prediction, so
  .compute_auc/.compute_pauc receive the biological ground truth
  (uncontaminated) for pROC.

## Usage

``` r
.classify_proteins(de_res_comp, ev_comp, alpha, lfc_thr, p_col)
```

## Arguments

- de_res_comp:

  DE results for one comparison

- ev_comp:

  Expected values for one comparison

- alpha:

  Significance threshold

- lfc_thr:

  Log fold-change threshold

- p_col:

  P-value column name

## Value

Data frame with added columns: truth, predicted, classification,
is_significant, direction_error

## Details

Two diagnostic columns are recorded alongside the classification and
take no part in any metric. `is_significant` is the significance rule
actually used here (p-value AND fold-change threshold), stored so that
tables and plots cannot drift from the classification by recomputing it.
`direction_error` separates the two kinds of FN: a change that was
missed altogether, and a change that was called significant with the
sign inverted.
