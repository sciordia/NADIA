# =============================================================================
# Regressions for the bugs fixed in the 2026-07 code reviews
# =============================================================================
#
# Every test here corresponds to a defect documented in one of the
# CODE_REVIEW_*.md reports. They are the tests that would have caught those
# bugs, and their job is to stop them coming back. Each one names the defect it
# guards, because a bare assertion is easy to "fix" by changing the expectation.

# --- Normalisation metrics: per-protein, not per-group ------------------------

test_that(".nm_pcv/.nm_pmad/.nm_pev return one value PER PROTEIN", {
    # Defect: the three used colMeans() and so returned one value per group
    # instead of one per protein. The bug changed which normalisation method the
    # benchmark recommended (Rlr instead of cycloess), so it was not cosmetic.
    set.seed(1)
    mat <- matrix(rnorm(20 * 6, mean = 20), nrow = 20,
                  dimnames = list(paste0("p", 1:20), paste0("s", 1:6)))
    groups <- rep(c("A", "B"), each = 3)

    for (f in list(NADIA:::.nm_pcv, NADIA:::.nm_pmad, NADIA:::.nm_pev)) {
        out <- f(mat, groups)
        expect_length(out, nrow(mat))          # 20, not 2
        expect_false(length(out) == length(unique(groups)))
    }
})

test_that(".nm_pcv responds to per-protein variability, not to the grand mean", {
    # A per-protein statistic must change when one protein's replicates spread
    # out, even if the overall mean of the matrix is unchanged.
    mat <- matrix(rep(10, 12), nrow = 3, dimnames = list(paste0("p", 1:3), NULL))
    groups <- rep(c("A", "B"), each = 2)

    tight <- NADIA:::.nm_pcv(mat, groups)

    loose <- mat
    loose[2, ] <- c(5, 15, 5, 15)   # same mean, much larger spread
    spread <- NADIA:::.nm_pcv(loose, groups)

    expect_gt(spread[2], tight[2])
    expect_equal(spread[1], tight[1])   # the other proteins are untouched
})


# --- Imputation metrics: ACC_OI on the masked cells only ----------------------

test_that(".im_acc_oi is computed only on the masked cells", {
    # Defect: it correlated whole rows, so the untouched values dominated and
    # every method scored close to 1. Restricting it to the masked cells is what
    # makes the metric discriminative.
    set.seed(2)
    truth <- matrix(rnorm(20 * 6), nrow = 20)
    imp   <- truth
    mask  <- matrix(FALSE, nrow = 20, ncol = 6)
    mask[cbind(1:20, sample(6, 20, replace = TRUE))] <- TRUE

    # Imputed values are pure noise, uncorrelated with the truth.
    imp[mask] <- rnorm(sum(mask))

    acc <- NADIA:::.im_acc_oi(truth, imp, mask)

    # If the whole row were used, the untouched cells would push this near 1.
    expect_lt(abs(acc), 0.6)
})

test_that(".im_acc_oi is NA for a constant imputation", {
    # A method like `min` writes the same value into every gap of a protein. The
    # correlation with a constant is undefined, and NA is the correct answer:
    # a method that cannot vary cannot be scored on how well it varies.
    truth <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2)
    imp   <- truth
    mask  <- matrix(c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE), nrow = 2)
    imp[mask] <- 0    # constant

    expect_true(is.na(NADIA:::.im_acc_oi(truth, imp, mask)))
})

test_that(".im_acc_oi is high when the imputation recovers the truth", {
    set.seed(3)
    truth <- matrix(rnorm(20 * 6), nrow = 20)
    imp   <- truth
    mask  <- matrix(FALSE, nrow = 20, ncol = 6)
    mask[cbind(1:20, sample(6, 20, replace = TRUE))] <- TRUE
    imp[mask] <- truth[mask] + rnorm(sum(mask), sd = 0.05)   # nearly perfect

    expect_gt(NADIA:::.im_acc_oi(truth, imp, mask), 0.9)
})


# --- Quantile.robust: no NA-induced bias -------------------------------------

test_that(".norm_quantile_robust bias does not grow with the missingness rate", {
    # Defect: the reference was built with sort(na.last = TRUE), mixing i/m and
    # i/n quantiles, so a column was shifted in proportion to how many NAs it
    # had. The fix interpolates every column onto a common quantile grid.
    #
    # A single run cannot show this: at 30 % NA the standard error of the median
    # is itself around 0.2, so any one difference is mostly noise. What
    # distinguishes the fix from the bug is that the MEAN difference stays flat
    # as missingness rises, instead of growing with it. Measured over the same
    # design, the old code gave -0.27, -0.81 and -1.35 at 10 %, 30 % and 50 %.

    bias_at <- function(na_frac, reps = 15) {
        d <- vapply(seq_len(reps), function(s) {
            set.seed(s)
            n <- 300
            x <- matrix(rnorm(n * 3, mean = 20, sd = 2), nrow = n)
            x[, 2] <- x[, 1]                                  # exact copy ...
            x[sample(n, round(na_frac * n)), 2] <- NA         # ... with gaps
            out  <- NADIA:::.norm_quantile_robust(x)
            both <- !is.na(out[, 1]) & !is.na(out[, 2])
            median(out[both, 2]) - median(out[both, 1])
        }, numeric(1))
        mean(d)
    }

    b10 <- bias_at(0.10)
    b50 <- bias_at(0.50)

    # Both must be small in absolute terms ...
    expect_lt(abs(b10), 0.15)
    expect_lt(abs(b50), 0.20)

    # ... and, crucially, the bias must not scale with the missingness. Under
    # the old code b50 was about five times b10.
    expect_lt(abs(b50), 4 * abs(b10) + 0.15)
})


# --- Hopkins statistic: reproducible -----------------------------------------

test_that(".nm_hopkins is reproducible", {
    # Defect: set.seed() was called AFTER the random points were generated, so
    # the statistic was not reproducible despite taking a seed.
    set.seed(5)
    mat <- matrix(rnorm(40 * 6), nrow = 40)

    a <- NADIA:::.nm_hopkins(mat, seed = 42L)
    b <- NADIA:::.nm_hopkins(mat, seed = 42L)
    expect_identical(a, b)

    # And a different seed should generally give a different value, otherwise
    # the seed is being ignored altogether.
    d <- NADIA:::.nm_hopkins(mat, seed = 7L)
    expect_false(isTRUE(all.equal(a, d)))
})


# --- Fold-change classification: sign-guarded --------------------------------

test_that("a logFC of exactly 0 is 'No Change' even at threshold 0", {
    # Defect: with logFC_threshold = 0 both the Up and the Down assignment
    # matched logFC == 0, so the second overwrote the first and the protein was
    # silently classified as Down. The fix requires a strict sign.
    de <- .nadia_fake_de(logFC = c(0, 0.5, -0.5), adj = c(0.001, 0.001, 0.001))
    expect_identical(as.character(de$Change), c("No Change", "Up", "Down"))
})

test_that("classification requires significance AND the threshold", {
    de <- .nadia_fake_de(logFC = c(2, 2, 0.1), adj = c(0.001, 0.5, 0.001),
                         lfc_thr = 1)
    expect_identical(as.character(de$Change),
                     c("Up",          # significant and above threshold
                       "No Change",   # above threshold but not significant
                       "No Change"))  # significant but below threshold
})

test_that("Change always has the three levels in a fixed order", {
    # Downstream code (plots, tables, colour maps) relies on this vocabulary.
    de <- .nadia_fake_de(logFC = 1, adj = 0.001)
    expect_true(is.factor(de$Change))
    expect_identical(levels(de$Change), c("Up", "Down", "No Change"))
})


# --- softHybrid: the sigmoid and the MNAR weight -----------------------------

test_that(".sigmoid is monotone, bounded and centred on x0", {
    f <- NADIA:::.sigmoid
    expect_equal(f(0, k = 1, x0 = 0), 0.5)
    expect_gt(f(5, k = 1, x0 = 0), f(-5, k = 1, x0 = 0))
    expect_lt(f(-50, k = 1, x0 = 0), 1e-6)
    expect_gt(f(50, k = 1, x0 = 0), 1 - 1e-6)

    # All outputs stay inside [0, 1], which is what lets the value be used as a
    # blending weight without further clamping.
    vals <- f(seq(-20, 20, by = 0.5), k = 2, x0 = 3)
    expect_true(all(vals >= 0 & vals <= 1))
})
