# =============================================================================
# Helpers shared by the test files
# =============================================================================

#' A small SummarizedExperiment with a controlled amount of missingness
#'
#' @param n_prot   Number of proteins.
#' @param n_per    Replicates per condition.
#' @param conds    Condition names.
#' @param na_prop  Fraction of cells set to NA.
#' @param seed     RNG seed.
.nadia_toy_se <- function(n_prot = 60, n_per = 4,
                          conds = c("A", "B", "D"),
                          na_prop = 0.1, seed = 1L) {
    set.seed(seed)
    samples <- paste0(rep(conds, each = n_per), "_", seq_len(n_per))
    n_col   <- length(samples)

    # A baseline abundance per protein plus a real effect in half of them, so
    # that a differential expression test has something to find.
    base <- matrix(rnorm(n_prot * n_col, mean = 20, sd = 0.4),
                   nrow = n_prot, ncol = n_col,
                   dimnames = list(paste0("p", seq_len(n_prot)), samples))
    effect <- rep(c(2, 0), length.out = n_prot)
    base[, grepl("^B_", samples)] <- base[, grepl("^B_", samples)] + effect

    if (na_prop > 0) {
        n_na <- round(na_prop * length(base))
        base[sample(length(base), n_na)] <- NA
    }

    SummarizedExperiment::SummarizedExperiment(
        assays  = list(cycloess = base),
        colData = S4Vectors::DataFrame(
            Column    = samples,
            Condition = factor(rep(conds, each = n_per), levels = conds),
            Replicate = rep(seq_len(n_per), times = length(conds)),
            row.names = samples))
}


#' Drive `.extract_limma_results()` into a chosen corner of its classifier
#'
#' The fold-change classification lives inside `.extract_limma_results()`, which
#' takes a limma fit. Some of the cases worth testing — notably a `logFC` of
#' exactly zero together with a significant p-value — cannot arise from real
#' data, since a zero coefficient gives t = 0 and p = 1.
#'
#' So a genuine fit is built first and its coefficients and p-values are then
#' overwritten with the values under test. The classifier itself is the real one;
#' only its input is constructed.
#'
#' @param logFC   Vector of fold changes to feed in.
#' @param adj     Vector of adjusted p-values to feed in.
#' @param lfc_thr Threshold passed as `logFC_up` / `-logFC_down`.
#' @param alpha   Significance threshold.
#' @return The data frame returned by `.extract_limma_results()`.
.nadia_fake_de <- function(logFC, adj, lfc_thr = 0, alpha = 0.05) {
    stopifnot(length(logFC) == length(adj))
    n <- length(logFC)

    y      <- matrix(rnorm(n * 6), nrow = n,
                     dimnames = list(paste0("p", seq_len(n)), NULL))
    design <- cbind(Intercept = 1, Group = rep(0:1, each = 3))
    fit    <- limma::eBayes(limma::lmFit(y, design))

    # Keep the fit's structure; replace the two quantities the classifier reads.
    fit$coefficients[, "Group"] <- logFC
    fit$p.value[,     "Group"]  <- adj
    fit <- fit[, "Group", drop = FALSE]

    NADIA:::.extract_limma_results(
        fit, comparisons = "B-A",
        logFC_up = lfc_thr, logFC_down = -lfc_thr,
        alpha = alpha, p_adj = TRUE)
}
