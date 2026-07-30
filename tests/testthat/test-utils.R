# =============================================================================
# Shared internal helpers (R/utils.R)
# =============================================================================
#
# These were duplicated across six modules until Phase 3 and were unified into a
# single definition. The tests below pin down the behaviour that the unification
# had to preserve, including the deliberate design decisions.

test_that("%||% only catches NULL, not NA or empty vectors", {
    `%||%` <- NADIA:::`%||%`

    expect_identical(NULL %||% "fallback", "fallback")
    expect_identical("value" %||% "fallback", "value")

    # These are the cases the operator must NOT catch: NA and a zero-length
    # vector are values, not absences.
    expect_true(is.na(NA %||% "fallback"))
    expect_identical(character(0) %||% "fallback", character(0))
    expect_identical(list() %||% "fallback", list())
    expect_identical(FALSE %||% "fallback", FALSE)
})


# --- Colour helpers ----------------------------------------------------------

test_that(".normalize_hex handles the accepted input forms", {
    f <- NADIA:::.normalize_hex

    expect_identical(f("#0e6655"), "#0E6655")
    expect_identical(f("0e6655"),  "#0E6655")   # without the leading '#'
    expect_identical(f("#0E6655"), "#0E6655")   # already normalised

    # An 8-digit RRGGBBAA colour drops its alpha channel.
    expect_identical(f("#0E665580"), "#0E6655")
    expect_identical(f("0E665580"),  "#0E6655")
})

test_that(".normalize_hex is not vectorised, and fails loudly when given a vector", {
    # Every call site passes a single colour. The nchar() == 8 branch uses if(),
    # so a vector of length > 1 raises an error rather than returning a quietly
    # wrong answer. That is the behaviour worth locking in: this helper feeds
    # colour strings straight into figures, where a silently mangled value would
    # be very hard to trace back.
    expect_error(NADIA:::.normalize_hex(c("#0E6655", "#D55E00")),
                 "condition has length")
})

test_that(".hex_to_rgba produces a Highcharts rgba() string", {
    f <- NADIA:::.hex_to_rgba

    expect_identical(f("#FFFFFF", 1),   "rgba(255, 255, 255, 1.00)")
    expect_identical(f("#000000", 0),   "rgba(0, 0, 0, 0.00)")
    expect_identical(f("#0E6655", 0.5), "rgba(14, 102, 85, 0.50)")

    # The '#' is optional and the alpha is formatted to two decimals.
    expect_identical(f("0E6655", 0.125), "rgba(14, 102, 85, 0.12)")
})

test_that(".hex_to_rgba and .darken_hex have no default for their 2nd argument", {
    # This is a design decision from the unification, not an oversight. The
    # duplicated copies disagreed on the default; since all 7 call sites always
    # pass the argument, no default is declared so that a call without it fails
    # visibly instead of silently picking one module's convention.
    expect_error(NADIA:::.hex_to_rgba("#0E6655"))
    expect_error(NADIA:::.darken_hex("#0E6655"))
})

test_that(".darken_hex darkens monotonically and clamps at black", {
    f <- NADIA:::.darken_hex

    expect_identical(f("#FFFFFF", 0),   "#FFFFFF")   # no darkening
    expect_identical(f("#FFFFFF", 1),   "#000000")   # fully dark
    expect_identical(f("#FFFFFF", 0.5), "#808080")

    # Never below zero, whatever the factor.
    expect_identical(f("#0E6655", 2), "#000000")

    # Monotone: a larger factor is never brighter.
    lum <- function(h) sum(grDevices::col2rgb(h))
    expect_lte(lum(f("#0E6655", 0.6)), lum(f("#0E6655", 0.3)))
})


# --- Feature filtering -------------------------------------------------------

test_that(".adjp_col builds the expected column name", {
    expect_identical(NADIA:::.adjp_col("B-A"), "adjP_B-A")
    expect_identical(NADIA:::.adjp_col("Treat.1-Ctrl"), "adjP_Treat.1-Ctrl")
})

test_that(".get_feature_ids accepts 'target' and 'specific' as synonyms", {
    # The two unified copies used different names for the same mode; both must
    # keep working, because each public module propagates its own vocabulary.
    d <- data.frame(FeatureID = c("p1", "p2", "p3"),
                    sig_any   = c(TRUE, FALSE, TRUE),
                    `adjP_B-A` = c(0.01, 0.20, 0.04),
                    check.names = FALSE)

    a <- NADIA:::.get_feature_ids(d, mode = "target",   comparison = "B-A")
    b <- NADIA:::.get_feature_ids(d, mode = "specific", comparison = "B-A")
    expect_identical(a, b)
    expect_identical(a, c("p1", "p3"))
})

test_that(".get_feature_ids modes select the right features", {
    d <- data.frame(FeatureID = c("p1", "p2", "p3"),
                    sig_any   = c(TRUE, FALSE, TRUE),
                    `adjP_B-A` = c(0.01, 0.20, 0.04),
                    check.names = FALSE)

    expect_identical(NADIA:::.get_feature_ids(d, "all"), c("p1", "p2", "p3"))
    expect_identical(NADIA:::.get_feature_ids(d, "any"), c("p1", "p3"))

    # alpha is respected in the targeted mode
    expect_identical(
        NADIA:::.get_feature_ids(d, "target", alpha = 0.02, comparison = "B-A"),
        "p1")
})

test_that(".get_feature_ids does not let NA FeatureIDs through (regression)", {
    # Regression: the PCA copy filtered with `sig_any == TRUE`, which returns NA
    # rows when sig_any contains NA and so leaked NA FeatureIDs downstream. The
    # which() from the heatmap copy is the one that was kept.
    d <- data.frame(FeatureID = c("p1", "p2", "p3"),
                    sig_any   = c(TRUE, NA, FALSE),
                    stringsAsFactors = FALSE)

    ids <- NADIA:::.get_feature_ids(d, mode = "any")
    expect_identical(ids, "p1")
    expect_false(anyNA(ids))
})

test_that(".get_feature_ids deduplicates by FeatureID", {
    # The long format repeats each feature once per sample; the helper must
    # return one ID per feature, not one per row.
    d <- data.frame(FeatureID = c("p1", "p1", "p2", "p2"),
                    sig_any   = c(TRUE, TRUE, FALSE, FALSE))
    expect_identical(NADIA:::.get_feature_ids(d, "all"), c("p1", "p2"))
})

test_that(".get_feature_ids errors are actionable", {
    d <- data.frame(FeatureID = "p1", sig_any = TRUE)

    expect_error(NADIA:::.get_feature_ids(d, mode = "target"),
                 "comparison", ignore.case = TRUE)
    expect_error(NADIA:::.get_feature_ids(d, "target", comparison = "X-Y"),
                 "adjP_X-Y", fixed = TRUE)
    expect_error(
        NADIA:::.get_feature_ids(data.frame(FeatureID = "p1"), mode = "any"),
        "sig_any", fixed = TRUE)
    expect_error(NADIA:::.get_feature_ids(d, mode = "nonsense"))
})


# --- Proteome Discoverer helpers ---------------------------------------------

test_that(".parse_gene_from_description extracts GN= and returns NA otherwise", {
    f <- NADIA:::.parse_gene_from_description

    x <- c("Serum albumin OS=Homo sapiens OX=9606 GN=ALB PE=1 SV=2",
           "Uncharacterised protein OS=Homo sapiens",
           "Protein X GN=TP53")

    expect_identical(f(x), c("ALB", NA, "TP53"))
    expect_true(is.na(f("no gene here")))
})

test_that(".parse_gene_from_description stops the gene at the first space", {
    # GN= runs to the next space; taking more would swallow the following tag.
    f <- NADIA:::.parse_gene_from_description
    expect_identical(f("x GN=ALB PE=1 SV=2"), "ALB")
})
