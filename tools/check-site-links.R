# =============================================================================
# Every local src= and href= in docs/ must resolve
# =============================================================================
#
# Run after pkgdown::build_site(). Two things it catches that the build itself
# does not:
#
#   * the vignette diagrams. They are referenced with knitr::include_graphics(),
#     which rmarkdown::find_external_resources() reports as neither `web` nor
#     `explicit` -- the pair pkgdown filters on -- so copy_article_images()
#     skips them and every diagram 404s.
#   * the <a href> around the README diagram. pkgdown rewrites the <img src>
#     inside it but leaves the link pointing at the repository path.
#
# Both are fixed by copying vignettes/figures into place; this asserts it worked.

files <- list.files("docs", pattern = "[.]html$", recursive = TRUE,
                    full.names = TRUE)
if (!length(files)) stop("No built site found under docs/.", call. = FALSE)

broken <- list()
n_checked <- 0L

for (f in files) {
    html <- paste(readLines(f, warn = FALSE), collapse = "\n")
    refs <- unlist(regmatches(
        html, gregexpr('(?:src|href)="[^"]+"', html)))
    refs <- sub('^(?:src|href)="', "", sub('"$', "", refs))

    # Absolute, protocol-relative, inline and in-page targets are not files.
    refs <- refs[!grepl("^(https?:|//|data:|mailto:|#|\\?)", refs)]
    refs <- sub("[#?].*$", "", refs)
    refs <- refs[nzchar(refs)]
    if (!length(refs)) next

    n_checked <- n_checked + length(refs)
    target <- file.path(dirname(f), refs)
    missing <- refs[!file.exists(target) & !dir.exists(target)]
    if (length(missing)) {
        broken[[f]] <- unique(missing)
    }
}

message(sprintf("Checked %d local references across %d pages.",
                n_checked, length(files)))

if (length(broken)) {
    for (f in names(broken)) {
        message("  ", sub("^docs/", "", f), ":")
        for (m in broken[[f]]) message("    -> ", m)
    }
    stop(sprintf("%d page(s) reference files that are not there.",
                 length(broken)), call. = FALSE)
}

message("No broken local references.")
