# =============================================================================
# Locating the agent skill shipped with the package
# =============================================================================

#' Locate the NADIA agent skill
#'
#' Returns the directory holding the agent skill that ships with NADIA. The
#' skill documents the analysis pipeline for a coding assistant: which function
#' to call at each stage, how to choose between the normalisation and imputation
#' methods, what to verify afterwards, and the mistakes the package's defaults
#' invite. It is plain Markdown, so any assistant can read it, and it follows
#' the layout Claude expects for a skill: a `SKILL.md` entry point that points
#' at the files under `reference/` and `scripts/`.
#'
#' To use it with Claude Code, copy the directory into the skills folder of a
#' project (`.claude/skills/nadia/`) or of the user (`~/.claude/skills/nadia/`).
#' The examples below print the command that does it.
#'
#' @param ... Character path segments below the skill directory, passed on to
#'   [base::system.file()]. With no arguments the directory itself is returned.
#'
#' @return A single character string: the absolute path, or `""` when the skill
#'   cannot be found, which is what `system.file()` returns for a missing file.
#'
#' @seealso `vignette("NADIA")` for the same pipeline written for a person.
#'
#' @examples
#' # Where the skill lives
#' nadia_skill_path()
#'
#' # Its entry point, and the files it refers to
#' nadia_skill_path("SKILL.md")
#' list.files(nadia_skill_path(), recursive = TRUE)
#'
#' # Install it for Claude Code, for the current user
#' cat(sprintf(
#'     "cp -R '%s' ~/.claude/skills/nadia",
#'     nadia_skill_path()
#' ))
#' @export
nadia_skill_path <- function(...) {
    system.file("skill", ..., package = "NADIA")
}
