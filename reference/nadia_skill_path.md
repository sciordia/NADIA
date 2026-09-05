# Locate the NADIA agent skill

Returns the directory holding the agent skill that ships with NADIA. The
skill documents the analysis pipeline for a coding assistant: which
function to call at each stage, how to choose between the normalisation
and imputation methods, what to verify afterwards, and the mistakes the
package's defaults invite. It is plain Markdown, so any assistant can
read it, and it follows the layout Claude expects for a skill: a
`SKILL.md` entry point that points at the files under `reference/` and
`scripts/`.

## Usage

``` r
nadia_skill_path(...)
```

## Arguments

- ...:

  Character path segments below the skill directory, passed on to
  [`base::system.file()`](https://rdrr.io/r/base/system.file.html). With
  no arguments the directory itself is returned.

## Value

A single character string: the absolute path, or `""` when the skill
cannot be found, which is what
[`system.file()`](https://rdrr.io/r/base/system.file.html) returns for a
missing file.

## Details

To use it with Claude Code, copy the directory into the skills folder of
a project (`.claude/skills/nadia/`) or of the user
(`~/.claude/skills/nadia/`). The examples below print the command that
does it.

## See also

[`vignette("NADIA")`](https://sciordia.github.io/NADIA/articles/NADIA.md)
for the same pipeline written for a person.

## Examples

``` r
# Where the skill lives
nadia_skill_path()
#> [1] "/home/runner/work/_temp/Library/NADIA/skill"

# Its entry point, and the files it refers to
nadia_skill_path("SKILL.md")
#> [1] "/home/runner/work/_temp/Library/NADIA/skill/SKILL.md"
list.files(nadia_skill_path(), recursive = TRUE)
#> [1] "SKILL.md"                   "reference/export.md"       
#> [3] "reference/import.md"        "reference/method-choice.md"
#> [5] "reference/processing.md"    "reference/validation.md"   
#> [7] "reference/visualisation.md" "scripts/run_pipeline.R"    

# Install it for Claude Code, for the current user
cat(sprintf(
    "cp -R '%s' ~/.claude/skills/nadia",
    nadia_skill_path()
))
#> cp -R '/home/runner/work/_temp/Library/NADIA/skill' ~/.claude/skills/nadia
```
