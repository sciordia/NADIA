# =============================================================================
# Documentation of the data shipped with the package
# =============================================================================

#' Example DIA experiment, already preprocessed
#'
#' Output of [preprocess_spectronaut()] on a trimmed Spectronaut report, ready to
#' feed [process_proteomics()] and the rest of the pipeline without having to
#' read the source file again.
#'
#' @format An S3 object of class `c("spectronaut_data", "proteomics_data",
#'   "list")` with three elements:
#'   \describe{
#'     \item{`metadata`}{Data frame of 12 rows and 7 columns, one row per
#'       injection: `R.FileName`, `R.Condition`, `R.Replicate`, `Coding` and the
#'       three identification counts reported by Spectronaut
#'       (`R.PrecursorsIdentified`, `R.StrippedSequencesIdentified`,
#'       `R.ProteinGroupsIdentified`).}
#'     \item{`protein_id`}{Data frame of 2,000 rows by 52 columns holding the
#'       identification information of every protein group.}
#'     \item{`protein_quant`}{Data frame of 2,000 rows by 44 columns: eight
#'       annotation columns plus, for each of the 12 samples, the number of
#'       precursors and of stripped sequences used for quantification and the
#'       intensity (`PG.Quantity_<condition>_<replicate>`).}
#'   }
#'
#' @details
#' The data come from a quantitative proteomics experiment acquired on an
#' Orbitrap Exploris at the Proteomics Facility of the Centro Nacional de
#' Biotecnologia (CNB-CSIC) and processed with Spectronaut v20. The original
#' experiment has 16 injections of the same digest split across four conditions
#' (A-D) of four replicates each, and 10,437 protein groups.
#'
#' So that the package examples run quickly, three conditions (A, B and D) and a
#' random sample of 2,000 protein groups are kept. Three and not two conditions
#' are retained because the soft-clustering in [pattern_profiler_analysis()]
#' groups profiles across conditions, and with only two points per profile there
#' is no pattern left to cluster.
#'
#' The sample is drawn at random, with a fixed seed, rather than picking the
#' best-covered proteins. That matters: **8.1 % of the intensities are missing
#' and 72.2 % of the proteins are complete**, closely matching the 7.9 % and
#' 72.1 % of the full experiment. Selecting the most-quantified proteins instead
#' would yield a dataset with no missing values at all, in which every imputation
#' method returns the input unchanged, and in which 87 % of the proteins come out
#' differentially expressed instead of the 31 % of the full data.
#'
#' The full recipe is in `system.file("scripts", "make_extdata.R", package =
#' "NADIA")`. The trimmed report this object is derived from is also
#' distributed as a file, in `system.file("extdata", "nadia_dia_report.tsv.gz",
#' package = "NADIA")`, so that [preprocess_spectronaut()] can be run end to
#' end.
#'
#' @source Proteomics Facility, Centro Nacional de Biotecnologia (CNB-CSIC).
#'
#' @examples
#' data(nadia_dia)
#' nadia_dia
#'
#' # The three conditions and their four replicates
#' table(nadia_dia$metadata$R.Condition)
#'
#' # Proteins and intensity columns
#' dim(nadia_dia$protein_quant)
#' grep("^PG.Quantity_", colnames(nadia_dia$protein_quant), value = TRUE)
#'
#' # Missing values, which is what the package is about. A Spectronaut report is
#' # long, one row per run and protein group, and where a protein was not
#' # quantified the row is simply absent; the pivot to wide format is what turns
#' # those absences into NA.
#' quant <- as.matrix(nadia_dia$protein_quant[
#'   grep("^PG.Quantity_", colnames(nadia_dia$protein_quant))])
#' round(100 * mean(is.na(quant)), 1)
"nadia_dia"
