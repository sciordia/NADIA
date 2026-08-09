# =============================================================================
# Provenance and preparation of inst/extdata/ and data/nadia_dia.rda
# =============================================================================
#
# Bioconductor requires documenting how the data shipped with a package were
# obtained and prepared. This script is the complete, reproducible recipe.
#
# The source files are NOT distributed, and are not part of the repository
# either: they are full reports of between 0.6 and 29 MB from experiments that
# are still unpublished, held by the maintainer. This script documents exactly
# how they were reduced to what the package does ship, so the provenance is on
# record even though the inputs cannot be handed out.
#
# To re-run it, place those reports in a directory and point RAW at it (below);
# the default assumes a `data-raw/` beside the package root. Anyone wanting the
# originals should contact the maintainer.
#
# Origin of the data
# ------------------
# All of it comes from quantitative proteomics experiments acquired at the
# Proteomics Facility of the Centro Nacional de Biotecnologia (CNB-CSIC).
#
# * DIA (Spectronaut): 16 injections of the same digest, split across four
#   conditions (A-D) of four replicates each, acquired on an Orbitrap Exploris
#   and processed with Spectronaut v20. Report at protein-group and run level
#   (long format).
# * DIA spike-in: the same report annotated with the species of every protein
#   group, for the benchmarking based on spiked-in organisms.
# * TMT (Proteome Discoverer): a TMTpro experiment of two mixes and ten
#   fractions, searched with three engines.
# * LFQ (Proteome Discoverer): a label-free experiment, with its sample ->
#   condition annotation file.
#
# Trimming applied
# ----------------
# Every report is cut down to 2,000 protein groups, drawn at random with a fixed
# seed. Two reasons, in this order:
#
# 1. The source experiments are unpublished. Shipping a report whole would
#    distribute the entire private dataset through a public package, which is not
#    what an example needs. A random subset illustrates the format and the
#    pipeline just as well.
# 2. It keeps `R CMD check` well under the 10-minute limit with more than a
#    hundred examples running.
#
# The species map is likewise restricted to the protein groups that survive in
# the DIA example, since nothing downstream looks up the rest.
#
# Three of the four conditions are kept (A, B and D) because with only two the
# soft clustering in the Pattern Profiler becomes degenerate -- it groups profiles
# across conditions, and with two points per profile there is no pattern left to
# group.
#
# The 2,000 proteins are drawn as a RANDOM SAMPLE with a fixed seed, and that
# choice matters. An earlier version took "the 2,000 proteins with the most
# quantified values", which sounds harmless but selects precisely the complete,
# high-abundance ones and destroys the two structures the package exists to deal
# with. Measured on this dataset:
#
#                        missing   Up     Down   No Change
#   full (10,437)          7.9 %   15.5 %  15.3 %   69.2 %
#   most-quantified        0.0 %   35.6 %  51.1 %   13.3 %   <- unrepresentative
#   random sample          8.1 %   16.4 %  15.6 %   67.9 %   <- faithful
#
# With 0 % missing values every imputation method returns the input matrix
# unchanged, so the imputation examples would demonstrate nothing -- in a package
# named after missing values. A random sample reproduces both the missingness rate
# and the differential-expression balance at the same size and speed.
# =============================================================================

CONDITIONS <- c("A", "B", "D")
N_PROTEINS <- 2000L
SEED       <- 42L

RAW      <- Sys.getenv("NADIA_RAW", "data-raw")   # see the note at the top
EXTDATA  <- "inst/extdata"
dir.create(EXTDATA, recursive = TRUE, showWarnings = FALSE)
dir.create("data",  recursive = TRUE, showWarnings = FALSE)

gz <- function(df, target) {
  con <- gzfile(target, "wt", compression = 9)
  on.exit(close(con), add = TRUE)
  utils::write.table(df, con, sep = "\t", row.names = FALSE, quote = FALSE,
                     na = "")
  message(sprintf("  %-32s %6.2f MB", basename(target),
                  file.size(target) / 1024^2))
}

# --- 1. Spectronaut DIA report ----------------------------------------------
message("1. nadia_dia_report.tsv.gz")
dia <- utils::read.delim(file.path(RAW, "Curso_Q24_DIA_Spectronaut_v20_Report.tsv"),
                         sep = "\t", header = TRUE, check.names = FALSE,
                         stringsAsFactors = FALSE)
dia <- dia[dia$R.Condition %in% CONDITIONS, , drop = FALSE]

# A random sample of protein groups, so that the missingness rate and the
# differential-expression balance of the full dataset are preserved (see the
# header for the numbers). The seed is fixed so the file is reproducible.
set.seed(SEED)
selected <- sample(unique(dia$PG.ProteinGroups), N_PROTEINS)
dia <- dia[dia$PG.ProteinGroups %in% selected, , drop = FALSE]
gz(dia, file.path(EXTDATA, "nadia_dia_report.tsv.gz"))

# --- 2. DIA report with species (spike-in, for benchmarking) -----------------
# Restricted to the protein groups kept above: the benchmarking examples only
# ever look up the 2,000 that are in the example dataset, so shipping the
# species of the whole experiment would expose data for no purpose.
message("2. nadia_dia_spikein.tsv.gz")
spk <- utils::read.delim(
  file.path(RAW, "Curso_Q24_DIA_Spectronaut_v20_Report_Species.tsv"),
  sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
spk <- spk[spk$PG.ProteinGroups %in% selected, , drop = FALSE]
gz(spk, file.path(EXTDATA, "nadia_dia_spikein.tsv.gz"))

# --- 3. Proteome Discoverer TMT report --------------------------------------
# TMT and LFQ come in wide format, one row per protein, so the sample is drawn
# straight from the rows. Same criterion as the DIA report: random, fixed seed.
message("3. nadia_tmt_report.tsv.gz")
tmt <- utils::read.delim(
  file.path(RAW, "20260527_Q25_TMTpro_TMT1y2_10Fr_Static_3engines_onlyRAW.tsv"),
  sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
set.seed(SEED)
tmt <- tmt[sort(sample(nrow(tmt), min(N_PROTEINS, nrow(tmt)))), , drop = FALSE]
gz(tmt, file.path(EXTDATA, "nadia_tmt_report.tsv.gz"))

# --- 4. Proteome Discoverer LFQ report and its annotation -------------------
message("4. nadia_lfq_report.tsv.gz + nadia_lfq_annotation.tsv")
lfq <- utils::read.delim(
  file.path(RAW, "20260710_AIturrate_2659_LFQ_QUANT_onlyRAW.tsv"),
  sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
set.seed(SEED)
lfq <- lfq[sort(sample(nrow(lfq), min(N_PROTEINS, nrow(lfq)))), , drop = FALSE]
gz(lfq, file.path(EXTDATA, "nadia_lfq_report.tsv.gz"))

file.copy(file.path(RAW, "20260710_AIturrate_2659_LFQ_QUANT_onlyRAW_Annot.tsv"),
          file.path(EXTDATA, "nadia_lfq_annotation.tsv"), overwrite = TRUE)

# --- 5. DIA-NN protein-group matrix -----------------------------------------
# The same 16 runs as the Spectronaut report in step 1, searched with DIA-NN
# instead. Wide format, so the sample is drawn from the rows as in steps 3-4.
#
# DIA-NN names each intensity column after the raw file, which says nothing
# about the design. The columns are renamed here to the
# "Abundance: <Condition>_<Replicate>" convention that preprocess_diann() reads,
# so that the example needs no annotation file and no extra argument. Reading
# the report exactly as DIA-NN wrote it is the `sample_names` route instead.
message("5. nadia_diann_report.tsv.gz")
dnn <- utils::read.delim(
  file.path(RAW, "20241014_CursoProtQ_DIA_DIANN_report.pg_matrix.tsv"),
  sep = "\t", header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)

annot_cols <- c("Protein.Group", "Protein.Ids", "Protein.Names", "Genes",
                "First.Protein.Description", "N.Sequences",
                "N.Proteotypic.Sequences")
run_cols <- setdiff(names(dnn), annot_cols)

# ...sample_A1.raw -> A_1
tag <- sub("^.*sample_([A-Za-z]+)([0-9]+)\\.raw$", "\\1_\\2", basename(run_cols))
stopifnot(!any(tag == basename(run_cols)))       # every column parsed
names(dnn)[match(run_cols, names(dnn))] <- paste0("Abundance: ", tag)

# The sample is drawn before dropping the unused conditions, exactly as in
# step 1, so that the 2,000 protein groups do not depend on which conditions
# are kept.
set.seed(SEED)
dnn <- dnn[sort(sample(nrow(dnn), min(N_PROTEINS, nrow(dnn)))), , drop = FALSE]
dnn <- dnn[, c(intersect(annot_cols, names(dnn)),
               paste0("Abundance: ", tag[sub("_.*$", "", tag) %in% CONDITIONS])),
           drop = FALSE]
gz(dnn, file.path(EXTDATA, "nadia_diann_report.tsv.gz"))

dnn_mat <- as.matrix(dnn[, grep("^Abundance: ", names(dnn))])
message(sprintf("  missing values: %.1f %% | complete proteins: %.1f %%",
                100 * mean(is.na(dnn_mat)),
                100 * mean(rowSums(is.na(dnn_mat)) == 0)))

# --- 6. Preprocessed object data/nadia_dia.rda ------------------------------
# Built from the already trimmed report in inst/extdata, so that the object and
# the example file describe exactly the same experiment.
message("6. data/nadia_dia.rda")
nadia_dia <- NADIA::preprocess_spectronaut(
  file_path       = file.path(EXTDATA, "nadia_dia_report.tsv.gz"),
  condition_order = CONDITIONS,
  verbose         = FALSE)

save(nadia_dia, file = "data/nadia_dia.rda",
     compress = "xz", compression_level = 9)
message(sprintf("  %-32s %6.2f MB", "nadia_dia.rda",
                file.size("data/nadia_dia.rda") / 1024^2))
message(sprintf("  protein_quant: %d x %d | %d samples | %d conditions",
                nrow(nadia_dia$protein_quant), ncol(nadia_dia$protein_quant),
                nrow(nadia_dia$metadata),
                length(unique(nadia_dia$metadata$R.Condition))))

# Report the missingness rate: it is the property that makes this dataset useful
# as an example, so a regression in the sampling should be visible here.
quant_cols <- grep("^PG.Quantity_", colnames(nadia_dia$protein_quant), value = TRUE)
mat <- as.matrix(nadia_dia$protein_quant[, quant_cols])
mat[mat == 0] <- NA
message(sprintf("  missing values: %.1f %% | complete proteins: %.1f %%",
                100 * mean(is.na(mat)),
                100 * mean(rowSums(is.na(mat)) == 0)))
