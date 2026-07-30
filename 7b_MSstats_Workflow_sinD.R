
# MSstats Pipeline ("noBIG")
raw <- read.csv("./data-raw/20260306_CursoProtQ_2024_DIA_Exploris_v20p5_sinImp_sinNorm_sinD_MSstats.tsv", sep = "\t")
annotation <-  read.csv("./data-raw/MSstats_Annotation_sinD.csv")

quant <-SpectronauttoMSstatsFormat(raw,
                                   annotation = annotation,
                                   intensity = "PeakArea",
                                   filter_with_Qvalue = TRUE,
                                   qvalue_cutoff = 0.01,
                                   useUniquePeptide = TRUE,
                                   removeFewMeasurements = TRUE,
                                   removeProtein_with1Feature = FALSE,
                                   summaryforMultipleRows = max
                                   )



processed.quant <- dataProcess(quant, normalization = "equalizeMedians")

# Define the comparison matrix
comparison <- matrix(c(
  -1,  1,  0,   # B - A
  -1,  0,  1,   # C - A
  0, -1,  1   # C - B
), ncol = 3, byrow = TRUE)
row.names(comparison) <- c("B-A", "C-A", "C-B")
colnames(comparison) <- c("A", "B", "C")

comparison_result <- groupComparison(contrast.matrix = comparison,
                                     data = processed.quant
                                     )

# Save DF
DEPs_results <- comparison_result$ComparisonResult
readr::write_tsv(DEPs_results, file = "./results/Q24_DIA_Spectronaut_20p5_MSstats_sinD//DEPs_results_MSstats_sinD.tsv")
