
# MSstats Pipeline ("noBIG")
raw <- read.csv("./data/20260306_CursoProtQ_2024_DIA_Exploris_v20p5_sinImp_sinNorm_sinD_Report.tsv", sep = "\t")
annotation <-  read.csv("./data/MSstats_Annotation_sinD.csv")

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
  -1,  1,  0,  0,   # B - A
  -1,  0,  1,  0,   # C - A
  -1,  0,  0,  1,   # D - A
  0, -1,  1,  0,   # C - B
  0, -1,  0,  1,   # D - B
  0,  0,  -1, 1    # D - C
), ncol = 4, byrow = TRUE)
row.names(comparison) <- c("B-A", "C-A", "D-A", "C-B", "D-B", "D-C")
colnames(comparison) <- c("A", "B", "C", "D")

comparison_result <- groupComparison(contrast.matrix = comparison,
                                     data = processed.quant
                                     )

# Save DF
DEPs_results <- comparison_result$ComparisonResult
readr::write_tsv(DEPs_results, file = "./results/Q24_DIA_Spectronaut_20p5_MSstats/DEPs_results_MSstats.tsv")
