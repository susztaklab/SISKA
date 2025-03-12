library(dplyr)
library(tidyr)
library(pheatmap)

#cell types to compare

cell_types <- c("Podo", "PEC", "PTS1", "PTS2", "PTS3", "injPT", "DTL_ATL", "MD", "TAL", "injTAL", "DCT", "DCT2", "CNT", "CD_PC", "injDCT_CNT", "ICA", "ICB",  "EC", "Stromal", "Immune")

all_data <- list()

for (cell_type in cell_types) {
  filename <- paste0(".../scSPECTRA/multispecies_marker_new/", cell_type, "_r2_results.csv")
  cell_data <- read.csv(filename, stringsAsFactors = FALSE)
  cell_data$Sample <- paste(cell_type, cell_data$Sample, sep = "_")

  # Add the cell type column
  cell_data$CellType <- cell_type

  all_data[[cell_type]] <- cell_data
}

combined_data <- do.call(rbind, all_data)

# Add CellType as a column before spreading
combined_data <- combined_data %>% mutate(CellType = factor(CellType, levels = cell_types))

final_data <- spread(combined_data, key = GeneSet, value = R2Value)
rownames(final_data) <- final_data$Sample
final_data$Sample <- NULL



##########we save marker statistic for all cell types#########

# Loop through each cell type
for (cell_type_of_interest in cell_types) {

  # Create a logical vector for the cell type of interest
  is_cell_type_of_interest <- grepl(paste0("^", cell_type_of_interest, "_"), rownames(final_data))

  # Subset the data
  data_of_interest <- final_data[is_cell_type_of_interest, ]
  data_other <- final_data[!is_cell_type_of_interest, ]

  # Initialize results dataframe
  results <- data.frame(GeneSet=colnames(final_data), MeanInGroup=NA, MeanOutOfGroup=NA, LogFoldChange=NA, PValue=NA)

  # Analysis for each gene set
  for (gene_set in colnames(final_data)) {
    in_group <- data_of_interest[, gene_set]
    out_group <- data_other[, gene_set]

    # Exclude NaN values
    in_group <- in_group[!is.nan(in_group)]
    out_group <- out_group[!is.nan(out_group)]

    # Wilcoxon test and calculations
    if (length(unique(in_group)) > 1 && length(unique(out_group)) > 1) {
      test_result <- wilcox.test(in_group, out_group, alternative = "greater")
      p_value <- test_result$p.value
    } else {
      p_value <- NA
    }

    log_fold_change <- log2(mean(in_group, na.rm = TRUE) + 0.01) - log2(mean(out_group, na.rm = TRUE) + 0.01)

    # Store results
    results[results$GeneSet == gene_set, ] <- c(gene_set, mean(in_group, na.rm = TRUE), mean(out_group, na.rm = TRUE), log_fold_change, p_value)
  }

  # Adjust P-values
  results$FDR <- p.adjust(results$PValue, method = "fdr")

  # Save the results for each cell type
  write.csv(results, paste0(".../scSPECTRA/multispecies_marker_new/statistic/", cell_type_of_interest, "_marker_features.csv"), row.names = FALSE)
}

