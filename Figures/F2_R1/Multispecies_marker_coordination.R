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




gene_sets <- c("GO.0097205.renal.filtration.BP", "GO.0008179.adenylate.cyclase.binding.MF", 
               "GO.0006525.arginine.metabolic.process.BP", "GO.0046835.carbohydrate.phosphorylation.BP", "GO.0005324.long.chain.fatty.acid.transporter.activity.MF",
               "GO.0019825.oxygen.binding.MF", "GO.0089718.amino.acid.import.across.plasma.membrane.BP",
               "GO.0001964.startle.response.BP", 
               "GO.0005005.transmembrane.ephrin.receptor.activity.MF",
               "GO.0006837.serotonin.transport.BP", "GO.0006972.hyperosmotic.response.BP", "GO.0001573.ganglioside.metabolic.process.BP", 
               "GO.0070266.necroptotic.process.BP", "GO.0098719.sodium.ion.import.across.plasma.membrane.BP", "GO.0006833.water.transport.BP", "GO.0034134.toll.like.receptor.2.signaling.pathway.BP" ,
               "GO.0090083.regulation.of.inclusion.body.assembly.BP", "GO.1902017.regulation.of.cilium.assembly.BP",  
               "GO.1905063.regulation.of.vascular.smooth.muscle.cell.differentiation.BP", "GO.0072606.interleukin.8.secretion.BP", "GO.0015301.anion.anion.antiporter.activity.MF", "GO.0015106.bicarbonate.transmembrane.transporter.activity.MF",  
               "GO.0006677.glycosylceramide.metabolic.process.BP", "GO.0033630.positive.regulation.of.cell.adhesion.mediated.by.integrin.BP", "GO.0000030.mannosyltransferase.activity.MF", 
               "GO.0003009.skeletal.muscle.contraction.BP", "GO.0002922.positive.regulation.of.humoral.immune.response.BP")



# Assuming final_data is your data frame ready for plotting
# Add a temporary column to final_data for ordering
final_data$Order <- match(final_data$CellType, cell_types)

# Reorder final_data based on the new Order column
final_data <- final_data[order(final_data$Order), ]

# Remove the temporary Order column
final_data$Order <- NULL

specific_gene_sets_data <- final_data[, gene_sets]



# Assume 'df' is your original dataframe

df <- final_data
# Create a vector of cell types
unique_cell_types <- unique(df$CellType)

# Create a vector of column names for gene sets (excluding 'CellType')
gene_set_columns <- setdiff(names(df), "CellType")

# Initialize a new dataframe for the averages
average_df <- data.frame(matrix(ncol = length(gene_set_columns), nrow = length(unique_cell_types)))
names(average_df) <- gene_set_columns
rownames(average_df) <- unique_cell_types



# Calculate the averages
for (gene_set in gene_set_columns) {
  for (cell_type in unique_cell_types) {
    subset_data <- df[df$CellType == cell_type, gene_set]
    average_df[cell_type, gene_set] <- mean(subset_data, na.rm = TRUE)
  }
}


average_df$CellType <- rownames(average_df)
# Now average_df contains the average of each gene set for each cell type
average_df$Order <- match(average_df$CellType, cell_types)

# Reorder final_data based on the new Order column
average_df <- average_df[order(average_df$Order), ]

# Remove the temporary Order column
average_df$Order <- NULL
average_df$CellType <- NULL

specific_gene_sets_data <- average_df[, gene_sets]

# Transpose the data
transposed_data <- t(specific_gene_sets_data)

# Plot the heatmap with transposed data
pheatmap(transposed_data, scale = "none", clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", clustering_method = "complete",
         show_rownames = TRUE, show_colnames = TRUE, cluster_rows = FALSE, cluster_cols = FALSE,
         fontsize = 10)

