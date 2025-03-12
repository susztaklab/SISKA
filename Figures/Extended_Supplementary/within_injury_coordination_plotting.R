library(ggplot2)
library(dplyr)
library(tidyr)

#we reload the data

# Read the data
filename <- paste0(".../Atlas/scSPECTRA/healthy_vs_injured_epithelial_coord_human/", "healthy_TAL_MD", "_r2_results.csv")
cell_data_healthy1 <- read.csv(filename, stringsAsFactors = FALSE)
cell_data_healthy1$condition <- "healthy"

# Read the data
filename <- paste0(".../Atlas/scSPECTRA/healthy_vs_injured_epithelial_coord_human/", "healthy_PT", "_r2_results.csv")
cell_data_healthy2 <- read.csv(filename, stringsAsFactors = FALSE)
cell_data_healthy2$condition <- "healthy"

# Read the data
filename <- paste0(".../Atlas/scSPECTRA/healthy_vs_injured_epithelial_coord_human/", "healthy_DCT_CNT_CD", "_r2_results.csv")
cell_data_healthy3 <- read.csv(filename, stringsAsFactors = FALSE)
cell_data_healthy3$condition <- "healthy"

# Read the data
filename <- paste0(".../Atlas/scSPECTRA/healthy_vs_injured_epithelial_coord_human/", "injured", "_r2_results.csv")
cell_data_injured <- read.csv(filename, stringsAsFactors = FALSE)
cell_data_injured$condition <- "injured"

# Combine the healthy datasets
cell_data_healthy_combined <- rbind(cell_data_healthy1, cell_data_healthy2, cell_data_healthy3)


all_data <- list()

all_data[["healthy"]] <- cell_data_healthy_combined
all_data[["injured"]] <- cell_data_injured

combined_data <- do.call(rbind, all_data)


# Pivot the data to have gene sets as columns and samples as rows

final_data <- spread(combined_data, key = GeneSet, value = R2Value)

# Set row names and remove the Sample column
rownames(final_data) <- final_data$Sample
final_data$Sample <- NULL

cell_type_of_interest <- "injured"

# Create a logical vector where TRUE represents rows with the cell type of interest
is_cell_type_of_interest <- final_data$condition == cell_type_of_interest

final_data$condition <- NULL

# Data for the cell type of interest
data_of_interest <- final_data[is_cell_type_of_interest, ]

# Data for all other cell types
data_other <- final_data[!is_cell_type_of_interest, ]

# Initialize a dataframe to store the results
results <- data.frame(GeneSet=colnames(final_data), MeanInGroup=NA, MeanOutOfGroup=NA, LogFoldChange=NA, PValue=NA)



# Perform Wilcoxon tests and calculate statistics for each gene set
for (gene_set in colnames(final_data)) {
  in_group <- data_of_interest[, gene_set]
  out_group <- data_other[, gene_set]
  
  # Exclude NaN values for the analysis
  in_group <- in_group[!is.nan(in_group)]
  out_group <- out_group[!is.nan(out_group)]
  
  # Perform Wilcoxon test, checking if both groups have more than one unique value
  if (length(unique(in_group)) > 1 && length(unique(out_group)) > 1) {
    test_result <- wilcox.test(in_group, out_group, alternative = "two.sided")
    p_value <- test_result$p.value
  } else {
    # Assign NA to p_value if one of the groups has one or no unique values
    p_value <- NA
  }
  
  log_fold_change <- log2(mean(in_group, na.rm = TRUE) + 0.01) - log2(mean(out_group, na.rm = TRUE) + 0.01)
  
  # Store the results
  results[results$GeneSet == gene_set, "MeanInGroup"] <- mean(in_group, na.rm = TRUE)
  results[results$GeneSet == gene_set, "MeanOutOfGroup"] <- mean(out_group, na.rm = TRUE)
  results[results$GeneSet == gene_set, "LogFoldChange"] <- log_fold_change
  results[results$GeneSet == gene_set, "PValue"] <- p_value
}

# Adjust P-values for multiple testing
results$FDR <- p.adjust(results$PValue, method = "fdr")

#filtered_results <- subset(results, LogFoldChange > 0)



gene_set <- "..."  

# Read the data
filename <- paste0(".../Atlas/scSPECTRA/healthy_vs_injured_epithelial_coord_human/", "healthy_TAL_MD", "_r2_results.csv")
cell_data_healthy1 <- read.csv(filename, stringsAsFactors = FALSE) %>% filter(GeneSet == gene_set)
cell_data_healthy1$condition <- "healthy"

# Read the data
filename <- paste0(".../Atlas/scSPECTRA/healthy_vs_injured_epithelial_coord_human/", "healthy_PT", "_r2_results.csv")
cell_data_healthy2 <- read.csv(filename, stringsAsFactors = FALSE) %>% filter(GeneSet == gene_set)
cell_data_healthy2$condition <- "healthy"

# Read the data
filename <- paste0(".../Atlas/scSPECTRA/healthy_vs_injured_epithelial_coord_human/", "healthy_DCT_CNT_CD", "_r2_results.csv")
cell_data_healthy3 <- read.csv(filename, stringsAsFactors = FALSE) %>% filter(GeneSet == gene_set)
cell_data_healthy3$condition <- "healthy"

# Read the data
filename <- paste0(".../Atlas/scSPECTRA/healthy_vs_injured_epithelial_coord_human/", "injured", "_r2_results.csv")
cell_data_injured <- read.csv(filename, stringsAsFactors = FALSE) %>% filter(GeneSet == gene_set)
cell_data_injured$condition <- "injured"

# Combine the healthy datasets
cell_data_healthy_combined <- rbind(cell_data_healthy1, cell_data_healthy2, cell_data_healthy3)

all_data <- list()

all_data[["injured"]] <- cell_data_injured
all_data[["healthy"]] <- cell_data_healthy_combined

combined_data <- do.call(rbind, all_data)

# Reorder factor levels
combined_data$condition <- factor(combined_data$condition, levels = c("injured", "healthy"))

# Create the density plot
ggplot(combined_data, aes(x = R2Value, fill = condition)) +
  geom_density(alpha = 0.5, adjust = 1) +
  scale_fill_manual(values = c("injured" = "red", "healthy" = "blue")) +
  scale_x_continuous(limits = c(0, 1)) +  # Set x-axis limits
  labs(title = gene_set, x = "", y = "") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(colour = "black", size = 12),
    axis.ticks = element_line(),
    legend.position = "none",  # Remove legend
    axis.line = element_line(colour = "black")
  )

