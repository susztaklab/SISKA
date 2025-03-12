library(dplyr)
library(stats)
library(ggplot2)


# Define the paths to the subfolders for healthy and diseased data
healthy_folder_path <- ".../scSPECTRA/celltype_coodination_species/healthy/PT_r2_results.csv"
diseased_folder_path <- ".../scSPECTRA/celltype_coodination_species/diseased/1B/injPT_r2_results.csv"

# Construct file paths
healthy_file_path <- paste0(healthy_folder_path)
diseased_file_path <- paste0(diseased_folder_path)

# Load R² results for healthy and diseased
r2_healthy <- read.csv(healthy_file_path)
r2_diseased <- read.csv(diseased_file_path)

# Preparing for analysis
results <- data.frame(GeneSet = character(), MeanHealthy = numeric(), MeanDiseased = numeric(), FoldChange = numeric(), PValue = numeric(), FDR = numeric())

# Get unique gene sets
gene_sets <- unique(c(r2_healthy$GeneSet, r2_diseased$GeneSet))

# Loop through each gene set to perform Wilcoxon test
for (gene_set in gene_sets) {
  healthy_values <- r2_healthy$R2Value[r2_healthy$GeneSet == gene_set]
  diseased_values <- r2_diseased$R2Value[r2_diseased$GeneSet == gene_set]
  
  # Wilcoxon test (two sided - both information relevant)
  test_result <- wilcox.test(healthy_values, diseased_values, alternative = "two.sided")
  
  # Calculate mean, fold change, and store results
  mean_healthy <- mean(healthy_values, na.rm = TRUE)
  mean_diseased <- mean(diseased_values, na.rm = TRUE)
  fold_change <- mean_diseased / mean_healthy
  
  results <- rbind(results, data.frame(GeneSet = gene_set, MeanHealthy = mean_healthy, MeanDiseased = mean_diseased, FoldChange = fold_change, PValue = test_result$p.value))
}

# Adjusting for multiple testing (FDR)
results$FDR <- p.adjust(results$PValue, method = "fdr")

# higher coordination in diseased cells
filtered_results <- results[results$FoldChange > 1.0 & results$FDR < 0.05, ]


#plotting the examples

gene_set <- "GO.0070102.interleukin.6.mediated.signaling.pathway.BP"  

# Load and filter R² results for healthy and diseased
r2_healthy <- read.csv(healthy_file_path) %>% filter(GeneSet == gene_set)
r2_diseased <- read.csv(diseased_file_path) %>% filter(GeneSet == gene_set)

# Label the data
r2_healthy$Condition <- "Healthy"
r2_diseased$Condition <- "Diseased"

# Combine the data
combined_data <- rbind(r2_healthy, r2_diseased)

# Create the density plot
ggplot(combined_data, aes(x = R2Value, fill = Condition)) +
  geom_density(alpha = 0.5, adjust = 1) +
  scale_fill_manual(values = c("Healthy" = "blue", "Diseased" = "red")) +
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
