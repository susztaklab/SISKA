library(ggplot2)
library(dplyr)

#working
# Read in the data for this sample
datH <- readRDS(file = paste0(sample_output_folder, "datH.rds"))
datD <- readRDS(file = paste0(sample_output_folder, "datD.rds"))

gene_neg <- "AKR1C1"
expression_H_neg <- datH[, gene_neg]
expression_D_neg <- datD[, gene_neg]

# Calculate density
density_H_neg <- density(expression_H_neg)

# Set x-axis and y-axis to start from 0, and add a slight margin to max values
xlim_neg <- c(0, max(density_H_neg$x, expression_D_neg) * 1.05)
ylim_neg <- c(0, max(density_H_neg$y) * 1.05)

# Plot the density with thicker blue line, ensuring axes meet at 0 precisely
plot(density_H_neg, main = paste(gene_neg), 
     xlab = "", ylab = "", col = "blue",
     xlim = xlim_neg, ylim = ylim_neg, lwd = 4, cex.axis = 2,
     xaxs = "i", yaxs = "i")

# Add red dashed line and points
abline(v = expression_D_neg, col = "red", lwd = 4, lty = 2)
points(expression_D_neg, 0, col = "red", pch = 19)




metadata = read.csv("/home/isilon/users/o_kloetzer/Atlas/scSPECTRA/R2_pval/Atlas_Extended_II_Albuminuria_gpt.csv")
metadata_subset = subset(metadata, subset = Disease_level1 != "Control" & Sex == "male")
diseased_samples = as.vector(metadata_subset$Sample)
diseased_samples

metadata = read.csv("/home/isilon/users/o_kloetzer/Atlas/scSPECTRA/R2_pval/Atlas_Extended_II_Albuminuria_gpt.csv")
metadata_subset = subset(metadata, subset = Disease_level1 != "Control" & Hypertension == "Yes")
diseased_samples = as.vector(metadata_subset$Sample)

diseased_samples = c('28-12510', '29-10010', '29-10393', '29-10400', '29-10406', '30-10018',
                     '31-10006', '31-10091', '31-10105', 'diabetic_3_humphreys_DKD',
                     'HK2558', 'HK2596', 'HK2851', 'KRAD70', 'PKD6_humphreys_ADPKD',
                     'PKD8_humphreys_ADPKD')

metadata_subset = subset(metadata, subset = Disease_level1 != "Control")
#diseased_samples = as.vector(metadata_subset$Sample)


gene_neg <- "DLG4"

# Initialize vectors to collect expression values across samples
all_expression_H_neg <- c()
all_expression_D_neg <- c()

file.exists(paste0(sample_output_folder, "datH.rds")) && file.exists(paste0(sample_output_folder, "datD.rds"))
paste0(sample_output_folder, "datH.rds")


# Loop over all diseased samples
for (sample in diseased_samples) {
  # Define the sample output folder
  sample_output_folder <- paste0(celltype_output_folder, "subfolder_sample_", sample, "/")
  
  # Check if both datH.rds and datD.rds files exist before proceeding
  if (file.exists(paste0(sample_output_folder, "datH.rds")) && file.exists(paste0(sample_output_folder, "datD.rds"))) {
    
    # Read in the data for this sample
    datH <- readRDS(file = paste0(sample_output_folder, "datH.rds"))
    datD <- readRDS(file = paste0(sample_output_folder, "datD.rds"))
    
    # Collect gene expression values for this sample
    expression_H_neg <- datH[, gene_neg]
    expression_D_neg <- datD[, gene_neg]
    
    # Extend the lists of expressions over multiple samples
    all_expression_H_neg <- c(all_expression_H_neg, expression_H_neg)
    all_expression_D_neg <- c(all_expression_D_neg, expression_D_neg)
    
  } else {
    # Print a message to inform that the sample folder was skipped
    message(paste("Skipping sample:", sample, "- Folder or files not found"))
  }
}

# Calculate densities for averaged expression values across all samples
density_H_neg <- density(all_expression_H_neg)
density_H_neg$y <- pmax(density_H_neg$y, 0)  # Ensure no negative values for y-axis

density_D_neg <- density(all_expression_D_neg)
density_D_neg$y <- pmax(density_D_neg$y, 0)

# Set x-axis and y-axis to start from 0, and add a slight margin to max values
xlim_neg <- c(0, max(density_H_neg$x, all_expression_D_neg) * 1.05)
ylim_neg <- c(0, max(density_H_neg$y) * 1.05)

# Plot the density for datH (in blue) with the same settings as before
plot(density_H_neg, main = paste(gene_neg), 
     xlab = "", ylab = "", col = "blue",
     xlim = xlim_neg, ylim = ylim_neg, lwd = 4, cex.axis = 2,
     xaxs = "i", yaxs = "i")

# Add red dashed line and points for datD values
abline(v = all_expression_D_neg, col = "red", lwd = 4, lty = 2)
points(all_expression_D_neg, rep(0, length(all_expression_D_neg)), col = "red", pch = 19)

# Optionally, add annotations if needed (commented out for now)
# text(all_expression_D_neg, rep(0, length(all_expression_D_neg)), labels = paste("D:", round(all_expression_D_neg, 2)), pos = 3, col = "red")

# Create two separate data sets: one for all H values and one for all D values
boxplot_data <- list(
  H = all_expression_H_neg,  # Gene expression from datH across multiple samples
  D = all_expression_D_neg   # Gene expression from datD across multiple samples
)

# Define the y-axis limit to start from 0 and go slightly above the maximum value
ylim_neg <- c(0, max(c(all_expression_H_neg, all_expression_D_neg)) * 1.05)
dev.off()
# Plot the boxplot with same styling principles
boxplot(boxplot_data,
        main = paste("Gene Expression for", gene_neg),
        xlab = "Samples", ylab = "Expression Level",
        col = c("blue", "red"),  # Blue for H and Red for D
        lwd = 2,                 # Line width for the box plot borders
        cex.axis = 2,
        ylim = ylim_neg)            # Larger font for axis labels














#comparing conditions

metadata_subset = subset(metadata, subset = Disease_level1 != "Control" & Sex == "female")
diseased_samples = as.vector(metadata_subset$Sample)

gene_neg <- "DLG4"

# Initialize vectors to collect expression values across samples
all_expression_D_neg1 <- c()
all_expression_D_neg2 <- c()

file.exists(paste0(sample_output_folder, "datH.rds")) && file.exists(paste0(sample_output_folder, "datD.rds"))
paste0(sample_output_folder, "datH.rds")

# Loop over all diseased samples
for (sample in diseased_samples) {
  # Define the sample output folder
  sample_output_folder <- paste0(celltype_output_folder, "subfolder_sample_", sample, "/")
  
  # Check if both datH.rds and datD.rds files exist before proceeding
  if (file.exists(paste0(sample_output_folder, "datH.rds")) && file.exists(paste0(sample_output_folder, "datD.rds"))) {
    
    # Read in the data for this sample
    
    datD <- readRDS(file = paste0(sample_output_folder, "datD.rds"))
    
    # Collect gene expression values for this sample
    
    expression_D_neg <- datD[, gene_neg]
    
    # Extend the lists of expressions over multiple samples
    
    all_expression_D_neg1 <- c(all_expression_D_neg1, expression_D_neg)
    
  } else {
    # Print a message to inform that the sample folder was skipped
    message(paste("Skipping sample:", sample, "- Folder or files not found"))
  }
}

metadata_subset = subset(metadata, subset = Disease_level1 != "Control" & Sex == "male")
diseased_samples = as.vector(metadata_subset$Sample)

file.exists(paste0(sample_output_folder, "datH.rds")) && file.exists(paste0(sample_output_folder, "datD.rds"))
paste0(sample_output_folder, "datH.rds")


# Loop over all diseased samples
for (sample in diseased_samples) {
  # Define the sample output folder
  sample_output_folder <- paste0(celltype_output_folder, "subfolder_sample_", sample, "/")
  
  # Check if both datH.rds and datD.rds files exist before proceeding
  if (file.exists(paste0(sample_output_folder, "datH.rds")) && file.exists(paste0(sample_output_folder, "datD.rds"))) {
    
    # Read in the data for this sample
    #datH <- readRDS(file = paste0(sample_output_folder, "datH.rds"))
    datD <- readRDS(file = paste0(sample_output_folder, "datD.rds"))
    
    # Collect gene expression values for this sample
    #expression_H_neg <- datH[, gene_neg]
    expression_D_neg <- datD[, gene_neg]
    
    # Extend the lists of expressions over multiple samples
    #all_expression_H_neg <- c(all_expression_H_neg, expression_H_neg)
    all_expression_D_neg2 <- c(all_expression_D_neg2, expression_D_neg)
    
  } else {
    # Print a message to inform that the sample folder was skipped
    message(paste("Skipping sample:", sample, "- Folder or files not found"))
  }
}

# Optionally, add annotations if needed (commented out for now)
# text(all_expression_D_neg, rep(0, length(all_expression_D_neg)), labels = paste("D:", round(all_expression_D_neg, 2)), pos = 3, col = "red")

# Create two separate data sets: one for all H values and one for all D values
boxplot_data <- list(
  H = all_expression_D_neg1,  # Gene expression from datH across multiple samples
  D = all_expression_D_neg2   # Gene expression from datD across multiple samples
)

# Define the y-axis limit to start from 0 and go slightly above the maximum value
ylim_neg <- c(0, max(c(all_expression_D_neg1, all_expression_D_neg2)) * 1.05)
dev.off()
# Plot the boxplot with same styling principles
boxplot(boxplot_data,
        main = paste("Gene Expression for", gene_neg),
        xlab = "Samples", ylab = "Expression Level",
        col = c("orange", "red"),  # Blue for H and Red for D
        lwd = 2,                 # Line width for the box plot borders
        cex.axis = 2,
        ylim = ylim_neg)            # Larger font for axis labels

