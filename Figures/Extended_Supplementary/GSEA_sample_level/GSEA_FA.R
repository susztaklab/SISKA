library(fgsea)
library(dplyr)  # For easier DataFrame manipulation

# Load the RDS file
gene_set_data <- readRDS("/home/isilon/users/o_kloetzer/Atlas/Revision/SISKA/output_SISKA_crossspecies_healthy/go_sets_modified_KEGG.rds")

# Assuming gene_set_data is a data frame or matrix with genes as row names and gene sets as column names
# Transpose the data if needed (if gene sets are rows and genes are columns)
if (all(rownames(gene_set_data) %in% colnames(gene_set_data))) {
  gene_set_data <- t(gene_set_data)
}

# Initialize an empty list to store GMT format data
gmt_list <- list()

# Iterate over each column (gene set)
for (gene_set in colnames(gene_set_data)) {
  # Extract genes that are part of the gene set
  genes_in_set <- rownames(gene_set_data)[gene_set_data[, gene_set] == 1]
  
  # Create a GMT entry: gene set name, description (empty here), and the list of genes
  gmt_list[[gene_set]] <- c(gene_set, "", genes_in_set)
}

# Convert the list into a character vector for writing to a file
gmt_lines <- sapply(gmt_list, function(x) paste(x, collapse = "\t"))

# Write to a GMT file
output_file <- "/home/isilon/users/o_kloetzer/Atlas/Revision/SISKA/output_SISKA_crossspecies_healthy/kegg_human.gmt"
writeLines(gmt_lines, con = output_file)



output_folder_base = "/home/isilon/users/o_kloetzer/Atlas/Revision/SISKA/output_SISKA_crossspecies_healthy/"
#output_folder_base = "/home/isilon/users/o_kloetzer/Atlas/Revision/lungatlas/lungatlas_clean/output_nonTumorRef_Primary/"

cell_types = "PT"
cell_type = cell_types

celltype_output_folder <- paste0(output_folder_base, cell_type, "/")

# List of samples
sample_list <- c("LYZ2_eYFP", "PDGFR-eYFP", "GGT_eYFP", "FA3")  

# Path to GMT file (gene set definitions)
gmt.file <- "/home/isilon/users/o_kloetzer/Atlas/Revision/SISKA/output_SISKA_crossspecies_healthy/kegg_human.gmt"
pathways <- gmtPathways(gmt.file)

# Extract all pathway names to set up the results DataFrame
pathway_names <- names(pathways)

# Initialize an empty DataFrame with pathways as columns and samples as rows, filled with NA
pvalue_df <- as.data.frame(matrix(NA, nrow = length(sample_list), ncol = length(pathway_names)))
ens_df <- as.data.frame(matrix(NA, nrow = length(sample_list), ncol = length(pathway_names)))
colnames(pvalue_df) <- pathway_names
rownames(pvalue_df) <- sample_list
colnames(ens_df) <- pathway_names
rownames(ens_df) <- sample_list

threshold = 0
# Loop through each sample and perform the analysis
for (sample in sample_list) {
  
  # Create pseudobulk for the query sample
  sample_output_folder <- paste0(celltype_output_folder, "subfolder_sample_", sample, "/")
  
  datH <- readRDS(file = paste0(sample_output_folder, "datH.rds"))
  datD <- readRDS(file = paste0(sample_output_folder, "datD.rds"))
  
  # Calculate the average expression for each gene in datH and datD
  avg_expression_H <- colMeans(datH)
  avg_expression_D <- colMeans(datD)
  
  # Filter genes based on a threshold in both datasets
  filtered_genes <- names(avg_expression_H)[avg_expression_H > threshold & avg_expression_D > threshold]
  
  # Subset data based on filtered genes
  dat_subsetH <- datH[, filtered_genes, drop = FALSE]
  dat_subsetD <- datD[, filtered_genes, drop = FALSE]
  
  # Calculate the mean and standard deviation for each gene in datH
  means_H <- apply(dat_subsetH, 2, mean)
  sds_H <- apply(dat_subsetH, 2, sd)
  
  # Transform the expression values in datD into z-scores
  dat_subsetD_z <- sweep(dat_subsetD, 2, means_H, "-")
  dat_subsetD_z <- sweep(dat_subsetD_z, 2, sds_H, "/")
  
  # Convert dat_subsetD_z to a data frame
  dat_subsetD_z_df <- as.data.frame(t(dat_subsetD_z))  # Transpose to get genes as rows
  
  # Create a named vector for z-scores with gene identifiers
  z_scores <- as.vector(dat_subsetD_z_df[, 1])
  ranks <- setNames(z_scores, rownames(dat_subsetD_z_df))
  
  # Sort the ranks in decreasing order (for GSEA)
  ranks <- sort(ranks, decreasing = TRUE)
  
  # Run fGSEA
  set.seed(42)
  fgseaRes <- fgsea(pathways = pathways, stats = ranks, minSize = 10, maxSize = 200)
  
  # Extract the padj values from fgseaRes
  padj_values <- fgseaRes$padj
  ens_values <- fgseaRes$NES
  
  # Map the pathways with padj values back to the main dataframe
  pathway_padjs <- setNames(padj_values, fgseaRes$pathway)
  pathway_ens <- setNames(ens_values, fgseaRes$pathway)
  
  
  # Fill in the corresponding row in pvalue_df for this sample
  pvalue_df[sample, names(pathway_padjs)] <- pathway_padjs
  ens_df[sample, names(pathway_ens)] <- pathway_ens
}

# Now pvalue_df will have the padj values for each pathway (columns) and each sample (rows)
# You can write this DataFrame to a CSV file for further analysis
write.csv(pvalue_df, file = "/home/isilon/users/o_kloetzer/Atlas/Revision/SISKA/output_SISKA_crossspecies_healthy_GSEA/Pval_PT.csv", row.names = TRUE)
write.csv(ens_df, file = "/home/isilon/users/o_kloetzer/Atlas/Revision/SISKA/output_SISKA_crossspecies_healthy_GSEA/NES_PT.csv", row.names = TRUE)















library(fgsea)
library(dplyr)  # For easier DataFrame manipulation

# Load the RDS file
gene_set_data <- readRDS("/home/isilon/users/o_kloetzer/Atlas/Revision/scSpectra_Revision/scSpectra_Mousekidney/output/go_sets_modified.rds")

# Assuming gene_set_data is a data frame or matrix with genes as row names and gene sets as column names
# Transpose the data if needed (if gene sets are rows and genes are columns)
if (all(rownames(gene_set_data) %in% colnames(gene_set_data))) {
  gene_set_data <- t(gene_set_data)
}

# Initialize an empty list to store GMT format data
gmt_list <- list()

# Iterate over each column (gene set)
for (gene_set in colnames(gene_set_data)) {
  # Extract genes that are part of the gene set
  genes_in_set <- rownames(gene_set_data)[gene_set_data[, gene_set] == 1]
  
  # Create a GMT entry: gene set name, description (empty here), and the list of genes
  gmt_list[[gene_set]] <- c(gene_set, "", genes_in_set)
}

# Convert the list into a character vector for writing to a file
gmt_lines <- sapply(gmt_list, function(x) paste(x, collapse = "\t"))

# Write to a GMT file
output_file <- "/home/isilon/users/o_kloetzer/Atlas/Revision/scSpectra_Revision/scSpectra_Mousekidney/output/kegg_mouse.gmt"
writeLines(gmt_lines, con = output_file)



output_folder_base = "/home/isilon/users/o_kloetzer/Atlas/Revision/scSpectra_Revision/scSpectra_Mousekidney/output/"
#output_folder_base = "/home/isilon/users/o_kloetzer/Atlas/Revision/lungatlas/lungatlas_clean/output_nonTumorRef_Primary/"

cell_types = "PT"
cell_type = cell_types

celltype_output_folder <- paste0(output_folder_base, cell_type, "/")

# List of samples
sample_list <- c("LYZ2_eYFP", "PDGFR-eYFP", "GGT_eYFP", "FA3")  

# Path to GMT file (gene set definitions)
gmt.file <- "/home/isilon/users/o_kloetzer/Atlas/Revision/scSpectra_Revision/scSpectra_Mousekidney/output/kegg_mouse.gmt"
pathways <- gmtPathways(gmt.file)

# Extract all pathway names to set up the results DataFrame
pathway_names <- names(pathways)

# Initialize an empty DataFrame with pathways as columns and samples as rows, filled with NA
pvalue_df <- as.data.frame(matrix(NA, nrow = length(sample_list), ncol = length(pathway_names)))
ens_df <- as.data.frame(matrix(NA, nrow = length(sample_list), ncol = length(pathway_names)))
colnames(pvalue_df) <- pathway_names
rownames(pvalue_df) <- sample_list
colnames(ens_df) <- pathway_names
rownames(ens_df) <- sample_list

threshold = 0
# Loop through each sample and perform the analysis
for (sample in sample_list) {
  
  # Create pseudobulk for the query sample
  sample_output_folder <- paste0(celltype_output_folder, "subfolder_sample_", sample, "/")
  
  datH <- readRDS(file = paste0(sample_output_folder, "datH.rds"))
  datD <- readRDS(file = paste0(sample_output_folder, "datD.rds"))
  
  # Calculate the average expression for each gene in datH and datD
  avg_expression_H <- colMeans(datH)
  avg_expression_D <- colMeans(datD)
  
  # Filter genes based on a threshold in both datasets
  filtered_genes <- names(avg_expression_H)[avg_expression_H > threshold & avg_expression_D > threshold]
  
  # Subset data based on filtered genes
  dat_subsetH <- datH[, filtered_genes, drop = FALSE]
  dat_subsetD <- datD[, filtered_genes, drop = FALSE]
  
  # Calculate the mean and standard deviation for each gene in datH
  means_H <- apply(dat_subsetH, 2, mean)
  sds_H <- apply(dat_subsetH, 2, sd)
  
  # Transform the expression values in datD into z-scores
  dat_subsetD_z <- sweep(dat_subsetD, 2, means_H, "-")
  dat_subsetD_z <- sweep(dat_subsetD_z, 2, sds_H, "/")
  
  # Convert dat_subsetD_z to a data frame
  dat_subsetD_z_df <- as.data.frame(t(dat_subsetD_z))  # Transpose to get genes as rows
  
  # Create a named vector for z-scores with gene identifiers
  z_scores <- as.vector(dat_subsetD_z_df[, 1])
  ranks <- setNames(z_scores, rownames(dat_subsetD_z_df))
  
  # Sort the ranks in decreasing order (for GSEA)
  ranks <- sort(ranks, decreasing = TRUE)
  
  # Run fGSEA
  set.seed(42)
  fgseaRes <- fgsea(pathways = pathways, stats = ranks, minSize = 10, maxSize = 200)
  
  # Extract the padj values from fgseaRes
  padj_values <- fgseaRes$padj
  ens_values <- fgseaRes$NES
  
  # Map the pathways with padj values back to the main dataframe
  pathway_padjs <- setNames(padj_values, fgseaRes$pathway)
  pathway_ens <- setNames(ens_values, fgseaRes$pathway)
  
  
  # Fill in the corresponding row in pvalue_df for this sample
  pvalue_df[sample, names(pathway_padjs)] <- pathway_padjs
  ens_df[sample, names(pathway_ens)] <- pathway_ens
}

# Now pvalue_df will have the padj values for each pathway (columns) and each sample (rows)
# You can write this DataFrame to a CSV file for further analysis
write.csv(pvalue_df, file = "/home/isilon/users/o_kloetzer/Atlas/Revision/scSpectra_Revision/scSpectra_Mousekidney/output/Pval_GSEA/Pval_PT.csv", row.names = TRUE)
write.csv(ens_df, file = "/home/isilon/users/o_kloetzer/Atlas/Revision/scSpectra_Revision/scSpectra_Mousekidney/output/R2_GSEA/NES_PT.csv", row.names = TRUE)

