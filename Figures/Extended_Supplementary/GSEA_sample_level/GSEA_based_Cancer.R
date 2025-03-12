library(fgsea)
library(dplyr)  # For easier DataFrame manipulation

output_folder_base = "/home/isilon/users/o_kloetzer/Atlas/Revision/lungatlas/lungatlas_clean/output_CancerRef_Primary/"
#output_folder_base = "/home/isilon/users/o_kloetzer/Atlas/Revision/lungatlas/lungatlas_clean/output_nonTumorRef_Primary/"

cell_types = "Cancer"

celltype_output_folder <- paste0(output_folder_base, cell_type, "/")

# List of samples
sample_list <- c("Goveia_Carmeliet_2020_patient_32",
                 "Goveia_Carmeliet_2020_patient_40",
                 "Goveia_Carmeliet_2020_patient_41",
                 "Goveia_Carmeliet_2020_patient_45",
                 "Goveia_Carmeliet_2020_patient_50",
                 "He_Fan_2021_P2",
                 "He_Fan_2021_P5",
                 "Kim_Lee_2020_P0028",
                 "Kim_Lee_2020_P0031",
                 "Kim_Lee_2020_P1006",
                 "Kim_Lee_2020_P1028",
                 "Kim_Lee_2020_P1049",
                 "Kim_Lee_2020_P1058",
                 "Lambrechts_Thienpont_2018_6149v1_1",
                 "Lambrechts_Thienpont_2018_6149v1_2",
                 "Lambrechts_Thienpont_2018_6149v2_3",
                 "Lambrechts_Thienpont_2018_6149v2_4",
                 "Lambrechts_Thienpont_2018_6653_6",
                 "Lambrechts_Thienpont_2018_6653_7",
                 "Laughney_Massague_2020_LX675",
                 "Laughney_Massague_2020_LX679")  

# Path to GMT file (gene set definitions)
gmt.file <- "/home/isilon/users/o_kloetzer/Atlas/Revision/lungatlas/lungatlas_clean/output_CancerRef_Primary/c2.cp.v2023.2.Hs.symbols.gmt.txt"
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
  fgseaRes <- fgsea(pathways = pathways, stats = ranks, minSize = 10, maxSize = 100)
  
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
write.csv(pvalue_df, file = "/home/isilon/users/o_kloetzer/Atlas/Revision/lungatlas/lungatlas_clean/CancerRef_GSEA/Pval/Pval_Cancer.csv", row.names = TRUE)
write.csv(ens_df, file = "/home/isilon/users/o_kloetzer/Atlas/Revision/lungatlas/lungatlas_clean/CancerRef_GSEA/Pval/NES_Cancer.csv", row.names = TRUE)






library(fgsea)
library(dplyr)  # For easier DataFrame manipulation

#output_folder_base = "/home/isilon/users/o_kloetzer/Atlas/Revision/lungatlas/lungatlas_clean/output_CancerRef_Primary/"
output_folder_base = "/home/isilon/users/o_kloetzer/Atlas/Revision/lungatlas/lungatlas_clean/output_nonTumorRef_Primary/"

cell_types = "Cancer"

celltype_output_folder <- paste0(output_folder_base, cell_type, "/")

# List of samples
sample_list <- c("Laughney_Massague_2020_LX653",
                 "Laughney_Massague_2020_LX675",
                 "Laughney_Massague_2020_LX676",
                 "Laughney_Massague_2020_LX679",
                 "Laughney_Massague_2020_LX680",
                 "Laughney_Massague_2020_LX682",
                 "Laughney_Massague_2020_LX684")  

# Path to GMT file (gene set definitions)
gmt.file <- "/home/isilon/users/o_kloetzer/Atlas/Revision/lungatlas/lungatlas_clean/output_CancerRef_Primary/c2.cp.v2023.2.Hs.symbols.gmt.txt"
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
  fgseaRes <- fgsea(pathways = pathways, stats = ranks, minSize = 10, maxSize = 100)
  
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
write.csv(pvalue_df, file = "/home/isilon/users/o_kloetzer/Atlas/Revision/lungatlas/lungatlas_clean/nonTumorref_GSEA/Pval/Pval_Cancer.csv", row.names = TRUE)
write.csv(ens_df, file = "/home/isilon/users/o_kloetzer/Atlas/Revision/lungatlas/lungatlas_clean/nonTumorref_GSEA/Pval/NES_Cancer.csv", row.names = TRUE)






