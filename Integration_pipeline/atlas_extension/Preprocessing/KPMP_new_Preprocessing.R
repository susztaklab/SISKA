library(DoubletFinder)
library(dplyr)
library(Seurat)
library(anndata)
library(reticulate)
library(sceasy)
library(patchwork)
library(SeuratDisk)

runDoubletFinder <- function(seurat_obj, PCs = 1:30, sct = FALSE, pN = 0.25, doublet_rate = 0.1, dir_path) {
  
  setwd(dir_path)
  
  obj_name <- deparse(substitute(seurat_obj))
  cat("Processing Seurat object:", obj_name, "\n")  # Print the object name
  
  # pK Identification
  sweep.res.list <- paramSweep_v3(seurat_obj, PCs = PCs, sct = sct)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  
  pdf("bcmvn.pdf", width = 12, height = 6)
  bcmvn <- find.pK(sweep.stats)
  dev.off()
  
  pK <- bcmvn$pK[which.max(bcmvn$BCmetric)]
  pK <- as.numeric(as.vector(pK))
  
  # Homotypic Doublet Proportion Estimate
  annotations <- seurat_obj@meta.data$RNA_snn_res.0.5
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(doublet_rate * nrow(seurat_obj@meta.data))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))
  
  # Run DoubletFinder
  seurat_obj <- doubletFinder_v3(seurat_obj, PCs = PCs, pN = pN, pK = pK, nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = sct)
  
  # Determine slotname and display results
  slotname <- paste0("DF.classifications_", pN, "_", pK, "_", nExp_poi.adj)
  print(table(seurat_obj[[slotname]]))
  print(pK)
  print(slotname)
  print(seurat_obj)
  
  # Store results in a new metadata column "doubletfinder"
  seurat_obj$doubletfinder <- seurat_obj[[slotname]]
  
  # Plot
  pdf("dimplot.pdf", width = 6, height = 6)
  print(DimPlot(seurat_obj, group.by = slotname))
  dev.off()
  
  return(seurat_obj)
}

#we use this code to spearate multiome samples (h5 matrix will be used) from standard samples (exactly 3 files in outs)

# Define the base directory
base_dir <- ".../Atlas/Atlas_Extension/KPMP_data_additional"

# List all sample directories
sample_dirs <- list.dirs(base_dir, full.names = TRUE, recursive = FALSE)

# Initialize vectors for the two categories
samples_with_3_files <- c()
samples_with_more_than_3_files <- c()

# Loop through each sample directory
for (sample_dir in sample_dirs) {
  # Construct the path to the "outs" subfolder
  outs_path <- file.path(sample_dir, "outs")
  
  # Count the number of files in the "outs" subfolder
  if (dir.exists(outs_path)) {
    file_count <- length(list.files(outs_path))
    
    # Classify the sample based on file count
    if (file_count == 3) {
      samples_with_3_files <- c(samples_with_3_files, basename(sample_dir))
    } else if (file_count > 3) {
      samples_with_more_than_3_files <- c(samples_with_more_than_3_files, basename(sample_dir))
    }
  }
}

# Print the sample lists
print(samples_with_3_files)
print(samples_with_more_than_3_files)

# Define the base directory
base_dir <- ".../Atlas/Atlas_Extension/KPMP_data_additional"

# Initialize a list to store Seurat objects
seurat_objects <- list()

# Process samples with more than 3 files
file_index <- 1
for (sample_name in samples_with_more_than_3_files) {
  outs_path <- file.path(base_dir, sample_name, "outs")
  h5_file <- file.path(outs_path, "filtered_feature_bc_matrix.h5")
  
  if (file.exists(h5_file)) {
    # Construct the DoubletFinder output directory for the sample
    doubletfinder_output_dir <- file.path(base_dir, sample_name, "doubletfinder")
    dir.create(doubletfinder_output_dir, recursive = TRUE, showWarnings = FALSE)
    counts <- Read10X_h5(h5_file)
    gene_expression_matrix <- counts$`Gene Expression`
    seurat_objects[[file_index]] <- CreateSeuratObject(counts = gene_expression_matrix, project = sample_name, assay = "RNA")
    
    # Calculate the percentage of mitochondrial genes
    seurat_objects[[file_index]][["percent.mt"]] <- PercentageFeatureSet(seurat_objects[[file_index]], pattern = "^MT-")
    
    # Quality Control
    seurat_objects[[file_index]] <- subset(seurat_objects[[file_index]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
    
    # Normalization and Feature Selection
    seurat_objects[[file_index]] <- NormalizeData(seurat_objects[[file_index]])
    seurat_objects[[file_index]] <- FindVariableFeatures(seurat_objects[[file_index]], selection.method = "vst", nfeatures = 2000)
    
    # Scaling
    seurat_objects[[file_index]] <- ScaleData(seurat_objects[[file_index]])
    
    # Dimensionality Reduction and Clustering
    seurat_objects[[file_index]] <- RunPCA(seurat_objects[[file_index]])
    seurat_objects[[file_index]] <- FindNeighbors(seurat_objects[[file_index]], dims = 1:30)
    seurat_objects[[file_index]] <- FindClusters(seurat_objects[[file_index]], resolution = 0.5)
    
    
    seurat_objects[[file_index]] <- runDoubletFinder(seurat_objects[[file_index]], PCs = 1:30, sct = FALSE, pN = 0.25, doublet_rate = 0.1, doubletfinder_output_dir)
    file_index <- file_index + 1
  }
}

# Process samples with exactly 3 files
for (sample_name in samples_with_3_files) {
  doubletfinder_output_dir <- file.path(base_dir, sample_name, "doubletfinder")
  dir.create(doubletfinder_output_dir, recursive = TRUE, showWarnings = FALSE)
  outs_path <- file.path(base_dir, sample_name, "outs")
  seurat_objects[[file_index]] <- CreateSeuratObject(counts = Read10X(outs_path), project = sample_name)
  
  # Calculate the percentage of mitochondrial genes
  seurat_objects[[file_index]][["percent.mt"]] <- PercentageFeatureSet(seurat_objects[[file_index]], pattern = "^MT-")
  
  # Quality Control
  seurat_objects[[file_index]] <- subset(seurat_objects[[file_index]], subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
  
  # Normalization and Feature Selection
  seurat_objects[[file_index]] <- NormalizeData(seurat_objects[[file_index]])
  seurat_objects[[file_index]] <- FindVariableFeatures(seurat_objects[[file_index]], selection.method = "vst", nfeatures = 2000)
  
  # Scaling
  seurat_objects[[file_index]] <- ScaleData(seurat_objects[[file_index]])
  
  # Dimensionality Reduction and Clustering
  seurat_objects[[file_index]] <- RunPCA(seurat_objects[[file_index]])
  seurat_objects[[file_index]] <- FindNeighbors(seurat_objects[[file_index]], dims = 1:30)
  seurat_objects[[file_index]] <- FindClusters(seurat_objects[[file_index]], resolution = 0.5)
  
  seurat_objects[[file_index]] <- runDoubletFinder(seurat_objects[[file_index]], PCs = 1:30, sct = FALSE, pN = 0.25, doublet_rate = 0.1, doubletfinder_output_dir)
  file_index <- file_index + 1
}


object <- merge(seurat_objects[[1]], y = seurat_objects[-1])

#percent_mt
object[["percent.mt"]] <- PercentageFeatureSet(object = object, pattern = "^MT-")

setwd(".../Atlas/Atlas_Extension/KPMP_data_additional")
saveRDS(object, file = 'KPMP_new_object_doubletfinder.rds')

#convert anndata
object_anndata <- convertFormat(object, from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)
write_h5ad(object_anndata, filename = "KPMP_new_object_doubletfinder.h5ad")
