#restart R
setwd(".../Envz/SoupX4")
library(renv)
renv::activate()

library(SoupX)
library(Seurat)


# Set the path to your directory
path <- ".../Atlas/Atlas_Extension/Cal_CKD_raw/"  

# List subdirectories
subfolders <- list.dirs(path = path, full.names = FALSE, recursive = FALSE)

# Print the list of subdirectories
print(subfolders)

# Loop through the subfolders
for (subfolder in subfolders) {
  
  setwd(paste0(path, subfolder))
  
  outputfile <- "_soupx_filt"
  
  toc = Seurat::Read10X_h5("filtered_feature_bc_matrix.h5")$`Gene Expression`
  tod = Seurat::Read10X_h5("raw_feature_bc_matrix.h5")$`Gene Expression`
  sc = SoupChannel(tod, toc)
  
  srat  <- CreateSeuratObject(counts = toc)
  
  #we need clustering information
  srat    <- SCTransform(srat, verbose = F)
  srat    <- RunPCA(srat, verbose = F)
  srat    <- RunUMAP(srat, dims = 1:30, verbose = F)
  srat    <- FindNeighbors(srat, dims = 1:30, verbose = F)
  srat    <- FindClusters(srat, verbose = T)
  
  #add information to soup object
  meta    <- srat@meta.data
  umap    <- srat@reductions$umap@cell.embeddings
  sc  <- setClusters(sc, setNames(meta$seurat_clusters, rownames(meta)))
  sc  <- setDR(sc, umap)
  head(meta)
  
  #run the soup function 
  sc  <- autoEstCont(sc)
  
  head(sc$soupProfile[order(sc$soupProfile$est, decreasing = T), ], n = 20)
  
  adj.matrix  <- adjustCounts(sc, roundToInt = T)
  DropletUtils:::write10xCounts(outputfile, adj.matrix)
  
}







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



# Initialize a list to store Seurat objects
seurat_objects <- list()


file_index <- 1

# Set the path to your directory
path <- ".../Atlas/Atlas_Extension/Cal_CKD_raw/"  

# List subdirectories
subfolders <- list.dirs(path = path, full.names = FALSE, recursive = FALSE)


for (sample_name in subfolders) {
  doubletfinder_output_dir <- file.path(path, sample_name, "doubletfinder")
  dir.create(doubletfinder_output_dir, recursive = TRUE, showWarnings = FALSE)
  outs_path <- file.path(path, sample_name, "_soupx_filt")
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

setwd(".../Atlas/Atlas_Extension/Cal_CKD_raw")
saveRDS(object, file = 'Cal_CKD_new_object.rds')

#convert anndata
object_anndata <- convertFormat(object, from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)
write_h5ad(object_anndata, filename = "Cal_CKD_new_object.h5ad")

sessionInfo()






> sessionInfo()
R version 4.0.2 (2020-06-22)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS:   /gpfs/fs02/apps/R/4.0.2/lib64/R/lib/libRblas.so
LAPACK: /gpfs/fs02/apps/R/4.0.2/lib64/R/lib/libRlapack.so

locale:
  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
[3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
[5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
[7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
[9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
  [1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
  [1] SeuratDisk_0.0.0.9020 patchwork_1.1.2       sceasy_0.0.7         
[4] reticulate_1.25       anndata_0.7.5.5       SeuratObject_4.1.3   
[7] Seurat_4.2.1          dplyr_1.0.10          DoubletFinder_2.0.3  
