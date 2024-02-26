#restart R
setwd(".../Envz/SoupX4")
library(renv)
renv::activate()
installed.packages() # to see if env is empty")

library(SoupX)
library(Seurat)


setwd(".../data/Humphreys_ADPKD/PKD7")
outputfile <- "GSM5627701_PKD7_soupx_filt"

toc = Seurat::Read10X_h5("GSM5627701_PKD7_filtered_feature_bc_matrix.h5")
tod = Seurat::Read10X_h5("GSM5627701_PKD7_raw_feature_bc_matrix.h5")
sc = SoupChannel(tod, toc)

srat  <- CreateSeuratObject(counts = toc)
srat

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











setwd(".../data/Humphreys_ADPKD/PKD8")
outputfile <- "GSM5627702_PKD8_soupx_filt"

toc = Seurat::Read10X_h5("GSM5627702_PKD8_filtered_feature_bc_matrix.h5")
tod = Seurat::Read10X_h5("GSM5627702_PKD8_raw_feature_bc_matrix.h5")
sc = SoupChannel(tod, toc)

srat  <- CreateSeuratObject(counts = toc)
srat

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











setwd(".../data/Humphreys_ADPKD/Ctrl1")
outputfile <- "GSM5627690_cont1_soupx_filt"

toc = Seurat::Read10X_h5("GSM5627690_cont1_filtered_feature_bc_matrix.h5")
tod = Seurat::Read10X_h5("GSM5627690_cont1_raw_feature_bc_matrix.h5")
sc = SoupChannel(tod, toc)

srat  <- CreateSeuratObject(counts = toc)
srat

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













setwd(".../data/Humphreys_ADPKD/Ctrl2")
outputfile <- "GSM5627691_cont2_soupx_filt"

toc = Seurat::Read10X_h5("GSM5627691_cont2_filtered_feature_bc_matrix.h5")
tod = Seurat::Read10X_h5("GSM5627691_cont2_raw_feature_bc_matrix.h5")
sc = SoupChannel(tod, toc)

srat  <- CreateSeuratObject(counts = toc)
srat

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












setwd(".../data/Humphreys_ADPKD/Ctrl3")
outputfile <- "GSM5627692_cont3_soupx_filt"

toc = Seurat::Read10X_h5("GSM5627692_cont3_filtered_feature_bc_matrix.h5")
tod = Seurat::Read10X_h5("GSM5627692_cont3_raw_feature_bc_matrix.h5")
sc = SoupChannel(tod, toc)

srat  <- CreateSeuratObject(counts = toc)
srat

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












setwd(".../data/Humphreys_ADPKD/Ctrl4")
outputfile <- "GSM5627693_cont4_soupx_filt"

toc = Seurat::Read10X_h5("GSM5627693_cont4_filtered_feature_bc_matrix.h5")
tod = Seurat::Read10X_h5("GSM5627693_cont4_raw_feature_bc_matrix.h5")
sc = SoupChannel(tod, toc)

srat  <- CreateSeuratObject(counts = toc)
srat

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













setwd(".../data/Humphreys_ADPKD/Ctrl5")
outputfile <- "GSM5627694_cont5_soupx_filt"

toc = Seurat::Read10X_h5("GSM5627694_cont5_filtered_feature_bc_matrix.h5")
tod = Seurat::Read10X_h5("GSM5627694_cont5_raw_feature_bc_matrix.h5")
sc = SoupChannel(tod, toc)

srat  <- CreateSeuratObject(counts = toc)
srat

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
  [1] SeuratObject_4.1.3 Seurat_4.3.0       SoupX_1.6.2        renv_0.16.0       




#new session


library(dplyr)
library(Seurat)
library(anndata)

library(reticulate)
library(sceasy)

library(patchwork)
library(SeuratDisk)

#set wd and meta
setwd(".../data/Humphreys_ADPKD")
meta_data <- read.csv(file = 'GSE185948_metadata_RNA.csv')
# Assume the following data frame

# Create new column
meta_data$barcodes_unique_new <- paste0(meta_data$patient, "_", meta_data$barcode)

# Print data frame
rownames(meta_data) <- meta_data$barcodes_unique_new


#####################continue with samples######################

setwd(".../data/Humphreys_ADPKD/PKD1/")
sample_name = "PKD1"

file <- CreateSeuratObject(counts = Read10X("GSM5627695_PKD1_soupx_filt"))
file <- RenameCells(file, add.cell.id = sample_name)
Idents(file) <- file@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$patient == "PKD1") 
file <- subset(x = file, idents = rownames(meta))

file1 <- AddMetaData(file, meta)




#####################continue with samples######################

setwd(".../data/Humphreys_ADPKD/PKD2/")
sample_name = "PKD2"

file <- CreateSeuratObject(counts = Read10X("GSM5627696_PKD2_soupx_filt"))
file <- RenameCells(file, add.cell.id = sample_name)
Idents(file) <- file@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$patient == sample_name) 
file <- subset(x = file, idents = rownames(meta))

file2 <- AddMetaData(file, meta)





#####################continue with samples######################

setwd(".../data/Humphreys_ADPKD/PKD3/")
sample_name = "PKD3"

file <- CreateSeuratObject(counts = Read10X("GSM5627697_PKD3_soupx_filt"))
file <- RenameCells(file, add.cell.id = sample_name)
Idents(file) <- file@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$patient == sample_name) 
file <- subset(x = file, idents = rownames(meta))

file3 <- AddMetaData(file, meta)







#####################continue with samples######################

setwd(".../data/Humphreys_ADPKD/PKD4/")
sample_name = "PKD4"

file <- CreateSeuratObject(counts = Read10X("GSM5627698_PKD4_soupx_filt"))
file <- RenameCells(file, add.cell.id = sample_name)
Idents(file) <- file@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$patient == sample_name) 
file <- subset(x = file, idents = rownames(meta))

file4 <- AddMetaData(file, meta)










#####################continue with samples######################

setwd(".../data/Humphreys_ADPKD/PKD5/")
sample_name = "PKD5"

file <- CreateSeuratObject(counts = Read10X("GSM5627699_PKD5_soupx_filt"))
file <- RenameCells(file, add.cell.id = sample_name)
Idents(file) <- file@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$patient == sample_name) 
file <- subset(x = file, idents = rownames(meta))

file5 <- AddMetaData(file, meta)










#####################continue with samples######################

setwd(".../data/Humphreys_ADPKD/PKD6/")
sample_name = "PKD6"

file <- CreateSeuratObject(counts = Read10X("GSM5627700_PKD6_soupx_filt"))
file <- RenameCells(file, add.cell.id = sample_name)
Idents(file) <- file@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$patient == sample_name) 
file <- subset(x = file, idents = rownames(meta))

file6 <- AddMetaData(file, meta)








#####################continue with samples######################

setwd(".../data/Humphreys_ADPKD/PKD7/")
sample_name = "PKD7"

file <- CreateSeuratObject(counts = Read10X("GSM5627701_PKD7_soupx_filt"))
file <- RenameCells(file, add.cell.id = sample_name)
Idents(file) <- file@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$patient == sample_name) 
file <- subset(x = file, idents = rownames(meta))

file7 <- AddMetaData(file, meta)









#####################continue with samples######################

setwd(".../data/Humphreys_ADPKD/PKD8/")
sample_name = "PKD8"

file <- CreateSeuratObject(counts = Read10X("GSM5627702_PKD8_soupx_filt"))
file <- RenameCells(file, add.cell.id = sample_name)
Idents(file) <- file@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$patient == sample_name) 
file <- subset(x = file, idents = rownames(meta))

file8 <- AddMetaData(file, meta)









#####################continue with samples######################

setwd(".../data/Humphreys_ADPKD/Ctrl1")
sample_name = "control1"

file <- CreateSeuratObject(counts = Read10X("GSM5627690_cont1_soupx_filt"))
file <- RenameCells(file, add.cell.id = sample_name)
Idents(file) <- file@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$patient == sample_name) 
file <- subset(x = file, idents = rownames(meta))

file9 <- AddMetaData(file, meta)





#####################continue with samples######################

setwd(".../data/Humphreys_ADPKD/Ctrl2")
sample_name = "control2"

file <- CreateSeuratObject(counts = Read10X("GSM5627691_cont2_soupx_filt"))
file <- RenameCells(file, add.cell.id = sample_name)
Idents(file) <- file@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$patient == sample_name) 
file <- subset(x = file, idents = rownames(meta))

file10 <- AddMetaData(file, meta)






#####################continue with samples######################

setwd(".../data/Humphreys_ADPKD/Ctrl3")
sample_name = "control3"

file <- CreateSeuratObject(counts = Read10X("GSM5627692_cont3_soupx_filt"))
file <- RenameCells(file, add.cell.id = sample_name)
Idents(file) <- file@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$patient == sample_name) 
file <- subset(x = file, idents = rownames(meta))

file11 <- AddMetaData(file, meta)








#####################continue with samples######################

setwd(".../data/Humphreys_ADPKD/Ctrl4")
sample_name = "control4"

file <- CreateSeuratObject(counts = Read10X("GSM5627693_cont4_soupx_filt"))
file <- RenameCells(file, add.cell.id = sample_name)
Idents(file) <- file@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$patient == sample_name) 
file <- subset(x = file, idents = rownames(meta))

file12 <- AddMetaData(file, meta)












#####################continue with samples######################

setwd(".../data/Humphreys_ADPKD/Ctrl5")
sample_name = "control5"

file <- CreateSeuratObject(counts = Read10X("GSM5627694_cont5_soupx_filt"))
file <- RenameCells(file, add.cell.id = sample_name)
Idents(file) <- file@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$patient == sample_name) 
file <- subset(x = file, idents = rownames(meta))

file13 <- AddMetaData(file, meta)







ADPKD_object <- merge(file1, y = c(file2, file3,
                                   file4, file5,
                                   file6, file7,
                                   file8, file9,
                                   file10, file11,
                                   file12, file13))

setwd(".../data/Humphreys_ADPKD")
saveRDS(ADPKD_object, file = 'ADPKD_object.rds')




#percent_mt
ADPKD_object[["percent.mt"]] <- PercentageFeatureSet(object = ADPKD_object, pattern = "^MT-")






object <- ADPKD_object
object_anndata <- convertFormat(object, from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)
write_h5ad(object_anndata, filename = "ADPKD_object.h5ad")

saveRDS(ADPKD_object, file = 'ADPKD_object.rds')



