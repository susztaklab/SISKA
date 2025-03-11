# our Standard Script for Preprocessing from Seurat Objects. Run for each sample and merge.

library (Seurat)
library (dplyr)
library (ggplot2)
library (harmony)
library (patchwork)
library (SoupX)
library (DoubletFinder)

#================================Ambient RNA Correction===============================#
Control1 = load10X('path/to/your/cellranger/outs/folder')
Control1 = autoEstCont(Control1)
Control1 = adjustCounts(Control1)
#Do it for all samples
#================================Doublet Removal===============================#
Control1.S <- CreateSeuratObject(counts = Control1, project = "Control1", min.cells = 3, min.features = 300) # create Seurat object
Control1.S <- NormalizeData(Control1.S)
Control1.S <- FindVariableFeatures(Control1.S, selection.method = "vst", nfeatures = 2000)
Control1.S <- ScaleData(Control1.S)
Control1.S <- RunPCA(Control1.S)
Control1.S <- RunUMAP(Control1.S, dims = 1:30)
sweep.res.list_kidney <- paramSweep_v3(Control1.S, PCs = 1:30, sct = FALSE)
sweep.stats_kidney <- summarizeSweep(sweep.res.list_kidney, GT = FALSE)
bcmvn_kidney <- find.pK(sweep.stats_kidney) #determine the pk
nExp_poi <- round(0.15*nrow(Control1.S@meta.data))
Control1.S <- doubletFinder_v3(seu_kidney, PCs = 1:30, pN = 0.25, pK = #depends on the previous sterp, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
Control1.S.Filtered <- subset(Control1.S, subset = #Doublet_related_column == "singlet")


#============================Merge Datasets======================================================#
object = merge (Control1.S.Filtered, y = c (#All filtered objects))
#============================Add percent.mt to dataset======================================================#
object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^mt-")

#============================Pre Processing======================================================#
object <- subset(object, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 15)


