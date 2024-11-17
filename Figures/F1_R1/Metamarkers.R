
library(MetaMarkers)
library(scRNAseq)
library(dplyr)

####################high resolution: Atlas6.6 - leiden 8.0

ob_mouse <- readRDS(".../Atlas6.6_seurat1_mouse_sce.rds")
ob_rat <- readRDS(".../Atlas6.6_seurat1_rat_sce.rds")
ob_human <- readRDS(".../Atlas6.6_seurat1_human_sce.rds")


assay(ob_mouse, "cpm") = convert_to_cpm(assay(ob_mouse))
assay(ob_rat, "cpm") = convert_to_cpm(assay(ob_rat))
assay(ob_human, "cpm") = convert_to_cpm(assay(ob_human))

markers_mouse = compute_markers(assay(ob_mouse, "cpm"), ob_mouse$leiden_scVI_8_0, ob_mouse$annotation_final_level1_harmonised)
markers_rat = compute_markers(assay(ob_rat, "cpm"), ob_rat$leiden_scVI_8_0, ob_rat$annotation_final_level1_harmonised)
markers_human = compute_markers(assay(ob_human, "cpm"), ob_human$leiden_scVI_8_0, ob_human$annotation_final_level1_harmonised)




setwd(".../MetaMarkers/Atlas6.6_Superresolution/Leiden_8.0")
export_markers(markers_mouse, "mouse_markers.csv")
export_markers(markers_human, "human_markers.csv")
export_markers(markers_rat, "rat_markers.csv")

markers = list(
  markers_human,
  markers_mouse,
  markers_rat
)
meta_markers = make_meta_markers(markers, detailed_stats = TRUE)
export_markers(meta_markers, "meta_markers.csv")


########## final global annotation 1B


library(MetaMarkers)
library(scRNAseq)
library(dplyr)



ob_mouse <- readRDS(".../Atlas6.6_seurat1_mouse_sce.rds")
ob_rat <- readRDS(".../Atlas6.6_seurat1_rat_sce.rds")
ob_human <- readRDS(".../Atlas6.6_seurat1_human_sce.rds")


assay(ob_mouse, "cpm") = convert_to_cpm(assay(ob_mouse))
assay(ob_rat, "cpm") = convert_to_cpm(assay(ob_rat))
assay(ob_human, "cpm") = convert_to_cpm(assay(ob_human))

markers_mouse = compute_markers(assay(ob_mouse, "cpm"), ob_mouse$annotation_final_level1B)
markers_rat = compute_markers(assay(ob_rat, "cpm"), ob_rat$annotation_final_level1B)
markers_human = compute_markers(assay(ob_human, "cpm"), ob_human$annotation_final_level1B)


setwd(".../MetaMarkers/Atlas6.6_annotation_final_1B")
export_markers(markers_mouse, "mouse_markers.csv")
export_markers(markers_human, "human_markers.csv")
export_markers(markers_rat, "rat_markers.csv")

markers = list(
  markers_human,
  markers_mouse,
  markers_rat
)
meta_markers = make_meta_markers(markers, detailed_stats = TRUE)
export_markers(meta_markers, "meta_markers.csv")





#############Gene set (DEGO) analysis including converting of object 

library(Seurat)

#converting tools could not be used for the DEGO anndata object. Instead, .csv files of the matrix were exported

counts <- read.csv(".../Atlas/DEGO/Matrix_export/matrix_rat.csv")
meta <- read.csv(".../Atlas/DEGO/Matrix_export/DEGO_R_meta.csv")
row.names(counts) <- counts$X  
row.names(meta) <- meta$X  
meta <- meta[,-1]
counts <- counts[,-1]

seurat_object <- CreateSeuratObject(
  t(counts),
  project = "rat",
  assay = "RNA",
)

seurat_object <- AddMetaData(seurat_object, meta, col.name = NULL)

seurat_object_sce <- as.SingleCellExperiment(seurat_object)

saveRDS(seurat_object_sce, file = ".../Atlas/DEGO/Matrix_export/DEGO_rat_sce.rds")

#same for other species


library(MetaMarkers)
library(scRNAseq)
library(dplyr)

#load the files you are interested in 
ob_mouse <- readRDS(".../Atlas/DEGO/Matrix_export/DEGO_mouse_sce.rds")
ob_rat <- readRDS(".../Atlas/DEGO/Matrix_export/DEGO_rat_sce.rds")
ob_human <- readRDS(".../Atlas/DEGO/Matrix_export/DEGO_human_sce.rds")

#we don't normailze since DEGOs are already based on averaging normalized values

markers_mouse = compute_markers(assay(ob_mouse), ob_mouse$annotation_final_level1B)
setwd(".../Atlas/DEGO/Matrix_export/")
export_markers(markers_mouse, "mouse_markers.csv")


markers_rat = compute_markers(assay(ob_rat), ob_rat$annotation_final_level1B)
export_markers(markers_rat, "rat_markers.csv")


markers_human = compute_markers(assay(ob_human), ob_human$annotation_final_level1B)
export_markers(markers_human, "human_markers.csv")


setwd(".../Atlas/DEGO/Matrix_export")
markers_human = read_meta_markers("human_markers.csv.gz")
markers_mouse = read_meta_markers("mouse_markers.csv.gz")
markers_rat = read_meta_markers("rat_markers.csv.gz")


markers = list(
  markers_human,
  markers_mouse,
  markers_rat
)

meta_markers = make_meta_markers(markers, fc_threshold = 1.1, detailed_stats = TRUE)
export_markers(meta_markers, "meta_markers.csv")

