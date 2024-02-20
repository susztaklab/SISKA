

library(Seurat)
library(dplyr)
library(reticulate)
library(sceasy)
library(anndata)
library(patchwork)
library(SeuratDisk)

#Changing gene symbols in Seurat objects to human one-to-one ortholog 

#import new ensembl gene list 
orthos <- read.csv(file = ".../Genelist_V2_subset.csv")

#import "raw" rds files with original gene symbols
## human
KPMP <- '.../KPMP_objectlarge.h5Seurat' #downloaded from KPMP website
S_HUMAN <- ".../Susztak_human.rds"

## mouse
HUMPHREYS_1M <- '.../GSE184652_RAW/Humphreys_1M_Atlas_final.rds'
HUMPHREYS_IRI <- ".../IRI_mouse.rds"

## rat
RAT_ABEDINI <- ".../Rat.All.Idents.9.12.2022.rds"
RAT_BALZER <- ".../ZSF1_comb_soupX_diet_object.rds"

#set WD
setwd(".../Atlas/objects")

################################
#change to human gene names (applied for each Seurat object)
mousekidney <- readRDS(file = HUMPHREYS_IRI)
DefaultAssay(mousekidney)<-"RNA"
mousegenes<-as.vector(unlist(orthos$Mouse.gene.name))
humangenes<-as.vector(unlist(orthos$Gene.name))
orthos$id <- paste0(orthos$Gene.name,"-",orthos$Mouse.gene.name)
MOUSE_orthos <- DietSeurat(
  mousekidney,
  counts = TRUE,
  data = TRUE,
  scale.data = FALSE,
  features = mousegenes,
  assays = NULL,
  dimreducs = NULL,
  graphs = NULL,
  misc = TRUE
)
MOUSE_orthos_copy <- MOUSE_orthos
genesMOUSE <- MOUSE_orthos@assays$RNA@counts@Dimnames[[1]]

##rename MOUSE
for (i in genesMOUSE){
  ie <-paste0("^",i,"$")  # using regex for exact match (if not it might annotate similar gene names)
  p0 <- which(grepl(ie, genesMOUSE))
  p <- which(grepl(ie,orthos$Mouse.gene.name)) # find position of the match
  p1<-p[1]  # keeping the first position from p
  f<-orthos$Gene.name[p1] # get the human gene name from the same position
  MOUSE_orthos@assays$RNA@counts@Dimnames[[1]][p0] <- as.character(f)
  MOUSE_orthos@assays$RNA@data@Dimnames[[1]][p0] <- as.character(f)
  rownames(MOUSE_orthos@assays$RNA@meta.features)[p0] <- as.character(f)}

## test renaming 
x <- length(MOUSE_orthos@assays$RNA@counts@Dimnames[[1]])
genes_genes_new1 <- c()
for (i in 1:x){
  genes_genes_new1[i] <- paste0(MOUSE_orthos@assays$RNA@counts@Dimnames[[1]][i],"-",MOUSE_orthos_copy@assays$RNA@counts@Dimnames[[1]][i])}  
genes_genes_new2 <- c()
for (i in 1:x){
  genes_genes_new2[i] <- paste0(MOUSE_orthos@assays$RNA@data@Dimnames[[1]][i],"-",MOUSE_orthos_copy@assays$RNA@data@Dimnames[[1]][i])}  
genes_genes_new3 <- c()
for (i in 1:x){
  genes_genes_new3[i] <- paste0(rownames(MOUSE_orthos@assays$RNA@meta.features)[i],"-",rownames(MOUSE_orthos_copy@assays$RNA@meta.features)[i])}  
test_vector1 <- c(which((genes_genes_new1 %in% orthos$id)))
test_vector2 <- c(which((genes_genes_new2 %in% orthos$id)))
test_vector3 <- c(which((genes_genes_new3 %in% orthos$id)))
length(test_vector1)
length(test_vector2)
length(test_vector3)

#save and convert
saveRDS(MOUSE_orthos, file = "mouseIRI_humangenenames.rds")
object <- MOUSE_orthos
object_anndata <- convertFormat(object, from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)
write_h5ad(object_anndata, filename = "mouseIRI_humangenenames.h5ad")
object <- mousekidney
object_anndata <- convertFormat(object, from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)
write_h5ad(object_anndata, filename = "mouseIRI.h5ad")

#continue



