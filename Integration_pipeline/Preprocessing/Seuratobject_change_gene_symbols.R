#New Preprocessing for Atlas3 with new Humphrey data (1M, created from h5 files)
#28.Nov 2022

library(Seurat)
library(dplyr)
library(reticulate)
library(sceasy)
library(anndata)
library(patchwork)
library(SeuratDisk)

#import new ensembl gene list 
orthos <- read.csv(file = "/home/kloetzer/Atlas/EnsemblGeneLists/Genelist_V2_subset.csv")

#import "raw" rds files with original gene names 
## human
KPMP <- '/home/kloetzer/data/KPMP/KPMP_objectlarge.h5Seurat'
S_HUMAN <- "/home/kloetzer/data/rat_integration/Amin/HK.SN.ForSteve.8.12.2022.rds"

## mouse
BERLIN1 <- '/home/kloetzer/data/mouse_berlin/mouseberlin.rds'
BERLIN2 <- '/home/kloetzer/data/mouse_berlin/mouseberlin2.rds'
HUMPHREYS_1M <- '/home/kloetzer/data/rat_integration/mouse_Humphrey/GSE184652_RAW/Humphreys_1M_Atlas_final.rds'
HUMPHREYS_IRI <- "/home/kloetzer/data/HumphreyIRIMouse/IRI_mouse.rds"

## rat
RAT_ABEDINI <- "/home/kloetzer/data/rat_integration/to_Konstantin/Rat.All.Idents.9.12.2022.rds"
RAT_BALZER <- "/home/kloetzer/data/Michael/ZSF1_comb_soupX_diet_object.rds"

#set WD
setwd("/home/kloetzer/Atlas/objects")

#change to human gene names
mousekidney <- readRDS(file = BERLIN1)
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
saveRDS(MOUSE_orthos, file = "mouseberlin1_humangenenames.rds")
object <- MOUSE_orthos
object_anndata <- convertFormat(object, from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)
write_h5ad(object_anndata, filename = "mouseberlin1_humangenenames.h5ad")
object <- mousekidney
object_anndata <- convertFormat(object, from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)
write_h5ad(object_anndata, filename = "mouseberlin1.h5ad")

####################################
#change to human gene names
mousekidney <- readRDS(file = BERLIN2)
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
saveRDS(MOUSE_orthos, file = "mouseberlin2_humangenenames.rds")
object <- MOUSE_orthos
object_anndata <- convertFormat(object, from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)
write_h5ad(object_anndata, filename = "mouseberlin2_humangenenames.h5ad")
object <- mousekidney
object_anndata <- convertFormat(object, from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)
write_h5ad(object_anndata, filename = "mouseberlin2.h5ad")


################################
#change to human gene names
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

##################################
#change to human gene names
mousekidney <- readRDS(file = HUMPHREYS_1M)
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
saveRDS(MOUSE_orthos, file = "Humphreys_1M_humangenes.rds")
object <- MOUSE_orthos
object_anndata <- convertFormat(object, from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)
write_h5ad(object_anndata, filename = "Humphreys_1M_humangenes.h5ad")
object <- mousekidney
object_anndata <- convertFormat(object, from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)
write_h5ad(object_anndata, filename = "Humphreys_1M.h5ad")

###############################Rat##########################
object <- readRDS(file = RAT_ABEDINI)
DefaultAssay(object)<-"RNA"
ratgenes<-as.vector(unlist(orthos$Rat.gene.name))
humangenes<-as.vector(unlist(orthos$Gene.name))
object_orthos <- DietSeurat(
  object,
  counts = TRUE,
  data = TRUE,
  scale.data = FALSE,
  features = ratgenes,
  assays = NULL,
  dimreducs = NULL,
  graphs = NULL,
  misc = TRUE
)
genes <- object_orthos@assays$RNA@counts@Dimnames[[1]]
object_orthos_copy <- object_orthos

for (i in genes){
  ie <-paste0("^",i,"$")  # using regex for exact match (if not it might annotate similar gene names)
  p0 <- which(grepl(ie, genes))
  p <- which(grepl(ie,orthos$Rat.gene.name)) # find position of the match
  p1<-p[1]  # keeping the first position from p
  f<-orthos$Gene.name[p1] # get the human gene name from the same position
  object_orthos@assays$RNA@counts@Dimnames[[1]][p0] <- as.character(f)
  object_orthos@assays$RNA@data@Dimnames[[1]][p0] <- as.character(f)
  rownames(object_orthos@assays$RNA@meta.features)[p0] <- as.character(f)}

## test renaming 
x <- length(object_orthos@assays$RNA@counts@Dimnames[[1]])
genes_genes_new1 <- c()
for (i in 1:x){
  genes_genes_new1[i] <- paste0(object_orthos@assays$RNA@counts@Dimnames[[1]][i],"-",object_orthos_copy@assays$RNA@counts@Dimnames[[1]][i])}  
genes_genes_new2 <- c()
for (i in 1:x){
  genes_genes_new2[i] <- paste0(object_orthos@assays$RNA@data@Dimnames[[1]][i],"-",object_orthos_copy@assays$RNA@data@Dimnames[[1]][i])}  
genes_genes_new3 <- c()
for (i in 1:x){
  genes_genes_new3[i] <- paste0(rownames(object_orthos@assays$RNA@meta.features)[i],"-",rownames(object_orthos_copy@assays$RNA@meta.features)[i])}  
test_vector1 <- c(which((genes_genes_new1 %in% orthos$id)))
test_vector2 <- c(which((genes_genes_new2 %in% orthos$id)))
test_vector3 <- c(which((genes_genes_new3 %in% orthos$id)))
length(test_vector1)
length(test_vector2)
length(test_vector3)

saveRDS(object_orthos, file = "RAT_ABEDINI_humangenes.rds")
object_anndata <- convertFormat(object_orthos, from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)
write_h5ad(object_anndata, filename = "RAT_ABEDINI_humangenes.h5ad")
object_anndata <- convertFormat(object, from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)
write_h5ad(object_anndata, filename = "RAT_ABEDINI.h5ad")

#########################
object <- readRDS(file = RAT_BALZER)
DefaultAssay(object)<-"RNA"
ratgenes<-as.vector(unlist(orthos$Rat.gene.name))
humangenes<-as.vector(unlist(orthos$Gene.name))
object_orthos <- DietSeurat(
  object,
  counts = TRUE,
  data = TRUE,
  scale.data = FALSE,
  features = ratgenes,
  assays = NULL,
  dimreducs = NULL,
  graphs = NULL,
  misc = TRUE
)
genes <- object_orthos@assays$RNA@counts@Dimnames[[1]]
object_orthos_copy <- object_orthos

for (i in genes){
  ie <-paste0("^",i,"$")  # using regex for exact match (if not it might annotate similar gene names)
  p0 <- which(grepl(ie, genes))
  p <- which(grepl(ie,orthos$Rat.gene.name)) # find position of the match
  p1<-p[1]  # keeping the first position from p
  f<-orthos$Gene.name[p1] # get the human gene name from the same position
  object_orthos@assays$RNA@counts@Dimnames[[1]][p0] <- as.character(f)
  object_orthos@assays$RNA@data@Dimnames[[1]][p0] <- as.character(f)
  rownames(object_orthos@assays$RNA@meta.features)[p0] <- as.character(f)}

## test renaming 
x <- length(object_orthos@assays$RNA@counts@Dimnames[[1]])
genes_genes_new1 <- c()
for (i in 1:x){
  genes_genes_new1[i] <- paste0(object_orthos@assays$RNA@counts@Dimnames[[1]][i],"-",object_orthos_copy@assays$RNA@counts@Dimnames[[1]][i])}  
genes_genes_new2 <- c()
for (i in 1:x){
  genes_genes_new2[i] <- paste0(object_orthos@assays$RNA@data@Dimnames[[1]][i],"-",object_orthos_copy@assays$RNA@data@Dimnames[[1]][i])}  
genes_genes_new3 <- c()
for (i in 1:x){
  genes_genes_new3[i] <- paste0(rownames(object_orthos@assays$RNA@meta.features)[i],"-",rownames(object_orthos_copy@assays$RNA@meta.features)[i])}  
test_vector1 <- c(which((genes_genes_new1 %in% orthos$id)))
test_vector2 <- c(which((genes_genes_new2 %in% orthos$id)))
test_vector3 <- c(which((genes_genes_new3 %in% orthos$id)))
length(test_vector1)
length(test_vector2)
length(test_vector3)

saveRDS(object_orthos, file = "RAT_BALZER_humangenes.rds")
object_anndata <- convertFormat(object_orthos, from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)
write_h5ad(object_anndata, filename = "RAT_BALZER_humangenes.h5ad")
object_anndata <- convertFormat(object, from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)
write_h5ad(object_anndata, filename = "RAT_BALZER.h5ad")

####################human###########################
object <- LoadH5Seurat(file = KPMP)
DefaultAssay(object)<-"RNA"
humangenes<-as.vector(unlist(orthos$Gene.name))
object_orthos <- DietSeurat(
  object,
  counts = TRUE,
  data = TRUE,
  scale.data = FALSE,
  features = humangenes,
  assays = NULL,
  dimreducs = NULL,
  graphs = NULL,
  misc = TRUE
)

saveRDS(object_orthos, file = "KPMP_humangenes_orthos.rds")
object_anndata <- convertFormat(object_orthos, from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)
write_h5ad(object_anndata, filename = "KPMP_humangenes_orthos.h5ad")
object_anndata <- convertFormat(object, from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)
write_h5ad(object_anndata, filename = "KPMP.h5ad")

#################################################
object <- readRDS(file = S_HUMAN)
DefaultAssay(object)<-"RNA"
humangenes<-as.vector(unlist(orthos$Gene.name))
object_orthos <- DietSeurat(
  object,
  counts = TRUE,
  data = TRUE,
  scale.data = FALSE,
  features = humangenes,
  assays = NULL,
  dimreducs = NULL,
  graphs = NULL,
  misc = TRUE
)

saveRDS(object_orthos, file = "Susztak_HUMAN_humangenes_orthos.rds")
object_anndata <- convertFormat(object_orthos, from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)
write_h5ad(object_anndata, filename = "Susztak_HUMAN_humangenes_orthos.h5ad")
object_anndata <- convertFormat(object, from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)
write_h5ad(object_anndata, filename = "Susztak_HUMAN.h5ad")

sessionInfo()