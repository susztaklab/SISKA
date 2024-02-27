setwd(".../Envz/hdWGCNA")
library(renv)
renv::activate()
installed.packages() # to see if env is empty")

# single-cell analysis package
library(Seurat)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)

# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

# optionally enable multithreading
enableWGCNAThreads(nThreads = 8)

# load the snRNA-seq dataset and subset
seurat_obj <- readRDS('.../Atlas/objects/Atlas6.4_seurat1.rds')
seurat_obj_subset <- subset(x = seurat_obj, subset = disease == "healthy")
seurat_obj <- subset(x = seurat_obj_subset, subset = annotation_final_level1 == "PT")

#Processing
all.genes <- rownames(seurat_obj)
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)

seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))




seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "PT" # the name of the hdWGCNA experiment
)


# construct metacells  in each group # maybe modify min.cells since 4 human samples were removed 
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("species", "orig_ident"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'pca', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'species' # set the Idents of the metacell seurat object
)

# normalize metacell expression matrix:
seurat_obj <- NormalizeMetacells(seurat_obj)




seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = "human", # the name of the group of interest in the group.by column
  group.by='species', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  slot = 'data' # using normalized data
)

setwd(".../Atlas/hdWGCNA/Atlas6.6_healthy_PT_Testrun1")

# Test different soft powers:
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)

# Update the fill scale to use the color palette
pdf(file = "soft_powers.pdf",  
    width = 12, 
    height = 12) 

# assemble with patchwork
wrap_plots(plot_list, ncol=2)
dev.off()


# construct co-expression network:
seurat_obj <- ConstructNetwork(
  seurat_obj, soft_power=6,
  setDatExpr=FALSE,
  tom_name = 'human' # name of the topoligical overlap matrix written to disk
)



# Update the fill scale to use the color palette
pdf(file = "Dendo.pdf",  
    width = 12, 
    height = 6)
PlotDendrogram(seurat_obj, main='human hdWGCNA Dendrogram')
dev.off()


# need to run ScaleData first or else throws an error:
seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))

# compute all MEs in the full single-cell dataset
seurat_obj <- ModuleEigengenes(
  seurat_obj,
  group.by.vars="orig.ident"
)

saveRDS(seurat_obj, file = ".../Atlas/hdWGCNA/Atlas6.6_healthy_PT_Testrun1/seurat_obj.rds")
seurat_obj <- readRDS(".../Atlas/hdWGCNA/Atlas6.6_healthy_PT_Testrun1/seurat_obj.rds")
setwd(".../Atlas/hdWGCNA/Atlas6.6_healthy_PT_Testrun1")
# harmonized module eigengenes:
hMEs <- GetMEs(seurat_obj)

# module eigengenes:
MEs <- GetMEs(seurat_obj, harmonized=FALSE)


# compute eigengene-based connectivity (kME):
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'species', group_name = 'human'
)

# rename the modules
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "Human_PT"
)

library(dplyr)
# plot genes ranked by kME for each module
p <- PlotKMEs(seurat_obj, ncol=5)

pdf(file = "PlotKMEs.pdf",  
    width = 12, 
    height = 12)
p
dev.off()

# get the module assignment table:
modules <- GetModules(seurat_obj)

# show the first 6 columns:
head(modules[,1:6])



# get hub genes
hub_df <- GetHubGenes(seurat_obj, n_hubs = 10)

head(hub_df)

# compute gene scoring for the top 25 hub genes by kME for each module
# with Seurat method
seurat_obj <- ModuleExprScore(
  seurat_obj,
  n_genes = 25,
  method='Seurat'
)

# get hMEs from seurat object
MEs <- GetMEs(seurat_obj, harmonized=TRUE)
mods <- colnames(MEs); mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)

# plot with Seurat's DotPlot function
p <- DotPlot(seurat_obj, features=mods, group.by = 'species')

# flip the x/y axes, rotate the axis labels, and change color scheme:
p <- p +
  coord_flip() +
  RotatedAxis() +
  scale_color_gradient2(high='red', mid='grey95', low='blue')

pdf(file = "dotplot.pdf",  
    width = 24, 
    height = 12)
p
dev.off()


# Plot INH-M4 hME using Seurat VlnPlot function
p <- VlnPlot(
  seurat_obj,
  features = 'Human_PT1',
  group.by = 'species',
  pt.size = 0 # don't show actual data points
)

# add box-and-whisker plots on top:
p <- p + geom_boxplot(width=.25, fill='white')

# change axis labels and remove legend:
p <- p + xlab('') + ylab('hME') + NoLegend()

pdf(file = "vln.pdf",  
    width = 12, 
    height = 12)
p
dev.off()
# plot output













sessionInfo()
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
  [1] hdWGCNA_0.2.17        WGCNA_1.72-1          fastcluster_1.2.3    
[4] dynamicTreeCut_1.63-1 patchwork_1.1.2       cowplot_1.1.1        
[7] SeuratObject_4.1.3    Seurat_4.3.0          renv_0.16.0          

loaded via a namespace (and not attached):
  [1] backports_1.4.1        Hmisc_4.7-2            plyr_1.8.8            
[4] igraph_1.4.1           lazyeval_0.2.2         sp_1.6-0              
[7] splines_4.0.2          listenv_0.9.0          scattermore_0.8       
[10] ggplot2_3.4.2          digest_0.6.31          foreach_1.5.2         
[13] htmltools_0.5.5        GO.db_3.12.1           fansi_1.0.4           
[16] magrittr_2.0.3         checkmate_2.1.0        memoise_2.0.1         
[19] tensor_1.5             cluster_2.1.0          doParallel_1.0.17     
[22] ROCR_1.0-11            globals_0.16.2         tester_0.1.7          
[25] matrixStats_0.63.0     spatstat.sparse_3.0-1  jpeg_0.1-8.1          
[28] colorspace_2.1-0       blob_1.2.4             ggrepel_0.9.3         
[31] xfun_0.38              dplyr_1.1.1            jsonlite_1.8.4        
[34] progressr_0.13.0       spatstat.data_3.0-1    impute_1.64.0         
[37] survival_3.1-12        zoo_1.8-11             iterators_1.0.14      
[40] glue_1.6.2             polyclip_1.10-4        gtable_0.3.3          
[43] leiden_0.4.3           future.apply_1.10.0    BiocGenerics_0.36.1   
[46] abind_1.4-5            scales_1.2.1           DBI_1.1.3             
[49] spatstat.random_3.1-4  miniUI_0.1.1.1         Rcpp_1.0.10           
[52] viridisLite_0.4.1      xtable_1.8-4           htmlTable_2.4.1       
[55] reticulate_1.28        foreign_0.8-80         bit_4.0.5             
[58] preprocessCore_1.52.1  Formula_1.2-5          stats4_4.0.2          
[61] htmlwidgets_1.6.2      httr_1.4.5             RColorBrewer_1.1-3    
[64] ellipsis_0.3.2         ica_1.0-3              pkgconfig_2.0.3       
[67] nnet_7.3-14            uwot_0.1.14            deldir_1.0-6          
[70] utf8_1.2.3             tidyselect_1.2.0       rlang_1.1.0           
[73] reshape2_1.4.4         later_1.3.0            AnnotationDbi_1.52.0  
[76] munsell_0.5.0          tools_4.0.2            cachem_1.0.7          
[79] cli_3.6.1              generics_0.1.3         RSQLite_2.3.1         
[82] ggridges_0.5.4         stringr_1.5.0          fastmap_1.1.1         
[85] goftest_1.2-3          knitr_1.42             bit64_4.0.5           
[88] fitdistrplus_1.1-8     purrr_1.0.1            RANN_2.6.1            
[91] pbapply_1.7-0          future_1.32.0          nlme_3.1-148          
[94] mime_0.12              compiler_4.0.2         rstudioapi_0.14       
[97] plotly_4.10.1          png_0.1-8              spatstat.utils_3.0-2  
[100] tibble_3.2.1           stringi_1.7.12         lattice_0.20-41       
[103] Matrix_1.5-3           vctrs_0.6.1            pillar_1.9.0          
[106] lifecycle_1.0.3        spatstat.geom_3.1-0    lmtest_0.9-40         
[109] RcppAnnoy_0.0.20       data.table_1.14.8      irlba_2.3.5.1         
[112] httpuv_1.6.9           R6_2.5.1               latticeExtra_0.6-29   
[115] promises_1.2.0.1       KernSmooth_2.23-17     gridExtra_2.3         
[118] IRanges_2.24.1         parallelly_1.35.0      codetools_0.2-16      
[121] MASS_7.3-51.6          sctransform_0.3.5      harmony_0.1.1         
[124] S4Vectors_0.28.1       parallel_4.0.2         grid_4.0.2            
[127] rpart_4.1-15           tidyr_1.3.0            Rtsne_0.16            
[130] spatstat.explore_3.1-0 Biobase_2.50.0         shiny_1.7.4           
[133] base64enc_0.1-3     
