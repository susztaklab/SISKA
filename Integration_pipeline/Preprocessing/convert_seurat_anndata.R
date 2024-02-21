library(Seurat)
library(dplyr)
library(reticulate)
library(sceasy)
library(anndata)
library(patchwork)
library(SeuratDisk)

#convert
setwd(".../")
object <- readRDS("object.rds")
object_anndata <- convertFormat(object, from="seurat", to="anndata", main_layer="counts", drop_single_values=FALSE)
write_h5ad(object_anndata, filename = "object.h5ad")

sessionInfo() # PMACS output below 

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
  [1] Signac_1.3.0          SeuratDisk_0.0.0.9020 patchwork_1.1.2      
[4] anndata_0.7.5.5       sceasy_0.0.7          reticulate_1.25      
[7] dplyr_1.0.10          SeuratObject_4.1.3    Seurat_4.2.1         

loaded via a namespace (and not attached):
  [1] fastmatch_1.1-3        plyr_1.8.8             igraph_1.3.5          
[4] lazyeval_0.2.2         sp_1.5-1               splines_4.0.2         
[7] BiocParallel_1.24.1    listenv_0.8.0          scattermore_0.8       
[10] SnowballC_0.7.0        GenomeInfoDb_1.26.7    ggplot2_3.4.0         
[13] digest_0.6.30          htmltools_0.5.3        fansi_1.0.3           
[16] magrittr_2.0.3         tensor_1.5             cluster_2.1.4         
[19] ROCR_1.0-11            globals_0.16.1         Biostrings_2.58.0     
[22] matrixStats_0.62.0     docopt_0.7.1           spatstat.sparse_3.0-0 
[25] colorspace_2.0-3       rappdirs_0.3.3         ggrepel_0.9.2         
[28] sparsesvd_0.2-2        crayon_1.5.2           RCurl_1.98-1.9        
[31] jsonlite_1.8.3         progressr_0.11.0       spatstat.data_3.0-0   
[34] survival_3.4-0         zoo_1.8-11             glue_1.6.2            
[37] polyclip_1.10-4        gtable_0.3.1           zlibbioc_1.36.0       
[40] XVector_0.30.0         leiden_0.4.3           future.apply_1.10.0   
[43] BiocGenerics_0.36.1    abind_1.4-5            scales_1.2.1          
[46] DBI_1.1.3              spatstat.random_3.0-1  miniUI_0.1.1.1        
[49] Rcpp_1.0.9             viridisLite_0.4.1      xtable_1.8-4          
[52] bit_4.0.4              stats4_4.0.2           htmlwidgets_1.5.4     
[55] httr_1.4.4             RColorBrewer_1.1-3     ellipsis_0.3.2        
[58] ica_1.0-3              pkgconfig_2.0.3        farver_2.1.1          
[61] ggseqlogo_0.1          uwot_0.1.14            deldir_1.0-6          
[64] utf8_1.2.2             here_1.0.1             tidyselect_1.2.0      
[67] rlang_1.0.6            reshape2_1.4.4         later_1.3.0           
[70] munsell_0.5.0          tools_4.0.2            cli_3.4.1             
[73] generics_0.1.3         ggridges_0.5.4         stringr_1.4.1         
[76] fastmap_1.1.0          goftest_1.2-3          bit64_4.0.5           
[79] fitdistrplus_1.1-8     purrr_0.3.5            RANN_2.6.1            
[82] pbapply_1.6-0          future_1.29.0          nlme_3.1-160          
[85] mime_0.12              slam_0.1-50            RcppRoll_0.3.0        
[88] hdf5r_1.3.7            compiler_4.0.2         plotly_4.10.1         
[91] png_0.1-7              spatstat.utils_3.0-1   tibble_3.1.8          
[94] tweenr_2.0.2           stringi_1.7.8          lattice_0.20-45       
[97] Matrix_1.5-3           vctrs_0.5.1            pillar_1.8.1          
[100] lifecycle_1.0.3        spatstat.geom_3.0-3    lmtest_0.9-40         
[103] RcppAnnoy_0.0.20       data.table_1.14.4      cowplot_1.1.1         
[106] bitops_1.0-7           irlba_2.3.5.1          httpuv_1.6.6          
[109] GenomicRanges_1.42.0   R6_2.5.1               promises_1.2.0.1      
[112] KernSmooth_2.23-20     gridExtra_2.3          lsa_0.73.3            
[115] IRanges_2.24.1         parallelly_1.32.1      codetools_0.2-18      
[118] MASS_7.3-58.1          assertthat_0.2.1       rprojroot_2.0.3       
[121] withr_2.5.0            qlcMatrix_0.9.7        sctransform_0.3.5     
[124] Rsamtools_2.6.0        S4Vectors_0.28.1       GenomeInfoDbData_1.2.4
[127] parallel_4.0.2         grid_4.0.2             tidyr_1.2.1           
[130] Rtsne_0.16             spatstat.explore_3.0-5 ggforce_0.4.1         
[133] shiny_1.7.3   