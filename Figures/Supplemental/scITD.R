###start here

setwd(".../Envz/scITD")
library(renv)
renv::activate()
library(scITD)
library(Matrix)
library(graphics)
library(dplyr)
##

meta <- read.csv(".../Atlas/Tables/Atlas_6.6_Metadata.csv")
rownames(meta) <- meta$X
meta <- meta[, -1]
colnames(meta)[colnames(meta) == "orig_ident"] <- "donors"
colnames(meta)[colnames(meta) == "annotation_final_level1"] <- "ctypes" 
counts <- readRDS('.../Atlas/Tables/Atlas6.4_matrix.rds')

#should be as factor 
meta <- meta %>% mutate_all(as.factor)




# 6.4
# set names for plot
setwd(".../Atlas/scITD/Atlas6.4_disease_injured")
PLOT1 <- "plotA_1.pdf"
PLOT2 <- "plotA_2.pdf"
PLOT3 <- "plotA_3.pdf"
PLOT4 <- "plotA_factor_control.pdf"

#we don't use prolif_tubule to not lose too many donors 
param_list <- initialize_params(ctypes_use = c('PT', 'TAL_MD', 'DCT_CNT_CD', 'DTL_ATL', 'PEC', 'Podo'),
                                ncores = 30, rand_seed = 10)

# create project container
pbmc_container <- make_new_container(count_data=counts, 
                                     meta_data=meta,
                                     params=param_list,
                                     label_donor_sex = FALSE)

# form the tensor from the data
pbmc_container <- form_tensor(pbmc_container, donor_min_cells=0,
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=.1,
                              scale_var = TRUE, var_scale_power = 2)

# number of genes included in the tensor
print(length(pbmc_container[["all_vargenes"]]))

# run the tensor decomposition
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(5,10),
                                 tucker_type = 'regular', rotation_type = 'hybrid')
#try sep of species - check if factors

# get donor scores-metadata associations
pbmc_container <- get_meta_associations(pbmc_container, vars_test=c('species', 'disease'), stat_use='pval')

# plot donor scores
pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('species', 'disease'),
                                    show_donor_ids = FALSE,
                                    add_meta_associations='pval')


pdf(PLOT1, width = 20, height = 30)
# show the donor scores heatmap
pbmc_container$plots$donor_matrix
dev.off()

# get significant genes
pbmc_container <- get_lm_pvals(pbmc_container)
# generate the loadings plots
pbmc_container <- get_all_lds_factor_plots(pbmc_container, 
                                           use_sig_only=TRUE,
                                           nonsig_to_zero=TRUE,
                                           sig_thresh=.02,
                                           display_genes=FALSE,
                                           gene_callouts = TRUE,
                                           callout_n_gene_per_ctype=3,
                                           show_var_explained = TRUE)

# arrange the plots into a figure and show the figure
myfig <- render_multi_plots(pbmc_container,data_type='loadings')
pdf(PLOT2, width = 40, height = 24)
# show the donor scores heatmap
myfig
dev.off()

############################# set factor ######################################
pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=4, method="fgsea", thresh=0.01, db_use=c("GO"), signed=TRUE, draw_plot = FALSE)
pdf(PLOT3, width = 20, height = 40)
pbmc_container$plots$gsea$`4` #set factor to plot
dev.off()

#full donor scores matrix 
ds_matrix <- pbmc_container[["tucker_results"]][[1]]


#unfolded loadings tensor
uf_loadings <- pbmc_container[["tucker_results"]][[2]]

write.csv(ds_matrix, file = "ds_matrix_A.csv")
write.csv(uf_loadings, file = "uf_loadings_A.csv")

#F1
f1_data <- get_one_factor(pbmc_container, factor_select=1)
f1_dscores <- f1_data[[1]]
f1_loadings <- f1_data[[2]]
write.csv(f1_dscores, file = "A_F1_dscores.csv")
write.csv(f1_loadings, file = "A_F1_loadings.csv")

#F4
f1_data <- get_one_factor(pbmc_container, factor_select=4)
f1_dscores <- f1_data[[1]]
f1_loadings <- f1_data[[2]]

write.csv(f1_dscores, file = "A_F4_dscores.csv")
write.csv(f1_loadings, file = "A_F4_loadings.csv")

matrix_pvalue_5 <- pbmc_container[["plots"]][["all_lds_plots"]][["4"]]@matrix
write.csv(matrix_pvalue_5, file = "A_F4_heatmap_matrix_pvalue.csv")

matrix_pvalue_5 <- pbmc_container[["plots"]][["all_lds_plots"]][["1"]]@matrix
write.csv(matrix_pvalue_5, file = "A_F1_heatmap_matrix_pvalue.csv")






library(dplyr)
#Back to local
#we created a matrix with loadings. Non significant genes are set to "0". We now want to know how many 
#0 each celltype has 
factor4_matrix <- read.csv(".../Atlas/scITD/Atlas6.4_epithelial_injury/A_F4_heatmap_matrix_pvalue.csv", 
                           row.names = 1)
head(factor4_matrix)

df <- factor4_matrix

df <- df %>% select("PT", "DCT_CNT_CD", "TAL_MD")

head(df)

# Remove rows (genes) with at least one '0'
df <- df %>% filter_all(all_vars(. != 0))
head(df)

# Calculate mean per row
df$mean <- rowMeans(df)
head(df)


write.csv(df, file = ".../Atlas/scITD/Atlas6.4_epithelial_injury/A_F4_heatmap_matrix_modified.csv")






















###start here
setwd(".../Envz/scITD")
library(renv)
renv::activate()
library(scITD)
library(Matrix)
library(graphics)
library(dplyr)

##
meta <- read.csv(".../Atlas/Tables/Atlas_6.4_healthy_scITD_metadata.csv")
rownames(meta) <- meta$X
meta <- meta[, -1]
colnames(meta)[colnames(meta) == "orig_ident"] <- "donors"
colnames(meta)[colnames(meta) == "annotation_final_level1"] <- "ctypes" 
counts <- readRDS('.../Atlas/Tables/Atlas6.4_healthy_matrix.rds')

#should be as factor 
meta <- meta %>% mutate_all(as.factor)

# set names for plot
setwd(".../Atlas/scITD/Atlas6.4_healthy")
PLOT1 <- "plotA_1.pdf"
PLOT2 <- "plotA_2.pdf"
PLOT3 <- "plotA_3.pdf"
PLOT4 <- "plotA_factor_control.pdf"

#we don't use prolif_tubule to not lose too many donors 
param_list <- initialize_params(ctypes_use = c('TAL_MD', 'DCT_CNT_CD', 'PT', 'IC', 'EC', 'Stromal', 'Immune', 'PEC', 'DTL_ATL', 'Podo'),
                                ncores = 30, rand_seed = 10)
# create project container
pbmc_container <- make_new_container(count_data=counts, 
                                     meta_data=meta,
                                     params=param_list,
                                     label_donor_sex = FALSE)

# form the tensor from the data
pbmc_container <- form_tensor(pbmc_container, donor_min_cells=0,
                              norm_method='trim', scale_factor=10000,
                              vargenes_method='norm_var_pvals', vargenes_thresh=.1,
                              scale_var = TRUE, var_scale_power = 2)

# number of genes included in the tensor
print(length(pbmc_container[["all_vargenes"]]))

# run the tensor decomposition
pbmc_container <- run_tucker_ica(pbmc_container, ranks=c(5,10),
                                 tucker_type = 'regular', rotation_type = 'hybrid')


# get donor scores-metadata associations
pbmc_container <- get_meta_associations(pbmc_container, vars_test=c('species', 'sex'), stat_use='pval')

# plot donor scores
pbmc_container <- plot_donor_matrix(pbmc_container, meta_vars=c('species', 'sex'),
                                    show_donor_ids = FALSE,
                                    add_meta_associations='pval')

pdf(PLOT1, width = 20, height = 30)
# show the donor scores heatmap
pbmc_container$plots$donor_matrix
dev.off()

# get significant genes
pbmc_container <- get_lm_pvals(pbmc_container)

# generate the loadings plots
pbmc_container <- get_all_lds_factor_plots(pbmc_container, 
                                           use_sig_only=TRUE,
                                           nonsig_to_zero=TRUE,
                                           sig_thresh=.02,
                                           display_genes=FALSE,
                                           gene_callouts = TRUE,
                                           callout_n_gene_per_ctype=3,
                                           show_var_explained = TRUE)

# arrange the plots into a figure and show the figure
myfig <- render_multi_plots(pbmc_container,data_type='loadings')
pdf(PLOT2, width = 40, height = 24)

# show the donor scores heatmap
myfig
dev.off()

############################# set factor ######################################
pbmc_container <- run_gsea_one_factor(pbmc_container, factor_select=1, method="fgsea", thresh=0.01, db_use=c("GO"), signed=TRUE, draw_plot = FALSE)
pdf(PLOT3, width = 20, height = 40)
pbmc_container$plots$gsea$`1` #set factor to plot
dev.off()

#full donor scores matrix 
ds_matrix <- pbmc_container[["tucker_results"]][[1]]


#unfolded loadings tensor
uf_loadings <- pbmc_container[["tucker_results"]][[2]]

write.csv(ds_matrix, file = "ds_matrix_A.csv")
write.csv(uf_loadings, file = "uf_loadings_A.csv")

#F1
f1_data <- get_one_factor(pbmc_container, factor_select=1)
f1_dscores <- f1_data[[1]]
f1_loadings <- f1_data[[2]]
write.csv(f1_dscores, file = "A_F1_dscores.csv")
write.csv(f1_loadings, file = "A_F1_loadings.csv")

#F2
f1_data <- get_one_factor(pbmc_container, factor_select=2)
f1_dscores <- f1_data[[1]]
f1_loadings <- f1_data[[2]]

write.csv(f1_dscores, file = "A_F2_dscores.csv")
write.csv(f1_loadings, file = "A_F2_loadings.csv")

#F1 significant genes 
matrix_pvalue_5 <- pbmc_container[["plots"]][["all_lds_plots"]][["1"]]@matrix
write.csv(matrix_pvalue_5, file = "A_F1_heatmap_matrix_pvalue.csv")

#F2 significant genes 
matrix_pvalue_5 <- pbmc_container[["plots"]][["all_lds_plots"]][["2"]]@matrix
write.csv(matrix_pvalue_5, file = "A_F2_heatmap_matrix_pvalue.csv")

library(dplyr)
#Back to local
#we created a matrix with loadings. Non significant genes are set to "0". We now want to know how many 
#0 each celltype has 
factor1_matrix <- read.csv(".../Atlas/scITD/Atlas6.4_healthy/A_F1_heatmap_matrix_pvalue.csv")

counts <- colSums(factor1_matrix != 0)
counts
factor1_matrix <- rbind(factor1_matrix, counts)
factor1_matrix[nrow(factor1_matrix),1] <- "Number of significant genes"

write.csv(factor1_matrix, file = ".../Atlas/scITD/Atlas6.4_healthy/A_F1_heatmap_matrix_pvalue_numbersig.csv")

