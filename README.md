# Cross-Species Integrated Single-cell Kidney Atlas (SISKA)

These are the codes related to our SISKA multi-species atlas. 

### Download the data here: 

SISKA, the extended human disease-centered atlas, or any related file from our Zenodo Repository

https://doi.org/10.5281/zenodo.15007208

Go to "Integration Pipeline" to find the code used to generate SISKA and the extended human atlas.

Go to "Figures" to browse custom codes used to generate F1-F6 and supplemental figures. 

Go to https://susztaklab.com/SISKA/ and select your species and gene of interest (beta).

Go to https://susztaklab.com/rhge/ to compare gene expression coordination of cellular functions between rodent models and human diseases (beta).

Go to https://github.com/kloetzerka/CellSpectra to learn about our R package CellSpectra 

## Input datasets 

Previously published input data can be accessed via the Kidney Precision Medicine Project website https://www.kpmp.org for KPMP human snRNA-seq (UCSD_HuBMAP_KPMP-Biopsy_10X-R_12032021.h5Seurat), and Gene Expression Omnibus accession codes GSE211785 (Susztak lab human snRNA-seq), GSE209821 (ZSF1 rat snRNA-seq), GSE183839 (DOCA rat snRNA-seq), GSE184652 (diabetic kidney disease mouse snRNA-seq) and GSE139107 (mouse AKI time course snRNA-seq). Unpublished raw and processed snRNA-seq mouse data have been deposited in Gene Expression Omnibus with the accession code GSE291551. 

Please check the original publications and respective data repositories for additional information regarding the input seurat objects. For details on package versions used for basic processing see our manuscript and Supplementary Tables.

## Integration

### Preprocessing

For our basic example script of intial processing see:
SoupX_Doubletfinder_Basic.R

1.) Changing to human gene symbols based on an Ensembl gene list

Script used to convert to Seurat objects to Human 1-to-1 orthologous genes:
Seuratobject_change_gene_symbols.R

2.) Converting to anndata objects

Code and environment information used to convert objects to anndata:
convert_seurat_anndata.R

3.) Merging and Preprocessing of anndata objects + pre-scVI processing (HVG selection, etc.)

Github_Preprocessing_1_Merge.ipynb

4.) Integration Run 1

Github_Preprocessing_2_Integration1.ipynb

5.) Cleaning

Github_Preprocessing_3_Cleaning1.ipynb

6.) Preparing Integration Run 2 (final)

Github_Preprocessing_4_HVG.ipynb

### scVI

Details on our final integration model and script

### Post integration

After integration processing included:

- Annotation

- Final Cleaning (doublet/ambient contamination with inconsitent marker expression, problematic outlier samples, etc.)

### Atlas extension

Additional published input data used for the extended human atlas was accessed via CELLxGENE (diabetic kidney disease, Wilson et al. 2022), GSE185948 (ADPKD), and the KPMP website (KPMP unpublished). An additional CKD dataset was recently published (Liu et al., Science 2025). Unpublished human snRNA-seq data from biopsies has been deposited in Zenodo (https://doi.org/10.5281/zenodo.15007208).

1.) Preprocessing

Basic preprocessing to generate input objects

2.) Integration and cleaning

Follow our Step-by-step pipeline to reproduce the extended atlas.


### Figures

Details about downstream analysis to generate the figures of our manuscript.
