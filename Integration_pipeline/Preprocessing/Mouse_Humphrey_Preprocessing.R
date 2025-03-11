#load and use Humphrey 1M data as h5 file
#28. Nov 2022

#$ module load hdf5/1.8.21

library(dplyr)
library(Seurat)
library(anndata)

#set wd and meta
setwd("/home/kloetzer/data/rat_integration/mouse_Humphrey")
meta_data <- read.csv(file = 'GSE184652_metadata.csv')
rownames(meta_data) <- meta_data$X
meta_data <- meta_data[,-1]
setwd("/home/kloetzer/data/rat_integration/mouse_Humphrey/GSE184652_RAW")

file1 <- CreateSeuratObject(counts = Read10X_h5("GSM5594468_E3019_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "E3019", min.cells = 3, min.features = 200)
file1 <- RenameCells(file1, add.cell.id = "E3019")
Idents(file1) <- file1@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "E3019")
file1 <- subset(x = file1, idents = rownames(meta))
file1 <- AddMetaData(file1, meta)

file2 <- CreateSeuratObject(counts = Read10X_h5("GSM5594469_A3020_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "A3020", min.cells = 3, min.features = 200)
file2 <- RenameCells(file2, add.cell.id = "A3020")
Idents(file2) <- file2@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "A3020")
file2 <- subset(x = file2, idents = rownames(meta))
file2 <- AddMetaData(file2, meta)

file3 <- CreateSeuratObject(counts = Read10X_h5("GSM5594470_F3021_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "F3021", min.cells = 3, min.features = 200)
file3 <- RenameCells(file3, add.cell.id = "F3021")
Idents(file3) <- file3@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "F3021")
file3 <- subset(x = file3, idents = rownames(meta))
file3 <- AddMetaData(file3, meta)

file4 <- CreateSeuratObject(counts = Read10X_h5("GSM5594471_G3022_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "G3022", min.cells = 3, min.features = 200)
file4 <- RenameCells(file4, add.cell.id = "G3022")
Idents(file4) <- file4@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "G3022")
file4 <- subset(x = file4, idents = rownames(meta))
file4 <- AddMetaData(file4, meta)

file5 <- CreateSeuratObject(counts = Read10X_h5("GSM5594472_H3023_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "H3023", min.cells = 3, min.features = 200)
file5 <- RenameCells(file5, add.cell.id = "H3023")
Idents(file5) <- file5@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "H3023")
file5 <- subset(x = file5, idents = rownames(meta))
file5 <- AddMetaData(file5, meta)

file6 <- CreateSeuratObject(counts = Read10X_h5("GSM5594473_N3024_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "N3024", min.cells = 3, min.features = 200)
file6 <- RenameCells(file6, add.cell.id = "N3024")
Idents(file6) <- file6@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "N3024")
file6 <- subset(x = file6, idents = rownames(meta))
file6 <- AddMetaData(file6, meta)

file7 <- CreateSeuratObject(counts = Read10X_h5("GSM5594474_B3025_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "B3025", min.cells = 3, min.features = 200)
file7 <- RenameCells(file7, add.cell.id = "B3025")
Idents(file7) <- file7@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "B3025")
file7 <- subset(x = file7, idents = rownames(meta))
file7 <- AddMetaData(file7, meta)

file8 <- CreateSeuratObject(counts = Read10X_h5("GSM5594475_A3026_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "A3026", min.cells = 3, min.features = 200)
file8 <- RenameCells(file8, add.cell.id = "A3026")
Idents(file8) <- file8@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "A3026")
file8 <- subset(x = file8, idents = rownames(meta))
file8 <- AddMetaData(file8, meta)

file9 <- CreateSeuratObject(counts = Read10X_h5("GSM5594476_B3027_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "B3027", min.cells = 3, min.features = 200)
file9 <- RenameCells(file9, add.cell.id = "B3027")
Idents(file9) <- file9@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "B3027")
file9 <- subset(x = file9, idents = rownames(meta))
file9 <- AddMetaData(file9, meta)

file10 <- CreateSeuratObject(counts = Read10X_h5("GSM5594477_C3028_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "C3028", min.cells = 3, min.features = 200)
file10 <- RenameCells(file10, add.cell.id = "C3028")
Idents(file10) <- file10@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "C3028")
file10 <- subset(x = file10, idents = rownames(meta))
file10 <- AddMetaData(file10, meta)


Atlas1 <- merge(file1, y = c(file2, file3,
                             file4, file5,
                             file6, file7,
                             file8, file9,
                             file10))

saveRDS(Atlas1, file = 'Mouse_Atlas_1_7_subset.rds')


file11 <- CreateSeuratObject(counts = Read10X_h5("GSM5594478_A3029_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "A3029", min.cells = 3, min.features = 200)
file11 <- RenameCells(file11, add.cell.id = "A3029")
Idents(file11) <- file11@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "A3029")
file11 <- subset(x = file11, idents = rownames(meta))
file11 <- AddMetaData(file11, meta)

file12 <- CreateSeuratObject(counts = Read10X_h5("GSM5594479_I3030_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "I3030", min.cells = 3, min.features = 200)
file12 <- RenameCells(file12, add.cell.id = "I3030")
Idents(file12) <- file12@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "I3030")
file12 <- subset(x = file12, idents = rownames(meta))
file12 <- AddMetaData(file12, meta)

file13 <- CreateSeuratObject(counts = Read10X_h5("GSM5594480_C3031_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "C3031", min.cells = 3, min.features = 200)
file13 <- RenameCells(file13, add.cell.id = "C3031")
Idents(file13) <- file13@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "C3031")
file13 <- subset(x = file13, idents = rownames(meta))
file13 <- AddMetaData(file13, meta)

file14 <- CreateSeuratObject(counts = Read10X_h5("GSM5594481_D3032_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "D3032", min.cells = 3, min.features = 200)
file14 <- RenameCells(file14, add.cell.id = "D3032")
Idents(file14) <- file14@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "D3032")
file14 <- subset(x = file14, idents = rownames(meta))
file14 <- AddMetaData(file14, meta)

file15 <- CreateSeuratObject(counts = Read10X_h5("GSM5594482_A3033_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "A3033", min.cells = 3, min.features = 200)
file15 <- RenameCells(file15, add.cell.id = "A3033")
Idents(file15) <- file15@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "A3033")
file15 <- subset(x = file15, idents = rownames(meta))
file15 <- AddMetaData(file15, meta)

file16 <- CreateSeuratObject(counts = Read10X_h5("GSM5594483_B3034_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "B3034", min.cells = 3, min.features = 200)
file16 <- RenameCells(file16, add.cell.id = "B3034")
Idents(file16) <- file16@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "B3034")
file16 <- subset(x = file16, idents = rownames(meta))
file16 <- AddMetaData(file16, meta)

file17 <- CreateSeuratObject(counts = Read10X_h5("GSM5594484_B3035_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "B3035", min.cells = 3, min.features = 200)
file17 <- RenameCells(file17, add.cell.id = "B3035")
Idents(file17) <- file17@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "B3035")
file17 <- subset(x = file17, idents = rownames(meta))
file17 <- AddMetaData(file17, meta)

file18 <- CreateSeuratObject(counts = Read10X_h5("GSM5594485_J3036_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "J3036", min.cells = 3, min.features = 200)
file18 <- RenameCells(file18, add.cell.id = "J3036")
Idents(file18) <- file18@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "J3036")
file18 <- subset(x = file18, idents = rownames(meta))
file18 <- AddMetaData(file18, meta)

file19 <- CreateSeuratObject(counts = Read10X_h5("GSM5594486_K3037_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "K3037", min.cells = 3, min.features = 200)
file19 <- RenameCells(file19, add.cell.id = "K3037")
Idents(file19) <- file19@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "K3037")
file19 <- subset(x = file19, idents = rownames(meta))
file19 <- AddMetaData(file19, meta)

file20 <- CreateSeuratObject(counts = Read10X_h5("GSM5594487_A3038_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "A3038", min.cells = 3, min.features = 200)
file20 <- RenameCells(file20, add.cell.id = "A3038")
Idents(file20) <- file20@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "A3038")
file20 <- subset(x = file20, idents = rownames(meta))
file20 <- AddMetaData(file20, meta)


Atlas2 <- merge(file11, y = c(file12, file13, file14, file15, file16,
                              file17, file18, file19, file20))

#saveRDS(Atlas2, file = 'Mouse_Atlas_Prototype_new2.rds')

saveRDS(Atlas2, file = 'Mouse_Atlas_2_7_subset.rds')


file21 <- CreateSeuratObject(counts = Read10X_h5("GSM5594488_B3039_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "B3039", min.cells = 3, min.features = 200)
file21 <- RenameCells(file21, add.cell.id = "B3039")
Idents(file21) <- file21@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "B3039")
file21 <- subset(x = file21, idents = rownames(meta))
file21 <- AddMetaData(file21, meta)

file22 <- CreateSeuratObject(counts = Read10X_h5("GSM5594489_D3060_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "D3060", min.cells = 3, min.features = 200)
file22 <- RenameCells(file22, add.cell.id = "D3060")
Idents(file22) <- file22@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "D3060")
file22 <- subset(x = file22, idents = rownames(meta))
file22 <- AddMetaData(file22, meta)

file23 <- CreateSeuratObject(counts = Read10X_h5("GSM5594490_I3061_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "I3061", min.cells = 3, min.features = 200)
file23 <- RenameCells(file23, add.cell.id = "I3061")
Idents(file23) <- file23@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "I3061")
file23 <- subset(x = file23, idents = rownames(meta))
file23 <- AddMetaData(file23, meta)

file24 <- CreateSeuratObject(counts = Read10X_h5("GSM5594491_J3062_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "J3062", min.cells = 3, min.features = 200)
file24 <- RenameCells(file24, add.cell.id = "J3062")
Idents(file24) <- file24@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "J3062")
file24 <- subset(x = file24, idents = rownames(meta))
file24 <- AddMetaData(file24, meta)

file25 <- CreateSeuratObject(counts = Read10X_h5("GSM5594492_G3063_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "G3063", min.cells = 3, min.features = 200)
file25 <- RenameCells(file25, add.cell.id = "G3063")
Idents(file25) <- file25@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "G3063")
file25 <- subset(x = file25, idents = rownames(meta))
file25 <- AddMetaData(file25, meta)

file26 <- CreateSeuratObject(counts = Read10X_h5("GSM5594493_H3064_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "H3064", min.cells = 3, min.features = 200)
file26 <- RenameCells(file26, add.cell.id = "H3064")
Idents(file26) <- file26@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "H3064")
file26 <- subset(x = file26, idents = rownames(meta))
file26 <- AddMetaData(file26, meta)

file27 <- CreateSeuratObject(counts = Read10X_h5("GSM5594494_I3065_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "I3065", min.cells = 3, min.features = 200)
file27 <- RenameCells(file27, add.cell.id = "I3065")
Idents(file27) <- file27@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "I3065")
file27 <- subset(x = file27, idents = rownames(meta))
file27 <- AddMetaData(file27, meta)

file28 <- CreateSeuratObject(counts = Read10X_h5("GSM5594495_P3066_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "P3066", min.cells = 3, min.features = 200)
file28 <- RenameCells(file28, add.cell.id = "P3066")
Idents(file28) <- file28@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "P3066")
file28 <- subset(x = file28, idents = rownames(meta))
file28 <- AddMetaData(file28, meta)

file29 <- CreateSeuratObject(counts = Read10X_h5("GSM5594496_Q3067_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "Q3067", min.cells = 3, min.features = 200)
file29 <- RenameCells(file29, add.cell.id = "Q3067")
Idents(file29) <- file29@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "Q3067")
file29 <- subset(x = file29, idents = rownames(meta))
file29 <- AddMetaData(file29, meta)

file30 <- CreateSeuratObject(counts = Read10X_h5("GSM5594497_G3068_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "G3068", min.cells = 3, min.features = 200)
file30 <- RenameCells(file30, add.cell.id = "G3068")
Idents(file30) <- file30@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "G3068")
file30 <- subset(x = file30, idents = rownames(meta))
file30 <- AddMetaData(file30, meta)


Atlas3 <- merge(file21, y = c(file22, file23, file24, file25, file26,
                              file27, file28, file29, file30))
saveRDS(Atlas3, file = 'Mouse_Atlas_3_7_subset.rds')


file31 <- CreateSeuratObject(counts = Read10X_h5("GSM5594498_H3069_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "H3069", min.cells = 3, min.features = 200)
file31 <- RenameCells(file31, add.cell.id = "H3069")
Idents(file31) <- file31@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "H3069")
file31 <- subset(x = file31, idents = rownames(meta))
file31 <- AddMetaData(file31, meta)

file32 <- CreateSeuratObject(counts = Read10X_h5("GSM5594499_E3070_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "E3070", min.cells = 3, min.features = 200)
file32 <- RenameCells(file32, add.cell.id = "E3070")
Idents(file32) <- file32@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "E3070")
file32 <- subset(x = file32, idents = rownames(meta))
file32 <- AddMetaData(file32, meta)

file33 <- CreateSeuratObject(counts = Read10X_h5("GSM5594500_K3071_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "K3071", min.cells = 3, min.features = 200)
file33 <- RenameCells(file33, add.cell.id = "K3071")
Idents(file33) <- file33@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "K3071")
file33 <- subset(x = file33, idents = rownames(meta))
file33 <- AddMetaData(file33, meta)

file34 <- CreateSeuratObject(counts = Read10X_h5("GSM5594501_L3072_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "L3072", min.cells = 3, min.features = 200)
file34 <- RenameCells(file34, add.cell.id = "L3072")
Idents(file34) <- file34@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "L3072")
file34 <- subset(x = file34, idents = rownames(meta))
file34 <- AddMetaData(file34, meta)

file35 <- CreateSeuratObject(counts = Read10X_h5("GSM5594502_I3073_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "I3073", min.cells = 3, min.features = 200)
file35 <- RenameCells(file35, add.cell.id = "I3073")
Idents(file35) <- file35@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "I3073")
file35 <- subset(x = file35, idents = rownames(meta))
file35 <- AddMetaData(file35, meta)

file36 <- CreateSeuratObject(counts = Read10X_h5("GSM5594503_J3074_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "J3074", min.cells = 3, min.features = 200)
file36 <- RenameCells(file36, add.cell.id = "J3074")
Idents(file36) <- file36@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "J3074")
file36 <- subset(x = file36, idents = rownames(meta))
file36 <- AddMetaData(file36, meta)

file37 <- CreateSeuratObject(counts = Read10X_h5("GSM5594504_E3075_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "E3075", min.cells = 3, min.features = 200)
file37 <- RenameCells(file37, add.cell.id = "E3075")
Idents(file37) <- file37@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "E3075")
file37 <- subset(x = file37, idents = rownames(meta))
file37 <- AddMetaData(file37, meta)

file38 <- CreateSeuratObject(counts = Read10X_h5("GSM5594505_K3076_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "K3076", min.cells = 3, min.features = 200)
file38 <- RenameCells(file38, add.cell.id = "K3076")
Idents(file38) <- file38@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "K3076")
file38 <- subset(x = file38, idents = rownames(meta))
file38 <- AddMetaData(file38, meta)

file39 <- CreateSeuratObject(counts = Read10X_h5("GSM5594506_R3077_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "R3077", min.cells = 3, min.features = 200)
file39 <- RenameCells(file39, add.cell.id = "R3077")
Idents(file39) <- file39@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "R3077")
file39 <- subset(x = file39, idents = rownames(meta))
file39 <- AddMetaData(file39, meta)

file40 <- CreateSeuratObject(counts = Read10X_h5("GSM5594507_I3078_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "I3078", min.cells = 3, min.features = 200)
file40 <- RenameCells(file40, add.cell.id = "I3078")
Idents(file40) <- file40@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "I3078")
file40 <- subset(x = file40, idents = rownames(meta))
file40 <- AddMetaData(file40, meta)


Atlas4 <- merge(file31, y = c(file32, file33, file34, file35, file36, file37,
                              file38, file39, file40))
saveRDS(Atlas4, file = 'Mouse_Atlas_4_7_subset.rds')

file41 <- CreateSeuratObject(counts = Read10X_h5("GSM5594508_J3079_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "J3079", min.cells = 3, min.features = 200)
file41 <- RenameCells(file41, add.cell.id = "J3079")
Idents(file41) <- file41@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "J3079")
file41 <- subset(x = file41, idents = rownames(meta))
file41 <- AddMetaData(file41, meta)

file42 <- CreateSeuratObject(counts = Read10X_h5("GSM5594509_F3080_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "F3080", min.cells = 3, min.features = 200)
file42 <- RenameCells(file42, add.cell.id = "F3080")
Idents(file42) <- file42@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "F3080")
file42 <- subset(x = file42, idents = rownames(meta))
file42 <- AddMetaData(file42, meta)

file43 <- CreateSeuratObject(counts = Read10X_h5("GSM5594510_M3081_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "M3081", min.cells = 3, min.features = 200)
file43 <- RenameCells(file43, add.cell.id = "M3081")
Idents(file43) <- file43@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "M3081")
file43 <- subset(x = file43, idents = rownames(meta))
file43 <- AddMetaData(file43, meta)

file44 <- CreateSeuratObject(counts = Read10X_h5("GSM5594511_N3082_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "N3082", min.cells = 3, min.features = 200)
file44 <- RenameCells(file44, add.cell.id = "N3082")
Idents(file44) <- file44@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "N3082")
file44 <- subset(x = file44, idents = rownames(meta))
file44 <- AddMetaData(file44, meta)

file45 <- CreateSeuratObject(counts = Read10X_h5("GSM5594512_L3083_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "L3083", min.cells = 3, min.features = 200)
file45 <- RenameCells(file45, add.cell.id = "L3083")
Idents(file45) <- file45@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "L3083")
file45 <- subset(x = file45, idents = rownames(meta))
file45 <- AddMetaData(file45, meta)

file46 <- CreateSeuratObject(counts = Read10X_h5("GSM5594513_K3084_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "K3084", min.cells = 3, min.features = 200)
file46 <- RenameCells(file46, add.cell.id = "K3084")
Idents(file46) <- file46@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "K3084")
file46 <- subset(x = file46, idents = rownames(meta))
file46 <- AddMetaData(file46, meta)

file47 <- CreateSeuratObject(counts = Read10X_h5("GSM5594514_F3085_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "F3085", min.cells = 3, min.features = 200)
file47 <- RenameCells(file47, add.cell.id = "F3085")
Idents(file47) <- file47@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "F3085")
file47 <- subset(x = file47, idents = rownames(meta))
file47 <- AddMetaData(file47, meta)

file48 <- CreateSeuratObject(counts = Read10X_h5("GSM5594515_L3086_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "L3086", min.cells = 3, min.features = 200)
file48 <- RenameCells(file48, add.cell.id = "L3086")
Idents(file48) <- file48@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "L3086")
file48 <- subset(x = file48, idents = rownames(meta))
file48 <- AddMetaData(file48, meta)

file49 <- CreateSeuratObject(counts = Read10X_h5("GSM5594516_K3087_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "K3087", min.cells = 3, min.features = 200)
file49 <- RenameCells(file49, add.cell.id = "K3087")
Idents(file49) <- file49@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "K3087")
file49 <- subset(x = file49, idents = rownames(meta))
file49 <- AddMetaData(file49, meta)

file50 <- CreateSeuratObject(counts = Read10X_h5("GSM5594517_L3089_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "L3089", min.cells = 3, min.features = 200)
file50 <- RenameCells(file50, add.cell.id = "L3089")
Idents(file50) <- file50@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "L3089")
file50 <- subset(x = file50, idents = rownames(meta))
file50 <- AddMetaData(file50, meta)


Atlas5 <- merge(file41, y = c(file42, file43, file44,
                              file45, file46, file47, file48, file49, file50))
saveRDS(Atlas5, file = 'Mouse_Atlas_5_7_subset.rds')

file51 <- CreateSeuratObject(counts = Read10X_h5("GSM5594518_G3090_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "G3090", min.cells = 3, min.features = 200)
file51 <- RenameCells(file51, add.cell.id = "G3090")
Idents(file51) <- file51@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "G3090")
file51 <- subset(x = file51, idents = rownames(meta))
file51 <- AddMetaData(file51, meta)

file52 <- CreateSeuratObject(counts = Read10X_h5("GSM5594519_O3091_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "O3091", min.cells = 3, min.features = 200)
file52 <- RenameCells(file52, add.cell.id = "O3091")
Idents(file52) <- file52@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "O3091")
file52 <- subset(x = file52, idents = rownames(meta))
file52 <- AddMetaData(file52, meta)

file53 <- CreateSeuratObject(counts = Read10X_h5("GSM5594520_P3092_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "P3092", min.cells = 3, min.features = 200)
file53 <- RenameCells(file53, add.cell.id = "P3092")
Idents(file53) <- file53@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "P3092")
file53 <- subset(x = file53, idents = rownames(meta))
file53 <- AddMetaData(file53, meta)

file54 <- CreateSeuratObject(counts = Read10X_h5("GSM5594521_M3093_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "M3093", min.cells = 3, min.features = 200)
file54 <- RenameCells(file54, add.cell.id = "M3093")
Idents(file54) <- file54@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "M3093")
file54 <- subset(x = file54, idents = rownames(meta))
file54 <- AddMetaData(file54, meta)

file55 <- CreateSeuratObject(counts = Read10X_h5("GSM5594522_N3094_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "N3094", min.cells = 3, min.features = 200)
file55 <- RenameCells(file55, add.cell.id = "N3094")
Idents(file55) <- file55@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "N3094")
file55 <- subset(x = file55, idents = rownames(meta))
file55 <- AddMetaData(file55, meta)

file56 <- CreateSeuratObject(counts = Read10X_h5("GSM5594523_I3095_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "I3095", min.cells = 3, min.features = 200)
file56 <- RenameCells(file56, add.cell.id = "I3095")
Idents(file56) <- file56@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "I3095")
file56 <- subset(x = file56, idents = rownames(meta))
file56 <- AddMetaData(file56, meta)

file57 <- CreateSeuratObject(counts = Read10X_h5("GSM5594524_M3096_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "M3096", min.cells = 3, min.features = 200)
file57 <- RenameCells(file57, add.cell.id = "M3096")
Idents(file57) <- file57@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "M3096")
file57 <- subset(x = file57, idents = rownames(meta))
file57 <- AddMetaData(file57, meta)

file58 <- CreateSeuratObject(counts = Read10X_h5("GSM5594525_T3097_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "T3097", min.cells = 3, min.features = 200)
file58 <- RenameCells(file58, add.cell.id = "T3097")
Idents(file58) <- file58@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "T3097")
file58 <- subset(x = file58, idents = rownames(meta))
file58 <- AddMetaData(file58, meta)

file59 <- CreateSeuratObject(counts = Read10X_h5("GSM5594526_M3098_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "M3098", min.cells = 3, min.features = 200)
file59 <- RenameCells(file59, add.cell.id = "M3098")
Idents(file59) <- file59@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "M3098")
file59 <- subset(x = file59, idents = rownames(meta))
file59 <- AddMetaData(file59, meta)

file60 <- CreateSeuratObject(counts = Read10X_h5("GSM5594527_N3099_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "N3099", min.cells = 3, min.features = 200)
file60 <- RenameCells(file60, add.cell.id = "N3099")
Idents(file60) <- file60@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "N3099")
file60 <- subset(x = file60, idents = rownames(meta))
file60 <- AddMetaData(file60, meta)

Atlas6 <- merge(file51, y = c(file52, file53, file54, file55, file56, file57, file58,
                              file59, file60))
saveRDS(Atlas6, file = 'Mouse_Atlas_6_7_subset.rds')

file61 <- CreateSeuratObject(counts = Read10X_h5("GSM5594528_A3100_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "A3100", min.cells = 3, min.features = 200)
file61 <- RenameCells(file61, add.cell.id = "A3100")
Idents(file61) <- file61@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "A3100")
file61 <- subset(x = file61, idents = rownames(meta))
file61 <- AddMetaData(file61, meta)

file62 <- CreateSeuratObject(counts = Read10X_h5("GSM5594529_G3101_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "G3101", min.cells = 3, min.features = 200)
file62 <- RenameCells(file62, add.cell.id = "G3101")
Idents(file62) <- file62@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "G3101")
file62 <- subset(x = file62, idents = rownames(meta))
file62 <- AddMetaData(file62, meta)

file63 <- CreateSeuratObject(counts = Read10X_h5("GSM5594530_H3102_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "H3102", min.cells = 3, min.features = 200)
file63 <- RenameCells(file63, add.cell.id = "H3102")
Idents(file63) <- file63@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "H3102")
file63 <- subset(x = file63, idents = rownames(meta))
file63 <- AddMetaData(file63, meta)

file64 <- CreateSeuratObject(counts = Read10X_h5("GSM5594531_O3103_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "O3103", min.cells = 3, min.features = 200)
file64 <- RenameCells(file64, add.cell.id = "O3103")
Idents(file64) <- file64@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "O3103")
file64 <- subset(x = file64, idents = rownames(meta))
file64 <- AddMetaData(file64, meta)

file65 <- CreateSeuratObject(counts = Read10X_h5("GSM5594532_P3104_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "P3104", min.cells = 3, min.features = 200)
file65 <- RenameCells(file65, add.cell.id = "P3104")
Idents(file65) <- file65@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "P3104")
file65 <- subset(x = file65, idents = rownames(meta))
file65 <- AddMetaData(file65, meta)

file66 <- CreateSeuratObject(counts = Read10X_h5("GSM5594533_J3105_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "J3105", min.cells = 3, min.features = 200)
file66 <- RenameCells(file66, add.cell.id = "J3105")
Idents(file66) <- file66@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "J3105")
file66 <- subset(x = file66, idents = rownames(meta))
file66 <- AddMetaData(file66, meta)

file67 <- CreateSeuratObject(counts = Read10X_h5("GSM5594534_U3106_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "U3106", min.cells = 3, min.features = 200)
file67 <- RenameCells(file67, add.cell.id = "U3106")
Idents(file67) <- file67@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "U3106")
file67 <- subset(x = file67, idents = rownames(meta))
file67 <- AddMetaData(file67, meta)

file68 <- CreateSeuratObject(counts = Read10X_h5("GSM5594535_V3107_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "V3107", min.cells = 3, min.features = 200)
file68 <- RenameCells(file68, add.cell.id = "V3107")
Idents(file68) <- file68@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "V3107")
file68 <- subset(x = file68, idents = rownames(meta))
file68 <- AddMetaData(file68, meta)

file69 <- CreateSeuratObject(counts = Read10X_h5("GSM5594536_O3108_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "O3108", min.cells = 3, min.features = 200)
file69 <- RenameCells(file69, add.cell.id = "O3108")
Idents(file69) <- file69@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "O3108")
file69 <- subset(x = file69, idents = rownames(meta))
file69 <- AddMetaData(file69, meta)

file70 <- CreateSeuratObject(counts = Read10X_h5("GSM5594537_P3109_raw_feature_bc_matrix.h5", use.names = TRUE, unique.features = TRUE), project = "P3109", min.cells = 3, min.features = 200)
file70 <- RenameCells(file70, add.cell.id = "P3109")
Idents(file70) <- file70@assays$RNA@counts@Dimnames[2]
meta <- subset(meta_data, meta_data$orig.ident == "P3109")
file70 <- subset(x = file70, idents = rownames(meta))
file70 <- AddMetaData(file70, meta)

Atlas7 <- merge(file61, y = c(file62, file63, file64, file65,
                              file66, file67, file68, file69, file70))
saveRDS(Atlas7, file = 'Mouse_Atlas_7_7_subset.rds')

Atlas <- merge(Atlas1, y = c(Atlas2, Atlas3, Atlas4, Atlas5, Atlas6, Atlas7))
saveRDS(Atlas, file = 'Humphreys_1M_Atlas_final.rds')



#IRImouse
mtx_control <- read.table(file = '/home/kloetzer/data/HumphreyIRIMouse/GSE139107_MouseIRI_control.dge.txt')
mtx_1 <- read.table(file = '/home/kloetzer/data/HumphreyIRIMouse/GSE139107_MouseIRI_4hours.dge.txt')
mtx_2 <- read.table(file = '/home/kloetzer/data/HumphreyIRIMouse/GSE139107_MouseIRI_12hours.dge.txt')
mtx_3 <- read.table(file = '/home/kloetzer/data/HumphreyIRIMouse/GSE139107_MouseIRI_2days.dge.txt')
mtx_4 <- read.table(file = '/home/kloetzer/data/HumphreyIRIMouse/GSE139107_MouseIRI_14days.dge.txt')
mtx_5 <- read.table(file = '/home/kloetzer/data/HumphreyIRIMouse/GSE139107_MouseIRI_6weeks.dge.txt')
meta_data <- read.table(file = '/home/kloetzer/data/HumphreyIRIMouse/GSE139107_MouseIRI.metadata.txt')

object_control <- CreateSeuratObject(counts = mtx_control, project = "IRIHumphrey")
object_1 <- CreateSeuratObject(counts = mtx_1, project = "IRIHumphrey")
object_2 <- CreateSeuratObject(counts = mtx_2, project = "IRIHumphrey")
object_3 <- CreateSeuratObject(counts = mtx_3, project = "IRIHumphrey")
object_4 <- CreateSeuratObject(counts = mtx_4, project = "IRIHumphrey")
Object_5 <- CreateSeuratObject(counts = mtx_5, project = "IRIHumphrey")
IRI_mice <- merge(object_control, y = c(object_1, object_2, object_3, object_4, Object_5), project = "IRIHumphrey")
IRI_mice <- AddMetaData(
  object = IRI_mice,
  metadata = meta_data)
IRI_mice[["percent.mt"]] <- PercentageFeatureSet(IRI_mice, pattern = "^mt-")
object <- IRI_mice
saveRDS(file = '/home/kloetzer/data/HumphreyIRIMouse/IRI_mouse.rds')



sessionInfo()
