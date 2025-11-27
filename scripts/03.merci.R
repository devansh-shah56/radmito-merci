#reading the coverage and variants file
library(MERCI)
library(Seurat)
library(dplyr)
library(ggplot2)
library(reticulate)

# location of the file; coverage and variants file
varFile <- "/home/devanshs/MERCI/MERCI_processed_files/kidwell/SRR15339529/out/SRR15339529.MT_variants.txt"
coverFile <- "/home/devanshs/MERCI/MERCI_processed_files/kidwell/SRR15339529/out/SRR15339529.Coverage_Cell.csv"

# reading variants file
data_test <- read.table(varFile, header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = NULL) 
mtSNVTable <- readMTvar_10x(varFile, minReads = 1000)


#identifying the cancer cells and macrophage cells, LoadH5Seurat: Load a saved 'Seurat' object from an h5Seurat file in R
#data <- read_h5ad("/home/devanshs/MERCI/MERCI_processed_files/kidwell/exp-info/GSE181410_mda-mb-231.h5ad", as ="HDF5AnnData")
#data <- CreateSeuratObject(counts = t(as.matrix(data$X)), meta.data = data$obs,min.features = 500, min.cells = 30)
#snhx = schard::h5ad2seurat('GSE181410_mda-mb-231.h5ad')
#saveRDS(snhx,"GSE181410_mda-mb-231.rds")

# reading the annotated RDS object
obj1 <- readRDS("/home/devanshs/MERCI/MERCI_processed_files/kidwell/SRR15339529/obj_1_annotated.rds")

# extracting the cell barcodes of macrophage and cancer cells
mac_cell <- colnames(obj1)[ which(!is.na(obj1$cell_type) & obj1$cell_type == "macrophage") ]
can_cell <- colnames(obj1)[ which(!is.na(obj1$cell_type) & obj1$cell_type == "mda-mb-231") ]

selected_Cells <- c(mac_cell, can_cell)

# reading coverage file
MTCoverage_info <- readCoverage_10x(coverFile, S.cells=cat(mac_cell, can_cell))

#next only ysing the selected cells
s.mtSNVTable <- mtSNVTable[mtSNVTable$Cell%in%selected_Cells, ]

###MERCI - LOO; generate the VAF matrix
mtSNV_ma3 <- MTmutMatrix_transform(MT_variants=s.mtSNVTable, MT_coverage=MTCoverage_info, donor_cells=mac_cell, mixed_cells=selected_Cells, min_d=5, min_observeRate= 0.1)

can_cell_sel <- intersect(can_cell, colnames(mtSNV_ma3))

# getting donor enriched mtSNVs
s.muts3 <- Denrich_mtMut_extr(varMatrix=mtSNV_ma3, donor_cells=mac_cell, mixed_cells=c(mac_cell, can_cell_sel), OR=2, Nmut_min=2)

# DNA rank
DNA_rank3 <- Cell_Neff_cal(varMatrix=mtSNV_ma3, MT_variants=s.mtSNVTable, MT_coverage=MTCoverage_info, donor_cells=mac_cell, mixed_cells=c(can_cell_sel, mac_cell), mutFeatures=s.muts3, adjust=FALSE) #calculating the Neff (DNA rank)

# getting the cell expression for selected cells
# create a subset Seurat object containing only the selected cells
subset_obj2 <- subset(obj1, cells = c(mac_cell, can_cell_sel))

expr_mat_subset2 <- Seurat::GetAssayData(subset_obj2, assay = "RNA", slot = "counts")

# RNA rank
RNA_rank3 <- MERCI_LOO_MT_est(expr_mat_subset2, reciever_cells=c(mac_cell, can_cell_sel), donor_cells=mac_cell)

# finding the transferred cells
CellN_stat3 <- CellNumber_test(DNA_rank3, RNA_rank3, Number_R=1000)

# predicting cells
MTreceiver_pre3 <- MERCI_ReceiverPre(DNA_rank3, RNA_rank3, top_rank=50)

#adding to the seurat object
# Keep only rows for cells that actually exist in subset_obj2
present_cells <- intersect(MTreceiver_pre3$cell, colnames(subset_obj2))
MTreceiver_present <- MTreceiver_pre3[MTreceiver_pre3$cell %in% present_cells, , drop = FALSE]

# Create a named vector where names are barcodes and values are predictions
flag_vec <- setNames(MTreceiver_present$prediction, MTreceiver_present$cell)

# Add to Seurat metadata; this will create a column "receiver_flag"
subset_obj2 <- Seurat::AddMetaData(subset_obj2, metadata = flag_vec, col.name = "receiver_flag")

# Normalize
subset_obj2 <- NormalizeData(subset_obj2, normalization.method = "LogNormalize", scale.factor = 1e4)
subset_obj2 <- FindVariableFeatures(subset_obj2, selection.method = "vst", nfeatures = 2000)
subset_obj2 <- ScaleData(subset_obj2, features = VariableFeatures(subset_obj2))
subset_obj2 <- RunPCA(subset_obj2, features = VariableFeatures(subset_obj2), verbose = FALSE)
subset_obj2 <- RunUMAP(subset_obj2, dims = 1:25, verbose = FALSE)

#saving object as RDS
saveRDS(subset_obj2, "/home/devanshs/MERCI/MERCI_processed_files/kidwell/SRR15339529/SRR15339529-mt.rds")

# Basic UMAP colored by the receiver_flag metadata
p1 <- DimPlot(subset_obj2, reduction = "umap", group.by = "receiver_flag", pt.size = 0.5) +
  ggplot2::ggtitle("UMAP colored by receiver_flag")

# UMAP with labels
p2 <- DimPlot(subset_obj2, reduction = "umap", group.by = "cell_type", pt.size = 0.5) +
  ggplot2::ggtitle("UMAP colored by cell_type")

# Print plots
print(p1)
print(p2)

p_side_by_side <- DimPlot(
  subset_obj2,
  reduction = "umap",
  group.by = "cell_type",     # color by cell type
  split.by = "receiver_flag", # create separate panels for receiver states
  pt.size = 0.6,
  label = FALSE
) + ggplot2::ggtitle("UMAP: cell_type split by receiver_flag")

# Print
print(p_side_by_side)

# 2. Single UMAP overlay: color by cell_type and use shape to indicate receiver_flag
#    (shape.by works when there are few unique values in receiver_flag)
p_overlay_shape <- DimPlot(
  subset_obj2,
  reduction = "umap",
  group.by = "cell_type",    # color = cell types
  shape.by = "receiver_flag",# shape = receiver vs non-receiver
  pt.size = 0.6,
  label = FALSE
) + ggplot2::ggtitle("UMAP: cell_type (color) with receiver_flag (shape)")

print(p_overlay_shape)
