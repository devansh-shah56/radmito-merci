# reading the coverage and variants file
library(MERCI)
library(Seurat)
library(dplyr)
library(ggplot2)
library(reticulate)

# location of the file; coverage and variants file
varFile_new <- "/home/devanshs/MERCI/MERCI_processed_files/kidwell/SRR15339531/out/SRR15339531.MT_variants.txt"
coverFile_new <- "/home/devanshs/MERCI/MERCI_processed_files/kidwell/SRR15339531/out/SRR15339531.Coverage_Cell.csv"

# reading variants file
data_test_new <- read.table(varFile_new, header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = NULL) 
mtSNVTable_new <- readMTvar_10x(varFile_new, minReads = 1000)

# reading the annotated RDS object
obj3_new <- readRDS("/home/devanshs/MERCI/MERCI_processed_files/kidwell/SRR15339531/obj_3_annotated.rds")

# extracting the cell barcodes of macrophage and cancer cells
mac_cell_new <- colnames(obj3_new)[ which(!is.na(obj3_new$cell_type) & obj3_new$cell_type == "macrophage") ]
can_cell_new <- colnames(obj3_new)[ which(!is.na(obj3_new$cell_type) & obj3_new$cell_type == "mda-mb-231") ]
selected_Cells_new <- c(mac_cell_new, can_cell_new)

# reading coverage file
MTCoverage_info_new <- readCoverage_10x(coverFile_new, S.cells = selected_Cells_new)

# next only using the selected cells
s_mtSNVTable_new <- mtSNVTable_new[mtSNVTable_new$Cell %in% selected_Cells_new, ]

### MERCI - LOO; generate the VAF matrix
mtSNV_ma3_new <- MTmutMatrix_transform(
  MT_variants   = s_mtSNVTable_new,
  MT_coverage   = MTCoverage_info_new,
  donor_cells   = can_cell_new,
  mixed_cells   = selected_Cells_new,
  min_d         = 5,
  min_observeRate = 0.1
)

can_cell_sel_new <- intersect(can_cell_new, colnames(mtSNV_ma3_new))
mac_cell_new <- intersect(mac_cell_new, colnames(mtSNV_ma3_new))

# getting donor enriched mtSNVs (donor = cancer; mixed = cancer + macrophage)
s_muts3_new <- Denrich_mtMut_extr(
  varMatrix   = mtSNV_ma3_new,
  donor_cells = can_cell_new,
  mixed_cells = c(can_cell_sel_new, mac_cell_new),
  OR          = 2,
  Nmut_min    = 2
)

# DNA rank (donor = cancer; mixed = cancer + macrophage)
DNA_rank3_new <- Cell_Neff_cal(
  varMatrix   = mtSNV_ma3_new,
  MT_variants = s_mtSNVTable_new,
  MT_coverage = MTCoverage_info_new,
  donor_cells = can_cell_new,
  mixed_cells = c(can_cell_sel_new, mac_cell_new),
  mutFeatures = s_muts3_new,
  adjust      = FALSE
)

# getting the cell expression for selected cells
# create a subset Seurat object containing only the selected cells
subset_obj2_new <- subset(obj3_new, cells = c(mac_cell_new, can_cell_sel_new))

expr_mat_subset2_new <- Seurat::GetAssayData(subset_obj2_new, assay = "RNA", slot = "counts")

# RNA rank (donor = cancer; receiver = cancer + macrophage)
RNA_rank3_new <- MERCI_LOO_MT_est(
  expr_mat_subset2_new,
  reciever_cells = c(mac_cell_new, can_cell_sel_new),
  donor_cells    = can_cell_new
)

# finding the transferred cells
CellN_stat3_new <- CellNumber_test(DNA_rank3_new, RNA_rank3_new, Number_R = 1000)

# predicting cells
MTreceiver_pre3_new <- MERCI_ReceiverPre(DNA_rank3_new, RNA_rank3_new, top_rank = 50)

# Keep only rows for cells that actually exist in subset_obj2_new
present_cells_new <- intersect(MTreceiver_pre3_new$cell, colnames(subset_obj2_new))
MTreceiver_present_new <- MTreceiver_pre3_new[MTreceiver_pre3_new$cell %in% present_cells_new, , drop = FALSE]

# Create a named vector where names are barcodes and values are predictions
flag_vec_new <- setNames(MTreceiver_present_new$prediction, MTreceiver_present_new$cell)

# Add to Seurat metadata; this will create a column "receiver_flag"
subset_obj2_new <- Seurat::AddMetaData(subset_obj2_new, metadata = flag_vec_new, col.name = "receiver_flag")

# Normalize
subset_obj2_new <- NormalizeData(subset_obj2_new, normalization.method = "LogNormalize", scale.factor = 1e4)
subset_obj2_new <- FindVariableFeatures(subset_obj2_new, selection.method = "vst", nfeatures = 2000)
subset_obj2_new <- ScaleData(subset_obj2_new, features = VariableFeatures(subset_obj2_new))
subset_obj2_new <- RunPCA(subset_obj2_new, features = VariableFeatures(subset_obj2_new), verbose = FALSE)
subset_obj2_new <- RunUMAP(subset_obj2_new, dims = 1:25, verbose = FALSE)

# saving object as RDS (SRR15339531)
saveRDS(subset_obj2_new, "/home/devanshs/MERCI/MERCI_processed_files/kidwell/SRR15339531/SRR15339531-mt.rds")

p1_new <- DimPlot(subset_obj2_new, reduction = "umap", group.by = "receiver_flag", pt.size = 0.5) +
  ggplot2::ggtitle("UMAP colored by receiver_flag")

# UMAP with labels
p2_new <- DimPlot(subset_obj2_new, reduction = "umap", group.by = "cell_type", pt.size = 0.5) +
  ggplot2::ggtitle("UMAP colored by cell_type")

# Print plots
print(p1_new)
print(p2_new)

# 2. Single UMAP overlay: color by cell_type and use shape to indicate receiver_flag
p_overlay_shape_new <- DimPlot(
  subset_obj2_new,
  reduction = "umap",
  group.by = "cell_type",     # color = cell types
  shape.by = "receiver_flag", # shape = receiver vs non-receiver
  pt.size = 0.6,
  label = FALSE
) + ggplot2::ggtitle("UMAP: cell_type (color) with receiver_flag (shape)")

print(p_overlay_shape_new)
