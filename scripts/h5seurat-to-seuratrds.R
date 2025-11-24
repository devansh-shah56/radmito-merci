# Author - Devansh Shah
# To convert the gene expression information (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181410) and extract barcodes to label cells as macrophages and cancer in the MERCI-LOO pipeline

# donwload necessary packages; SeuratDisk and schard
devtools::install_github("cellgeni/schard")
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")

#import libraries
library(SeuratDisk)
library(schard)

#converting h5seurat to h5ad
Convert("GSE181410_macrophage.h5seurat", "GSE181410_macrophage.h5ad")
Convert("GSE181410_mda-mb-231.h5seurat", "GSE181410_mda-mb-231.h5ad")

#converting h5ad to seurat RDS
snhx = schard::h5ad2seurat('GSE181410_mda-mb-231.h5ad')
saveRDS(snhx,"GSE181410_mda-mb-231.rds")

snhx1 = schard::h5ad2seurat('GSE181410_macrophage.h5ad')
saveRDS(snhx1,"GSE181410_macrophage.rds")
