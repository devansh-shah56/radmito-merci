# RADMito-MERCI

A computational pipeline for tracing mitochondrial transfer between cancer cells and immune cells using single-cell sequencing data.

## Project Overview

This repository implements **MERCI** (Mitochondrial-Enabled Reconstruction of Cellular Interactions), a statistical deconvolution method designed to trace and quantify mitochondrial trafficking between cancer cells and T cells at single-cell resolution. The project builds upon the original MERCI framework published by Zhang et al. (2023) in *Cancer Cell*.

### Research Context

Recent studies have revealed that cancer cells can hijack mitochondria from nearby T cells through tunneling nanotubes or extracellular vesicles. This process metabolically empowers cancer cells while depleting immune cells, representing a novel mechanism of immune evasion. Understanding this mitochondrial transfer dynamic is important for cancer biology and potential therapeutic interventions.

## Key Features

- **Mitochondrial variant calling** from single-cell RNA-seq data
- **Dual-ranking system** combining DNA mutations and gene expression profiles
- **Receiver cell identification** based on donor-derived mitochondrial signatures
- **Visualization tools** for UMAP plotting and heteroplasmy analysis
- **Support for both 10X Genomics and other single-cell platforms**

## Pipeline Architecture

The pipeline consists of two main modules:

### 1. MERCI-mtSNP (Python)
Calls mitochondrial single nucleotide variants from aligned sequencing data.

**Input Requirements:**
- BAM/BAI files (output from CellRanger, STAR, or similar alignment tools)
- Reference genome FASTA file
- Cell barcodes file (for 10X data)

**Supported Data Types:**
- 10X scRNA-seq
- 10X mtscATAC-seq
- Smart-seq2
- Bulk RNA-seq/ATAC-seq
- BD-Rhapsody scRNA-seq

**Output:**
- Per-cell mitochondrial coverage CSV file
- Mitochondrial variants TXT file

### 2. MERCI R Package
Predicts mitochondrial recipient cells and quantifies transferred mitochondria composition.

**Core Functions:**
- `readMTvar_10x()` - Read mitochondrial variant files
- `readCoverage_10x()` - Process coverage data
- `MTmutMatrix_transform()` - Generate variant-by-cell matrix
- `Denrich_mtMut_extr()` - Identify donor-enriched mutations
- `Cell_Neff_cal()` - Calculate DNA rank scores
- `MERCI_LOO_MT_est()` - Estimate RNA rank and mitochondrial fractions
- `MERCI_ReceiverPre()` - Predict receiver cells

## Installation

### Python Dependencies

```bash
pip install pandas numpy pysam matplotlib
```

### R Dependencies

```r
# Install MERCI package
devtools::install_github("shyhihihi/MERCI")

# Required R packages
install.packages(c("Seurat", "dplyr", "ggplot2", "Matrix"))
```

## Usage

### Data Preprocessing

#### 1. Generate Required Input Files

For 10X Genomics data, you need three files:
- `matrix.mtx` - Sparse count matrix
- `features.tsv` - Gene annotations (must include MT genes)
- `barcodes.tsv` - Cell barcodes

Create Seurat-compatible RDS objects:

```r
library(Seurat)

# Read 10X data
data_dir <- "/path/to/10x_dir"
counts <- Read10X(data.dir = data_dir)

# Create sparse matrix
cell_exp <- as(counts, "dgCMatrix")

# Create cell metadata
cell_names <- colnames(cell_exp)
cell_info <- data.frame(
  cell_name = cell_names,
  cell_type = NA_character_,
  MTtransfer = NA_character_,
  culture_history = NA_character_,
  stringsAsFactors = FALSE
)
rownames(cell_info) <- cell_info$cell_name

# Save objects
saveRDS(cell_exp, "cell_exp.RDS")
saveRDS(cell_info, "cell_info.RDS")
```

#### 2. Call Mitochondrial Variants

```bash
python MERCI-mtSNP.py \
  -D 10x_scRNA-seq \
  -o output_directory \
  -S sample_name \
  -b path/to/bamfile.bam \
  -f path/to/reference.fa \
  -c path/to/barcodes_dir \
  -M 20 \
  -s human
```

### MERCI Analysis Pipeline

```r
library(MERCI)
library(Seurat)

# 1. Read mitochondrial variants
varFile <- "path/to/mtSNV_profile.txt"
mtSNV <- readMTvar_10x(varFile, min_reads = 1000)

# 2. Read coverage data
coverFile <- "path/to/coverage.csv"
MT_coverage <- readCoverage_10x(coverFile, S.cells = selected_cells)

# 3. Generate variant matrix
mtSNV_matrix <- MTmutMatrix_transform(mtSNV, threshold = 0.01)

# 4. Identify donor-enriched mutations
donor_muts <- Denrich_mtMut_extr(
  varMatrix = mtSNV_matrix,
  donor_cells = donor_cell_names,
  receiver_cells = receiver_cell_names,
  OddsRatio = 5,
  pvalue = 0.01,
  qvalue = 0.05
)

# 5. Calculate DNA rank
DNA_rank <- Cell_Neff_cal(
  varMatrix = mtSNV_matrix,
  MT_variants = mtSNV,
  MT_coverage = MT_coverage,
  donor_cells = donor_cell_names,
  mixed_cells = receiver_cell_names,
  mutFeatures = donor_muts,
  adjust = FALSE
)

# 6. Load expression data
cell_exp <- readRDS("cell_exp.RDS")

# 7. Calculate RNA rank
RNA_rank <- MERCI_LOO_MT_est(
  cell_exp,
  receiver_cells = receiver_cell_names,
  donor_cells = donor_cell_names,
  organism = 'human'
)

# 8. Predict receiver cells
receiver_pred <- MERCI_ReceiverPre(
  DNA_rank,
  RNA_rank,
  cutoff = 0.5
)
```

### Visualization

```r
# Create Seurat object with MERCI results
seurat_obj <- CreateSeuratObject(counts = cell_exp)
seurat_obj$mito_receiver <- receiver_pred$prediction

# Standard Seurat workflow
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- RunPCA(seurat_obj)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)

# Plot receiver cells
DimPlot(seurat_obj, group.by = "mito_receiver")
ggsave("receiver_cells_umap.png", width = 8, height = 6, dpi = 300)
```

## Project Structure

```
radmito-merci/
├── scripts/
│   ├── 01_preprocessing.R
│   ├── 02_variant_calling.sh
│   └── 03_merci_analysis.R
├── data/
│   ├── raw/
│   └── processed/
├── output/
│   ├── figures/
│   └── results/
├── docs/
└── README.md
```

## Key Concepts

### DNA Rank Score
Quantifies the presence of donor-derived mitochondrial mutations in receiver cells. Calculated using an effective count statistic based on donor-enriched variants.

### RNA Rank Score
Estimates the relative abundance of transferred mitochondria using support vector regression on mitochondrial gene expression profiles.

### Receiver Cell Identification
Cells are classified as receivers when both DNA and RNA rank scores exceed predefined thresholds (typically top 10-50%, depending on expected receiver proportion).

## Biological Interpretation

### Mitochondrial Transfer Phenotype
MERCI identifies cells that have acquired mitochondria from donor cells. In cancer contexts, this often represents:

- **T cell to cancer cell transfer**: Metabolic empowerment of cancer cells
- **Cancer cell to T cell transfer**: T cell dysfunction and immune evasion

### Signature Genes
The mitochondrial transfer phenotype is associated with:
- Cytoskeleton remodeling genes (PIM1, MYO1B, PFN1, ABI1)
- Energy production pathways
- TNF-alpha signaling
- Cell adhesion molecules (PVRL2)

## Technical Considerations

### Data Quality Requirements
- Minimum mitochondrial read coverage: 1000 reads per cell
- Recommended sequencing depth: 500K-1M reads per cell for scRNA-seq
- Must retain mitochondrial genes in feature set (check with `grep '^MT-'`)

### Normalization
- Use median nCount_RNA as scaling factor for normalization
- Avoid adding +1 to sparse count matrices
- Apply log transformation after scaling

### Parameter Selection
- **DNA/RNA rank cutoff**: 10% for stringent, 50% for sensitive detection
- **Odds Ratio threshold**: 5-10 for donor enrichment
- **Coverage threshold**: 10X minimum per variant site

## Troubleshooting

### Common Issues

1. **Missing cell reads in variant file**
   - Ensure mtSNV file includes all required columns
   - Check barcode compatibility between files

2. **Low receiver detection**
   - Verify receiver proportion is above 7% of sample
   - Adjust cutoff thresholds
   - Check data quality and coverage

3. **Memory issues with ScaleData**
   - Use SCTransform as alternative to NormalizeData
   - Process data in batches
   - Increase available RAM

## Project Goals

This project aims to:

1. Establish a reproducible pipeline for MERCI analysis
2. Apply the method to diverse cancer datasets
3. Compare mitochondrial transfer patterns across:
   - Different cancer types
   - Metastatic vs non-metastatic tumors
   - Treatment-naive vs treated samples
4. Integrate with heteroplasmy analysis
5. Validate findings through experimental datasets

## Course Information

This project was developed as part of the DH607 (Introduction to Computational Multi-omics) course at IIT Bombay, focusing on single-cell genomics and mitochondrial biology.

## References

1. Zhang H, et al. (2023). Systematic investigation of mitochondrial transfer between cancer cells and T cells at single-cell resolution. *Cancer Cell*, 41(10):1788-1802.

2. MERCI Software: https://github.com/shyhihihi/MERCI

3. Immune evasion through mitochondrial transfer paper (2024). *Nature*, 636:447-455.

## Contributors

- Devansh Shah
- Ravi Kant

## License

MIT License

## Contact

For questions or collaboration inquiries, please open an issue on this repository.

---

**Note**: This pipeline is under active development. The methods and parameters may be updated as the project progresses.