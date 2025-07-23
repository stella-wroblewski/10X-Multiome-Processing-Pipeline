# Multiome and Xenium Processing Pipeline

**Author**: Stella Wroblewski
**Date**: 07/23/2025  

## Overview

This repository contains an R pipeline for processing 10x Genomics Multiome (scRNA + scATAC) and Xenium spatial transcriptomics data. The pipeline performs quality control, normalization, dimensionality reduction, clustering, and multimodal integration using Seurat and Signac.

## Key Features

1. **Dual Processing**: Handles both Multiome (RNA + ATAC) and Xenium spatial data
2. **Memory Management**: Configurable settings for processing large datasets
3. **Quality Control**: Modality-specific QC thresholds for RNA, ATAC, and spatial data
4. **Integration Methods**: Weighted Nearest Neighbors (WNN) for multimodal data integration
5. **Export Formats**: Outputs h5ad files for Python compatibility
6. **Organism Support**: Works with mouse and human data

## Repository Contents

- **`multiome_xenium_pipeline.R`**
  - Main pipeline script with:
    - User configuration section
    - Automated package installation
    - Processing functions for each modality
    - Integration workflows
    - Visualization and export utilities

- **`README.md`**
  - This documentation file

## Getting Started

### 1. System Requirements

- R version ≥ 4.0.0
- 16GB RAM minimum (32GB recommended)
- 10GB free disk space per sample
- 4+ CPU cores recommended

### 2. Dependencies

The pipeline will automatically install required packages:

**Core Packages:**
- Seurat (v4+)
- Signac
- SeuratDisk
- Matrix, future

**Visualization:**
- ggplot2, patchwork, viridis, RColorBrewer

**Bioconductor:**
- EnsDb.Mmusculus.v79 (mouse) or EnsDb.Hsapiens.v86 (human)
- BSgenome packages

### 3. Input Data Requirements

#### Multiome Data
- `.h5` file: 10x filtered feature matrix
- `.tsv.gz` file: ATAC fragments (must be indexed with .tbi)
- `.bed` file: Peak calls (optional)

#### Xenium Data  
- `cell_feature_matrix.h5`: Expression matrix
- `cells.csv.gz`: Cell metadata with spatial coordinates

### 4. Configuration

Edit the configuration section at the top of the script:

```r
# Set input file paths
INPUT_CONFIG <- list(
  multiome_rna_h5 = "path/to/sample.h5",
  multiome_atac_fragments = "path/to/fragments.tsv.gz",
  xenium_matrix_h5 = "path/to/cell_feature_matrix.h5",
  xenium_cells_csv = "path/to/cells.csv.gz",
  sample_name = "my_sample",
  organism = "mouse"  # or "human"
)

# Adjust memory settings
MEMORY_CONFIG <- list(
  max_cores = 4,
  max_memory_gb = 32,
  save_intermediates = TRUE,
  downsample_plots = TRUE,
  max_cells_plot = 20000
)

# Set QC thresholds
QC_CONFIG <- list(
  rna = list(
    min_features = 200,
    max_features = 10000,
    max_mt_percent = 20
  ),
  atac = list(
    min_fragments = 1000,
    max_fragments = 100000,
    tss_enrichment_min = 2
  ),
  xenium = list(
    min_features = 100,
    max_features = 5000
  )
)
```

### 5. Running the Pipeline

```bash
# Run the pipeline
Rscript multiome_xenium_pipeline.R

# Monitor progress
tail -f processed_data/preprocessing_log.txt
```

## Pipeline Workflow

### 1. RNA Processing (Multiome)
- Load 10x h5 file
- Calculate QC metrics (% mitochondrial, % ribosomal)
- Filter cells based on QC thresholds
- Normalize data (log normalization)
- Find variable features (2000 genes)
- Scale data, run PCA
- Run UMAP and find clusters

### 2. ATAC Processing (Multiome)
- Load fragment files
- Create genomic bins
- Calculate ATAC QC metrics (TSS enrichment, nucleosome signal)
- Filter cells
- Run TF-IDF normalization
- Perform LSI dimensionality reduction
- Run UMAP and find clusters

### 3. Multimodal Integration
- Find cells present in both RNA and ATAC
- Run WNN integration
- Generate integrated UMAP
- Find multimodal clusters
- Calculate modality weights

### 4. Xenium Processing
- Load expression matrix and spatial coordinates
- Calculate QC metrics
- Filter cells
- Standard scRNA-seq workflow (normalize, scale, PCA, UMAP)
- Generate spatial visualizations

### 5. Output Generation
- Save Seurat objects (.rds files)
- Export to h5ad format
- Generate QC plots and UMAPs
- Save cluster markers
- Create processing summary

## Output Structure

```
processed_data/
├── figures/              # QC plots, UMAPs, spatial plots
├── qc/                   # Metrics CSVs, marker genes
├── objects/              # Seurat objects (if save_intermediates=TRUE)
├── h5ad/                 # Python-compatible files
└── preprocessing_log.txt # Detailed log with timestamps
```

## Key Functions

**`log_message()`**: Logs messages with timestamps to console and file

**`save_qc_metrics()`**: Exports cell metadata and QC metrics to CSV

**`create_qc_plots()`**: Generates violin plots for QC metrics

**`save_plot()`**: Saves plots with consistent settings

**`convert_to_h5ad()`**: Converts Seurat objects to h5ad format

## Troubleshooting

### Memory Issues
- Reduce `max_cores` in MEMORY_CONFIG
- Set `save_intermediates = TRUE` to save progress
- Process samples individually

### File Not Found
- Use absolute paths in INPUT_CONFIG
- Ensure fragment files have .tbi index
- Check file permissions

### Package Installation
- Run `BiocManager::install()` for Bioconductor packages
- Check R version compatibility
- Install system dependencies (hdf5, etc.)

## Example Usage

```r
# Process a mouse Multiome + Xenium dataset
INPUT_CONFIG <- list(
  multiome_rna_h5 = "/data/multiome/filtered_feature_bc_matrix.h5",
  multiome_atac_fragments = "/data/multiome/atac_fragments.tsv.gz",
  xenium_matrix_h5 = "/data/xenium/cell_feature_matrix.h5",
  xenium_cells_csv = "/data/xenium/cells.csv.gz",
  sample_name = "mouse_brain_001",
  organism = "mouse"
)

# Run pipeline
source("multiome_xenium_pipeline.R")
```

## References

- [Seurat](https://satijalab.org/seurat/)
- [Signac](https://stuartlab.org/signac/)
- [10x Genomics Multiome](https://www.10xgenomics.com/products/single-cell-multiome-atac-plus-gene-expression)
- [10x Genomics Xenium](https://www.10xgenomics.com/platforms/xenium)

