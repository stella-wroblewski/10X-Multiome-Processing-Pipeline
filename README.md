# Multiome Processing Pipeline (scRNA + scATAC)

**Author**: Stella Wroblewski
**Date**: 03/04/2025

## Overview

This repository provides a generalized pipeline for processing 10x multiome data (simultaneous single-cell RNA sequencing and single-cell ATAC sequencing). The code in `multiome_pipeline_publish_ready.R` demonstrates how to load, filter, normalize, integrate, and analyze multiome datasets, enabling downstream analyses such as dimensional reduction, clustering, annotation, differential expression (RNA), and differential accessibility (ATAC).

## Key Features
1. **Modular Design**: The pipeline is organized into sections, making it easy to modify or reuse specific functions.
2. **Seurat + Signac**: Uses Seurat (for scRNA) and Signac (for scATAC) to handle multi-modal single-cell data.
3. **Customizable QC Steps**: Allows user-defined thresholds for basic quality control on RNA and ATAC data.
4. **Integration Methods**: Demonstrates how to integrate data with Harmony and Weighted Nearest Neighbors (WNN).
5. **Cell Type Annotation**: Provides templates for sc-type, SingleR, or other annotation approaches.
6. **Differential Expression/Accessibility**: Shows example workflows for discovering DE genes/peaks.
7. **Flexible**: You can adjust the genome, annotation database, and thresholds for your specific organism or experiment.

## Repository Contents

- **multiome_pipeline_publish_ready.R**
  - A script containing:
    - **`process_sample()`** function: A custom function to read, filter, and process scRNA + scATAC from a single 10x multiome sample.
    - Example usage for multiple samples: how to load data from multiple samples, integrate them, and run downstream analyses.
    - Sections on QC, dimensional reduction, integration (Harmony + WNN), differential expression, and visualization.

- **README.md**
  - The document you are reading now, describing how to set up and run the pipeline, including key dependencies.

## Getting Started

### 1. Dependencies and Installation

This pipeline heavily relies on the following R packages:

- **[Seurat](https://satijalab.org/seurat/)** (v4+ recommended)
- **[Signac](https://stuartlab.org/signac/)** (for scATAC)
- **[Harmony](https://github.com/immunogenomics/harmony)** (for data integration)
- **[dplyr](https://dplyr.tidyverse.org/)** (data manipulation)
- **[patchwork](https://patchwork.data-imaginist.com/)** (plot arrangement)
- **[EnsDb.Mmusculus.v79](https://bioconductor.org/packages/release/data/annotation/html/EnsDb.Mmusculus.v79.html)** (mouse annotation database; modify if working with another organism)
- **[hdf5r](https://cran.r-project.org/web/packages/hdf5r/index.html)** (for reading 10x .h5 files)

Additionally, you may find these packages useful:
- **[ggplot2](https://ggplot2.tidyverse.org/)** (data visualization)
- **[clusterProfiler](https://yulab-smu.top/clusterProfiler-book/)**, **[ReactomePA](https://www.bioconductor.org/packages/release/bioc/html/ReactomePA.html)**, **[org.Mm.eg.db](https://bioconductor.org/packages/org.Mm.eg.db/)** (for enrichment analyses)
- **[openxlsx](https://cran.r-project.org/web/packages/openxlsx/index.html)**, **[HGNChelper](https://cran.r-project.org/web/packages/HGNChelper/index.html)**, etc.

Make sure to install all packages before running the pipeline:
```r
install.packages("devtools")
install.packages("remotes")
# Example for installing Signac
remotes::install_github("timoast/signac")
# Install Harmony
remotes::install_github("immunogenomics/harmony")
# For Bioconductor packages
if(!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
BiocManager::install("EnsDb.Mmusculus.v79")
BiocManager::install("org.Mm.eg.db")
BiocManager::install("ReactomePA")
# etc.
```

### 2. Data Requirements

1. **10x HDF5 File** (`.h5`): Typically named `filtered_feature_bc_matrix.h5`.
2. **ATAC Fragments File** (`.tsv.gz`): Contains chromatin accessibility info.
3. **Annotation**: Provided via an EnsDb database or other annotation resource suitable for your organism.

### 3. Usage Instructions

1. **Clone or download** the repository.
2. **Open `multiome_pipeline_publish_ready.R`** in R or RStudio.
3. **Set file paths**:
   - Replace the placeholder `file_paths` and `frag_files` with paths to your own data.
   - Keep a consistent naming scheme in `sample_ids` to facilitate merging.
4. **Run the `process_sample()`** function for each sample, producing a list of Seurat objects.
5. **Integration** (optional but recommended if you have multiple samples or batches):
   - Use SCT, Harmony, or WNN integration as described.
6. **Cell Type Annotation**: Insert your method of choice (sc-type, SingleR, or manual marker-based annotation).
7. **Differential Expression / Accessibility**: Adjust groupings and run `FindMarkers()` to identify features that differ between conditions.
8. **Visualization**: Generate UMAPs, FeaturePlots, DotPlots, heatmaps, or your own specialized visualizations.

### 4. Pipeline Outline

1. **Data Ingestion**
   - Read HDF5 with `Read10X_h5()`
   - Create RNA and ATAC assays
2. **Quality Control**
   - Filter cells by `nCount_RNA`, `nCount_ATAC`, and `% mitochondrial reads`
3. **Normalization & Dimensional Reduction**
   - RNA: `SCTransform`, `RunPCA`, `RunUMAP`
   - ATAC: `RunTFIDF`, `RunSVD`, `RunUMAP`
4. **Merging & Integration** (Optional)
   - Merge objects, run Harmony (batch correction), Weighted Nearest Neighbor (multi-modal integration)
5. **Downstream Analysis**
   - Find clusters, cell-type annotation, differential expression/accessibility
6. **Visualization**
   - Plot UMAP projections, gene expression patterns, etc.

### 5. Tips & Best Practices

1. **Parameter Tuning**: QC thresholds are data dependent; experiment with them.
2. **Annotation Sources**: Always double-check your annotation (e.g., for human vs. mouse) and convert gene IDs where needed.
3. **Storage & Memory**: Multiome data can be large. Consider HPC or cloud resources if running into memory constraints.
4. **Version Control**: Keep track of package versions using a lock file or environment manager (e.g., `renv`).

### 6. Example Code Snippet

Below is a minimal example to show how you might set up your environment:

```r
# Step 1: Load libraries
library(Seurat)
library(Signac)
library(EnsDb.Mmusculus.v79)

# Step 2: Define your file paths
file_paths <- list(
  sample1 = "~/path/to/sample1_filtered_feature_bc_matrix.h5",
  sample2 = "~/path/to/sample2_filtered_feature_bc_matrix.h5"
)
frag_files <- list(
  sample1 = "~/path/to/sample1_atac_fragments.tsv.gz",
  sample2 = "~/path/to/sample2_atac_fragments.tsv.gz"
)
sample_ids <- names(file_paths)

# Step 3: Process each sample
raw_samples <- mapply(
  FUN = process_sample,
  file_path = file_paths,
  frag_file = frag_files,
  sample_id = sample_ids,
  SIMPLIFY = FALSE
)

# Step 4 (optional): Integrate multiple samples
# ... SCT, Harmony, or WNN steps ...

# Step 5: Analyze, visualize, etc.

```

### 7. Known Issues or Caveats

- **Large Datasets**: Processing very large datasets may require more memory or HPC resources.
- **Organism-Specific**: The script uses mouse (`EnsDb.Mmusculus.v79`) by default. For other organisms, change the database or remove the annotation block.
- **Paths & Filenames**: The function `process_sample()` expects 10x `.h5` files with `Gene Expression` and `Peaks` keys, plus the corresponding `.tsv.gz` fragments.

### 8. References

- [Seurat v4](https://satijalab.org/seurat/)
- [Signac](https://stuartlab.org/signac/)
- [Harmony Integration](https://github.com/immunogenomics/harmony)
- [sc-type](https://github.com/IanevskiAleksandr/sc-type)
- [SingleR](https://bioconductor.org/packages/SingleR)

### 9. Contact

For questions, bug reports, or improvements, please open an issue or contact Stella Wroblewski at swroblewski@tulane.edu.

---

**Enjoy your single-cell multiome analysis!**

