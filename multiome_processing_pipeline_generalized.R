#!/usr/bin/env Rscript

#' ============================================================================
#' MULTIOME AND XENIUM PREPROCESSING PIPELINE
#' ============================================================================
#' 
#' Description:
#'   A comprehensive pipeline for preprocessing 10x Multiome (RNA + ATAC) and 
#'   10x Xenium spatial transcriptomics data. This pipeline performs quality 
#'   control, normalization, dimensionality reduction, clustering, and 
#'   integration analysis.
#' 
#' Features:
#'   - Memory-optimized processing for large datasets
#'   - Publication-ready visualizations
#'   - Multimodal integration using Weighted Nearest Neighbor (WNN) analysis
#'   - Export to h5ad format for downstream analysis
#'   - Comprehensive QC metrics and reports
#' 
#' Author: [Stella Wroblewski]

#' 
#' ============================================================================

# ============================================================================
# USER CONFIGURATION - MODIFY THIS SECTION FOR YOUR DATA
# ============================================================================

# Set your input file paths here
INPUT_CONFIG <- list(
  # Multiome input files
  multiome_rna_h5 = "path/to/your/multiome_sample.h5",           # 10x h5 file with RNA counts
  multiome_atac_fragments = "path/to/your/fragments.tsv.gz",     # ATAC fragments file
  multiome_atac_peaks = "path/to/your/peaks.bed",                # Peak file (optional)
  
  # Xenium input files  
  xenium_matrix_h5 = "path/to/your/cell_feature_matrix.h5",      # Xenium expression matrix
  xenium_cells_csv = "path/to/your/cells.csv.gz",                # Xenium cell metadata
  
  # Sample metadata
  sample_name = "sample1",                                        # Name for your sample
  organism = "mouse"                                              # Options: "mouse" or "human"
)

# Output directory configuration
OUTPUT_DIR <- "processed_data"

# Memory configuration (adjust based on your system)
MEMORY_CONFIG <- list(
  max_cores = 4,              # Number of CPU cores to use
  max_memory_gb = 32,         # Maximum memory in GB
  save_intermediates = TRUE,  # Save intermediate objects to disk
  downsample_plots = TRUE,    # Downsample large datasets for plotting
  max_cells_plot = 20000      # Maximum cells to use in plots
)

# Quality control thresholds
QC_CONFIG <- list(
  # RNA QC thresholds
  rna = list(
    min_features = 200,       # Minimum genes per cell
    max_features = 10000,     # Maximum genes per cell
    max_mt_percent = 20,      # Maximum mitochondrial percentage
    max_ribo_percent = 30,    # Maximum ribosomal percentage
    min_cells = 3             # Minimum cells expressing a gene
  ),
  
  # ATAC QC thresholds
  atac = list(
    min_fragments = 1000,           # Minimum fragments per cell
    max_fragments = 100000,         # Maximum fragments per cell
    min_frip = 0.3,                 # Minimum fraction of reads in peaks
    max_blacklist_ratio = 0.05,     # Maximum blacklist region ratio
    nucleosome_signal_max = 10,     # Maximum nucleosome signal
    tss_enrichment_min = 2          # Minimum TSS enrichment
  ),
  
  # Xenium QC thresholds
  xenium = list(
    min_features = 100,       # Minimum genes per cell
    max_features = 5000,      # Maximum genes per cell
    max_mt_percent = 20,      # Maximum mitochondrial percentage
    min_cells = 3             # Minimum cells expressing a gene
  )
)

# ============================================================================
# SETUP AND INITIALIZATION
# ============================================================================

# Create output directory structure
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUTPUT_DIR, "figures"), showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "qc"), showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "objects"), showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "h5ad"), showWarnings = FALSE)

# Initialize logging
log_file <- file.path(OUTPUT_DIR, "preprocessing_log.txt")
log_conn <- file(log_file, open = "w")

#' Log messages to both console and file
#' @param message Character string to log
log_message <- function(message) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  full_message <- paste0("[", timestamp, "] ", message)
  cat(full_message, "\n")
  cat(full_message, "\n", file = log_conn, append = TRUE)
}

log_message("Starting preprocessing pipeline")
log_message(paste("Output directory:", OUTPUT_DIR))
log_message(paste("Sample name:", INPUT_CONFIG$sample_name))

# ============================================================================
# PACKAGE INSTALLATION AND LOADING
# ============================================================================

# Define required packages
required_packages <- c(
  # Core packages
  "Seurat", "Signac", "SeuratDisk", "Matrix", "future",
  
  # Visualization
  "ggplot2", "patchwork", "viridis", "RColorBrewer", "cowplot",
  
  # Data manipulation
  "dplyr", "tidyr", "stringr",
  
  # File I/O
  "rhdf5", "hdf5r", "rtracklayer",
  
  # Bioconductor packages - will be specified separately
  "BiocManager"
)

# Bioconductor packages based on organism
bioc_packages <- if (INPUT_CONFIG$organism == "mouse") {
  c("EnsDb.Mmusculus.v79", "BSgenome.Mmusculus.UCSC.mm10")
} else {
  c("EnsDb.Hsapiens.v86", "BSgenome.Hsapiens.UCSC.hg38")
}

# Install missing packages
log_message("Checking and installing required packages...")
missing_packages <- required_packages[!required_packages %in% installed.packages()[,"Package"]]

if (length(missing_packages) > 0) {
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

# Install Bioconductor packages if needed
if (!all(bioc_packages %in% installed.packages()[,"Package"])) {
  if (!"BiocManager" %in% installed.packages()[,"Package"]) {
    install.packages("BiocManager")
  }
  BiocManager::install(bioc_packages, update = FALSE, ask = FALSE)
}

# Load all packages
log_message("Loading required packages...")
suppressPackageStartupMessages({
  # Core analysis
  library(Seurat)
  library(Signac)
  library(SeuratDisk)
  library(Matrix)
  library(future)
  
  # Visualization
  library(ggplot2)
  library(patchwork)
  library(viridis)
  library(RColorBrewer)
  library(cowplot)
  
  # Data manipulation
  library(dplyr)
  library(tidyr)
  library(stringr)
  
  # File I/O
  library(rhdf5)
  library(hdf5r)
  library(rtracklayer)
  
  # Load organism-specific packages
  if (INPUT_CONFIG$organism == "mouse") {
    library(EnsDb.Mmusculus.v79)
    library(BSgenome.Mmusculus.UCSC.mm10)
  } else {
    library(EnsDb.Hsapiens.v86)
    library(BSgenome.Hsapiens.UCSC.hg38)
  }
})

# ============================================================================
# SYSTEM CONFIGURATION
# ============================================================================

# Configure parallel processing
n_cores <- min(MEMORY_CONFIG$max_cores, parallel::detectCores() - 1)
plan("multicore", workers = n_cores)
options(future.globals.maxSize = MEMORY_CONFIG$max_memory_gb * 1024^3)

# Set random seed for reproducibility
set.seed(42)

# Configure Seurat options
options(Seurat.object.assay.version = "v5")

log_message(paste("Using", n_cores, "cores for processing"))
log_message(paste("Maximum memory set to", MEMORY_CONFIG$max_memory_gb, "GB"))

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

#' Save QC metrics to CSV file
#' @param object Seurat object
#' @param stage Character string describing the processing stage
save_qc_metrics <- function(object, stage) {
  qc_data <- object[[]]
  qc_file <- file.path(OUTPUT_DIR, "qc", paste0(INPUT_CONFIG$sample_name, "_", stage, "_metrics.csv"))
  write.csv(qc_data, qc_file, row.names = TRUE)
  log_message(paste("Saved QC metrics for", stage))
  return(qc_data)
}

#' Create standardized violin plots for QC metrics
#' @param object Seurat object
#' @param features Vector of features to plot
#' @param title Plot title
create_qc_plots <- function(object, features, title) {
  # Downsample if needed
  if (MEMORY_CONFIG$downsample_plots && ncol(object) > MEMORY_CONFIG$max_cells_plot) {
    cells_use <- sample(colnames(object), MEMORY_CONFIG$max_cells_plot)
    object <- subset(object, cells = cells_use)
  }
  
  p <- VlnPlot(object, features = features, pt.size = 0, ncol = 2) & 
    theme_minimal() &
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  p <- p + plot_annotation(title = title)
  return(p)
}

#' Save plot with consistent settings
#' @param plot ggplot object
#' @param filename Output filename
#' @param width Plot width in inches
#' @param height Plot height in inches
save_plot <- function(plot, filename, width = 10, height = 8) {
  output_path <- file.path(OUTPUT_DIR, "figures", filename)
  ggsave(output_path, plot, width = width, height = height, dpi = 300, bg = "white")
  log_message(paste("Saved plot:", filename))
}

#' Clean up memory after large operations
clean_memory <- function() {
  gc()
  invisible()
}

# ============================================================================
# PROCESS MULTIOME RNA DATA
# ============================================================================

log_message("Processing Multiome RNA data")

# Check if input file exists
if (!file.exists(INPUT_CONFIG$multiome_rna_h5)) {
  stop("Multiome RNA file not found: ", INPUT_CONFIG$multiome_rna_h5)
}

# Load RNA data
log_message("Loading RNA count matrix...")
rna_counts <- Read10X_h5(INPUT_CONFIG$multiome_rna_h5)

# Handle both single and multi-genome cases
if (is.list(rna_counts)) {
  # Multi-genome data - take the first genome
  rna_counts <- rna_counts[[1]]
}

# Create Seurat object
rna_obj <- CreateSeuratObject(
  counts = rna_counts,
  project = INPUT_CONFIG$sample_name,
  min.cells = QC_CONFIG$rna$min_cells,
  min.features = 0  # Filter after QC metrics calculation
)

# Add sample metadata
rna_obj$orig.ident <- INPUT_CONFIG$sample_name
rna_obj$technology <- "Multiome"

# Calculate QC metrics
log_message("Calculating RNA QC metrics...")
# Mitochondrial genes - handle both mouse and human
if (INPUT_CONFIG$organism == "mouse") {
  rna_obj[["percent.mt"]] <- PercentageFeatureSet(rna_obj, pattern = "^mt-|^Mt-")
  rna_obj[["percent.ribo"]] <- PercentageFeatureSet(rna_obj, pattern = "^Rp[sl]")
} else {
  rna_obj[["percent.mt"]] <- PercentageFeatureSet(rna_obj, pattern = "^MT-")
  rna_obj[["percent.ribo"]] <- PercentageFeatureSet(rna_obj, pattern = "^RP[SL]")
}

# Save pre-QC metrics
pre_qc_metrics <- save_qc_metrics(rna_obj, "rna_pre_qc")

# Create QC plots
qc_plot <- create_qc_plots(
  rna_obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
  title = "RNA Pre-QC Metrics"
)
save_plot(qc_plot, "rna_pre_qc.png", width = 12, height = 10)

# Apply QC filters
log_message("Filtering cells based on QC thresholds...")
rna_obj <- subset(
  rna_obj,
  subset = nFeature_RNA > QC_CONFIG$rna$min_features &
           nFeature_RNA < QC_CONFIG$rna$max_features &
           percent.mt < QC_CONFIG$rna$max_mt_percent &
           percent.ribo < QC_CONFIG$rna$max_ribo_percent
)

# Save post-QC metrics
post_qc_metrics <- save_qc_metrics(rna_obj, "rna_post_qc")

# Create post-QC plots
qc_plot_post <- create_qc_plots(
  rna_obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
  title = "RNA Post-QC Metrics"
)
save_plot(qc_plot_post, "rna_post_qc.png", width = 12, height = 10)

# Standard RNA processing workflow
log_message("Normalizing RNA data...")
rna_obj <- NormalizeData(rna_obj)

log_message("Finding variable features...")
rna_obj <- FindVariableFeatures(rna_obj, selection.method = "vst", nfeatures = 2000)

# Plot variable features
var_plot <- VariableFeaturePlot(rna_obj)
var_plot <- LabelPoints(plot = var_plot, points = head(VariableFeatures(rna_obj), 10), repel = TRUE)
save_plot(var_plot, "rna_variable_features.png")

log_message("Scaling data...")
rna_obj <- ScaleData(rna_obj)

log_message("Running PCA...")
rna_obj <- RunPCA(rna_obj, npcs = 50, verbose = FALSE)

# Save elbow plot
elbow_plot <- ElbowPlot(rna_obj, ndims = 50) + 
  theme_minimal() + 
  ggtitle("RNA PCA Elbow Plot")
save_plot(elbow_plot, "rna_elbow_plot.png", width = 8, height = 6)

# Run UMAP and clustering
log_message("Running UMAP and clustering...")
rna_obj <- RunUMAP(rna_obj, dims = 1:30)
rna_obj <- FindNeighbors(rna_obj, dims = 1:30)
rna_obj <- FindClusters(rna_obj, resolution = 0.8)

# Create UMAP plot
umap_plot <- DimPlot(rna_obj, reduction = "umap", label = TRUE, repel = TRUE) +
  theme_minimal() +
  ggtitle("RNA UMAP - Clusters")
save_plot(umap_plot, "rna_umap_clusters.png")

# Find cluster markers
log_message("Finding RNA cluster markers...")
rna_markers <- FindAllMarkers(
  rna_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

# Save markers
write.csv(
  rna_markers,
  file.path(OUTPUT_DIR, "qc", paste0(INPUT_CONFIG$sample_name, "_rna_markers.csv")),
  row.names = FALSE
)

# Save RNA object
if (MEMORY_CONFIG$save_intermediates) {
  saveRDS(rna_obj, file.path(OUTPUT_DIR, "objects", paste0(INPUT_CONFIG$sample_name, "_rna.rds")))
  log_message("Saved RNA object")
}

# Log summary statistics
log_message(sprintf(
  "RNA processing complete: %d cells, %d genes, %d clusters",
  ncol(rna_obj),
  nrow(rna_obj),
  length(unique(rna_obj$seurat_clusters))
))

clean_memory()

# ============================================================================
# PROCESS MULTIOME ATAC DATA
# ============================================================================

log_message("Processing Multiome ATAC data")

# Check if input files exist
if (!file.exists(INPUT_CONFIG$multiome_atac_fragments)) {
  stop("ATAC fragments file not found: ", INPUT_CONFIG$multiome_atac_fragments)
}

# Get genome information
if (INPUT_CONFIG$organism == "mouse") {
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
  seqlevelsStyle(annotations) <- "UCSC"
  genome_name <- "mm10"
} else {
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
  seqlevelsStyle(annotations) <- "UCSC"
  genome_name <- "hg38"
}

# Create fragment object using cells from RNA
log_message("Creating ATAC fragment object...")
fragments <- CreateFragmentObject(
  path = INPUT_CONFIG$multiome_atac_fragments,
  cells = colnames(rna_obj)
)

# Create bins for initial quantification
log_message("Creating genomic bins...")
genome_bins <- GenomeBinMatrix(
  fragments = fragments,
  genome = seqinfo(annotations),
  binsize = 5000  # 5kb bins
)

# Create ATAC assay
atac_assay <- CreateChromatinAssay(
  counts = genome_bins,
  fragments = fragments,
  genome = genome_name,
  annotation = annotations
)

# Create ATAC Seurat object
atac_obj <- CreateSeuratObject(
  counts = atac_assay,
  assay = "ATAC",
  project = INPUT_CONFIG$sample_name
)

# Add metadata
atac_obj$orig.ident <- INPUT_CONFIG$sample_name
atac_obj$technology <- "Multiome"

# Calculate ATAC QC metrics
log_message("Calculating ATAC QC metrics...")
atac_obj <- NucleosomeSignal(atac_obj)
atac_obj <- TSSEnrichment(atac_obj, fast = FALSE)

# Additional QC metrics
total_fragments <- CountFragments(INPUT_CONFIG$multiome_atac_fragments)
rownames(total_fragments) <- total_fragments$CB
atac_obj$fragments <- total_fragments[colnames(atac_obj), "frequency_count"]
atac_obj$pct_reads_in_peaks <- atac_obj$peak_region_fragments / atac_obj$fragments * 100

# Save pre-QC metrics
pre_qc_metrics_atac <- save_qc_metrics(atac_obj, "atac_pre_qc")

# Create ATAC QC plots
atac_qc_features <- c("fragments", "TSS.enrichment", "nucleosome_signal", "pct_reads_in_peaks")
atac_qc_features <- atac_qc_features[atac_qc_features %in% colnames(atac_obj[[]])]

if (length(atac_qc_features) > 0) {
  qc_plot_atac <- create_qc_plots(
    atac_obj,
    features = atac_qc_features,
    title = "ATAC Pre-QC Metrics"
  )
  save_plot(qc_plot_atac, "atac_pre_qc.png", width = 12, height = 10)
}

# Apply ATAC QC filters
log_message("Filtering cells based on ATAC QC thresholds...")
atac_obj <- subset(
  atac_obj,
  subset = fragments > QC_CONFIG$atac$min_fragments &
           fragments < QC_CONFIG$atac$max_fragments &
           TSS.enrichment > QC_CONFIG$atac$tss_enrichment_min &
           nucleosome_signal < QC_CONFIG$atac$nucleosome_signal_max
)

# Save post-QC metrics
post_qc_metrics_atac <- save_qc_metrics(atac_obj, "atac_post_qc")

# Standard ATAC processing workflow
log_message("Processing ATAC data...")
atac_obj <- RunTFIDF(atac_obj)
atac_obj <- FindTopFeatures(atac_obj, min.cutoff = 'q0')
atac_obj <- RunSVD(atac_obj)

# Plot correlation with sequencing depth
depth_cor <- DepthCor(atac_obj)
depth_plot <- data.frame(
  Component = 1:length(depth_cor),
  Correlation = depth_cor
) %>%
  ggplot(aes(x = Component, y = Correlation)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  ggtitle("ATAC LSI Component Depth Correlation")
save_plot(depth_plot, "atac_depth_correlation.png", width = 8, height = 6)

# Run UMAP and clustering (excluding first component)
log_message("Running UMAP and clustering for ATAC...")
atac_obj <- RunUMAP(atac_obj, reduction = 'lsi', dims = 2:30)
atac_obj <- FindNeighbors(atac_obj, reduction = 'lsi', dims = 2:30)
atac_obj <- FindClusters(atac_obj, resolution = 0.8)

# Create ATAC UMAP plot
umap_plot_atac <- DimPlot(atac_obj, reduction = "umap", label = TRUE, repel = TRUE) +
  theme_minimal() +
  ggtitle("ATAC UMAP - Clusters")
save_plot(umap_plot_atac, "atac_umap_clusters.png")

# Save ATAC object
if (MEMORY_CONFIG$save_intermediates) {
  saveRDS(atac_obj, file.path(OUTPUT_DIR, "objects", paste0(INPUT_CONFIG$sample_name, "_atac.rds")))
  log_message("Saved ATAC object")
}

# Log summary statistics
log_message(sprintf(
  "ATAC processing complete: %d cells, %d features, %d clusters",
  ncol(atac_obj),
  nrow(atac_obj),
  length(unique(atac_obj$seurat_clusters))
))

clean_memory()

# ============================================================================
# MULTIMODAL INTEGRATION (WNN)
# ============================================================================

log_message("Performing multimodal integration with WNN")

# Find common cells between RNA and ATAC
common_cells <- intersect(colnames(rna_obj), colnames(atac_obj))
log_message(paste("Found", length(common_cells), "common cells between RNA and ATAC"))

# Subset to common cells
rna_obj <- subset(rna_obj, cells = common_cells)
atac_obj <- subset(atac_obj, cells = common_cells)

# Create multimodal object
multimodal_obj <- rna_obj
multimodal_obj[["ATAC"]] <- atac_obj[["ATAC"]]

# Run WNN analysis
log_message("Running WNN analysis...")
multimodal_obj <- FindMultiModalNeighbors(
  multimodal_obj,
  reduction.list = list("pca", "lsi"),
  dims.list = list(1:30, 2:30),
  modality.weight.name = "RNA.weight"
)

# Run UMAP on WNN
multimodal_obj <- RunUMAP(
  multimodal_obj,
  nn.name = "weighted.nn",
  reduction.name = "wnn.umap",
  reduction.key = "wnnUMAP_"
)

# Find clusters on WNN
multimodal_obj <- FindClusters(
  multimodal_obj,
  graph.name = "wsnn",
  algorithm = 3,
  resolution = 0.8
)

# Create WNN visualization
wnn_plot <- DimPlot(
  multimodal_obj,
  reduction = "wnn.umap",
  label = TRUE,
  repel = TRUE
) +
  theme_minimal() +
  ggtitle("WNN UMAP - Integrated Clusters")
save_plot(wnn_plot, "wnn_umap_clusters.png")

# Plot modality weights
weight_plot <- VlnPlot(
  multimodal_obj,
  features = "RNA.weight",
  group.by = "seurat_clusters",
  pt.size = 0
) +
  theme_minimal() +
  ggtitle("RNA Weight by Cluster")
save_plot(weight_plot, "wnn_modality_weights.png", width = 10, height = 6)

# Save multimodal object
if (MEMORY_CONFIG$save_intermediates) {
  saveRDS(multimodal_obj, file.path(OUTPUT_DIR, "objects", paste0(INPUT_CONFIG$sample_name, "_multimodal.rds")))
  log_message("Saved multimodal object")
}

# Find WNN markers
log_message("Finding WNN cluster markers...")
DefaultAssay(multimodal_obj) <- "RNA"
wnn_markers <- FindAllMarkers(
  multimodal_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

write.csv(
  wnn_markers,
  file.path(OUTPUT_DIR, "qc", paste0(INPUT_CONFIG$sample_name, "_wnn_markers.csv")),
  row.names = FALSE
)

clean_memory()

# ============================================================================
# PROCESS XENIUM DATA
# ============================================================================

log_message("Processing Xenium spatial data")

# Check if input files exist
if (!file.exists(INPUT_CONFIG$xenium_matrix_h5)) {
  stop("Xenium matrix file not found: ", INPUT_CONFIG$xenium_matrix_h5)
}
if (!file.exists(INPUT_CONFIG$xenium_cells_csv)) {
  stop("Xenium cells file not found: ", INPUT_CONFIG$xenium_cells_csv)
}

# Load Xenium data
log_message("Loading Xenium data...")
xenium_counts <- Read10X_h5(INPUT_CONFIG$xenium_matrix_h5)

# Handle both single and multi-genome cases
if (is.list(xenium_counts)) {
  xenium_counts <- xenium_counts[[1]]
}

# Load cell metadata
cell_metadata <- read.csv(INPUT_CONFIG$xenium_cells_csv)
rownames(cell_metadata) <- cell_metadata$cell_id

# Create Seurat object
xenium_obj <- CreateSeuratObject(
  counts = xenium_counts,
  project = paste0(INPUT_CONFIG$sample_name, "_xenium"),
  meta.data = cell_metadata,
  min.cells = QC_CONFIG$xenium$min_cells
)

# Add metadata
xenium_obj$orig.ident <- INPUT_CONFIG$sample_name
xenium_obj$technology <- "Xenium"

# Calculate QC metrics
log_message("Calculating Xenium QC metrics...")
if (INPUT_CONFIG$organism == "mouse") {
  xenium_obj[["percent.mt"]] <- PercentageFeatureSet(xenium_obj, pattern = "^mt-|^Mt-")
} else {
  xenium_obj[["percent.mt"]] <- PercentageFeatureSet(xenium_obj, pattern = "^MT-")
}

# Save pre-QC metrics
pre_qc_xenium <- save_qc_metrics(xenium_obj, "xenium_pre_qc")

# Create QC plots
qc_plot_xenium <- create_qc_plots(
  xenium_obj,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  title = "Xenium Pre-QC Metrics"
)
save_plot(qc_plot_xenium, "xenium_pre_qc.png", width = 12, height = 8)

# Create spatial QC plot if coordinates exist
if (all(c("x_centroid", "y_centroid") %in% colnames(xenium_obj[[]]))) {
  spatial_qc <- ggplot(
    xenium_obj[[]],
    aes(x = x_centroid, y = y_centroid, color = nFeature_RNA)
  ) +
    geom_point(size = 0.5, alpha = 0.8) +
    scale_color_viridis() +
    theme_minimal() +
    coord_fixed() +
    ggtitle("Spatial Distribution of Gene Count")
  save_plot(spatial_qc, "xenium_spatial_qc.png", width = 10, height = 10)
}

# Apply QC filters
log_message("Filtering cells based on QC thresholds...")
xenium_obj <- subset(
  xenium_obj,
  subset = nFeature_RNA > QC_CONFIG$xenium$min_features &
           nFeature_RNA < QC_CONFIG$xenium$max_features &
           percent.mt < QC_CONFIG$xenium$max_mt_percent
)

# Save post-QC metrics
post_qc_xenium <- save_qc_metrics(xenium_obj, "xenium_post_qc")

# Standard processing workflow
log_message("Processing Xenium data...")
xenium_obj <- NormalizeData(xenium_obj)
xenium_obj <- FindVariableFeatures(xenium_obj, selection.method = "vst", nfeatures = 2000)
xenium_obj <- ScaleData(xenium_obj)
xenium_obj <- RunPCA(xenium_obj, npcs = 50, verbose = FALSE)

# Run UMAP and clustering
log_message("Running UMAP and clustering for Xenium...")
xenium_obj <- RunUMAP(xenium_obj, dims = 1:30)
xenium_obj <- FindNeighbors(xenium_obj, dims = 1:30)
xenium_obj <- FindClusters(xenium_obj, resolution = 0.8)

# Create UMAP plot
umap_xenium <- DimPlot(xenium_obj, reduction = "umap", label = TRUE, repel = TRUE) +
  theme_minimal() +
  ggtitle("Xenium UMAP - Clusters")
save_plot(umap_xenium, "xenium_umap_clusters.png")

# Create spatial plot if coordinates exist
if (all(c("x_centroid", "y_centroid") %in% colnames(xenium_obj[[]]))) {
  spatial_clusters <- ggplot(
    xenium_obj[[]],
    aes(x = x_centroid, y = y_centroid, color = factor(seurat_clusters))
  ) +
    geom_point(size = 0.5, alpha = 0.8) +
    theme_minimal() +
    coord_fixed() +
    ggtitle("Spatial Distribution of Clusters") +
    theme(legend.position = "right")
  save_plot(spatial_clusters, "xenium_spatial_clusters.png", width = 12, height = 10)
}

# Find markers
log_message("Finding Xenium cluster markers...")
xenium_markers <- FindAllMarkers(
  xenium_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

write.csv(
  xenium_markers,
  file.path(OUTPUT_DIR, "qc", paste0(INPUT_CONFIG$sample_name, "_xenium_markers.csv")),
  row.names = FALSE
)

# Save Xenium object
if (MEMORY_CONFIG$save_intermediates) {
  saveRDS(xenium_obj, file.path(OUTPUT_DIR, "objects", paste0(INPUT_CONFIG$sample_name, "_xenium.rds")))
  log_message("Saved Xenium object")
}

# Log summary statistics
log_message(sprintf(
  "Xenium processing complete: %d cells, %d genes, %d clusters",
  ncol(xenium_obj),
  nrow(xenium_obj),
  length(unique(xenium_obj$seurat_clusters))
))

clean_memory()

# ============================================================================
# EXPORT TO H5AD FORMAT
# ============================================================================

log_message("Converting objects to h5ad format")

#' Convert Seurat object to h5ad format
#' @param seurat_obj Seurat object to convert
#' @param output_name Name for the output file (without extension)
convert_to_h5ad <- function(seurat_obj, output_name) {
  tryCatch({
    # Save as h5Seurat first
    h5seurat_file <- file.path(OUTPUT_DIR, "h5ad", paste0(output_name, ".h5Seurat"))
    SaveH5Seurat(seurat_obj, filename = h5seurat_file, overwrite = TRUE)
    
    # Convert to h5ad
    h5ad_file <- file.path(OUTPUT_DIR, "h5ad", paste0(output_name, ".h5ad"))
    Convert(h5seurat_file, dest = h5ad_file, overwrite = TRUE)
    
    # Remove intermediate file
    file.remove(h5seurat_file)
    
    log_message(paste("Successfully converted", output_name, "to h5ad format"))
  }, error = function(e) {
    log_message(paste("Error converting", output_name, "to h5ad:", e$message))
  })
}

# Convert all objects
convert_to_h5ad(rna_obj, paste0(INPUT_CONFIG$sample_name, "_rna"))
convert_to_h5ad(atac_obj, paste0(INPUT_CONFIG$sample_name, "_atac"))
convert_to_h5ad(multimodal_obj, paste0(INPUT_CONFIG$sample_name, "_multimodal"))
convert_to_h5ad(xenium_obj, paste0(INPUT_CONFIG$sample_name, "_xenium"))

# ============================================================================
# GENERATE SUMMARY REPORT
# ============================================================================

log_message("Generating summary report")

# Create summary data frame
summary_data <- data.frame(
  Dataset = c("Multiome RNA", "Multiome ATAC", "Multimodal WNN", "Xenium"),
  Cells_Input = c(
    nrow(pre_qc_metrics),
    nrow(pre_qc_metrics_atac),
    length(common_cells),
    nrow(pre_qc_xenium)
  ),
  Cells_Output = c(
    ncol(rna_obj),
    ncol(atac_obj),
    ncol(multimodal_obj),
    ncol(xenium_obj)
  ),
  Features = c(
    nrow(rna_obj),
    nrow(atac_obj),
    nrow(multimodal_obj),
    nrow(xenium_obj)
  ),
  Clusters = c(
    length(unique(rna_obj$seurat_clusters)),
    length(unique(atac_obj$seurat_clusters)),
    length(unique(multimodal_obj$seurat_clusters)),
    length(unique(xenium_obj$seurat_clusters))
  ),
  stringsAsFactors = FALSE
)

# Calculate filtering rate
summary_data$Filter_Rate <- round(
  (summary_data$Cells_Input - summary_data$Cells_Output) / summary_data$Cells_Input * 100, 
  2
)

# Save summary
write.csv(
  summary_data,
  file.path(OUTPUT_DIR, "processing_summary.csv"),
  row.names = FALSE
)

# Print summary
print(summary_data)

# Create summary plot
summary_plot <- summary_data %>%
  tidyr::pivot_longer(cols = c("Cells_Input", "Cells_Output"), 
                      names_to = "Stage", 
                      values_to = "Cells") %>%
  ggplot(aes(x = Dataset, y = Cells, fill = Stage)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle("Cell Count Summary") +
  scale_fill_manual(values = c("Cells_Input" = "#3498db", "Cells_Output" = "#2ecc71"))

save_plot(summary_plot, "processing_summary.png", width = 10, height = 6)

# ============================================================================
# CLEANUP AND COMPLETION
# ============================================================================

log_message("Pipeline completed successfully!")
log_message(paste("All results saved to:", OUTPUT_DIR))
log_message("Summary of outputs:")
log_message("  - Figures: ./figures/")
log_message("  - QC metrics: ./qc/")
log_message("  - Seurat objects: ./objects/")
log_message("  - H5AD files: ./h5ad/")

# Close log file
close(log_conn)

# Print final message
cat("\n")
cat("========================================\n")
cat("PREPROCESSING COMPLETE!\n")
cat(paste("Results saved to:", OUTPUT_DIR, "\n"))
cat("========================================\n")
