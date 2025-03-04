#' Title: Multiome (10x) Data Processing Pipeline for Combined scRNA-seq and scATAC-seq
#' Author: Stella Wroblewski
#' Date: 03/04/2025
#' Description:
#' This script provides a generalized framework for processing 10x multiome (scRNA + scATAC) data.
#' It showcases an end-to-end analysis pipeline using Seurat and Signac. Steps include data loading,
#' quality control, data integration (using SCTransform, Harmony, and Weighted Nearest Neighbor),
#' cell type annotation, differential expression, and various downstream visualizations.

######################################################################
# 1. Load Necessary Libraries
######################################################################
library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(EnsDb.Mmusculus.v79)
library(hdf5r)
library(patchwork)
library(harmony)

######################################################################
# 2. Define a Function to Process Individual Samples
######################################################################

process_sample <- function(file_path, frag_file, sample_id) {
  # Description: Reads 10x multiome data (RNA + ATAC), performs basic QC,
  # and returns a Seurat object with RNA and ATAC assays.

  # 1. Read the 10x hdf5 file
  inputdata.10x <- Read10X_h5(file_path)

  # 2. Extract RNA and ATAC data
  rna_counts <- inputdata.10x$`Gene Expression`
  atac_counts <- inputdata.10x$Peaks

  # 3. Create Seurat object from RNA data
  sample_obj <- CreateSeuratObject(counts = rna_counts)

  # 4. Calculate mitochondrial percentage (modify pattern for different organisms)
  sample_obj[["percent.mt"]] <- PercentageFeatureSet(sample_obj, pattern = "^MT-")

  # 5. Prepare ATAC data
  grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
  grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
  atac_counts <- atac_counts[as.vector(grange.use), ]

  # 6. Obtain gene annotations for the relevant genome
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
  seqlevelsStyle(annotations) <- 'UCSC'
  genome(annotations) <- "mm10"  # Adjust if using a different organism

  # 7. Create ChromatinAssay and add to the Seurat object
  chrom_assay <- CreateChromatinAssay(
    counts = atac_counts,
    sep = c(":", "-"),
    genome = 'mm10',  # Adjust genome as needed
    fragments = frag_file,
    min.cells = 10,
    annotation = annotations
  )
  sample_obj[["ATAC"]] <- chrom_assay

  # 8. Basic QC filtering (adjust thresholds as appropriate)
  sample_obj <- subset(
    x = sample_obj,
    subset = nCount_ATAC < 7e4 &
      nCount_ATAC > 5e3 &
      nCount_RNA < 25000 &
      nCount_RNA > 1000 &
      percent.mt < 20
  )

  # 9. RNA analysis (SCT, PCA, UMAP)
  DefaultAssay(sample_obj) <- "RNA"
  sample_obj <- SCTransform(sample_obj, verbose = FALSE) %>%
    RunPCA() %>%
    RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

  # 10. ATAC analysis (TF-IDF, SVD, UMAP)
  DefaultAssay(sample_obj) <- "ATAC"
  sample_obj <- RunTFIDF(sample_obj)
  sample_obj <- FindTopFeatures(sample_obj, min.cutoff = 'q0')
  sample_obj <- RunSVD(sample_obj)
  sample_obj <- RunUMAP(
    sample_obj,
    reduction = 'lsi',
    dims = 2:50,
    reduction.name = "umap.atac",
    reduction.key = "atacUMAP_"
  )

  # 11. Add sample ID
  sample_obj$orig.ident <- sample_id
  return(sample_obj)
}

######################################################################
# 3. Example: Processing Multiple Samples
######################################################################
# Below is a skeleton example of how you might process multiple samples.
# Replace the placeholders with your own file paths and sample IDs.

# file_paths <- list(
#   sample1 = "path/to/sample1_filtered_feature_bc_matrix.h5",
#   sample2 = "path/to/sample2_filtered_feature_bc_matrix.h5"
# )
# frag_files <- list(
#   sample1 = "path/to/sample1_atac_fragments.tsv.gz",
#   sample2 = "path/to/sample2_atac_fragments.tsv.gz"
# )

# sample_ids <- names(file_paths)
# raw_samples <- mapply(process_sample, file_paths, frag_files, sample_ids, SIMPLIFY = FALSE)
# names(raw_samples) <- sample_ids

######################################################################
# 4. Integration Workflow
######################################################################
# After processing each sample, you can integrate them.
# Below is an example using SCT integration + Harmony, followed by WNN.

# 4A. Ensure Variable Features Identified on Each Sample
# for (i in seq_along(raw_samples)) {
#   DefaultAssay(raw_samples[[i]]) <- "RNA"
#   raw_samples[[i]] <- SCTransform(raw_samples[[i]], verbose = FALSE, variable.features.n = 2000)
# }

# 4B. Rename Cells to Keep Them Unique
# for (i in seq_along(raw_samples)) {
#   raw_samples[[i]] <- RenameCells(raw_samples[[i]], add.cell.id = sample_ids[i])
# }

# 4C. Select Integration Features (SCT)
# DefaultAssay(raw_samples[[1]]) <- "SCT"
# features <- SelectIntegrationFeatures(object.list = raw_samples, nfeatures = 3000)

# 4D. Run PCA on Each Sample Prior to Merging
# for (i in seq_along(raw_samples)) {
#   DefaultAssay(raw_samples[[i]]) <- "SCT"
#   raw_samples[[i]] <- RunPCA(raw_samples[[i]], features = features, verbose = FALSE)
# }

# 4E. Merge into a Single Object
# combined <- merge(raw_samples[[1]], y = raw_samples[-1], add.cell.id = sample_ids, merge.data = TRUE)

# 4F. Run PCA Again on Combined Object
# DefaultAssay(combined) <- "SCT"
# combined <- RunPCA(combined, features = features, verbose = FALSE)

# 4G. Run Harmony Integration
# combined <- RunHarmony(
#   object = combined,
#   group.by.vars = "orig.ident",
#   dims = 1:30
# )

# 4H. UMAP on Harmony Embeddings
# combined <- RunUMAP(
#   combined,
#   reduction = "harmony",
#   dims = 1:30,
#   reduction.name = "umap.rna",
#   reduction.key = "rnaUMAP_",
#   verbose = FALSE
# )

# 4I. ATAC Analysis on Combined Object
# DefaultAssay(combined) <- "ATAC"
# combined <- RunTFIDF(combined)
# combined <- FindTopFeatures(combined, min.cutoff = 'q0')
# combined <- RunSVD(combined)
# combined <- RunUMAP(
#   combined,
#   dims = 2:50,
#   reduction = "lsi",
#   reduction.name = "umap.atac",
#   reduction.key = "atacUMAP_",
#   verbose = FALSE
# )

# 4J. Weighted Nearest Neighbors (WNN)
# combined <- FindMultiModalNeighbors(
#   combined,
#   reduction.list = list("harmony", "lsi"),
#   dims.list = list(1:30, 2:50)
# )
# combined <- RunUMAP(
#   combined,
#   nn.name = "weighted.nn",
#   reduction.name = "wnn.umap",
#   reduction.key = "wnnUMAP_"
# )

######################################################################
# 5. Example: Cell Type Annotation
######################################################################
# You can use your preferred cell type annotation method (sc-type, SingleR, etc.).
# For example:

# library(openxlsx)
# library(HGNChelper)
# source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# gs_list <- gene_sets_prepare("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx", "Brain")
# DefaultAssay(combined) <- "SCT"
# scRNAseqData <- as.matrix(combined[["SCT"]]@scale.data)
# es.max <- sctype_score(
#   scRNAseqData = scRNAseqData,
#   scaled = TRUE,
#   gs = gs_list$gs_positive,
#   gs2 = gs_list$gs_negative
# )
# combined@meta.data$sctype_classification <- NA
# for (i in colnames(es.max)) {
#   best_cell_type <- rownames(es.max)[which.max(es.max[, i])]
#   combined@meta.data$sctype_classification[match(i, rownames(combined@meta.data))] <- best_cell_type
# }

######################################################################
# 6. Differential Expression & Accessibility
######################################################################
# Suppose you split samples into groups. Then:

# combined$group <- ifelse(combined$orig.ident %in% c("sample1", "sample2"), "GroupA", "GroupB")

# RNA DE (using SCT assay)
# DefaultAssay(combined) <- "SCT"
# combined <- PrepSCTFindMarkers(combined)
# de_genes <- FindMarkers(
#   combined,
#   ident.1 = "GroupA",
#   ident.2 = "GroupB",
#   assay = "SCT"
# )

# ATAC DE (Peak Accessibility)
# DefaultAssay(combined) <- "ATAC"
# de_accessibility <- FindMarkers(
#   combined,
#   ident.1 = "GroupA",
#   ident.2 = "GroupB",
#   assay = "ATAC"
# )

######################################################################
# 7. Visualization Examples
######################################################################
# 7A. Dimensional Reductions
# p1 <- DimPlot(combined, reduction = "umap.rna", group.by = "orig.ident")
# p2 <- DimPlot(combined, reduction = "umap.atac", group.by = "orig.ident")
# p3 <- DimPlot(combined, reduction = "wnn.umap", group.by = "orig.ident")
# p1 + p2 + p3

# 7B. FeaturePlot on SCT or a Custom RNA_activity Assay
# DefaultAssay(combined) <- "SCT"
# FeaturePlot(combined, features = c("Gene1", "Gene2"), reduction = "wnn.umap")

######################################################################
# End of Pipeline
######################################################################
# This script is intended as a template for multiome data analysis.
# Adjust parameters, thresholds, and references for your specific project.

