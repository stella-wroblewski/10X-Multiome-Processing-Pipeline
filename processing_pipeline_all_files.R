# Load necessary libraries
library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(EnsDb.Mmusculus.v79)
library(hdf5r)
library(patchwork)
library(harmony)

# Function to process individual samples
process_sample <- function(file_path, frag_file, sample_id) {
  # Read the 10x hdf5 file
  inputdata.10x <- Read10X_h5(file_path)
  
  # Extract RNA and ATAC data
  rna_counts <- inputdata.10x$`Gene Expression`
  atac_counts <- inputdata.10x$Peaks
  
  # Create Seurat object for RNA data
  pbmc <- CreateSeuratObject(counts = rna_counts)
  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
  
  # Process ATAC data
  grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
  grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
  atac_counts <- atac_counts[as.vector(grange.use), ]
  annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
  seqlevelsStyle(annotations) <- 'UCSC'
  genome(annotations) <- "mm10"
  
  # Create ChromatinAssay
  chrom_assay <- CreateChromatinAssay(
    counts = atac_counts,
    sep = c(":", "-"),
    genome = 'mm10',
    fragments = frag_file,
    min.cells = 10,
    annotation = annotations
  )
  pbmc[["ATAC"]] <- chrom_assay
  
  #Basic QC
  pbmc <- subset(
    x = pbmc,
    subset = nCount_ATAC < 7e4 &
    nCount_ATAC > 5e3 &
    nCount_RNA < 25000 &
    nCount_RNA > 1000 &
    percent.mt < 20
  )
  
  # RNA analysis
  DefaultAssay(pbmc) <- "RNA"
  pbmc <- SCTransform(pbmc, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')

  # ATAC analysis
   DefaultAssay(pbmc) <- "ATAC"
   pbmc <- RunTFIDF(pbmc)
   pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
   pbmc <- RunSVD(pbmc)
   pbmc <- RunUMAP(pbmc, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
   pbmc$orig.ident <- sample_id
   return(pbmc)
}

# File paths and fragment file paths for each sample
file_paths <- list(
  i1 = "/Users/stellawroblewski/Library/CloudStorage/Box-Box/Stella Wroblewski/Projects/Multiome Processing Pipeline/10x sample files/i1_filtered_feature_bc_matrix.h5",
  i5 = "/Users/stellawroblewski/Library/CloudStorage/Box-Box/Stella Wroblewski/Projects/Multiome Processing Pipeline/10x sample files/I5_filtered_feature_bc_matrix.h5",
  m1 = "/Users/stellawroblewski/Library/CloudStorage/Box-Box/Stella Wroblewski/Projects/Multiome Processing Pipeline/10x sample files/m1_filtered_feature_bc_matrix.h5",
  m4 = "/Users/stellawroblewski/Library/CloudStorage/Box-Box/Stella Wroblewski/Projects/Multiome Processing Pipeline/10x sample files/M4_filtered_feature_bc_matrix.h5"
)
frag_files <- list(
  i1 = "/Users/stellawroblewski/Library/CloudStorage/Box-Box/Stella Wroblewski/Projects/Multiome Processing Pipeline/10x sample files/i1_atac_fragments.tsv.gz",
  i5 = "/Users/stellawroblewski/Library/CloudStorage/Box-Box/Stella Wroblewski/Projects/Multiome Processing Pipeline/10x sample files/i5_atac_fragments.tsv.gz",
  m1 = "/Users/stellawroblewski/Library/CloudStorage/Box-Box/Stella Wroblewski/Projects/Multiome Processing Pipeline/10x sample files/m1_atac_fragments.tsv.gz",
  m4 = "/Users/stellawroblewski/Library/CloudStorage/Box-Box/Stella Wroblewski/Projects/Multiome Processing Pipeline/10x sample files/m4_atac_fragments.tsv.gz"
)

# Step 1: Process each sample using mapply to pass both file paths and fragment files
sample_ids <- names(file_paths)
first_raw_samples <- mapply(process_sample, file_paths, frag_files, sample_ids, SIMPLIFY = FALSE)
names(first_raw_samples) <- sample_ids

raw_samples <- first_raw_samples

# Step 2: Ensure variable features are found for each sample in the RNA assay using SCTransform
for (i in seq_along(raw_samples)) {
  DefaultAssay(raw_samples[[i]]) <- "RNA"
  raw_samples[[i]] <- SCTransform(raw_samples[[i]], verbose = FALSE, variable.features.n = 2000)  # SCTransform will identify and save variable features
}

# Step 3: Add unique cell identifiers to ensure unique cell names
for (i in seq_along(raw_samples)) {
  raw_samples[[i]] <- RenameCells(raw_samples[[i]], add.cell.id = sample_ids[i])
}

# Step 4: Select integration features based on SCT assay
DefaultAssay(raw_samples[[1]]) <- "SCT"
features <- SelectIntegrationFeatures(object.list = raw_samples, nfeatures = 3000)

# Step 5: Prepare the datasets for Harmony integration (PCA is based on SCT)
for (i in seq_along(raw_samples)) {
  DefaultAssay(raw_samples[[i]]) <- "SCT"
  raw_samples[[i]] <- RunPCA(raw_samples[[i]], features = features, verbose = FALSE)
}

# Step 6: Merge the PCA results into a single Seurat object while ensuring cell names are unique
combined <- merge(raw_samples[[1]], y = raw_samples[-1], add.cell.id = sample_ids, merge.data = TRUE)

# Step 7: Run PCA on the combined object using SCT assay
DefaultAssay(combined) <- "SCT"
combined <- RunPCA(combined, features = features, verbose = FALSE)

# Step 8: Run Harmony for integration based on SCT assay PCA
combined <- RunHarmony(
  object = combined,               
  group.by.vars = "orig.ident",    
  dims = 1:30                     
)

# Step 9: Proceed with further analysis using Harmony reduction
combined <- RunUMAP(combined, reduction = "harmony", dims = 1:30, reduction.name = "umap.rna", reduction.key = "rnaUMAP_", verbose = FALSE)

# Step 10: Run TF-IDF and SVD for ATAC assay
DefaultAssay(combined) <- "ATAC"
combined <- RunTFIDF(combined)
combined <- FindTopFeatures(combined, min.cutoff = 'q0')
combined <- RunSVD(combined)

# Step 11: Run UMAP for ATAC using LSI reduction
combined <- RunUMAP(combined, dims = 2:50, reduction = "lsi", reduction.name = "umap.atac", reduction.key = "atacUMAP_", verbose = FALSE)

# Step 12: Run WNN (Weighted Nearest Neighbor) analysis using Harmony and LSI embeddings
combined <- FindMultiModalNeighbors(combined, reduction.list = list("harmony", "lsi"), dims.list = list(1:30, 2:50))
final_combined <- RunUMAP(combined, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

# Set the identities to a meaningful metadata variable, e.g., sample ID or cell type
Idents(final_combined) <- "orig.ident"

# Plotting UMAPs
p1 <- DimPlot(final_combined, reduction = "umap.rna", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(final_combined, reduction = "umap.atac", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(final_combined, reduction = "wnn.umap", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("WNN")
combined_plot <- p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
print(combined_plot)

# FeaturePlot with correct default assay set back to "SCT"
DefaultAssay(final_combined) <- "SCT"
p6 <- FeaturePlot(final_combined, features = c("Apoe", "Gfap", "Acat1"), reduction = 'wnn.umap', max.cutoff = 3, ncol = 3)
print(p6)

# Check if any variable features were identified in the SCT assay
cat("Number of Variable Features in Combined Object (SCT):\n")
print(length(VariableFeatures(final_combined)))

# If no variable features, manually set them
if (length(VariableFeatures(final_combined)) == 0) {
  cat("No variable features identified. Aggregating manually...\n")
  
  # Manually aggregate variable features
  sct_variable_features <- unique(unlist(lapply(final_combined@assays$SCT@SCTModel.list, VariableFeatures)))
  
  # Set the aggregated variable features to the SCT assay in the combined object
  VariableFeatures(final_combined) <- sct_variable_features
  
  cat("Number of Variable Features after manual aggregation:\n")
  print(length(VariableFeatures(final_combined)))
}

# Identify the top 10 variable features in SCT assay for visualization
top10 <- head(VariableFeatures(final_combined), 10)

# Generate heatmap with valid features if available
if (length(top10) > 0) {
  DoHeatmap(final_combined, features = top10, slot = "data") + NoLegend()
} else {
  cat("No valid features found for heatmap generation.\n")
}






######### SINGLE-R INTEGRATION ###########

library(SingleR)
library(celldex)

# Load necessary libraries
library(Seurat)
library(SingleR)
library(celldex)
library(harmony)

# Assuming 'final_combined' is the integrated dataset
# Load the reference dataset for mouse cells
ref <- MouseRNAseqData()

# Step 1: Prepare data for SingleR
DefaultAssay(final_combined) <- "SCT"  # Or your assay of interest (e.g., RNA)
sct_data <- GetAssayData(final_combined, slot = "data", assay = "SCT")

# Step 2: Run SingleR for cell type identification
singler_results <- SingleR(
  test = sct_data, 
  ref = ref, 
  labels = ref$label.main
)

# Step 3: Add SingleR labels to metadata
final_combined$SingleR_labels <- singler_results$labels

# Optionally, add fine labels if available
if (!is.null(singler_results$labels.fine)) {
  final_combined$SingleR_fine_labels <- singler_results$labels.fine
}

# Step 4: Visualization using SingleR labels
Idents(final_combined) <- "SingleR_labels"
p1 <- DimPlot(final_combined, reduction = "wnn.umap", group.by = "SingleR_labels", label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("WNN UMAP with SingleR Cell Types")
print(p1)

library(SeuratDisk)

# Convert Seurat object to AnnData format
SaveH5Seurat(final_combined, filename = "final_combined.h5Seurat", overwrite = TRUE)

# Convert H5Seurat to AnnData (H5AD) format
Convert("final_combined.h5Seurat", dest = "h5ad")






######### DIFFERENTIAL EXPRESSION BETWEEN MOCK AND INFECTED GROUPS (RNA AND ATAC) ###########

next_combined <- final_combined

# Step 1: Verify and Set Identifiers
# Check the unique values in `orig.ident` to ensure they match the sample identifiers
unique(next_combined$orig.ident)

# Set identities based on `orig.ident`
Idents(next_combined) <- "orig.ident"

# Create a new metadata column for group assignment
next_combined$group <- ifelse(next_combined$orig.ident %in% c("m1", "m4"), "Mock", "Infected")

# Verify the new `group` column
table(next_combined$group)

# Step 2: Update Idents for Differential Expression Analysis (RNA)
# Ensure that `group` column is used for differential expression
Idents(next_combined) <- "group"

# Prepare SCT data for differential expression
next_combined <- PrepSCTFindMarkers(next_combined)

# Perform Differential Expression Analysis (RNA)
de_genes <- FindMarkers(next_combined, ident.1 = "Mock", ident.2 = "Infected", assay = "SCT")

# View the top differentially expressed genes
head(de_genes)

# Export RNA differential expression results to CSV
write.csv(de_genes, file = "/Users/stellawroblewski/Desktop/new_differential_expression_results_mockvsinfected.csv", row.names = TRUE)

# Step 3: Update Idents for ATAC Analysis
DefaultAssay(next_combined) <- "ATAC"
Idents(next_combined) <- "group"

# Perform Differential Accessibility Analysis (ATAC)
de_accessibility <- FindMarkers(next_combined, ident.1 = "Mock", ident.2 = "Infected", assay = "ATAC")

# View the top differentially accessible regions
head(de_accessibility)

# Export ATAC differential accessibility results to CSV
write.csv(de_accessibility, file = "/Users/stellawroblewski/Desktop/differential_accessibility_results_mockvsinfected.csv", row.names = TRUE)







########## HEATMAP of top 500 DEGs (differentially expressed between mock and infected)#############

########## HEATMAP of top 500 DEGs (differentially expressed between mock and infected)#############

# Set the default assay to RNA
DefaultAssay(next_combined) <- "RNA"

# Check if 'orig.ident' represents the sample information (i1, i5, m1, m4)
unique(next_combined$orig.ident)

# Use 'orig.ident' as the sample identifier
next_combined$sample <- next_combined$orig.ident

# Ensure the identities are set based on 'sample'
next_combined <- SetIdent(next_combined, value = "sample")

# Identify the top 500 differentially expressed genes based on adjusted p-value
top500_genes <- de_genes %>%
  arrange(p_val_adj) %>%
  head(500) %>%
  pull(gene)

# Calculate average expression for each sample and each of the top 500 genes using 'counts' layer
average_expression <- sapply(unique(next_combined$sample), function(sample) {
  cells <- WhichCells(next_combined, idents = sample)
  expr_data <- FetchData(next_combined, vars = top500_genes, cells = cells, layer = "counts")
  colMeans(expr_data)
})

# Set row names and column names for average_expression
rownames(average_expression) <- top500_genes
colnames(average_expression) <- unique(next_combined$sample)

# Reorder columns to be in the specified order (m1, m4, i1, i5)
average_expression <- average_expression[, c("m1", "m4", "i1", "i5")]

# Rename columns to the specified labels
colnames(average_expression) <- c("Mock Female", "Mock Male", "Infected Male", "Infected Female")

# Standardize the data
average_expression <- t(scale(t(average_expression)))

# Print average expression values to verify
print("Average expression for each of the top 500 genes for each group:")
print(average_expression)

# Create and display the heatmap with enhanced color scale and clustering
pheatmap(
  average_expression,
  cluster_rows = TRUE,  # Cluster rows (genes) to see patterns
  cluster_cols = TRUE,  # Cluster columns (samples) to see group differences
  show_colnames = TRUE,  # Show sample names
  show_rownames = FALSE,  # Do not show gene names
  color = colorRampPalette(c("blue", "white", "red"))(100),  # Enhanced color scale
  main = "Heatmap of Top 500 Differentially Expressed Genes"  # Title of the heatmap
)

# Save the heatmap as a PNG
png("heatmap_top500_genes.png", width = 1200, height = 800)
pheatmap(
  average_expression,
  cluster_rows = TRUE,  
  cluster_cols = TRUE,  
  show_colnames = TRUE,  
  show_rownames = FALSE,  
  color = colorRampPalette(c("blue", "white", "red"))(100),  
  main = "Heatmap of Top 500 Differentially Expressed Genes"  
)
dev.off()






### NEEDS FIXING
########## HEATMAP of top 500 DEGs (differentially expressed between all samples)#############

# Perform differential expression analysis between all samples (i1, i5, m1, m4)
de_genes_all_samples <- FindMarkers(next_combined, ident.1 = "i1", ident.2 = "i5", ident.3 = "m1", ident.4 = "m4", assay = "RNA")

# Ensure 'gene' is a column in de_genes_all_samples
de_genes_all_samples <- de_genes_all_samples %>%
  rownames_to_column(var = "gene")

# Identify the top 500 differentially expressed genes based on adjusted p-value
top500_genes_all_samples <- de_genes_all_samples %>%
  arrange(p_val_adj) %>%
  head(500) %>%
  pull(gene)

# Check if the genes exist in the dataset
existing_genes <- intersect(top500_genes_all_samples, rownames(next_combined))
missing_genes <- setdiff(top500_genes_all_samples, existing_genes)
print(paste("Missing genes:", paste(missing_genes, collapse = ", ")))

# Calculate average expression for each sample and each of the top 500 genes
average_expression <- sapply(unique(next_combined$sample), function(sample) {
  cells <- WhichCells(next_combined, idents = sample)
  expr_data <- FetchData(next_combined, vars = existing_genes, cells = cells, layer = "counts") # Use 'layer' instead of 'slot'
  colMeans(expr_data)
})

# Set row names and column names for average_expression
rownames(average_expression) <- existing_genes
colnames(average_expression) <- unique(next_combined$sample)

# Reorder columns to be in the specified order (m1, m4, i1, i5)
average_expression <- average_expression[, c("m1", "m4", "i1", "i5")]

# Rename columns to the specified labels
colnames(average_expression) <- c("Mock Female", "Mock Male", "Infected Male", "Infected Female")

# Standardize the data
average_expression <- t(scale(t(average_expression)))

# Print average expression values to verify
print("Average expression for each of the top 500 genes for each group:")
print(average_expression)

# Print the p-values for the top 500 genes
print("P-values for the top 500 genes:")
print(de_genes_all_samples %>% filter(gene %in% existing_genes) %>% select(gene, p_val_adj))

# Create and save the heatmap with enhanced color scale and clustering
# Save the heatmap as a PNG file

pheatmap(
  average_expression,
  cluster_rows = TRUE,  # Cluster rows (genes) to see patterns
  cluster_cols = TRUE,  # Cluster columns (samples) to see group differences
  show_colnames = TRUE,  # Show sample names
  show_rownames = FALSE,  # Do not show gene names
  color = colorRampPalette(c("blue", "white", "red"))(100),  # Enhanced color scale
  main = "Heatmap of Top 500 Differentially Expressed Genes (All Samples)"  # Title of the heatmap
)








########### UMAPS #############


# Set the identities to 'orig.ident' so we can subset based on sample IDs
Idents(final_combined) <- "orig.ident"

# Plot UMAPs for each sample individually
p_i1 <- DimPlot(final_combined, reduction = "wnn.umap", cells = WhichCells(final_combined, idents = "i1"), 
                label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("i1 Sample") + NoLegend()

p_i5 <- DimPlot(final_combined, reduction = "wnn.umap", cells = WhichCells(final_combined, idents = "i5"), 
                label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("i5 Sample") + NoLegend()

p_m1 <- DimPlot(final_combined, reduction = "wnn.umap", cells = WhichCells(final_combined, idents = "m1"), 
                label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("m1 Sample") + NoLegend()

p_m4 <- DimPlot(final_combined, reduction = "wnn.umap", cells = WhichCells(final_combined, idents = "m4"), 
                label = TRUE, label.size = 3.5, repel = TRUE) + ggtitle("m4 Sample") + NoLegend()

# Arrange the individual plots side by side
p_all_samples <- (p_i1 | p_i5 | p_m1 | p_m4)
print(p_all_samples)


# Set identities to 'orig.ident' to access individual sample IDs
Idents(final_combined) <- "orig.ident"

# UMAP for Mock samples (m1 and m4) without labels
p_mock <- DimPlot(final_combined, reduction = "wnn.umap", cells.highlight = WhichCells(final_combined, idents = c("m1", "m4")), 
                  cols.highlight = "blue", label = FALSE, repel = TRUE) + 
  ggtitle("Mock Samples (m1, m4) UMAP") + NoLegend()

# UMAP for Infected samples (i1 and i5) without labels
p_infected <- DimPlot(final_combined, reduction = "wnn.umap", cells.highlight = WhichCells(final_combined, idents = c("i1", "i5")), 
                      cols.highlight = "red", label = FALSE, repel = TRUE) + 
  ggtitle("Infected Samples (i1, i5) UMAP") + NoLegend()

# Arrange the UMAPs for Mock and Infected side by side using patchwork
combined_overlayed_plot <- (p_mock | p_infected) & theme(plot.title = element_text(hjust = 0.5))
print(combined_overlayed_plot)






########### KEGG Pathway Analysis ###########


# Load the required libraries
library(clusterProfiler)
library(org.Mm.eg.db)  # Annotation database for mouse
library(enrichplot)
library(dplyr)

# Convert gene symbols to Entrez IDs for pathway analysis
de_genes$gene <- rownames(de_genes)
gene_list <- de_genes %>% arrange(p_val_adj) %>% pull(gene)
gene_list_entrez <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# Run KEGG pathway enrichment analysis
kegg_enrich <- enrichKEGG(gene = gene_list_entrez$ENTREZID, organism = "mmu", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
write.csv(as.data.frame(kegg_enrich), file = "kegg_pathway_analysis_results.csv", row.names = FALSE)

# Plot KEGG pathway enrichment results ordered by gene ratio 
dotplot(kegg_enrich, showCategory = 20) + ggtitle("KEGG Pathway Enrichment Analysis - Ordered by Gene Ratio")







########### Reactome and GO Pathway Analysis ###########

# Load the necessary libraries for Reactome pathway enrichment
library(ReactomePA)

# GO Enrichment Analysis (Biological Process)
go_bp_enrich <- enrichGO(gene = gene_list_entrez$ENTREZID, OrgDb = org.Mm.eg.db, ont = "BP", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
write.csv(as.data.frame(go_bp_enrich), file = "go_bp_enrichment_results.csv", row.names = FALSE)

# Reactome Pathway Enrichment Analysis
reactome_enrich <- enrichPathway(gene = gene_list_entrez$ENTREZID, organism = "mouse", pvalueCutoff = 0.05, qvalueCutoff = 0.05)
write.csv(as.data.frame(reactome_enrich), file = "reactome_pathway_enrichment_results.csv", row.names = FALSE)

# GO Enrichment Analysis (Biological Process) - ordered by Gene Ratio
dotplot(go_bp_enrich, showCategory = 20) +
  aes(x = GeneRatio, y = reorder(Description, GeneRatio), fill = pvalue) +  # Reorder by GeneRatio and set fill to pvalue
  scale_fill_continuous(low = "red", high = "blue") +  # Color based on p-value
  ggtitle("GO Biological Process Enrichment - Ordered by Gene Ratio") +
  xlab("Gene Ratio") +
  ylab("GO Biological Process") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20),  # Increase title size
    axis.title = element_text(size = 16),  # Increase axis title size
    axis.text = element_text(size = 14),   # Increase axis text size
    legend.title = element_text(size = 14),  # Increase legend title size
    legend.text = element_text(size = 12)    # Increase legend text size
  )

# Reactome Pathway Enrichment Analysis - ordered by Gene Ratio
dotplot(reactome_enrich, showCategory = 20) +
  aes(x = GeneRatio, y = reorder(Description, GeneRatio), fill = pvalue) +  # Reorder by GeneRatio and set fill to pvalue
  scale_fill_continuous(low = "red", high = "blue") +  # Color based on p-value
  ggtitle("Reactome Pathway Enrichment - Ordered by Gene Ratio") +
  xlab("Gene Ratio") +
  ylab("Reactome Pathways") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20),  # Increase title size
    axis.title = element_text(size = 16),  # Increase axis title size
    axis.text = element_text(size = 14),   # Increase axis text size
    legend.title = element_text(size = 14),  # Increase legend title size
    legend.text = element_text(size = 12)    # Increase legend text size
  )
















