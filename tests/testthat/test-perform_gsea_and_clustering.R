# Load libraries
library(Seurat)
library(ggplot2)
library(dplyr)

# Mock RunGSEA function
RunGSEA <- function(seurat_object, category, verbose = TRUE) {
  seurat_object@misc$gsea <- list(enrichment_scores = rnorm(ncol(seurat_object)))
  return(seurat_object)
}

# Mock perform_gsea_and_clustering function
perform_gsea_and_clustering <- function(seurat_object, gsea_category, umap_dims, clustering_resolution, max_terms_per_factor) {
  # Run GSEA (mock)
  seurat_object <- RunGSEA(seurat_object, category = gsea_category)

  # Perform clustering (mock UMAP embedding from NMF reduction)
  nmf_embedding <- matrix(rnorm(ncol(seurat_object) * length(umap_dims)), ncol = length(umap_dims))
  rownames(nmf_embedding) <- colnames(seurat_object)

  # Add UMAP reduction
  seurat_object@reductions$umap <- CreateDimReducObject(embeddings = nmf_embedding, key = "UMAP_", assay = "RNA")

  # Perform clustering
  seurat_object <- FindNeighbors(seurat_object, reduction = "umap", dims = 1:length(umap_dims))
  seurat_object <- FindClusters(seurat_object, resolution = clustering_resolution)

  return(seurat_object)
}

# Create a mock Seurat object
set.seed(42)
genes <- paste0("Gene", 1:50)
cells <- paste0("Cell", 1:10)
expr_matrix <- matrix(rnorm(50 * 10), nrow = 50, ncol = 10, dimnames = list(genes, cells))
seurat_obj <- CreateSeuratObject(counts = expr_matrix)

# Add mock NMF reduction
nmf_embedding <- matrix(rnorm(ncol(seurat_obj) * 5), nrow = ncol(seurat_obj), ncol = 5)
rownames(nmf_embedding) <- colnames(seurat_obj)
seurat_obj@reductions$nmf <- CreateDimReducObject(embeddings = nmf_embedding, key = "NMF_", assay = "RNA")

# Run the perform_gsea_and_clustering function
seurat_obj <- perform_gsea_and_clustering(
  seurat_object = seurat_obj,
  gsea_category = "C7",
  umap_dims = 1:5,
  clustering_resolution = 0.6,
  max_terms_per_factor = 5
)

# Check the clustering and dimensionality reduction results
print(seurat_obj)





