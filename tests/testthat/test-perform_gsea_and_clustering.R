library(testthat)
library(Seurat)
library(ggplot2)
library(viridis)
library(dplyr)

# Mock function for RunGSEA (as we don't have the actual function)
RunGSEA <- function(seurat_object, category, verbose = FALSE) {
  # Simulate adding some GSEA-related metadata
  seurat_object@meta.data$GSEA_Dummy <- sample(letters, ncol(seurat_object), replace = TRUE)
  return(seurat_object)
}

# Mock function for GSEAHeatmap (returns a dummy plot)
GSEAHeatmap <- function(seurat_object, reduction, max.terms.per.factor) {
  ggplot() + ggtitle("Mock GSEA Heatmap")
}

# Function for performing GSEA and clustering
perform_gsea_and_clustering <- function(seurat_object, gsea_category, umap_dims, clustering_resolution, max_terms_per_factor) {
  # Run GSEA
  seurat_object <- RunGSEA(seurat_object, category = gsea_category)

  # Perform clustering based on UMAP dimensions
  seurat_object <- FindNeighbors(seurat_object, dims = umap_dims)
  seurat_object <- FindClusters(seurat_object, resolution = clustering_resolution)

  # Calculate UMAP embedding for visualization
  seurat_object <- RunUMAP(seurat_object, dims = umap_dims)

  # Generate GSEA heatmap and save it in the 'misc' slot
  heatmap <- GSEAHeatmap(seurat_object, reduction = "umap", max.terms.per.factor = max_terms_per_factor)
  if (is.null(seurat_object@misc)) {
    seurat_object@misc <- list()
  }
  seurat_object@misc$plots <- list(heatmap)

  return(seurat_object)
}

# Test function for perform_gsea_and_clustering
test_that("perform_gsea_and_clustering works as expected", {
  # Create a mock expression matrix
  genes <- paste0("Gene", 1:50)
  cells <- paste0("Cell", 1:10)
  expr_matrix <- matrix(rnorm(500), nrow = 50, ncol = 10, dimnames = list(genes, cells))

  # Create a mock Seurat object
  seurat_obj <- CreateSeuratObject(counts = expr_matrix)

  # Add a mock 'nmf' reduction (assuming it is an NMF-based reduction)
  seurat_obj@reductions$nmf <- CreateDimReducObject(
    embeddings = matrix(rnorm(10 * 20), nrow = 10, ncol = 20),
    key = "NMF_", assay = "RNA"
  )

  # Add mock metadata for 'cell_type' column
  seurat_obj@meta.data$cell_type <- sample(c("TypeA", "TypeB", "TypeC"), 10, replace = TRUE)

  # Run the function perform_gsea_and_clustering
  updated_seurat <- perform_gsea_and_clustering(
    seurat_object = seurat_obj,
    gsea_category = "C7",
    umap_dims = 1:10,
    clustering_resolution = 0.6,
    max_terms_per_factor = 5
  )

  # Check if the output is a Seurat object
  expect_s3_class(updated_seurat, "Seurat")

  # Check if the GSEA column was added to metadata
  expect_true("GSEA_Dummy" %in% colnames(updated_seurat@meta.data))

  # Check if clusters were identified (seurat_clusters should exist)
  expect_true("seurat_clusters" %in% colnames(updated_seurat@meta.data))

  # Check if UMAP was calculated (umap reduction should exist)
  expect_true("umap" %in% names(updated_seurat@reductions))

  # Check if the cluster names were properly set (based on PCA comparison)
  expect_true(!is.null(levels(updated_seurat@meta.data$seurat_clusters)))

  # Check if the cluster visualization (e.g., FeaturePlot) was produced (mocked plot should exist)
  expect_true(inherits(updated_seurat@misc$plots[[1]], "gg")) # Assuming the plot was saved in `misc$plots`
})



