library(Seurat)
library(ggplot2)
library(dplyr)
library(viridis)

# Mock RunGSEA function (to simulate the behavior of actual GSEA)
RunGSEA <- function(seurat_object, category, verbose = TRUE) {
  # Add a mock GSEA result (enrichment scores or similar results)
  seurat_object[["gsea"]] <- list(enrichment_scores = rnorm(ncol(seurat_object)))
  return(seurat_object)
}

# Mock GSEAHeatmap function (to simulate the behavior of GSEA heatmap visualization)
GSEAHeatmap <- function(seurat_object, reduction = "nmf", max.terms.per.factor = 3) {
  # Create a simple mock heatmap plot for testing purposes
  p <- ggplot(data.frame(x = 1:10, y = rnorm(10)), aes(x = x, y = y)) +
    geom_tile(aes(fill = y)) +
    scale_fill_viridis_c()
  return(p)
}

# Mock FeaturePlot function (to simulate plotting of features)
FeaturePlot <- function(seurat_object, features, raster = FALSE) {
  # Create a simple mock plot of features
  p <- ggplot(data.frame(x = 1:10, y = rnorm(10)), aes(x = x, y = y)) +
    geom_point() +
    theme_bw() +
    labs(title = "Feature Plot")
  return(p)
}

# Create a mock Seurat object
set.seed(42)
genes <- paste0("Gene", 1:50)
cells <- paste0("Cell", 1:10)
expr_matrix <- matrix(rnorm(50 * 10), nrow = 50, ncol = 10, dimnames = list(genes, cells))
seurat_obj <- CreateSeuratObject(counts = expr_matrix)

# Manually create the 'layers' slot
seurat_obj@assays$RNA@layers <- list()
seurat_obj@assays$RNA@layers[["counts"]] <- expr_matrix
seurat_obj@assays$RNA@layers[["scaled_counts"]] <- expr_matrix * 1.5  # Mock scaled data
seurat_obj@assays$RNA@layers[["normalized_counts"]] <- expr_matrix / rowSums(expr_matrix)

# Add NMF results (simulating with mock data)
seurat_obj@reductions$nmf <- list(
  cell.embeddings = matrix(rnorm(10 * 5), nrow = 10, ncol = 5),  # 10 cells, 5 NMF dimensions
  loadings = matrix(rnorm(50 * 5), nrow = 50, ncol = 5)  # 50 genes, 5 NMF dimensions
)

# Add mock clustering and PCA-related metadata
seurat_obj$seurat_clusters <- factor(sample(1:3, size = 10, replace = TRUE))
seurat_obj$cell_type <- factor(sample(1:3, size = 10, replace = TRUE))  # Mock cell types

# Run the function perform_gsea_and_clustering
seurat_obj <- perform_gsea_and_clustering(
  seurat_object = seurat_obj,
  gsea_category = "C7",
  umap_dims = 1:5,  # Use the first 5 NMF dimensions
  clustering_resolution = 0.6,
  max_terms_per_factor = 5
)

# Check the output to verify the results (e.g., clusters, UMAP)
print(seurat_obj)




