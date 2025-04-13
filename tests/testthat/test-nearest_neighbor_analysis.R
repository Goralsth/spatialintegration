library(testthat)
library(Seurat)
library(stringr)

# Create mock Seurat object
create_mock_seurat <- function(expr_matrix) {
  seurat_obj <- CreateSeuratObject(counts = expr_matrix)
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj)
  return(seurat_obj)
}

test_that("nearest_neighbor_analysis correctly identifies neighbors", {
  # Create common features
  genes <- c("GENEA", "GENEB", "GENEC", "GENED")
  cells <- c("Cell1", "Cell2", "Cell3", "Cell4")

  # Create an expression matrix
  expr_matrix <- matrix(runif(16, min = 1, max = 10),
                        nrow = 4,
                        dimnames = list(genes, cells))

  # Create Seurat object
  seurat_object <- create_mock_seurat(expr_matrix)

  # Run PCA and find neighbors
  seurat_object <- ScaleData(seurat_object)
  seurat_object <- RunPCA(seurat_object, features = VariableFeatures(seurat_object))
  seurat_object <- FindNeighbors(seurat_object, dims = 1:2)

  # Check that nearest neighbors are computed
  expect_true(!is.null(seurat_object@graphs$RNA_snn))
  expect_true(ncol(seurat_object@graphs$RNA_snn) == length(cells))
  expect_true(nrow(seurat_object@graphs$RNA_snn) == length(cells))

  # Check that PCA was run
  expect_true(!is.null(seurat_object@reductions$pca))
})

