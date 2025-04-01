# Load required libraries
library(testthat)
library(Seurat)
library(ggplot2)
library(reshape2)

# Define the test
test_that("plot_nmf_correlation_heatmap works correctly", {

  # Create mock expression matrix
  genes <- paste0("Gene", 1:50)
  cells <- paste0("Cell", 1:10)
  expr_matrix <- matrix(rnorm(500), nrow = 50, ncol = 10, dimnames = list(genes, cells))

  # Create a mock Seurat object
  seurat_obj <- CreateSeuratObject(counts = expr_matrix)

  # Add mock NMF reduction to the Seurat object
  seurat_obj@reductions$nmf <- CreateDimReducObject(
    embeddings = matrix(rnorm(10 * 5), nrow = 10, ncol = 5),  # 10 cells, 5 NMF factors
    key = "NMF_", assay = "RNA"
  )

  # Add mock metadata for 'subclass_label' and 'segment'
  seurat_obj@meta.data$subclass_label <- sample(c("TypeA", "TypeB", "TypeC"), 10, replace = TRUE)
  seurat_obj@meta.data$segment <- sample(c("Segment1", "Segment2"), 10, replace = TRUE)

  # Test the plot_nmf_correlation_heatmap function
  output_file <- tempfile(fileext = ".png")

  p <- plot_nmf_correlation_heatmap(seurat_object = seurat_obj, output_file = output_file)

  # Check that the result is a ggplot object
  expect_s3_class(p, "gg")  # Check that p is a ggplot object

  # Check if the file is created
  expect_true(file.exists(output_file))

  # Clean up: remove the generated file after the test
  unlink(output_file)
})

