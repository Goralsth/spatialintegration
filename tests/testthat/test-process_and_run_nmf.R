# Load required libraries
library(testthat)
library(Seurat)
library(Matrix)

# Define the test
test_that("process_and_run_nmf works correctly", {

  # Create a mock expression matrix for the Seurat object
  genes <- paste0("Gene", 1:50)
  cells <- paste0("Cell", 1:10)

  # Create a random expression matrix for the Seurat object
  expr_matrix <- matrix(rnorm(500), nrow = 50, ncol = 10, dimnames = list(genes, cells))

  # Create a mock Seurat object
  seurat_obj <- CreateSeuratObject(counts = expr_matrix)

  # Test the process_and_run_nmf function with a specified rank for NMF
  seurat_obj_processed <- process_and_run_nmf(seurat_obj, nmf_rank = 5)

  # Check that NMF results are assigned to the reductions slot
  expect_true("nmf" %in% names(seurat_obj_processed@reductions))  # Ensure 'nmf' reduction exists

  # Check if the NMF reduction contains embeddings and loadings
  nmf_reduction <- seurat_obj_processed@reductions[["nmf"]]
  expect_true("embeddings" %in% names(nmf_reduction))  # Ensure 'embeddings' exist in NMF
  expect_true("loadings" %in% names(nmf_reduction))    # Ensure 'loadings' exist in NMF

  # Check that the embeddings and loadings are not empty
  expect_gt(nrow(nmf_reduction$embeddings), 0)  # Ensure embeddings have rows
  expect_gt(ncol(nmf_reduction$embeddings), 0)  # Ensure embeddings have columns
  expect_gt(nrow(nmf_reduction$loadings), 0)    # Ensure loadings have rows
  expect_gt(ncol(nmf_reduction$loadings), 0)    # Ensure loadings have columns

  # Check if the correct number of factors (k) were used
  expect_equal(ncol(nmf_reduction$embeddings), 5)  # Ensure 5 factors (k = 5)

  # Check for absence of errors during the process
  expect_silent(process_and_run_nmf(seurat_obj, nmf_rank = 5))
})

