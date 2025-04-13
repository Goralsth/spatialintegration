# Load required libraries
library(testthat)
library(Seurat)
library(Matrix)

# Define the test for the process_and_run_nmf function
test_that("process_and_run_nmf works correctly", {

  # Create a mock expression matrix for testing
  genes <- paste0("Gene", 1:50)  # 50 genes
  cells <- paste0("Cell", 1:10)  # 10 cells
  expr_matrix <- matrix(rnorm(50 * 10), nrow = 50, ncol = 10, dimnames = list(genes, cells))

  # Create a Seurat object with the expression matrix
  seurat_obj <- CreateSeuratObject(counts = expr_matrix)

  # Manually create the 'layers' slot if it doesn't exist
  if (is.null(seurat_obj@assays$RNA@layers)) {
    seurat_obj@assays$RNA@layers <- list()  # Initialize layers slot
  }

  # Add the counts matrix to the 'counts' layer
  seurat_obj@assays$RNA@layers[["counts"]] <- expr_matrix

  # Create a scaled version of the data (for this example, we multiply by a scaling factor)
  scaled_expr_matrix <- expr_matrix * 1.5  # This would typically be scaled data
  seurat_obj@assays$RNA@layers[["scaled_counts"]] <- scaled_expr_matrix

  # Add a normalized version of the data (here, we're dividing by row sums as an example)
  normalized_expr_matrix <- expr_matrix / rowSums(expr_matrix)
  seurat_obj@assays$RNA@layers[["normalized_counts"]] <- normalized_expr_matrix

  # Run the process_and_run_nmf function with specific rank for NMF
  seurat_obj_processed <- process_and_run_nmf(seurat_obj, nmf_rank = 3)

  # Check that NMF results are assigned to the reductions slot
  expect_true("nmf" %in% names(seurat_obj_processed@reductions))

  # Check if the NMF reduction contains embeddings and loadings
  nmf_reduction <- seurat_obj_processed@reductions[["nmf"]]
  expect_true("embeddings" %in% names(nmf_reduction))
  expect_true("loadings" %in% names(nmf_reduction))

  # Check that the embeddings and loadings are not empty
  expect_gt(nrow(nmf_reduction$embeddings), 0)
  expect_gt(ncol(nmf_reduction$embeddings), 0)
  expect_gt(nrow(nmf_reduction$loadings), 0)
  expect_gt(ncol(nmf_reduction$loadings), 0)

  # Check if the correct number of factors (k) were used in the NMF model
  expect_equal(ncol(nmf_reduction$embeddings), 3)

  # Verify that the scaled counts are assigned to the correct layer
  expect_true("scaled_counts" %in% names(seurat_obj_processed@assays$RNA@layers))

  # Check if the data layer contains the scaled counts
  expect_identical(seurat_obj_processed@assays$RNA@layers[["data"]],
                   seurat_obj_processed@assays$RNA@layers[["scaled_counts"]])

  # Check if no errors are thrown during the process
  expect_silent(process_and_run_nmf(seurat_obj, nmf_rank = 3))
})
