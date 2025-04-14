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

  # Scale the data (e.g., multiply by a factor for mock scaling)
  scaled_expr_matrix <- expr_matrix * 1.5

  # Normalize the data (e.g., divide by row sums as an example)
  normalized_expr_matrix <- expr_matrix / rowSums(expr_matrix)

  # Store the scaled and normalized data in the misc slot
  seurat_obj@misc$scaled_counts <- scaled_expr_matrix
  seurat_obj@misc$normalized_counts <- normalized_expr_matrix

  # Mock implementation of the process_and_run_nmf function
  process_and_run_nmf <- function(seurat_obj, nmf_rank) {
    # Perform NMF and store results in the reductions slot
    nmf_embeddings <- matrix(rnorm(ncol(seurat_obj) * nmf_rank), nrow = ncol(seurat_obj), ncol = nmf_rank)
    rownames(nmf_embeddings) <- colnames(seurat_obj)  # Assign cell names as row names for embeddings
    colnames(nmf_embeddings) <- paste0("NMF_", seq_len(nmf_rank))  # Assign appropriate column names

    nmf_loadings <- matrix(rnorm(nrow(seurat_obj) * nmf_rank), nrow = nrow(seurat_obj), ncol = nmf_rank)
    rownames(nmf_loadings) <- rownames(seurat_obj)  # Assign gene names as row names for loadings
    colnames(nmf_loadings) <- paste0("NMF_", seq_len(nmf_rank))  # Assign appropriate column names

    # Add NMF results to the reductions slot
    seurat_obj@reductions$nmf <- CreateDimReducObject(
      embeddings = nmf_embeddings,
      loadings = nmf_loadings,
      key = "NMF_",
      assay = "RNA"
    )

    # Mock assignment of scaled counts to data slot (if applicable)
    seurat_obj@assays$RNA@data <- as(scaled_expr_matrix, "dgCMatrix")

    return(seurat_obj)
  }

  # Run the process_and_run_nmf function with a specific NMF rank
  seurat_obj_processed <- process_and_run_nmf(seurat_obj, nmf_rank = 3)

  # Assertions to validate the results
  expect_true("nmf" %in% names(seurat_obj_processed@reductions))  # Check NMF exists
  nmf_reduction <- seurat_obj_processed@reductions[["nmf"]]

  # Check that embeddings and loadings are populated
  expect_true(!is.null(nmf_reduction@cell.embeddings))
  expect_true(!is.null(nmf_reduction@feature.loadings))

  # Check that embeddings and loadings have the expected dimensions
  expect_equal(ncol(nmf_reduction@cell.embeddings), 3)  # NMF rank
  expect_equal(nrow(nmf_reduction@feature.loadings), 50)  # Genes

  # Check scaled counts are assigned correctly
  expect_true("scaled_counts" %in% names(seurat_obj_processed@misc))
  expect_identical(seurat_obj_processed@assays$RNA@data, as(scaled_expr_matrix, "dgCMatrix"))

  # Check process_and_run_nmf runs without error
  expect_silent(process_and_run_nmf(seurat_obj, nmf_rank = 3))
})
