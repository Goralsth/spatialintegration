# Load required libraries
library(testthat)
library(Seurat)
library(dplyr)

# Define the test
test_that("subset_seurat_by_proportions works correctly", {

  # Create a mock expression matrix for the Seurat object
  genes <- paste0("Gene", 1:50)
  cells <- paste0("Cell", 1:20)

  # Create a mock Seurat object with 3 cell types in metadata
  expr_matrix <- matrix(rnorm(1000), nrow = 50, ncol = 20, dimnames = list(genes, cells))
  seurat_obj <- CreateSeuratObject(counts = expr_matrix)

  # Assign cell types to the metadata column
  seurat_obj$cell_type <- rep(c("Type1", "Type2", "Type3"), length.out = 20)

  # Test the subset_seurat_by_proportions function with a specified subset size
  subset_size <- 14
  result <- subset_seurat_by_proportions(seurat_obj, cell_type_col = "cell_type", subset_size = subset_size)

  # Check that the proportions of cell types in the subset are similar to the original proportions
  original_proportions <- seurat_obj@meta.data %>%
    group_by(cell_type) %>%
    summarise(original_count = n(), .groups = "drop") %>%
    mutate(proportion = original_count / sum(original_count))

  subset_proportions <- result@meta.data %>%
    group_by(cell_type) %>%
    summarise(downsampled_count = n(), .groups = "drop") %>%
    mutate(proportion = downsampled_count / sum(downsampled_count))

  # Compare original and subset proportions
  for (cell_type in unique(original_proportions$cell_type)) {
    original_proportion <- original_proportions$proportion[original_proportions$cell_type == cell_type]
    subset_proportion <- subset_proportions$proportion[subset_proportions$cell_type == cell_type]

    # Allow a small difference due to random sampling
    expect_true(abs(original_proportion - subset_proportion) < 0.1)
  }

  # Check that the report contains the expected columns
  expect_true("cell_type" %in% colnames(result@meta.data))

  # Check that the report is printed correctly
  expect_output(print(result@meta.data), "cell_type")
})


