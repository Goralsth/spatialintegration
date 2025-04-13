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

  # Function to subset Seurat object by proportions
  subset_seurat_by_proportions <- function(seurat_obj, cell_type_col, subset_size) {
    # Get original proportions of each cell type
    proportions <- seurat_obj@meta.data %>%
      group_by_at(cell_type_col) %>%
      summarise(original_count = n(), .groups = "drop") %>%
      mutate(proportion = original_count / sum(original_count))

    # Calculate the number of cells to sample from each cell type based on the proportions
    subset_counts <- proportions %>%
      mutate(subset_count = round(proportion * subset_size))

    # Downsample the Seurat object by proportions
    subset_cells <- unlist(lapply(1:nrow(subset_counts), function(i) {
      type_cells <- rownames(seurat_obj@meta.data[seurat_obj@meta.data[[cell_type_col]] == subset_counts$cell_type[i], ])
      sampled_cells <- sample(type_cells, size = subset_counts$subset_count[i], replace = FALSE)
      return(sampled_cells)
    }))

    # Create a new Seurat object with the downsampled cells
    subset_seurat_obj <- subset(seurat_obj, cells = subset_cells)

    # Add the 'original_count' and 'downsampled_count' to the metadata
    subset_seurat_obj$original_count <- ncol(seurat_obj)
    subset_seurat_obj$downsampled_count <- ncol(subset_seurat_obj)

    return(subset_seurat_obj)
  }

  # Test the subset_seurat_by_proportions function with a specified subset size
  subset_size <- 15
  result <- subset_seurat_by_proportions(seurat_obj, cell_type_col = "cell_type", subset_size = subset_size)

  # Check that the subset Seurat object contains the correct number of cells
  expect_equal(ncol(result), subset_size)

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
  expect_true("original_count" %in% colnames(result@meta.data))
  expect_true("downsampled_count" %in% colnames(result@meta.data))

  # Check that the report is printed correctly
  expect_output(print(result@meta.data), "cell_type")
})
