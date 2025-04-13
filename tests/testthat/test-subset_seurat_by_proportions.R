library(testthat)
library(Seurat)
library(dplyr)
library(progress)

# Create a mock Seurat object for testing
create_mock_seurat <- function(num_cells = 100, num_genes = 50, cell_types = c("TypeA", "TypeB")) {
  # Generate random expression data
  expr_matrix <- matrix(rnorm(num_cells * num_genes), nrow = num_genes, ncol = num_cells)
  rownames(expr_matrix) <- paste0("Gene", 1:num_genes)
  colnames(expr_matrix) <- paste0("Cell", 1:num_cells)

  # Generate random cell type metadata
  cell_type_col <- sample(cell_types, num_cells, replace = TRUE)
  meta_data <- data.frame(cell_type = cell_type_col, row.names = colnames(expr_matrix))

  # Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = expr_matrix, meta.data = meta_data)
  return(seurat_obj)
}

test_that("subset_seurat_by_proportions correctly subsets Seurat object while maintaining cell type proportions", {
  # Create a mock Seurat object with 100 cells and 2 cell types (TypeA and TypeB)
  seurat_obj <- create_mock_seurat(num_cells = 100, num_genes = 50, cell_types = c("TypeA", "TypeB"))

  # Set the desired subset size (e.g., 50 cells)
  subset_size <- 50

  # Call the function to subset the Seurat object
  result <- subset_seurat_by_proportions(
    seurat_object = seurat_obj,
    cell_type_col = "cell_type",
    subset_size = subset_size
  )

  # Extract the resulting subsetted Seurat object and the report
  seurat_subset <- result
  report <- result@meta.data %>%
    group_by(cell_type) %>%
    summarise(count = n())

  # Check that the number of cells in the subset equals the desired subset size
  expect_equal(ncol(seurat_subset), subset_size)

  # Check that the proportions of the cell types are approximately preserved in the subset
  original_proportions <- table(seurat_obj@meta.data$cell_type) / ncol(seurat_obj)
  subset_proportions <- table(seurat_subset@meta.data$cell_type) / ncol(seurat_subset)

  for (cell_type in unique(seurat_obj@meta.data$cell_type)) {
    expect_equal(subset_proportions[cell_type], original_proportions[cell_type], tolerance = 0.1)
  }

  # Check that the report shows the original and downsampled counts for each cell type
  expect_true("original_count" %in% colnames(report))
  expect_true("downsampled_count" %in% colnames(report))

  # Check that the total downsampled count matches the subset size
  total_downsampled <- sum(report$downsampled_count)
  expect_equal(total_downsampled, subset_size)
})

