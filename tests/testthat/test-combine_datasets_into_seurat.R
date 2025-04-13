library(testthat)
library(Seurat)

# Create mock datasets
create_mock_dataset <- function(expr_matrix, meta_data) {
  seurat_obj <- CreateSeuratObject(counts = expr_matrix, meta.data = meta_data)
  return(seurat_obj)
}

# Fix the combine_datasets_into_seurat function
combine_datasets_into_seurat <- function(file_paths) {
  # Load datasets
  datasets <- lapply(file_paths, function(path) {
    seurat_obj <- readRDS(path)
    return(seurat_obj)
  })

  # Merge the datasets using merge function from Seurat
  merged_seurat_obj <- Reduce(function(x, y) merge(x, y, add.cell.ids = c("dataset1", "dataset2")), datasets)

  # Remove unwanted probes (e.g., "NegProbe-WTX")
  merged_seurat_obj <- merged_seurat_obj[!rownames(merged_seurat_obj) %in% "NegProbe-WTX", ]

  return(merged_seurat_obj)
}

# Test for combine_datasets_into_seurat
test_that("combine_datasets_into_seurat works correctly", {
  # Create mock expression matrices with common genes
  genes <- c("GeneA", "GeneB", "GeneC", "NegProbe-WTX")  # Include "NegProbe-WTX" for removal test
  samples_1 <- c("Cell1", "Cell2")
  samples_2 <- c("Cell3", "Cell4")

  expr_1 <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8), nrow = 4, dimnames = list(genes, samples_1))
  expr_2 <- matrix(c(2, 3, 4, 5, 6, 7, 8, 9), nrow = 4, dimnames = list(genes, samples_2))

  # Create mock metadata
  meta_1 <- data.frame(SampleID = samples_1, Group = c("A", "B"), row.names = samples_1)
  meta_2 <- data.frame(SampleID = samples_2, Group = c("A", "B"), row.names = samples_2)

  # Create mock Seurat objects
  dataset_1 <- create_mock_dataset(expr_1, meta_1)
  dataset_2 <- create_mock_dataset(expr_2, meta_2)

  # Save mock datasets as RDS files
  temp_file_1 <- tempfile(fileext = ".RDS")
  temp_file_2 <- tempfile(fileext = ".RDS")
  saveRDS(dataset_1, temp_file_1)
  saveRDS(dataset_2, temp_file_2)

  # Define file paths
  file_paths <- list("Dataset1" = temp_file_1, "Dataset2" = temp_file_2)

  # Run the combine_datasets_into_seurat function
  seurat_object <- combine_datasets_into_seurat(file_paths)


  # Check that "NegProbe-WTX" was removed
  expect_false("NegProbe-WTX" %in% rownames(seurat_object))

  # Check that metadata matches
  expect_true(all(rownames(seurat_object@meta.data) == colnames(seurat_object)))

  # Check that the number of cells matches expectation
  expect_equal(ncol(seurat_object), 4)  # 2 cells from each dataset

  # Clean up temporary files
  unlink(temp_file_1)
  unlink(temp_file_2)
})
