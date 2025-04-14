library(Seurat)
library(stringr)

# Create mock Seurat object
create_mock_seurat <- function(expr_matrix) {
  seurat_obj <- CreateSeuratObject(counts = expr_matrix)
  return(seurat_obj)
}

# Mock implementation of merge_allen_geomx
merge_allen_geomx <- function(seurat_obj1, seurat_obj2) {
  # Merge the two Seurat objects
  merged_seurat <- merge(seurat_obj1, y = seurat_obj2, merge.data = TRUE)

  # Retrieve the merged count matrix
  merged_counts <- GetAssayData(merged_seurat, slot = "counts")

  # Optionally store the merged count matrix in the 'misc' slot
  merged_seurat@misc$merged_counts <- merged_counts

  return(merged_seurat)
}

# Test the merge function
test_that("merge_allen_geomx correctly merges two Seurat objects", {
  # Create common features
  common_genes <- c("GENEA", "GENEB", "GENEC")
  extra_genes_1 <- c("GENEX")
  extra_genes_2 <- c("GENEY")

  # Create expression matrices with some unique and common genes
  expr_1 <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8),
                   nrow = 4,
                   dimnames = list(c(common_genes, extra_genes_1), c("Cell1", "Cell2")))

  expr_2 <- matrix(c(2, 3, 4, 5, 6, 7, 8, 9),
                   nrow = 4,
                   dimnames = list(c(common_genes, extra_genes_2), c("Cell3", "Cell4")))

  # Create Seurat objects
  seurat_object <- create_mock_seurat(expr_1)
  allen_object <- create_mock_seurat(expr_2)

  # Run the merge function
  merged_seurat <- merge_allen_geomx(seurat_object, allen_object)


  # Check that the features include both common and unique genes
  all_genes <- unique(c(rownames(seurat_object), rownames(allen_object)))
  expect_true(all(all_genes %in% rownames(merged_seurat)))

  # Check that the cell count is correct after merging
  expect_equal(ncol(merged_seurat), 4)

  # Check that the merged count matrix matches the expected dimensions
  merged_counts <- GetAssayData(merged_seurat, slot = "counts")
  expect_equal(dim(merged_counts), c(length(all_genes), 4))
})

