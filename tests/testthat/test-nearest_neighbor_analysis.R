library(testthat)
library(Seurat)
library(stringr)

# Create mock Seurat object
create_mock_seurat <- function(expr_matrix) {
  seurat_obj <- CreateSeuratObject(counts = expr_matrix)
  return(seurat_obj)
}

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

  # Check that the output is a Seurat object
  expect_s3_class(merged_seurat, "Seurat")

  # Check that only common features are retained
  expect_equal(rownames(merged_seurat), common_genes)

  # Check that cell numbers are correct after merging
  expect_equal(ncol(merged_seurat), 4)

  # Check that a new merged count matrix exists in the RNA assay layers
  expect_true("counts_merge" %in% names(merged_seurat@assays$RNA@layers))

  # Check that the merged count matrix has expected dimensions
  expect_equal(dim(merged_seurat@assays$RNA@layers$counts_merge), c(length(common_genes), 4))
})

