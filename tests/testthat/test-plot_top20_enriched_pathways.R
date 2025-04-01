# Load required libraries
library(testthat)
library(Seurat)
library(ggplot2)
library(dplyr)

# Define the test
test_that("plot_top20_enriched_pathways works correctly", {

  # Create mock expression matrix for Seurat and Allen objects
  genes <- paste0("Gene", 1:50)
  pathways <- paste0("Pathway", 1:50)
  cells <- paste0("Cell", 1:10)

  # Create a random expression matrix for both Seurat and Allen
  expr_matrix <- matrix(rnorm(500), nrow = 50, ncol = 10, dimnames = list(genes, cells))

  # Create mock Seurat and Allen Seurat objects
  seurat_obj <- CreateSeuratObject(counts = expr_matrix)
  allen_obj <- CreateSeuratObject(counts = expr_matrix)

  # Add mock NMF reduction with GSEA results
  nmf_seurat <- matrix(rnorm(10 * 5), nrow = 10, ncol = 5)  # 10 cells, 5 factors
  nmf_allen <- matrix(rnorm(10 * 5), nrow = 10, ncol = 5)

  seurat_obj@reductions$nmf <- CreateDimReducObject(embeddings = nmf_seurat, key = "NMF_", assay = "RNA")
  allen_obj@reductions$nmf <- CreateDimReducObject(embeddings = nmf_allen, key = "NMF_", assay = "RNA")

  # Mock GSEA NES data (Normalized Enrichment Scores)
  seurat_obj@reductions$nmf@misc$gsea$nes <- matrix(rnorm(50 * 5), nrow = 50, ncol = 5)
  allen_obj@reductions$nmf@misc$gsea$nes <- matrix(rnorm(50 * 5), nrow = 50, ncol = 5)

  # Set rownames for NES matrices to be the pathway names
  rownames(seurat_obj@reductions$nmf@misc$gsea$nes) <- pathways
  rownames(allen_obj@reductions$nmf@misc$gsea$nes) <- pathways

  # Test the plot_top20_enriched_pathways function
  output_file <- tempfile(fileext = ".pdf")

  p <- plot_top20_enriched_pathways(
    seurat_object = seurat_obj,
    allen_object = allen_obj,
    nmf_factor = 1,
    output_file = output_file
  )

  # Check that the result is a ggplot object
  expect_s3_class(p, "gg")  # Check that p is a ggplot object

  # Check if the file is created
  expect_true(file.exists(output_file))

  # Clean up: remove the generated file after the test
  unlink(output_file)
})


