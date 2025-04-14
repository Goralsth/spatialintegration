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
  nmf_seurat <- matrix(rnorm(10 * 5), nrow = 10, ncol = 5, dimnames = list(cells, paste0("NMF_", 1:5)))  # Add rownames for cells
  nmf_allen <- matrix(rnorm(10 * 5), nrow = 10, ncol = 5, dimnames = list(cells, paste0("NMF_", 1:5)))  # Add rownames for cells

  seurat_obj@reductions$nmf <- CreateDimReducObject(embeddings = nmf_seurat, key = "NMF_", assay = "RNA")
  allen_obj@reductions$nmf <- CreateDimReducObject(embeddings = nmf_allen, key = "NMF_", assay = "RNA")

  # Mock GSEA NES data (Normalized Enrichment Scores)
  seurat_obj@reductions$nmf@misc$gsea$nes <- matrix(rnorm(50 * 5), nrow = 50, ncol = 5, dimnames = list(pathways, paste0("NMF_", 1:5)))
  allen_obj@reductions$nmf@misc$gsea$nes <- matrix(rnorm(50 * 5), nrow = 50, ncol = 5, dimnames = list(pathways, paste0("NMF_", 1:5)))

  # Define the plot_top20_enriched_pathways function (mock for this test)
  plot_top20_enriched_pathways <- function(seurat_object, allen_object, nmf_factor, output_file) {
    # Extract GSEA NES data for the selected NMF factor
    nes_seurat <- seurat_object@reductions$nmf@misc$gsea$nes[, nmf_factor]
    nes_allen <- allen_object@reductions$nmf@misc$gsea$nes[, nmf_factor]

    # Combine NES data into a data frame
    gsea_results <- data.frame(
      Pathway = rownames(seurat_object@reductions$nmf@misc$gsea$nes),
      NES_Seurat = nes_seurat,
      NES_Allen = nes_allen
    )

    # Filter top 20 pathways by NES_Seurat
    top_pathways <- gsea_results %>%
      arrange(desc(NES_Seurat)) %>%
      head(20)

    # Create the ggplot object
    p <- ggplot(top_pathways, aes(x = reorder(Pathway, NES_Seurat), y = NES_Seurat)) +
      geom_bar(stat = "identity", aes(fill = NES_Allen)) +
      coord_flip() +
      scale_fill_viridis_c() +
      theme_minimal() +
      labs(title = paste("Top 20 Enriched Pathways for Factor", nmf_factor),
           x = "Pathways", y = "NES (Seurat)")

    # Save the plot to a file
    ggsave(output_file, plot = p, width = 8, height = 6)

    return(p)
  }

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

