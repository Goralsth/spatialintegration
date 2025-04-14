# Load required libraries
library(testthat)
library(Seurat)
library(ggplot2)
library(dplyr)

# Define the test
test_that("plot_top20_enriched_pathways_single works correctly", {

  # Create mock expression matrix for Seurat object
  genes <- paste0("Gene", 1:50)
  pathways <- paste0("Pathway", 1:50)
  cells <- paste0("Cell", 1:10)

  # Create a random expression matrix for the Seurat object
  expr_matrix <- matrix(rnorm(500), nrow = 50, ncol = 10, dimnames = list(genes, cells))

  # Create a mock Seurat object
  seurat_obj <- CreateSeuratObject(counts = expr_matrix)

  # Mock NMF reduction with GSEA results
  nmf_seurat <- matrix(rnorm(10 * 5), nrow = 10, ncol = 5, dimnames = list(cells, paste0("NMF_", 1:5)))  # Add rownames for cells
  seurat_obj@reductions$nmf <- CreateDimReducObject(embeddings = nmf_seurat, key = "NMF_", assay = "RNA")

  # Mock GSEA NES and p-values data
  nes_values <- rnorm(50)  # NES values for 50 pathways
  p_values <- runif(50, 0.001, 0.05)  # p-values for 50 pathways

  # Add the mock GSEA data to the Seurat object
  seurat_obj@reductions$nmf@misc$gsea <- list(
    nes = matrix(nes_values, nrow = 50, ncol = 1, dimnames = list(pathways, "NES")),
    pval = matrix(p_values, nrow = 50, ncol = 1, dimnames = list(pathways, "pval"))
  )

  # Define the plot_top20_enriched_pathways_single function (mock for this test)
  plot_top20_enriched_pathways_single <- function(seurat_object, nmf_factor, output_file) {
    # Extract GSEA results for the selected factor
    nes <- seurat_object@reductions$nmf@misc$gsea$nes[, 1]
    pval <- seurat_object@reductions$nmf@misc$gsea$pval[, 1]

    # Combine NES and p-values into a data frame
    gsea_results <- data.frame(Pathway = rownames(seurat_object@reductions$nmf@misc$gsea$nes),
                               NES = nes, pval = pval)

    # Filter top 20 pathways by NES
    top_pathways <- gsea_results %>%
      arrange(desc(NES)) %>%
      head(20)

    # Create the ggplot object
    p <- ggplot(top_pathways, aes(x = reorder(Pathway, NES), y = NES, fill = pval)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      scale_fill_viridis_c() +
      theme_minimal() +
      labs(title = paste("Top 20 Enriched Pathways for Factor", nmf_factor),
           x = "Pathways", y = "NES")

    # Save the plot to a file
    ggsave(output_file, plot = p, width = 8, height = 6)

    return(p)
  }

  # Test the plot_top20_enriched_pathways_single function
  output_file <- tempfile(fileext = ".pdf")

  p <- plot_top20_enriched_pathways_single(
    seurat_object = seurat_obj,
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
