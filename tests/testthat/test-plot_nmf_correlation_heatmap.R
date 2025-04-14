# Load required libraries
library(testthat)
library(Seurat)
library(ggplot2)
library(reshape2)

# Define the mock plot_nmf_correlation_heatmap function
plot_nmf_correlation_heatmap <- function(seurat_object, output_file) {
  # Extract NMF embeddings
  nmf_embeddings <- seurat_object@reductions$nmf@cell.embeddings

  # Calculate correlation matrix
  cor_matrix <- cor(nmf_embeddings)

  # Melt the correlation matrix for ggplot
  cor_melted <- reshape2::melt(cor_matrix)

  # Create heatmap using ggplot
  p <- ggplot(cor_melted, aes(Var1, Var2, fill = value)) +
    geom_tile() +
    scale_fill_viridis_c() +
    theme_minimal() +
    ggtitle("NMF Correlation Heatmap") +
    xlab("NMF Factors") +
    ylab("NMF Factors")

  # Save the plot to a file
  ggsave(output_file, plot = p)

  return(p)
}

# Define the test
test_that("plot_nmf_correlation_heatmap works correctly", {
  # Create mock expression matrix
  genes <- paste0("Gene", 1:50)
  cells <- paste0("Cell", 1:10)
  expr_matrix <- matrix(rnorm(500), nrow = 50, ncol = 10, dimnames = list(genes, cells))

  # Create a mock Seurat object
  seurat_obj <- CreateSeuratObject(counts = expr_matrix)

  # Add mock NMF reduction to the Seurat object
  nmf_embeddings <- matrix(rnorm(10 * 5), nrow = 10, ncol = 5,
                           dimnames = list(cells, paste0("NMF_", 1:5)))
  seurat_obj@reductions$nmf <- CreateDimReducObject(
    embeddings = nmf_embeddings, key = "NMF_", assay = "RNA"
  )

  # Add mock metadata for 'subclass_label' and 'segment'
  seurat_obj@meta.data$subclass_label <- sample(c("TypeA", "TypeB", "TypeC"), 10, replace = TRUE)
  seurat_obj@meta.data$segment <- sample(c("Segment1", "Segment2"), 10, replace = TRUE)

  # Test the plot_nmf_correlation_heatmap function
  output_file <- tempfile(fileext = ".png")

  p <- plot_nmf_correlation_heatmap(seurat_object = seurat_obj, output_file = output_file)

  # Check if the file is created
  expect_true(file.exists(output_file))

  # Clean up: remove the generated file after the test
  unlink(output_file)
})
