#' Correlation Heatmap of NMF Factors and Cell Types
#'
#' This function computes the Pearson correlation between NMF factor loadings
#' and one-hot encoded cell type annotations (subclass labels and segments) and
#' generates a heatmap of the correlations.
#'
#' @param seurat_object A Seurat object containing NMF factor loadings and metadata.
#' @param output_file Path to save the heatmap plot as a file (default: "cell_type_nmfFactor_corr.png").
#' @return A ggplot object representing the heatmap of correlations.
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradient2 theme_minimal element_text labs ggsave
#' @examples
#' library(Seurat)
#' library(ggplot2)
#' p <- plot_nmf_correlation_heatmap(seurat_object = seurat_merge, output_file = "heatmap.png")
#' print(p)
#' @export
plot_nmf_correlation_heatmap <- function(seurat_object, output_file = "cell_type_nmfFactor_corr.png") {
  # Extract NMF factor loadings
  nmf_factors <- seurat_object@reductions$nmf@cell.embeddings

  # Align rownames between NMF factors and metadata
  common_cells <- intersect(rownames(nmf_factors), rownames(seurat_object@meta.data))
  nmf_factors <- nmf_factors[common_cells, ]

  # Extract cell type annotations
  cell_types <- data.frame(
    subclass_label = seurat_object@meta.data$subclass_label,
    segment = seurat_object@meta.data$segment,
    row.names = rownames(seurat_object@meta.data)
  )

  # Handle missing values
  cell_types$subclass_label[is.na(cell_types$subclass_label)] <- "NA_subclass"
  cell_types$segment[is.na(cell_types$segment)] <- "NA_segment"

  # Convert to factors
  cell_types$subclass_label <- factor(cell_types$subclass_label)
  cell_types$segment <- factor(cell_types$segment)

  # One-hot encode the 'subclass_label' and 'segment' factors
  subclass_one_hot <- model.matrix(~ 0 + subclass_label, data = cell_types)
  segment_one_hot <- model.matrix(~ 0 + segment, data = cell_types)

  # Combine the one-hot encoded matrices
  one_hot_encoded <- cbind(subclass_one_hot, segment_one_hot)

  # Compute correlations for each NMF factor
  cor_results <- apply(nmf_factors, 2, function(factor) {
    apply(one_hot_encoded, 2, function(cell_vector) {
      cor(factor, cell_vector, method = "pearson")
    })
  })

  # Convert the correlation results into a data frame
  cor_df <- as.data.frame(cor_results)
  cor_df$NMF_factor <- rownames(cor_df)

  # Reshape the data for heatmap visualization
  cor_df_melted <- reshape2::melt(cor_df, id.vars = "NMF_factor", variable.name = "Cell_Type", value.name = "Correlation")

  # Create a heatmap using ggplot2
  p <- ggplot(cor_df_melted, aes(x = Cell_Type, y = NMF_factor, fill = Correlation)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    theme_minimal() +
    labs(x = "Cell Type", y = "NMF Factor", fill = "Correlation") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

  # Save the plot
  ggsave(output_file, p)

  return(p)
}

