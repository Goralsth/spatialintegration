#' Identify and Plot Top 20 Enriched Pathways for Correlated NMF Factor
#'
#' This function identifies the top 20 gene sets enriched in the specified NMF factor
#' for a given dataset, weighted by NES and p-values, and generates a dot plot
#' of the top 20 enriched pathways.
#'
#' @param seurat_object A Seurat object containing GSEA results in the NMF reduction.
#' @param nmf_factor The NMF factor to analyze for enrichment (default is 9).
#' @param output_file The path to save the output dot plot as a PDF (default is "top20_dotplot.pdf").
#' @return A ggplot object of the dot plot for the top 20 enriched pathways.
#' @importFrom ggplot2 ggplot aes geom_point scale_size_continuous scale_color_gradient coord_flip labs theme_minimal theme element_text ggsave
#' @examples
#' library(Seurat)
#' library(ggplot2)
#' p <- plot_top20_enriched_pathways_dotplot(
#'   seurat_object = seurat_object,
#'   nmf_factor = 9,
#'   output_file = "top20_dotplot.pdf"
#' )
#' print(p)
#' @export
plot_top20_enriched_pathways_single <- function(seurat_object,
                                                nmf_factor = 9,
                                                output_file = "top20_dotplot.pdf") {
  # Step 1: Extract NES values and p-values for the specified NMF factor
  nes_values <- seurat_object@reductions[["nmf"]]@misc[["gsea"]][["nes"]][, nmf_factor]
  p_values <- seurat_object@reductions[["nmf"]]@misc[["gsea"]][["pval"]][, nmf_factor]

  # Check if NES and p-values are correctly extracted
  if (is.null(nes_values) || is.null(p_values)) {
    stop("NES or p-values could not be extracted from the Seurat object.")
  }

  nes_values <- as.numeric(nes_values)  # Ensure that NES values are numeric
  p_values <- as.numeric(p_values)     # Ensure that p-values are numeric

  # Handle zero p-values by replacing them with a small value
  p_values[p_values == 0] <- 1e-10

  # Step 2: Get pathway names
  pathways <- rownames(seurat_object@reductions[["nmf"]]@misc[["gsea"]][["nes"]])

  # Step 3: Create a data frame with a combined score (e.g., NES weighted by p-value)
  # Using -log10(p-value) to emphasize significant pathways
  pathway_df <- data.frame(
    pathway = pathways,
    nes = nes_values,
    pval = p_values,
    combined_score = nes_values * -log10(p_values)
  )

  # Remove rows with missing or invalid values
  pathway_df <- na.omit(pathway_df)
  pathway_df <- pathway_df[is.finite(pathway_df$combined_score), ]

  # Check if there are valid rows left
  if (nrow(pathway_df) == 0) {
    stop("No valid data available after filtering missing or invalid values.")
  }

  # Step 4: Sort by combined score and get the top 20 pathways
  top20_pathways <- pathway_df[order(pathway_df$combined_score, decreasing = TRUE), ][1:20, ]

  # Step 5: Create the ggplot dot plot
  p <- ggplot(top20_pathways, aes(y = reorder(pathway, combined_score), x = nes)) +
    geom_point(aes(size = -log10(pval), color = combined_score)) +
    scale_size_continuous(
      name = "-log10(p-value)",  # Legend title for size
      range = c(2, 10)          # Adjust point size range as desired
    ) +
    scale_color_gradient(
      low = "blue",
      high = "red",
      name = "Combined Score"
    ) +
    labs(
      x = "Normalized Enrichment Score (NES)",  # Swap x and y labels if necessary
      y = "Pathway",
      title = paste("Top 20 Enriched Pathways in NMF (Factor", nmf_factor, ")")
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 10),  # Customize axis text
      axis.text.y = element_text(size = 10),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 10)
    )


  # Step 6: Save the plot as a PDF
  ggsave(output_file, p, width = 10, height = 8, units = "cm")

  return(p)
}

