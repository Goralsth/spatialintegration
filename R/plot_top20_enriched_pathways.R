#' Identify and Plot Top 20 Enriched Pathways for Correlated NMF Factors
#'
#' This function identifies the top 20 gene sets enriched in the most correlated NMF factors
#' for spatial and single-cell datasets and generates a comparative bar chart.
#'
#' @param seurat_object A Seurat object containing GSEA results in the NMF reduction for the Seurat dataset.
#' @param allen_object A Seurat object containing GSEA results in the NMF reduction for the Allen dataset.
#' @param nmf_factor The NMF factor to analyze for enrichment (default is 9).
#' @param output_file The path to save the output bar chart as a PDF (default is "top20.pdf").
#' @return A ggplot object of the bar chart for the top 20 enriched pathways.
#' @importFrom ggplot2 ggplot aes geom_bar coord_flip labs theme_minimal theme element_text ggsave
#' @examples
#' library(Seurat)
#' library(ggplot2)
#' p <- plot_top20_enriched_pathways(
#'   seurat_object = seurat_object,
#'   allen_object = allen_object,
#'   nmf_factor = 9,
#'   output_file = "top20_enriched.pdf"
#' )
#' print(p)
#' @export
plot_top20_enriched_pathways <- function(seurat_object,
                                         allen_object,
                                         nmf_factor = 9,
                                         output_file = "top20.pdf") {
  # Step 1: Extract NES values for the specified NMF factor
  seurat_nes <- seurat_object@reductions[["nmf"]]@misc[["gsea"]][["nes"]][, nmf_factor]
  allen_nes <- allen_object@reductions[["nmf"]]@misc[["gsea"]][["nes"]][, nmf_factor]

  # Step 2: Get pathway names
  seurat_pathways <- rownames(seurat_object@reductions[["nmf"]]@misc[["gsea"]][["nes"]])
  allen_pathways <- rownames(allen_object@reductions[["nmf"]]@misc[["gsea"]][["nes"]])

  # Step 3: Create data frames for both datasets
  seurat_df <- data.frame(pathway = seurat_pathways, nes = seurat_nes)
  allen_df <- data.frame(pathway = allen_pathways, nes = allen_nes)

  # Step 4: Sort by NES and get the top 20 pathways
  seurat_top20 <- seurat_df[order(seurat_df$nes, decreasing = TRUE), ][1:20, ]
  allen_top20 <- allen_df[order(allen_df$nes, decreasing = TRUE), ][1:20, ]

  # Add dataset labels
  seurat_top20$dataset <- "Seurat"
  allen_top20$dataset <- "Allen"

  # Combine both datasets into one data frame
  combined_df <- rbind(seurat_top20, allen_top20)

  # Step 5: Create the ggplot bar chart
  p <- ggplot(combined_df, aes(x = reorder(pathway, nes), y = nes, fill = dataset)) +
    geom_bar(stat = "identity", position = "dodge") +
    coord_flip() +
    labs(x = "Pathway", y = "Normalized Enrichment Score (NES)",
         title = "Top 20 Enriched Pathways in NMF (Seurat vs Allen)") +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 5),
      legend.key.size = unit(0.00001, 'cm'),
      legend.title = element_text(size = 2),
      legend.text = element_text(size = 2),
      legend.key.height = unit(0.3, 'cm'),
      legend.key.width = unit(3, 'cm')
    )

  # Step 6: Save the plot as a PDF
  ggsave(output_file, p, width = 40, height = 30, units = "cm")

  return(p)
}
