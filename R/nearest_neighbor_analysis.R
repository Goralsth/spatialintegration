#' Nearest Neighbor Matching and Analysis for Seurat and Allen Data
#'
#' This function finds the nearest neighbors in the Allen data for each Seurat cell based on NMF factors,
#' computes correlations between matched NMF factors, and extracts NES values and pathways for NMF factors
#' from both datasets.
#'
#' @param nmf_seurat A matrix of NMF factors from the Seurat dataset.
#' @param nmf_allen A matrix of NMF factors from the Allen dataset.
#' @param seurat_object A Seurat object containing GSEA results in the NMF reduction.
#' @param allen_object A Seurat object containing GSEA results in the NMF reduction for the Allen dataset.
#' @param nmf_factor The NMF factor to analyze NES values for. Default is 9.
#' @return A list containing the correlation matrix, pathway data frames for Seurat and Allen, and the matched Allen NMF matrix.
#' @importFrom FNN get.knnx
#' @importFrom stats cor
#' @importFrom pheatmap pheatmap
#' @examples
#' library(Seurat)
#' library(FNN)
#' results <- nearest_neighbor_analysis(
#'   nmf_seurat = nmf_seurat_matrix,
#'   nmf_allen = nmf_allen_matrix,
#'   seurat_object = seurat_object,
#'   allen_object = allen_object
#' )
#' results$correlation_matrix
#' @export
nearest_neighbor_analysis <- function(nmf_seurat,
                                      nmf_allen,
                                      seurat_object,
                                      allen_object,
                                      nmf_factor = 9) {
  # Load necessary library
  library(FNN)

  # Find nearest neighbors in Allen for each Seurat cell
  nn_search <- get.knnx(nmf_allen, nmf_seurat, k = 1)  # k=1 gets the closest match

  # Extract matched cells
  nmf_allen_matched <- nmf_allen[nn_search$nn.index, , drop = FALSE]

  # Compute correlation
  cor_matrix <- cor(nmf_seurat, nmf_allen_matched, method = "pearson")
  pheatmap(cor_matrix, cluster_rows = TRUE, cluster_cols = TRUE,
           main = "Correlation of NMF Factors (Seurat vs. Nearest Allen Cells)")

  # Extract NES values for the specified NMF factor
  seurat_nmf_nes <- seurat_object@reductions[["nmf"]]@misc[["gsea"]][["nes"]][, nmf_factor]
  allen_nmf_nes <- allen_object@reductions[["nmf"]]@misc[["gsea"]][["nes"]][, nmf_factor]

  # Get pathway names
  seurat_pathways <- rownames(seurat_object@reductions[["nmf"]]@misc[["gsea"]][["nes"]])
  allen_pathways <- rownames(allen_object@reductions[["nmf"]]@misc[["gsea"]][["nes"]])

  # Create data frames for both Seurat and Allen NMF pathways and NES
  seurat_df <- data.frame(pathway = seurat_pathways, nes = seurat_nmf_nes)
  allen_df <- data.frame(pathway = allen_pathways, nes = allen_nmf_nes)

  return(list(
    correlation_matrix = cor_matrix,
    seurat_pathway_data = seurat_df,
    allen_pathway_data = allen_df,
    nmf_allen_matched = nmf_allen_matched
  ))
}

