#' Perform GSEA and Clustering on Seurat Object
#'
#' This function performs Gene Set Enrichment Analysis (GSEA), dimensionality reduction,
#' clustering, and visualization on a Seurat object using specified parameters.
#'
#' @param seurat_object A Seurat object to perform analysis on.
#' @param gsea_category Character. The GSEA category to use for enrichment analysis.
#' Default is "C7" (immunologic signatures).
#' @param umap_dims Numeric vector. The dimensions to use for UMAP and clustering. Default is 1:30.
#' @param clustering_resolution Numeric. The resolution parameter for clustering. Default is 0.5.
#' @param max_terms_per_factor Numeric. The maximum number of terms to display per factor in the GSEA heatmap. Default is 3.
#'
#' @return A Seurat object with updated metadata, clustering, and visualization results.
#'
#' @examples
#' # Example usage:
#' seurat_object <- perform_gsea_and_clustering(
#'   seurat_object = seurat_obj,
#'   gsea_category = "C7",
#'   umap_dims = 1:20,
#'   clustering_resolution = 0.6,
#'   max_terms_per_factor = 5
#' )
#'
#' @import Seurat ggplot2 viridis dplyr
#' @export
perform_gsea_and_clustering <- function(seurat_object,
                                        gsea_category = "C7",
                                        umap_dims = 1:30,
                                        clustering_resolution = 0.5,
                                        max_terms_per_factor = 3) {
  # Run GSEA
  seurat_object <- RunGSEA(seurat_object, category = gsea_category, verbose = FALSE)

  # Generate GSEA heatmap
  p<-GSEAHeatmap(seurat_object, reduction = "nmf", max.terms.per.factor = max_terms_per_factor)
print(p)

# Adjust umap_dims based on available dimensions in "nmf"
available_dims <- ncol(seurat_object@reductions$nmf@cell.embeddings)
umap_dims <- umap_dims[umap_dims <= available_dims]
if (length(umap_dims) == 0) {
  stop("No valid dimensions available for UMAP.")
}

# Find neighbors, clusters, and run UMAP
seurat_object <- seurat_object %>%
  FindNeighbors(dims = umap_dims, reduction = "nmf") %>%
  FindClusters(resolution = clustering_resolution, verbose = FALSE) %>%
  RunUMAP(reduction = "nmf", dims = umap_dims, verbose = FALSE)

  # Plot NMF features
  p<-FeaturePlot(seurat_object, features = paste0("NMF_", 1:6), raster = FALSE)
print(p)
  # Compare NMF clusters with PCA clusters if cell_type metadata exists
  if ("cell_type" %in% colnames(seurat_object@meta.data)) {
    df <- data.frame(
      "nmf_clusters" = seurat_object@meta.data$seurat_clusters,
      "pca_clusters" = seurat_object@meta.data$cell_type
    )

    df <- df[!is.na(df$pca_clusters), ] %>%
      group_by(nmf_clusters) %>%
      count(pca_clusters) %>%
      mutate(freq = n / sum(n))

    comparison_plot <- ggplot(df, aes(x = nmf_clusters, y = pca_clusters, size = freq, color = n)) +
      geom_point() +
      theme_bw() +
      labs(
        x = "NMF cluster",
        y = "PCA cluster",
        size = "Proportion\nof cluster",
        color = "Cells in\nNMF cluster"
      ) +
      scale_color_viridis_c(option = "D")
    print(comparison_plot)

    # Rename clusters based on PCA comparison
    cluster_names <- df %>%
      slice(which.max(n)) %>%
      pull(pca_clusters)

    levels(seurat_object@meta.data$seurat_clusters) <- make.unique(as.vector(cluster_names))

    # Generate UMAP plot with updated cluster names
    DimPlot(seurat_object, reduction = "umap", label = TRUE, group.by = "seurat_clusters")
  } else {
    message("No 'cell_type' column found in metadata for PCA comparison.")
  }

  return(seurat_object)
}

