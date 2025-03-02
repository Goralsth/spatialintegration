#' Perform GSEA, Clustering, and Visualization with NMF on a Seurat Object
#'
#' This function runs GSEA using NMF factors, clusters the data, and generates visualizations
#' including heatmaps, UMAP, and feature plots. It also compares NMF clusters with PCA clusters.
#'
#' @param seurat_object A Seurat object with NMF reduction calculated.
#' @param gsea_category The category for GSEA (e.g., "C7"). Default is "C7".
#' @param umap_dims The number of dimensions to use for UMAP. Default is 1:30.
#' @param clustering_resolution The resolution for clustering. Default is 0.5.
#' @param max_terms_per_factor The maximum number of terms to display per factor in GSEA heatmap. Default is 3.
#' @return A list containing the updated Seurat object and a ggplot comparing NMF and PCA clusters.
#' @importFrom Seurat RunGSEA GSEAHeatmap FindNeighbors FindClusters RunUMAP FeaturePlot DimPlot
#' @importFrom dplyr group_by count mutate slice pull
#' @importFrom ggplot2 ggplot aes geom_point theme_bw labs scale_color_viridis_c
#' @examples
#' library(Seurat)
#' library(ggplot2)
#' seurat_object <- process_and_run_nmf(seurat_object, nmf_rank = 5)
#' results <- perform_gsea_and_clustering(seurat_object)
#' results$comparison_plot
#' @export
perform_gsea_and_clustering <- function(seurat_object,
                                        gsea_category = "C7",
                                        umap_dims = 1:30,
                                        clustering_resolution = 0.5,
                                        max_terms_per_factor = 3) {
  # Run GSEA
  seurat_object <- RunGSEA(seurat_object, category = gsea_category, verbose = FALSE)

  # Generate GSEA heatmap
  GSEAHeatmap(seurat_object, reduction = "nmf", max.terms.per.factor = max_terms_per_factor)

  # Find neighbors, clusters, and run UMAP
  seurat_object <- seurat_object %>%
    FindNeighbors(dims = umap_dims, reduction = "nmf") %>%
    FindClusters(resolution = clustering_resolution, verbose = FALSE) %>%
    RunUMAP(reduction = "nmf", dims = umap_dims, verbose = FALSE)

  # Plot NMF features
  FeaturePlot(seurat_object, features = paste0("NMF_", 1:6), raster = FALSE)

  # Compare NMF clusters with PCA clusters
  df <- data.frame(
    "nmf_clusters" = seurat_object@meta.data$seurat_clusters,
    "pca_clusters" = seurat_object@meta.data$cell_type
  )

  df <- df[!is.na(df$pca_clusters), ]
  df <- df %>%
    group_by(nmf_clusters) %>%
    count(pca_clusters) %>%
    mutate(freq = n / sum(n))

  comparison_plot <- ggplot(df, aes(nmf_clusters, pca_clusters, size = freq, color = n)) +
    geom_point() +
    theme_bw() +
    labs(
      x = "NMF cluster",
      y = "PCA cluster",
      size = "Proportion\nof cluster",
      color = "Cells in\nNMF cluster"
    ) +
    scale_color_viridis_c(option = "D")

  # Rename clusters based on PCA comparison
  cluster_names <- df %>%
    slice(which.max(n)) %>%
    pull(pca_clusters)

  levels(seurat_object@meta.data$seurat_clusters) <- make.unique(as.vector(cluster_names))

  # Generate UMAP plot with updated cluster names
  DimPlot(
    seurat_object,
    reduction = "umap",
    label = TRUE,
    group
