% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/perform_gsea_and_clustering.R
\name{perform_gsea_and_clustering}
\alias{perform_gsea_and_clustering}
\title{Perform GSEA and Clustering on Seurat Object}
\usage{
perform_gsea_and_clustering(
  seurat_object,
  gsea_category = "C7",
  umap_dims = 1:30,
  clustering_resolution = 0.5,
  max_terms_per_factor = 3
)
}
\arguments{
\item{seurat_object}{A Seurat object to perform analysis on.}

\item{gsea_category}{Character. The GSEA category to use for enrichment analysis.
Default is "C7" (immunologic signatures).}

\item{umap_dims}{Numeric vector. The dimensions to use for UMAP and clustering. Default is 1:30.}

\item{clustering_resolution}{Numeric. The resolution parameter for clustering. Default is 0.5.}

\item{max_terms_per_factor}{Numeric. The maximum number of terms to display per factor in the GSEA heatmap. Default is 3.}
}
\value{
A Seurat object with updated metadata, clustering, and visualization results.
}
\description{
This function performs Gene Set Enrichment Analysis (GSEA), dimensionality reduction,
clustering, and visualization on a Seurat object using specified parameters.
}
\examples{
# Example usage:
seurat_object <- perform_gsea_and_clustering(
  seurat_object = seurat_obj,
  gsea_category = "C7",
  umap_dims = 1:20,
  clustering_resolution = 0.6,
  max_terms_per_factor = 5
)

}
