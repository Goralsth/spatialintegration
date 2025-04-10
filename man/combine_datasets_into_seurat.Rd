% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/combine_datasets_into_seurat.R
\name{combine_datasets_into_seurat}
\alias{combine_datasets_into_seurat}
\title{Load and Combine Datasets into a Seurat Object}
\usage{
combine_datasets_into_seurat(file_paths)
}
\arguments{
\item{file_paths}{A named list of file paths to the RDS files. The names should correspond
to dataset labels (e.g., `list("ACA_MOS_MOp_6mpi" = "path/to/file1.RDS", ...)`).}
}
\value{
A Seurat object containing combined expression data and metadata.
}
\description{
This function loads multiple datasets, extracts and normalizes expression data,
aligns them based on common features and metadata, and combines them into a Seurat object.
}
\examples{
#' @examples
file_paths <- list.files("/path/to/", pattern = "target_data_*.RDS", full.names = TRUE)
if (length(file_paths) > 0) {
  seurat_obj <- combine_datasets_into_seurat(file_paths)
}
file_paths <- list(
  "ACA_MOS_MOp_6mpi" = "/path/to/target_data_ACA_MOS_MOp_6mpi.RDS",
  "ACA_MOS_MOp_9mpi" = "/path/to/target_data_ACA_MOS_MOp_9mpi.RDS",
  "LAMDA" = "/path/to/target_data_LAMDA.RDS",
  "MultiRegion" = "/path/to/target_data_MultiRegion.RDS",
  "SNcAmygNew" = "/path/to/target_data_SNcAmygNew.RDS"
)
seurat_object <- combine_datasets_into_seurat(file_paths)
}
