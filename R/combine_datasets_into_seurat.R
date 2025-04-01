#' Load and Combine Datasets into a Seurat Object
#'
#' This function loads multiple datasets, extracts and normalizes expression data,
#' aligns them based on common features and metadata, and combines them into a Seurat object.
#'
#' @param file_paths A named list of file paths to the RDS files. The names should correspond
#'   to dataset labels (e.g., `list("ACA_MOS_MOp_6mpi" = "path/to/file1.RDS", ...)`).
#' @return A Seurat object containing combined expression data and metadata.
#' @importFrom Seurat CreateSeuratObject
#' @examples
#' #' @examples
#' file_paths <- list.files("/path/to/", pattern = "target_data_*.RDS", full.names = TRUE)
#' if (length(file_paths) > 0) {
#'   seurat_obj <- combine_datasets_into_seurat(file_paths)
#' }

#' file_paths <- list(
#'   "ACA_MOS_MOp_6mpi" = "/path/to/target_data_ACA_MOS_MOp_6mpi.RDS",
#'   "ACA_MOS_MOp_9mpi" = "/path/to/target_data_ACA_MOS_MOp_9mpi.RDS",
#'   "LAMDA" = "/path/to/target_data_LAMDA.RDS",
#'   "MultiRegion" = "/path/to/target_data_MultiRegion.RDS",
#'   "SNcAmygNew" = "/path/to/target_data_SNcAmygNew.RDS"
#' )
#' seurat_object <- combine_datasets_into_seurat(file_paths)
#' @export
combine_datasets_into_seurat <- function(file_paths) {
  # Load datasets
  datasets <- lapply(file_paths, readRDS)

  # Extract q_norm data and find common row names
  q_norm_list <- lapply(datasets, function(data) data@assayData$q_norm)
  common_row_names <- Reduce(intersect, lapply(q_norm_list, rownames))

  # Subset q_norm data to common row names
  q_norm_list <- lapply(q_norm_list, function(q_norm) q_norm[common_row_names, , drop = FALSE])

  # Combine all q_norm datasets into one
  q_norm <- do.call(cbind, q_norm_list)

  # Remove the row named "NegProbe-WTX"
  q_norm <- q_norm[rownames(q_norm) != "NegProbe-WTX", ]

  # Extract phenoData and find common column names
  pheno_list <- lapply(datasets, function(data) data@phenoData@data)
  common_pheno_col_names <- Reduce(intersect, lapply(pheno_list, colnames))

  # Subset phenoData to common column names
  pheno_list <- lapply(pheno_list, function(pheno) pheno[, common_pheno_col_names, drop = FALSE])

  # Combine all phenoData into one
  pheno_combined <- do.call(rbind, pheno_list)

  # Ensure column names of q_norm match row names of pheno_combined
  if (!all(colnames(q_norm) %in% rownames(pheno_combined))) {
    stop("Column names of q_norm must match row names of pheno_combined.")
  }

  # Subset pheno_combined to match q_norm
  pheno_combined <- pheno_combined[colnames(q_norm), , drop = FALSE]

  # Create and return Seurat object
  seurat_object <- CreateSeuratObject(
    counts = q_norm,                 # Use q_norm as expression data
    meta.data = pheno_combined       # Add metadata
  )

  return(seurat_object)
}
