#' Merge Seurat Object with Allen SEA AD Atlas
#'
#' This function merges a Seurat object with the Allen SEA AD Atlas after finding and subsetting
#' both datasets to their common features (row names).
#'
#' @param seurat_object A Seurat object to be merged.
#' @param allen_object A Seurat object representing the Allen SEA AD Atlas.
#' @return A merged Seurat object containing only the common features between the two datasets.
#' @examples
#' # Assuming seurat_object and Allen_SEA_AD_Atlas are pre-loaded Seurat objects
#' merged_seurat <- merge_allen_geomx(seurat_object, Allen_SEA_AD_Atlas)
#' @export
merge_allen_geomx <- function(seurat_object, allen_object) {
  # Convert rownames of the Seurat object to uppercase
  rownames(seurat_object) <- str_to_upper(rownames(seurat_object))

  # Create a list of datasets
  data_list <- list(seurat_object, allen_object)

  # Find common row names across the datasets
  common_row_names <- Reduce(intersect, lapply(data_list, rownames))

  # Subset both datasets by the common row names
  seurat_object <- subset(seurat_object, features = common_row_names)
  allen_object <- subset(allen_object, features = common_row_names)

  # Merge the two Seurat objects
  seurat_merge <- merge(seurat_object, allen_object)


  # Extract and bind count matrices by columns
  counts_1 <- GetAssayData(seurat_object, slot = "counts")
  counts_2 <- GetAssayData(allen_object, slot = "counts")
  merged_counts <- cbind(counts_1, counts_2) # column-binding coulumns for the merged count matrix


  seurat_merge@assays$RNA@layers$counts_merge<-merged_counts

  # Return the merged Seurat object
  return(seurat_merge)
}
