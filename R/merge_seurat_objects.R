#' Merge Two Seurat Objects
#'
#' This function merges two Seurat objects by combining their count matrices and metadata.
#' It ensures compatibility between the two datasets and creates a new Seurat object.
#'
#' @param seurat_object A Seurat object to be merged.
#' @param Allen_SEA_AD_Atlas Another Seurat object to merge with `seurat_object`.
#' @return A new Seurat object containing the merged counts and metadata.
#' @importFrom Seurat CreateSeuratObject
#' @examples
#' # Example usage:
#' seurat_merge <- merge_seurat_objects(seurat_object, Allen_SEA_AD_Atlas)
#' @export
merge_seurat_objects <- function(seurat_object, Allen_SEA_AD_Atlas) {

  rownames(seurat_object)<-str_to_upper(rownames(seurat_object))
  data_list <- list(seurat_object,Allen_SEA_AD_Atlas)
  common_row_names <- Reduce(intersect, lapply(data_list, rownames))

  seurat_object <- subset(seurat_object, features = common_row_names)
  Allen_SEA_AD_Atlas <- subset(Allen_SEA_AD_Atlas, features = common_row_names)

  # Extract count matrices
  express_geo <- seurat_object@assays$RNA@layers$counts
  express_allen <- Allen_SEA_AD_Atlas@assays$RNA@counts
  express_merge <- cbind(express_geo, express_allen)

  # Extract metadata
  meta_geo <- seurat_object@meta.data
  meta_allen <- Allen_SEA_AD_Atlas@meta.data

  # Ensure rownames in meta_geo match column names of express_geo
  rownames(meta_geo) <- colnames(express_geo)

  # Ensure rownames in meta_allen match column names of express_allen
  rownames(meta_allen) <- colnames(express_allen)

  # Get the union of column names
  all_columns <- union(colnames(meta_geo), colnames(meta_allen))

  # Add missing columns with NA values to meta_geo
  missing_columns_geo <- setdiff(all_columns, colnames(meta_geo))
  meta_geo[missing_columns_geo] <- NA

  # Add missing columns with NA values to meta_allen
  missing_columns_allen <- setdiff(all_columns, colnames(meta_allen))
  meta_allen[missing_columns_allen] <- NA

  # Align columns
  meta_geo <- meta_geo[, all_columns, drop = FALSE]
  meta_allen <- meta_allen[, all_columns, drop = FALSE]

  # Combine metadata
  meta_merge <- rbind(meta_geo, meta_allen)

  # Ensure rownames are unique and non-empty
  rownames(meta_merge) <- make.unique(rownames(meta_merge), sep = "_")
  rownames(meta_merge) <- ifelse(rownames(meta_merge) == "", paste0("sample_", seq_len(nrow(meta_merge))), rownames(meta_merge))

  # Ensure colnames of counts matrix match rownames of metadata
  colnames(express_merge) <- rownames(meta_merge)

  # Create Seurat object
  seurat_merge <- CreateSeuratObject(
    counts = express_merge,
    meta.data = meta_merge
  )

  # Return the merged Seurat object
  return(seurat_merge)
}
