#' Subset Seurat Object While Maintaining Proportions of Cell Types
#'
#' This function randomly subsets a Seurat object, preserving the proportions of specified cell types.
#' The cell type information is obtained from a column in the metadata. A report of the original
#' and downsampled counts for each cell type is returned. A progress bar is displayed during the sampling process.
#'
#' @param seurat_object A Seurat object to subset.
#' @param cell_type_col Character. The name of the metadata column containing cell type annotations.
#' @param subset_size Numeric. The desired total number of cells in the subset.
#' @return A list containing:
#'   \item{seurat_subset}{A Seurat object containing the subset of cells with preserved cell type proportions.}
#'   \item{report}{A data frame with the original and downsampled counts for each cell type.}
#'
#' @examples
#' # Example usage:
#' result <- subset_seurat_by_proportions(
#'   seurat_object = seurat_obj,
#'   cell_type_col = "cell_type",
#'   subset_size = 1000
#' )
#' seurat_subset <- result$seurat_subset
#' report <- result$report
#'
#' @import dplyr Seurat progress
#' @export
subset_seurat_by_proportions <- function(seurat_object, cell_type_col, subset_size) {
  # Check if the specified column exists in the metadata
  if (!cell_type_col %in% colnames(seurat_object@meta.data)) {
    stop("The specified column '", cell_type_col, "' does not exist in the metadata.")
  }

  # Extract metadata and calculate original counts and proportions for each cell type
  meta_data <- seurat_object@meta.data
  proportions <- meta_data %>%
    group_by(.data[[cell_type_col]]) %>%
    summarise(
      original_count = n(),
      .groups = "drop"
    ) %>%
    mutate(proportion = original_count / sum(original_count))

  # Determine the number of cells to sample for each cell type
  proportions <- proportions %>%
    mutate(downsampled_count = round(proportion * subset_size))

  # Initialize an empty vector to hold selected cell IDs
  selected_cells <- c()

  # Initialize the progress bar
  pb <- progress::progress_bar$new(
    total = nrow(proportions),
    format = "  Downsampling [:bar] :percent :current/:total"
  )

  # For each cell type, randomly sample the appropriate number of cell IDs
  for (i in seq_len(nrow(proportions))) {
    pb$tick() # Update progress bar

    cell_type <- proportions[[cell_type_col]][i]
    downsampled_count <- proportions$downsampled_count[i]

    # Get cell IDs for the current cell type
    cell_ids <- rownames(meta_data[meta_data[[cell_type_col]] == cell_type, ])

    # Randomly sample cell IDs
    selected_cells <- c(selected_cells, sample(cell_ids, size = downsampled_count))
  }

  # Subset the Seurat object using the selected cell IDs
  seurat_subset <- subset(seurat_object, cells = selected_cells)

  # Create the report
  report <- proportions %>%
    select(cell_type = .data[[cell_type_col]], original_count, downsampled_count)

  return(list(
    seurat_subset = seurat_subset,
    report = report
  ))
}

