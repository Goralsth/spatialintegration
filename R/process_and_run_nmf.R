#' Perform Preprocessing, PCA, and NMF on a Seurat Object
#'
#' This function performs data normalization, scaling, PCA, and NMF on a Seurat object.
#' It assigns the results to the appropriate slots and provides tools for visualization.
#'
#' @param seurat_object A Seurat object to be processed.
#' @param nmf_rank The number of factors (k) for NMF. If NULL, automatic rank determination is performed.
#' @param nmf_reduction_name The name of the reduction slot to assign the NMF results. Default is "nmf".
#' @param nmf_reduction_key The prefix for NMF dimensions in the reduction slot. Default is "NMF_".
#' @param seed The random seed for reproducibility. Default is 16.
#' @return A processed Seurat object with PCA and NMF results.
#' @importFrom Seurat NormalizeData ScaleData FindVariableFeatures RunPCA CreateDimReducObject
#' @importFrom SeuratObject DefaultAssay GetAssayData VariableFeatures
#' @examples
#' library(Seurat)
#' seurat_object <- CreateSeuratObject(counts = matrix(rnorm(1000), nrow = 10))
#' seurat_object <- process_and_run_nmf(seurat_object, nmf_rank = 5)
#' RankPlot(seurat_object, detail.level = 2)
#' @export
process_and_run_nmf <- function(seurat_object,
                                nmf_rank = NULL,
                                nmf_reduction_name = "nmf",
                                nmf_reduction_key = "NMF_",
                                seed = 16) {
  # Normalize, scale, and perform PCA
  seurat_object <- NormalizeData(seurat_object)
  seurat_object <- ScaleData(seurat_object)
  seurat_object <- FindVariableFeatures(seurat_object)
  seurat_object <- Seurat::RunPCA(seurat_object)

  # Assign data slot to "layers"
  seurat_object@assays[["RNA"]]@layers[["data"]] <- seurat_object@assays[["RNA"]]@layers[["counts"]]

  # Run NMF
  RunNMF.Seurat <- function(object,
                            k = NULL,
                            reduction.name = "nmf",
                            reduction.key = "NMF_",
                            ...) {
    # Default assay and data validation
    assay <- DefaultAssay(object)
    v <- GetAssayData(object, assay = assay, slot = "data")@x
    if (sum(as.integer(v)) == sum(v)) {
      object <- PreprocessData(object, assay = assay)
    }
    A <- GetAssayData(object, assay = assay, slot = "data")
    rnames <- rownames(A)
    cnames <- colnames(A)

    # Run NMF with a single rank or automatic determination
    if (is.null(k)) {
      nmf_model <- ard_nmf(A, k_max = 10, tol = 1e-5, maxit = 100, verbose = 2)
    } else {
      nmf_model <- run_nmf(A, k, tol = 1e-5, maxit = 100, verbose = 2)
    }

    rownames(nmf_model$h) <- colnames(nmf_model$w) <- paste0(reduction.key, 1:nrow(nmf_model$h))
    rownames(nmf_model$w) <- rnames
    colnames(nmf_model$h) <- cnames

    # Assign results to reductions
    object[[reduction.name]] <- CreateDimReducObject(
      embeddings = t(nmf_model$h),
      loadings = nmf_model$w,
      assay = assay,
      key = reduction.key,
      misc = list("cv_data" = nmf_model$cv_data)
    )
    return(object)
  }

  # Set seed and run NMF
  set.seed(seed)
  seurat_object <- RunNMF.Seurat(
    seurat_object,
    k = nmf_rank,
    reduction.name = nmf_reduction_name,
    reduction.key = nmf_reduction_key
  )

  # Return processed Seurat object
  return(seurat_object)
}
