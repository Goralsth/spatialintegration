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
#' @import Seurat SeuratObject
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
  #seurat_object <- NormalizeData(seurat_object)
  seurat_object@assays$RNA$data <- seurat_object@assays$RNA$counts
  seurat_object <- ScaleData(seurat_object)
  seurat_object <- FindVariableFeatures(seurat_object)
  seurat_object <- Seurat::RunPCA(seurat_object)

  # --- SCALE TO EQUAL CONTRIBUTION OF TRANSCRIPTS BETWEEN DATASETS ---
  geo_ind <- is.na(seurat_object@meta.data[["specimen_name"]])
  geo_sum <- sum(seurat_object@assays[["RNA"]]@layers[["counts"]][, geo_ind])
  allen_sum <- sum(seurat_object@assays[["RNA"]]@layers[["counts"]][, !geo_ind])
  scale_factor_allen <- geo_sum / allen_sum

  counts <- seurat_object@assays[["RNA"]]@layers[["counts"]]
  counts_scaled <- counts
  counts_scaled[, !geo_ind] <- counts_scaled[, !geo_ind] * scale_factor_allen
  seurat_object@assays[["RNA"]]@layers[["scaled_counts"]] <- counts_scaled
  # ---------------------------------------------------------------------

  # Assign data slot to "layers"
  seurat_object@assays[["RNA"]]@layers[["data"]] <- seurat_object@assays[["RNA"]]@layers[["scaled_counts"]]

  RunNMF.Seurat <- function(object,
                            split.by = NULL,
                            k = NULL,
                            assay = NULL,
                            reps = 3,
                            tol = 1e-5,
                            L1 = 0.01,
                            L2 = 0,
                            verbose = 2,
                            reduction.name = "nmf",
                            reduction.key = "NMF_",
                            maxit = 100,
                            test.set.density = 0.05,
                            learning.rate = 0.8,
                            tol.overfit = 1e-4,
                            trace.test.mse = 5,
                            threads = 0,
                            features = NULL,
                            ...) {
    if (is.null(assay)) {
      assay <- DefaultAssay(object)
    }

    # Check if data has been normalized
    v <- GetAssayData(object, assay = assay, slot = "data")@x
    if (sum(as.integer(v)) == sum(v)) {
      object <- PreprocessData(object, assay = assay)
    }
    A <- GetAssayData(object, assay = assay, slot = "data")

    if (!is.null(features)) {
      if (features[[1]] == "var.features") {
        A <- A[VariableFeatures(object, assay = assay), ]
      } else if (is.integer(features) || is.character(features)) {
        A <- A[features, ]
      } else {
        stop("'features' vector was invalid.")
      }
    }

    rnames <- rownames(A)
    cnames <- colnames(A)

    if (!is.null(split.by)) {
      split.by <- as.integer(as.numeric(as.factor(object[[split.by]][, 1]))) - 1
      if (any(sapply(split.by, is.na))) {
        stop("'split.by' cannot contain NA values")
      }
      A@x <- A@x + 1
      A@x <- A@x - 1
      A <- weight_by_split(A, split.by, length(unique(split.by)))
    }

    At <- Matrix::t(A)
    seed.use <- abs(.Random.seed[[3]])

    if (!is.null(k) && length(k) > 1) {
      cv_data <- cross_validate_nmf(
        A = A,
        ranks = k,
        n_replicates = reps,
        tol = tol * 10,
        maxit = maxit,
        verbose = verbose,
        L1 = L1,
        L2 = L2,
        threads = threads,
        test_density = test.set.density,
        tol_overfit = tol.overfit,
        trace_test_mse = trace.test.mse
      )
      best_rank <- GetBestRank(cv_data, tol.overfit)
      if (verbose >= 1) {
        cat("best rank: ", best_rank, "\n")
      }
      cat("\nfitting final model:\n")
      nmf_model <- run_nmf(A, best_rank, tol, maxit, verbose > 1, L1, L2, threads)
    } else if (is.null(k)) {
      nmf_model <- ard_nmf(
        A = A,
        k_init = k,
        k_max = 1e4,
        k_min = 2,
        n_replicates = reps,
        tol = tol,
        maxit = maxit,
        verbose = verbose,
        L1 = L1,
        L2 = L2,
        threads = threads,
        test_density = test.set.density,
        learning_rate = learning.rate,
        tol_overfit = tol.overfit,
        trace_test_mse = trace.test.mse
      )
      cv_data <- nmf_model$cv_data
    } else if (length(k) == 1) {
      nmf_model <- run_nmf(A, k, tol, maxit, verbose > 1, L1, L2, threads)
      cv_data <- list()
    } else {
      stop("value for 'k' was invalid")
    }

    rownames(nmf_model$h) <- colnames(nmf_model$w) <- paste0(reduction.key, 1:nrow(nmf_model$h))
    rownames(nmf_model$w) <- rnames
    colnames(nmf_model$h) <- cnames

    object[[reduction.name]] <- CreateDimReducObject(
      embeddings = t(nmf_model$h),
      loadings = nmf_model$w,
      assay = assay,
      key = reduction.key,
      misc = list("cv_data" = cv_data)
    )

    return(object)
  }

  set.seed(seed)
  seurat_object <- RunNMF.Seurat(
    seurat_object,
    k = nmf_rank,
    reduction.name = nmf_reduction_name,
    reduction.key = nmf_reduction_key
  )

  return(seurat_object)
}
