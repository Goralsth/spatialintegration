## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----load geomx data and create suerat object, eval=FALSE---------------------
# 
# #define file paths
# file_paths <- list(
#   "ACA_MOS_MOp_6mpi" = "/varidata/research/projects/henderson/THOMAS/Machine_Learning_For_Bioinformatics_Project/target_data_ACA_MOS_MOp_6mpi.RDS",
#   "ACA_MOS_MOp_9mpi" = "/varidata/research/projects/henderson/THOMAS/Machine_Learning_For_Bioinformatics_Project/target_data_ACA_MOS_MOp_9mpi.RDS",
#   "LAMDA" = "/varidata/research/projects/henderson/THOMAS/Machine_Learning_For_Bioinformatics_Project/target_data_LAMDA.RDS",
#   "MultiRegion" = "/varidata/research/projects/henderson/THOMAS/Machine_Learning_For_Bioinformatics_Project/target_data_MultiRegion.RDS",
#   "SNcAmygNew" = "/varidata/research/projects/henderson/THOMAS/Machine_Learning_For_Bioinformatics_Project/target_data_SNcAmygNew.RDS"
# )
# 
# 
# 
# #make suerat object
# seurat_object <- combine_datasets_into_seurat(file_paths)
# 
# 

## ----get nmf embeddings, eval=FALSE-------------------------------------------
# seurat_object <- process_and_run_nmf(seurat_object, nmf_rank = 5)
# RankPlot(seurat_object, detail.level = 2)
# 

## ----load reference Allen data, eval=FALSE------------------------------------
# 
# Allen_SEA_AD_Atlas <- readRDS("/varidata/research/projects/henderson/THOMAS/Machine_Learning_For_Bioinformatics_Project/Reference_MTG_RNAseq_all-nuclei.2022-06-07.rds")
# 

## ----run NMF on allen data, eval=FALSE----------------------------------------
# Allen_SEA_AD_Atlas <- process_and_run_nmf(Allen_SEA_AD_Atlas, nmf_rank = 5)
# 

## ----save objects, eval=FALSE-------------------------------------------------
# saveRDS(seurat_object,"seurat_object.rds")
# saveRDS(Allen_SEA_AD_Atlas,"Allen_SEA_AD_Atlas.rds")

## ----load in nmf embedded data, eval=FALSE------------------------------------
# 
#  seurat_object <- readRDS("/varidata/research/projects/henderson/THOMAS/SpatialIntegration/SpatialIntegration/seurat_object.rds")
# Allen_SEA_AD_Atlas <- readRDS("/varidata/research/projects/henderson/THOMAS/SpatialIntegration/SpatialIntegration/Allen_SEA_AD_Atlas.rds")
# 
# # Replace NA values in the "subclass_label" column of @meta.data with "NA"
# Allen_SEA_AD_Atlas@meta.data[["subclass_label"]][is.na(Allen_SEA_AD_Atlas@meta.data[["subclass_label"]])] <- "NotLabled"
# 

## ----subset seurat object for upload to data file in package, eval=FALSE------
# 
# seurat_object_sub<-subset_seurat_by_proportions(seurat_object = seurat_object, cell_type_col="segment", 500)
# allen_sub<-subset_seurat_by_proportions(seurat_object = Allen_SEA_AD_Atlas, cell_type_col="subclass_label", 6000)
# 

## ----eval=FALSE---------------------------------------------------------------
# 
#   seurat_merge<-merge_seurat_objects(seurat_object = seurat_object_sub, allen_sub)
# 
#   #save data to data file
#   use_data(seurat_merge, compress = "bzip2", overwrite = TRUE)
# 

## -----------------------------------------------------------------------------

library(devtools)


## ----set lib path and install package-----------------------------------------
#cahnge below libpath to directory with sufficient storage for large packages on an HPC environment
#.libPaths(c("/varidata/research/projects/henderson/THOMAS/R/4.4/library", .libPaths()))

install_github("Goralsth/spatialintegration", build_vignettes = TRUE)

library(SpatialIntegration)


## -----------------------------------------------------------------------------

#load the data
data(seurat_merge)


## -----------------------------------------------------------------------------
is(seurat_merge, "Seurat")

seurat_merge

## -----------------------------------------------------------------------------
#get indiceds of geomx data
geo_ind<-is.na(seurat_merge@meta.data[["specimen_name"]])

#get sum of geomx data
sum(seurat_merge@assays[["RNA"]]@layers[["counts"]][geo_ind])
#get sum of allen data
sum(seurat_merge@assays[["RNA"]]@layers[["counts"]][-geo_ind])


## -----------------------------------------------------------------------------

#for whatever reason, the singlet package will not load as a dependency, so load it here
library(singlet)

## -----------------------------------------------------------------------------

seurat_merge<-process_and_run_nmf(seurat_merge, nmf_rank = 15)


## -----------------------------------------------------------------------------


seurat_merge<-perform_gsea_and_clustering(seurat_merge)



## -----------------------------------------------------------------------------
#View umap clusters
print(DimPlot(seurat_merge, reduction = "umap"))
#View GeoMx segment lables
print(DimPlot(seurat_merge, reduction = "umap", group.by ="segment"))
#View Cell subclass labels from allen
print(DimPlot(seurat_merge, reduction = "umap", group.by ="subclass_label"))

#Veiw Cells labled NA in original Allen data
print(DimPlot(seurat_merge, reduction = "umap", group.by ="class_label"))

## -----------------------------------------------------------------------------
#View each cell type label and spatial annotations correlation to each nmf factor
plot_nmf_correlation_heatmap(seurat_merge, "cell_type_nmfFactor_corr.png")

## -----------------------------------------------------------------------------
#View nmf loadings on umaps
FeaturePlot(seurat_merge, features = paste0("NMF_", c(2,3,7,8,c(12:15))), raster = FALSE)


## -----------------------------------------------------------------------------

#get genesets that load with NeuN related factor
plot_top20_enriched_pathways_single(seurat_object = seurat_merge,
  nmf_factor = 13,
  output_file = "top20_enriched_nmf6.pdf"
)

#get genesets that load with pSyn related factor
plot_top20_enriched_pathways_single(seurat_object = seurat_merge,
  nmf_factor = 2,
  output_file = "top20_enriched_nmf15.pdf"
)



## -----------------------------------------------------------------------------

GSEAHeatmap(seurat_merge, reduction = "nmf", max.terms.per.factor = 50)
FeaturePlot(seurat_merge, features = paste0("NMF_", 1:10), raster = FALSE)
DimPlot(seurat_merge, reduction = "umap", group.by ="segment")

