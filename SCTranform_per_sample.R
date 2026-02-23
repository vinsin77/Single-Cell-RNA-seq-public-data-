#SCTransform step (individual run for per sample)
#3 normals and 8 patients: 11 samples in total

#first refreshing the memory:
rm(list = ls())
gc()

library(Seurat)
library(dplyr)
getwd() ##must be the work directory
data_dir <- "/mnt/bctl/public/chen_PMID33033240/cellranger"

#loading all QC results
sample_list <- readRDS("rds_objects/all_samples_QC.rds")

#summary
length(sample_list)
sapply(sample_list, ncol)

#using 8 cores:
#library(future)
#plan(sequential)         
#plan("multicore", workers = 8)

#SCTransform-normalisation for each sample (e.g. normal2)
sample_list[["normal2"]] <- SCTransform(sample_list[["normal2"]], vars.to.regress = "percent.mt", verbose = TRUE)
#and save it each time
saveRDS(sample_list, "rds_objects/all_samples_SCT_progress.rds")

#after completing the normalisation, restart R environment
library(Seurat)
sample_list <- readRDS("rds_objects/all_samples_SCT_progress.rds")
length(sample_list)

features <- SelectIntegrationFeatures(
  object.list = sample_list,
  nfeatures = 3000
)
length(features)
####################instead just run this:
sample_list <- readRDS("rds_objects/all_samples_QC.rds")                                                       
sample_list <- lapply(sample_list, function(x) {
  SCTransform(
    x,
    assay = "RNA",
    vars.to.regress = "percent.mt",
    verbose = TRUE
  )
})

sapply(sample_list, DefaultAssay)
sapply(sample_list, function(x) class(x[["SCT"]]))

saveRDS(sample_list, "all_samples_SCT_clean.rds")       
features <- SelectIntegrationFeatures(sample_list, nfeatures = 3000)
length(features)
sample_list <- PrepSCTIntegration(
  object.list = sample_list,
  anchor.features = features
)

anchors <- FindIntegrationAnchors(
  object.list = sample_list,
  normalization.method = "SCT",
  anchor.features = features
)

combined <- IntegrateData(
  anchors,
  normalization.method = "SCT"
)

