#Overall Strategy
#Step 1 — Load each sample separately +
#Step 2 — QC per sample +
#Step 3 — SCTransform normalization
#Step 4 — Integration across 11 samples
#Step 5 — PCA + UMAP + clustering
#Step 6 — Identify major compartments
#Step 7 — Compare tumor vs normal
#Step 8 — Subcluster immune cells

#connecting remote host (bctl-gpu)
#working on VS-code
#module spider R/4.4.1-gfbf-2023b
module spider Seurat/5.1.0
module load R-bundle-Bioconductor/3.19-foss-2023b-R-4.4.1

library(Seurat)
library(dplyr)

#directories:
setwd("/main/my_folder/test_out") #working path/outputs
getwd()
data_dir <- "/main/public/paper1/cellranger" #data location, there are normal and cancer samples

sample_dirs <- list.dirs(data_dir, recursive = FALSE)
sample_dirs
sample_list <- list()

#testing for one sample only:
test_counts <- Read10X(data.dir = file.path(data_dir, "normal1"))
str(test_counts)
#> Formal class 'dgCMatrix' [package "Matrix"] with 6 slot...

#Use parallelization for SCTransform
library(future)
plan("multicore", workers = 8) # or however many cores you have

#running for all 11 samples:
for (sample_path in sample_dirs) {
   sample_name <- basename(sample_path)
  
  counts <- Read10X(data.dir = sample_path)
  
 
  seurat_obj <- CreateSeuratObject(counts = counts,
                                   project = sample_name,
                                   min.cells = 3,
                                   min.features = 200)
  
  seurat_obj$sample <- sample_name
  
  # Add condition metadata
  if (sample_name %in% c("normal1","normal2","normal3")) {
    seurat_obj$condition <- "Normal"
  } else {
    seurat_obj$condition <- "Tumor"
  }
  
  sample_list[[sample_name]] <- seurat_obj
}
##check cell numbers
sapply(sample_list, ncol)
#>normal1 normal2 normal3 patient1 patient2 patient3 patient4 patient5 9946 8250 6671 10013 12284 7730 8514 8377
sum(sapply(sample_list, ncol))
#>[1] 106115


for (sample_name in names(sample_list)) {
  #add mitochondrial percentage
  sample_list[[sample_name]][["percent.mt"]] <- PercentageFeatureSet(
    sample_list[[sample_name]],
    pattern = "^MT-"  #human
  )
  #QC filtering
  sample_list[[sample_name]] <- subset(
    sample_list[[sample_name]],
    subset = nFeature_RNA > 300 &
             nFeature_RNA < 6000 &
             percent.mt < 15
  )
   cat("Finished QC for:", sample_name, "\n")
}
#saving
#dir.create("rds_objects", showWarnings = FALSE)
#saveRDS(sample_list, file = "rds_objects/all_samples_QC.rds")


# normalisation stage using SCTransform
#install.packages('BiocManager')
#BiocManager::install('glmGamPoi') #for much faster implementation
library(glmGamPoi)

for (sample_name in names(sample_list)) {
  sample_list[[sample_name]] <- SCTransform(sample_list[[sample_name]],
                                  vars.to.regress = "percent.mt",
                                  verbose = FALSE)
}
##SCTransform takes a while ~ afew hours!!
features <- SelectIntegrationFeatures(object.list = sample_list, nfeatures = 3000)

sample_list <- PrepSCTIntegration(object.list = sample_list,
                                  anchor.features = features)

anchors <- FindIntegrationAnchors(object.list = sample_list,
                                  normalization.method = "SCT",
                                  anchor.features = features)

integrated <- IntegrateData(anchorset = anchors,
                            normalization.method = "SCT")
integrated <- RunPCA(integrated)
ElbowPlot(integrated)

integrated <- RunUMAP(integrated, dims = 1:30)

integrated <- FindNeighbors(integrated, dims = 1:30)
integrated <- FindClusters(integrated, resolution = 0.6)

DimPlot(integrated, label = TRUE)

DimPlot(integrated, group.by = "sample")
DimPlot(integrated, group.by = "condition")
FeaturePlot(integrated, features = c(
  "CD3D",      # T cells
  "MS4A1",     # B cells
  "LYZ",       # Myeloid
  "EPCAM",     # Epithelial (tumor likely here)
  "COL1A1",    # Fibroblast
  "PECAM1"     # Endothelial
))

markers <- FindAllMarkers(integrated, only.pos = TRUE)
epi <- subset(integrated, idents = "Epithelial")

Idents(epi) <- "condition"

deg_epi <- FindMarkers(epi,
                       ident.1 = "Tumor",
                       ident.2 = "Normal")
immune <- subset(integrated, idents = c("T cells","B cells","Myeloid","NK"))

immune <- RunPCA(immune)
immune <- RunUMAP(immune, dims = 1:20)
immune <- FindNeighbors(immune, dims = 1:20)
immune <- FindClusters(immune, resolution = 0.6)

DimPlot(immune, label = TRUE)
