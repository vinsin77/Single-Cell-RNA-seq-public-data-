#Overall Strategy
#Step 1 — Load each sample separately
#Step 2 — QC per sample
#Step 3 — SCTransform normalization
#Step 4 — Integration across 11 samples
#Step 5 — PCA + UMAP + clustering
#Step 6 — Identify major compartments
#Step 7 — Compare tumor vs normal
#Step 8 — Subcluster immune cells

library(Seurat)
library(dplyr)

sample_dirs <- list.dirs("data", recursive = FALSE)

sample_list <- list()

for (i in seq_along(sample_dirs)) {
  
  counts <- Read10X(data.dir = file.path(sample_dirs[i], "filtered_feature_bc_matrix"))
  
  sample_name <- basename(sample_dirs[i])
  
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
  
  sample_list[[i]] <- seurat_obj
}

for (i in 1:length(sample_list)) {
  
  sample_list[[i]][["percent.mt"]] <- PercentageFeatureSet(
    sample_list[[i]],
    pattern = "^MT-"
  )
  
  sample_list[[i]] <- subset(
    sample_list[[i]],
    subset = nFeature_RNA > 300 &
             nFeature_RNA < 6000 &
             percent.mt < 15
  )
}
for (i in 1:length(sample_list)) {
  
  sample_list[[i]][["percent.mt"]] <- PercentageFeatureSet(
    sample_list[[i]],
    pattern = "^MT-"
  )
  
  sample_list[[i]] <- subset(
    sample_list[[i]],
    subset = nFeature_RNA > 300 &
             nFeature_RNA < 6000 &
             percent.mt < 15
  )
}
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
