
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(ggplot2)
library(cowplot)
library(dplyr)
library(RColorBrewer)

#Analysis after QC 
load("seurat_filtered.RData")

# Splitting Seurat object by original identity and processing each subset
obj.list_f <- SplitObject(filtered_seurat, split.by = "orig.ident")
Rh.list <- obj.list_f
new_list <- list()
for (i in 1:length(Rh.list)) {
    obj.list    = Rh.list[[i]]
    obj.list[["RNA"]]@meta.features <- data.frame(row.names = rownames(obj.list[["RNA"]]))
    obj.list <- NormalizeData(obj.list, verbose = FALSE)
    obj.list <- FindVariableFeatures(obj.list, selection.method = "vst", 
    nfeatures = 3000, verbose = FALSE)
    new_list[[i]] <-  obj.list
}

# Assigning names to processed subsets
names(new_list) <- names(Rh.list)

# Selecting integration features and performing data integration
features <- SelectIntegrationFeatures(object.list = new_list)
new_list <- lapply(X = new_list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})

# Finding integration anchors, integrating data, and performing dimensionality reduction and cluster identification 
Rh.anchors <- FindIntegrationAnchors(object.list = new_list, anchor.features = features, reduction = "rpca" ,   dims = 1:50)
Rh.integrated <- IntegrateData(anchorset = Rh.anchors, dims = 1:50)
Rh.integrated <- ScaleData(Rh.integrated, verbose = FALSE)
Rh.integrated <- RunPCA(Rh.integrated, verbose = FALSE)
Rh.integrated <- RunUMAP(Rh.integrated, dims = 1:50)
Rh.integrated <- FindNeighbors(Rh.integrated, reduction = "pca", dims = 1:30)
Rh.integrated  <- FindClusters(Rh.integrated, resolution = 0.5)

# Plotting UMAP visualization
tiffc("umap.tiff", units="in", width=9, height=9, res=300)
DimPlot(Rh.integrated, group.by = c("orig.ident","cluster"),raster=FALSE)
dev.off()

# Save the integrated Seurat object
saveRDS(Rh.integrated, "Rh_seurat_integrated.RData")

