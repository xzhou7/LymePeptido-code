#Recluster monocytes

#set up working directory
setwd(dir = "C:/Users/zhoux/Box/Xin.Chen.Shareable/R7_NR/")
setwd(dir = "~/Library/CloudStorage/Box-Box/Xin.Chen.Shareable/R7_NR")
getwd()
load("./Data/022723_Step2.3_Subset5000_RunUMAP.RData")

#load necessary package
library(data.table)
library(dplyr)
library(Seurat)
library(ggplot2)
library(stringr)
library(patchwork)
library(metap)
library(multtest)
library(Rcpp)
library(EnhancedVolcano)

#Select Myeloid for reclustering
table(pbmc.integrated.5000$seurat_clusters)

Myeloid <- subset(pbmc.integrated.5000, seurat_clusters %in% c("27", "18", "32", "22", "34", "16", "14", "24", "2", "28",
                                                               "9", "6", "30", "19", "37", "35", "27", "21"))
DefaultAssay(Myeloid) <- "integrated"
Myeloid_variable <- FindVariableFeatures(Myeloid, selection.method = "vst", nfeatures = 100,verbose = FALSE)
Myeloid_scale <- ScaleData(Myeloid_variable, verbose = T)
Myeloid_PCA <- RunPCA(Myeloid_scale, npcs = 50, verbose = T)
ElbowPlot(Myeloid_PCA)
Myeloid_UMAP <- RunUMAP(Myeloid_PCA, reduction = "pca", dims = 1:10, verbose = FALSE)
Myeloid_Nei <-FindNeighbors(Myeloid_UMAP, reduction = "pca", dims= 1:10)
Myeloid.integrated <-FindClusters(Myeloid_Nei, resolution = 0.3)

p1 <- DimPlot(Myeloid.integrated, label = T, split.by = "condition")
p1
ggsave(filename = "./Results/022823_Meyloid with nfeTURE100_dim1:10_res_0.3.pdf", p1, width = 8, height = 3, dpi = 300)

DefaultAssay(Myeloid) <- "antibody"
antibody.list <- row.names(Myeloid)
marker.plot <- FeaturePlot(Myeloid.integrated,features = "CD14.MPHIP9.CD14.AHS0037.pAbO",min.cutoff = 0, max.cutoff = 1000) + scale_colour_gradient(low = "#D3D3D3",high = "#8C1515")
marker.plot 
marker.plot <- FeaturePlot(Myeloid.integrated,features = "CD16.3G8.FCGR3A.AHS0053.pAbO",min.cutoff = 0, max.cutoff = 1000) + scale_colour_gradient(low = "#D3D3D3",high = "#8C1515")
marker.plot

#DefaultAssay(mono.integrated) <- "antibody"
#all.genes.rna <- rownames(mono.integrated)
#all.genes.rna[all.genes.rna %like% "CD14.MPHIP9.CD14.AHS0037.pAbO"]

#find markers for cluster
pbmc.markers <- FindAllMarkers(Myeloid.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
write.csv(pbmc.markers,"./Results/Myeloid.integrated gene markerswith nfeTURE100_dim1:10_res_0.3.csv", row.names = FALSE)

save(Myeloid.integrated, file = "./Data/Step4.1_022823_recluster_Myeloid.integrated_final.RData")
