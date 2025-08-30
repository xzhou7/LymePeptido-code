#Single cell Lyme disease post merge
#Author: Xin Chen, Ph.D.
#Date Created: 010623
#Date Updated: 01/14/23
#After Find clusters On Cluster, run UMAP

#Intergrate Trim data
#link: https://satijalab.org/seurat/articles/integration_introduction.html
#set up working directory
#setwd("/Users/xzhou7/Library/CloudStorage/Box-Box/Xin.Chen.Shareable/R7_NR")

setwd(dir = "~/Library/CloudStorage/Box-Box/Xin.Chen.Shareable/R7_NR")
setwd(dir = "C:/Users/zhoux/Box/Xin.Chen.Shareable/R7_NR/")
#load necessary package
library(data.table)
library(Seurat)
library(ggplot2)
library(stringr)
library(patchwork)

getwd()
load("./Data/011423_R6_Step2.2_FindCluster.RData")
source("./00_colorKey.R")
#pbmc.integrated <- IntegrateData(anchorset = pbmc.anchors, dims = 1:30)

#pbmc.integrated <- ScaleData(pbmc.integrated, verbose = F)

#pbmc.integrated <- RunPCA(pbmc.integrated, npcs = 50, verbose = F)

#pbmc.integrated <-FindNeighbors(pbmc.integrated, dims= 1:10)

#pbmc.integrated <-FindClusters(pbmc.integrated, resolution = 1.2)

#according to the set idents, randomly subset each group to 3000 cells
# table(pbmc.integrated$subject)
# Idents(pbmc.integrated) <- "subject"
# pbmc.integrated.3000 <- subset(pbmc.integrated, downsample = 3000)
# table(pbmc.integrated.3000$subject)
# 
# Idents(pbmc.integrated.3000) <- "seurat_clusters"
# pbmc.integrated.3000 <- RunUMAP(pbmc.integrated.3000, dims = 1:10, verbose = T)
# 
# DimPlot(pbmc.integrated.3000)
FeaturePlot(pbmc.integrated.3000, "FCGR3A")

#downsampling to 5000 cells
Idents(pbmc.integrated) <- "subject"
pbmc.integrated.5000 <- subset(pbmc.integrated, downsample = 5000)
table(pbmc.integrated.5000$subject)

ElbowPlot(pbmc.integrated.5000)

Idents(pbmc.integrated.5000) <- "seurat_clusters"
pbmc.integrated.5000 <- RunUMAP(pbmc.integrated.5000, dims = 1:15, verbose = T)

DimPlot(pbmc.integrated.5000)

FeaturePlot(pbmc.integrated.5000, "FOXP3",min.cutoff = 0, sort.cell = T)

FeaturePlot(pbmc.integrated.5000, "CCL5",min.cutoff = 0, sort.cell = T)


FeaturePlot(pbmc.integrated.5000, "CCL5", min.cutoff = 0, sort.cell = T, split.by = "condition")


#save(pbmc.integrated, file = "./Data/012423_Step2.3_RunUMAP.RData")

save(pbmc.integrated.5000, file = "./Data/022123_Step2.3_Subset5000_RunUMAP.RData")


