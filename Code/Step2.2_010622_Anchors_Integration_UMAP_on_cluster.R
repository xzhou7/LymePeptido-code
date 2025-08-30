#Single cell Lyme disease post merge
#Author: Xin Chen, Ph.D.
#Date Created: Nov.16.2021
#Date Updated: 11/16/21
#Run On Cluster

#Intergrate Trim data
#link: https://satijalab.org/seurat/articles/integration_introduction.html
#set up working directory
setwd("/oak/stanford/scg/lab_mmdavis/xx213/cluster_run/")

#load necessary package
library(data.table)
library(Seurat)
library(ggplot2)
library(stringr)
library(patchwork)

getwd()

load("/oak/stanford/scg/lab_mmdavis/xx213/cluster_run/Step2.010623_beforeMerge.RData")

pbmc.anchors <- FindIntegrationAnchors(object.list = pbmc.list, anchor.features = 350, dims = 1:30)

save(pbmc.anchors, file = "./Step2.2_Intergrated_data_Anchors.RData")

pbmc.anchors

pbmc.integrated <- IntegrateData(anchorset = pbmc.anchors, dims = 1:30)

pbmc.integrated <- ScaleData(pbmc.integrated, verbose = F)

pbmc.integrated <- RunPCA(pbmc.integrated, npcs = 50, verbose = F)

pbmc.integrated <-FindNeighbors(pbmc.integrated, dims= 1:10)

pbmc.integrated <-FindClusters(pbmc.integrated, resolution = 1.2)

pbmc.integrated <- RunUMAP(pbmc.integrated, dims = 1:10, verbose = F)


save(pbmc.integrated, file = "./072522_R6_Step3.1_subset.id.RData")
