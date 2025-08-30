#Single cell Lyme disease post merge
#Author: Xin Chen, Ph.D.
#Date Created: 010623
#Date Updated: 01/14/23
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
load("/oak/stanford/scg/lab_mmdavis/xx213/cluster_run/Step2.2_Intergrated_data_Anchors.RData")
pbmc.integrated <- IntegrateData(anchorset = pbmc.anchors, dims = 1:30)

pbmc.integrated <- ScaleData(pbmc.integrated, verbose = F)

pbmc.integrated <- RunPCA(pbmc.integrated, npcs = 50, verbose = F)

pbmc.integrated <-FindNeighbors(pbmc.integrated, dims= 1:10)

pbmc.integrated <-FindClusters(pbmc.integrated, resolution = 1.2)

save(pbmc.integrated, file = "./011423_R6_Step2.2_FindCluster.RData")
