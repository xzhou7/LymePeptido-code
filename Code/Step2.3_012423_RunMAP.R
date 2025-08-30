#Single cell Lyme disease post merge
#Author: Xin Chen, Ph.D.
#Date Created: 010623
#Date Updated: 01/24/23
#After Find clusters On Cluster, run UMAP

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
setwd("/oak/stanford/scg/lab_mmdavis/xx213/cluster_run/")
load("/oak/stanford/scg/lab_mmdavis/xx213/cluster_run/011423_R6_Step2.2_FindCluster.RData")

pbmc.integrated <- RunUMAP(pbmc.integrated, dims = 1:10, verbose = F)

save(pbmc.integrated, file = "/oak/stanford/scg/lab_mmdavis/xx213/cluster_run/012423_Step2.3_RunUMAP.RData")