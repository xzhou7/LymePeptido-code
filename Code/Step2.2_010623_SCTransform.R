#Single cell Lyme disease Batch setup
#Author: Xin Chen, Ph.D.
#Date Created: 01623
# SCTransform

#set up working directory
setwd("/Users/xinchen/Library/CloudStorage/Box-Box/Xin.Chen.Shareable/R7_NR")
getwd()
load("./Data/Step2.1_010623_QC_trim_by_batch.RData")

#load necessary package
library(data.table)
library(Seurat)
library(ggplot2)
library(stringr)
library(patchwork)


#Split the dataset into a list of 5 seurat batches
pbmc.list <- SplitObject(pbmc,split.by = "batch")
pbmc.list

#Normalize and identify variable features for each dataset independently
for (i in 1:length(pbmc.list)) {
  pbmc.list[[i]] <- SCTransform(pbmc.list[[i]], verbose = T)
  pbmc.list[[i]] <- FindVariableFeatures(pbmc.list[[i]], selection.method = "vst", nfeatures = 350,verbose = FALSE)
}

save(pbmc.list, file = "./Data/Step2.010623_beforeMerge.RData")


#Perform Integration (no memory for here, move to cluster)
scp ./Desktop/Step2.010623_beforeMerge.RData xchen7@scg.stanford.edu:/labs/mmdavis/xx213/cluster_run/
pbmc.anchors <- FindIntegrationAnchors(object.list = pbmc.list, anchor.features = 350, dims = 1:30)

save(pbmc.anchors, file = "./Data/Step2.2_Intergrated_data_Anchors.RData")

#Creats an 'integrated' data assay
pbmc.integrated <- IntegrateData(anchorset = pbmc.anchors, dims = 1:30)

save(pbmc.integrated, file = "./Data/Step2.2_Intergrated_data_Integrated.RData")

