#Peptidoglycan_spleen Data Analysis 
#Apr 8, 2024
#Updated Apr 15, 2025

#Overall Tutorial 
#https://bioconductor.org/books/release/OSCA/

#UMAP
#https://mp.weixin.qq.com/s/bZ9DYkKOKYo2w1OaoBxc1w 

#ClusterProfile
#https://mp.weixin.qq.com/s?__biz=MzUzMTEwODk0Ng==&mid=2247513807&idx=1&sn=677cdfe6558dc6b1406daa1d56f2b9e6&chksm=fa4573f2cd32fae4333f310732bbd46b76981f2c0e2379d4121c470c2cbe4539a821305be06e&scene=132#wechat_redirect

#Monocle3
#http://htmlpreview.github.io/?https://github.com/satijalab/seurat-wrappers/blob/master/docs/monocle3.html

#Palo
#https://winnie09.github.io/Wenpin_Hou/pages/Palo.html

#Vison
#https://yoseflab.github.io/VISION/articles/VISION-vignette.html 

library(dplyr)
library(Seurat)
library(patchwork)
library(data.table)
library(dplyr)
library(ggplot2)
library(stringr)
library(metap)
library(multtest)
library(Rcpp)
library(sctransform)
library(glmGamPoi)
library(cowplot)
library(ggunchull)
library(tidydr)
library(ggsci)
library(Cairo)
library(DoubletFinder)
library(harmony)

#mac
setwd("/Users/xzhou7/Library/CloudStorage/Dropbox/lyme_disease/Pepti_Spleen_singlecell/")

#windows
setwd("D:/Dropbox/lyme_disease/Pepti_Spleen_singlecell/")

getwd()

#load control cells
pbmc.data <- Read10X(data.dir = "./Raw_Data/control/filtered_feature_bc_matrix/")
pbmc.data
pbmc <- CreateSeuratObject(counts = pbmc.data$`Gene Expression`, project = "spleen-control")
pbmc[["ADT"]]  <-CreateAssayObject(counts = pbmc.data$`Antibody Capture`, project="spleen-control")
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc[["percent.RP"]] <- PercentageFeatureSet(pbmc, pattern = "^RP")
#control <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 20 & nCount_RNA > 300 & nCount_RNA < 40000)
#VlnPlot(control, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

DefaultAssay(pbmc) <- "ADT"
VlnPlot(pbmc, c("Hasgtag 1-TotalSeqC","Hasgtag 2-TotalSeqC","Hasgtag 3-TotalSeqC"))

pbmc[["percent.1"]] <- PercentageFeatureSet(pbmc, pattern = "^Hasgtag 1")
pbmc[["percent.2"]] <- PercentageFeatureSet(pbmc, pattern = "^Hasgtag 2")
pbmc[["percent.3"]] <- PercentageFeatureSet(pbmc, pattern = "^Hasgtag 3")

pbmc$Hashtag <- ifelse(pbmc$percent.1 > pbmc$percent.2 & pbmc$percent.1 > pbmc$percent.3, "Hashtag1",
                       ifelse(pbmc$percent.2 > pbmc$percent.1 & pbmc$percent.2 > pbmc$percent.3, "Hashtag2", "Hashtag3"))

#donor1
donor1 <- subset(pbmc, Hashtag == "Hashtag1")
DefaultAssay(donor1) <- "RNA"

donor1[["percent.mt"]] <- PercentageFeatureSet(donor1, pattern = "^MT-")
donor1[["percent.RP"]] <- PercentageFeatureSet(donor1, pattern = "^RP")
donor1[[]]
plot1 <- FeatureScatter(donor1, feature1 = "nCount_RNA", feature2 = "percent.mt") + geom_hline(yintercept = 10)
plot2 <- FeatureScatter(donor1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+ geom_hline(yintercept = 8000)
plot1 + plot2

donor1 <- subset(donor1, subset = nFeature_RNA > 1000 & nFeature_RNA < 8000 & percent.mt < 7.5 & nCount_RNA > 200 & nCount_RNA < 100000)
donor1


#donor 2
donor2 <- subset(pbmc, Hashtag == "Hashtag2")
DefaultAssay(donor2) <- "RNA"

donor2[["percent.mt"]] <- PercentageFeatureSet(donor2, pattern = "^MT-")
donor2[["percent.RP"]] <- PercentageFeatureSet(donor2, pattern = "^RP")

plot1 <- FeatureScatter(donor2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(donor2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

donor2 <- subset(donor2, subset = nFeature_RNA > 1000 & nFeature_RNA < 8000 & percent.mt < 7.5 & nCount_RNA > 200 & nCount_RNA < 100000)

#donor 3
donor3 <- subset(pbmc, Hashtag == "Hashtag3")
DefaultAssay(donor3) <- "RNA"

donor3[["percent.mt"]] <- PercentageFeatureSet(donor3, pattern = "^MT-")
donor3[["percent.RP"]] <- PercentageFeatureSet(donor3, pattern = "^RP")

plot1 <- FeatureScatter(donor3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(donor3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

donor3 <- subset(donor3, subset = nFeature_RNA > 1000 & nFeature_RNA < 8000 & percent.mt < 7.5 & nCount_RNA > 200 & nCount_RNA < 100000)

donor1 <- SCTransform(donor1, vst.flavor = "v2", verbose = FALSE) %>% RunPCA(npcs = 15, verbose = T)
donor1 <- RunUMAP(donor1, dims = 1:15)
donor1 <- FindNeighbors(donor1, dims = 1:15, verbose = FALSE)
donor1 <- FindClusters(donor1, resolution = 0.8)
# individual clusters
DimPlot(donor1, reduction = "umap")

donor2 <- SCTransform(donor2, vst.flavor = "v2", verbose = FALSE) %>% RunPCA(npcs = 15, verbose = T)
donor2 <- RunUMAP(donor2, dims = 1:15)
donor2 <- FindNeighbors(donor2, dims = 1:15, verbose = FALSE)
donor2 <- FindClusters(donor2, resolution = 0.8)
DimPlot(donor2, reduction = "umap")

donor3 <- SCTransform(donor3, vst.flavor = "v2", verbose = FALSE) %>% RunPCA(npcs = 15, verbose = T)
donor3 <- RunUMAP(donor3, dims = 1:15)
donor3 <- FindNeighbors(donor3, dims = 1:15, verbose = FALSE)
donor3 <- FindClusters(donor3, resolution = 0.8)
DimPlot(donor3, reduction = "umap")

#load peptidoglycan cells
pbmc.data.pep <- Read10X(data.dir = "./Raw_Data/peptidoglycan/filtered_feature_bc_matrix/")
pbmc.data.pep
pbmc.pep <- CreateSeuratObject(counts = pbmc.data.pep$`Gene Expression`, project = "spleen_peptidoglycan")
pbmc.pep[["ADT"]]  <-CreateAssayObject(counts = pbmc.data.pep$`Antibody Capture`, project="spleen_peptidoglycan")
pbmc.pep[["percent.mt"]] <- PercentageFeatureSet(pbmc.pep, pattern = "^MT-")
pbmc.pep[["percent.RP"]] <- PercentageFeatureSet(pbmc.pep, pattern = "^RP")
#peptidoglycan <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 20 & nCount_RNA > 300 & nCount_RNA < 40000)
#VlnPlot(peptidoglycan, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

pbmc.pep
DefaultAssay(pbmc.pep) <- "ADT"
VlnPlot(pbmc.pep, c("Hasgtag 1-TotalSeqC","Hasgtag 2-TotalSeqC","Hasgtag 3-TotalSeqC"))

pbmc.pep[["percent.1"]] <- PercentageFeatureSet(pbmc.pep, pattern = "^Hasgtag 1")
pbmc.pep[["percent.2"]] <- PercentageFeatureSet(pbmc.pep, pattern = "^Hasgtag 2")
pbmc.pep[["percent.3"]] <- PercentageFeatureSet(pbmc.pep, pattern = "^Hasgtag 3")

pbmc.pep$Hashtag <- ifelse(pbmc.pep$percent.1 > pbmc.pep$percent.2 & pbmc.pep$percent.1 > pbmc.pep$percent.3, "Hashtag1",
                         ifelse(pbmc.pep$percent.2 > pbmc.pep$percent.1 & pbmc.pep$percent.2 > pbmc.pep$percent.3, "Hashtag2", "Hashtag3"))
  
table(pbmc.pep$Hashtag)
pbmc.pep[[]]

pep_donor1 <- subset(pbmc.pep, Hashtag == "Hashtag1")
pep_donor2 <- subset(pbmc.pep, Hashtag == "Hashtag2")
pep_donor3 <- subset(pbmc.pep, Hashtag == "Hashtag3")

DefaultAssay(pep_donor1) <- "RNA"

pep_donor1[["percent.mt"]] <- PercentageFeatureSet(pep_donor1, pattern = "^MT-")
pep_donor1[["percent.RP"]] <- PercentageFeatureSet(pep_donor1, pattern = "^RP")
pep_donor1[[]]
plot1 <- FeatureScatter(pep_donor1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pep_donor1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

pep_donor1 <- subset(pep_donor1, subset = nFeature_RNA > 1000 & nFeature_RNA < 8000 & percent.mt < 7.5 & nCount_RNA > 200 & nCount_RNA < 100000)
pep_donor1

#donor 2
DefaultAssay(pep_donor2) <- "RNA"

pep_donor2[["percent.mt"]] <- PercentageFeatureSet(pep_donor2, pattern = "^MT-")
pep_donor2[["percent.RP"]] <- PercentageFeatureSet(pep_donor2, pattern = "^RP")

plot1 <- FeatureScatter(pep_donor2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pep_donor2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

pep_donor2 <- subset(pep_donor2, subset = nFeature_RNA > 1000 & nFeature_RNA < 8000 & percent.mt < 7.5 & nCount_RNA > 200 & nCount_RNA < 100000)

#donor 3
DefaultAssay(pep_donor3) <- "RNA"

pep_donor3[["percent.mt"]] <- PercentageFeatureSet(pep_donor3, pattern = "^MT-")
pep_donor3[["percent.RP"]] <- PercentageFeatureSet(pep_donor3, pattern = "^RP")

plot1 <- FeatureScatter(pep_donor3, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pep_donor3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

pep_donor3 <- subset(pep_donor3, subset = nFeature_RNA > 1000 & nFeature_RNA < 8000 & percent.mt < 7.5 & nCount_RNA > 200 & nCount_RNA < 100000)

pep_donor1 <- SCTransform(pep_donor1, vst.flavor = "v2", verbose = FALSE) %>% RunPCA(npcs = 15, verbose = T)
pep_donor1 <- RunUMAP(pep_donor1, dims = 1:15)
pep_donor1 <- FindNeighbors(pep_donor1, dims = 1:15, verbose = FALSE)
pep_donor1 <- FindClusters(pep_donor1, resolution = 0.8)
DimPlot(pep_donor1, reduction = "umap")

ElbowPlot(pep_donor1)

#FeaturePlot(pep_donor2, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
#                               "CD8A"))

pep_donor2 <- SCTransform(pep_donor2, vst.flavor = "v2", verbose = FALSE) %>% RunPCA(npcs = 15, verbose = T)
pep_donor2 <- RunUMAP(pep_donor2, dims = 1:15)
pep_donor2 <- FindNeighbors(pep_donor2, dims = 1:15, verbose = FALSE)
pep_donor2 <- FindClusters(pep_donor2, resolution = 0.8)
DimPlot(pep_donor2, reduction = "umap")
#FeaturePlot(pep_donor2, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP",
#                                     "CD8A"))

pep_donor3 <- SCTransform(pep_donor3, vst.flavor = "v2", verbose = FALSE) %>% RunPCA(npcs = 15, verbose = T)
pep_donor3 <- RunUMAP(pep_donor3, dims = 1:15)
pep_donor3 <- FindNeighbors(pep_donor3, dims = 1:15, verbose = FALSE)
pep_donor3 <- FindClusters(pep_donor3, resolution = 0.8)
DimPlot(pep_donor3, reduction = "umap")

#########################################
#rm(pbmc.data)
spleen.list <- list(control_1 = donor1, control_2= donor2, control_3=donor3, pep_1 = pep_donor1, pep_2= pep_donor2, pep_3=pep_donor3)
spleen.list

# 
# ########################################
# #DoubletFinder
# # --- Process Donor 1 (Control) ---
# 
# print("Processing Donor 1 (Control)")
# 
# # 1. Define the expected doublet rate - ADJUST THIS VALUE!
# # Example: Assuming ~4% doublet rate for this specific sample size.
# # Calculate based on 10x recommendations (e.g., 0.075 for 10k cells)
# # doubletRate_d1 <- (ncol(donor1) / 10000) * 0.075 # Example calculation
# doubletRate_d1 <- 0.04 # <<< Replace with your estimated rate for donor1
# 
# # 2. Optimize pK value
# # Note: This step can be computationally intensive
# sweep.res.list_donor1 <- paramSweep(donor1, PCs = 1:15, sct = TRUE) # Use sct = TRUE because you used SCTransform
# sweep.stats_donor1 <- summarizeSweep(sweep.res.list_donor1, GT = FALSE)
# bcmvn_donor1 <- find.pK(sweep.stats_donor1)
# 
# # Find the pK value associated with the maximum BCmetric
# pK_optimal_d1 <- as.numeric(as.character(bcmvn_donor1$pK[which.max(bcmvn_donor1$BCmetric)]))
# print(paste("Optimal pK for Donor 1:", pK_optimal_d1))
# 
# # 3. Calculate the expected number of doublets
# nExp_poi_d1 <- round(doubletRate_d1 * ncol(donor1))
# print(paste("Expected doublets for Donor 1:", nExp_poi_d1))
# 
# # 4. Run DoubletFinder
# # Ensure your cluster resolution ('resolution = 0.8') generated a column named 'sct_snn_res.0.8'
# # Check with: colnames(donor1@meta.data)
# # If the name is different, adjust 'annotations_col' accordingly.
# annotations_col_d1 <- "SCT_snn_res.0.8"  # Change if your cluster column name is different
# donor1 <- doubletFinder(donor1,
#                            PCs = 1:15,
#                            pN = 0.25, # Default parameter, usually doesn't need changing
#                            pK = pK_optimal_d1,
#                            nExp = nExp_poi_d1,
#                            reuse.pANN = NULL, # Set to NULL for the first run
#                            sct = TRUE) # Use sct = TRUE
# 
# # 5. Check the results and find the classification column name
# # DoubletFinder adds metadata columns. The classification column usually looks like:
# # DF.classifications_pN_pK_nExp
# print("Metadata columns after DoubletFinder:")
# print(tail(colnames(donor1@meta.data)))
# 
# # Find the classification column name programmatically (adjust pattern if needed)
# df_classification_col_d1 <- grep("DF.classifications", colnames(donor1@meta.data), value = TRUE)[1]
# if (is.na(df_classification_col_d1)) {
#   stop("Could not find DoubletFinder classification column for Donor 1!")
# }
# print(paste("Using classification column:", df_classification_col_d1))
# 
# # Optional: Visualize doublets on UMAP
# DimPlot(donor1, group.by = df_classification_col_d1, reduction = "umap")
# FeaturePlot(donor1, "CD14")
# 
# # Optional: Check distribution of classifications
# print("DoubletFinder Classification Summary for Donor 1:")
# print(table(donor1[[df_classification_col_d1]]))
# 
# # 6. Remove doublets
# print(paste("Cells before doublet removal (Donor 1):", ncol(donor1)))
# donor1_singlets <- subset(donor1, subset = !!sym(df_classification_col_d1) == "Singlet")
# print(paste("Cells after doublet removal (Donor 1):", ncol(donor1_singlets)))
# 
# # Replace the original object with the singlet-only object
# donor1 <- donor1_singlets
# rm(donor1_singlets) # Clean up temporary object
# 
# #Harmony#
# for (name in names(spleen.list)) {
#   spleen.list[[name]]$orig.ident <- name
# }
# 
# spleen.merged <- merge(spleen.list[[1]], y = spleen.list[-1], add.cell.ids = names(spleen.list), project = "Spleen")
# spleen.merged <- SCTransform(spleen.merged, verbose = T)
# spleen.merged <- RunPCA(spleen.merged, verbose = T)
# 
# spleen.merged <- RunHarmony(spleen.merged, group.by.vars = "orig.ident")
# 
# spleen.merged <- RunHarmony(
#   object = spleen.merged,
#   group.by.vars = "orig.ident",
#   reduction.use = "pca",         
#   dims.use = 1:30,               
#   theta = 2,                     
#   plot_convergence = TRUE        
# )
# 
# spleen.merged <- RunUMAP(spleen.merged, reduction = "harmony", min.dist = 0.8,dims = 1:20)
# spleen.merged <- FindNeighbors(spleen.merged, reduction = "harmony", dims = 1:20)
# spleen.merged <- FindClusters(spleen.merged, resolution = 0.5)
# 
# DimPlot(spleen.merged, reduction = "umap", group.by = "orig.ident")  # batch check
# DimPlot(spleen.merged, reduction = "umap", label = TRUE)             # cluster labels
# 
# #Harmony End#
 

features <- SelectIntegrationFeatures(object.list = spleen.list, nfeatures = 4000)
spleen.list <- PrepSCTIntegration(object.list = spleen.list, anchor.features = features)
#this step take time
spleen.anchors <- FindIntegrationAnchors(object.list = spleen.list, normalization.method = "SCT", anchor.features = features,reference = c(1, 2, 3))
spleen.combined.sct <- IntegrateData(anchorset = spleen.anchors, normalization.method = "SCT")
spleen.combined.sct <- FindVariableFeatures(spleen.combined.sct, selection.method = "vst", nfeatures = 4000)
spleen.combined.sct
save(spleen.combined.sct, file = "./Integreted_4000Feature.RData")

features <- SelectIntegrationFeatures(object.list = spleen.list, nfeatures = 5000)
spleen.list <- PrepSCTIntegration(object.list = spleen.list, anchor.features = features)
#this step take time
spleen.anchors <- FindIntegrationAnchors(object.list = spleen.list, normalization.method = "SCT", anchor.features = features,reference = c(1, 2, 3))
spleen.combined.sct2 <- IntegrateData(anchorset = spleen.anchors, normalization.method = "SCT")
spleen.combined.sct2 <- FindVariableFeatures(spleen.combined.sct, selection.method = "vst", nfeatures = 5000)
spleen.combined.sct2
save(spleen.combined.sct2, file = "./Integreted_5000Feature.RData")


load("./Integreted_4000Feature.RData")
load("./Integreted_5000Feature.RData")
spleen.combined.sct.1 <- spleen.combined.sct
rm(spleen.combined.sct)
spleen.combined.sct.2 <- spleen.combined.sct2
rm(spleen.combined.sct2)

load("./RObject/Integreted_archive.RData")


spleen.combined.sct <- RunPCA(spleen.combined.sct, verbose = FALSE)
spleen.combined.sct <- RunUMAP(spleen.combined.sct, dims = 1:30, min.dist = 0.7, verbose = FALSE)
spleen.combined.sct <- FindNeighbors(spleen.combined.sct, dims = 1:30, verbose = FALSE)
spleen.combined.sct <- FindClusters(spleen.combined.sct,resolution = 0.6, verbose = FALSE)

spleen.combined.sct.1 <- RunPCA(spleen.combined.sct.1, verbose = FALSE)
spleen.combined.sct.1 <- RunUMAP(spleen.combined.sct.1, dims = 1:30, min.dist = 0.7, verbose = FALSE)
spleen.combined.sct.1 <- FindNeighbors(spleen.combined.sct.1, dims = 1:30, verbose = FALSE)
spleen.combined.sct.1 <- FindClusters(spleen.combined.sct.1,resolution = 0.7, verbose = FALSE)

spleen.combined.sct.2 <- RunPCA(spleen.combined.sct.2, verbose = FALSE)
spleen.combined.sct.2 <- RunUMAP(spleen.combined.sct.2, dims = 1:30,  verbose = FALSE)
spleen.combined.sct.2 <- FindNeighbors(spleen.combined.sct.2, dims = 1:30, verbose = FALSE)
spleen.combined.sct.2 <- FindClusters(spleen.combined.sct.2,resolution = 0.8, verbose = FALSE)
# 
# #check which integration works better
# library(kBET)
# embedding <- Embeddings(spleen.combined.sct, reduction = "pca")[, 1:10]
# batch <- spleen.combined.sct$Hashtag
# #note: the next step takes time
# kbet_res <- kBET(embedding, batch, testSize = 500)
# mean(kbet_res$stats$kBET.observed)  # Closer to 1 is better
# 
# embedding_1 <- Embeddings(spleen.combined.sct.1, reduction = "pca")[, 1:10]
# batch_1 <- spleen.combined.sct.1$Hashtag  # or use $sample if that's your batch indicator
# kbet_res_1 <- kBET(embedding_1, batch_1, testSize = 500)
# mean(kbet_res_1$stats$kBET.observed)
# 
# embedding_2 <- Embeddings(spleen.combined.sct.2, reduction = "pca")[, 1:10]
# batch_2 <- spleen.combined.sct.2$Hashtag  # or use $sample if appropriate
# kbet_res_2 <- kBET(embedding_2, batch_2, testSize = 500)
# mean(kbet_res_2$stats$kBET.observed)

DimPlot(spleen.combined.sct, label = T) + NoLegend() + coord_fixed()
spleen.combined.sct$sample <- paste(spleen.combined.sct$orig.ident, spleen.combined.sct$Hashtag, sep = "")
DimPlot(spleen.combined.sct, label = F, split.by = "sample", ncol=3) + NoLegend() + coord_fixed()
spleen.combined.sct[[]]

spleen.combined.sct$sample <- as.character(spleen.combined.sct$sample)

sample_rename <- c(
  "spleen_controlHashtag1"       = "control_1",
  "spleen_controlHashtag2"       = "control_2",
  "spleen_controlHashtag3"       = "control_3",
  "spleen_peptidoglycanHashtag1" = "Pepti_1",
  "spleen_peptidoglycanHashtag2" = "Pepti_2",
  "spleen_peptidoglycanHashtag3" = "Pepti_3")

spleen.combined.sct@meta.data$sample <- sample_rename[as.character(spleen.combined.sct@meta.data$sample)]

VlnPlot(spleen.combined.sct, features = c("CD3D", "CD14", "FOXP3"), group.by ="sample",ncol = 3, pt.size = 0)
housekeeping_genes <- c("ACTB", "GAPDH", "B2M", "sct_RPL13A", "RPS18", "EEF1A1", "HPRT1", "PPIA")

VlnPlot(spleen.combined.sct, features = housekeeping_genes, group.by = "sample", pt.size = 0, ncol = 4)

#save(spleen.combined.sct, file = "./Integreted_RunUMAP_041625.RData")

DimPlot(spleen.combined.sct.1, label = T) + NoLegend() + coord_fixed()
spleen.combined.sct.1$sample <- paste(spleen.combined.sct.1$orig.ident, spleen.combined.sct.1$Hashtag, sep = "")
sample_rename2 <- c(
  "spleen-controlHashtag1"       = "control_1",
  "spleen-controlHashtag2"       = "control_2",
  "spleen-controlHashtag3"       = "control_3",
  "spleen_peptidoglycanHashtag1" = "Pepti_1",
  "spleen_peptidoglycanHashtag2" = "Pepti_2",
  "spleen_peptidoglycanHashtag3" = "Pepti_3")
spleen.combined.sct.1@meta.data$sample <- sample_rename2[as.character(spleen.combined.sct.1@meta.data$sample)]

DimPlot(spleen.combined.sct.1, label = F, split.by = "sample", ncol=3) + NoLegend() + coord_fixed()
VlnPlot(spleen.combined.sct.1, features = c("CD3D", "CD14", "FOXP3"), group.by ="sample", ncol = 3, pt.size = 0)
VlnPlot(spleen.combined.sct.1, features = housekeeping_genes, group.by = "sample", pt.size = 0, ncol = 4)

###
DimPlot(spleen.combined.sct.2, label = T) + NoLegend() + coord_fixed()
spleen.combined.sct.2$sample <- paste(spleen.combined.sct.2$orig.ident, spleen.combined.sct.2$Hashtag, sep = "")
spleen.combined.sct.2@meta.data$sample <- sample_rename2[as.character(spleen.combined.sct.2@meta.data$sample)]
DimPlot(spleen.combined.sct.2, label = F, split.by = "sample", ncol=3) + NoLegend() + coord_fixed()

FeatureTarget <- "FCGR3A"
p1 <- FeaturePlot(spleen.combined.sct,FeatureTarget, order = T,min.cutoff = 0) + NoLegend() + coord_fixed()
p1.1 <- DimPlot(spleen.combined.sct, label = T) + NoLegend() + coord_fixed() + ggtitle("3000")
p2 <- FeaturePlot(spleen.combined.sct.1,FeatureTarget, order = T,min.cutoff = 0) + NoLegend() + coord_fixed()
p2.1 <- DimPlot(spleen.combined.sct.1, label = T) + NoLegend() + coord_fixed()+ ggtitle("4000")
p3 <- FeaturePlot(spleen.combined.sct.2,FeatureTarget, order = T,min.cutoff = 0) + NoLegend() + coord_fixed()
p3.1 <- DimPlot(spleen.combined.sct.2, label = T) + NoLegend() + coord_fixed() + ggtitle("5000")

p.compare <- (p1 + p1.1) | (p2 + p2.1) | (p3 + p3.1)
p.compare

#decision: use 5000 anchor, for the downsteam analysis 
rm(spleen.combined.sct)
rm(spleen.combined.sct.1)
spleen.combined.sct <- spleen.combined.sct.2
rm(spleen.combined.sct.2)


VlnPlot(spleen.combined.sct, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

DimPlot(spleen.combined.sct, label = TRUE) + NoLegend()

# Visualize canonical marker genes on the sctransform embedding.
FeaturePlot(spleen.combined.sct, features = c("CD8A", "CD3D", "MS4A1", "CD79A", "NKG7","CD4", "CD14","SELL","FOXP3"), pt.size = 0.2,
            ncol = 3, min.cutoff = 0)

FeaturePlot(spleen.combined.sct, features = c("CD1C","FCGR3A","CD14","SELL"), pt.size = 0.2,
            ncol = 3, min.cutoff = 0, order = T)


FeaturePlot(spleen.combined.sct, features = c("CXCL2","CCL3","NAMPT"), pt.size = 0.2,
            ncol = 3, min.cutoff = 0, order = T)

FeaturePlot(spleen.combined.sct, features = c("CD27"), pt.size = 0.2,
            ncol = 1, min.cutoff = 0, order = T) + coord_equal()

FeaturePlot(spleen.combined.sct, features = c("CD8A", "CD3D", "MS4A1", "CD79A", "NKG7","CD4", "CD14","SELL","FOXP3","CD1C","FCGR3A", 
                                              "CD3E", "CCR7", "S100A4","TBX21","RORC", "IFNG", "IL7R", "NCAM1", "ITGAM","CD38", "TIGIT","ITGAL", "SELL","CTLA4","ANXA5","CCL5","KLRB1","TRDC",
                                              "GZMA", "GZMK","GZMB", "PRF1", "FOXP3", "IL2RA","IFI6", "STAT1","MX1","GATA3","FCER1G","FGFBP2",
                                              "ITGB1","KIR2DL1","CD160","CD3D","KLRC2","GZMH"), pt.size = 0.2,ncol = 3, min.cutoff = 0)

FeaturePlot(spleen.combined.sct, features = c("MKI67", "TOP2A", "CCNB1", "CDK1"), pt.size = 0, min.cutoff = 0)

Marker.list <- c("CD8A", "CD3D", "MS4A1", "CD79A", "NKG7","CD4", "CD14","SELL","FOXP3","IL10"."CD1C","FCGR3A", "HLA-DRA","IGHG1","IGHD","BCL6", "LAG3", "HAVCR2", "CCR5", 
                 "CD3E", "CCR7", "S100A4","TBX21","RORC", "IFNG", "IL7R", "NCAM1", "ITGAM","CD38", "TIGIT","ITGAL", "SELL","CTLA4","ANXA5","CCL5","KLRB1","TRDC","IGHV3-23",
                 "GZMA", "GZMK","GZMB", "PRF1", "FOXP3", "IL2RA","IFI6", "STAT1","MX1","GATA3","FCER1G","FGFBP2","GZMH","CCL4","FSCN1","AICDA","CD27","IL21",
                 "ITGB1","KIR2DL1","CD160","CD3D","KLRC2", "IGLV3-1", "sct_KLF4", "MKI67", "TOP2A", "CCNB1", "CDK1","MZB1","XBP1","IRF8","PAX5", "CXCR4")

DefaultAssay(spleen.combined.sct) <- "integrated"

FeaturePlot(spleen.combined.sct, features = c("CD1C","FCGR3A"), pt.size = 0.2,
            ncol = 3, min.cutoff = 0)

FeaturePlot(spleen.combined.sct, features = c("IGHV2-26"),order=T, pt.size = 0.2, min.cutoff = 0)

#Integrated Markers
for (i in 1:length(Marker.list)){
  marker <- Marker.list[i]
  print(marker)
  marker.plot <- FeaturePlot(spleen.combined.sct,features = marker, min.cutoff = 0, order=T) + scale_colour_gradient(low = "#D3D3D3",high = "#8C1515")
  path <- paste0("./Markers/Inte.",marker,".pdf")
  ggsave(path, marker.plot, width = 4, height = 3.5, dpi = 300)
}

DimPlot(spleen.combined.sct, label = T, split.by = "orig.ident") + NoLegend()

DefaultAssay(spleen.combined.sct) <- "SCT"
Idents(spleen.combined.sct) <- "orig.ident"

gene.list <- row.names(spleen.combined.sct)
gene.list[str_detect(gene.list, "IL17")]

#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2976349/
#有变???
FeaturePlot(spleen.combined.sct, features = "sct_CXCL1", pt.size = 0.2, ncol = 3, min.cutoff = 0, max.cutoff = 1, split.by = "orig.ident", order = T) + coord_equal()
FeaturePlot(spleen.combined.sct, features = "sct_CXCL2", pt.size = 0.2, ncol = 3, min.cutoff = 0, max.cutoff = 1, split.by = "orig.ident", order = T) + coord_equal()
FeaturePlot(spleen.combined.sct, features = "CXCL3", pt.size = 0.2, ncol = 3, min.cutoff = 0, max.cutoff = 1, split.by = "orig.ident", order = T) + coord_equal()
FeaturePlot(spleen.combined.sct, features = "CXCL8", pt.size = 0.2, ncol = 3, min.cutoff = 0, max.cutoff = 1, split.by = "orig.ident", order = T) + coord_equal()
FeaturePlot(spleen.combined.sct, features = "sct_IL17F", pt.size = 0.2, ncol = 3, min.cutoff = 0, max.cutoff = 1, split.by = "orig.ident", order = T) + coord_equal()
FeaturePlot(spleen.combined.sct, features = "IL1A", pt.size = 0.2, ncol = 3, min.cutoff = 0, max.cutoff = 3, split.by = "orig.ident", order = T) + coord_equal()
FeaturePlot(spleen.combined.sct, features = "PTGS2", pt.size = 0.2, ncol = 3, min.cutoff = 0,max.cutoff = 1, split.by = "orig.ident", order = T) + coord_equal()
FeaturePlot(spleen.combined.sct, features = "sct_IL1B", pt.size = 0.2, ncol = 3, min.cutoff = 0,max.cutoff = 3, split.by = "orig.ident", order = T) + coord_equal()
FeaturePlot(spleen.combined.sct, features = "GZMK", pt.size = 0.2, ncol = 3, min.cutoff = 0, split.by = "orig.ident", order = T) + coord_equal()

FeaturePlot(spleen.combined.sct, features = "sct_IKZF2", pt.size = 0.2, ncol = 3, min.cutoff = 0, max.cutoff = 5, split.by = "orig.ident", order = T) + coord_equal()

FeaturePlot(spleen.combined.sct, features = "IL4", pt.size = 0.2, ncol = 3, min.cutoff = 0, max.cutoff = 1, split.by = "orig.ident", order = T) + coord_equal()
FeaturePlot(spleen.combined.sct, features = "IL13", pt.size = 0.2, ncol = 3, min.cutoff = 0, max.cutoff = 1, split.by = "orig.ident", order = T) + coord_equal()

FeaturePlot(spleen.combined.sct, features = "sct_CCL22", pt.size = 0.2, ncol = 3, min.cutoff = 0,max.cutoff = 4, split.by = "orig.ident", order = T) + coord_equal()

#没变???
FeaturePlot(spleen.combined.sct, features = "CXCL3", pt.size = 0.2, ncol = 3, min.cutoff = 0,max.cutoff = 3, split.by = "orig.ident", order = T) + coord_equal()
FeaturePlot(spleen.combined.sct, features = "CXCL4", pt.size = 0.2, ncol = 3, min.cutoff = 0,max.cutoff = 3, split.by = "orig.ident", order = T) + coord_equal()
FeaturePlot(spleen.combined.sct, features = "IL6", pt.size = 0.2, ncol = 3, min.cutoff = 0,max.cutoff = 3, split.by = "orig.ident", order = T) + coord_equal()
FeaturePlot(spleen.combined.sct, features = "TNF", pt.size = 0.2, ncol = 3, min.cutoff = 0,max.cutoff = 3, split.by = "orig.ident", order = T) + coord_equal()
FeaturePlot(spleen.combined.sct, features = "IL18", pt.size = 0.2, ncol = 3, min.cutoff = 0,max.cutoff = 3, split.by = "orig.ident", order = T) + coord_equal()
FeaturePlot(spleen.combined.sct, features = "FSCN1", pt.size = 0.2, ncol = 3, min.cutoff = 0,max.cutoff = 5, split.by = "orig.ident", order = T) + coord_equal()
kir_genes <- c("KIR3DL1", "KIR2DL3", "KIR2DL1", "KIR2DL4", "KIR3DL2", "KIR2DS4")
spleen.combined.sct[["KIR_score"]] <- Matrix::colMeans(GetAssayData(spleen.combined.sct, slot = "data")[kir_genes, ])
FeaturePlot(spleen.combined.sct, features = "KIR_score",pt.size = 0.2, ncol = 3, min.cutoff = 0,max.cutoff = 5, split.by = "orig.ident", order = T) + coord_equal()
FeaturePlot(spleen.combined.sct, features = "RORC", pt.size = 0.2, ncol = 3, min.cutoff = 0, max.cutoff = 3, split.by = "orig.ident", order = T) + coord_equal()

FeaturePlot(spleen.combined.sct, features = "sct_IL1B", pt.size = 0.2, ncol = 3, min.cutoff = 0,max.cutoff = 3, split.by = "orig.ident", order = T) + coord_equal()

VlnPlot(spleen.combined.sct, features =c("KIR3DL1", "KIR2DL3", "KIR2DL1", "KIR2DL4", "KIR3DL2", "KIR2DS4"), ncol = 3, pt.size = 0, group.by = "sample")
table(spleen.combined.sct$Hashtag,spleen.combined.sct$orig.ident)

#Hashtag three catch less monocytes, likely experiment variation
DimPlot(spleen.combined.sct, label = T) + NoLegend()

##Mark Cell Type
#Note：no ASC found in this

DefaultAssay(spleen.combined.sct)

FeaturePlot(spleen.combined.sct, features = "CD8A", pt.size = 0.2, ncol = 1, min.cutoff = 0, max.cutoff = 1, order = T, label = T)
FeaturePlot(spleen.combined.sct, features = "CD4", pt.size = 0.2, ncol = 1, min.cutoff = 0, max.cutoff = 1,  order = T, label = T)
FeaturePlot(spleen.combined.sct, features = "IL17F", pt.size = 0.2, ncol = 1, min.cutoff = 0, max.cutoff = 1, order = T, label = T)
FeaturePlot(spleen.combined.sct, features = "FOXP3", pt.size = 0.2, ncol = 1, min.cutoff = 0,  order = T, label = T)
FeaturePlot(spleen.combined.sct, features = "FCGR3A", pt.size = 0.2, ncol = 1, min.cutoff = 0, order = T, label = T)
FeaturePlot(spleen.combined.sct, features = "MS4A1", pt.size = 0.2, ncol = 1, min.cutoff = 0, order = T, label = T)

#make seurat_clusters larger cell type
#spleen.combined.sct$seurat_clusters <- spleen.combined.sct$integrated_snn_res.0.8
#Idents(spleen.combined.sct) <- spleen.combined.sct$seurat_clusters
Idents(spleen.combined.sct)
all.markers.spleen <- FindAllMarkers(spleen.combined.sct, only.pos = T)
write.csv(file = "./table/all.markers.spleen.csv", all.markers.spleen)
# 
# #B cell
# spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '2' = "B_Cell")
# spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '20' = "B_Cell")
# spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '0' = "B_Cell")
# spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '5' = "B_Cell")
# spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '24' = "B_Cell")
# spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '11' = "B_Cell")
# spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '3' = "B_Cell")
# spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '16' = "B_Cell")
# spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '7' = "B_Cell")
# spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '9' = "B_Cell")
# spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '25' = "B_Cell")
# 
# #Monocytes
# spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '1' = "Monocytes/DC")
# 
# #NK
# spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '14' = "NK")
# spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '10' = "NK")
# spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '17' = "NK")
# 
# #CD8
# spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '6' = "CD8")
# spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '4' = "CD8")
# spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '12' = "CD8")
# spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '19' = "CD8")
# spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '13' = "CD8")
# spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '23' = "CD8")
# 
# #CD4
# spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '8' = "CD4")
# spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '21' = "CD4")
# spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '15' = "Treg")
# 
# #not determined
# spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '18' = "NotDetermined")
# spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '22' = "NotDetermined")
# spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '26' = "NotDetermined")
# 

DimPlot(spleen.combined.sct, label = T) + coord_equal()
FeaturePlot(spleen.combined.sct, features = c("IFNG"),order=T, 
            split.by="orig.ident", pt.size = 0.2, min.cutoff = 0)+coord_equal()

#identical(pbmc.clean$seurat_clusters, pbmc.clean$integrated_snn_res.1.2)
Idents(spleen.combined.sct) <- spleen.combined.sct$seurat_clusters
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '0' = "X00_Mem_B")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '1' = "X01_Macrophage_DC")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '2' = "X02_NaiveB")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '3' = "X03_Mem_B")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '4' = "X04_Mem_B")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '5' = "X05_GZMK_CD8T")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '6' = "X06_Mem_B")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '7' = "X07_Naive_CD4T")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '8' = "X08_PreGC_lightZone_B")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '9' = "X09_IGKV3_20_Mem_B")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '10' = "X10_Rest_NK")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '11' = "X11_IGHV3_23_Mem_B")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '12' = "X12_GZMB_CD8T")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '13' = "X13_Th17_MAIT")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '14' = "X14_GC_DarkZone_B")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '15' = "X15_Tissue_Resident_CD8T")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '16' = "X16_Cytotoxic_regulatory_T_cells")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '17' = "X17_Activated_NK")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '18' = "X18_Plasmablast_Precursors_B")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '19' = "X19_Undetermined")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '20' = "X20_gamma_delta_T ")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '21' = "X21_Treg")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '22' = "X22_PlasmaBlast")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '23' = "X23_IGKV3_20_Naive_B")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '24' = "X24_Naive_CD8T")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '25' = "X25_IGHV3_23_Naive_B")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '26' = "X26_moDC_like_FSCN1+")

spleen.combined.sct$celltype <- Idents(spleen.combined.sct)
table(spleen.combined.sct$orig.ident)

spleen.combined.sct.16 <- subset(spleen.combined.sct, idents = "X16_Cytotoxic_regulatory_T_cells")
Idents(spleen.combined.sct.16) <- spleen.combined.sct.16$orig.ident

spleen.combined.sct.16 <- PrepSCTFindMarkers(spleen.combined.sct.16)
library(MAST)
markers_cluster16 <- FindMarkers(
  object = spleen.combined.sct.16,
  ident.1 = "spleen-control",
  ident.2 = "spleen_peptidoglycan",
  group.by = "orig.ident",
  test.use = "wilcox",
  assay = "SCT",
  slot = "data")

markers_cluster16$gene <- row.names(markers_cluster16)

DimPlot(spleen.combined.sct, label =T) + NoLegend() + coord_fixed()

library(EnhancedVolcano)

# Generate the volcano plot
EnhancedVolcano(markers_cluster16,
                lab = markers_cluster16$gene,
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,             # Adjust threshold as needed
                FCcutoff = 0.25,            # Log2 fold-change threshold (adjustable)
                pointSize = 2.5,
                labSize = 4.0,
                title = 'DEGs in Cluster 16 (Control vs Peptidoglycan)',
                subtitle = 'MAST, SCT assay',
                caption = 'Log2FC threshold = 0.25, adj.P < 0.05',
                legendPosition = 'right',
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.5,
                col = c('grey30', 'forestgreen', 'royalblue', 'red2')  # optional color scheme
)

meta <- spleen.combined.sct@meta.data

# Create a count table of celltypes by sample
celltype_counts <- table(meta$sample, meta$celltype)

# Convert to a data frame
celltype_df <- as.data.frame(celltype_counts)

# Rename for clarity
colnames(celltype_df) <- c("Sample", "CellType", "Count")

# Calculate percentages
celltype_percent <- celltype_df %>%
  group_by(Sample) %>%
  mutate(Percent = Count / sum(Count) * 100)

# View the result
head(celltype_percent)


celltype_percent <- celltype_percent %>%
  mutate(Group = ifelse(grepl("^control", Sample), "Control", "Peptidoglycan"))

ggplot(celltype_percent, aes(x = Group, y = Percent, fill = Group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.7) +
  facet_wrap(~ CellType, scales = "free_y", ncol = 6) +
  theme_minimal() +
  labs(title = "Cell Type Proportions by Condition",
       x = "Condition",
       y = "Cell Type Percentage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 8),
        legend.position = "none")
library(ggpubr)

p <- ggplot(celltype_percent, aes(x = Group, y = Percent, fill = Group)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.7) +
  facet_wrap(~ CellType, scales = "free_y", ncol = 6) +
  stat_compare_means(method = "t.test", label = "p.format") +  # Add p-values
  theme_minimal() +
  labs(title = "Cell Type Proportions by Condition",
       x = "Condition", y = "Percentage") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 8),
        legend.position = "none")

p

ggsave(filename = "./Result.Percent.pdf", p, width = 10, height = 15,dpi=300)

save(spleen.combined.sct, file = "./RObject/02.Celltype_labelled_UMAP.RData")


