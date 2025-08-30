#Single cell Lyme disease post merge
#Author: Xin Chen, Ph.D.
#Date Created: 022723
#Date Updated: 02/27/23
#After Find clusters On Cluster, run UMAP

#set up working directory
#setwd("/Users/xzhou7/Library/CloudStorage/Box-Box/Xin.Chen.Shareable/R7_NR")

setwd(dir = "~/Library/CloudStorage/Box-Box/Xin.Chen.Shareable/R7_NR")
setwd(dir = "C:/Users/zhoux/Box/Xin.Chen.Shareable/R7_NR/")
#load necessary package
library(data.table)
library(dplyr)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(stringr)
library(patchwork)
library(metap)
library(multtest)
library(Rcpp)
library(ggsci)
library(cowplot)
getwd()
load("./Data/022723_Step2.3_Subset5000_RunUMAP.RData")

DefaultAssay(pbmc.integrated.5000) <- "integrated"

pbmc.integrated.5000[[]]
p1 <- DimPlot(pbmc.integrated.5000, reduction = "umap", group.by = "seurat_clusters",label=T) + 
  ggtitle("UMAP of PBMC from Lyme Disease Patients/Controls")
p1
ggsave(filename = "./Results/UMAP.pdf", p1, width=8.5, height = 7, dpi = 300)

p1 <- DimPlot(pbmc.integrated.5000, reduction = "umap", split.by = "condition",label=T) + 
  ggtitle("UMAP of PBMC from Lyme Disease Patients/Controls")
p1
ggsave(filename = "./Results/UMAP.by condition B4 cleanup.pdf", p1, width=20, height = 7, dpi = 300)

#Check different conditions before clean up/by subject
count2.2 <- data.frame(table(pbmc.integrated.5000$subject, pbmc.integrated.5000$seurat_clusters))
total.bysubject <- data.frame(table(pbmc.integrated.5000$subject))
colnames(total.bysubject)[2] <- "Total"
count2.2$conditon <- "HD"
count2.2$conditon[str_detect(count2.2$Var1, "PTLD")] <- "PTLD"
count2.2$conditon[str_detect(count2.2$Var1, "PTLDN")] <- "PTLDN"
count2.2$conditon[str_detect(count2.2$Var1, "RTH")] <- "RTH"
count2.2 <- merge(count2.2, total.bysubject, by="Var1")
count2.2$Percent <- count2.2$Freq/count2.2$Total
p2.2 <- ggplot(count2.2, aes(x=conditon, y=Percent, color=conditon)) + geom_boxplot(alpha=0.9)
p2.2 <- p2.2 + facet_wrap(.~Var2, scales = "free") + theme_cowplot() + scale_color_aaas()
#p2 <- p2 + geom_jitter(color="black", size=0.4, alpha=0.9)
p2.2

ggsave(filename = "./Results/022823_cluster percentage between condition b4 cleanup.pdf", p2.2, width = 20, height = 10, dpi = 300)

pbmc.integrated.5000[[]]

#Check markers for CL32

Idents(pbmc.integrated.5000) <- pbmc.integrated.5000$integrated_snn_res.1.2

Idents(pbmc.integrated.5000) <- "integrated_snn_res.1.2"
marker1 <- FindMarkers(pbmc.integrated.5000,grouping.var="integrated_snn_res.1.2",ident.1='32',assay="integrated")
marker1

#Check antibody for CL32
DefaultAssay(pbmc.integrated.5000) <- "antibody"
antibody.list <- row.names(pbmc.integrated.5000)
for (i in 1:length(antibody.list)){
  marker <- antibody.list[i]
  print(marker)
  marker.plot <- FeaturePlot(pbmc.integrated.5000,features = marker,min.cutoff = 0, max.cutoff = 1000) + scale_colour_gradient(low = "#D3D3D3",high = "#8C1515")
  path <- paste0("./Markers/Before_cleanup/ANTIBODY.",marker,".pdf")
  ggsave(path, marker.plot,width = 4, height = 3.5, dpi = 300)
}

#Clean cluster
CL0 <- subset(pbmc.integrated.5000, seurat_clusters == 0)
cell.id.0 <- CellSelector(DimPlot(CL0))

CL1 <- subset(pbmc.integrated.5000, seurat_clusters == 1)
cell.id.1 <- CellSelector(DimPlot(CL1))

CL2 <- subset(pbmc.integrated.5000, seurat_clusters == 2)
cell.id.2 <- CellSelector(DimPlot(CL2))

CL3 <- subset(pbmc.integrated.5000, seurat_clusters == 3)
cell.id.3 <- CellSelector(DimPlot(CL3))

CL4 <- subset(pbmc.integrated.5000, seurat_clusters == 4)
cell.id.4 <- CellSelector(DimPlot(CL4))

CL5 <- subset(pbmc.integrated.5000, seurat_clusters == 5)
cell.id.5 <- CellSelector(DimPlot(CL5))

CL6 <- subset(pbmc.integrated.5000, seurat_clusters == 6)
cell.id.6 <- CellSelector(DimPlot(CL6))

CL7 <- subset(pbmc.integrated.5000, seurat_clusters == 7)
cell.id.7 <- CellSelector(DimPlot(CL7))

CL8 <- subset(pbmc.integrated.5000, seurat_clusters == 8)
cell.id.8 <- CellSelector(DimPlot(CL8))

CL9 <- subset(pbmc.integrated.5000, seurat_clusters == 9)
cell.id.9 <- CellSelector(DimPlot(CL9))

CL10 <- subset(pbmc.integrated.5000, seurat_clusters == 10)
cell.id.10 <- CellSelector(DimPlot(CL10))

CL11 <- subset(pbmc.integrated.5000, seurat_clusters == 11)
cell.id.11 <- CellSelector(DimPlot(CL11))

CL12 <- subset(pbmc.integrated.5000, seurat_clusters == 12)
cell.id.12 <- CellSelector(DimPlot(CL12))

CL13 <- subset(pbmc.integrated.5000, seurat_clusters == 13)
cell.id.13 <- CellSelector(DimPlot(CL13))

CL14 <- subset(pbmc.integrated.5000, seurat_clusters == 14)
cell.id.14 <- CellSelector(DimPlot(CL14))

CL15 <- subset(pbmc.integrated.5000, seurat_clusters == 15)
cell.id.15 <- CellSelector(DimPlot(CL15))

CL16 <- subset(pbmc.integrated.5000, seurat_clusters == 16)
cell.id.16 <- CellSelector(DimPlot(CL16))

CL17 <- subset(pbmc.integrated.5000, seurat_clusters == 17)
cell.id.17 <- CellSelector(DimPlot(CL17))

CL18 <- subset(pbmc.integrated.5000, seurat_clusters == 18)
cell.id.18 <- CellSelector(DimPlot(CL18))

CL19 <- subset(pbmc.integrated.5000, seurat_clusters == 19)
cell.id.19 <- CellSelector(DimPlot(CL19))

CL20 <- subset(pbmc.integrated.5000, seurat_clusters == 20)
cell.id.20 <- CellSelector(DimPlot(CL20))

CL21 <- subset(pbmc.integrated.5000, seurat_clusters == 21)
cell.id.21 <- CellSelector(DimPlot(CL21))

CL22 <- subset(pbmc.integrated.5000, seurat_clusters == 22)
cell.id.22 <- CellSelector(DimPlot(CL22))

CL23 <- subset(pbmc.integrated.5000, seurat_clusters == 23)
cell.id.23 <- CellSelector(DimPlot(CL23))

CL24 <- subset(pbmc.integrated.5000, seurat_clusters == 24)
cell.id.24 <- CellSelector(DimPlot(CL24))

CL25 <- subset(pbmc.integrated.5000, seurat_clusters == 25)
cell.id.25 <- CellSelector(DimPlot(CL25))

CL26 <- subset(pbmc.integrated.5000, seurat_clusters == 26)
cell.id.26 <- CellSelector(DimPlot(CL26))

CL27 <- subset(pbmc.integrated.5000, seurat_clusters == 27)
cell.id.27 <- CellSelector(DimPlot(CL27))

CL28 <- subset(pbmc.integrated.5000, seurat_clusters == 28)
cell.id.28 <- CellSelector(DimPlot(CL28))

CL29 <- subset(pbmc.integrated.5000, seurat_clusters == 29)
cell.id.29 <- CellSelector(DimPlot(CL29))

#Remove cluster after 30 dur to low cell number
table(pbmc.integrated.5000$seurat_clusters)

pbmc.integrated[[]]

subset.id <- c(cell.id.0,cell.id.1, cell.id.2, cell.id.3, cell.id.4, cell.id.5,
               cell.id.6, cell.id.7, cell.id.8, cell.id.9, cell.id.10, cell.id.11, 
               cell.id.12, cell.id.13, cell.id.14, cell.id.15, cell.id.16, cell.id.17, 
               cell.id.18, cell.id.19, cell.id.20, cell.id.21, cell.id.22, cell.id.23,
               cell.id.24, cell.id.25, cell.id.26, cell.id.27, cell.id.28, cell.id.29)

save(subset.id, file = "./Data/Step4_022723_subset.id.RData")

#upload to cluster:
#scp ./Desktop/Step4_022723_subset.id.RData xchen7@scg.stanford.edu:/labs/mmdavis/xx213/cluster_run/ 

module load R/4.0

setwd("/oak/stanford/scg/lab_mmdavis/xx213/cluster_run/")

#load necessary package
library(data.table)
library(Seurat)
library(ggplot2)
library(stringr)
library(patchwork)
library(dplyr)
library(metap)
library(multtest)
library(Rcpp)

getwd()
load("/oak/stanford/scg/lab_mmdavis/xx213/cluster_run/Step4_022723_subset.id.RData")
load("/oak/stanford/scg/lab_mmdavis/xx213/cluster_run/022123_Step2.3_Subset5000_RunUMAP.RData")

pbmc.clean <- subset(pbmc.integrated.5000, cells = subset.id)
save(pbmc.clean, file = "./Step4_022723_Clean_cluster.RData")

DimPlot(pbmc.clean, group.by = "seurat_clusters", label = T)
?DimPlot
pbmc.integrated[[]]

save(pbmc.clean, file = "./Step4_022723_Clean_cluster.RData")

#transfer back to desktop
DimPlot(pbmc.clean, group.by = "seurat_clusters", label = T)
?DimPlot
pbmc.integrated[[]]

save(pbmc.clean, file = "./Data/012022_R5_Step3.1_Clean_cluster.RData")
