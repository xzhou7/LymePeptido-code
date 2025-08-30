#Single cell Lyme disease 
#Author: Xin Chen, Ph.D.
#Date Created: N022723
#After clean cluster, rename the cluster for each subsets
#Link: https://satijalab.org/seurat/articles/integration_introduction.html


#set up working directory
#setwd("/Users/xinchen/Library/CloudStorage/Box-Box/Xin.Chen.Shareable/R7_NR/")
setwd("~/Library/CloudStorage/Dropbox/lyme_disease/R7_NR/")
#setwd("/Users/xzhou7/Library/CloudStorage/Box-Box/Xin.Chen.Shareable/R7_NR/")

#setwd("C:/Users/zhoux/Box/Xin.Chen.Shareable/R7_NR/")
getwd()

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
library(cowplot)
source("./Code/00_colorKey_XC.R")
load("./Data/Step4_022723_Clean_cluster.RData")

pbmc.clean
DefaultAssay(pbmc.clean) <- "integrated"
rownames(pbmc.clean)[str_detect(rownames(pbmc.clean), "FCER1A")]
p1 <- DimPlot(pbmc.clean, reduction = "umap", label = T,raster=FALSE) + coord_equal()
p1
#ggsave(filename = "./Results/UMAP_clean_cluster.pdf", p1, width=8.5, height = 7, dpi = 300)

rownames(pbmc.clean)[str_detect(rownames(pbmc.clean), "CD14")]

all.list <- c("CD4", "CD8A", "NKG7", "CD1C","FCGR3A", "CD3E", "CCR7", "S100A4","TBX21","RORC", "IFNG", "IL7R", "NCAM1", "ITGAM",
              "CD38", "TIGIT","ITGAL", "SELL","CTLA4","ANXA5","CCL5","KLRB1","CASP1", "PIK3CA", "UTS2", "PDZK1IP1","TRDC",
              "GZMA", "GZMK","GZMB", "PRF1", "FOXP3", "IL2RA", "LAG3", "ISG15", "IFI6", "STAT1", "LY6E", "MX1", "MX2", "GATA3","FCER1G","FGFBP2","SPON2",
              "CLIC3","IGFBP7","MYOM2","LTB","ITGB1","KIR2DL4","CCL3","CD160","SPTSSB","CD3D","KLRC2","LAMP1","KIR2DL3")

Marker.list <- c("CD4", "CD8A", "NKG7", "CD1C","FCGR3A", "CD3E", "CCR7", "S100A4","TBX21","RORC", "IFNG", "IL7R", "NCAM1", "ITGAM",
                 "CD38", "TIGIT","ITGAL", "SELL","CTLA4","ANXA5","CCL5","KLRB1","TRDC",
                 "GZMA", "GZMK","GZMB", "PRF1", "FOXP3", "IL2RA",  "IFI6", "STAT1","MX1","GATA3","FCER1G","FGFBP2",
                 "ITGB1","KIR2DL1","CD160","CD3D","KLRC2", "FCER1A")

Marker.list[(! Marker.list %in% rownames(pbmc.clean))]

DefaultAssay(pbmc.clean) <- "integrated"
#Integrated Markers
for (i in 1:length(Marker.list)){
  marker <- Marker.list[i]
  print(marker)
  marker.plot <- FeaturePlot(pbmc.clean,features = marker,min.cutoff = 0) + scale_colour_gradient(low = "#D3D3D3",high = "#8C1515")
  path <- paste0("./Markers/Inte.",marker,".pdf")
  ggsave(path, marker.plot,width = 4, height = 3.5, dpi = 300)
}

Marker.list_SCT <- c("CD4", "CD8A", "NKG7", "CD14","CD1C","FCGR3A", "CD3E", "CCR7", "S100A4","TBX21","RORC", "IFNG", "IL7R", "NCAM1", "ITGAM",
                     "CD38", "TIGIT","ITGAL", "SELL","CTLA4","ANXA5","CCL5","KLRB1","TRDC",
                     "GZMA", "GZMK","GZMB", "PRF1", "FOXP3", "IL2RA",  "IFI6", "STAT1","MX1","GATA3","FCER1G","FGFBP2",
                     "ITGB1","KIR2DL1","CD160","CD3D","KLRC2","IGHA1-secreted")


DefaultAssay(pbmc.clean) <- "SCT"
#SCT Markers
for (i in 1:length(Marker.list_SCT)){
  marker <- Marker.list_SCT[i]
  print(marker)
  marker.plot <- FeaturePlot(pbmc.clean,features = marker,min.cutoff = 0) + scale_colour_gradient(low = "#D3D3D3",high = "#8C1515")
  path <- paste0("./Markers/SCT.",marker,".pdf")
  ggsave(path, marker.plot,width = 4, height = 3.5, dpi = 300)
}

DefaultAssay(pbmc.clean) <- "antibody"
antibody.list <- row.names(pbmc.clean)
for (i in 1:length(antibody.list)){
  marker <- antibody.list[i]
  print(marker)
  marker.plot <- FeaturePlot(pbmc.clean,features = marker,min.cutoff = 0, max.cutoff = 1000) + scale_colour_gradient(low = "#D3D3D3",high = "#8C1515")
  path <- paste0("./Markers/ANTIBODY.",marker,".pdf")
  ggsave(path, marker.plot,width = 4, height = 3.5, dpi = 300)
}

rownames(pbmc.clean)[str_detect(rownames(pbmc.clean), "IGHA")]
#plot single marker


#pbmc.clean$integrated_snn_res.1.5 <- factor
pbmc.clean[[]]

DimPlot(pbmc.clean,label = T,raster=FALSE)

#identical(pbmc.clean$seurat_clusters, pbmc.clean$integrated_snn_res.1.2)
Idents(pbmc.clean) <- pbmc.clean$integrated_snn_res.1.2
pbmc.clean <- RenameIdents(object = pbmc.clean, '0' = "X00_NavieCD4T")
pbmc.clean <- RenameIdents(object = pbmc.clean, '1' = "X01_NaiveCD4T")
pbmc.clean <- RenameIdents(object = pbmc.clean, '2' = "X02_CD14mono")
pbmc.clean <- RenameIdents(object = pbmc.clean, '3' = "X03_CD4T")
pbmc.clean <- RenameIdents(object = pbmc.clean, '4' = "X04_NK")
pbmc.clean <- RenameIdents(object = pbmc.clean, '5' = "X05_B")
pbmc.clean <- RenameIdents(object = pbmc.clean, '6' = "X06_CD14mono")
pbmc.clean <- RenameIdents(object = pbmc.clean, '7' = "X07_GZMB_CD8T")
pbmc.clean <- RenameIdents(object = pbmc.clean, '8' = "X08_CD4T")
pbmc.clean <- RenameIdents(object = pbmc.clean, '9' = "X09_CD14mono")
pbmc.clean <- RenameIdents(object = pbmc.clean, '10' = "X10_TH17_CD4T")
pbmc.clean <- RenameIdents(object = pbmc.clean, '11' = "X11_GZMK_CD8T")
pbmc.clean <- RenameIdents(object = pbmc.clean, '12' = "X12_NaiveCD8T")
pbmc.clean <- RenameIdents(object = pbmc.clean, '13' = "X13_IgD_B")
pbmc.clean <- RenameIdents(object = pbmc.clean, '14' = "X14_CD14mono")
pbmc.clean <- RenameIdents(object = pbmc.clean, '15' = "X15_CD4T")
pbmc.clean <- RenameIdents(object = pbmc.clean, '16' = "X16_CD14mono")
pbmc.clean <- RenameIdents(object = pbmc.clean, '17' = "X17_NK")
pbmc.clean <- RenameIdents(object = pbmc.clean, '18' = "X18_CD16mono")
pbmc.clean <- RenameIdents(object = pbmc.clean, '19' = "X19_myeloid")
pbmc.clean <- RenameIdents(object = pbmc.clean, '20' = "X20_CD4T")
pbmc.clean <- RenameIdents(object = pbmc.clean, '21' = "X21_DC")
pbmc.clean <- RenameIdents(object = pbmc.clean, '22' = "X22_IntermediateMono")
pbmc.clean <- RenameIdents(object = pbmc.clean, '23' = "X23_T")
pbmc.clean <- RenameIdents(object = pbmc.clean, '24' = "X24_CD14mono")
pbmc.clean <- RenameIdents(object = pbmc.clean, '25' = "X25_B")
pbmc.clean <- RenameIdents(object = pbmc.clean, '26' = "X26_NKT")
pbmc.clean <- RenameIdents(object = pbmc.clean, '27' = "X27_DR+GZMB+DC")
pbmc.clean <- RenameIdents(object = pbmc.clean, '28' = "X28_CD14mono")
pbmc.clean <- RenameIdents(object = pbmc.clean, '29' = "X29_NKT")

pbmc.clean$CellType.1 <- Idents(pbmc.clean)
pbmc.clean$CellType.1 <- factor(pbmc.clean$CellType.1, levels = sort(unique(pbmc.clean$CellType.1),decreasing=T))
Idents(pbmc.clean) <- pbmc.clean$CellType.1
p.label <- DimPlot(pbmc.clean,label = F,raster=FALSE) + coord_equal()
p.label
#ggsave(filename = "./Results/UMAP.LABEL .pdf", p.label, width = 18, height = 12, dpi=300)
#ggsave(filename = "../lyme_disease/Manuscript/Figures/FigureS1/S1A.pdf", p.label, width = 8, height = 5, dpi = 300)

Idents(pbmc.clean) <- pbmc.clean$integrated_snn_res.1.2
Idents(pbmc.clean) <- pbmc.clean$CellType.1
DimPlot(pbmc.clean, label=T)
pbmc.integrated <- RenameIdents(object = pbmc.integrated, "X22_IgA.B" = "X22_IgA.B_New")
Idents(pbmc.clean) <- "integrated_snn_res.1.2"
#marker1 <- FindMarkers(pbmc.clean,grouping.var="integrated_snn_res.1.2",ident.1='20',assay="integrated")
#marker1

#Highlight specific cluster (need to start from the first of the three line, do not excute from middle)
pbmc.pbmc.clean$Target.Cluster <- "others"
pbmc.pbmc.clean$Target.Cluster[pbmc.pbmc.clean$integrated_snn_res.1.5 == "17" ] <- "17"
DimPlot(pbmc.pbmc.clean, group.by = "Target.Cluster") + scale_color_manual(values=c("#8C1515", "#D3D3D3"))

# find markers for every cluster compared to all remaining cells, report only the positive
# ones
# pbmc.clean
# DefaultAssay(pbmc.clean) <- "integrated"
# pbmc.markers <- FindAllMarkers(pbmc.clean, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# pbmc.markers %>%
#   group_by(cluster) %>%
#   slice_max(n = 10, order_by = avg_log2FC)
#write.csv(pbmc.markers,"./Results/T cell gene markers.csv", row.names = FALSE)
pbmc.markers <- read.csv("./Results/T cell gene markers.csv")

rownames(pbmc.integrated)[str_detect(rownames(pbmc.integrated), "S100A4")]
FeaturePlot(pbmc.integrated, "GZMB",min.cutoff = 0)
FeaturePlot(pbmc.integrated, "CD1C",min.cutoff = 0)

DimPlot(pbmc.clean,split.by = "condition",raster=FALSE)

table(pbmc.clean$condition)

save(pbmc.clean, file = "./Data/Step4.1_022823_rename_clean_cluster.RData")

#Check percentage cluster between condition after clean up
setwd(dir = "~/Library/CloudStorage/Dropbox/lyme_disease/R7_NR/")
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
library(ggrepel)

getwd()
load("./Data/Step4.1_022823_rename_clean_cluster.RData")

#by subject
count2.2 <- data.frame(table(pbmc.clean$subject, pbmc.clean$CellType.1))
total.bysubject <- data.frame(table(pbmc.clean$subject))
colnames(total.bysubject)[2] <- "Total"
count2.2$conditon <- "HD"
count2.2$conditon[str_detect(count2.2$Var1, "PTLD")] <- "PTLD"
count2.2$conditon[str_detect(count2.2$Var1, "PTLDN")] <- "PTLDN"
count2.2$conditon[str_detect(count2.2$Var1, "RTH")] <- "RTH"
count2.2 <- merge(count2.2, total.bysubject, by="Var1")
count2.2$Percent <- count2.2$Freq/count2.2$Total

count2.2$SLICE <- "SLICE_I"
count2.2$SLICE[str_detect(count2.2$conditon, "PTLDN")] <- "SLICE_III"
count2.2$SLICE[str_detect(count2.2$conditon, "HD")] <- "SLICE_III"

count2.2$conditon <- factor(count2.2$conditon, levels = c("HD", "PTLDN", "RTH", "PTLD"))

comparasions <- list(c("HD", "PTLDN"), c("RTH","PTLD"))

p2.2 <- ggplot(count2.2, aes(x=conditon, y=Percent, color=conditon)) + geom_boxplot(alpha=0.9)
p2.2 <- p2.2 + facet_wrap(.~Var2, scales = "free") + theme_cowplot() + scale_color_aaas()
#p2 <- p2 + geom_jitter(color="black", size=0.4, alpha=0.9)
p2.2 <- p2.2 + stat_compare_means(method="t.test", comparisons = comparasions)
p2.2

#ggsave(filename = "./Results/061625_cluster percentage between condition after cleanup.pdf", p2.2, width = 20, height = 20, dpi = 300)


celltype.list <- unique(count2.2$Var2) %>% as.character()
 
# comparasions <- list(c("HD", "PTLDN"), c("RTH","PTLD"))
# for (i in 1:30){
#   print(celltype.list[i])
#   p2.3 <- subset(count2.2, Var2 %in% celltype.list[i]) %>% ggplot(aes(x=conditon, y=Percent, color=conditon)) + geom_boxplot(alpha=0.9)+stat_compare_means(method="t.test", comparisons = comparasions)
#   p2.3
#   path = paste0("../lyme_disease/Manuscript/Figures/FigureS1/Proportion_Whole/", celltype.list[i], ".pdf")
#   ggsave(filename = path, p2.3, width = 5, height = 5, dpi=300)
# }

#for single comparasion
celltype.list[1]
comparasions <- list(c("HD", "PTLDN"), c("RTH","PTLD"))
p2.3 <- subset(count2.2, Var2 %in% celltype.list[1]) %>% ggplot(aes(x=conditon, y=Percent, color=conditon)) + geom_boxplot(alpha=0.9)+stat_compare_means(method="t.test", comparisons = comparasions)
p2.3

##########For Figure 1C
Idents(pbmc.clean) <- pbmc.clean$integrated_snn_res.1.2
pbmc.clean <- RenameIdents(object = pbmc.clean, '0' = "NaiveCD4T")
pbmc.clean <- RenameIdents(object = pbmc.clean, '1' = "NaiveCD4T")
pbmc.clean <- RenameIdents(object = pbmc.clean, '2' = "Monocytes")
pbmc.clean <- RenameIdents(object = pbmc.clean, '3' = "CD4T")
pbmc.clean <- RenameIdents(object = pbmc.clean, '4' = "NK")
pbmc.clean <- RenameIdents(object = pbmc.clean, '5' = "B")
pbmc.clean <- RenameIdents(object = pbmc.clean, '6' = "Monocytes")
pbmc.clean <- RenameIdents(object = pbmc.clean, '7' = "CD8T")
pbmc.clean <- RenameIdents(object = pbmc.clean, '8' = "CD4T")
pbmc.clean <- RenameIdents(object = pbmc.clean, '9' = "Monocytes")
pbmc.clean <- RenameIdents(object = pbmc.clean, '10' = "CD4T")
pbmc.clean <- RenameIdents(object = pbmc.clean, '11' = "CD8T")
pbmc.clean <- RenameIdents(object = pbmc.clean, '12' = "NaiveCD8T")
pbmc.clean <- RenameIdents(object = pbmc.clean, '13' = "B")
pbmc.clean <- RenameIdents(object = pbmc.clean, '14' = "Monocytes")
pbmc.clean <- RenameIdents(object = pbmc.clean, '15' = "CD4T")
pbmc.clean <- RenameIdents(object = pbmc.clean, '16' = "Monocytes")
pbmc.clean <- RenameIdents(object = pbmc.clean, '17' = "NK")
pbmc.clean <- RenameIdents(object = pbmc.clean, '18' = "Monocytes")
pbmc.clean <- RenameIdents(object = pbmc.clean, '19' = "Monocytes")
pbmc.clean <- RenameIdents(object = pbmc.clean, '20' = "CD4T")
pbmc.clean <- RenameIdents(object = pbmc.clean, '21' = "DC")
pbmc.clean <- RenameIdents(object = pbmc.clean, '22' = "Monocytes")
pbmc.clean <- RenameIdents(object = pbmc.clean, '23' = "CD8T")
pbmc.clean <- RenameIdents(object = pbmc.clean, '24' = "Monocytes")
pbmc.clean <- RenameIdents(object = pbmc.clean, '25' = "B")
pbmc.clean <- RenameIdents(object = pbmc.clean, '26' = "CD8T")
pbmc.clean <- RenameIdents(object = pbmc.clean, '27' = "DC")
pbmc.clean <- RenameIdents(object = pbmc.clean, '28' = "Monocytes")
pbmc.clean <- RenameIdents(object = pbmc.clean, '29' = "NK")

pbmc.clean$CellType.1C <- Idents(pbmc.clean)
pbmc.clean$CellType.1C <- factor(pbmc.clean$CellType.1C, levels = sort(unique(pbmc.clean$CellType.1C),decreasing=T))
Idents(pbmc.clean) <- pbmc.clean$CellType.1C
p.label.1C <- DimPlot(pbmc.clean,label = T,raster=FALSE) + coord_equal() + scale_color_manual(values=my_palette)
p.label.1C

#ggsave(filename = "../lyme_disease/Manuscript/Figures/Figure1/Figure1B.pdf",p.label.1C, width = 6, height = 5, dpi = 300)
DotPlot(pbmc.clean, features = c("CD4", "CD8A", "NKG7", "CD14","CD1C","FCGR3A", "CD3E", "CCR7", "S100A4","TBX21","RORC", "IFNG", "IL7R", "NCAM1", "ITGAM",
                                 "CD38", "TIGIT","ITGAL", "SELL","CTLA4","ANXA5","CCL5","KLRB1","TRDC",
                                 "GZMA", "GZMK","GZMB", "PRF1", "FOXP3", "IL2RA",  "IFI6", "STAT1","MX1","GATA3","FCER1G","FGFBP2",
                                 "ITGB1","KIR2DL1","CD160","CD3D","KLRC2","IGHA1-secreted"))  
rownames(pbmc.clean)


pbmc.clean$CellType.1C <- factor(pbmc.clean$CellType.1C, levels = c("NaiveCD4T", "CD4T", "NaiveCD8T", "CD8T","NK", "DC", "Monocytes",  "B"))

Idents(pbmc.clean)

cells_to_keep <- sample(colnames(pbmc.clean), 30000)
pbmc.subset <- subset(pbmc.clean, cells = cells_to_keep)
p.label.1C <- DimPlot(pbmc.subset,label = F,raster=FALSE) + coord_equal()# + scale_color_manual(values=my_palette) + scale_fill_manual(values=my_palette_alpha_0.7)
p.label.1C
celltype_position_1=pbmc.clean@reductions$umap@cell.embeddings %>% as.data.frame() %>%
  cbind(celltype=pbmc.clean@meta.data$CellType.1C) %>%
  group_by(celltype) %>% dplyr::summarise(UMAP_1 = median(UMAP_1),UMAP_2 = median(UMAP_2))

p.label.1Cp2 <- p.label.1C+ geom_label_repel(data = celltype_position_1,aes(x=UMAP_1, y=UMAP_2, label=celltype, color=celltype),
                                    fontface="bold", point.padding=unit(0.1, "lines"), alpha=0.95)+theme(legend.position = "none") +
  scale_color_manual(values = my_palette)
p.label.1Cp2

#ggsave(filename = "../lyme_disease/Manuscript/Figures/Figure1/Figure1B_30Kcell.pdf",p.label.1Cp2, width = 6, height = 5, dpi = 300)

Idents(pbmc.clean) <- "CellType.1C"
#dp.1C <- DotPlot(pbmc.clean, features = c("CD4", "CCR7","CD8A","CD3E","GZMK","GZMB","NKG7","NCAM1", "CD1C", "FCGR3A", 
#                                 "FCER1G", "PECAM1","ICAM1", "CD83", "MS4A1", "CD79A"))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), legend.position = "top")
dp.1C <- DotPlot(pbmc.clean, features = c("CD4", "CCR7","CD8A","CD3E","GZMK","GZMB","NKG7","NCAM1", "CD1C", "FCGR3A", 
                                          "FCER1G", "PECAM1","ICAM1", "CD83", "MS4A1", "CD79A"))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

dp.1C <- dp.1C + scale_color_gradient2(low = "#D3D3D3", high = "#8C1515")
dp.1C
#ggsave(filename = "../Manuscript/Figures/Figure1-cohorts and PBMCs overview/Figure1C_v2.pdf", dp.1C, width = 7, height = 4, dpi = 300)

count2.3 <- data.frame(table(pbmc.clean$subject, pbmc.clean$CellType.1C))
total.bysubject <- data.frame(table(pbmc.clean$subject))
colnames(total.bysubject)[2] <- "Total"
count2.3$conditon <- "HD"
count2.3$conditon[str_detect(count2.3$Var1, "PTLD")] <- "PTLD"
count2.3$conditon[str_detect(count2.3$Var1, "PTLDN")] <- "PTLDN"
count2.3$conditon[str_detect(count2.3$Var1, "RTH")] <- "RTH"
count2.3 <- merge(count2.3, total.bysubject, by="Var1")
count2.3$Percent <- count2.3$Freq/count2.3$Total

count2.3$SLICE <- "SLICE_I"
count2.3$SLICE[str_detect(count2.3$conditon, "PTLDN")] <- "SLICE_III"
count2.3$SLICE[str_detect(count2.3$conditon, "HD")] <- "SLICE_III"

count2.3$conditon <- factor(count2.3$conditon, levels = c("HD", "PTLDN", "RTH", "PTLD"))

count2.3 <- filter(count2.3, ! Var1 %in% c("PTLD2V1","RTH1V7b5","RTH2V5b5","RTH5V3b2"))     

p2.3 <- ggplot(count2.3, aes(x=conditon, y=Percent, color=conditon)) + geom_boxplot(alpha=0.9)
p2.3 <- p2.3 + facet_wrap(.~Var2, scales = "free") + theme_cowplot() + scale_color_aaas()
#p2 <- p2 + geom_jitter(color="black", size=0.4, alpha=0.9)
p2.3 <- p2.3 + stat_compare_means(comparisons = comparasions)
p2.3

#######modify##########070725
# Compute CD8ratio
CD8ratio <- cbind(count2.3 %>% filter(Var2 == "CD8T"), count2.3 %>% filter(Var2 == "NaiveCD8T"))
colnames(CD8ratio)[1:7] <- paste(colnames(CD8ratio)[1:7], "1", sep = "_")
CD8ratio$CD8ratio <- CD8ratio$Percent_1 / CD8ratio$Percent

# Plot with jittered points and boxplot
pCD8ratio <- ggplot(CD8ratio, aes(x = conditon, y = CD8ratio, color = conditon)) +
  geom_jitter(width = 0.2, size = 2) +   # Add jittered data points
  geom_boxplot(alpha = 0.7, outlier.alpha = 0, width = 0.6) +  # Match Plot 2 style
  stat_compare_means(comparisons = comparasions) +  # Add p-values
  scale_y_log10() +  # Log-scale Y axis
  scale_color_manual(values = condition_color) +  # Match color scheme
  theme_cowplot(font_size = 14) +
  XZ_flip_x_label()  # Custom function for rotated x-axis labels

# Display plot
pCD8ratio

# Save plot
ggsave(
  filename = "../Manuscript/Figures/Figure1-cohorts and PBMCs overview/Figure1EA_v2.pdf",
  plot = pCD8ratio,
  width = 4,
  height = 5,
  dpi = 300
)

##################CD4####################
# Compute CD4ratio
CD4ratio <- cbind(count2.3 %>% filter(Var2 == "CD4T"), count2.3 %>% filter(Var2 == "NaiveCD4T"))
colnames(CD4ratio)[1:7] <- paste(colnames(CD4ratio)[1:7], "1", sep = "_")
CD4ratio$CD4ratio <- CD4ratio$Percent_1 / CD4ratio$Percent

# Plot with jittered points and boxplot
pCD4ratio <- ggplot(CD4ratio, aes(x = conditon, y = CD4ratio, color = conditon)) +
  geom_jitter(width = 0.2, size = 2) +   # Jittered dots
  geom_boxplot(alpha = 0.7, outlier.alpha = 0, width = 0.6) +  # Match Plot 1 style
  stat_compare_means(comparisons = comparasions) +  # p-values
  scale_y_log10() +  # Log y-axis
  scale_color_manual(values = condition_color) +
  theme_cowplot(font_size = 14) +
  XZ_flip_x_label()

# Display plot
pCD4ratio

# Save plot
ggsave(
  filename = "../Manuscript/Figures/Figure1-cohorts and PBMCs overview/Figure1EB_v2.pdf",
  plot = pCD4ratio,
  width = 4,
  height = 5,
  dpi = 300
)

#######

# Load cowplot (if not already loaded)
library(cowplot)
# Remove legend from CD8 plot only
pCD8ratio_nolegend <- pCD8ratio + theme(legend.position = "none") + ggtitle("CD8 Memory/Naïve Ratio")
pCD4ratio <- pCD4ratio + ggtitle("CD4 Memory/Naïve Ratio")

# Combine plots
combined_plot <- plot_grid(
  pCD8ratio_nolegend,
  pCD4ratio,
  labels = c("A", "B"),
  ncol = 2,
  align = "hv",
  axis = "tb"
)

# Display plot
combined_plot

# Save to file
ggsave(
  filename = "../Manuscript/Figures/Figure1-cohorts and PBMCs overview/Figure1E_combined_v2.pdf",
  plot = combined_plot,
  width = 6.5,
  height = 5,
  dpi = 300
)







