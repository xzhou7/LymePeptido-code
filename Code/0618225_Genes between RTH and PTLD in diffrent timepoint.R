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
pbmc.clean <- RenameIdents(object = pbmc.clean, '3' = "X03_Treg")
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
p.label <- DimPlot(pbmc.clean,label = T,raster=FALSE) + coord_equal()
p.label

ggsave(filename = "./Results/061725Dimplot_label1.pdf",
       plot = p.label,
       width = 12, height = 12, dpi = 300)


FeaturePlot(pbmc.clean, features = "FOXP3", raster = FALSE, min.cutoff = 0) +
  scale_color_gradient(low = "#D3D3D3", high = "#8C1515") +
  coord_equal()

FeaturePlot(pbmc.clean, features = "RORC", raster = FALSE, min.cutoff = 0) +
  scale_color_gradient(low = "#D3D3D3", high = "#8C1515") +
  coord_equal()

###############X_03CD4 check marker#######
#Idents(pbmc.clean) <- "CellType.1"  # or "CellType.1C" if that’s where your labels are
#levels(pbmc.clean)

#markers_X03 <- FindMarkers(pbmc.clean, ident.1 = "X03_CD4T", only.pos = TRUE, 
#                           min.pct = 0.25, logfc.threshold = 0.25)

#head(markers_X03[order(-markers_X03$avg_log2FC), ])
######Check Tregs between HD and PTLDN
Idents(pbmc.clean) <- "CellType.1"
X03_CD4T <- subset(pbmc.clean, idents = "X03_CD4T")
table(X03_CD4T$condition)
X03_CD4T <- subset(X03_CD4T, condition %in% c("HD", "PTLDN"))
Idents(X03_CD4T) <- "condition"
de_HD_PTLDN <- FindMarkers(X03_CD4T, ident.1 = "PTLDN", ident.2 = "HD",
                           logfc.threshold = 0.1, min.pct = 0.1)
de_HD_PTLDN$gene <- rownames(de_HD_PTLDN)

library(ggplot2)

de_HD_PTLDN$significance <- "NS"
de_HD_PTLDN$significance[de_HD_PTLDN$avg_log2FC > 0.25 & de_HD_PTLDN$p_val_adj < 0.05] <- "Up"
de_HD_PTLDN$significance[de_HD_PTLDN$avg_log2FC < -0.25 & de_HD_PTLDN$p_val_adj < 0.05] <- "Down"

volcano_plot <- ggplot(de_HD_PTLDN, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significance)) +
  geom_point(alpha = 0.8, size = 1.2) +
  geom_text_repel(data = top_genes, aes(label = gene), size = 3, max.overlaps = Inf) +
  scale_color_manual(values = c("Up" = "#8C1515", "Down" = "#1f78b4", "NS" = "gray")) +
  theme_minimal(base_size = 12) +
  labs(
    x = "log2 Fold Change (PTLDN vs HD)",
    y = "-log10 adjusted p-value",
    title = "X03_CD4T: PTLDN vs HD"
  )

# Display
volcano_plot
