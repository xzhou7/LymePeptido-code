#Last Updated 04/09/2024

library(Seurat)
library(MAST)
library(monocle3)
library(dplyr)
library(stringr)
library(Scillus)
library(ggplot2)
library(ComplexHeatmap)
library(ggrepel)
library(ggpubr)
library(DESeq2)
library(EnhancedVolcano)
library(cowplot)
library(reshape2)

#mac
setwd("~/Library/CloudStorage/Box-Box/Xin.Chen.Shareable/R7_NR/")

#windows
setwd("D:/Dropbox/lyme_disease/R7_NR/")
getwd()

source("./Code/xztools.R")
source("./Code/00_colorKey.R")
patient_data <- read.csv("./Results/patient_data_addPTLDN.csv", header = T)

# #initial analysis (note that reanalysis is below after line 200)
# load("./Data/Step5.1_032723_T_NK_clean.RData")
# XXX.clean <- T_NK.clean
# DefaultAssay(XXX.clean) <- "integrated"
# XXX.clean <- FindVariableFeatures(XXX.clean, selection.method = "vst", nfeatures = 350)
# all.genes <- rownames(XXX.clean)
# XXX.clean <- ScaleData(XXX.clean, features = all.genes)
# DimPlot(XXX.clean, reduction = "pca")
# #JackStraw(XXX.clean, num.replicate = 100)
# ElbowPlot(XXX.clean)
# XXX.clean <- FindNeighbors(XXX.clean, dims = 1:15)
# XXX.clean <- FindClusters(XXX.clean, resolution = 0.7)
# XXX.clean <- RunUMAP(XXX.clean, dims = 1:15)
# 
# DimPlot(XXX.clean, raster=FALSE)
# 
# all.genes[str_detect(all.genes, "KIR")]
# 
# DefaultAssay(XXX.clean) <- "integrated"
# FeaturePlot(XXX.clean,"CD8A", min.cutoff = 0, order = T,raster=FALSE)
# FeaturePlot(XXX.clean,"CD4", min.cutoff = 0, order = T ,raster=FALSE)
# 
# FeaturePlot(XXX.clean,"CD37", min.cutoff = 0, order = T,raster=FALSE)
# FeaturePlot(XXX.clean,"KLRC1", min.cutoff = 0, order = T ,raster=FALSE)
# FeaturePlot(XXX.clean,"CD160", min.cutoff = 0, order = T ,raster=FALSE)
# FeaturePlot(XXX.clean,"CCR6", min.cutoff = 0, order = T ,raster=FALSE)
# FeaturePlot(XXX.clean,"CXCR4", min.cutoff = 0, order = T ,raster=FALSE)
# FeaturePlot(XXX.clean,"FCN1", min.cutoff = 0, order = T ,raster=FALSE)
# FeaturePlot(XXX.clean,"IFIT3", min.cutoff = 0, order = T ,raster=FALSE)     
# 
# FeaturePlot(XXX.clean,"KIR2DL1", min.cutoff = 0, order = T ,raster=FALSE)
# FeaturePlot(XXX.clean,"NCAM1", min.cutoff = 0, order = T ,raster=FALSE)
# FeaturePlot(XXX.clean,"GZMH", min.cutoff = 0, order = T ,raster=FALSE)
# 
# FeaturePlot(XXX.clean,"GZMB", min.cutoff = 0, order = T ,raster=FALSE)
# FeaturePlot(XXX.clean,"GZMK", min.cutoff = 0, order = T ,raster=FALSE)
# FeaturePlot(XXX.clean,"KLRB1", min.cutoff = 0, order = T, raster=FALSE)
# FeaturePlot(XXX.clean,"SELL", min.cutoff = 0, order = T, raster=FALSE)
# 
# FeaturePlot(XXX.clean,"PRF1", min.cutoff = 0, order = T, raster=FALSE)
# FeaturePlot(XXX.clean,"LEF1", min.cutoff = 0, order = T, raster=FALSE)
# FeaturePlot(XXX.clean,"LEF1", min.cutoff = 0, order = T, raster=FALSE)
# 
# FeaturePlot(XXX.clean,"FOXP3", min.cutoff = 0, order = T, raster=FALSE)
# FeaturePlot(XXX.clean,"RORC", min.cutoff = 0, order = T, raster=FALSE)
# FeaturePlot(XXX.clean,"TBX21", min.cutoff = 0, order = T, raster=FALSE)
# FeaturePlot(XXX.clean,"GATA3", min.cutoff = 0, order = T, raster=FALSE)
# FeaturePlot(XXX.clean,"TRDC", min.cutoff = 0, order = T, raster=FALSE)
# 
# DefaultAssay(XXX.clean) <- "antibody"
# rownames(XXX.clean)
# 
# FeaturePlot(XXX.clean,"CCR7.CCR7.AHS0273.pAbO", min.cutoff = 1, max.cutoff = 500, order = T, raster=FALSE)
# RidgePlot(XXX.clean,"CCR7.CCR7.AHS0273.pAbO", group.by = "orig.ident") & xlim (0,15)
# 
# FeaturePlot(XXX.clean,"CD45RA.HI100.PTPRC.AHS0009.pAbO", min.cutoff = 1, max.cutoff = 500, order = T, raster=FALSE)
# FeaturePlot(XXX.clean,"TCR.gamma-delta.11F2.TRG-TRD.AHS0142.pAbO", min.cutoff = 1, max.cutoff = 500, order = T, raster=FALSE)
# 
# # All.Markers <- FindAllMarkers(XXX.clean, logfc.threshold = 0.25, only.pos = T,test.use = "MAST")
# # All.Markers
# # write.csv(file = "./Results/All_Markers.T.NK.csv",All.Markers)
# 
# All.Markers <- read.csv(file = "./Results/All_Markers.T.NK.csv", header = T)
# 
# #filter(All.Markers, cluster == 12)
# 
# # heatmap.genemarker <- plot_heatmap(dataset = XXX.clean, 
# #              markers = All.Markers,
# #              sort_var = c("seurat_clusters","condition"),
# #              anno_var = c("seurat_clusters","condition"),
# #              anno_colors = list("Set2",                                             # RColorBrewer palette
# #                                 c("red","orange","yellow","purple") # color vector
# #                                ))
# # 
# # pdf(file="./Results/heatmap.T.marker.pdf",width = 7, height = 15)
# # draw(heatmap.genemarker)
# # dev.off()
# 
# all.genes
# 
# FindMarkers(XXX.clean, 2, 0, logfc.threshold = 0.25, only.pos = F,test.use = "MAST")
# FindMarkers(XXX.clean, 2, 3, logfc.threshold = 0.25, only.pos = F,test.use = "MAST")
# FindMarkers(XXX.clean, 12, 14, logfc.threshold = 0.25, only.pos = F,test.use = "MAST")
# 
# DimPlot(XXX.clean,label = T, raster=FALSE) & coord_equal()
# 
# XXX.clean$celltype <- "Other"
# XXX.clean$celltype[XXX.clean$integrated_snn_res.0.7 == "0"] <- "Naive_CD4"
# XXX.clean$celltype[XXX.clean$integrated_snn_res.0.7 == "1"] <- "CD56dim_NK"
# XXX.clean$celltype[XXX.clean$integrated_snn_res.0.7 == "2"] <- "Activated_CD4"
# XXX.clean$celltype[XXX.clean$integrated_snn_res.0.7 == "3"] <- "Activated_CD4_2"
# XXX.clean$celltype[XXX.clean$integrated_snn_res.0.7 == "4"] <- "GZMK_CD8"
# XXX.clean$celltype[XXX.clean$integrated_snn_res.0.7 == "5"] <- "Naive_CD8"
# XXX.clean$celltype[XXX.clean$integrated_snn_res.0.7 == "6"] <- "MAIT_T"
# XXX.clean$celltype[XXX.clean$integrated_snn_res.0.7 == "7"] <- "GZMB_CD8"
# XXX.clean$celltype[XXX.clean$integrated_snn_res.0.7 == "8"] <- "GammaDelta_T"
# XXX.clean$celltype[XXX.clean$integrated_snn_res.0.7 == "9"] <- "CD56bright_NK"
# XXX.clean$celltype[XXX.clean$integrated_snn_res.0.7 == "10"] <- "TFH_CD4"
# XXX.clean$celltype[XXX.clean$integrated_snn_res.0.7 == "11"] <- "Treg"
# XXX.clean$celltype[XXX.clean$integrated_snn_res.0.7 == "12"] <- "LEF1_NK"
# XXX.clean$celltype[XXX.clean$integrated_snn_res.0.7 == "13"] <- "Proliferating_Cell"
# XXX.clean$celltype[XXX.clean$integrated_snn_res.0.7 == "14"] <- "IFN_NK"
# 
# # Cluster 0
# CL0 <- subset(XXX.clean, integrated_snn_res.0.7 == 0)
# cell.id.0 <- CellSelector(DimPlot(CL0))
# 
# # Cluster 1
# CL1 <- subset(XXX.clean, integrated_snn_res.0.7 == 1)
# cell.id.1 <- CellSelector(DimPlot(CL1))
# 
# # Cluster 2
# CL2 <- subset(XXX.clean, integrated_snn_res.0.7 == 2)
# cell.id.2 <- CellSelector(DimPlot(CL2))
# 
# # Cluster 3
# CL3 <- subset(XXX.clean, integrated_snn_res.0.7 == 3)
# cell.id.3 <- CellSelector(DimPlot(CL3))
# 
# # Cluster 4
# CL4 <- subset(XXX.clean, integrated_snn_res.0.7 == 4)
# cell.id.4 <- CellSelector(DimPlot(CL4))
# 
# # Cluster 5
# CL5 <- subset(XXX.clean, integrated_snn_res.0.7 == 5)
# cell.id.5 <- CellSelector(DimPlot(CL5))
# 
# # Cluster 6
# CL6 <- subset(XXX.clean, integrated_snn_res.0.7 == 6)
# cell.id.6 <- CellSelector(DimPlot(CL6))
# 
# # Cluster 7
# CL7 <- subset(XXX.clean, integrated_snn_res.0.7 == 7)
# cell.id.7 <- CellSelector(DimPlot(CL7))
# 
# # CLueter 8
# CL8 <- subset(XXX.clean, integrated_snn_res.0.7 == 8)
# cell.id.8 <- CellSelector(DimPlot(CL8))
# 
# # Cluster 9
# CL9 <- subset(XXX.clean, integrated_snn_res.0.7 == 9)
# cell.id.9 <- CellSelector(DimPlot(CL9))
# 
# # Cluster 10
# CL10 <- subset(XXX.clean, integrated_snn_res.0.7 == 10)
# cell.id.10 <- CellSelector(DimPlot(CL10))
# 
# # Cluster 11
# CL11 <- subset(XXX.clean, integrated_snn_res.0.7 == 11)
# cell.id.11 <- CellSelector(DimPlot(CL11))
# 
# # Cluster 12
# CL12 <- subset(XXX.clean, integrated_snn_res.0.7 == 12)
# cell.id.12 <- CellSelector(DimPlot(CL12))
# 
# # Cluster 13
# CL13 <- subset(XXX.clean, integrated_snn_res.0.7 == 13)
# cell.id.13 <- CellSelector(DimPlot(CL13))
# 
# # Cluster 14
# CL14 <- subset(XXX.clean, integrated_snn_res.0.7 == 14)
# cell.id.14 <- CellSelector(DimPlot(CL14))
# 
# # subset.id <- c(cell.id.0, cell.id.1, cell.id.2, cell.id.3, cell.id.4, cell.id.5,
# #                cell.id.6, cell.id.7, cell.id.8, cell.id.9, cell.id.10, cell.id.11,
# #                cell.id.12, cell.id.13, cell.id.14)
# # save(subset.id, file = "./Data/Step4_070523_NK_T_subset.id.RData")
# load("./Data/Step4_070523_NK_T_subset.id.RData")
# 
# T.NK.clean_Annotated <- subset(XXX.clean, cells = subset.id)
#save(T.NK.clean_Annotated, file = "./Data/Step5_070523_T.NK.clean_Annotated.RData")

#################################################################################################
#reanalysis start here
#################################################################################################
load("./Data/Step5_070523_T.NK.clean_Annotated.RData")
T.NK.clean_Annotated$celltype[T.NK.clean_Annotated$integrated_snn_res.0.7 == "8"] <- "GammaDelta_T"
T.NK.clean_Annotated$celltype[T.NK.clean_Annotated$integrated_snn_res.0.7 == "10"] <- "TFH_CD4"
T.NK.clean_Annotated$celltype[T.NK.clean_Annotated$integrated_snn_res.0.7 == "13"] <- "Proliferating_Cell"

FindMarkers(T.NK.clean_Annotated, ident.1 = "IFN_CD4", ident.2 = "Activated_CD4", only.pos = T)

Idents(T.NK.clean_Annotated) <- "celltype"

T.NK.clean_Annotated_keep <- sample(colnames(T.NK.clean_Annotated), 40000)
T.NK.clean_Annotated.subset <- subset(T.NK.clean_Annotated, cells = T.NK.clean_Annotated_keep)

p_umap <- DimPlot(T.NK.clean_Annotated.subset, label =F, raster=FALSE) & coord_equal() 
#p_umap <- p_umap + scale_color_manual(values = my_palette) + scale_fill_manual(values = my_palette_alpha_0.7)
p_umap

# #annotate label position
celltype_position=T.NK.clean_Annotated@reductions$umap@cell.embeddings %>% as.data.frame() %>%
  cbind(celltype=T.NK.clean_Annotated@meta.data$celltype) %>%
  group_by(celltype) %>% dplyr::summarise(UMAP_1 = median(UMAP_1),UMAP_2 = median(UMAP_2))

p_umap2 <- p_umap+ geom_label_repel(data = celltype_position,aes(x=UMAP_1, y=UMAP_2, label=celltype, color=celltype),
                   fontface="bold", point.padding=unit(0.1, "lines"), alpha=0.95)+theme(legend.position = "none") +
  scale_color_manual(values = T.NK.celltype_colors) 
p_umap2
#ggsave(filename = "./Results/UMAP_T_NK.pdf",p_umap2, width = 5, height = 5)
#ggsave(filename = "../lyme_disease/Manuscript/Figures/Figure2/Figure2A.pdf",p_umap2, width = 5, height = 5)

all.gene <- row.names(T.NK.clean_Annotated)
all.gene[str_detect(all.gene, "IL")]

T.NK.clean_Annotated_dot <- T.NK.clean_Annotated
T.NK.clean_Annotated_dot$celltype <- factor(T.NK.clean_Annotated_dot$celltype, levels = c("Proliferating_Cell","CD56bright_NK", "CD56dim_NK", "IFN_NK", "LEF1_NK",
                                                "GammaDelta_T", "MAIT_T", "Naive_CD8", "GZMB_CD8", "GZMK_CD8", 
                                                "Naive_CD4", "Activated_CD4", "Activated_CD4_2", "Treg", "TFH_CD4"))

Markers_TNK <- c("TOP2A", "PCNA", "XCL1", "TNF", "IFIT3", "CTSW","PRF1", "CCL4", "CST7","CCL5", "GZMA", "GZMH",
                 "CD3E","IL7R", "LEF1", "CCR7","SELL","TRDC","TRGC2", "CD8A", "RORC","GZMK","CD4",
                 "FOXP3","IL2RA", "RORA", "PRDM1", "CTLA4",  "BCL6", "BTG1","RGS1","CXCR5")

Idents(T.NK.clean_Annotated_dot) <-"celltype"

pdot <- DotPlot(T.NK.clean_Annotated_dot, features = Markers_TNK)
pdot <- pdot + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  
pdot <- pdot + scale_color_gradient2(low = "#D3D3D3", high = "#8C1515")
pdot
#ggsave(filename = "./Paper_Figures/dot_marker_NKT2.pdf", pdot, width = 9, height = 5, dpi = 300)
#ggsave(filename = "../lyme_disease/Manuscript/Figures/Figure2/Figure2B.pdf", pdot, width = 9, height = 5, dpi = 300)

rm(T.NK.clean_Annotated_dot)

DefaultAssay(T.NK.clean_Annotated) <- "SCT"
# FeaturePlot(T.NK.clean_Annotated, "TNF", order = T, split.by = "condition")
# FeaturePlot(T.NK.clean_Annotated, "IL10", order = T, split.by = "condition")
# FeaturePlot(T.NK.clean_Annotated, "IL32", order = T, split.by = "condition")
# FeaturePlot(T.NK.clean_Annotated, "RORC", order = T, split.by = "condition")
# FeaturePlot(T.NK.clean_Annotated, "FOXP3", order = T, split.by = "condition", min.cutoff = 0)
# FeaturePlot(T.NK.clean_Annotated, "RORC", order = T, min.cutoff = 0, raster=FALSE)
# FeaturePlot(T.NK.clean_Annotated, "GATA3", order = T, min.cutoff = 0, raster=FALSE)
# FeaturePlot(T.NK.clean_Annotated, "CXCR5", order = T, min.cutoff = 0,raster=FALSE)
# FeaturePlot(T.NK.clean_Annotated, "IFIT3", order = T,min.cutoff = 0,split.by = "condition", raster=FALSE)
# FeaturePlot(T.NK.clean_Annotated, "IL2RA", order = T,min.cutoff = 0, raster=FALSE)
# FeaturePlot(T.NK.clean_Annotated, "CCR6", order = T, min.cutoff = 0, raster=FALSE)
# FeaturePlot(T.NK.clean_Annotated, "TBX21", order = T, min.cutoff = 0, raster=FALSE)
# FeaturePlot(T.NK.clean_Annotated, "IL21", order = T, min.cutoff = 0, raster=FALSE)
# 
# FeaturePlot(T.NK.clean_Annotated, "S100A4", order = T,min.cutoff = 0,raster=FALSE)
# FeaturePlot(T.NK.clean_Annotated, "KLRB1", order = T, split.by = "batch", min.cutoff = 0)
# FeaturePlot(T.NK.clean_Annotated, "FOXP3", order = T, split.by = "batch", min.cutoff = 0)

FeaturePlot(T.NK.clean_Annotated, "AHR", order = T, min.cutoff = 0)


#by subject
count2.2 <- data.frame(table(T.NK.clean_Annotated$subject, T.NK.clean_Annotated$celltype))
total.bysubject <- data.frame(table(T.NK.clean_Annotated$subject))
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

p2.2 <- ggplot(count2.2, aes(x=conditon, y=Percent, color=conditon)) + geom_boxplot(alpha=0.7, outlier.alpha = 0)
p2.2 <- p2.2 + facet_wrap(.~Var2, scales = "free") + theme_cowplot()
p2.2 <- p2.2 + geom_jitter(color="black", size=0.4, alpha=0.9)
p2.2 <- p2.2 + stat_compare_means(method="t.test", comparisons = comparasions) + scale_color_manual(values=condition_color)
p2.2

#ggsave(filename = "../lyme_disease/Manuscript/Figures/FigureS2/Percentage.pdf",p2.2, width = 12, height = 10, dpi = 300 )

#####Treg
T.NK.clean_Annotated_Treg <- subset(T.NK.clean_Annotated, idents = "Treg")
#DimPlot(T.NK.clean_Annotated_Treg,split.by = "condition")
Idents(T.NK.clean_Annotated_Treg) <- "condition"
#FeaturePlot(T.NK.clean_Annotated_Treg, "S100A4")
DefaultAssay(T.NK.clean_Annotated_Treg) <- "integrated"

#T.NK.clean_Annotated_Treg <- PrepSCTFindMarkers(T.NK.clean_Annotated_Treg, verbose = T)

PTLD_RTH <- FindMarkers(T.NK.clean_Annotated_Treg, ident.1 = "PTLD", ident.2 = "RTH", logfc.threshold = log2(1.1), min.pct = 0.3,
                        test.use = "MAST", verbose = TRUE)

PTLDN_HD <- FindMarkers(T.NK.clean_Annotated_Treg, ident.1 = "PTLDN", ident.2 = "HD", logfc.threshold = log2(1.1), min.pct = 0.3,
                        test.use = "MAST", verbose = TRUE)

PTLD_RTH$gene <- rownames(PTLD_RTH)
PTLDN_HD$gene <- rownames(PTLDN_HD)

merged_df <- merge(PTLD_RTH,PTLDN_HD, by="gene") %>% filter(gene != "IL32")
#write.csv(file = "D:/OneDrive - Stanford/Desktop/test.csv", merged_df)

x_max <- max(abs(min(merged_df$avg_log2FC.x, na.rm = TRUE)), max(merged_df$avg_log2FC.x, na.rm = TRUE))
y_max <- max(abs(min(merged_df$avg_log2FC.y, na.rm = TRUE)), max(merged_df$avg_log2FC.y, na.rm = TRUE))

merged_df_sig <- filter(merged_df, p_val_adj.x < 0.05 | p_val_adj.y < 0.05)

filter(merged_df_sig, avg_log2FC.x < 0 & avg_log2FC.y < 0)$gene
filter(merged_df_sig, avg_log2FC.x > 0 & avg_log2FC.y > 0)$gene

filter(merged_df, gene == "S100A4")
filter(merged_df, gene == "KLRB1")
filter(merged_df, gene == "FOXP3")
filter(merged_df, gene == "CD69")

filter(merged_df_sig, gene == "TNF")
filter(merged_df_sig, gene == "CTLA4")
filter(merged_df_sig, gene == "KLRB1")
filter(merged_df_sig, gene == "RGS1")

p.treg <- merged_df_sig %>%
  ggplot(aes(x=avg_log2FC.x, y=avg_log2FC.y, label=gene)) +
  geom_point() +
  coord_cartesian(xlim = c(-x_max, x_max), ylim = c(-y_max, y_max)) +
  scale_x_continuous(sec.axis = sec_axis(~., name = "")) +
  scale_y_continuous(sec.axis = sec_axis(~., name = "")) +
  geom_vline(xintercept = 0, linetype="dashed", color = "black") +
  geom_hline(yintercept = 0, linetype="dashed", color = "black") +
  geom_label_repel() +
  geom_label_repel(data = merged_df_sig %>% filter(gene %in% c("CTLA4","S100A4","KLRB1")), aes(x=avg_log2FC.x, y=avg_log2FC.y, label=gene)) +
  theme_classic() + xlab("RTH <-- | --> PTLD") + ylab("HD <-- | --> PTLDN")
p.treg

#ggsave(filename = "../lyme_disease/Manuscript/Figures/Figure4/Figure4B.pdf",p.treg, width = 5, height = 5, dpi = 300)
#write.csv(file = "../lyme_disease/Manuscript/Figures/Figure4/Figure4B.csv",merged_df)

slice1.VP <- EnhancedVolcano(PTLD_RTH,
                lab = PTLD_RTH$gene,
                FCcutoff=0.15,
                xlim = c(-0.7, 0.7),
                x = 'avg_log2FC',
                y = 'p_val_adj')
slice1.VP

#ggsave(filename = "../lyme_disease/Manuscript/Figures/FigureS4/PTLD_RTH.pdf",slice1.VP, width = 7, height = 7, dpi = 300)

slice3.VP <- EnhancedVolcano(PTLDN_HD,
                lab = PTLDN_HD$gene,
                FCcutoff=0.15,
                xlim = c(-2, 2),
                x = 'avg_log2FC',
                y = 'p_val_adj')
slice3.VP

#ggsave(filename = "../lyme_disease/Manuscript/Figures/FigureS4/PTLDN_HD.pdf",slice3.VP, width = 7, height = 7, dpi = 300)

slice1.VP + slice3.VP

Treg_PTLD_RTH <- PTLD_RTH
Treg_PTLDN_HD <- PTLDN_HD

slice3.VP.treg <- EnhancedVolcano(Treg_PTLDN_HD,
                             lab = Treg_PTLDN_HD$gene,
                             FCcutoff=0.15,
                             xlim = c(-2, 2),
                             x = 'avg_log2FC',
                             y = 'p_val_adj',
                             title = "Treg")
slice3.VP.treg

ggsave(filename = "./DEG.Treg.Treg_PTLDN_HD.pdf", slice3.VP.treg, width = 15, height = 10)

######################################
DefaultAssay(T.NK.clean_Annotated) <- "antibody"
antibody.list <- rownames(T.NK.clean_Annotated)
antibody.list

FeaturePlot(T.NK.clean_Annotated,"CD158e1.KIR3DL1.AHS0211.pAbO",split.by = "condition", raster=T, max.cutoff = 100)
FeaturePlot(T.NK.clean_Annotated,"KIR.NKAT2.KIR2DL2-KIR2DL3-KIR2DS2.AHS0209.pAbO",split.by = "condition", raster=T, max.cutoff = 100)
FeaturePlot(T.NK.clean_Annotated, features = c("CD158e1.KIR3DL1.AHS0211.pAbO", "KIR.NKAT2.KIR2DL2-KIR2DL3-KIR2DS2.AHS0209.pAbO"),
            split.by = "condition", raster=T, max.cutoff = 100, blend = TRUE)

Idents(T.NK.clean_Annotated) <- "celltype"
T.NK.clean_Annotated$CMV <- "No"

T.NK.clean_Annotated$CMV[T.NK.clean_Annotated$subject %in% c("HD6","HD9", "HD10")] <- "Yes"
table(T.NK.clean_Annotated$subject, T.NK.clean_Annotated$CMV)

T.NK.clean_AnnotatedCMV <- subset(T.NK.clean_Annotated, subset = CMV == "No")
table(T.NK.clean_AnnotatedCMV$CMV)

FeaturePlot(T.NK.clean_AnnotatedCMV,"CD158e1.KIR3DL1.AHS0211.pAbO", split.by = "condition", order = T)

table(T.NK.clean_AnnotatedCMV$subject)

EFC_CD8.clean_noCMV <- subset(T.NK.clean_AnnotatedCMV, idents = c("GZMB_CD8","GammaDelta_T"))
EFC_CD8.plot <- DimPlot(EFC_CD8.clean_noCMV)
DimPlot(EFC_CD8.clean_noCMV)

total.EFC_CD8_number <- table(EFC_CD8.clean_noCMV$subject) %>% data.frame()
colnames(total.EFC_CD8_number)[2] <- "Total_Number"

DefaultAssay(EFC_CD8.clean_noCMV) <- "antibody"
antibody.list
p <- RidgePlot(EFC_CD8.clean_noCMV, features = "CD158e1.KIR3DL1.AHS0211.pAbO")
p <- p + xlim(0, 25) + geom_vline(xintercept = 20, linetype=2)
p

#ggsave("../lyme_disease/Manuscript/Figures/Figure4/antibody.expression.distri.pdf",p, width = 7, height = 4, dpi = 300)

#set the antibody name that you want statistics
# EFC_CD8.KIR <- subset(EFC_CD8.clean_noCMV, `CD158e1.KIR3DL1.AHS0211.pAbO` > 20)
# EFC_CD8.KIR <- subset(EFC_CD8.clean_noCMV, `KIR.NKAT2.KIR2DL2-KIR2DL3-KIR2DS2.AHS0209.pAbO` > 20)
EFC_CD8.KIR <- subset(EFC_CD8.clean_noCMV,`CD158e1.KIR3DL1.AHS0211.pAbO` > 20 & `KIR.NKAT2.KIR2DL2-KIR2DL3-KIR2DS2.AHS0209.pAbO` > 20)

KIR_CD8_number <- table(EFC_CD8.KIR$subject) %>% data.frame()
KIR_CD8_number
colnames(KIR_CD8_number)[2] <- "KIR_CD8"

KIR.Percent <- merge(total.EFC_CD8_number,KIR_CD8_number, by = "Var1", all.x = T)
KIR.Percent$KIR_CD8[is.na(KIR.Percent$KIR_CD8)] <- 0
KIR.Percent <-KIR.Percent %>%  mutate(KIR_Percentage = KIR_CD8/Total_Number)
KIR.Percent

KIR.Percent$condition <- "RTH"
KIR.Percent$condition[str_detect(KIR.Percent$Var1,"HD")] <- "HD"
KIR.Percent$condition[str_detect(KIR.Percent$Var1,"PTLD")] <- "PTLD"
KIR.Percent$condition[str_detect(KIR.Percent$Var1,"PTLDN")] <- "PTLDN"

batchlist <- table(EFC_CD8.clean_noCMV$subject, EFC_CD8.clean_noCMV$batch) %>% data.frame()
batchlist <- batchlist %>% filter(Freq >0) %>% select(Var1, Var2) %>% unique
batchlist

KIR.Percent <- merge(KIR.Percent,batchlist,by="Var1")
KIR.Percent

comparasion_HD_PTLDN_RTH_PTLD <- list(c("HD", "PTLDN"), c("RTH", "PTLD"))
KIR.Percent$condition <- factor(KIR.Percent$condition, levels = c("HD", "PTLDN", "RTH", "PTLD"))

#remove the first two batch
p.percent.kir.cd8 <- ggplot(filter(KIR.Percent, ! Var2 %in% c("batch1", "batch2")), aes(x=condition, y=KIR_Percentage)) + geom_jitter() + geom_boxplot(alpha=0.5)
p.percent.kir.cd8 <- p.percent.kir.cd8 + stat_compare_means(comparisons = comparasion_HD_PTLDN_RTH_PTLD,  method = "t.test")
p.percent.kir.cd8

#remove last datapoint
p.percent.kir.cd8_2 <- filter(KIR.Percent, !str_detect(KIR.Percent$Var1,"V7")) %>% 
  filter(! Var2 %in% c("batch1", "batch2")) %>% 
  ggplot(aes(x=condition, y=KIR_Percentage)) + geom_jitter() + geom_boxplot(alpha=0.5)
p.percent.kir.cd8_2 <- p.percent.kir.cd8_2 + stat_compare_means(comparisons = comparasion_HD_PTLDN_RTH_PTLD,method = "t.test")
p.percent.kir.cd8_2

compare <-  filter(KIR.Percent, !str_detect(KIR.Percent$Var1,"V7")) %>% filter(! Var2 %in% c("batch1", "batch2"))
compare$SubjectID <- str_extract(compare$Var1, "\\D+\\d+")

# lme4 <- lmer(KIR_Percentage  ~ condition   + (1|SubjectID), data =compare)
# summary(lme4)

#compare by individual
compare_byindi <- compare %>%
  group_by(SubjectID) %>%
  summarise(mean_KIR_Percentage = mean(KIR_Percentage, na.rm = TRUE)) %>% 
  mutate(condition = str_extract(SubjectID, "\\D+"))

ggplot(compare_byindi, aes(x=condition, y=mean_KIR_Percentage)) + geom_boxplot() + geom_jitter() +  
  stat_compare_means(comparisons = comparasion_HD_PTLDN_RTH_PTLD,method = "t.test")

t.test(filter(compare_byindi, condition == "PTLDN")$mean_KIR_Percentage, 
  filter(compare_byindi, condition == "HD")$mean_KIR_Percentage)

t.test(filter(compare_byindi, condition == "PTLD")$mean_KIR_Percentage, 
       filter(compare_byindi, condition == "RTH")$mean_KIR_Percentage)

patient_data$Var1 <- paste(patient_data$RTH_Patient,patient_data$visit, sep="V")

patient_data$Var1 <- ifelse(patient_data$group == "PTLDN", 
                            patient_data$RTH_Patient, 
                            patient_data$Var1)

compare$Var1 <- as.character(compare$Var1)
compare$Var1[compare$Var1 == "RTH2V5b5"] <- "RTH2V5"
compare$Var1[compare$Var1 == "RTH5V3b5"] <- "RTH5V3"

compare2 <- merge(compare,patient_data, by= "Var1") %>% filter(condition != "PTLDN")

compareKIRCD8 <- compare2 %>% ggplot(aes(x=KIR_Percentage, y = numb_sx, color= condition)) + geom_point() + geom_smooth(method = "lm") #+ facet_wrap(.~condition, scales = "free")
compareKIRCD8 <- compareKIRCD8 + theme_classic() + ggtitle("KIR+ CD8 and Lyme Symptoms Between Groups") + xlab ("KIR+ CD8 Ratio") + ylab("Number of Lyme Symptoms")
compareKIRCD8 <- compareKIRCD8 + scale_color_manual(values = condition_color) + coord_fixed(ratio = 1/500)
compareKIRCD8
#write.csv(file = "./Results/KIR_T_Percent.csv",KIR.Percent)
#write.csv(file = "../lyme_disease/Manuscript/Figures/Figure4/KIR_T_Percent.csv",KIR.Percent)
#ggsave(filename = "../lyme_disease/Manuscript/Figures/Figure4/Figure4C.pdf", compareKIRCD8, height = 5, width = 5, dpi=300)

lm.data <- lm(numb_sx ~ KIR_Percentage : condition, data= compare2)
summary(lm.data)

#find the effector CD8 DEG
EFC_CD8.clean <- subset(T.NK.clean_Annotated, idents= c("GammaDelta_T"))

Idents(EFC_CD8.clean) <- "condition"
DimPlot(EFC_CD8.clean)

DefaultAssay(EFC_CD8.clean) <- "SCT"
DefaultAssay(EFC_CD8.clean) <- "integrated"

EFC_CD8.clean <- PrepSCTFindMarkers(EFC_CD8.clean)

DefaultAssay(EFC_CD8.clean)
PTLD_RTH_EFFCD8 <- FindMarkers(EFC_CD8.clean, ident.1 = "PTLD", ident.2 = "RTH", logfc.threshold = log(1.1),test.use = "MAST", min.pct = 0.3)
PTLDN_HD_EFFCD8 <- FindMarkers(EFC_CD8.clean, ident.1 = "PTLDN", ident.2 = "HD", logfc.threshold = log(1.1),test.use = "MAST", min.pct = 0.3)

PTLD_RTH_EFFCD8$gene <- rownames(PTLD_RTH_EFFCD8)
PTLDN_HD_EFFCD8$gene <- rownames(PTLDN_HD_EFFCD8)

merged_df.EFFCD8 <- merge(PTLD_RTH_EFFCD8,PTLDN_HD_EFFCD8, by="gene") %>% filter(gene != "IL32")

x_max <- max(abs(min(merged_df.EFFCD8$avg_log2FC.x, na.rm = TRUE)), max(merged_df.EFFCD8$avg_log2FC.x, na.rm = TRUE))
y_max <- max(abs(min(merged_df.EFFCD8$avg_log2FC.y, na.rm = TRUE)), max(merged_df.EFFCD8$avg_log2FC.y, na.rm = TRUE))

merged_df.EFFCD8_sig <- filter(merged_df.EFFCD8, p_val_adj.x < 0.05 | p_val_adj.y < 0.05)
  
filter(merged_df.EFFCD8_sig, avg_log2FC.x < 0 & avg_log2FC.y < 0)$gene
filter(merged_df.EFFCD8_sig, avg_log2FC.x > 0 & avg_log2FC.y > 0)$gene

filter(merged_df.EFFCD8_sig, gene == "TRAC")
filter(merged_df.EFFCD8_sig, gene == "EOMES")

p.effcd8 <- merged_df.EFFCD8_sig %>%
  ggplot(aes(x=avg_log2FC.x, y=avg_log2FC.y, label=gene)) + 
  geom_point() +
  coord_cartesian(xlim = c(-x_max, x_max), ylim = c(-y_max, y_max)) +
  scale_x_continuous(sec.axis = sec_axis(~., name = "")) +
  scale_y_continuous(sec.axis = sec_axis(~., name = "")) + 
  geom_vline(xintercept = 0, linetype="dashed", color = "black") +
  geom_hline(yintercept = 0, linetype="dashed", color = "black") +
  geom_label_repel(data = merged_df.EFFCD8_sig %>% filter(gene %in% c("TRAC","TRDC","TRGC2", "IL7R", "MYC", "LAIR2", "FASLG", "GZMK","TBX21","IL2","GATA3", "EOMES")), aes(x=avg_log2FC.x, y=avg_log2FC.y, label=gene)) + 
  theme_classic() + xlab("RTH <-- | --> PTLD") + ylab("HD <-- | --> PTLDN")
p.effcd8

#ggsave(filename = "../lyme_disease/Manuscript/Figures/Figure4/gammadeltaT.pdf",p.effcd8, width = 5, height = 5, dpi = 300)
#write.csv(file = "../lyme_disease/Manuscript/Figures/Figure4/gammadeltaT.csv",merged_df)

slice1.VP <- EnhancedVolcano(PTLD_RTH_EFFCD8,
                             lab = PTLD_RTH_EFFCD8$gene,
                             FCcutoff=0.15,
                             xlim = c(-2, 2),
                             x = 'avg_log2FC',
                             y = 'p_val_adj', 
                             title = "GammaDeltaT")
slice1.VP

#ggsave(filename = "../lyme_disease/Manuscript/Figures/FigureS4/PTLD_RTH_GDT.pdf",slice1.VP, width = 7, height = 7, dpi = 300)

slice3.VP_GDT <- EnhancedVolcano(PTLDN_HD_EFFCD8,
                             lab = PTLDN_HD_EFFCD8$gene,
                             FCcutoff=0.15,
                             xlim = c(-2.3, 2.3),
                             x = 'avg_log2FC',
                             y = 'p_val_adj', 
                             title = "GammaDeltaT")
slice3.VP_GDT

#ggsave(filename = "../lyme_disease/Manuscript/Figures/FigureS4/PTLDN_HD_GDT.pdf",slice3.VP_GDT, width = 7, height = 7, dpi = 300)
#ggsave(filename = "./DEG.GDT.pdf",slice3.VP_GDT,width = 15, height = 10)

GammaDelta_T_PTLD_RTH <- PTLD_RTH_EFFCD8
GammaDelta_T_PTLDN_HD <- PTLDN_HD_EFFCD8

#find the effector Activated CD4 DEG
AC_CD4.clean <- subset(T.NK.clean_Annotated, idents= c("Activated_CD4", "Activated_CD4_2"))

Idents(AC_CD4.clean) <- "condition"
DimPlot(AC_CD4.clean)

DefaultAssay(AC_CD4.clean) <- "SCT"
AC_CD4.clean <- PrepSCTFindMarkers(AC_CD4.clean)
DefaultAssay(AC_CD4.clean) <- "integrated"

SCTResults(object=AC_CD4.clean, slot="umi.assay")

FeaturePlot(T.NK.clean_Annotated, features = "CCR6")

DefaultAssay(AC_CD4.clean) 
PTLD_RTH_ACTCD4 <- FindMarkers(AC_CD4.clean, ident.1 = "PTLD", ident.2 = "RTH", logfc.threshold = log(1.1), test.use = "MAST", min.pct = 0.2)
PTLDN_HD_ACTCD4 <- FindMarkers(AC_CD4.clean, ident.1 = "PTLDN", ident.2 = "HD", logfc.threshold = log(1.1), test.use = "MAST", min.pct = 0.2)

PTLD_RTH_ACTCD4$gene <- rownames(PTLD_RTH_ACTCD4)
PTLDN_HD_ACTCD4$gene <- rownames(PTLDN_HD_ACTCD4)

merged_df.ACTCD4 <- merge(PTLD_RTH_ACTCD4,PTLDN_HD_ACTCD4, by="gene") %>% filter(gene != "IL32")

x_max <- max(abs(min(merged_df.ACTCD4$avg_log2FC.x, na.rm = TRUE)), max(merged_df.ACTCD4$avg_log2FC.x, na.rm = TRUE))
y_max <- max(abs(min(merged_df.ACTCD4$avg_log2FC.y, na.rm = TRUE)), max(merged_df.ACTCD4$avg_log2FC.y, na.rm = TRUE))

merged_df.ACTCD4_sig <- filter(merged_df.ACTCD4, p_val_adj.x < 0.05 | p_val_adj.y < 0.05)

filter(merged_df.ACTCD4_sig, avg_log2FC.x < 0 & avg_log2FC.y < 0)$gene
filter(merged_df.ACTCD4_sig, avg_log2FC.x > 0 & avg_log2FC.y > 0)$gene

filter(merged_df.EFFCD8_sig, gene == "TRAC")
filter(merged_df.ACTCD4_sig, gene == "RORC")
filter(merged_df.ACTCD4_sig, gene == "CD52")

p.actcd4 <- merged_df.ACTCD4_sig %>%
  ggplot(aes(x=avg_log2FC.x, y=avg_log2FC.y, label=gene)) + 
  geom_point() +
  coord_cartesian(xlim = c(-x_max, x_max), ylim = c(-y_max, y_max)) +
  scale_x_continuous(sec.axis = sec_axis(~., name = "")) +
  scale_y_continuous(sec.axis = sec_axis(~., name = "")) + 
  geom_vline(xintercept = 0, linetype="dashed", color = "black") +
  geom_hline(yintercept = 0, linetype="dashed", color = "black") +
  geom_label_repel() + theme_classic() + xlab("RTH <-- | --> PTLD") + ylab("HD <-- | --> PTLDN") + ggtitle("Mature CD4")
p.actcd4

#write.csv(file = "../lyme_disease/Manuscript/Figures/Figure4/CD4T_DEG.csv")
#ggsave(filename = "../lyme_disease/Manuscript/Figures/Figure2/Figure2F.pdf", p.actcd4, width = 20, height = 20, dpi = 300)

CD4_PTLD_RTH <- PTLD_RTH_ACTCD4
CD4_PTLDN_HD <- PTLDN_HD_ACTCD4

slice3.VP_CD4 <- EnhancedVolcano(CD4_PTLDN_HD,
                             lab = CD4_PTLDN_HD$gene,
                             FCcutoff=0.15,
                             xlim = c(-2.3, 2.3),
                             x = 'avg_log2FC',
                             y = 'p_val_adj', 
                             title = "CD4_activated")
slice3.VP_CD4

#ggsave(filename = "./DEG.CD4.pdf",slice3.VP_CD4,width = 15, height = 10)

############
NK.clean <- subset(T.NK.clean_Annotated, idents= c("CD56dim_NK", "CD56bright_NK","LEF1_NK","IFN_NK"))
DimPlot(NK.clean)

Idents(NK.clean) <- "condition"

FeaturePlot(NK.clean, "B3GAT1", order = T, split.by = "condition")

DefaultAssay(NK.clean) <- "integrated"

DefaultAssay(NK.clean) 
PTLD_RTH_NK <- FindMarkers(NK.clean, ident.1 = "PTLD", ident.2 = "RTH", logfc.threshold = log(1.1), test.use = "MAST", min.pct = 0.2)
PTLDN_HD_NK <- FindMarkers(NK.clean, ident.1 = "PTLDN", ident.2 = "HD", logfc.threshold = log(1.1), test.use = "MAST", min.pct = 0.2)

PTLD_RTH_NK$gene <- rownames(PTLD_RTH_NK)
PTLDN_HD_NK$gene <- rownames(PTLDN_HD_NK)

merged_df.NK <- merge(PTLD_RTH_NK,PTLDN_HD_NK, by="gene") %>% filter(gene != "IL32")

x_max <- max(abs(min(merged_df.NK$avg_log2FC.x, na.rm = TRUE)), max(merged_df.NK$avg_log2FC.x, na.rm = TRUE))
y_max <- max(abs(min(merged_df.NK$avg_log2FC.y, na.rm = TRUE)), max(merged_df.NK$avg_log2FC.y, na.rm = TRUE))

merged_df.NK_sig <- filter(merged_df.NK, p_val_adj.x < 0.05 | p_val_adj.y < 0.05)

filter(merged_df.NK_sig, avg_log2FC.x < 0 & avg_log2FC.y < 0)$gene
filter(merged_df.NK_sig, avg_log2FC.x > 0 & avg_log2FC.y > 0)$gene

p.nk <- merged_df.NK_sig %>%
  ggplot(aes(x=avg_log2FC.x, y=avg_log2FC.y, label=gene)) + 
  geom_point() +
  coord_cartesian(xlim = c(-x_max, x_max), ylim = c(-y_max, y_max)) +
  scale_x_continuous(sec.axis = sec_axis(~., name = "")) +
  scale_y_continuous(sec.axis = sec_axis(~., name = "")) + 
  geom_vline(xintercept = 0, linetype="dashed", color = "black") +
  geom_hline(yintercept = 0, linetype="dashed", color = "black") +
  geom_label_repel() + theme_classic() + xlab("RTH <-- | --> PTLD") + ylab("HD <-- | --> PTLDN") + ggtitle("NK")
p.nk

NK_PTLD_RTH <- PTLD_RTH_NK
NK_PTLDN_HD <- PTLDN_HD_NK


slice3.VP.NK <- EnhancedVolcano(PTLDN_HD_NK,
                             lab = PTLDN_HD_NK$gene,
                             FCcutoff=0.15,
                             xlim = c(-2.3, 2.3),
                             x = 'avg_log2FC',
                             y = 'p_val_adj', 
                             title = "NK cell")
slice3.VP.NK

#ggsave(filename = "./DEG.NK.pdf",slice3.VP.NK,width = 15, height = 10)

####################################
CD8.clean <- subset(T.NK.clean_Annotated, idents= c("GZMB_CD8", "GZMK_CD8"))

Idents(CD8.clean) <- "condition"
DimPlot(CD8.clean)
FeaturePlot(CD8.clean, "GZMK", order = T, split.by = "condition")
DefaultAssay(CD8.clean) <- "integrated"

DefaultAssay(CD8.clean) 
PTLD_RTH_CD8 <- FindMarkers(CD8.clean, ident.1 = "PTLD", ident.2 = "RTH", logfc.threshold = log(1.1), test.use = "MAST", min.pct = 0.2)
PTLDN_HD_CD8 <- FindMarkers(CD8.clean, ident.1 = "PTLDN", ident.2 = "HD", logfc.threshold = log(1.1), test.use = "MAST", min.pct = 0.2)

PTLD_RTH_CD8$gene <- rownames(PTLD_RTH_CD8)
PTLDN_HD_CD8$gene <- rownames(PTLDN_HD_CD8)

merged_df.CD8 <- merge(PTLD_RTH_CD8,PTLDN_HD_CD8, by="gene") %>% filter(gene != "IL32")

x_max <- max(abs(min(merged_df.CD8$avg_log2FC.x, na.rm = TRUE)), max(merged_df.CD8$avg_log2FC.x, na.rm = TRUE))
y_max <- max(abs(min(merged_df.CD8$avg_log2FC.y, na.rm = TRUE)), max(merged_df.CD8$avg_log2FC.y, na.rm = TRUE))

merged_df.CD8_sig <- filter(merged_df.CD8, p_val_adj.x < 0.05 | p_val_adj.y < 0.05)

filter(merged_df.CD8_sig, avg_log2FC.x < 0 & avg_log2FC.y < 0)$gene
filter(merged_df.CD8_sig, avg_log2FC.x > 0 & avg_log2FC.y > 0)$gene

p.cd8 <- merged_df.CD8_sig %>%
  ggplot(aes(x=avg_log2FC.x, y=avg_log2FC.y, label=gene)) + 
  geom_point() +
  coord_cartesian(xlim = c(-x_max, x_max), ylim = c(-y_max, y_max)) +
  scale_x_continuous(sec.axis = sec_axis(~., name = "")) +
  scale_y_continuous(sec.axis = sec_axis(~., name = "")) + 
  geom_vline(xintercept = 0, linetype="dashed", color = "black") +
  geom_hline(yintercept = 0, linetype="dashed", color = "black") +
  geom_label_repel() + theme_classic() + xlab("RTH <-- | --> PTLD") + ylab("HD <-- | --> PTLDN") + ggtitle("CD8")
p.cd8

CD8_PTLD_RTH <- PTLD_RTH_CD8
CD8_PTLDN_HD <- PTLDN_HD_CD8

slice3.VP_CD8 <- EnhancedVolcano(CD8_PTLDN_HD,
                             lab = CD8_PTLDN_HD$gene,
                             FCcutoff=0.15,
                             xlim = c(-2.3, 2.3),
                             x = 'avg_log2FC',
                             y = 'p_val_adj', 
                             title = "CD8 cell")
slice3.VP_CD8
#ggsave(filename = "./DEG.CD8.pdf",slice3.VP_CD8,width = 15, height = 10)

#TFH
TFH.clean <- subset(T.NK.clean_Annotated, idents= c("TFH_CD4"))

Idents(TFH.clean) <- "condition"
FeaturePlot(TFH.clean, "CXCL5", order = T, split.by = "condition")
DefaultAssay(TFH.clean) <- "integrated"

DefaultAssay(TFH.clean) 
PTLD_RTH_TFH <- FindMarkers(TFH.clean, ident.1 = "PTLD", ident.2 = "RTH", logfc.threshold = log(1.1), test.use = "MAST", min.pct = 0.3)
PTLDN_HD_TFH <- FindMarkers(TFH.clean, ident.1 = "PTLDN", ident.2 = "HD", logfc.threshold = log(1.1), test.use = "MAST", min.pct = 0.3)

PTLD_RTH_TFH$gene <- rownames(PTLD_RTH_TFH)
PTLDN_HD_TFH$gene <- rownames(PTLDN_HD_TFH)

merged_df.TFH <- merge(PTLD_RTH_TFH,PTLDN_HD_TFH, by="gene") %>% filter(gene != "IL32")

x_max <- max(abs(min(merged_df.TFH$avg_log2FC.x, na.rm = TRUE)), max(merged_df.TFH$avg_log2FC.x, na.rm = TRUE))
y_max <- max(abs(min(merged_df.TFH$avg_log2FC.y, na.rm = TRUE)), max(merged_df.TFH$avg_log2FC.y, na.rm = TRUE))

merged_df.TFH_sig <- filter(merged_df.TFH, p_val_adj.x < 0.05 | p_val_adj.y < 0.05)

filter(merged_df.TFH_sig, avg_log2FC.x < 0 & avg_log2FC.y < 0)$gene
filter(merged_df.TFH_sig, avg_log2FC.x > 0 & avg_log2FC.y > 0)$gene

p.tfh <- merged_df.TFH_sig %>%
  ggplot(aes(x=avg_log2FC.x, y=avg_log2FC.y, label=gene)) + 
  geom_point() +
  coord_cartesian(xlim = c(-x_max, x_max), ylim = c(-y_max, y_max)) +
  scale_x_continuous(sec.axis = sec_axis(~., name = "")) +
  scale_y_continuous(sec.axis = sec_axis(~., name = "")) + 
  geom_vline(xintercept = 0, linetype="dashed", color = "black") +
  geom_hline(yintercept = 0, linetype="dashed", color = "black") +
  geom_label_repel() + theme_classic() + xlab("RTH <-- | --> PTLD") + ylab("HD <-- | --> PTLDN") + ggtitle("TFH")
p.tfh

TFH_PTLD_RTH <- PTLD_RTH_TFH
TFH_PTLDN_HD <- PTLDN_HD_TFH

slice3.VP_TFH <- EnhancedVolcano(TFH_PTLDN_HD,
                             lab = TFH_PTLDN_HD$gene,
                             FCcutoff=0.15,
                             xlim = c(-2.3, 2.3),
                             x = 'avg_log2FC',
                             y = 'p_val_adj', 
                             title = "TFH cell")
slice3.VP_TFH
#ggsave(filename = "./DEG.TFH.pdf",slice3.VP_TFH,width = 15, height = 10)

#MAIT
MAIT.clean <- subset(T.NK.clean_Annotated, idents= c("MAIT_T"))

Idents(MAIT.clean) <- "condition"
FeaturePlot(MAIT.clean, "RORC", order = T, split.by = "condition")
DefaultAssay(MAIT.clean) <- "integrated"

DefaultAssay(MAIT.clean) 
PTLD_RTH_MAIT <- FindMarkers(MAIT.clean, ident.1 = "PTLD", ident.2 = "RTH", logfc.threshold = log(1.1), test.use = "MAST", min.pct = 0.3)
PTLDN_HD_MAIT <- FindMarkers(MAIT.clean, ident.1 = "PTLDN", ident.2 = "HD", logfc.threshold = log(1.1), test.use = "MAST", min.pct = 0.3)

PTLD_RTH_MAIT$gene <- rownames(PTLD_RTH_MAIT)
PTLDN_HD_MAIT$gene <- rownames(PTLDN_HD_MAIT)

merged_df.MAIT <- merge(PTLD_RTH_MAIT,PTLDN_HD_MAIT, by="gene") %>% filter(gene != "IL32")

x_max <- max(abs(min(merged_df.MAIT$avg_log2FC.x, na.rm = TRUE)), max(merged_df.MAIT$avg_log2FC.x, na.rm = TRUE))
y_max <- max(abs(min(merged_df.MAIT$avg_log2FC.y, na.rm = TRUE)), max(merged_df.MAIT$avg_log2FC.y, na.rm = TRUE))

merged_df.MAIT_sig <- filter(merged_df.MAIT, p_val_adj.x < 0.05 | p_val_adj.y < 0.05)

filter(merged_df.MAIT_sig, avg_log2FC.x < 0 & avg_log2FC.y < 0)$gene
filter(merged_df.MAIT_sig, avg_log2FC.x > 0 & avg_log2FC.y > 0)$gene

p.mait <- merged_df.MAIT_sig %>%
  ggplot(aes(x=avg_log2FC.x, y=avg_log2FC.y, label=gene)) + 
  geom_point() +
  coord_cartesian(xlim = c(-x_max, x_max), ylim = c(-y_max, y_max)) +
  scale_x_continuous(sec.axis = sec_axis(~., name = "")) +
  scale_y_continuous(sec.axis = sec_axis(~., name = "")) + 
  geom_vline(xintercept = 0, linetype="dashed", color = "black") +
  geom_hline(yintercept = 0, linetype="dashed", color = "black") +
  geom_label_repel() + theme_classic() + xlab("RTH <-- | --> PTLD") + ylab("HD <-- | --> PTLDN") + ggtitle("MAIT")
p.mait

MAIT_PTLD_RTH <- PTLD_RTH_MAIT
MAIT_PTLDN_HD <- PTLDN_HD_MAIT

slice3.VP_MAIT <- EnhancedVolcano(MAIT_PTLDN_HD,
                             lab = MAIT_PTLDN_HD$gene,
                             FCcutoff=0.15,
                             xlim = c(-2.3, 2.3),
                             x = 'avg_log2FC',
                             y = 'p_val_adj', 
                             title = "MAIT cell")
slice3.VP_MAIT
#ggsave(filename = "./DEG.MAIT.pdf",slice3.VP_MAIT,width = 15, height = 10)

################################
Treg_PTLDN_HD$celltype <- "Treg"
GammaDelta_T_PTLDN_HD$celltype <- "GammaDelta_T"
CD4_PTLDN_HD$celltype <- "Memory_CD4"
NK_PTLDN_HD$celltype <- "NK"
CD8_PTLDN_HD$celltype <- "Memory_CD8"
TFH_PTLDN_HD$celltype <- "TFH"
MAIT_PTLDN_HD$celltype <- "MAIT"

PTLDN_HD_ALL <- rbind(Treg_PTLDN_HD,GammaDelta_T_PTLDN_HD,CD4_PTLDN_HD,NK_PTLDN_HD,CD8_PTLDN_HD,TFH_PTLDN_HD,MAIT_PTLDN_HD)

#write.csv(file = "../lyme_disease/Manuscript/Figures/Figure2/All_lympho_PTLDN_HD_DEG.csv", PTLDN_HD_ALL)
#PTLDN_HD_ALL <- read.csv("../lyme_disease/Manuscript/Figures/Figure2/All_lympho_PTLDN_HD_DEG.csv")

#find markers for actcd4
Idents(T.NK.clean_Annotated) <- "celltype"
PrepSCTFindMarkers(T.NK.clean_Annotated)
cluster0.markers <- FindMarkers(T.NK.clean_Annotated, ident.1 = "Activated_CD4", logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
T.NK.clean_Annotated[[]]

DefaultAssay(T.NK.clean_Annotated.subset) <- "SCT"
RORC_Feature <- FeaturePlot(T.NK.clean_Annotated.subset, features = "RORC" , order = T) + scale_color_gradient(low = "#D3D3D3", high = "#8C1515") + coord_equal()
CCR6_Feature <- FeaturePlot(T.NK.clean_Annotated.subset, features = "CCR6" ,order = T) + scale_color_gradient(low = "#D3D3D3", high = "#8C1515")+ coord_equal()
FeaturePlot(T.NK.clean_Annotated.subset, features = "IL17E" ,order = T) + scale_color_gradient(low = "#D3D3D3", high = "#8C1515")+ coord_equal()
CD161_Feature <- FeaturePlot(T.NK.clean_Annotated.subset, features = "KLRB1" ,order = T) + scale_color_gradient(low = "#D3D3D3", high = "#8C1515")+ coord_equal()
CCR6_Feature <- FeaturePlot(T.NK.clean_Annotated.subset, features = "CCR6" ,order = T) + scale_color_gradient(low = "#D3D3D3", high = "#8C1515")+ coord_equal()

RORC_Feature
CCR6_Feature
CD161_Feature
#ggsave(filename = "../lyme_disease/Manuscript/Figures/Figure2/Figure2D.1.pdf", RORC_Feature, width = 5, height = 5, dpi = 300)  
#ggsave(filename = "../lyme_disease/Manuscript/Figures/Figure2/Figure2D.2.pdf", CCR6_Feature, width = 5, height = 5, dpi = 300)

######Figure2 还有点问题
gene.for.anno <- read_excel("./Results/genelistfor annotation.xlsx")

TFH.gene <- filter(gene.for.anno, celltype == "TFH")$gene
CD4_T_genes <- filter(gene.for.anno, celltype == "CD4 T")$gene
CD8_T_genes <- filter(gene.for.anno, celltype == "CD8 T")$gene
gammaDelta_T_genes <- filter(gene.for.anno, celltype == "gammaDelta T")$gene
MAIT_genes <- filter(gene.for.anno, celltype == "MAIT")$gene
NK_genes <- filter(gene.for.anno, celltype == "NK")$gene
TFH_genes <- filter(gene.for.anno, celltype == "TFH")$gene
Treg_genes <- filter(gene.for.anno, celltype == "Treg")$gene

p <- filter(PTLDN_HD_ALL,p_val_adj < 0.05) %>%  ggplot(aes(x = celltype, y = avg_log2FC, label = gene)) +
  geom_jitter(aes(color = (gene %in% Treg_genes & celltype == "Treg") | 
                    (gene %in% CD8_T_genes & celltype == "Memory_CD8")|
                    (gene %in% gammaDelta_T_genes & celltype == "GammaDelta_T")|
                    (gene %in% MAIT_genes & celltype == "MAIT")|
                    (gene %in% CD4_T_genes & celltype == "Memory_CD4")|
                    (gene %in% NK_genes & celltype == "NK")|
                    (gene %in% TFH.gene & celltype == "TFH"))) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
  
  geom_text_repel(data = subset(PTLDN_HD_ALL, (gene %in% Treg_genes & celltype == "Treg")), vjust = -1, hjust = 0.5, color = "red") +
  geom_text_repel(data = subset(PTLDN_HD_ALL, (gene %in% CD8_T_genes  & celltype == "Memory_CD8")), vjust = -1, hjust = 0.5, color = "red") +
  geom_text_repel(data = subset(PTLDN_HD_ALL, (gene %in% gammaDelta_T_genes  & celltype == "GammaDelta_T")), vjust = -1, hjust = 0.5, color = "red") +
  geom_text_repel(data = subset(PTLDN_HD_ALL, (gene %in% MAIT_genes  & celltype == "MAIT")), vjust = -1, hjust = 0.5, color = "red") +
  geom_text_repel(data = subset(PTLDN_HD_ALL, (gene %in% CD4_T_genes  & celltype == "Memory_CD4")), vjust = -1, hjust = 0.5, color = "red") +
  geom_text_repel(data = subset(PTLDN_HD_ALL, (gene %in% NK_genes  & celltype == "NK")), vjust = -1, hjust = 0.5, color = "red") +
  geom_text_repel(data = subset(PTLDN_HD_ALL, (gene %in% TFH.gene  & celltype == "TFH")), vjust = -1, hjust = 0.5, color = "red") +
  
  theme_minimal() +  
  labs(title = "Differentially Expressed Genes By Cell Type",
       x = "Cell Type",
       y = "Log2 Fold Change") + XZ_flip_x_label() +  guides(color = FALSE)
p  

#ggsave(filename = "../lyme_disease/Manuscript/Figures/Figure2/Gene.List.pdf", p, width = 12, height = 8, dpi=300)

##########
count2.2
wide_data <- dcast(count2.2, Var1 + conditon + SLICE ~ Var2, value.var = "Percent")
patient_data


