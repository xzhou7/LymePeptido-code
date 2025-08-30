#Single cell Lyme disease post merge
#Author: Xin Chen, Ph.D.
#Date Created: 032723
#T cell analysis

#set up working directory
setwd(dir = "~/Library/CloudStorage/Box-Box/Xin.Chen.Shareable/R7_NR")
getwd()

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

load("./Data/022723_Step2.3_Subset5000_RunUMAP.RData")
p1 <- DimPlot(pbmc.integrated.5000, group.by = "seurat_clusters", label = T)
p1

T_NK <- subset(pbmc.integrated.5000, seurat_clusters %in% c("4", "29", "17", "26", "23", "7", "11", "10", "12", "20",
                                                               "15", "3", "1", "0", "8", "33"))
p1 <- DimPlot(T_NK, group.by = "seurat_clusters", label = T) + coord_fixed(ratio=1)
ggsave(filename = "./Results/T_NK/T_NK_UMAP.pdf", p1, height = 8, width = 9, dpi = 300)

FeaturePlot(T_NK, "FOXP3",min.cutoff = 0, sort.cell = T)
FeaturePlot(T_NK, "CCL5",min.cutoff = 0, sort.cell = T)
FeaturePlot(T_NK, "CD8A", min.cutoff = 0, sort.cell = T)
FeaturePlot(T_NK, "CD4", min.cutoff = 0, sort.cell = T)

FeaturePlot(T_NK, "GNLY", min.cutoff = 0, sort.cell = T)
FeaturePlot(T_NK, "GZMB", min.cutoff = 0, sort.cell = T)
FeaturePlot(T_NK, "GZMK", min.cutoff = 0, sort.cell = T)
FeaturePlot(T_NK, "CCR7", min.cutoff = 0, sort.cell = T)
FeaturePlot(T_NK, "CCR6", min.cutoff = 0, sort.cell = T)
FeaturePlot(T_NK, "TRDV2", min.cutoff = 0, sort.cell = T)

#Clean-up

CL0 <- subset(T_NK, seurat_clusters == 0)
cell.id.0 <- CellSelector(DimPlot(CL0))

CL1 <- subset(T_NK, seurat_clusters == 1)
cell.id.1 <- CellSelector(DimPlot(CL1))

CL3 <- subset(T_NK, seurat_clusters == 3)
cell.id.3 <- CellSelector(DimPlot(CL3))

CL4 <- subset(T_NK, seurat_clusters == 4)
cell.id.4 <- CellSelector(DimPlot(CL4))

CL7 <- subset(T_NK, seurat_clusters == 7)
cell.id.7 <- CellSelector(DimPlot(CL7))

CL8 <- subset(T_NK, seurat_clusters == 8)
cell.id.8 <- CellSelector(DimPlot(CL8))

CL10 <- subset(T_NK, seurat_clusters == 10)
cell.id.10 <- CellSelector(DimPlot(CL10))

CL11 <- subset(T_NK, seurat_clusters == 11)
cell.id.11 <- CellSelector(DimPlot(CL11))

CL12 <- subset(T_NK, seurat_clusters == 12)
cell.id.12 <- CellSelector(DimPlot(CL12))

CL15 <- subset(T_NK, seurat_clusters == 15)
cell.id.15 <- CellSelector(DimPlot(CL15))

CL17 <- subset(T_NK, seurat_clusters == 17)
cell.id.17 <- CellSelector(DimPlot(CL17))

CL20 <- subset(T_NK, seurat_clusters == 20)
cell.id.20 <- CellSelector(DimPlot(CL20))

CL23 <- subset(T_NK, seurat_clusters == 23)
cell.id.23 <- CellSelector(DimPlot(CL23))

CL26 <- subset(T_NK, seurat_clusters == 26)
cell.id.26 <- CellSelector(DimPlot(CL26))

CL29 <- subset(T_NK, seurat_clusters == 29)
cell.id.29 <- CellSelector(DimPlot(CL29))

CL33 <- subset(T_NK, seurat_clusters == 33)
cell.id.33 <- CellSelector(DimPlot(CL33))


subset.id <- c(cell.id.0,cell.id.1, cell.id.3, cell.id.4, cell.id.7, cell.id.8,
               cell.id.10, cell.id.11, cell.id.12, cell.id.15, cell.id.17, cell.id.20, cell.id.23,cell.id.26, cell.id.29, cell.id.33 )
 
save(subset.id, file = "./Data/Step5_032723_T_NK_subset.id.RData")

T_NK.clean <- subset(pbmc.integrated.5000, cells = subset.id)
save(T_NK.clean, file = "./Data/Step5.1_032723_T_NK_clean.RData")



#################################
#set up working directory
setwd(dir = "~/Library/CloudStorage/Box-Box/Xin.Chen.Shareable/R7_NR")
getwd()

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
library(RColorBrewer)

load("./Data/Step5.1_032723_T_NK_clean.RData")
p1 <- DimPlot(T_NK.clean, group.by = "seurat_clusters", label = T)
p1
ggsave(filename = "./Results/032923_T_NK_UMAP.pdf", p1, width = 5, height = 4, dpi = 300)
FeaturePlot(T_NK.clean, "CD8A", min.cutoff = 0, sort.cell = T)

DefaultAssay(T_NK.clean) <- "antibody"
P = rownames(T_NK.clean)
P[P %like% "KIR"]

DefaultAssay(T_NK.clean) <- "integrated"
P = rownames(T_NK.clean)
P[P %like% "KIR"]

A <- subset(T_NK.clean, seurat_clusters %in% c("33"))
DefaultAssay(A) <- "integrated"
FeaturePlot(A, "CD8A", min.cutoff = 0, sort.cell = T)
FeaturePlot(A, "CD4", min.cutoff = 0, sort.cell = T)
FeaturePlot(A, "CCR6", min.cutoff = 0, sort.cell = T)
FeaturePlot(A, "NKG7", min.cutoff = 0, sort.cell = T)
FeaturePlot(A, "CCR6", min.cutoff = 0, sort.cell = T)


FeaturePlot(T_NK.clean, "CD8A", min.cutoff = 0, sort.cell = T)
FeaturePlot(T_NK.clean, "CD4", min.cutoff = 0, sort.cell = T)
FeaturePlot(T_NK.clean, "GNLY", min.cutoff = 0, sort.cell = T)
FeaturePlot(T_NK.clean, "GZMB", min.cutoff = 0, sort.cell = T)
FeaturePlot(T_NK.clean, "GZMK", min.cutoff = 0, sort.cell = T)
FeaturePlot(T_NK.clean, "CCR7", min.cutoff = 0, sort.cell = T)
FeaturePlot(T_NK.clean, "CCR6", min.cutoff = 0, sort.cell = T)
FeaturePlot(T_NK.clean, "CD158e1.KIR3DL1.AHS0211.pAbO", max.cutoff = 100, sort.cell = T)
FeaturePlot(T_NK.clean, "KIR.NKAT2.KIR2DL2-KIR2DL3-KIR2DS2.AHS0209.pAbO", max.cutoff = 100, sort.cell = T)
FeaturePlot(T_NK.clean, "KIR2DL1", max.cutoff = 2, sort.cell = T)

#######subset CD8T and CD4T#########
Tcell <- subset(pbmc.integrated.5000, seurat_clusters %in% c("23", "7", "11", "10", "12", "20","15", "3", "1", "0", "8"))
DefaultAssay(Tcell) <- "integrated"
Tcell_variable <- FindVariableFeatures(Tcell, selection.method = "vst", nfeatures = 100,verbose = FALSE)
Tcell_scale <- ScaleData(Tcell_variable, verbose = T)
Tcell_PCA <- RunPCA(Tcell_scale, npcs = 50, verbose = T)
ElbowPlot(Tcell_PCA)
Tcell_UMAP <- RunUMAP(Tcell_PCA, reduction = "pca", dims = 1:10, verbose = FALSE)
Tcell_Nei <-FindNeighbors(Tcell_UMAP, reduction = "pca", dims= 1:10)
Tcell.integrated <-FindClusters(Tcell_Nei, resolution = 1)
p1 <- DimPlot(Tcell.integrated, label = T)
p1
ggsave(filename = "./Results/032923_Tcell with nfeTURE100_dim1:10_res_1.pdf", p1, width = 5, height = 4, dpi = 300)
save(Tcell.integrated, file = "./Data/032923_Step5_Tcell.Rdata")

FeaturePlot(Tcell.integrated, "CD8A", min.cutoff = 0, sort.cell = T)
FeaturePlot(Tcell.integrated, "CD4", min.cutoff = 0, sort.cell = T)
FeaturePlot(Tcell.integrated, "FOXP3", min.cutoff = 0, sort.cell = T)
FeaturePlot(Tcell.integrated, "GNLY", min.cutoff = 0, sort.cell = T)
FeaturePlot(Tcell.integrated, "GZMB", min.cutoff = 0, sort.cell = T)
FeaturePlot(Tcell.integrated, "GZMK", min.cutoff = 0, sort.cell = T)
FeaturePlot(Tcell.integrated, "CCR7", min.cutoff = 0, sort.cell = T)
FeaturePlot(Tcell.integrated, "CCR6", min.cutoff = 0, sort.cell = T)
FeaturePlot(Tcell.integrated, "CD158e1.KIR3DL1.AHS0211.pAbO", max.cutoff = 100, sort.cell = T)
FeaturePlot(Tcell.integrated, "KIR.NKAT2.KIR2DL2-KIR2DL3-KIR2DS2.AHS0209.pAbO", max.cutoff = 100, sort.cell = T)
FeaturePlot(Tcell.integrated, "KIR2DL1", max.cutoff = 2, sort.cell = T)

######subset CD8 T cells########################
CD8T <- subset(Tcell.integrated, seurat_clusters %in% c("0", "6", "7",  "12", "15", "2"))
p1 <- DimPlot(CD8T, label = T)
p1
ggsave(filename = "./Results/032923_CD8T_b4_cleanup.pdf", p1, width = 5, height = 4, dpi = 300)
#####Cleanup CD8 T ########
CL0 <- subset(CD8T, seurat_clusters == 0)
cell.id.0 <- CellSelector(DimPlot(CL0))

CL7 <- subset(CD8T, seurat_clusters == 7)
cell.id.7 <- CellSelector(DimPlot(CL7))

CL6 <- subset(CD8T, seurat_clusters == 6)
cell.id.6 <- CellSelector(DimPlot(CL6))

CL12 <- subset(CD8T, seurat_clusters == 12)
cell.id.12 <- CellSelector(DimPlot(CL12))

CL15 <- subset(CD8T, seurat_clusters == 15)
cell.id.15 <- CellSelector(DimPlot(CL15))

CL2 <- subset(CD8T, seurat_clusters == 2)
cell.id.2 <- CellSelector(DimPlot(CL2))

subset.id <- c(cell.id.0, cell.id.7, cell.id.6, cell.id.12, cell.id.15, cell.id.2 )

save(subset.id, file = "./Data/Step5_032723_CD8T_subset.id.RData")

CD8T.clean <- subset(Tcell.integrated, cells = subset.id)
save(CD8T.clean, file = "./Data/Step5.1_032923_CD8T_clean.RData")

##########################CD8T######################
#set up working directory
setwd(dir = "~/Library/CloudStorage/Box-Box/Xin.Chen.Shareable/R7_NR")
getwd()

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

load("./Data/Step5.1_032923_CD8T_clean.RData")
p1 <- DimPlot(CD8T.clean, group.by = "seurat_clusters", label = T, split.by = "condition")
p1
ggsave(filename = "./Results/CD8 T/032923_CD8T_clean_UMAP.pdf", p1, width = 13, height = 4, dpi = 300)
FeaturePlot(CD8T.clean, "CD158e1.KIR3DL1.AHS0211.pAbO", max.cutoff = 100, sort.cell = T, split.by = "condition")

#find markers for cluster
pbmc.markers <- FindAllMarkers(CD8T.clean, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>%
  group_by(cluster) %>%
  slice_max(n = 10, order_by = avg_log2FC)

write.csv(pbmc.markers,"./Results/CD8T.clean_gene_markers.csv", row.names = F)
table(CD8T.clean$integrated_snn_res.1)

CD8T.clean$cluster <- "Naive"
CD8T.clean$cluster[CD8T.clean$integrated_snn_res.1 == 2] <- "GZMB+KIR+"
CD8T.clean$cluster[CD8T.clean$integrated_snn_res.1 == 6] <- "GZMK+"
CD8T.clean$cluster[CD8T.clean$integrated_snn_res.1 == 7] <- "CCR6+CD161+"
CD8T.clean$cluster[CD8T.clean$integrated_snn_res.1 == 12] <- "TCRgd+"
CD8T.clean$cluster[CD8T.clean$integrated_snn_res.1 == 15] <- "GNLY+"

cluster_color = c(
  "GZMB+KIR+" = ggsci::pal_d3()(n=3)[1], 
  "GZMK+" = ggsci::pal_d3()(n=3)[2], 
  "CCR6+CD161+" = ggsci::pal_d3()(n=3)[2],
  "TCRgd+" = ggsci::pal_d3()(n=3)[2],
  "GNLY+" = ggsci::pal_d3()(n=3)[2],
  "Naive" = ggsci::pal_d3()(n=3)[3]
)

p2 <- DimPlot(CD8T.clean, group.by = "cluster", split.by = "condition") + scale_color_manual(values = cluster_color) + coord_fixed(ratio=1)
p2

p1 <- DimPlot(CD8T.clean, group.by = "cluster", label = T, split.by = "condition")
p1
save(CD8T.clean, file = "./Data/032923_CD8T_Clean.Rdata")
ggsave(filename = "./Results/032923_CD8T_clean_UMAP_marker.pdf", p1, width = 10, height = 4, dpi = 300)


FeaturePlot(CD8T.clean, c("GNLY", "GZMB","GZMK", "CCR7", "CCR6", "TRDC", "CD158e1.KIR3DL1.AHS0211.pAbO", "KIR2DL1"), min.cutoff = 0, sort.cell = T)
FeaturePlot(CD8T.clean, "CD4", min.cutoff = 0, sort.cell = T)
FeaturePlot(CD8T.clean, "GNLY", min.cutoff = 0, sort.cell = T)
FeaturePlot(CD8T.clean, "GZMB", min.cutoff = 0, sort.cell = T)
FeaturePlot(CD8T.clean, "GZMK", min.cutoff = 0, sort.cell = T)
FeaturePlot(CD8T.clean, "CCR7", min.cutoff = 0, sort.cell = T)
FeaturePlot(CD8T.clean, "CCR6", min.cutoff = 0, sort.cell = T)
FeaturePlot(CD8T.clean, "KIR2DL1", min.cutoff = 0, sort.cell = T, split.by = "condition")
FeaturePlot(CD8T.clean, "CD158e1.KIR3DL1.AHS0211.pAbO", min.cutoff = 0, sort.cell = T, split.by = "condition")


FeaturePlot(CD8T.clean, "CD158e1.KIR3DL1.AHS0211.pAbO", min.cutoff = 0, max.cutoff = 20, split.by = "condition") #+ scale_color_gradient(low = "#D3D3D3",high = "#8C1515")
KIR_cluster <- subset(CD8T.clean, cluster == "GZMB+KIR+")
FeaturePlot(KIR_cluster, "CD158e1.KIR3DL1.AHS0211.pAbO", min.cutoff = 0, max.cutoff = 20, split.by = "condition") + scale_color_gradient(low = "#D3D3D3",high = "#8C1515")

FeaturePlot(KIR_cluster, "CD158e1.KIR3DL1.AHS0211.pAbO",max.cutoff = 20, split.by = "condition")+ scale_color_gradient(low = "#D3D3D3",high = "#8C1515")

table(CD8T.clean$cluster)
mypalette <- brewer.pal(7,"RdBu")

FeaturePlot(CD8T.clean, "CD158e1.KIR3DL1.AHS0211.pAbO", min.cutoff = 0, split.by = "condition", cols = "orange")


ggsave(filename = "./Paper_Figures/Myeloid/Myeloid_annotated.pdf",p2, width = 5, height = 4, dpi = 300)
table(CD8T.clean$cluster)

A <- subset(CD8T.clean, cluster %in% c("GZMB+KIR+"))
FeaturePlot(A, "CD158e1.KIR3DL1.AHS0211.pAbO", min.cutoff = 0, sort.cell = T, split.by = "condition")
FeaturePlot(A, "KIR2DL1", min.cutoff = 0, sort.cell = T, split.by = "condition")

A <- subset(CD8T.clean, cluster %in% c("GNLY+"))
FeaturePlot(A, "CD158e1.KIR3DL1.AHS0211.pAbO", min.cutoff = 0, sort.cell = T, split.by = "condition")
FeaturePlot(A, "KIR2DL1", min.cutoff = 0, sort.cell = T, split.by = "condition")


AB <- rownames(CD8T.clean)
AB[AB %like% "KIR"]

all.genes.rna <- rownames(pbmc)
all.genes.rna[all.genes.rna %like% "GAPDH"]

DefaultAssay(AB) <- "integrated"
FeaturePlot(A, "CD8A", min.cutoff = 0, sort.cell = T)
FeaturePlot(A, "CD4", min.cutoff = 0, sort.cell = T)
FeaturePlot(A, "CCR6", min.cutoff = 0, sort.cell = T)
FeaturePlot(A, "NKG7", min.cutoff = 0, sort.cell = T)
FeaturePlot(A, "CCR6", min.cutoff = 0, sort.cell = T)




#check CD8+KIR percetge
DefaultAssay(CD8T.clean) <- "antibody"
rownames(CD8T.clean)

antibody_number <- 31
rownames(CD8T.clean)[antibody_number]

CD8T.clean@assays$antibody[antibody_number] %>% as.numeric() %>% data.frame() %>% ggplot(aes(x=.)) + geom_histogram(bins = 50) + xlim(0.1,50) + ggtitle(rownames(CD8T.clean)[antibody_number])
p.KIR_ALL <- FeaturePlot(CD8T.clean, features = rownames(CD8T.clean)[antibody_number], max.cutoff = 100, min.cutoff = 20, order = T, split.by = "condition")&coord_fixed(ratio=1) 
p.KIR_ALL <- p.KIR_ALL + patchwork::plot_layout(ncol = 2, nrow = 2)
p.KIR_ALL
ggsave("./Paper_Figures/KIR/KIR3/CD8KIR3.UMAP.pdf",p.KIR_ALL, width = 7, height = 6, dpi=300)

#FeaturePlot(Myeloid.clean, features = rownames(Myeloid.clean)[antibody_number], max.cutoff = 100, min.cutoff = 10, order = T, split.by = "batch")&coord_fixed(ratio=1) 

antibody_number <- 29
rownames(CD8T.clean)[antibody_number]

CD8T.clean@assays$antibody[antibody_number] %>% as.numeric() %>% data.frame() %>% ggplot(aes(x=.)) + geom_histogram(bins = 50) + xlim(0.1,50) + ggtitle(rownames(CD8T.clean)[antibody_number])
p.CD158e <- FeaturePlot(CD8T.clean, features = rownames(CD8T.clean)[antibody_number], max.cutoff = 100, min.cutoff = 20, order = T, split.by = "condition")&coord_fixed(ratio=1) 
p.CD158e <- p.CD158e + patchwork::plot_layout(ncol = 2, nrow = 2)
p.CD158e

ggsave("./Paper_Figures/KIR/KIR2/CD8TKIR2.UMAP.pdf",p.CD158e, width = 7, height = 6, dpi=300)

batchtable <- table(CD8T.clean$subject, CD8T.clean$batch) %>% data.frame() %>% filter(Freq>0) %>% select(Var1, Var2) %>% unique()

barplot.df <- table(CD8T.clean$cluster, CD8T.clean$subject) %>% data.frame() 

barplot.df$condition <- "PTLD"
barplot.df$condition[str_detect(barplot.df$Var2, "HD", negate = FALSE)] <- "HD"
barplot.df$condition[str_detect(barplot.df$Var2, "RTH", negate = FALSE)] <- "RTH"
barplot.df$condition[str_detect(barplot.df$Var2, "PTLDN", negate = FALSE)] <- "PTLDN"
barplot.df$condition[str_detect(barplot.df$Var2, "HD9", negate = FALSE)] <- "CMV Infection"
barplot.df$condition[str_detect(barplot.df$Var2, "HD6", negate = FALSE)] <- "CMV Infection"
barplot.df$condition[str_detect(barplot.df$Var2, "HD10", negate = FALSE)] <- "CMV Infection"
barplot.df$timepoint <- "V1V5"
barplot.df$timepoint[str_detect(barplot.df$Var2, "V7", negate = FALSE)] <- "V7"

########Draw cluster percentage##########
p3 <- filter(barplot.df, timepoint != "V7")
p3 <- filter(p3,Var2 != "RTH2V5b1") %>% 
  filter(Var2 != "RTH5V3b5") %>%
  filter(Var2 != "PTLD2V1") %>%
  filter(condition != "PTLDN")


cluster_color = c(
  "GZMB+KIR+" = ggsci::pal_d3()(n=6)[1], 
  "GZMK+" = ggsci::pal_d3()(n=6)[2], 
  "CCR6+CD161+" = ggsci::pal_d3()(n=6)[3],
  "TCRgd+" = ggsci::pal_d3()(n=6)[4],
  "GNLY+" = ggsci::pal_d3()(n=6)[5],
  "Naive" = ggsci::pal_d3()(n=6)[6]
)
p4 <- ggplot(p3, aes(x=Var2,y=Freq, fill=Var1)) + geom_bar(position = "fill", stat = "identity") + facet_grid(.~ condition, scales = "free") + theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_y_continuous(labels = scales::percent)
p4 <- p4 + ylab("Cluster Percentage") + xlab("Individual and Timepoint")
p4 <- p4+ scale_y_continuous(labels = scales::percent)
p4
ggsave("./Paper_Figures/barplot_bycluster_CD8Tcluster.pdf",p4, width = 10, height = 7, dpi=300)

########################################

barplot.df <- barplot.df %>% group_by(Var2) %>% mutate(Percent = Freq/sum(Freq))
barplot.df <- merge(barplot.df, batchtable, by.x ="Var2", by.y="Var1")

CD8T.KIR2D<- subset(CD8T.clean, subset = `KIR.NKAT2.KIR2DL2-KIR2DL3-KIR2DS2.AHS0209.pAbO` > 20)
CD8T.KIR2D
CD8T.KIR3DL1 <- subset(CD8T.clean, subset = `CD158e1.KIR3DL1.AHS0211.pAbO` > 20)
CD8T.KIR3DL1

KIR2D.number <- table(CD8T.KIR2D$subject,CD8T.KIR2D$cluster) %>% data.frame()
KIR3DL1.number <- table(CD8T.KIR3DL1$subject,CD8T.KIR3DL1$cluster) %>% data.frame()

colnames(KIR2D.number)[3] <- "KIR2D"
colnames(KIR3DL1.number)[3] <- "KIR3DL1"

KIRpercent <- merge(barplot.df,KIR2D.number, by.x=c("Var2", "Var1"), by.y = c("Var1", "Var2")) %>% merge(KIR3DL1.number,by.x=c("Var2", "Var1"), by.y = c("Var1", "Var2"))
KIRpercent

KIRpercent <- KIRpercent %>% mutate(KIR2D_perc = KIR2D/Freq) %>% mutate(KIR3DL1_perc = KIR3DL1/Freq)

my_comparisons <- list( c("HD", "PTLD"), c("HD", "RTH"), c("PTLD", "RTH") )

KIRpercent$condition <- factor(KIRpercent$condition, levels = c("HD", "RTH","PTLD"))

barplot.df.wide <- reshape2::dcast(barplot.df, Var2 + condition + timepoint ~ Var1, value.var = "Percent")
barplot.df.wide
colnames(barplot.df.wide)[5] = "GNLY"
colnames(barplot.df.wide)[6] = "GZMB_KIR"
ggplot(barplot.df.wide, aes(x= GZMB_KIR, y= Naive)) + geom_point() + geom_smooth(method = "lm") + xlim(0,1) + ylim(0,0.3)


#KIR3D
#GZMB+KIR+
p.KIR3per.1 <- KIRpercent %>% filter(Var1 %in% c("GZMB+KIR+")) %>% filter(Freq > 10) %>% filter(condition !="PTLDN") %>% ggplot(aes(x=condition, y=KIR3DL1_perc)) + geom_jitter(pt.size=0.1) + geom_boxplot(alpha=0.3, outlier.alpha = 0)
p.KIR3per.1 <- p.KIR3per.1  + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method = "t.test",hide.ns = T) + theme_cowplot() + ggtitle("KIR3D in GZMB+KIR_CLUSTER")
p.KIR3per.1

t.test((filter(KIRpercent, Var1 == "GZMB+KIR+", condition =="RTH", Freq >10)$KIR3DL1_perc), (filter(KIRpercent, Var1 == "GZMB+KIR+", condition =="PTLD",Freq >10)$KIR3DL1_perc))

#GZMK+
p.KIR3per.2 <- KIRpercent %>% filter(Var1 != "GZMK+") %>% filter(Freq > 300) %>% filter(condition !="PTLDN") %>% ggplot(aes(x=condition, y=KIR3DL1_perc)) + geom_jitter(pt.size=0.1) + geom_boxplot(alpha=0.3,outlier.alpha = 0)
p.KIR3per.2 <- p.KIR3per.2  + stat_compare_means(comparisons = my_comparisons,label = "p.signif", method = "t.test", hide.ns = T) + theme_cowplot() + ggtitle("KIR3D in GZMK_CLUSTER")
p.KIR3per.2

#GNLY+
p.KIR3per.3 <- KIRpercent %>% filter(Var1 != "GNLY+") %>% filter(Freq > 300) %>% filter(condition !="PTLDN") %>% ggplot(aes(x=condition, y=KIR3DL1_perc)) + geom_jitter(pt.size=0.1) + geom_boxplot(alpha=0.3,outlier.alpha = 0)
p.KIR3per.3 <- p.KIR3per.3  + stat_compare_means(comparisons = my_comparisons,label = "p.signif", method = "t.test", hide.ns = T) + theme_cowplot() + ggtitle("KIR3D in GNLY_CLUSTER")
p.KIR3per.3

#TCRgd+
p.KIR3per.4 <- KIRpercent %>% filter(Var1 != "TCRgd+") %>% filter(Freq > 300) %>% filter(condition !="PTLDN") %>% ggplot(aes(x=condition, y=KIR3DL1_perc)) + geom_jitter(pt.size=0.1) + geom_boxplot(alpha=0.3,outlier.alpha = 0)
p.KIR3per.4 <- p.KIR3per.4  + stat_compare_means(comparisons = my_comparisons,label = "p.signif", method = "t.test", hide.ns = T) + theme_cowplot() + ggtitle("KIR3D in TCRgd_CLUSTER")
p.KIR3per.4

#CCR6+CD161+
p.KIR3per.5 <- KIRpercent %>% filter(Var1 != "CCR6+CD161+") %>% filter(Freq > 300) %>% filter(condition !="PTLDN") %>% ggplot(aes(x=condition, y=KIR3DL1_perc)) + geom_jitter(pt.size=0.1) + geom_boxplot(alpha=0.3,outlier.alpha = 0)
p.KIR3per.5 <- p.KIR3per.5  + stat_compare_means(comparisons = my_comparisons,label = "p.signif", method = "t.test", hide.ns = T) + theme_cowplot() + ggtitle("KIR3D in CCR6+CD161_CLUSTER")
p.KIR3per.5

#Naive
p.KIR3per.6 <- KIRpercent %>% filter(Var1 != "Naive") %>% filter(Freq > 300) %>% filter(condition !="PTLDN") %>% ggplot(aes(x=condition, y=KIR3DL1_perc)) + geom_jitter(pt.size=0.1) + geom_boxplot(alpha=0.3,outlier.alpha = 0)
p.KIR3per.6 <- p.KIR3per.6  + stat_compare_means(comparisons = my_comparisons,label = "p.signif", method = "t.test", hide.ns = T) + theme_cowplot() + ggtitle("KIR3D in Naive_CLUSTER")
p.KIR3per.6

p.KIR3 <- p.KIR3per.1 + p.KIR3per.2 + p.KIR3per.3 + p.KIR3per.4 + p.KIR3per.5 + p.KIR3per.6
p.KIR3
ggsave(filename = "./Paper_Figures/KIR/KIR3/CD8T_percent.compare.pdf", p.KIR3, width = 14, height = 8,dpi = 300)

#Combined
p.KIR3per.7 <- KIRpercent %>% filter(Var1 %in% c("GZMB+KIR+", "CCR6+CD161+", "GNLY+","TCRgd+", "GZMK+")) %>% filter(Freq > 10) %>% filter(condition !="PTLDN") %>% ggplot(aes(x=condition, y=KIR3DL1_perc)) + geom_jitter(pt.size=0.1) + geom_boxplot(alpha=0.3, outlier.alpha = 0)
p.KIR3per.7 <- p.KIR3per.7  + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method = "t.test",hide.ns = T) + theme_cowplot() + ggtitle("KIR3D in GZMB+KIR_CLUSTER")
p.KIR3per.7

t.test((filter(KIRpercent, Var1 == "GZMB+KIR+", condition =="RTH", Freq >10)$KIR3DL1_perc), (filter(KIRpercent, Var1 == "GZMB+KIR+", condition =="PTLD",Freq >10)$KIR3DL1_perc))
t.test((filter(KIRpercent, Var1 == "GZMB+KIR+", condition =="HD", Freq >10)$KIR3DL1_perc), (filter(KIRpercent, Var1 == "GZMB+KIR+", condition =="PTLD",Freq >10)$KIR3DL1_perc))
t.test((filter(KIRpercent, Var1 == "GZMB+KIR+", condition =="HD", Freq >10)$KIR3DL1_perc), (filter(KIRpercent, Var1 == "GZMB+KIR+", condition =="RTH",Freq >10)$KIR3DL1_perc))


#KIR2D
#GZMB+KIR+
p.KIR2per.1 <- KIRpercent %>% filter(Var1 %in% c("GZMB+KIR+")) %>% filter(Freq > 10) %>% filter(condition !="PTLDN") %>% ggplot(aes(x=condition, y=KIR2D_perc)) + geom_jitter(pt.size=0.1) + geom_boxplot(alpha=0.3, outlier.alpha = 0)
p.KIR2per.1 <- p.KIR2per.1  + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method = "t.test",hide.ns = T) + theme_cowplot() + ggtitle("KIR2D in GZMB+KIR_CLUSTER")
p.KIR2per.1

t.test((filter(KIRpercent, Var1 == "GZMB+KIR+", condition =="RTH", Freq >10)$KIR2D_perc), (filter(KIRpercent, Var1 == "GZMB+KIR+", condition =="PTLD",Freq >10)$KIR2D_perc))

#GZMK+
p.KIR2per.2 <- KIRpercent %>% filter(Var1 != "GZMB+KIR+") %>% filter(Freq > 300) %>% filter(condition !="PTLDN") %>% ggplot(aes(x=condition, y=KIR2D_perc)) + geom_jitter(pt.size=0.1) + geom_boxplot(alpha=0.3,outlier.alpha = 0)
p.KIR2per.2 <- p.KIR2per.2  + stat_compare_means(comparisons = my_comparisons,label = "p.signif", method = "t.test", hide.ns = T) + theme_cowplot() + ggtitle("KIR2D in GZMK_CLUSTER")
p.KIR2per.2

#TCRgd+
p.KIR2per.3 <- KIRpercent %>% filter(Var1 != "TCRgd+") %>% filter(Freq > 300) %>% filter(condition !="PTLDN") %>% ggplot(aes(x=condition, y=KIR2D_perc)) + geom_jitter(pt.size=0.1) + geom_boxplot(alpha=0.3,outlier.alpha = 0)
p.KIR2per.3 <- p.KIR2per.3  + stat_compare_means(comparisons = my_comparisons,label = "p.signif", method = "t.test", hide.ns = T) + theme_cowplot() + ggtitle("KIR2D in TCRgd_CLUSTER")
p.KIR2per.3

#GNLY+
p.KIR2per.4 <- KIRpercent %>% filter(Var1 != "GNLY+") %>% filter(Freq > 300) %>% filter(condition !="PTLDN") %>% ggplot(aes(x=condition, y=KIR2D_perc)) + geom_jitter(pt.size=0.1) + geom_boxplot(alpha=0.3,outlier.alpha = 0)
p.KIR2per.4 <- p.KIR2per.4  + stat_compare_means(comparisons = my_comparisons,label = "p.signif", method = "t.test", hide.ns = T) + theme_cowplot() + ggtitle("KIR2D in GNLY_CLUSTER")
p.KIR2per.4

#CCR6+CD161+
p.KIR2per.5 <- KIRpercent %>% filter(Var1 != "CCR6+CD161+") %>% filter(Freq > 300) %>% filter(condition !="PTLDN") %>% ggplot(aes(x=condition, y=KIR2D_perc)) + geom_jitter(pt.size=0.1) + geom_boxplot(alpha=0.3,outlier.alpha = 0)
p.KIR2per.5 <- p.KIR2per.5  + stat_compare_means(comparisons = my_comparisons,label = "p.signif", method = "t.test", hide.ns = T) + theme_cowplot() + ggtitle("KIR2D in CCR6+CD161_CLUSTER")
p.KIR2per.5

#Naive
p.KIR2per.6 <- KIRpercent %>% filter(Var1 != "Naive") %>% filter(Freq > 300) %>% filter(condition !="PTLDN") %>% ggplot(aes(x=condition, y=KIR2D_perc)) + geom_jitter(pt.size=0.1) + geom_boxplot(alpha=0.3,outlier.alpha = 0)
p.KIR2per.6 <- p.KIR2per.6  + stat_compare_means(comparisons = my_comparisons,label = "p.signif", method = "t.test", hide.ns = T) + theme_cowplot() + ggtitle("KIR2D in Naive_CLUSTER")
p.KIR2per.6

p.KIR2 <- p.KIR2per.1 + p.KIR2per.2 + p.KIR2per.3 + p.KIR2per.4 + p.KIR2per.5 + p.KIR2per.6
p.KIR2
ggsave(filename = "./Paper_Figures/KIR/KIR2/CD8T_percent.compare.pdf", p.KIR2, width = 10, height = 8,dpi = 300)


###for total CD8T###
KIR.CD8T.subset <- subset(CD8T.clean, CD158e1.KIR3DL1.AHS0211.pAbO > 20)
CD8T.clean$condition <- factor(CD8T.clean$condition,levels=c("HD", "RTH", "PTLD"))
kir.cd8.number <- data.frame(table(KIR.CD8T.subset$subject, KIR.CD8T.subset$cluster))
cd8.number <- data.frame(table(CD8T.clean$subject,CD8T.clean$cluster))

kirratio <- merge(kir.cd8.number,cd8.number, by=c("Var1","Var2"))
kirratio$Condition <- "HD"
kirratio$Condition[str_detect(kirratio$Var1, pattern = "PTLDN")] <- "PTLDN"
kirratio$Condition[str_detect(kirratio$Var1, pattern = "PTLD")] <- "PTLD"
kirratio$Condition[str_detect(kirratio$Var1, pattern = "RTH")] <- "RTH"

kirratio$Condition[str_detect(kirratio$Var1, pattern = "HD9")] <- "CMV positive"
kirratio$Condition[str_detect(kirratio$Var1, pattern = "HD6")] <- "CMV positive"
kirratio$Condition[str_detect(kirratio$Var1, pattern = "HD10")] <- "CMV positive"

kirratio$ratio <- kirratio$Freq.x/kirratio$Freq.y
p.kircd8 <- filter(kirratio, Condition != "PTLDN") %>% ggplot(aes(x=Condition, y=ratio)) + geom_boxplot() + ggtitle("% KIR3DL1 Positive CD8 T cells")
#p.kircd8
p.kircd9 <- filter(kirratio, Condition != "PTLDN") %>% ggplot(aes(x=Condition, y=ratio)) + geom_point(size=2) + ggtitle("% KIR3DL1 Positive CD8 T cells")
p.kircd9 <- p.kircd9 + scale_y_continuous(labels = scales::percent) #+ geom_errorbar(aes(ymin=(mean(ratio)-sd(ratio)), ymax=(mean(ratio)+sd(ratio))), color="red",width=.2)
p.kircd9 <- p.kircd9 + stat_summary(fun = median, fun.min = median, fun.max = median,geom = "crossbar", width = 0.5) + theme_cowplot() + ylab("Percentage")
#p.kircd9 <- p.kircd9 + stat_summary(fun = median, fun.min = median, fun.max = median, width = 0.5) + theme_cowplot() + ylab("Percentage")
p.kircd9 <- p.kircd9 + theme(axis.text.x = element_text(angle = 45, hjust=1)) + scientific_theme
p.kircd9
