#Single cell Lyme disease post merge
#Author: Xin Chen, Ph.D.
#Date Created: 022823
#Cleanup myeloid cluster after reclusters

#set up working directory
#setwd("/Users/xzhou7/Library/CloudStorage/Box-Box/Xin.Chen.Shareable/R7_NR")

setwd(dir = "/Users/xzhou7/Library/CloudStorage/Dropbox/lyme_disease/R7_NR")
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

# load("./Data/Step4.1_022823_recluster_Myeloid.integrated_final.RData")
# 
# table(Myeloid.integrated$integrated_snn_res.0.3)
# 
# CL0 <- subset(Myeloid.integrated, seurat_clusters == 0)
# cell.id.0 <- CellSelector(DimPlot(CL0))
# 
# CL1 <- subset(Myeloid.integrated, seurat_clusters == 1)
# cell.id.1 <- CellSelector(DimPlot(CL1))
# 
# CL2 <- subset(Myeloid.integrated, seurat_clusters == 2)
# cell.id.2 <- CellSelector(DimPlot(CL2))
# 
# CL3 <- subset(Myeloid.integrated, seurat_clusters == 3)
# cell.id.3 <- CellSelector(DimPlot(CL3))
# 
# CL4 <- subset(Myeloid.integrated, seurat_clusters == 4)
# cell.id.4 <- CellSelector(DimPlot(CL4))
# 
# CL5 <- subset(Myeloid.integrated, seurat_clusters == 5)
# cell.id.5 <- CellSelector(DimPlot(CL5))
# 
# CL6 <- subset(Myeloid.integrated, seurat_clusters == 6)
# cell.id.6 <- CellSelector(DimPlot(CL6))
# 
# CL7 <- subset(Myeloid.integrated, seurat_clusters == 7)
# cell.id.7 <- CellSelector(DimPlot(CL7))
# 
# CL8 <- subset(Myeloid.integrated, seurat_clusters == 8)
# cell.id.8 <- CellSelector(DimPlot(CL8))
# 
# 
# subset.id <- c(cell.id.0,cell.id.1, cell.id.2, cell.id.3, cell.id.4, cell.id.5,
#                cell.id.6, cell.id.7, cell.id.8)
# 
# save(subset.id, file = "./Data/Step4_022823_myeloid_subset.id.RData")
# 
# Myeloid.clean <- subset(Myeloid.integrated, cells = subset.id)
# save(Myeloid.clean, file = "./Data/Step4_022823_Myeloid.cluster_clean.RData")

load(file = "./Data/Step4_022823_Myeloid.cluster_clean.RData")
Myeloid.clean

p1 <- DimPlot(Myeloid.clean, group.by = "seurat_clusters", label = T)
p1

#ggsave(filename = "./Paper_Figures/Myeloid/Myeloid.pdf",p1, width = 5, height = 4, dpi = 300)

#find markers for cluster
pbmc.markers <- FindAllMarkers(Myeloid.clean, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)

#write.csv(pbmc.markers,"./Results/Myeloid.clean_gene_markers.csv", row.names = F)

#pbmc.markers <-read.csv("./Results/Myeloid.clean_gene_markers.csv", header = T)

Myeloid.clean$cluster <- "S100A9_Mono"
Myeloid.clean$cluster[Myeloid.clean$integrated_snn_res.0.3 == 1] <- "HLA_Mono"
Myeloid.clean$cluster[Myeloid.clean$integrated_snn_res.0.3 == 2] <- "HypoInflam_Mono"
Myeloid.clean$cluster[Myeloid.clean$integrated_snn_res.0.3 == 3] <- "NAMPT_Mono"
Myeloid.clean$cluster[Myeloid.clean$integrated_snn_res.0.3 == 4] <- "CD16_Mono"
Myeloid.clean$cluster[Myeloid.clean$integrated_snn_res.0.3 == 5] <- "X5_Mono"
Myeloid.clean$cluster[Myeloid.clean$integrated_snn_res.0.3 == 6] <- "X6_Mono"
Myeloid.clean$cluster[Myeloid.clean$integrated_snn_res.0.3 == 7] <- "cDC_Mono"
Myeloid.clean$cluster[Myeloid.clean$integrated_snn_res.0.3 == 8] <- "pDC_Mono"

cluster_color = c(
  "S100A9_Mono" = ggsci::pal_d3()(n=9)[9], 
  "HLA_Mono" = ggsci::pal_d3()(n=9)[3], 
  "HypoInflam_Mono" = ggsci::pal_d3()(n=9)[4],
  "NAMPT_Mono" = ggsci::pal_d3()(n=9)[2],
  "CD16_Mono" = ggsci::pal_d3()(n=9)[6],
  "X5_Mono" = ggsci::pal_d3()(n=9)[5],
  "X6_Mono" = ggsci::pal_d3()(n=9)[8],
  "cDC_Mono" = ggsci::pal_d3()(n=9)[7],
  "pDC_Mono" = ggsci::pal_d3()(n=9)[1]
)

# FindMarkers(Myeloid.clean, ident.1= 5, ident.2 =6)
# 
# FeaturePlot(Myeloid.clean, "IL32", min=0, order=T, split.by = "condition")
# 
# FeaturePlot(Myeloid.clean, "NKG7", min=0, order=T, split.by = "condition")
# 
# FeaturePlot(Myeloid.clean, "GNLY", min=0, order=T, split.by = "condition")
# 
# FeaturePlot(Myeloid.clean, "CCL4", min=0, order=T, split.by = "condition")
# FeaturePlot(Myeloid.clean, "SELL", min=0, order=T, split.by = "batch")

p2 <- DimPlot(Myeloid.clean, group.by = "cluster", label = T) + scale_color_manual(values = cluster_color) + coord_fixed(ratio=1)
p2

#ggsave(filename = "./Paper_Figures/Myeloid/Myeloid_annotated.pdf",p2, width = 5, height = 4, dpi = 300)

barplot.df <- table(Myeloid.clean$cluster, Myeloid.clean$subject) %>% data.frame() 

barplot.df$condition <- "PTLD"
barplot.df$condition[str_detect(barplot.df$Var2, "HD", negate = FALSE)] <- "HD"
barplot.df$condition[str_detect(barplot.df$Var2, "RTH", negate = FALSE)] <- "RTH"
barplot.df$condition[str_detect(barplot.df$Var2, "PTLDN", negate = FALSE)] <- "PTLDN"

barplot.df$timepoint <- "V1V5"
barplot.df$timepoint[str_detect(barplot.df$Var2, "V7", negate = FALSE)] <- "V7"

barplot.df$Var1 <- factor(barplot.df$Var1, levels = c("CD16_Mono","S100A9_Mono", "HLA_Mono", "NAMPT_Mono", "cDC_Mono","pDC_Mono", "X6_Mono", "X5_Mono", "HypoInflam_Mono" ))

p3 <- filter(barplot.df, timepoint != "V7") %>% 
  filter(condition != "PTLDN") %>% 
  ggplot(aes(x=Var2,y=Freq, fill=Var1)) + geom_bar(position = "fill", stat = "identity") + scale_fill_manual(values=cluster_color) + facet_grid(.~ condition, scales = "free") + theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
p3 <- p3 + ylab("Cluster Ratio") + xlab("Individual and Timepoint")
p3

#ggsave("./Paper_Figures/barplot_bycluster.pdf",p3, width = 10, height = 7, dpi=300)

barplot.df <- barplot.df %>% group_by(Var2) %>% mutate(Percent = Freq/sum(Freq))

barplot.df.wide <- reshape2::dcast(barplot.df, Var2 + condition + timepoint ~ Var1, value.var = "Percent")
barplot.df.wide

ggplot(barplot.df.wide, aes(x=HypoInflam_Mono, y= CD16_Mono)) + geom_point() + geom_smooth(method = "lm") + xlim(0,1) + ylim(0,0.3)

# #https://briatte.github.io/ggcorr/
# barplot.df.wide %>% select(-Var2, -condition, -timepoint) %>% ggcorr(geom = "circle",low = "blue", mid = "white", high = "red")

library(ggcorrplot)
#http://www.sthda.com/english/wiki/ggcorrplot-visualization-of-a-correlation-matrix-using-ggplot2 
p4.1 <- barplot.df.wide %>% select(-Var2, -condition, -timepoint) %>% cor(method="spearman") %>% ggcorrplot(hc.order = TRUE, outline.col = "white", type="lower",colors = c("#006B81", "#F4F4F4", "#E04F39")) + ggtitle("correlation by all condition")
p4.1

p4.2 <- barplot.df.wide %>% filter(condition == "HD") %>% select(-Var2, -condition, -timepoint) %>% cor(method="spearman") %>% ggcorrplot(hc.order = TRUE, outline.col = "white", type="lower",colors = c("#006B81", "#F4F4F4", "#E04F39")) + ggtitle("correlation by HD")
p4.2

p4.3 <- barplot.df.wide %>% filter(condition == "RTH") %>% select(-Var2, -condition, -timepoint) %>% cor(method="spearman") %>% ggcorrplot(hc.order = TRUE, outline.col = "white", type="lower",colors = c("#006B81", "#F4F4F4", "#E04F39")) + ggtitle("correlation by RTH")
p4.3

p4.4 <- barplot.df.wide %>% filter(condition == "PTLD") %>% select(-Var2, -condition, -timepoint) %>% cor(method="spearman") %>% ggcorrplot(hc.order = TRUE, outline.col = "white", type="lower",colors = c("#006B81", "#F4F4F4", "#E04F39")) + ggtitle("correlation by PTLD")
p4.4

p4 <- p4.2 + p4.3 + p4.4 +  plot_layout(guides = "collect")
#ggsave(filename = "./Paper_Figures/Mono_Internal_Correlation.pdf", p4, width = 11, height = 5, dpi=300)

Idents(Myeloid.clean) <- "cluster"

#https://www.frontiersin.org/articles/10.3389/fimmu.2019.02035/full 
#https://www.cell.com/immunity/pdf/S1074-7613(19)30334-6.pdf 
#cd48: https://pubmed.ncbi.nlm.nih.gov/16148114/ 
#https://www-science-org.stanford.idm.oclc.org/doi/10.1126/science.aah4573 

#https://www.pnas.org/doi/pdf/10.1073/pnas.2002476117
# NK-Like Monocytes:
# "GZMA", "GNLY", "KLRB1", "NKG7"

pdot <- DotPlot(Myeloid.clean, features = c("S100A9","S100A10","S100A12", "SELL","HLA.DRA","HLA.DPA1","HLA.DMA","CD163", "CD83","CCL3","CCL4", 
                                            "IL1B","CXCL8","IFIT3","NAMPT", "SLC2A3","CLEC4E","CD48","ITGB1","FCGR3A","CX3CR1","IL32", "IL7R","CCL5","CD68","CD63","FCER1A","CD74",
                                            "GZMB","CD4") ,cols = c("lightgrey", "#820000"))
pdot <- pdot + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))  
pdot

VariableFeatures(Myeloid.clean)
#ggsave(filename = "./Paper_Figures/dot_marker_mono.pdf",pdot, width = 9, height = 5, dpi = 300)

#antibodies
DefaultAssay(Myeloid.clean) <- "antibody"

#known antibody that are affected by batch: CCR7.CCR7.AHS0273.pAbO, CD127.IL7R.AHS0028.pAbO
rownames(Myeloid.clean)

antibody_number <- 31
rownames(Myeloid.clean)[antibody_number]

Myeloid.clean@assays$antibody[antibody_number] %>% as.numeric() %>% data.frame() %>% ggplot(aes(x=.)) + geom_histogram(bins = 50) + xlim(0.1,50) + ggtitle(rownames(Myeloid.clean)[antibody_number])
p.KIR_ALL <- FeaturePlot(Myeloid.clean, features = rownames(Myeloid.clean)[antibody_number], max.cutoff = 100, min.cutoff = 20, order = T, split.by = "condition")&coord_fixed(ratio=1) 
p.KIR_ALL <- p.KIR_ALL + patchwork::plot_layout(ncol = 2, nrow = 2)
p.KIR_ALL
#ggsave("./Paper_Figures/KIR/KIR3/KIR3.UMAP.pdf",p.KIR_ALL, width = 7, height = 6, dpi=300)

#FeaturePlot(Myeloid.clean, features = rownames(Myeloid.clean)[antibody_number], max.cutoff = 100, min.cutoff = 10, order = T, split.by = "batch")&coord_fixed(ratio=1) 

antibody_number <- 29
rownames(Myeloid.clean)[antibody_number]

Myeloid.clean@assays$antibody[antibody_number] %>% as.numeric() %>% data.frame() %>% ggplot(aes(x=.)) + geom_histogram(bins = 50) + xlim(0.1,50) + ggtitle(rownames(Myeloid.clean)[antibody_number])
p.CD158e <- FeaturePlot(Myeloid.clean, features = rownames(Myeloid.clean)[antibody_number], max.cutoff = 100, min.cutoff = 20, order = T, split.by = "condition")&coord_fixed(ratio=1) 
p.CD158e <- p.CD158e + patchwork::plot_layout(ncol = 2, nrow = 2)
p.CD158e

ggsave("./Paper_Figures/KIR/KIR2/KIR2.UMAP.pdf",p.CD158e, width = 7, height = 6, dpi=300)

batchtable <- table(Myeloid.clean$subject, Myeloid.clean$batch) %>% data.frame() %>% filter(Freq>0) %>% select(Var1, Var2) %>% unique()

barplot.df <- merge(barplot.df, batchtable, by.x ="Var2", by.y="Var1")

Myeloid.KIR2D<- subset(Myeloid.clean, subset = `KIR.NKAT2.KIR2DL2-KIR2DL3-KIR2DS2.AHS0209.pAbO` > 20)
Myeloid.KIR2D
Myeloid.KIR3DL1 <- subset(Myeloid.clean, subset = `CD158e1.KIR3DL1.AHS0211.pAbO` > 20)
Myeloid.KIR3DL1

KIR2D.number <- table(Myeloid.KIR2D$subject,Myeloid.KIR2D$cluster) %>% data.frame()
KIR3DL1.number <- table(Myeloid.KIR3DL1$subject,Myeloid.KIR3DL1$cluster) %>% data.frame()

colnames(KIR2D.number)[3] <- "KIR2D"
colnames(KIR3DL1.number)[3] <- "KIR3DL1"

KIRpercent <- merge(barplot.df,KIR2D.number, by.x=c("Var2", "Var1"), by.y = c("Var1", "Var2")) %>% merge(KIR3DL1.number,by.x=c("Var2", "Var1"), by.y = c("Var1", "Var2"))
KIRpercent

KIRpercent <- KIRpercent %>% mutate(KIR2D_perc = KIR2D/Freq) %>% mutate(KIR3DL1_perc = KIR3DL1/Freq)

my_comparisons <- list( c("HD", "PTLD"), c("HD", "RTH"), c("PTLD", "RTH") )

KIRpercent$condition <- factor(KIRpercent$condition, levels = c("HD", "RTH","PTLD"))

#KIR3D
#HypoInflam_Mono
p.KIR3per.1 <- KIRpercent %>% filter(Var1 %in% c("HypoInflam_Mono")) %>% filter(Freq > 10) %>% filter(condition !="PTLDN") %>% ggplot(aes(x=condition, y=KIR3DL1_perc)) + geom_jitter(pt.size=0.1) + geom_boxplot(alpha=0.3, outlier.alpha = 0)
p.KIR3per.1 <- p.KIR3per.1  + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method = "t.test",hide.ns = T) + theme_cowplot() + ggtitle("KIR3D in HyperInflam Mono")
p.KIR3per.1

t.test((filter(KIRpercent, Var1 == "HypoInflam_Mono", condition =="RTH", Freq >10)$KIR3DL1_perc), (filter(KIRpercent, Var1 == "HypoInflam_Mono", condition =="PTLD",Freq >10)$KIR3DL1_perc))

#HLA_Mono
p.KIR3per.2 <- KIRpercent %>% filter(Var1 != "HypoInflam_Mono") %>% filter(Freq > 300) %>% filter(condition !="PTLDN") %>% ggplot(aes(x=condition, y=KIR3DL1_perc)) + geom_jitter(pt.size=0.1) + geom_boxplot(alpha=0.3,outlier.alpha = 0)
p.KIR3per.2 <- p.KIR3per.2  + stat_compare_means(comparisons = my_comparisons,label = "p.signif", method = "t.test", hide.ns = T) + theme_cowplot() + ggtitle("KIR3D in none-HyperInflam Mono")
p.KIR3per.2

p.KIR3 <- p.KIR3per.1 + p.KIR3per.2
p.KIR3
ggsave(filename = "./Paper_Figures/KIR/KIR3/percent.compare.pdf", p.KIR3, width = 6, height = 5,dpi = 300)


#KIR2D
#HypoInflam_Mono
p.KIR2per.1 <- KIRpercent %>% filter(Var1 %in% c("HypoInflam_Mono")) %>% filter(Freq > 10) %>% filter(condition !="PTLDN") %>% ggplot(aes(x=condition, y=KIR2D_perc)) + geom_jitter(pt.size=0.1) + geom_boxplot(alpha=0.3, outlier.alpha = 0)
p.KIR2per.1 <- p.KIR2per.1  + stat_compare_means(comparisons = my_comparisons,label = "p.signif",method = "t.test",hide.ns = T) + theme_cowplot() + ggtitle("KIR2D in HyperInflam Mono")
p.KIR2per.1

t.test((filter(KIRpercent, Var1 == "HypoInflam_Mono", condition =="RTH", Freq >10)$KIR2D_perc), (filter(KIRpercent, Var1 == "HypoInflam_Mono", condition =="PTLD",Freq >10)$KIR2D_perc))

#HLA_Mono
p.KIR2per.2 <- KIRpercent %>% filter(Var1 != "HypoInflam_Mono") %>% filter(Freq > 300) %>% filter(condition !="PTLDN") %>% ggplot(aes(x=condition, y=KIR2D_perc)) + geom_jitter(pt.size=0.1) + geom_boxplot(alpha=0.3,outlier.alpha = 0)
p.KIR2per.2 <- p.KIR2per.2  + stat_compare_means(comparisons = my_comparisons,label = "p.signif", method = "t.test", hide.ns = T) + theme_cowplot() + ggtitle("KIR2D in none-HyperInflam Mono")
p.KIR2per.2

p.KIR2 <- p.KIR2per.1 + p.KIR2per.2
p.KIR2
ggsave(filename = "./Paper_Figures/KIR/KIR2/percent.compare.pdf", p.KIR2, width = 6, height = 5,dpi = 300)
# 
# filter(KIRpercent, Var1 == "HypoInflam_Mono", condition =="RTH")
# 
# chisq.test(cbind(c(22,2572), c(82,4450)))
# 
# sum(filter(KIRpercent, Var1 == "HypoInflam_Mono", condition =="HD", Freq >10)$KIR3DL1)
# sum(filter(KIRpercent, Var1 == "HypoInflam_Mono", condition =="HD", Freq >100)$Freq)
# 
# sum(filter(KIRpercent, Var1 == "HypoInflam_Mono", condition =="RTH", Freq >10)$KIR3DL1)
# sum(filter(KIRpercent, Var1 == "HypoInflam_Mono", condition =="RTH", Freq >100)$Freq)
# 
# sum(filter(KIRpercent, Var1 == "HypoInflam_Mono", condition =="PTLD",Freq >10)$KIR3DL1)
# sum(filter(KIRpercent, Var1 == "HypoInflam_Mono", condition =="PTLD",Freq >100)$Freq)
# 
# sum(filter(KIRpercent, Var1 != "HypoInflam_Mono", condition =="HD", Freq >10)$KIR3DL1)
# sum(filter(KIRpercent, Var1 != "HypoInflam_Mono", condition =="HD", Freq >100)$Freq)
# 
# sum(filter(KIRpercent, Var1 != "HypoInflam_Mono", condition =="RTH", Freq >10)$KIR3DL1)
# sum(filter(KIRpercent, Var1 != "HypoInflam_Mono", condition =="RTH", Freq >100)$Freq)
# 
# sum(filter(KIRpercent, Var1 != "HypoInflam_Mono", condition =="PTLD",Freq >10)$KIR3DL1)
# sum(filter(KIRpercent, Var1 != "HypoInflam_Mono", condition =="PTLD",Freq >100)$Freq)
# 
# 
# p.rth.IM <- data.frame(group=c("KIR3DL1+", "KIR3DL1-"),value=c(22, 2572-22)) %>% ggplot(aes(x="", y=value,fill=group)) + geom_bar(stat="identity", width=1) + coord_polar("y", start=0) + theme_cowplot() + scale_fill_manual(values = c("grey", "blue"))
# p.rth.IM<- p.rth.IM + ggtitle("RTH HyperInflam Mono KIR3DL1 %")
# p.rth.IM
# 
# p.ptld.IM <- data.frame(group=c("KIR3DL1+", "KIR3DL1-"),value=c(82, 4450-82)) %>% ggplot(aes(x="", y=value,fill=group)) + geom_bar(stat="identity", width=1) + coord_polar("y", start=0) + theme_cowplot() + scale_fill_manual(values = c("grey", "blue"))
# p.ptld.IM<- p.ptld.IM + ggtitle("PTLD HyperInflam Mono KIR3DL1 %")
# p.ptld.IM
# 
# p.hd.NIM <- data.frame(group=c("KIR3DL1+", "KIR3DL1-"),value=c(162, 8292-162)) %>% ggplot(aes(x="", y=value,fill=group)) + geom_bar(stat="identity", width=1) + coord_polar("y", start=0) + theme_cowplot() + scale_fill_manual(values = c("grey", "blue"))
# p.hd.NIM<- p.hd.NIM + ggtitle("HD NONE HyperInflam Mono KIR3DL1 %")
# p.hd.NIM
# 
# p.rth.NIM <- data.frame(group=c("KIR3DL1+", "KIR3DL1-"),value=c(266, 22375-266)) %>% ggplot(aes(x="", y=value,fill=group)) + geom_bar(stat="identity", width=1) + coord_polar("y", start=0) + theme_cowplot() + scale_fill_manual(values = c("grey", "blue"))
# p.rth.NIM<- p.rth.NIM + ggtitle("RTH NONE HyperInflam Mono KIR3DL1 %")
# p.rth.NIM
# 
# p.ptld.NIM <- data.frame(group=c("KIR3DL1+", "KIR3DL1-"),value=c(76, 6566-76)) %>% ggplot(aes(x="", y=value,fill=group)) + geom_bar(stat="identity", width=1) + coord_polar("y", start=0) + theme_cowplot() + scale_fill_manual(values = c("grey", "blue"))
# p.ptld.NIM<- p.ptld.NIM + ggtitle("PTLD NONE HyperInflam Mono KIR3DL1 %")
# p.ptld.NIM
# 
# pKIR.bar <- p.rth.IM + p.ptld.IM + p.hd.NIM + p.rth.NIM + p.ptld.NIM
# ggsave(filename = "Paper_Figures/KIR/KIR3/KIR.Pie.pdf",pKIR.bar, width = 12, height = 8, dpi=300)

FeaturePlot(Myeloid.clean, features = "KIR.NKAT2.KIR2DL2-KIR2DL3-KIR2DS2.AHS0209.pAbO", max.cutoff = 500, min.cutoff = 10, order = T, split.by = "condition")&coord_fixed(ratio=1) 


DefaultAssay(Myeloid.clean) <- "RNA"
FeatureScatter(Myeloid.clean, feature1 = "nCount_RNA",feature2 = "nFeature_RNA")





