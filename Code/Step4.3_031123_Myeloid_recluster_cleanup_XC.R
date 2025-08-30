#Single cell Lyme disease post merge
#Author: Xin Chen, Ph.D.
#Date Created: 022823
#Cleanup myeloid cluster after reclusters
#Updated 120624

#set up working directory
#setwd("C:/Users/zhoux/Box/Xin.Chen.Shareable/R7_NR")

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
library(stringr)
library(EnhancedVolcano)
library(tidygraph)
library(igraph)
library(ggraph)

#mac
setwd(dir = "~/Library/CloudStorage/Dropbox/lyme_disease/R7_NR/")



load(file = "./Data/Step4_022823_Myeloid.cluster_clean.RData")

#Myeloid.clean
#DefaultAssay(Myeloid.clean) <- "antibody"
#rownames(Myeloid.clean[["antibody"]])
#DefaultAssay(Myeloid.clean) <- "antibody"

#DefaultAssay(Myeloid.clean) <- "antibody"
#library(ggplot2)
# Generate FeaturePlot for CD14
#DefaultAssay(Myeloid.clean) <- "antibody"
#p <- FeaturePlot(
#  Myeloid.clean,
#  features = "CD14.MPHIP9.CD14.AHS0037.pAbO", 
#  reduction = "umap",
#  cols = c("lightgray", "#8C1515"),
#  min.cutoff = "0",   # sets lower bound of color scale
#  max.cutoff = "1000"    # set to the upper range you want
#)

# Save as PDF
#ggsave("CD14_UMAP.pdf", plot = p, width = 6, height = 5)


p1 <- DimPlot(Myeloid.clean, group.by = "seurat_clusters", label = T, split.by = "condition")
p1

#ggsave(filename = "./Paper_Figures/Myeloid/Myeloid.pdf",p1, width = 5, height = 4, dpi = 300)



FeaturePlot(Myeloid.clean, "CD14", min.cutoff = 0)

Myeloid.clean$cluster <- "S100A9_Mono"
Myeloid.clean$cluster[Myeloid.clean$integrated_snn_res.0.3 == 1] <- "HLA_Mono"
Myeloid.clean$cluster[Myeloid.clean$integrated_snn_res.0.3 == 2] <- "HyperInflam_Mono"
Myeloid.clean$cluster[Myeloid.clean$integrated_snn_res.0.3 == 3] <- "NAMPT_Mono"
Myeloid.clean$cluster[Myeloid.clean$integrated_snn_res.0.3 == 4] <- "CD16_Mono"
Myeloid.clean$cluster[Myeloid.clean$integrated_snn_res.0.3 == 5] <- "IL32_Mono"
Myeloid.clean$cluster[Myeloid.clean$integrated_snn_res.0.3 == 6] <- "LEF1_Mono"
Myeloid.clean$cluster[Myeloid.clean$integrated_snn_res.0.3 == 7] <- "cDC"
Myeloid.clean$cluster[Myeloid.clean$integrated_snn_res.0.3 == 8] <- "pDC"

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

p2 <- DimPlot(Myeloid.clean, group.by = "cluster", label = T) + scale_color_manual(values = Mono_color) + coord_fixed(ratio=1)
p2

Idents(Myeloid.clean) <- "cluster"

Myeloid.clean_keep <- sample(colnames(Myeloid.clean), 20000)
Myeloid.clean.subset <- subset(Myeloid.clean, cells = Myeloid.clean_keep)

p_umap <- DimPlot(Myeloid.clean.subset, label =F, raster=FALSE) & coord_equal() 
#p_umap <- p_umap + scale_color_manual(values = my_palette) + scale_fill_manual(values = my_palette_alpha_0.7)
p_umap

# #annotate label position
celltype_position=Myeloid.clean@reductions$umap@cell.embeddings %>% as.data.frame() %>%
  cbind(celltype=Myeloid.clean@meta.data$cluster) %>%
  group_by(celltype) %>% dplyr::summarise(UMAP_1 = median(UMAP_1),UMAP_2 = median(UMAP_2))

p_umap2 <- p_umap+ geom_label_repel(data = celltype_position,aes(x=UMAP_1, y=UMAP_2, label=celltype, color=celltype),
                                    fontface="bold", point.padding=unit(0.1, "lines"), alpha=0.95)+theme(legend.position = "none") +
  scale_color_manual(values = Mono_color)

p_umap2 <- p_umap+ theme(legend.position = "none") + scale_color_manual(values = Mono_color)

p_umap2
#######if want a white background for the label######

# Calculate median UMAP position per cluster
celltype_position <- Myeloid.clean@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  cbind(celltype = Myeloid.clean@meta.data$cluster) %>%
  group_by(celltype) %>%
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")

# Build UMAP plot with white-background cluster labels
p_umap2 <- p_umap +
  geom_label_repel(
    data = celltype_position,
    aes(x = UMAP_1, y = UMAP_2, label = celltype, color = celltype),
    fontface = "bold",
    fill = "white",           # White background for label
    point.padding = unit(0.1, "lines"),
    alpha = 0.95,
    show.legend = FALSE
  ) +
  theme(legend.position = "none") +
  scale_color_manual(values = Mono_color)

# Show the plot
p_umap2

#ggsave(filename = "/Users/xinchen/Library/CloudStorage/Dropbox/lyme_disease/Manuscript/Figures/Figure3-autoimmune monocytes/Umap_label_XC.pdf",p_umap2, width = 5, height = 4, dpi = 300)
####################################





table(Myeloid.clean$cluster)


#ggsave(filename = "/Users/xinchen/Library/CloudStorage/Dropbox/lyme_disease/Manuscript/Figures/Figure3-autoimmune monocytes/Umap_no_label_XC.pdf",p_umap2, width = 5, height = 4, dpi = 300)
#ggsave(filename = "/Users/xinchen/Library/CloudStorage/Dropbox/lyme_disease/Manuscript/Figures/Figure3-autoimmune monocytes/sig_Myeloid_annotated_bycondition.pdf", p2.2, width = 10, height = 5, dpi = 300)
### by count
count2.2 <- data.frame(table(Myeloid.clean$subject, Myeloid.clean$cluster))
total.bysubject <- data.frame(table(Myeloid.clean$subject))
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

#count2.2 <- subset(count2.2, Conditon %in% c("PTLD", "RTH") )
count2.2 <- subset(count2.2, Var2 %in% c("CD16_Mono","HLA_Mono","NAMP_Mono","S100A9_Mono","HyperInflam_Mono"))
count2.2$Var2 <- factor(count2.2$Var2, levels = c("CD16_Mono","HyperInflam_Mono","HLA_Mono","NAMP_Mono","S100A9_Mono"))


count2.2$conditon <- factor(count2.2$conditon, levels = c("HD", "PTLDN", "RTH", "PTLD"))

comparasions <- list(c("HD", "PTLDN"), c("RTH","PTLD"))

p2.2 <- ggplot(count2.2, aes(x=conditon, y=Percent, color=conditon)) + geom_boxplot(alpha=0.7, outlier.alpha = 0, size = 0.7)
p2.2 <- p2.2 + facet_wrap(.~Var2, scales = "free", nrow = 1) + theme_cowplot()
p2.2 <- p2.2 + geom_jitter(color="black", size=0.4, alpha=0.9)
p2.2 <- p2.2 + stat_compare_means(method="t.test", comparisons = comparasions) + scale_color_manual(values=condition_color) + XZ_flip_x_label()
p2.2 <- ggplot(count2.2, aes(x = conditon, y = Percent, color = conditon)) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.9) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0, width = 0.6) +
  facet_wrap(. ~ Var2, scales = "free", nrow = 1) +
  stat_compare_means(method = "t.test", comparisons = comparasions) +
  scale_color_manual(values = condition_color) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12),
    axis.text.x = element_text(size = 10, angle = 0, hjust = 0.5),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm")
  )

p2.2

ggsave(filename = "/Users/xinchen/Library/CloudStorage/Dropbox/lyme_disease/Manuscript/Figures/Figure3-autoimmune monocytes/sig_Myeloid_annotated_bycondition.pdf", p2.2, width = 10, height = 4, dpi = 300)


 ############# end
Myeloid.clean$condition_level <- factor(Myeloid.clean$condition,levels = c("HD", "PTLDN", "RTH","PTLD"))

# Setting a seed for reproducibility
set.seed(123) 
sampled_cells <- list()
for (i in levels(Myeloid.clean$condition_level)) {
  print(i)
  # Identify cells in the current condition
  cells_in_condition <- WhichCells(Myeloid.clean, expression = condition_level == i)
  sampled_cells[[i]] <- sample(cells_in_condition, size = 10000, replace = FALSE)
}
all_sampled_cells <- unlist(sampled_cells)
Myeloid.clean.sampled15000 <- subset(Myeloid.clean, cells = all_sampled_cells)
table(Myeloid.clean.sampled15000$condition_level)

p2.5 <- DimPlot(Myeloid.clean, group.by = "cluster", split.by = "condition_level") + scale_color_manual(values = Mono_color) + coord_fixed(ratio=1)
p2.5

ggsave(filename = "/Users/xinchen/Library/CloudStorage/Dropbox/lyme_disease/Manuscript/Figures/Figure3-autoimmune monocytes/sig_Myeloid_annotated_bycondition.pdf", p2.5, width = 30, height = 4, dpi = 300)

#by condition
table(Myeloid.clean$condition)

p2 <- DimPlot(Myeloid.clean, group.by = "cluster", label = T, split.by = "condition") + scale_color_manual(values = Mono_color) + coord_fixed(ratio=1)
p2


#annaotated Myeloid.clean
#save(Myeloid.clean, file = "./Data/Step4.3_032523_Myeloid.cluster_clean_rename_cluster.RData")

###############################################################################################
setwd(dir = "/Users/xinchen/Library/CloudStorage/Dropbox/lyme_disease/R7_NR/")
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


#source("./00_colorKey_XC.R")
#source("./Code/00_colorKey.R")
source("./Code/xztools.R")
source("./Code/00_colorKey_XC.R")
#reload
load(file = "./Data/Step4.3_032523_Myeloid.cluster_clean_rename_cluster.RData")
Myeloid.clean[[]]
table(Myeloid.clean$subject)
table(Myeloid.clean$cluster)
barplot.df <- table(Myeloid.clean$cluster, Myeloid.clean$subject) %>% data.frame() 

barplot.df$condition <- "PTLD"
barplot.df$condition[str_detect(barplot.df$Var2, "HD", negate = FALSE)] <- "HD"
barplot.df$condition[str_detect(barplot.df$Var2, "RTH", negate = FALSE)] <- "RTH"
barplot.df$condition[str_detect(barplot.df$Var2, "PTLDN", negate = FALSE)] <- "PTLDN"
barplot.df$condition <- factor(barplot.df$condition, levels = c("HD", "PTLDN", "RTH", "PTLD"))

table(barplot.df$condition)

barplot.df$timepoint <- "V1V5"
barplot.df$timepoint[str_detect(barplot.df$Var2, "V7", negate = FALSE)] <- "V7"

barplot.df$Var1 <- factor(barplot.df$Var1, levels = c("CD16_Mono","S100A9_Mono", "HLA_Mono", "NAMPT_Mono", "cDC_Mono","pDC_Mono", "X5_Mono", "X6_Mono", "HyperInflam_Mono"))
#barplot.df$Var2 <- factor(barplot.df$Var2, levels = c("HD1","HD2", "HD3","HD4", "HD5","HD6", "HD7","HD8","HD9","HD10", "HD11","HD12","PTLDN1", "PTLDN2", "PTLDN3","PTLDN4","PTLDN5","PTLDN6",
#                                                      "PTLDN7", "PTLDN8", "PTLD1V2","PTLD1V3", "PTLD1V5","PTLD2V2", "PTLD2V3","PTLD2V5", "PTLD2V7","PTLD3V2","PTLD3V3","PTLD3V5", "PTLD3V7","PTLD4V3","PTLD4V7",
#                                                     "RTH1V2","RTH1V3", "RTH1V5","RTH1V7b1","RTH1V7b5", "RTH2V2", "RTH2V5b1", "RTH2V5b5","RTH2V7","RTH3V3","RTH3V5","RTH4V2","RTH4V3","RTH4V5",
#                                                                                                    "RTH4V7","RTH5V2","RTH5V3b2","RTH5V3b5","RTH5V7","RTH6V3","RTH6V5","RTH6V7"))


barplot.df$Var2 <- factor(barplot.df$Var2, levels = c("HD1","HD2", "HD3","HD4", "HD5","HD6", "HD7","HD8","HD9","HD10", "HD11","HD12","PTLDN1", "PTLDN2", "PTLDN3","PTLDN4","PTLDN5","PTLDN6",
                                                      "PTLDN7", "PTLDN8", "PTLD1V2","PTLD1V3", "PTLD1V5","PTLD2V2", "PTLD2V3","PTLD2V5", "PTLD2V7","PTLD3V2","PTLD3V3","PTLD3V5", "PTLD3V7","PTLD4V3","PTLD4V7",
                                                      "RTH1V2","RTH1V3", "RTH1V5","RTH1V7b1","RTH1V7b5", "RTH2V2","RTH2V5b5","RTH2V7","RTH3V3","RTH3V5","RTH4V2","RTH4V3","RTH4V5",
                                                      "RTH4V7","RTH5V2","RTH5V3b2","RTH5V7","RTH6V3","RTH6V5","RTH6V7"))

table(barplot.df$condition)  # Does "PTLD" have any counts?
unique(barplot.df$condition)
levels(barplot.df$condition)
setdiff(levels(barplot.df$Var1), names(cluster_color))


#table(Myeloid.clean$subject)
p3 <- filter(barplot.df, timepoint != "V7") %>%
  filter(Var2 != "RTH2V5b1") %>% 
  filter(Var2 != "RTH5V3b5") %>%
  filter(Var2 != "PTLD2V1") %>%
  ggplot(aes(x=Var2,y=Freq, fill=Var1)) + geom_bar(position = "fill", stat = "identity") + scale_fill_manual(values=cluster_color) + facet_grid(.~ condition, scales = "free") + theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_y_continuous(labels = scales::percent)
p3 <- p3 + ylab("Cluster Percentage") + xlab("Individual and Timepoint")
p3 <- p3+ scale_y_continuous(labels = scales::percent)
p3

ggsave("../Manuscript/Figures/Figure3-autoimmune monocytes/barplot_bycluster_removeduplicate_noV7.pdf",p3, width = 8, height = 5, dpi=300)

#table(Myeloid.clean$subject)
p3 <- filter(barplot.df, Var2 != "RTH2V5b1") %>% 
  filter(Var2 != "RTH5V3b5") %>%
  filter(Var2 != "PTLD2V1") %>%
  ggplot(aes(x=Var2,y=Freq, fill=Var1)) + geom_bar(position = "fill", stat = "identity") + scale_fill_manual(values=cluster_color) + facet_grid(.~ condition, scales = "free") + theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_y_continuous(labels = scales::percent)
p3 <- p3 + ylab("Cluster Percentage") + xlab("Individual and Timepoint")
p3 <- p3+ scale_y_continuous(labels = scales::percent)
p3
#ggsave("../lyme_disease/Manuscript/Figures/Figure3/barplot_bycluster_removeduplicate.pdf",p3, width = 10, height = 7, dpi=300)


#barplot for HD and PTLDN
# barplot.df <- table(Myeloid.clean$cluster, Myeloid.clean$subject) %>% data.frame() 
# 
# barplot.df$condition <- "PTLD"
# barplot.df$condition[str_detect(barplot.df$Var2, "HD", negate = FALSE)] <- "HD"
# barplot.df$condition[str_detect(barplot.df$Var2, "RTH", negate = FALSE)] <- "RTH"
# barplot.df$condition[str_detect(barplot.df$Var2, "PTLDN", negate = FALSE)] <- "PTLDN"
# 
# barplot.df$timepoint <- "V1V5"
# barplot.df$timepoint[str_detect(barplot.df$Var2, "V7", negate = FALSE)] <- "V7"
# 
# table(Myeloid.clean$cluster)
# barplot.df <- filter(barplot.df, timepoint != "V7")
# 
# barplot.df <- filter(barplot.df, Var2 != "RTH2V5b1")
# barplot.df <- filter(barplot.df, Var2 != "RTH5V3b5")
# barplot.df <- filter(barplot.df,Var2 != "PTLD2V1")
# 
# barplot.df <- filter(barplot.df, condition != "PTLD")
# barplot.df <- filter(barplot.df, condition != "RTH")

#barplot.df$Var1 <- factor(barplot.df$Var1, levels = c("CD16_Mono","S100A9_Mono", "HLA_Mono", "NAMPT_Mono", "cDC_Mono","pDC_Mono", "X6_Mono", "X5_Mono", "HypoInflam_Mono" ))
#barplot.df$Var2 <- factor(barplot.df$Var2, levels = c("HD1","HD2", "HD3","HD4", "HD5","HD6", "HD7","HD8","HD9","HD10", "HD11","HD12",
#                                                      "PTLD1V2","PTLD1V3", "PTLD1V5","PTLD2V2", "PTLD2V3","PTLD2V5", "PTLD2V7","PTLD3V2","PTLD3V3","PTLD3V5", "PTLD3V7","PTLD4V3","PTLD4V7",
#                                                      "RTH1V2","RTH1V3", "RTH1V5","RTH1V7b1","RTH1V7b5", "RTH2V2", "RTH2V5b1", "RTH2V5b5","RTH2V7","RTH3V3","RTH3V5","RTH4V2","RTH4V3","RTH4V5","RTH4V7","RTH5V2","RTH5V3b2","RTH5V3b5","RTH5V7","RTH6V3","RTH6V5","RTH6V7"))
#table(Myeloid.clean$subject)
# p3 <- filter(barplot.df, timepoint != "V7") %>%
#   filter(Var2 != "RTH2V5b1") %>% 
#   filter(Var2 != "RTH5V3b5") %>%
#   filter(Var2 != "PTLD2V1") %>%
#   ggplot(aes(x=Var2,y=Freq, fill=Var1)) + geom_bar(position = "fill", stat = "identity") + scale_fill_manual(values=cluster_color) + facet_grid(.~ condition, scales = "free") + theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_y_continuous(labels = scales::percent)
# p3 <- p3 + ylab("Cluster Percentage") + xlab("Individual and Timepoint")
# p3 <- p3+ scale_y_continuous(labels = scales::percent)
# p3
# ggsave("./Paper_Figures/barplot_bycluster_removeduplicate.pdf",p3, width = 10, height = 7, dpi=300)

barplot.df <- barplot.df %>% group_by(Var2) %>% mutate(Percent = Freq/sum(Freq))

barplot.df.wide <- reshape2::dcast(barplot.df, Var2 + condition + timepoint ~ Var1, value.var = "Percent")
barplot.df.wide

ggplot(barplot.df.wide, aes(x=HyperInflam_Mono, y= CD16_Mono)) + geom_point() + geom_smooth(method = "lm") + xlim(0,1) + ylim(0,0.3)

# #https://briatte.github.io/ggcorr/
# barplot.df.wide %>% select(-Var2, -condition, -timepoint) %>% ggcorr(geom = "circle",low = "blue", mid = "white", high = "red")

library(ggcorrplot)
#install.packages("ggcorrplot")

#http://www.sthda.com/english/wiki/ggcorrplot-visualization-of-a-correlation-matrix-using-ggplot2 
p4.1 <- barplot.df.wide %>% select(-Var2, -condition, -timepoint) %>% cor(method="spearman") %>% ggcorrplot(hc.order = TRUE, outline.col = "white", type="lower",colors = c("#006B81", "#F4F4F4", "#E04F39")) + ggtitle("correlation by all condition")
p4.1

p4.2 <- barplot.df.wide %>% filter(condition == "HD") %>% select(-Var2, -condition, -timepoint) %>% cor(method="pearson") %>% ggcorrplot(hc.order = TRUE, outline.col = "white", type="lower",colors = c("#006B81", "#F4F4F4", "#E04F39")) + ggtitle("correlation by HD")
p4.2

p4.3 <- barplot.df.wide %>% filter(condition == "RTH") %>% select(-Var2, -condition, -timepoint) %>% cor(method="pearson") %>% ggcorrplot(hc.order = TRUE, outline.col = "white", type="lower",colors = c("#006B81", "#F4F4F4", "#E04F39")) + ggtitle("correlation by RTH")
p4.3

p4.4 <- barplot.df.wide %>% filter(condition == "PTLD") %>% select(-Var2, -condition, -timepoint) %>% cor(method="pearson") %>% ggcorrplot(hc.order = TRUE, outline.col = "white", type="lower",colors = c("#006B81", "#F4F4F4", "#E04F39")) + ggtitle("correlation by PTLD")
p4.4

p4.5 <- barplot.df.wide %>% filter(condition == "PTLDN") %>% select(-Var2, -condition, -timepoint) %>% cor(method="pearson") %>% ggcorrplot(hc.order = TRUE, outline.col = "white", type="lower",colors = c("#006B81", "#F4F4F4", "#E04F39")) + ggtitle("correlation by PTLDN")
p4.5

p4 <- (p4.2 + p4.3) / (p4.4 + p4.5) + plot_layout(guides = "collect")
p4
#ggsave(filename = "../lyme_disease/Manuscript/Figures/Figure3/Mono_Internal_Correlation.pdf", p4, width = 7, height = 7, dpi=300)


Idents(Myeloid.clean) <- "cluster"

#https://www.frontiersin.org/articles/10.3389/fimmu.2019.02035/full 
#https://www.cell.com/immunity/pdf/S1074-7613(19)30334-6.pdf 
#cd48: https://pubmed.ncbi.nlm.nih.gov/16148114/ 
#https://www-science-org.stanford.idm.oclc.org/doi/10.1126/science.aah4573 

#https://www.pnas.org/doi/pdf/10.1073/pnas.2002476117
# NK-Like Monocytes:
# "GZMA", "GNLY", "KLRB1", "NKG7"
Myeloid.clean.universe <- row.names(Myeloid.clean)
Myeloid.clean.universe[str_detect(Myeloid.clean.universe,"CD2")]

FeaturePlot(Myeloid.clean, "IL6", min.cutoff = 0)

pdot <- DotPlot(Myeloid.clean, features = c("S100A9","S100A10","S100A12", "SELL","HLA.DRA","HLA.DPA1","HLA.DMA","CD163", "IL6", "CD83","CCL3","CCL4", 
                                            "IL1B","CXCL8","IFIT3","NAMPT", "SLC2A3","CLEC4E","CD48","ITGB1","FCGR3A","CX3CR1","CD160", "CD28","CTLA4",
                                            "IL32","LEF1", "IL7R","CCL5","CD1C", "FCER1A","CD74","GZMB","CD4"))
pdot <- pdot + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_color_gradient2(low = "#D3D3D3", high = "#8C1515")
pdot

VariableFeatures(Myeloid.clean)
#ggsave(filename = "../Manuscript/Figures/Figure3-autoimmune monocytes/dot_marker_mono.pdf",pdot, width = 10, height = 4, dpi = 300)

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

#ggsave("./Paper_Figures/KIR/KIR2/KIR2.UMAP.pdf",p.CD158e, width = 7, height = 6, dpi=300)

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
#ggsave(filename = "./Paper_Figures/KIR/KIR3/percent.compare.pdf", p.KIR3, width = 6, height = 5,dpi = 300)

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

###################################Identify cluster markers#####################################
Idents(Myeloid.clean) <- "cluster"
table(Myeloid.clean$cluster)

DefaultAssay(Myeloid.clean) <- "SCT"
Myeloid.clean$cluster <- factor(Myeloid.clean$cluster, levels = c("pDC", "cDC", "X5_Mono", "X6_Mono", "CD16_Mono", "NAMPT_Mono", "HyperInflam_Mono", "HLA_Mono", "S100A9_Mono"))
p1.vln <- VlnPlot(Myeloid.clean,  features = c("CCL4", "CCL3", "IL1B","IL6", "TNF", "IL1RN","CXCL8", "CXCL2", "CD83"), group.by = "cluster", pt.size = 0.0, same.y.lims = T) #& stat_compare_means(ref.group = "HyperInflam_Mono", method = "t.test")
p1.vln <- p1.vln & scale_fill_manual(values = Mono_color)
p1.vln

#ggsave(filename = "/Users/xinchen/Library/CloudStorage/Dropbox/lyme_disease/Manuscript/Figures/Figure3-autoimmune monocytes/VlnPlot_Hypoinflam_cluster_Markers_XC.pdf", p1.vln, width = 8, height = 10,dpi = 300)


###### Modyfy by XC

library(patchwork)

# Create list of plots — one per feature
features <- c("CCL4", "CCL3", "IL1B", "IL6", "TNF", "IL1RN", "CXCL8", "CXCL2", "CD83")

vln_list <- lapply(features, function(gene) {
  VlnPlot(Myeloid.clean, features = gene, group.by = "cluster", pt.size = 0.0) +
    scale_fill_manual(values = Mono_color) +
    theme(legend.position = "none",
  axis.title.y = element_blank()
    )
})
vln_list[[1]] <- vln_list[[1]] + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
vln_list[[2]] <- vln_list[[2]] + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
vln_list[[3]] <- vln_list[[3]] + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
vln_list[[4]] <- vln_list[[4]] + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
vln_list[[5]] <- vln_list[[5]] + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
vln_list[[6]] <- vln_list[[6]] + theme(axis.text.x = element_blank(), axis.title.x = element_blank())

library(patchwork)
p1.vln <- wrap_plots(vln_list, ncol = 3)
#ggsave("/Users/xinchen/Library/CloudStorage/Dropbox/lyme_disease/Manuscript/Figures/Figure3-autoimmune monocytes/VlnPlot_Hypoinflam_cluster_Markers_modified.pdf", p1.vln, width = 10, height = 7, dpi = 300)


#####




#########Feature plot##########
VlnPlot(Myeloid.clean, features = c("IFIT3"), group.by = "cluster", pt.size = 0.1) & scale_fill_manual(values = Mono_color)

FeaturePlot(Myeloid.clean, features = c("CXCL10"), min.cutoff = 0)
FeaturePlot(Myeloid.clean, features = c("IFIT3"), min.cutoff = 0)
################################################################################################
#Differential Abundant Analysis
################################################################################################
table(Myeloid.clean$cluster)

DefaultAssay(Myeloid.clean) <- "integrated"
Idents(Myeloid.clean) <- "condition"

S100A9_Mono <- Myeloid.clean %>% subset(cluster == "S100A9_Mono")
DimPlot(S100A9_Mono)
FeaturePlot(S100A9_Mono, "CD14")

PTLD_RTH_S100A9_Mono <- FindMarkers(S100A9_Mono, ident.1 = "PTLD", ident.2 = "RTH", logfc.threshold = log(1.1),test.use = "MAST", min.pct = 0.3)
PTLDN_HD_S100A9_Mono <- FindMarkers(S100A9_Mono, ident.1 = "PTLDN", ident.2 = "HD", logfc.threshold = log(1.1),test.use = "MAST", min.pct = 0.3)

PTLD_RTH_S100A9_Mono$gene <- rownames(PTLD_RTH_S100A9_Mono)
PTLDN_HD_S100A9_Mono$gene <- rownames(PTLDN_HD_S100A9_Mono)

merged_df.S100A9 <- merge(PTLD_RTH_S100A9_Mono,PTLDN_HD_S100A9_Mono, by="gene") %>% filter(gene != "IL32")

x_max <- max(abs(min(merged_df.S100A9$avg_log2FC.x, na.rm = TRUE)), max(merged_df.S100A9$avg_log2FC.x, na.rm = TRUE))
y_max <- max(abs(min(merged_df.S100A9$avg_log2FC.y, na.rm = TRUE)), max(merged_df.S100A9$avg_log2FC.y, na.rm = TRUE))

merged_df.S100A9_sig <- filter(merged_df.S100A9, p_val_adj.x < 0.05 | p_val_adj.y < 0.05)

filter(merged_df.S100A9_sig, avg_log2FC.x < 0 & avg_log2FC.y < 0)$gene
filter(merged_df.S100A9_sig, avg_log2FC.x > 0 & avg_log2FC.y > 0)$gene

p.S100A9 <- merged_df.S100A9_sig %>%
  ggplot(aes(x=avg_log2FC.x, y=avg_log2FC.y, label=gene)) + 
  geom_point() +
  coord_cartesian(xlim = c(-x_max, x_max), ylim = c(-y_max, y_max)) +
  scale_x_continuous(sec.axis = sec_axis(~., name = "")) +
  scale_y_continuous(sec.axis = sec_axis(~., name = "")) + 
  geom_vline(xintercept = 0, linetype="dashed", color = "black") +
  geom_hline(yintercept = 0, linetype="dashed", color = "black") +
  geom_label_repel() + theme_classic() + xlab("RTH <-- | --> PTLD") + ylab("HD <-- | --> PTLDN") + ggtitle("S100A9 Mono")
p.S100A9

slice3.VP_S100A9 <- EnhancedVolcano(PTLDN_HD_S100A9_Mono,
                                  lab = PTLDN_HD_S100A9_Mono$gene,
                                  FCcutoff=0.15,
                                  xlim = c(-0.7, 0.7),
                                  x = 'avg_log2FC',
                                  y = 'p_val_adj', 
                                  title = "S100A9 Mono")
slice3.VP_S100A9

##################################################################### Subset CD16 Monocytes
CD16_Mono <- Myeloid.clean %>% subset(cluster == "CD16_Mono")

DimPlot(CD16_Mono)
FeaturePlot(CD16_Mono, "FCGR3A")

PTLD_RTH_CD16_Mono <- FindMarkers(CD16_Mono, ident.1 = "PTLD", ident.2 = "RTH", logfc.threshold = log(1.1), test.use = "MAST", min.pct = 0.3)
PTLDN_HD_CD16_Mono <- FindMarkers(CD16_Mono, ident.1 = "PTLDN", ident.2 = "HD", logfc.threshold = log(1.1), test.use = "MAST", min.pct = 0.3)

PTLD_RTH_CD16_Mono$gene <- rownames(PTLD_RTH_CD16_Mono)
PTLDN_HD_CD16_Mono$gene <- rownames(PTLDN_HD_CD16_Mono)

merged_df.CD16 <- merge(PTLD_RTH_CD16_Mono, PTLDN_HD_CD16_Mono, by = "gene") %>% filter(gene != "IL32")

x_max <- max(abs(min(merged_df.CD16$avg_log2FC.x, na.rm = TRUE)), max(merged_df.CD16$avg_log2FC.x, na.rm = TRUE))
y_max <- max(abs(min(merged_df.CD16$avg_log2FC.y, na.rm = TRUE)), max(merged_df.CD16$avg_log2FC.y, na.rm = TRUE))

merged_df.CD16_sig <- filter(merged_df.CD16, p_val_adj.x < 0.05 | p_val_adj.y < 0.05)

filter(merged_df.CD16_sig, avg_log2FC.x < 0 & avg_log2FC.y < 0)$gene
filter(merged_df.CD16_sig, avg_log2FC.x > 0 & avg_log2FC.y > 0)$gene

p.CD16 <- merged_df.CD16_sig %>%
  ggplot(aes(x = avg_log2FC.x, y = avg_log2FC.y, label = gene)) + 
  geom_point() +
  coord_cartesian(xlim = c(-x_max, x_max), ylim = c(-y_max, y_max)) +
  scale_x_continuous(sec.axis = sec_axis(~., name = "")) +
  scale_y_continuous(sec.axis = sec_axis(~., name = "")) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_label_repel() + theme_classic() + xlab("RTH <-- | --> PTLD") + ylab("HD <-- | --> PTLDN") + ggtitle("CD16 Mono")
p.CD16

slice3.VP_CD16 <- EnhancedVolcano(PTLDN_HD_CD16_Mono,
                                  lab = PTLDN_HD_CD16_Mono$gene,
                                  FCcutoff = 0.15,
                                  xlim = c(-0.7, 0.7),
                                  x = 'avg_log2FC',
                                  y = 'p_val_adj', 
                                  title = "CD16 Mono")
slice3.VP_CD16

########################################## Subset HyperInflam Monocytes
HyperInflam_Mono <- Myeloid.clean %>% subset(cluster == "HyperInflam_Mono")

DimPlot(HyperInflam_Mono)
FeaturePlot(HyperInflam_Mono, "CCL3", split.by = "condition")

PTLD_RTH_HyperInflam_Mono <- FindMarkers(HyperInflam_Mono, ident.1 = "PTLD", ident.2 = "RTH", logfc.threshold = log(1.1), test.use = "MAST", min.pct = 0.3)
PTLDN_HD_HyperInflam_Mono <- FindMarkers(HyperInflam_Mono, ident.1 = "PTLDN", ident.2 = "HD", logfc.threshold = log(1.1), test.use = "MAST", min.pct = 0.3)

PTLD_RTH_HyperInflam_Mono$gene <- rownames(PTLD_RTH_HyperInflam_Mono)
PTLDN_HD_HyperInflam_Mono$gene <- rownames(PTLDN_HD_HyperInflam_Mono)

merged_df.HyperInflam <- merge(PTLD_RTH_HyperInflam_Mono, PTLDN_HD_HyperInflam_Mono, by = "gene") %>% filter(gene != "IL32")

x_max <- max(abs(min(merged_df.HyperInflam$avg_log2FC.x, na.rm = TRUE)), max(merged_df.HyperInflam$avg_log2FC.x, na.rm = TRUE))
y_max <- max(abs(min(merged_df.HyperInflam$avg_log2FC.y, na.rm = TRUE)), max(merged_df.HyperInflam$avg_log2FC.y, na.rm = TRUE))

merged_df.HyperInflam_sig <- filter(merged_df.HyperInflam, p_val_adj.x < 0.05 | p_val_adj.y < 0.05)

filter(merged_df.HyperInflam_sig, avg_log2FC.x < 0 & avg_log2FC.y < 0)$gene
filter(merged_df.HyperInflam_sig, avg_log2FC.x > 0 & avg_log2FC.y > 0)$gene

p.HyperInflam <- merged_df.HyperInflam_sig %>%
  ggplot(aes(x = avg_log2FC.x, y = avg_log2FC.y, label = gene)) + 
  geom_point() +
  coord_cartesian(xlim = c(-x_max, x_max), ylim = c(-y_max, y_max)) +
  scale_x_continuous(sec.axis = sec_axis(~., name = "")) +
  scale_y_continuous(sec.axis = sec_axis(~., name = "")) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_label_repel() + theme_classic() + xlab("RTH <-- | --> PTLD") + ylab("HD <-- | --> PTLDN") + ggtitle("HyperInflam Mono")
p.HyperInflam

slice3.VP_HyperInflam <- EnhancedVolcano(PTLDN_HD_HyperInflam_Mono,
                                         lab = PTLDN_HD_HyperInflam_Mono$gene,
                                         FCcutoff = 0.15,
                                         xlim = c(-2.3, 2.3),
                                         x = 'avg_log2FC',
                                         y = 'p_val_adj', 
                                         title = "HyperInflam Mono")

# slice3.VP_HyperInflam <- EnhancedVolcano(PTLD_RTH_HyperInflam_Mono,
#                                          lab = PTLD_RTH_HyperInflam_Mono$gene,
#                                          FCcutoff = 0.15,
#                                          xlim = c(-2.3, 2.3),
#                                          x = 'avg_log2FC',
#                                          y = 'p_val_adj', 
#                                          title = "HyperInflam Mono")
slice3.VP_HyperInflam

########################################## Subset NAMPT Monocytes
NAMPT_Mono <- Myeloid.clean %>% subset(cluster == "NAMPT_Mono")

DimPlot(NAMPT_Mono)
FeaturePlot(NAMPT_Mono, "NAMPT")

PTLD_RTH_NAMPT_Mono <- FindMarkers(NAMPT_Mono, ident.1 = "PTLD", ident.2 = "RTH", logfc.threshold = log(1.1), test.use = "MAST", min.pct = 0.3)
PTLDN_HD_NAMPT_Mono <- FindMarkers(NAMPT_Mono, ident.1 = "PTLDN", ident.2 = "HD", logfc.threshold = log(1.1), test.use = "MAST", min.pct = 0.3)

PTLD_RTH_NAMPT_Mono$gene <- rownames(PTLD_RTH_NAMPT_Mono)
PTLDN_HD_NAMPT_Mono$gene <- rownames(PTLDN_HD_NAMPT_Mono)

merged_df.NAMPT <- merge(PTLD_RTH_NAMPT_Mono, PTLDN_HD_NAMPT_Mono, by = "gene") %>% filter(gene != "IL32")

x_max <- max(abs(min(merged_df.NAMPT$avg_log2FC.x, na.rm = TRUE)), max(merged_df.NAMPT$avg_log2FC.x, na.rm = TRUE))
y_max <- max(abs(min(merged_df.NAMPT$avg_log2FC.y, na.rm = TRUE)), max(merged_df.NAMPT$avg_log2FC.y, na.rm = TRUE))

merged_df.NAMPT_sig <- filter(merged_df.NAMPT, p_val_adj.x < 0.05 | p_val_adj.y < 0.05)

filter(merged_df.NAMPT_sig, avg_log2FC.x < 0 & avg_log2FC.y < 0)$gene
filter(merged_df.NAMPT_sig, avg_log2FC.x > 0 & avg_log2FC.y > 0)$gene

p.NAMPT <- merged_df.NAMPT_sig %>%
  ggplot(aes(x = avg_log2FC.x, y = avg_log2FC.y, label = gene)) + 
  geom_point() +
  coord_cartesian(xlim = c(-x_max, x_max), ylim = c(-y_max, y_max)) +
  scale_x_continuous(sec.axis = sec_axis(~., name = "")) +
  scale_y_continuous(sec.axis = sec_axis(~., name = "")) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_label_repel() + theme_classic() + xlab("RTH <-- | --> PTLD") + ylab("HD <-- | --> PTLDN") + ggtitle("NAMPT Mono")
p.NAMPT

slice3.VP_NAMPT <- EnhancedVolcano(PTLDN_HD_NAMPT_Mono,
                                   lab = PTLDN_HD_NAMPT_Mono$gene,
                                   FCcutoff = 0.15,
                                   xlim = c(-1, 1),
                                   x = 'avg_log2FC',
                                   y = 'p_val_adj', 
                                   title = "NAMPT Mono")
slice3.VP_NAMPT

# Subset HLA Monocytes
HLA_Mono <- Myeloid.clean %>% subset(cluster == "HLA_Mono")

DimPlot(HLA_Mono)
FeaturePlot(HLA_Mono, "CD4")

PTLD_RTH_HLA_Mono <- FindMarkers(HLA_Mono, ident.1 = "PTLD", ident.2 = "RTH", logfc.threshold = log(1.1), test.use = "MAST", min.pct = 0.3)
PTLDN_HD_HLA_Mono <- FindMarkers(HLA_Mono, ident.1 = "PTLDN", ident.2 = "HD", logfc.threshold = log(1.1), test.use = "MAST", min.pct = 0.3)

PTLD_RTH_HLA_Mono$gene <- rownames(PTLD_RTH_HLA_Mono)
PTLDN_HD_HLA_Mono$gene <- rownames(PTLDN_HD_HLA_Mono)

merged_df.HLA <- merge(PTLD_RTH_HLA_Mono, PTLDN_HD_HLA_Mono, by = "gene") %>% filter(gene != "IL32")

x_max <- max(abs(min(merged_df.HLA$avg_log2FC.x, na.rm = TRUE)), max(merged_df.HLA$avg_log2FC.x, na.rm = TRUE))
y_max <- max(abs(min(merged_df.HLA$avg_log2FC.y, na.rm = TRUE)), max(merged_df.HLA$avg_log2FC.y, na.rm = TRUE))

merged_df.HLA_sig <- filter(merged_df.HLA, p_val_adj.x < 0.05 | p_val_adj.y < 0.05)

filter(merged_df.HLA_sig, avg_log2FC.x < 0 & avg_log2FC.y < 0)$gene
filter(merged_df.HLA_sig, avg_log2FC.x > 0 & avg_log2FC.y > 0)$gene

p.HLA <- merged_df.HLA_sig %>%
  ggplot(aes(x = avg_log2FC.x, y = avg_log2FC.y, label = gene)) + 
  geom_point() +
  coord_cartesian(xlim = c(-x_max, x_max), ylim = c(-y_max, y_max)) +
  scale_x_continuous(sec.axis = sec_axis(~., name = "")) +
  scale_y_continuous(sec.axis = sec_axis(~., name = "")) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_label_repel() + theme_classic() + xlab("RTH <-- | --> PTLD") + ylab("HD <-- | --> PTLDN") + ggtitle("HLA Mono")
p.HLA

slice3.VP_HLA <- EnhancedVolcano(PTLDN_HD_HLA_Mono,
                                 lab = PTLDN_HD_HLA_Mono$gene,
                                 FCcutoff = 0.15,
                                 xlim = c(-0.75, 0.75),
                                 x = 'avg_log2FC',
                                 y = 'p_val_adj', 
                                 title = "HLA Mono")
slice3.VP_HLA

############# Subset X5 Monocytes############
X5_Mono <- Myeloid.clean %>% subset(cluster == "X5_Mono")

DimPlot(X5_Mono)
FeaturePlot(X5_Mono, "CCL3")

PTLD_RTH_X5_Mono <- FindMarkers(X5_Mono, ident.1 = "PTLD", ident.2 = "RTH", logfc.threshold = log(1.1), test.use = "MAST", min.pct = 0.3)
PTLDN_HD_X5_Mono <- FindMarkers(X5_Mono, ident.1 = "PTLDN", ident.2 = "HD", logfc.threshold = log(1.1), test.use = "MAST", min.pct = 0.3)

PTLD_RTH_X5_Mono$gene <- rownames(PTLD_RTH_X5_Mono)
PTLDN_HD_X5_Mono$gene <- rownames(PTLDN_HD_X5_Mono)

merged_df.X5 <- merge(PTLD_RTH_X5_Mono, PTLDN_HD_X5_Mono, by = "gene") %>% filter(gene != "IL32")

x_max <- max(abs(min(merged_df.X5$avg_log2FC.x, na.rm = TRUE)), max(merged_df.X5$avg_log2FC.x, na.rm = TRUE))
y_max <- max(abs(min(merged_df.X5$avg_log2FC.y, na.rm = TRUE)), max(merged_df.X5$avg_log2FC.y, na.rm = TRUE))

merged_df.X5_sig <- filter(merged_df.X5, p_val_adj.x < 0.05 | p_val_adj.y < 0.05)

filter(merged_df.X5_sig, avg_log2FC.x < 0 & avg_log2FC.y < 0)$gene
filter(merged_df.X5_sig, avg_log2FC.x > 0 & avg_log2FC.y > 0)$gene

p.X5 <- merged_df.X5_sig %>%
  ggplot(aes(x = avg_log2FC.x, y = avg_log2FC.y, label = gene)) + 
  geom_point() +
  coord_cartesian(xlim = c(-x_max, x_max), ylim = c(-y_max, y_max)) +
  scale_x_continuous(sec.axis = sec_axis(~., name = "")) +
  scale_y_continuous(sec.axis = sec_axis(~., name = "")) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_label_repel() + theme_classic() + xlab("RTH <-- | --> PTLD") + ylab("HD <-- | --> PTLDN") + ggtitle("X5 Mono")
p.X5

slice3.VP_X5 <- EnhancedVolcano(PTLDN_HD_X5_Mono,
                                lab = PTLDN_HD_X5_Mono$gene,
                                FCcutoff = 0.15,
                                xlim = c(-1, 1),
                                x = 'avg_log2FC',
                                y = 'p_val_adj', 
                                title = "X5 Mono")
slice3.VP_X5

################# Subset X6 Monocytes
X6_Mono <- Myeloid.clean %>% subset(cluster == "X6_Mono")

DimPlot(X6_Mono)
FeaturePlot(X6_Mono, "F13A1")

PTLD_RTH_X6_Mono <- FindMarkers(X6_Mono, ident.1 = "PTLD", ident.2 = "RTH", logfc.threshold = log(1.1), test.use = "MAST", min.pct = 0.3)
PTLDN_HD_X6_Mono <- FindMarkers(X6_Mono, ident.1 = "PTLDN", ident.2 = "HD", logfc.threshold = log(1.1), test.use = "MAST", min.pct = 0.3)

PTLD_RTH_X6_Mono$gene <- rownames(PTLD_RTH_X6_Mono)
PTLDN_HD_X6_Mono$gene <- rownames(PTLDN_HD_X6_Mono)

merged_df.X6 <- merge(PTLD_RTH_X6_Mono, PTLDN_HD_X6_Mono, by = "gene") %>% filter(gene != "IL32")

x_max <- max(abs(min(merged_df.X6$avg_log2FC.x, na.rm = TRUE)), max(merged_df.X6$avg_log2FC.x, na.rm = TRUE))
y_max <- max(abs(min(merged_df.X6$avg_log2FC.y, na.rm = TRUE)), max(merged_df.X6$avg_log2FC.y, na.rm = TRUE))

merged_df.X6_sig <- filter(merged_df.X6, p_val_adj.x < 0.05 | p_val_adj.y < 0.05)

filter(merged_df.X6_sig, avg_log2FC.x < 0 & avg_log2FC.y < 0)$gene
filter(merged_df.X6_sig, avg_log2FC.x > 0 & avg_log2FC.y > 0)$gene

p.X6 <- merged_df.X6_sig %>%
  ggplot(aes(x = avg_log2FC.x, y = avg_log2FC.y, label = gene)) + 
  geom_point() +
  coord_cartesian(xlim = c(-x_max, x_max), ylim = c(-y_max, y_max)) +
  scale_x_continuous(sec.axis = sec_axis(~., name = "")) +
  scale_y_continuous(sec.axis = sec_axis(~., name = "")) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_label_repel() + theme_classic() + xlab("RTH <-- | --> PTLD") + ylab("HD <-- | --> PTLDN") + ggtitle("X6 Mono")
p.X6

slice3.VP_X6 <- EnhancedVolcano(PTLDN_HD_X6_Mono,
                                lab = PTLDN_HD_X6_Mono$gene,
                                FCcutoff = 0.15,
                                xlim = c(-1, 1),
                                x = 'avg_log2FC',
                                y = 'p_val_adj', 
                                title = "X6 Mono")
slice3.VP_X6

######################### Subset cDC (conventional Dendritic Cells)
cDC <- Myeloid.clean %>% subset(cluster == "cDC")

DimPlot(cDC)
FeaturePlot(cDC, "CD1C")

PTLD_RTH_cDC <- FindMarkers(cDC, ident.1 = "PTLD", ident.2 = "RTH", logfc.threshold = log(1.1), test.use = "MAST", min.pct = 0.3)
PTLDN_HD_cDC <- FindMarkers(cDC, ident.1 = "PTLDN", ident.2 = "HD", logfc.threshold = log(1.1), test.use = "MAST", min.pct = 0.3)

PTLD_RTH_cDC$gene <- rownames(PTLD_RTH_cDC)
PTLDN_HD_cDC$gene <- rownames(PTLDN_HD_cDC)

merged_df.cDC <- merge(PTLD_RTH_cDC, PTLDN_HD_cDC, by = "gene") %>% filter(gene != "IL32")

x_max <- max(abs(min(merged_df.cDC$avg_log2FC.x, na.rm = TRUE)), max(merged_df.cDC$avg_log2FC.x, na.rm = TRUE))
y_max <- max(abs(min(merged_df.cDC$avg_log2FC.y, na.rm = TRUE)), max(merged_df.cDC$avg_log2FC.y, na.rm = TRUE))

merged_df.cDC_sig <- filter(merged_df.cDC, p_val_adj.x < 0.05 | p_val_adj.y < 0.05)

filter(merged_df.cDC_sig, avg_log2FC.x < 0 & avg_log2FC.y < 0)$gene
filter(merged_df.cDC_sig, avg_log2FC.x > 0 & avg_log2FC.y > 0)$gene

p.cDC <- merged_df.cDC_sig %>%
  ggplot(aes(x = avg_log2FC.x, y = avg_log2FC.y, label = gene)) + 
  geom_point() +
  coord_cartesian(xlim = c(-x_max, x_max), ylim = c(-y_max, y_max)) +
  scale_x_continuous(sec.axis = sec_axis(~., name = "")) +
  scale_y_continuous(sec.axis = sec_axis(~., name = "")) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_label_repel() + theme_classic() + xlab("RTH <-- | --> PTLD") + ylab("HD <-- | --> PTLDN") + ggtitle("cDC")
p.cDC

slice3.VP_cDC <- EnhancedVolcano(PTLDN_HD_cDC,
                                 lab = PTLDN_HD_cDC$gene,
                                 FCcutoff = 0.15,
                                 xlim = c(-1, 1),
                                 x = 'avg_log2FC',
                                 y = 'p_val_adj', 
                                 title = "cDC")
slice3.VP_cDC

# Subset pDC (plasmacytoid Dendritic Cells)
pDC <- Myeloid.clean %>% subset(cluster == "pDC")

# Visualize the subsetLEF1
DimPlot(pDC)
FeaturePlot(Myeloid.clean, "IL32")

PTLD_RTH_pDC <- FindMarkers(pDC, ident.1 = "PTLD", ident.2 = "RTH", logfc.threshold = log(1.1), test.use = "MAST", min.pct = 0.3)
PTLDN_HD_pDC <- FindMarkers(pDC, ident.1 = "PTLDN", ident.2 = "HD", logfc.threshold = log(1.1), test.use = "MAST", min.pct = 0.3)

PTLD_RTH_pDC$gene <- rownames(PTLD_RTH_pDC)
PTLDN_HD_pDC$gene <- rownames(PTLDN_HD_pDC)

merged_df.pDC <- merge(PTLD_RTH_pDC, PTLDN_HD_pDC, by = "gene") %>% filter(gene != "IL32")

x_max <- max(abs(min(merged_df.pDC$avg_log2FC.x, na.rm = TRUE)), max(merged_df.pDC$avg_log2FC.x, na.rm = TRUE))
y_max <- max(abs(min(merged_df.pDC$avg_log2FC.y, na.rm = TRUE)), max(merged_df.pDC$avg_log2FC.y, na.rm = TRUE))

merged_df.pDC_sig <- filter(merged_df.pDC, p_val_adj.x < 0.05 | p_val_adj.y < 0.05)

filter(merged_df.pDC_sig, avg_log2FC.x < 0 & avg_log2FC.y < 0)$gene
filter(merged_df.pDC_sig, avg_log2FC.x > 0 & avg_log2FC.y > 0)$gene

p.pDC <- merged_df.pDC_sig %>%
  ggplot(aes(x = avg_log2FC.x, y = avg_log2FC.y, label = gene)) + 
  geom_point() +
  coord_cartesian(xlim = c(-x_max, x_max), ylim = c(-y_max, y_max)) +
  scale_x_continuous(sec.axis = sec_axis(~., name = "")) +
  scale_y_continuous(sec.axis = sec_axis(~., name = "")) + 
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_label_repel() + theme_classic() + xlab("RTH <-- | --> PTLD") + ylab("HD <-- | --> PTLDN") + ggtitle("pDC")
p.pDC

slice3.VP_pDC <- EnhancedVolcano(PTLDN_HD_pDC,
                                 lab = PTLDN_HD_pDC$gene,
                                 FCcutoff = 0.15,
                                 xlim = c(-1, 1),
                                 x = 'avg_log2FC',
                                 y = 'p_val_adj', 
                                 title = "pDC")
slice3.VP_pDC

ggsave(filename = "../lyme_disease/Manuscript/Figures/FigureS3/S100A9_double.pdf", p.S100A9, width = 10, height = 10, dpi = 300)
ggsave(filename = "../lyme_disease/Manuscript/Figures/FigureS3/CD16_double.pdf", p.CD16, width = 10, height = 10, dpi = 300)
ggsave(filename = "../lyme_disease/Manuscript/Figures/FigureS3/cDC_double.pdf", p.cDC, width = 10, height = 10, dpi = 300)
ggsave(filename = "../lyme_disease/Manuscript/Figures/FigureS3/HLA_double.pdf", p.HLA, width = 10, height = 10, dpi = 300)
ggsave(filename = "../lyme_disease/Manuscript/Figures/FigureS3/HyperInflam_double.pdf", p.HyperInflam, width = 10, height = 10, dpi = 300)
ggsave(filename = "../lyme_disease/Manuscript/Figures/FigureS3/NAMPT_double.pdf", p.NAMPT, width = 10, height = 10, dpi = 300)
ggsave(filename = "../lyme_disease/Manuscript/Figures/FigureS3/pDC_double.pdf", p.pDC, width = 10, height = 10, dpi = 300)
ggsave(filename = "../lyme_disease/Manuscript/Figures/FigureS3/X5_double.pdf", p.X5, width = 10, height = 10, dpi = 300)
ggsave(filename = "../lyme_disease/Manuscript/Figures/FigureS3/X6_double.pdf", p.X6, width = 10, height = 10, dpi = 300)

ggsave(filename = "../lyme_disease/Manuscript/Figures/FigureS3/S100A9_PTLDN_HD.pdf", slice3.VP_S100A9, width = 10, height = 10, dpi = 300)
ggsave(filename = "../lyme_disease/Manuscript/Figures/FigureS3/CD16_PTLDN_HD.pdf", slice3.VP_CD16, width = 10, height = 10, dpi = 300)
ggsave(filename = "../lyme_disease/Manuscript/Figures/FigureS3/cDC_PTLDN_HD.pdf", slice3.VP_cDC, width = 10, height = 10, dpi = 300)
ggsave(filename = "../lyme_disease/Manuscript/Figures/FigureS3/HLA_PTLDN_HD.pdf", slice3.VP_HLA, width = 10, height = 10, dpi = 300)
ggsave(filename = "../lyme_disease/Manuscript/Figures/FigureS3/HyperInflam_PTLDN_HD.pdf", slice3.VP_HyperInflam, width = 10, height = 10, dpi = 300)
ggsave(filename = "../lyme_disease/Manuscript/Figures/FigureS3/NAMPT_PTLDN_HD.pdf", slice3.VP_NAMPT, width = 10, height = 10, dpi = 300)
ggsave(filename = "../lyme_disease/Manuscript/Figures/FigureS3/pDC_PTLDN_HD.pdf", slice3.VP_pDC, width = 10, height = 10, dpi = 300)
ggsave(filename = "../lyme_disease/Manuscript/Figures/FigureS3/X5_PTLDN_HD.pdf", slice3.VP_X5, width = 10, height = 10, dpi = 300)
ggsave(filename = "../lyme_disease/Manuscript/Figures/FigureS3/X6_PTLDN_HD.pdf", slice3.VP_X6, width = 10, height = 10, dpi = 300)

################################merge
PTLDN_HD_S100A9_Mono$celltype <- "S100A9_Mono"
PTLDN_HD_CD16_Mono$celltype <- "CD16_Mono"
PTLDN_HD_HyperInflam_Mono$celltype <- "HyperInflam_Mono"
PTLDN_HD_NAMPT_Mono$celltype <- "NAMPT_Mono"
PTLDN_HD_HLA_Mono$celltype <- "HLA_Mono"
PTLDN_HD_X5_Mono$celltype <- "X5_Mono"
PTLDN_HD_X6_Mono$celltype <- "X6_Mono"
PTLDN_HD_cDC$celltype <- "cDC"
PTLDN_HD_pDC$celltype <- "pDC"

PTLDN_HD_MONO_ALL <- rbind(PTLDN_HD_S100A9_Mono,PTLDN_HD_CD16_Mono,
                           PTLDN_HD_HyperInflam_Mono,PTLDN_HD_NAMPT_Mono,
                           PTLDN_HD_HLA_Mono,PTLDN_HD_X5_Mono,PTLDN_HD_X6_Mono,
                           PTLDN_HD_cDC,PTLDN_HD_pDC)

#write.csv(file = "../lyme_disease/Manuscript/Figures/Figure3/All_lympho_PTLDN_HD_DEG.csv", PTLDN_HD_MONO_ALL)
#PTLDN_HD_MONO_ALL <- read.csv("../lyme_disease/Manuscript/Figures/Figure3/All_lympho_PTLDN_HD_DEG.csv")

gene.for.anno <- c("CCL3","CCL4","EGR1","CX3CR1", "TNF", "IL1B", "CCL3L1", "BCL2A1",
                   "JUNB", "HLA.DPA1","TYROBP", "FCN1")

CD16_Mono_genes <- filter(gene.for.anno, celltype == "CD16_Mono")$gene
cDC_genes <- filter(gene.for.anno, celltype == "cDC")$gene
HLA_Mono_genes <- filter(gene.for.anno, celltype == "HLA_Mono")$gene
HyperInflam_Mono_genes <- filter(gene.for.anno, celltype == "HyperInflam_Mono")$gene
NAMPT_Mono_genes <- filter(gene.for.anno, celltype == "NAMPT_Mono")$gene
pDC_genes <- filter(gene.for.anno, celltype == "pDC")$gene
S100A9_Mono_genes <- filter(gene.for.anno, celltype == "S100A9_Mono")$gene
X5_Mono_genes <- filter(gene.for.anno, celltype == "X5_Mono")$gene
X6_Mono_genes <- filter(gene.for.anno, celltype == "X6_Mono")$gene

p_monocytes <- filter(PTLDN_HD_MONO_ALL, p_val_adj < 0.05) %>%
  ggplot(aes(x = celltype, y = avg_log2FC, label = gene)) +
  geom_jitter(aes(color = 
                    (gene %in% CD16_Mono_genes & celltype == "CD16_Mono") |
                    (gene %in% cDC_genes & celltype == "cDC") |
                    (gene %in% HLA_Mono_genes & celltype == "HLA_Mono") |
                    (gene %in% HyperInflam_Mono_genes & celltype == "HyperInflam_Mono") |
                    (gene %in% NAMPT_Mono_genes & celltype == "NAMPT_Mono") |
                    (gene %in% pDC_genes & celltype == "pDC") |
                    (gene %in% S100A9_Mono_genes & celltype == "S100A9_Mono") |
                    (gene %in% X5_Mono_genes & celltype == "X5_Mono") |
                    (gene %in% X6_Mono_genes & celltype == "X6_Mono"))) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey")) +
  
  geom_text_repel(data = subset(PTLDN_HD_MONO_ALL, (gene %in% CD16_Mono_genes & celltype == "CD16_Mono")), vjust = -1, hjust = 0.5, color = "red") +
  geom_text_repel(data = subset(PTLDN_HD_MONO_ALL, (gene %in% cDC_genes & celltype == "cDC")), vjust = -1, hjust = 0.5, color = "red") +
  geom_text_repel(data = subset(PTLDN_HD_MONO_ALL, (gene %in% HLA_Mono_genes & celltype == "HLA_Mono")), vjust = -1, hjust = 0.5, color = "red") +
  geom_text_repel(data = subset(PTLDN_HD_MONO_ALL, (gene %in% HyperInflam_Mono_genes & celltype == "HyperInflam_Mono")), vjust = -1, hjust = 0.5, color = "red") +
  geom_text_repel(data = subset(PTLDN_HD_MONO_ALL, (gene %in% NAMPT_Mono_genes & celltype == "NAMPT_Mono")), vjust = -1, hjust = 0.5, color = "red") +
  geom_text_repel(data = subset(PTLDN_HD_MONO_ALL, (gene %in% pDC_genes & celltype == "pDC")), vjust = -1, hjust = 0.5, color = "red") +
  geom_text_repel(data = subset(PTLDN_HD_MONO_ALL, (gene %in% S100A9_Mono_genes & celltype == "S100A9_Mono")), vjust = -1, hjust = 0.5, color = "red") +
  geom_text_repel(data = subset(PTLDN_HD_MONO_ALL, (gene %in% X5_Mono_genes & celltype == "X5_Mono")), vjust = -1, hjust = 0.5, color = "red") +
  geom_text_repel(data = subset(PTLDN_HD_MONO_ALL, (gene %in% X6_Mono_genes & celltype == "X6_Mono")), vjust = -1, hjust = 0.5, color = "red") +
  
  theme_minimal() +
  labs(title = "Differentially Expressed Genes By Cell Type (Monocytes/DC)",
       x = "Cell Type",
       y = "Log2 Fold Change") +
  XZ_flip_x_label() +
  guides(color = FALSE)

p_monocytes


