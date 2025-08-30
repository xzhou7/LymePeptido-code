#DEGs and Pathway
devtools::install_github("sajuukLyu/ggunchull", type = "source")
devtools::install_github('chris-mcginnis-ucsf/DoubletFinder')
devtools::install_github("immunogenomics/harmony")
# Install Harmony from CRAN
install.packages("harmony")

library(dplyr)
library(tidyr)
library(Seurat)
library(patchwork)
library(data.table)
library(ggplot2)
library(stringr)
library(metap)
library(multtest)
library(Rcpp)
library(sctransform)
library(glmGamPoi)
library(cowplot)
library(ggrepel)
library(tidydr)
library(ggsci)
library(Cairo)
library(DoubletFinder)
library(harmony)
library(ggpubr)
library(scales)
#library(ggunchull)
library(EnhancedVolcano)
#install.packages("ggrepel")
library(ggrepel)
#install.packages("devtools")


#mac
setwd("~/Library/CloudStorage/Dropbox/lyme_disease/Pepti_Spleen_singlecell/")
getwd()

#load data
load(file = "./RObject/02.Celltype_labelled_UMAP.RData")
source("Code/00_colorKey_XC.R")

# Set default assay if needed
DefaultAssay(spleen.combined.sct) <- "integrated"  # or "SCT" if using integrated data

# Define your genes of interest
genes.of.interest <- c("CCL3", "CCL4", "CCL8", "CXCL2")  # example genes

# Set condition identity
Idents(spleen.combined.sct) <- "seurat_clusters"
DefaultAssay(spleen.combined.sct) <- "integrated"
# Plot
VlnPlot(spleen.combined.sct, features = genes.of.interest, group.by = "orig.ident", pt.size = 0)


table(spleen.combined.sct$seurat_clusters)
DimPlot(spleen.combined.sct, label = TRUE)
Idents(spleen.combined.sct)

#load
DimPlot(spleen.combined.sct) + coord_fixed()

Idents(spleen.combined.sct) <- spleen.combined.sct$seurat_clusters
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '0' = "Mem_B")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '1' = "Myeloid_cells")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '2' = "Naive_B")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '3' = "Mem_B")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '4' = "Mem_B")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '5' = "GZMK_CD8T")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '6' = "Mem_B")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '7' = "Naive_CD4T")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '8' = "GerminalCenter_B")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '9' = "Mem_B")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '10' = "Rest_NK")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '11' = "Mem_B")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '12' = "GZMB_CD8T")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '13' = "Mem_CD4")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '14' = "GerminalCenter_B")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '15' = "Mem_CD8T")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '16' = "Mem_CD8T")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '17' = "Activated_NK")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '18' = "Plasmablast")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '19' = "Undetermined")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '20' = "Gamma_Delta_T")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '21' = "Treg")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '22' = "Plasmablast")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '23' = "Naive_B")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '24' = "Naive_CD8T")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '25' = "Naive_B")
spleen.combined.sct <- RenameIdents(object = spleen.combined.sct, '26' = "Myeloid_cells")
spleen.combined.sct$celltype2 <- Idents(spleen.combined.sct)

DimPlot(spleen.combined.sct) + coord_fixed() + scale_color_manual(values = pep_cluster2_color)

FeaturePlot(spleen.combined.sct, "rna_IL17F",split.by = "orig.ident", order = T, min.cutoff = 0, max.cutoff = 1) & coord_fixed()
FeaturePlot(spleen.combined.sct, "CCL22",split.by = "orig.ident", order = T, min.cutoff = 0, max.cutoff = 1) & coord_fixed()
FeaturePlot(spleen.combined.sct, "IL1B",split.by = "orig.ident", order = T, min.cutoff = 0, max.cutoff = 1) & coord_fixed()


FeaturePlot(spleen.combined.sct, "CXCL1",split.by = "orig.ident", order = T, min.cutoff = 0, max.cutoff = 1) & coord_fixed()
FeaturePlot(spleen.combined.sct, "CXCL2",split.by = "orig.ident", order = T, min.cutoff = 0, max.cutoff = 1) & coord_fixed()
FeaturePlot(spleen.combined.sct, "CXCL3",split.by = "orig.ident", order = T, min.cutoff = 0, max.cutoff = 1) & coord_fixed()
FeaturePlot(spleen.combined.sct, "CXCL5",split.by = "orig.ident", order = T, min.cutoff = 0, max.cutoff = 1) & coord_fixed()
FeaturePlot(spleen.combined.sct, "CXCL8",split.by = "orig.ident", order = T, min.cutoff = 0, max.cutoff = 1) & coord_fixed()
FeaturePlot(spleen.combined.sct, "CXCL10",split.by = "orig.ident", order = T, min.cutoff = 0, max.cutoff = 1) & coord_fixed()
FeaturePlot(spleen.combined.sct, "CCL3",split.by = "orig.ident", order = T, min.cutoff = 0, max.cutoff = 1) & coord_fixed()
FeaturePlot(spleen.combined.sct, "CCL7",split.by = "orig.ident", order = T, min.cutoff = 0, max.cutoff = 1) & coord_fixed()

FeaturePlot(spleen.combined.sct, "CDH5", order = T, min.cutoff = 0, max.cutoff = 10)

DefaultAssay(spleen.combined.sct) <- "integrated"
stromal_genes <- c("VIM", "SPARC", "CD44", "COL6A1", "MMP9", "TIMP1")
spleen.combined.sct <- AddModuleScore(spleen.combined.sct, features = list(stromal_genes), name = "StromalScore")
stromalscore <- FeaturePlot(spleen.combined.sct, features = "StromalScore1",min.cutoff = 0, max.cutoff = 4, order = T) + coord_equal()
stromalscore
#ggsave(filename = "../Manuscript/Figures/Figure 5-immune organoids and autoreactive B cells/Stromal_UMAP.pdf",stromalscore,width = 5, height = 6, dpi = 300)

stromal_vln <- VlnPlot(
  spleen.combined.sct,
  features = "StromalScore1",
  group.by = "celltype2",
  pt.size = 0  # Remove all dots
) +
  NoLegend() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = pep_cluster2_color) +
  geom_hline(yintercept = c(1, 2, 3, 4), linetype = "dashed", color = "black")

stromal_vln

#ggsave(filename = "../Manuscript/Figures/Figure 5-immune organoids and autoreactive B cells/Stromal_vln.pdf",stromal_vln,width = 10, height = 6, dpi = 300)

celltype_position <- spleen.combined.sct@reductions$umap@cell.embeddings %>% as.data.frame() %>%
  cbind(celltype = spleen.combined.sct@meta.data$celltype2) %>%   # or celltype
  group_by(celltype) %>% dplyr::summarise(UMAP_1 = median(umap_1),UMAP_2 = median(umap_2))

p_umap <- DimPlot(spleen.combined.sct, group.by = "celltype2", label = FALSE, cols = pep_cluster2_color) +
  coord_fixed() + theme(legend.position = "none") + ggtitle("Spleen Organoid") + scale_color_manual(values = pep_cluster2_color)
p_umap

p_umap2 <- p_umap + geom_label_repel(
    data = celltype_position, aes(x = UMAP_1, y = UMAP_2, label = celltype, color = celltype),
    fontface = "bold",#size = 6, 
    point.padding = unit(0.1, "lines"),alpha = 0.95,show.legend = FALSE) + scale_color_manual(values = pep_cluster2_color)
p_umap2

ggsave(filename = "../Manuscript/Figures/Figure 5-immune organoids and autoreactive B cells/UMAP.pdf",p_umap2, width = 6.5, height = 6.5, dpi=300 )
#ggsave(filename = "../Manuscript/Figures/Figure 5-immune organoids and autoreactive B cells/UMAP.png",p_umap2, width = 6.5, height = 6.5, dpi=300 )


Idents(spleen.combined.sct) <- spleen.combined.sct$celltype2
# This will find DEGs for each cell type compared to all other cells
#all_markers <- FindAllMarkers(
#  spleen.combined.sct,
#  assay = "RNA",         # Use "RNA" for normalized counts, "SCT" for SCTransform
#  only.pos = TRUE,       # Keep only upregulated genes
#  min.pct = 0.25,        # Expressed in at least 25% of cells in either group
#  logfc.threshold = 0.25 # Log2 fold change threshold
#)

table(spleen.combined.sct$celltype2)
table(spleen.combined.sct$celltype)


Idents(spleen.combined.sct) <- spleen.combined.sct$celltype
DimPlot(spleen.combined.sct, label = T)


colnames(spleen.combined.sct[[]]) 

spleen.combined.sct$celltype <- as.character(spleen.combined.sct$celltype)
spleen.combined.sct$celltype[spleen.combined.sct$celltype == "X13_Th17/MAIT"] <- "X13_Th17_MAIT"
spleen.combined.sct$celltype <- factor(spleen.combined.sct$celltype)

DefaultAssay(spleen.combined.sct) <- "RNA"
spleen.combined.sct <- NormalizeData(spleen.combined.sct,normalization.method = "LogNormalize",scale.factor = 10000)
DefaultAssay(spleen.combined.sct) <- "integrated"

#######################################################################
#Calculate Sample By Cluster-percentagewide celltype no difference, so use celltype2
meta_data <- spleen.combined.sct@meta.data 

cell_counts <- meta_data %>%
  group_by(sample, celltype2) %>%
  summarise(cell_number = n()) %>%
  ungroup()

total_counts <- meta_data %>%
  group_by(sample) %>%
  summarise(total_cells = n())

cell_percentages <- left_join(cell_counts, total_counts, by = "sample") %>%
  mutate(percent = 100 * cell_number / total_cells)

cell_percentages <- left_join(cell_percentages, 
                              meta_data %>% distinct(sample, orig.ident, Hashtag), 
                              by = "sample")
head(cell_percentages)
cell_percentages <- cell_percentages %>%
  mutate(orig.ident = factor(orig.ident, 
                             levels = c("spleen-control", "spleen_peptidoglycan"), 
                             labels = c("Control", "Peptidoglycan")))

table(cell_percentages$celltype2)
celltype2.list <- c("Myeloid_cells", "Naive_B", "Naive_CD8T", "Plasmablast", "Mem_CD4", "Gamma_Delta_T",
                    "Rest_NK", "Mem_CD8T", "GerminalCenter_B", "GZMB_CD8T", "Mem_B", "Undetermined",
                    "Activated_NK", "Naive_CD4T", "GZMK_CD8T","Treg")

for (i in 1:length(celltype2.list)){
  target_cell <- celltype2.list[i]
  print(target_cell)
  p <- ggpaired(cell_percentages %>% filter(celltype2 == target_cell),
                x = "orig.ident",  # you might need to create a column for treatment group
                y = "percent",
                id = "Hashtag",
                line.color = "gray",
                line.size = NA,
                color = "orig.ident",
                palette = "d3") + ggtitle(target_cell) +
    ylab("Percentage") + xlab("") + scale_y_continuous(labels = label_percent(scale = 1))
  p <- p + stat_compare_means(paired = TRUE, method = "t.test")
  p
  path_cellid <- paste0("../Manuscript/Figures/Figure 5-immune organoids and autoreactive B cells/Peptidoglycan_Cell_Percent-XC/Celltype2.",
                        target_cell, ".pdf", sep="")
  
  ggsave(filename = path_cellid, p, width = 3, height = 4, dpi=300)}


ident()
DimPlot(spleen.combined.sct)

singleTcalculate <- cell_percentages %>% filter(celltype2 == "Plasmablast")
singleTcalculate_wide <- singleTcalculate %>%select(Hashtag, orig.ident, percent) %>% pivot_wider(names_from = orig.ident, values_from = percent)
t_test_result <- t.test(singleTcalculate_wide$Peptidoglycan, singleTcalculate_wide$Control, paired = TRUE)
t_test_result

###Dot Plot
spleen.combined.sct$celltype2 <- factor(spleen.combined.sct$celltype2, levels = c("Naive_B", "Mem_B", "GerminalCenter_B", "Plasmablast",
                                                                                  "Naive_CD4T","Mem_CD4","Treg", 
                                                                                  "Naive_CD8T","Mem_CD8T","GZMK_CD8T","GZMB_CD8T","Gamma_Delta_T",
                                                                                  "Rest_NK","Activated_NK",
                                                                                  "Myeloid_cells", "Undetermined"))
Idents(spleen.combined.sct) <- "celltype2"
pdot <- DotPlot(spleen.combined.sct, features = c("MS4A1","CD79A","IGHD","CD38","TNFRSF13B","BCL6", "CD27","AICDA","MKI67", "CCNB1","MZB1","XBP1","IRF8",
                                                   "SELL","CCR7","IL7R",  "CD3D",
                                                  "CD4","GATA3", "FOXP3","RORC","IL2RA","CTLA4","KLRB1","TRDC",
                                                  "CD8A","TIGIT","GZMA", "GZMK","GZMB","FCGR3A","NCAM1", "CD69","TNF",
                                                  "CD14","FCER1G",  "ITGAM", "ITGAL", "ANXA5",  "ITGB1","KIR2DL1","CD160","KLRC2","FSCN1","VIM"))




pdot <- pdot + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_color_gradient(low = "#D3D3D3", high = "#8C1515")
pdot
getwd()
#ggsave(filename = "../Manuscript/Figures/Figure 5-immune organoids and autoreactive B cells/DotPlot.pdf",pdot,width = 12, height = 5, dpi=300 )

FeaturePlot(spleen.combined.sct,"NCAM1",min.cutoff = 0)

#Subset All SubCluster
Spleen_B <- subset(spleen.combined.sct, idents = c("Naive_B", "Mem_B", "GerminalCenter_B", "Plasmablast"))
DimPlot(Spleen_B) + coord_fixed()
umap_dataB <- Embeddings(Spleen_B, reduction = "umap")
cells_to_keepB <- rownames(umap_dataB)[(umap_dataB[,2] >= -5) & (umap_dataB[,1] <= 7.5)]
Spleen_B <- subset(Spleen_B, cells = cells_to_keepB)
DimPlot(Spleen_B) + coord_fixed()

Spleen_CD4T <- subset(spleen.combined.sct, idents = c("Naive_CD4T","Mem_CD4", "Treg"))
DimPlot(Spleen_CD4T) + coord_fixed()
umap_data4T <- Embeddings(Spleen_CD4T, reduction = "umap")
cells_to_keep4T <- rownames(umap_data4T)[(umap_data4T[,2] <= -7) & (umap_data4T[,1] <= 4)]
Spleen_CD4T <- subset(Spleen_CD4T, cells = cells_to_keep4T)
DimPlot(Spleen_CD4T) + coord_fixed()


##########check Treg
levels(Idents(spleen.combined.sct))

Spleen_Treg <- subset(spleen.combined.sct, idents = "X21_Treg")
DimPlot(Spleen_Treg) + coord_fixed()
umap_data4Treg <- Embeddings(Spleen_Treg, reduction = "umap")
cells_to_keep4Teg <- rownames(umap_data4Treg)[(umap_data4Treg[,2] <= -7) & (umap_data4Treg[,1] <= 4)]
Spleen_Treg <- subset(Spleen_Treg, cells = cells_to_keep4Teg)
DimPlot(Spleen_Treg) + coord_fixed()

##########check Th17
Spleen_TH17 <- subset(spleen.combined.sct, idents = "X13_Th17/MAIT")
DimPlot(Spleen_TH17) + coord_fixed()
umap_data4TH17 <- Embeddings(Spleen_TH17, reduction = "umap")
cells_to_keep4TH17 <- rownames(umap_data4TH17)[(umap_data4TH17[,2] <= -7) & (umap_data4TH17[,1] <= 4)]
Spleen_TH17 <- subset(Spleen_TH17, cells = cells_to_keep4TH17)
DimPlot(Spleen_TH17) + coord_fixed()


#####################
Spleen_CD8T <- subset(spleen.combined.sct, idents = c("GZMK_CD8T","GZMB_CD8T","Gamma_Delta_T","Naive_CD8T", "Mem_CD8T"))
DimPlot(Spleen_CD8T) + coord_fixed()
umap_data8T <- Embeddings(Spleen_CD8T, reduction = "umap")
cells_to_keep8T <- rownames(umap_data8T)[(umap_data8T[,2] <= -1.5) & (umap_data8T[,1] <= 7)]
Spleen_CD8T <- subset(Spleen_CD8T, cells = cells_to_keep8T)
DimPlot(Spleen_CD8T) + coord_fixed()

#levels(Idents(spleen.combined.sct))

Spleen_macrophage <- subset(spleen.combined.sct, idents = c("Myeloid_cells"))
DimPlot(Spleen_macrophage) + coord_fixed()
umap_dataM <- Embeddings(Spleen_macrophage, reduction = "umap")
cells_to_keepM <- rownames(umap_dataM)[(umap_dataM[,2] >= 0) & (umap_dataM[,1] >= 7)]
Spleen_macrophage <- subset(Spleen_macrophage, cells = cells_to_keepM)
DimPlot(Spleen_macrophage) + coord_fixed()

Spleen_NK <- subset(spleen.combined.sct, idents = c("Activated_NK","Rest_NK"))
DimPlot(Spleen_NK)

# ###########################################
# perform analysis of percentage on each of the sub-cell type
# just chance the target_obj to anything 
# Result not look very good
#
# target_obj <- Spleen_NK
# celltype_var <- "celltype"  # can be adapted later if your T cell uses another label
# 
# meta_data <- target_obj@meta.data
# 
# # 2. Calculate % of each celltype per sample
# cell_counts <- meta_data %>%
#   group_by(sample, celltype) %>%
#   summarise(cell_number = n()) %>%
#   ungroup()
# 
# total_counts <- meta_data %>%
#   group_by(sample) %>%
#   summarise(total_cells = n())
# 
# cell_percentages <- left_join(cell_counts, total_counts, by = "sample") %>%
#   mutate(percent = 100 * cell_number / total_cells)
# 
# cell_percentages <- left_join(
#   cell_percentages, 
#   meta_data %>% distinct(sample, orig.ident, Hashtag), 
#   by = "sample"
# )
# 
# cell_percentages <- cell_percentages %>%
#   mutate(orig.ident = factor(orig.ident, 
#                              levels = c("spleen-control", "spleen_peptidoglycan"), 
#                              labels = c("Control", "Peptidoglycan")))
# 
# table(cell_percentages$celltype)
# 
# # 3. Get unique cell types
# celltype_list <- unique(meta_data[[celltype_var]]) %>% as.character()
# 
# # 4. Loop through each cell type and generate paired plots
# for (target_cell in celltype_list) {
#   print(target_cell)
#   
#   p <- ggpaired(
#     cell_percentages %>% filter(.data[[celltype_var]] == target_cell),
#     x = "orig.ident",
#     y = "percent",
#     id = "Hashtag",
#     line.color = "gray",
#     line.size = 0.9,
#     color = "orig.ident",
#     palette = "d3"
#   ) +
#     ggtitle(target_cell) +
#     ylab("Percentage") +
#     xlab("") +
#     scale_y_continuous(labels = label_percent(scale = 1)) +
#     stat_compare_means(paired = TRUE, method = "t.test")
#   
#   # Optional: adapt the output directory based on object name
#   obj_label <- deparse(substitute(target_obj))
#   
#   path_out <- paste0(
#     "../Manuscript/Figures/Figure 5-immune organoids and autoreactive B cells/Peptidoglycan_Cell_Percent/",
#     obj_label, ".", target_cell, ".pdf"
#   )
#   
#   ggsave(filename = path_out, plot = p, width = 3, height = 3.5, dpi = 300)
# }

#Setup DEG Analysis
Spleen_B$celltype2
Assays(Spleen_B)
DefaultAssay(Spleen_B) <- "RNA"

Idents(Spleen_B) <-"orig.ident" 
B_deg_results <- FindMarkers(Spleen_B,ident.2 = "spleen-control",ident.1 = "spleen_peptidoglycan",
  test.use = "MAST",logfc.threshold = log(1.2), min.pct = 0.05, verbose = TRUE)

B_deg_results$gene <- rownames(B_deg_results)
B_deg_results
filter(B_deg_results, gene %in% c("CD83"))

P.B.deg <- EnhancedVolcano(B_deg_results,
                                    lab = B_deg_results$gene,
                                    FCcutoff=0.15,
                                    xlim = c(-5.5, 5.5),
                                    x = 'avg_log2FC',
                                    y = 'p_val_adj', 
                                    title = "B cell")
P.B.deg


DefaultAssay(Spleen_B) <- "SCT"
FeaturePlot(Spleen_B, "CCL22", split.by="orig.ident",min.cutoff = 0, max.cutoff = 3)
FeaturePlot(Spleen_B, "IL1B", split.by="orig.ident",min.cutoff = 0, max.cutoff = 3)
FeaturePlot(Spleen_B, "EBI3", split.by="orig.ident",min.cutoff = 0, max.cutoff = 3)
FeaturePlot(spleen.combined.sct, "EBI3", split.by="orig.ident", min.cutoff = 0)

B_deg_results

#################################################################
Spleen_macrophage$celltype2
Assays(Spleen_macrophage)
DefaultAssay(Spleen_macrophage) <- "RNA"
Idents(Spleen_macrophage) <- "orig.ident"

M_deg_results <- FindMarkers(Spleen_macrophage, ident.2 = "spleen-control",ident.1 = "spleen_peptidoglycan",
  test.use = "MAST", logfc.threshold = log(1.2),min.pct = 0.05,verbose = TRUE)

M_deg_results$gene <- rownames(M_deg_results)

filter(M_deg_results, gene %in% c("CD83", "IL1B", "CCL22", "EBI3","CXCL2","CCL3","CCL4","IL6","TNF", "CD86", "MRC1", "CD163"))
filter(M_deg_results, gene %in% c( "CD86", "MRC1", "CD163","CD74","HLA-DRA","CXCL8"))

P.M.deg <- EnhancedVolcano(M_deg_results,lab = M_deg_results$gene,
  FCcutoff = 0.15, x = 'avg_log2FC',
  y = 'p_val_adj',
  xlim = c(-8, 8),
  title = "Macrophage")

P.M.deg

FeaturePlot(Spleen_macrophage, "CXCL5",split.by = "orig.ident", min.cutoff = 0, max.cutoff = 3, order = T)
VlnPlot(Spleen_macrophage, "CXCL5",split.by = "sample")
VlnPlot(Spleen_macrophage, "FBP1",split.by = "sample")

FeaturePlot(Spleen_macrophage, "MRC1",split.by = "orig.ident", min.cutoff = 0, max.cutoff = 5, order = T)
FeaturePlot(Spleen_macrophage, "CD86",split.by = "orig.ident", min.cutoff = 0, max.cutoff = 5, order = T)
FeaturePlot(spleen.combined.sct, features = c("F13A1", "IL1B", "TNF", "CD86", "MRC1", "CD163"), min.cutoff = 0)

#################################################################
Spleen_CD4T$celltype2
DimPlot(Spleen_CD4T)
FeaturePlot(Spleen_CD8T, "sct_IL22", max.cutoff = 1)
Assays(Spleen_CD4T)
DefaultAssay(Spleen_CD4T) <- "RNA"
Idents(Spleen_CD4T) <- "orig.ident"

T4_deg_results <- FindMarkers(Spleen_CD4T, ident.2 = "spleen-control",ident.1 = "spleen_peptidoglycan",
  test.use = "MAST", logfc.threshold = log(1.2), min.pct = 0.05,verbose = TRUE)

T4_deg_results$gene <- rownames(T4_deg_results)
T4_deg_results

filter(T4_deg_results, gene %in% c("CD40LG", "IFNG", "IL2", "FOXP3", "TNFRSF4", "IL17F", "IL22"))

P.T4.deg <- EnhancedVolcano(
  T4_deg_results,
  lab = T4_deg_results$gene,
  FCcutoff = 0.15,
  xlim = c(-5.5, 5.5),
  x = 'avg_log2FC',
  y = 'p_val_adj',
  title = "CD4 T Cells")

P.T4.deg

DefaultAssay(Spleen_CD4T) <- "SCT"

FeaturePlot(Spleen_CD4T, "sct_IL17F",   split.by = "orig.ident", min.cutoff = 0, max.cutoff = 3)
FeaturePlot(Spleen_CD4T, "IL2",    split.by = "orig.ident", min.cutoff = 0, max.cutoff = 3)
FeaturePlot(Spleen_CD4T, "FOXP3",  split.by = "orig.ident", min.cutoff = 0, max.cutoff = 3)
FeaturePlot(Spleen_CD4T, "TNFRSF4",split.by = "orig.ident", min.cutoff = 0, max.cutoff = 3)

##################################################################xc-Check Tregs an th17 from celltype, but the sig genes are both IL1B and CXCL8
##########check Treg
levels(Idents(spleen.combined.sct))

Spleen_Treg <- subset(spleen.combined.sct, idents = "X21_Treg")
DimPlot(Spleen_Treg) + coord_fixed()
umap_data4Treg <- Embeddings(Spleen_Treg, reduction = "umap")
cells_to_keep4Teg <- rownames(umap_data4Treg)[(umap_data4Treg[,2] <= -7) & (umap_data4Treg[,1] <= 4)]
Spleen_Treg <- subset(Spleen_Treg, cells = cells_to_keep4Teg)
DimPlot(Spleen_Treg) + coord_fixed()

##########check Th17
Spleen_TH17 <- subset(spleen.combined.sct, idents = "X13_Th17/MAIT")
DimPlot(Spleen_TH17) + coord_fixed()
umap_data4TH17 <- Embeddings(Spleen_TH17, reduction = "umap")
cells_to_keep4TH17 <- rownames(umap_data4TH17)[(umap_data4TH17[,2] <= -7) & (umap_data4TH17[,1] <= 4)]
Spleen_TH17 <- subset(Spleen_TH17, cells = cells_to_keep4TH17)
DimPlot(Spleen_TH17) + coord_fixed()

##########
Spleen_CD4T$celltype
Spleen_CD4T$orig.ident
DimPlot(Spleen_CD4T)

DefaultAssay(Spleen_Treg) <- "RNA"
Idents(Spleen_Treg) <- "orig.ident"

X21_Treg_results <- FindMarkers(Spleen_Treg, ident.2 = "spleen-control",ident.1 = "spleen_peptidoglycan",
                                test.use = "MAST", logfc.threshold = log(1.2), min.pct = 0.05,verbose = TRUE)

X21_Treg_results$gene <- rownames(X21_Treg_results)
X21_Treg_results

filter(X21_Treg_results, gene %in% c("CD40LG", "IFNG", "IL2", "FOXP3", "TNFRSF4", "IL17F", "IL22"))

X21_Treg.deg <- EnhancedVolcano(
  X21_Treg_results,
  lab = X21_Treg_results$gene,
  x = 'avg_log2FC',
  y = 'p_val_adj',
  FCcutoff = 0.15,
  pCutoff = 0.05,
  xlim = c(-5.5, 5.5),
  title = "X21 Treg (CD4+ Tregs)"
)


########For Th17 cells

DefaultAssay(Spleen_TH17) <- "RNA"
Idents(Spleen_TH17) <- "orig.ident"

TH17_results <- FindMarkers(Spleen_TH17, ident.2 = "spleen-control",ident.1 = "spleen_peptidoglycan",
                            test.use = "MAST", logfc.threshold = log(1.2), min.pct = 0.05,verbose = TRUE)

TH17_results$gene <- rownames(TH17_results)
TH17_results

filter(TH17_results, gene %in% c("CD40LG", "IFNG", "IL2", "FOXP3", "TNFRSF4", "IL17F", "IL22"))

TH17.deg <- EnhancedVolcano(
  TH17_results,
  lab = TH17_results$gene,
  x = 'avg_log2FC',
  y = 'p_val_adj',
  FCcutoff = 0.15,
  pCutoff = 0.05,
  xlim = c(-5.5, 5.5),
  title = "TH17"
)

TH17.deg




##################################################################xc
#################################################################
Spleen_CD8T$celltype2
Assays(Spleen_CD8T)
DefaultAssay(Spleen_CD8T) <- "RNA"
Idents(Spleen_CD8T) <- "orig.ident"

T8_deg_results <- FindMarkers(
  Spleen_CD8T,
  ident.2 = "spleen-control",
  ident.1 = "spleen_peptidoglycan",
  test.use = "MAST",
  logfc.threshold = log(1.2),  # log2(1.1) ~ 0.095
  min.pct = 0.05,
  verbose = TRUE)

T8_deg_results$gene <- rownames(T8_deg_results)

filter(T8_deg_results, gene %in% c("GZMB", "GZMK", "PRF1", "IFNG", "CXCR6", "IL22"))

P.T8.deg <- EnhancedVolcano(
  T8_deg_results,
  lab = T8_deg_results$gene,
  FCcutoff = 0.15,
  x = 'avg_log2FC',
  y = 'p_val_adj',
  xlim = c(-5.5, 5.5),
  title = "CD8 T Cells")

P.T8.deg

DefaultAssay(Spleen_CD8T) <- "SCT"

FeaturePlot(Spleen_CD8T, "GZMB",   split.by = "orig.ident", min.cutoff = 0, max.cutoff = 3)
FeaturePlot(Spleen_CD8T, "GZMK",   split.by = "orig.ident", min.cutoff = 0, max.cutoff = 3)
FeaturePlot(Spleen_CD8T, "PRF1",   split.by = "orig.ident", min.cutoff = 0, max.cutoff = 3)
FeaturePlot(Spleen_CD8T, "CXCR6",  split.by = "orig.ident", min.cutoff = 0, max.cutoff = 3)
FeaturePlot(Spleen_CD8T, "IFNG",   split.by = "orig.ident", min.cutoff = 0, max.cutoff = 3)

#################################################################
Spleen_NK$celltype2
Assays(Spleen_NK)
DefaultAssay(Spleen_NK) <- "RNA"
Idents(Spleen_NK) <- "orig.ident"

NK_deg_results <- FindMarkers(
  Spleen_NK,
  ident.2 = "spleen-control",
  ident.1 = "spleen_peptidoglycan",
  test.use = "MAST",
  logfc.threshold = log(1.2),
  min.pct = 0.05,
  verbose = TRUE)

NK_deg_results$gene <- rownames(NK_deg_results)

filter(NK_deg_results, gene %in% c("GZMB", "GZMK", "PRF1", "NKG7", "IFNG", "KLRD1", "IL22"))

P.NK.deg <- EnhancedVolcano(
  NK_deg_results,
  lab = NK_deg_results$gene,
  FCcutoff = 0.15,
  x = 'avg_log2FC',
  y = 'p_val_adj',
  xlim = c(-8, 8),
  title = "NK Cells")

P.NK.deg

DefaultAssay(Spleen_NK) <-"SCT"
FeaturePlot(Spleen_NK, "DUSP4",   split.by = "orig.ident", min.cutoff = 0, max.cutoff = 3)
FeaturePlot(Spleen_NK, "STAT1",   split.by = "orig.ident", min.cutoff = 0, max.cutoff = 3)
FeaturePlot(Spleen_NK, "CD160",   split.by = "orig.ident", min.cutoff = 0, max.cutoff = 3)
FeaturePlot(Spleen_NK, "IL22",   split.by = "orig.ident", min.cutoff = 0, max.cutoff = 3, order = T)
FeaturePlot(Spleen_NK, "NCAM1",   split.by = "orig.ident", min.cutoff = 0, max.cutoff = 3, order = T)
FeaturePlot(Spleen_NK, "TNF",   split.by = "orig.ident", min.cutoff = 0, max.cutoff = 3, order = T)

#################################################################
DefaultAssay(spleen.combined.sct) <- "RNA"

Idents(spleen.combined.sct) <- "orig.ident"

#This step take 8 hrs
# ALL_deg_results <- FindMarkers(
#   spleen.combined.sct,
#   ident.1 = "spleen-control",
#   ident.2 = "spleen_peptidoglycan",
#   test.use = "MAST",
#   logfc.threshold = log(1.1),
#   min.pct = 0.01,
#   verbose = TRUE)
# 
# ALL_deg_results$gene <- rownames(ALL_deg_results)

filter(ALL_deg_results, gene == "IL22")
VlnPlot(spleen.combined.sct, "sct_IL22", split.by = "sample")
FeaturePlot(spleen.combined.sct, "IL22", split.by = "orig.ident", max.cutoff = 1)

FeaturePlot(spleen.combined.sct, "ERVK3-1", split.by = "orig.ident", max.cutoff = 3)


All_deg_vocano <- EnhancedVolcano(
  ALL_deg_results,
  lab = ALL_deg_results$gene,
  FCcutoff = 0.15,
  #xlim = c(-5.5, 5.5),
  x = 'avg_log2FC',
  y = 'p_val_adj',
  title = "All")

All_deg_vocano

#write.csv(ALL_deg_results, file = "../Manuscript/Figures/Figure 5-immune organoids and autoreactive B cells/ALL_DEG.csv")

B_deg_results$celltype <- "B_Cell"
T4_deg_results$celltype <- "CD4T_Cell"
T8_deg_results$celltype <- "CD8T_Cell"
M_deg_results$celltype <- "Macrophage_DC"
NK_deg_results$celltype <- "NK_Cell"

spleen_DEG <- rbind(B_deg_results,T4_deg_results,T8_deg_results,M_deg_results,NK_deg_results)

#save(spleen.combined.sct,spleen_DEG, file = "./RObject/03.DEGfinal.RData")
#write.csv(file = "../Manuscript/Figures/Figure 5-immune organoids and autoreactive B cells/Spleen_DEG.csv", spleen_DEG)

spleen_DEG[str_detect(spleen_DEG$gene, "ERV"),]
spleen_DEG[str_detect(spleen_DEG$gene, "MDA5"),]
spleen_DEG[str_detect(spleen_DEG$gene, "IFIH1"),]

ALL_DEG[str_detect(ALL_DEG$gene, "ERV"),]
ALL_DEG[str_detect(ALL_DEG$gene, "MDA5"),]
ALL_DEG[str_detect(ALL_DEG$gene, "IFIH1"),]

# DNA damage and repair genes (example subset)
ddr_genes <- c("TP53", "CDKN1A", "GADD45A", "ATM", "ATR", "BRCA1", "BRCA2", "RAD51", 
               "MSH2", "MSH6", "MLH1", "PMS2", "PARP1", "XRCC1")

# APOBEC family
apobec_genes <- c("APOBEC3A", "APOBEC3B", "APOBEC3C", "APOBEC3D", "APOBEC3F", "APOBEC3G")

# Create named lists
mutation_sets <- list(DDR = ddr_genes, APOBEC = apobec_genes)

spleen.combined.sct <- AddModuleScore(
  spleen.combined.sct,
  features = mutation_sets,
  name = "MutationScore")

VlnPlot(spleen.combined.sct,
        features = c("MutationScore1", "MutationScore2"),
        group.by = "orig.ident",
        pt.size = 0) + stat_compare_means(label = "p.format", method = "wilcox.test")

meta_data <- spleen.combined.sct@meta.data
wilcox.test(MutationScore1 ~ orig.ident, data = meta_data)
wilcox.test(MutationScore2 ~ orig.ident, data = meta_data)

celltypes <- unique(meta_data$celltype2)

results <- lapply(celltypes, function(ct) {
  df <- meta_data %>% filter(celltype2 == ct)
  test <- wilcox.test(MutationScore1 ~ orig.ident, data = df)
  data.frame(celltype = ct, p_value = test$p.value)
})

do.call(rbind, results)

results <- lapply(celltypes, function(ct) {
  df <- meta_data %>% filter(celltype2 == ct)
  test <- wilcox.test(MutationScore2 ~ orig.ident, data = df)
  data.frame(celltype = ct, p_value = test$p.value)
})

do.call(rbind, results)


#calculate somatic mutation score
celltypes <- unique(meta_data$celltype2)

# For MutationScore1
results1 <- lapply(celltypes, function(ct) {
  df <- meta_data %>% filter(celltype2 == ct)
  
  if (length(unique(df$orig.ident)) > 1) {
    test <- wilcox.test(MutationScore1 ~ orig.ident, data = df)
    mean_values <- df %>%
      group_by(orig.ident) %>%
      summarise(mean_score = mean(MutationScore1, na.rm = TRUE))
    
    # Ensure ordering: control vs treated
    control_mean <- mean_values$mean_score[mean_values$orig.ident == "spleen-control"]
    treated_mean <- mean_values$mean_score[mean_values$orig.ident == "spleen_peptidoglycan"]
    
    log2_fc <- log2(treated_mean / control_mean)
    
    data.frame(celltype = ct,
               p_value = test$p.value,
               log2FC = log2_fc)
  } else {
    data.frame(celltype = ct, p_value = NA, log2FC = NA)
  }
})

# For MutationScore2
results2 <- lapply(celltypes, function(ct) {
  df <- meta_data %>% filter(celltype2 == ct)
  
  if (length(unique(df$orig.ident)) > 1) {
    test <- wilcox.test(MutationScore2 ~ orig.ident, data = df)
    mean_values <- df %>%
      group_by(orig.ident) %>%
      summarise(mean_score = mean(MutationScore2, na.rm = TRUE))
    
    control_mean <- mean_values$mean_score[mean_values$orig.ident == "spleen-control"]
    treated_mean <- mean_values$mean_score[mean_values$orig.ident == "spleen_peptidoglycan"]
    
    log2_fc <- log2(treated_mean / control_mean)
    
    data.frame(celltype = ct,
               p_value = test$p.value,
               log2FC = log2_fc)
  } else {
    data.frame(celltype = ct, p_value = NA, log2FC = NA)
  }
})

# Combine and label results
MutationScore1_result <- do.call(rbind, results1) %>% mutate(score = "MutationScore1")
MutationScore2_result <- do.call(rbind, results2) %>% mutate(score = "MutationScore2")

combined_results <- rbind(MutationScore1_result, MutationScore2_result)
combined_results


meta_data <- spleen.combined.sct@meta.data

# Define score types and cell types
score_types <- c("MutationScore1", "MutationScore2")
celltypes <- unique(meta_data$celltype2)

# Run Wilcoxon test and compute fold change for each score type and cell type
results <- lapply(score_types, function(score) {
  lapply(celltypes, function(ct) {
    df <- meta_data %>% filter(celltype2 == ct)
    
    # Perform Wilcoxon test
    test <- wilcox.test(as.formula(paste(score, "~ orig.ident")), data = df)
    
    # Compute fold change (mean difference between groups)
    fc <- df %>%
      group_by(orig.ident) %>%
      summarise(mean = mean(.data[[score]], na.rm = TRUE)) %>%
      summarise(log2FC = log2(mean[1] / mean[2])) %>%
      pull(log2FC)
    
    data.frame(score = score, celltype = ct, p_value = test$p.value, log2FC = fc)
  }) %>% bind_rows()
}) %>% bind_rows()

# Convert to heatmap-like plot
results$negLog10P <- -log10(results$p_value)

results_adj <- results %>%
  mutate(log2FC = ifelse(p_value > 0.05 | is.na(log2FC), 0, log2FC))

ggplot(results_adj, aes(x = score, y = celltype, fill = log2FC, size = -log10(p_value))) +
  geom_point(shape = 21, color = "black") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0, na.value = "grey90") +
  scale_size_continuous(range = c(1, 8)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(
    x = "Mutation Score",
    y = "Cell Type",
    fill = "log2FC",
    size = "-log10(p-value)",
    title = "Mutation Score Heatmap")

VlnPlot(spleen.combined.sct, features = "PMS2", group.by = "celltype2")

shm_genes <- list(c("AICDA", "UNG", "MSH2", "MSH6", "EXO1", "PMS2"))
spleen.combined.sct <- AddModuleScore(spleen.combined.sct, features = shm_genes, name = "SHM_Score")
VlnPlot(spleen.combined.sct, features = "SHM_Score1", group.by = "orig.ident") + facet_wrap(~ celltype2, scales = "free_y")

FetchData(spleen.combined.sct, vars = c("SHM_Score1", "orig.ident", "celltype2")) %>% 
  ggplot(aes(x = orig.ident, y = SHM_Score1, fill = celltype2)) +
  geom_violin(trim = FALSE) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.1) +
  facet_wrap(~ celltype2, scales = "free_y") +
  theme_minimal() +
  labs(x = "Sample", y = "SHM Score", title = "SHM Score by Sample and Cell Type") + scale_fill_manual(values = pep_cluster2_color)

SHM_meta_data <- spleen.combined.sct@meta.data %>% dplyr::select(orig.ident, celltype2, SHM_Score1)

celltypes <- unique(SHM_meta_data$celltype2)

shm_results <- lapply(celltypes, function(ct) {
  df <- SHM_meta_data %>% filter(celltype2 == ct)
  
  test <- wilcox.test(SHM_Score1 ~ orig.ident, data = df)
  
  fc <- df %>%
    group_by(orig.ident) %>%
    summarise(mean_score = mean(SHM_Score1), .groups = "drop") %>%
    summarise(log2FC = log2(mean_score[2]/mean_score[1])) %>%
    pull(log2FC)
  
  data.frame(celltype = ct, p_value = test$p.value, log2FC = fc)
})

shm_results
shm_df <- do.call(rbind, lapply(shm_results, as.data.frame))
shm_df$p_adj <- p.adjust(shm_df$p_value, method = "BH")
shm_df %>% arrange(p_adj)

FeaturePlot(spleen.combined.sct, features = "CXCL8", split.by = "orig.ident")







