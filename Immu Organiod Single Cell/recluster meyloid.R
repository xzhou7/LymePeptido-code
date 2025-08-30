



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
library(dplyr)
library(tidyr)


library(Seurat)

#source("./00_colorKey_XC.R")
#mac
setwd("~/Library/CloudStorage/Dropbox/lyme_disease/Pepti_Spleen_singlecell/")
getwd()
#load data
load(file = "./RObject/02.Celltype_labelled_UMAP.RData")

#DimPlot(spleen.combined.sct,label = T)
#Idents(spleen.combined.sct) <- "seurat_clusters"
#levels(Idents(spleen.combined.sct))
#myeloid <- subset(spleen.combined.sct, idents = c("26", "1"))

#try to figure out what is x26
#markers_X26 <- FindMarkers(myeloid, ident.1 = "26", logfc.threshold = 0.25, test.use = "wilcox")
#markers_X26 <- markers_X26[markers_X26$avg_log2FC > 0, ]
#write.csv(markers_X26, file = "./Result/markers_X26.csv")
#up_sig_26 <- markers_X26 %>%
#  dplyr::filter(avg_log2FC > 0, p_val_adj < 0.05) %>%
#  dplyr::arrange(desc(avg_log2FC)) %>%
#  head(50)  # Top 30 by log2FC
#write.csv(up_sig_26, file = "./Result/markers_X26.csv")


setwd("~/Library/CloudStorage/Dropbox/lyme_disease/Pepti_Spleen_singlecell/")
getwd()
library(Seurat)
source("./Code/00_colorKey_XC.R")
load("./RObject/04.myeloid_cells.RData")
levels(Idents(myeloid))

library(Seurat)
library(ggplot2)
library(dplyr)

# Set ADT as default for demuxing
DefaultAssay(myeloid) <- "ADT"
VlnPlot(myeloid, c("Hasgtag 1-TotalSeqC", "Hasgtag 2-TotalSeqC", "Hasgtag 3-TotalSeqC"))

# Assign hashtags
myeloid[["percent.1"]] <- PercentageFeatureSet(myeloid, pattern = "^Hasgtag 1")
myeloid[["percent.2"]] <- PercentageFeatureSet(myeloid, pattern = "^Hasgtag 2")
myeloid[["percent.3"]] <- PercentageFeatureSet(myeloid, pattern = "^Hasgtag 3")

myeloid$Hashtag <- ifelse(myeloid$percent.1 > myeloid$percent.2 & myeloid$percent.1 > myeloid$percent.3, "Hashtag1",
                          ifelse(myeloid$percent.2 > myeloid$percent.1 & myeloid$percent.2 > myeloid$percent.3, "Hashtag2", "Hashtag3"))
table(myeloid$Hashtag)

# Function to process each donor (start from SCTransform)
process_donor <- function(seurat_obj) {
  DefaultAssay(seurat_obj) <- "RNA"
  
  seurat_obj <- SCTransform(seurat_obj, vst.flavor = "v2", verbose = FALSE)
  seurat_obj <- RunPCA(seurat_obj, npcs = 15, verbose = FALSE)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:15)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:15)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.8)
  
  return(seurat_obj)
}

# Subset by donor and process
pep_donor1 <- process_donor(subset(myeloid, Hashtag == "Hashtag1"))
pep_donor2 <- process_donor(subset(myeloid, Hashtag == "Hashtag2"))
pep_donor3 <- process_donor(subset(myeloid, Hashtag == "Hashtag3"))

# Combine into list
myeloid.list <- list(donor1 = pep_donor1,
                     donor2 = pep_donor2,
                     donor3 = pep_donor3)

features <- SelectIntegrationFeatures(object.list = myeloid.list, nfeatures = 4000)

myeloid.list <- PrepSCTIntegration(object.list = myeloid.list, anchor.features = features)

# This step may take time
myeloid.anchors <- FindIntegrationAnchors(
  object.list = myeloid.list,
  normalization.method = "SCT",
  anchor.features = features,
  reference = c(1, 2, 3)
)

myeloid.combined.sct <- IntegrateData(
  anchorset = myeloid.anchors,
  normalization.method = "SCT"
)

myeloid.combined.sct <- FindVariableFeatures(
  myeloid.combined.sct,
  selection.method = "vst",
  nfeatures = 4000
)

myeloid.combined.sct

save(myeloid.combined.sct, file = "./RObject/Integrated_4000Feature_myeloid.RData")

setwd("~/Library/CloudStorage/Dropbox/lyme_disease/Pepti_Spleen_singlecell/")
getwd()
library(Seurat)
library(ggplot2)
library(dplyr)
source("./Code/00_colorKey_XC.R")

load("./RObject/Integrated_4000Feature_myeloid.RData")
myeloid.combined.sct <- RunPCA(myeloid.combined.sct, verbose = FALSE)
myeloid.combined.sct <- RunUMAP(myeloid.combined.sct, dims = 1:30, min.dist = 0.7, verbose = FALSE)
myeloid.combined.sct <- FindNeighbors(myeloid.combined.sct, dims = 1:30, verbose = FALSE)
myeloid.combined.sct <- FindClusters(myeloid.combined.sct,resolution = 0.6, verbose = FALSE)

myeloid.combined.sct$orig.ident <- factor(myeloid.combined.sct$orig.ident, levels = c("spleen-control", "spleen_peptidoglycan"))
#table(myeloid.combined.sct$orig.ident)

p1 <- DimPlot(myeloid.combined.sct, label = T, split.by = "orig.ident") + NoLegend() + coord_fixed()
p2 <- DimPlot(myeloid.combined.sct, label = T, split.by = "sample") + NoLegend() + coord_fixed()
p2
ggsave("./Result/myeloid/myeloid_dimplot_conditions.pdf", plot = p1, width = 5, height = 3)
ggsave("./Result/myeloid/myeloid_dimplot_samples.pdf", plot = p2, width = 10, height = 3)

# Set the identity to clusters if not already
#Idents(myeloid.combined.sct) <- "seurat_clusters"

# Run marker detection for all clusters
#all_markers <- FindAllMarkers(
#  object = myeloid.combined.sct,
#  only.pos = TRUE,             # Only return upregulated genes
#  min.pct = 0.25,              # Expressed in at least 25% of cells in either cluster
#  logfc.threshold = 0.25,      # Minimum log2FC
#  test.use = "wilcox"          # Wilcoxon test
#)

#write.csv(all_markers, file = "./Result/myeloid/Myeloid_AllCluster_Markers.csv")
#top30 <- all_markers %>% 
#  group_by(cluster) %>% 
#  top_n(n = 30, wt = avg_log2FC)
#write.csv(top30, file = "./Result/myeloid/Myeloid_AllCluster_Markers_top30.csv")


levels(myeloid.combined.sct)


######rename cluster
new_cluster_ids <- c(
  "0" = "memory B cells",
  "1" = "cytotoxic T cells",
  "2" = "monocyte-derived inflam macrophages",
  "3" = "inflam tissue macrophages",
  "4" = "anti-inflam tissue macrophages",
  "5" = "differentiating monocytes",
  "6" = "dendritic cells",
  "7" = "hematopoietic stem and progenitor cells",
  "8" = "proliferating myeloid cells"
)


#myeloid.combined.sct <- RenameIdents(myeloid.combined.sct, new_cluster_ids)
#myeloid.combined.sct$celltype <- Idents(myeloid.combined.sct)
#DimPlot(myeloid.combined.sct, label = TRUE, group.by = "celltype") + NoLegend()
#save(myeloid.combined.sct, file = "./RObject/Integrated_myeloid_renamed.RData")

library(ggplot2)
library(ggrepel)
library(dplyr)
source("./Code/00_colorKey_XC.R")
p_umap <- DimPlot(myeloid.combined.sct, label =F, raster=FALSE) & coord_equal() 
p_umap
umap_df <- as.data.frame(Embeddings(myeloid.combined.sct, "umap"))
umap_df$celltype <- myeloid.combined.sct@meta.data$celltype

# Use correct column names: umap_1 and umap_2
celltype_position <- umap_df %>%
  group_by(celltype) %>%
  summarise(
    umap_1 = median(umap_1),
    umap_2 = median(umap_2),
    .groups = "drop"
  )

p_umap2 <- ggplot(umap_df, aes(x = umap_1, y = umap_2, color = celltype)) +
  geom_point(size = 0.5, alpha = 0.8) +
  geom_label_repel(
    data = celltype_position,
    aes(x = umap_1, y = umap_2, label = celltype, color = celltype),
    fontface = "bold",
    fill = "white",
    point.padding = unit(0.1, "lines"),
    alpha = 0.95
    # Removed: show.legend = FALSE
  ) +
  theme_void() +
  theme(legend.position = "right") +
  scale_color_manual(values = spleen_myeloid_color)

p_umap2 <- ggplot(umap_df, aes(x = umap_1, y = umap_2, color = celltype)) +
  geom_point(size = 0.5, alpha = 0.8) +
  geom_label_repel(
    data = celltype_position,
    aes(x = umap_1, y = umap_2, label = celltype, color = celltype),
    fontface = "bold",
    fill = "white",
    point.padding = unit(0.1, "lines"),
    alpha = 0.95,
    show.legend = FALSE
  ) +
  theme_void() +
  theme(legend.position = "none") +
  scale_color_manual(values = spleen_myeloid_color)

p_umap2

pdot<-DotPlot(myeloid.combined.sct, features = c("BANK1", "NAPSB","CD19", "CD79B","PAX5", "CD247", "CD7","GZMA","TRBC1",
  "ATOX1","FTH1","LY96","TLR4","SEMA4A","TRPM2","DUSP3", 
  "WARS1","S100A11", "VIM", "FCGR1A", 
  "CST3", "NPC2", "MT-ATP6", "TGFBI", 
  "DPYD","LRMDA", "ITPR2","ZEB2", 
  "LAMP3", "CRACD","RARRES2", "PNRC1",
  "GFI1B", "TPSAB1","GATA2", "GCSAML", 
 "FTL","MYL9","YBX1","SPI1")) + RotatedAxis()

pdot <- pdot + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_color_gradient2(low = "#D3D3D3", high = "#8C1515")
pdot
ggsave(filename = "../Manuscript/Figures/Figure3-autoimmune monocytes/dot_marker_mono.pdf",pdot, width = 10, height = 4, dpi = 300)





# Set identity classes to cluster assignments if not already set
Idents(myeloid.combined.sct) <- "seurat_clusters"

# Find markers for cluster 4 vs all other clusters
cluster4.markers <- FindMarkers(myeloid.combined.sct, 
                                ident.1 = "4", 
                                only.pos = TRUE, 
                                logfc.threshold = 0.25, 
                                test.use = "MAST")

# View top markers
head(cluster4.markers, 20)

# Find markers for cluster 2 vs all other clusters
cluster2.markers <- FindMarkers(myeloid.combined.sct, 
                                ident.1 = "2", 
                                only.pos = TRUE, 
                                logfc.threshold = 0.25, 
                                test.use = "MAST")

# View top markers
head(cluster2.markers, 20)

# Save to file
write.csv(cluster2.markers, file = "./Result/Integrated_4000Feature_myeloid_Cluster2_Markers.csv")

VlnPlot(myeloid.combined.sct, features = c("CD14", "PSAP", "CTSD", "MRC1", "CD163", "CD68", "CXCL8"))

VlnPlot(myeloid.combined.sct, features = c("CD14", "LY96", "FCGR3A", "S100A8", "S100A9", "CD68", "CSF1R", "CD163", "MRC1", "FOLR2", "TREM2"))


monocyte_genes <- c("CD14", "LY96", "S100A8", "S100A9", "CCR2")
macrophage_genes <- c("CD68", "CD163", "MRC1", "TREM2", "APOE", "FOLR2")

myeloid.combined.sct <- AddModuleScore(myeloid.combined.sct, features = list(monocyte_genes), name = "Monocyte_Score")
myeloid.combined.sct <- AddModuleScore(myeloid.combined.sct, features = list(macrophage_genes), name = "Macrophage_Score")

DefaultAssay(myeloid.combined.sct) <- "RNA"
DefaultAssay(myeloid.combined.sct) <- "integrated"

FeaturePlot(myeloid.combined.sct, features = "LAMP3")
FeaturePlot(myeloid.combined.sct, features = "EBF1", split.by = "orig.ident")

myeloid.combined.sct$orig.ident
DotPlot(myeloid.combined.sct, features = list(
  Monocyte = c("CD14", "LY96", "S100A8", "S100A9", "CCR2"),
  Macrophage = c("CD68", "CD163", "MRC1", "TREM2", "APOE")
)) + RotatedAxis()






















# Step 1: Normalize if not already done (using SCTransform is recommended)
myeloid <- SCTransform(myeloid, verbose = FALSE, method = "glmGamPoi")
# Step 2: Run PCA
myeloid <- RunPCA(myeloid, verbose = FALSE)
# Step 4: Find neighbors and clusters
myeloid <- FindNeighbors(myeloid, dims = 1:30)
myeloid <- FindClusters(myeloid, resolution = 0.5)  # adjust resolution if needed
myeloid_markers <- FindAllMarkers(myeloid,
                                  assay = "SCT",            # or "RNA" if you're using raw normalized data
                                  only.pos = TRUE,          # only keep upregulated markers
                                  min.pct = 0.25,           # expressed in at least 25% of cells
                                  logfc.threshold = 0.25)   # minimum log fold change
# View top 5 markers per cluster
#top20_meyloid_marker <- myeloid_markers %>%
#  group_by(cluster) %>%
#  top_n(n = 20, wt = avg_log2FC)
#write.csv(top20_meyloid_marker, "./Result/top20_meyloid_marker.csv", row.names = FALSE)
# Step 5: Run UMAP or tSNE
#myeloid <- RunUMAP(myeloid, dims = 1:30)  # or use harmony_dims if Harmony used
# Step 6: Visualize
#DimPlot(myeloid, label = TRUE, group.by = "seurat_clusters") + NoLegend()
#save(myeloid, file = "./RObject/04.myeloid_cells_annota.RData")



#load
load("./RObject/04.myeloid_cells_annota.RData")
DimPlot(myeloid, split.by = "orig.ident", label = T)
DimPlot(myeloid, split.by = "sample", label = T)
table(myeloid$sample)
DimPlot(myeloid,
        reduction = "umap",
        group.by = "celltype",       # or "seurat_clusters" if you prefer cluster numbers
        split.by = "orig.ident",      # the column indicating treatment
        label = TRUE,
        repel = TRUE) +
  NoLegend()

#remove cluster 0 as it is B Cells
# Check current cluster identities
# Ensure identities are set correctly
Idents(myeloid) <- "seurat_clusters"

# Remove clusters "0" and "3"
myeloid.filtered <- subset(myeloid, idents = setdiff(levels(Idents(myeloid)), c("0", "3")))

levels(Idents(myeloid.filtered))
#Step 1: Normalize if not already done (using SCTransform is recommended)
myeloid.filtered <- SCTransform(myeloid.filtered, verbose = FALSE, method = "glmGamPoi")
# Step 2: Run PCA
myeloid.filtered <- RunPCA(myeloid.filtered, verbose = FALSE)
# Step 3: (Optional) Use Harmony to correct for batch effects if applicable
library(harmony)

myeloid.filtered <- RunHarmony(myeloid.filtered, group.by.vars = "sample")  # replace with your batch variable
myeloid.filtered@reductions$pca <- myeloid.filtered@reductions$harmony  # overwrite PCA with Harmony
# Step 4: Find neighbors and clusters
myeloid.filtered <- FindNeighbors(myeloid.filtered, dims = 1:30)
myeloid.filtered <- FindClusters(myeloid.filtered, resolution = 0.5)  # adjust resolution if needed
myeloid.filtered_markers <- FindAllMarkers(myeloid.filtered,
                                  assay = "SCT",            # or "RNA" if you're using raw normalized data
                                  only.pos = TRUE,          # only keep upregulated markers
                                  min.pct = 0.25,           # expressed in at least 25% of cells
                                  logfc.threshold = 0.25)   # minimum log fold change

write.csv(myeloid.filtered_markers,
          file = "/Users/xinchen/Library/CloudStorage/Dropbox/lyme_disease/Pepti_Spleen_singlecell/Result/myeloid_filtered_markers.csv",
          row.names = FALSE)

Myeloid.clean <-myeloid.filtered
count2.2 <- data.frame(table(Myeloid.clean$sample, Myeloid.clean$seurat_clusters))
total.bysubject <- data.frame(table(Myeloid.clean$sample))
colnames(total.bysubject)[2] <- "Total"
count2.2$conditon <- "Control"
count2.2$conditon[str_detect(count2.2$Var1, "Pepti")] <- "Peptidoglycan"
count2.2 <- merge(count2.2, total.bysubject, by="Var1")
count2.2$Percent <- count2.2$Freq/count2.2$Total

count2.2$conditon <- factor(count2.2$conditon, levels = c("Control", "Peptidoglycan"))

comparasions <- list("Control", "Peptidoglycan")
library(ggplot2)


p2.2 <- ggplot(count2.2, aes(x=conditon, y=Percent, color=conditon)) + geom_boxplot(alpha=0.7, outlier.alpha = 0, size = 0.7)
library(cowplot)
p2.2 <- p2.2 + facet_wrap(.~Var2, scales = "free", nrow = 1) + theme_cowplot()
p2.2 <- p2.2 + geom_jitter(color="black", size=0.4, alpha=0.9)
library(ggpubr)

p2.2 <- p2.2 + stat_compare_means(method="t.test", comparisons = comparasions) + scale_color_manual(values=condition_color) + XZ_flip_x_label()
p2.2 <- ggplot(count2.2, aes(x = conditon, y = Percent, color = conditon)) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.9) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0, width = 0.6) +
  facet_wrap(. ~ Var2, scales = "free", nrow = 1) +
  stat_compare_means(method = "t.test", comparisons = comparasions) +
  scale_color_manual(values = condition_color_spleen) +
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












##Idents(myeloid)  # Show current identity class

#Idents(myeloid) <- "seurat_clusters"  # Ensure identities are set from metadata
#levels(Idents(myeloid))  # Should list cluster IDs like "0", "1", ...
#myeloid.filtered <- subset(myeloid, idents = setdiff(levels(Idents(myeloid)), "0"))
DimPlot(myeloid.filtered, split.by = "orig.ident", label = T)
DimPlot(myeloid.filtered, split.by = "sample", label = T)
table(myeloid.filtered$sample)

######rename cluster
# Define a named vector with new cluster names
#new_cluster_ids <- c(
#  "1" = "Activated monocytes/macrophages1",
#  "2" = "Homeostatic Macrophages",
#  "4" = "Activated monocytes/macrophages2",
#  "5" = "Activated monocytes/macrophages3",
#  "6" = "mature LAMP3⁺ DCs",
#  "7" = "megakaryocyte–erythroid progenitors",
#  "8" = "Classical Monocytes"
#)

# Apply the renaming
#myeloid <- RenameIdents(myeloid, new_cluster_ids)
#myeloid$celltype <- Idents(myeloid)

#levels(Idents(myeloid))

#DimPlot(myeloid, label = TRUE, repel = TRUE) + NoLegend()
#DimPlot(myeloid, label = TRUE, repel = TRUE, seperate) + NoLegend()
#######need a dorplot

# Load marker list (if not already loaded)
# Marker genes per cluster
#marker_genes <- c(
#  "S100A8", "S100A9", "LYZ",            # Classical Monocytes
#  "FCGR3A", "TREM2", "SIGLEC10",        # Non-Classical Monocytes
#  "HLA-DRA", "CD74",                    # Antigen-Presenting Monocytes
#  "FSCN1", "CLEC9A", "XCR1",            # cDC1
#  "CD1C", "CLEC10A",                    # cDC2
#  "IL1B", "CXCL8", "TNF",               # Inflammatory Monocytes
#  "MKI67", "CENPF", "TYMS",             # Proliferating Myeloid Cells
#  "LAMP3", "CCR7", "IDO1",              # Mature/Migratory DCs
#  "PHLDA2", "S100A8"                    #activation-stressed
#)

myeloid$celltype <- Idents(myeloid)  # Save renamed identities as metadata column


DotPlot(myeloid, 
        features = marker_genes, 
        group.by = "celltype") +
  RotatedAxis() +
  theme(axis.text.x = element_text(size = 8)) +
  ggtitle("Dot Plot of Key Myeloid Markers")


ggsave("Myeloid_Marker_DotPlot.pdf", width = 10, height = 6)

#####

###################barplot-not significant
#mac
setwd("~/Library/CloudStorage/Dropbox/lyme_disease/Pepti_Spleen_singlecell/")
getwd()
library(Seurat)
library(dplyr)

source("./Code/00_colorKey_XC.R")

#load
load("./RObject/04.myeloid_cells_annota.RData")
Idents(myeloid) <- "seurat_clusters"

# Remove clusters "0" and "3"
Myeloid.clean <- subset(myeloid, idents = setdiff(levels(Idents(myeloid)), c("0", "3")))

library(stringr)

count2.2 <- data.frame(table(Myeloid.clean$sample, Myeloid.clean$seurat_clusters))
total.bysubject <- data.frame(table(Myeloid.clean$sample))
colnames(total.bysubject)[2] <- "Total"
count2.2$conditon <- "Control"
count2.2$conditon[str_detect(count2.2$Var1, "Pepti")] <- "Peptidoglycan"
count2.2 <- merge(count2.2, total.bysubject, by="Var1")
count2.2$Percent <- count2.2$Freq/count2.2$Total

count2.2$conditon <- factor(count2.2$conditon, levels = c("Control", "Peptidoglycan"))

comparasions <- list("Control", "Peptidoglycan")
library(ggplot2)


p2.2 <- ggplot(count2.2, aes(x=conditon, y=Percent, color=conditon)) + geom_boxplot(alpha=0.7, outlier.alpha = 0, size = 0.7)
library(cowplot)
p2.2 <- p2.2 + facet_wrap(.~Var2, scales = "free", nrow = 1) + theme_cowplot()
p2.2 <- p2.2 + geom_jitter(color="black", size=0.4, alpha=0.9)
library(ggpubr)

p2.2 <- p2.2 + stat_compare_means(method="t.test", comparisons = comparasions) + scale_color_manual(values=condition_color) + XZ_flip_x_label()
p2.2 <- ggplot(count2.2, aes(x = conditon, y = Percent, color = conditon)) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.9) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0, width = 0.6) +
  facet_wrap(. ~ Var2, scales = "free", nrow = 1) +
  stat_compare_means(method = "t.test", comparisons = comparasions) +
  scale_color_manual(values = condition_color_spleen) +
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

ggsave(filename = "./Result/Figure_Myeloid_Cluster_Percent.pdf",
       plot = p2.2,
       width = 10,         # adjust as needed
       height = 3,         # adjust as needed
       units = "in")       # can also use "cm"

############ X1
# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)

# Subset Cluster 1
X1 <- subset(Myeloid.clean, seurat_clusters == "1")
DimPlot(X1)  # Optional visualization

# Set condition as identity for comparison
Idents(X1) <- "orig.ident"
table(Idents(X1))  # Check identities: should include 'spleen-control' and 'spleen_peptidoglycan'

# Differential expression: Peptidoglycan vs Control using MAST
Peptidoglycan_X1 <- FindMarkers(
  object = X1,
  ident.1 = "spleen_peptidoglycan",
  ident.2 = "spleen-control",
  logfc.threshold = log(1.1),
  test.use = "MAST",
  min.pct = 0.3
)

# Prepare DEG result
Peptidoglycan_X1$gene <- rownames(Peptidoglycan_X1)
Peptidoglycan_X1$log2FC <- Peptidoglycan_X1$avg_log2FC  # For convenience

# Annotate significance
Peptidoglycan_X1 <- Peptidoglycan_X1 %>%
  mutate(
    sig = case_when(
      p_val_adj < 0.05 & log2FC > 0 ~ "Up",
      p_val_adj < 0.05 & log2FC < 0 ~ "Down",
      TRUE ~ "NS"
    )
  )

# Set color palette
volcano_colors <- c("Up" = "#E41A1C", "Down" = "#377EB8", "NS" = "gray80")

# Volcano plot with gene labels for significant DEGs
volcano_plot <- ggplot(Peptidoglycan_X1, aes(x = log2FC, y = -log10(p_val_adj), color = sig)) +
  geom_point(alpha = 0.8, size = 1.5) +
  geom_text_repel(
    data = filter(Peptidoglycan_X1, sig != "NS"),
    aes(label = gene),
    size = 3,
    max.overlaps = 100
  ) +
  scale_color_manual(values = volcano_colors) +
  geom_vline(xintercept = c(-log(1.1), log(1.1)), linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  labs(
    title = "Volcano Plot: Peptidoglycan vs Control (Cluster X1)",
    x = "log2 Fold Change",
    y = "-log10 Adjusted p-value"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.title = element_blank())

# Display the plot
print(volcano_plot)

ggsave("./Result/Volcano_X1_Peptidoglycan_vs_Control.pdf", plot = volcano_plot, width = 6, height = 5)


############
############
library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)

# Define function to generate and save volcano plots and DEG tables
make_volcano_plot <- function(seurat_obj, cluster_id, output_dir = ".") {
  # Subset cluster
  cluster_obj <- subset(seurat_obj, seurat_clusters == cluster_id)
  Idents(cluster_obj) <- "orig.ident"
  
  # Run differential expression
  deg <- FindMarkers(
    object = cluster_obj,
    ident.1 = "spleen_peptidoglycan",
    ident.2 = "spleen-control",
    logfc.threshold = log(1.1),
    test.use = "MAST",
    min.pct = 0.3
  )
  
  # Annotate and label
  deg$gene <- rownames(deg)
  deg$log2FC <- deg$avg_log2FC
  deg <- deg %>%
    mutate(
      sig = case_when(
        p_val_adj < 0.05 & log2FC > 0 ~ "Up",
        p_val_adj < 0.05 & log2FC < 0 ~ "Down",
        TRUE ~ "NS"
      )
    )
  
  # Save DEG table to CSV (including gene, avg_log2FC, p_val_adj, pct.1, pct.2)
  csv_path <- file.path(output_dir, paste0("DEG_Cluster", cluster_id, "_Peptidoglycan_vs_Control.csv"))
  write.csv(deg[, c("gene", "avg_log2FC", "p_val_adj", "pct.1", "pct.2")], file = csv_path, row.names = FALSE)
  
  # Define color scheme
  volcano_colors <- c("Up" = "#E41A1C", "Down" = "#377EB8", "NS" = "gray80")
  
  # Create volcano plot
  p <- ggplot(deg, aes(x = log2FC, y = -log10(p_val_adj), color = sig)) +
    geom_point(alpha = 0.8, size = 1.5) +
    geom_text_repel(
      data = filter(deg, sig != "NS"),
      aes(label = gene),
      size = 3,
      max.overlaps = 100
    ) +
    scale_color_manual(values = volcano_colors) +
    geom_vline(xintercept = c(-log(1.1), log(1.1)), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
    labs(
      title = paste0("Volcano Plot: Peptidoglycan vs Control (Cluster ", cluster_id, ")"),
      x = "log2 Fold Change",
      y = "-log10 Adjusted p-value"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.title = element_blank())
  
  # Save to PDF
  pdf_path <- file.path(output_dir, paste0("Volcano_Cluster", cluster_id, "_Peptidoglycan_vs_Control.pdf"))
  ggsave(pdf_path, plot = p, width = 6, height = 5)
  
  return(p)
}

# Run for each cluster
clusters_to_plot <- c("1", "2", "4", "5", "6", "7", "8")
plots <- lapply(clusters_to_plot, function(clust) {
  make_volcano_plot(Myeloid.clean, clust, output_dir = "/Users/xinchen/Library/CloudStorage/Dropbox/lyme_disease/Pepti_Spleen_singlecell/Result/")
})


#######
# Load required libraries
library(dplyr)
library(readr)
library(stringr)

# Define the function
search_gene_across_clusters <- function(gene_list, deg_dir) {
  # Get list of all DEG CSV files in the directory
  deg_files <- list.files(deg_dir, pattern = "^DEG_Cluster.*\\.csv$", full.names = TRUE)
  
  # Initialize result list
  all_hits <- lapply(deg_files, function(file) {
    # Read the file
    df <- read_csv(file)
    
    # Extract cluster number from filename
    cluster_id <- str_extract(basename(file), "(?<=Cluster)\\d+")
    
    # Filter for genes of interest
    df_filtered <- df %>% 
      filter(gene %in% gene_list) %>%
      mutate(cluster = cluster_id)
    
    return(df_filtered)
  })
  
  # Combine all results
  result_df <- bind_rows(all_hits)
  
  return(result_df)
}

# ==== Example usage ====
# Set the directory containing the DEG_Cluster*.csv files
deg_dir <- "/Users/xinchen/Library/CloudStorage/Dropbox/lyme_disease/Pepti_Spleen_singlecell/Result/"

# Define the genes you want to search (can be one or many)
genes_of_interest <- c("CHI3L1", "IRG1", "IFNB1", "CCL7", "EDN1", "CCL8", "CXCL11", 
                       "CCL5", "IL1RN", "IL18", "CD80", "SOCS1", "IL7", "CD83", "CCL4", "CCL3")

# Search across all clusters
gene_hits <- search_gene_across_clusters(genes_of_interest, deg_dir)

# View result
print(gene_hits)

# Optional: save to CSV
write.csv(gene_hits, file = file.path(deg_dir, "GeneSearch_Results.csv"), row.names = FALSE)


######search gene from Myeloid_Cluster_Markers.csv
# Load necessary libraries
library(dplyr)
library(readr)

# Define function
search_genes_in_myeloid_markers <- function(gene_list, file_path) {
  # Read in the marker file
  markers <- read_csv(file_path)
  
  # Filter for genes of interest
  hits <- markers %>%
    filter(gene %in% gene_list)
  
  return(hits)
}

# ==== Example usage ====
# Full path to the marker file
marker_file <- "/Users/xinchen/Library/CloudStorage/Dropbox/lyme_disease/Pepti_Spleen_singlecell/Result/Myeloid_Cluster_Markers.csv"

# Define the genes you're interested in
genes_of_interest <- c("CHI3L1", "IRG1", "IFNB1", "CCL7", "EDN1", "CCL8", "CXCL11", 
                       "CCL5", "IL1RN", "IL18", "CD80", "SOCS1", "IL7", "CD83", "CCL4", "CCL3")

# Search
myeloid_hits <- search_genes_in_myeloid_markers(genes_of_interest, marker_file)

# View the result
print(myeloid_hits)

# Optional: save results to CSV
write.csv(myeloid_hits,
          file = "/Users/xinchen/Library/CloudStorage/Dropbox/lyme_disease/Pepti_Spleen_singlecell/Result/Myeloid_GeneSearch_Results.csv",
          row.names = FALSE)







