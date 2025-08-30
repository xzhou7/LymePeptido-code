#Single cell Lyme disease post merge
#Author: Xin Chen, Ph.D.
#Date Created: 052424
#trajactory
#https://personal.sron.nl/~pault/#sec:qualitative

#set up working directory
setwd("C:/Users/zhoux/Box/Xin.Chen.Shareable/R7_NR")

#setwd(dir = "~/Library/CloudStorage/Box-Box/Xin.Chen.Shareable/R7_NR")

getwd()
#load necessary package
library(data.table)
library(dplyr)
library(Seurat)
library(SeuratWrappers)
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
library(monocle3)

source("./Code/00_colorKey.R")
load(file = "./Data/Step4_022823_Myeloid.cluster_clean.RData")
pbmc.markers <- read.csv("./Results/Myeloid.clean_gene_markers.csv", header = T)

Myeloid.clean

p1 <- DimPlot(Myeloid.clean, group.by = "seurat_clusters", label = T, split.by = "condition")
p1

filter(pbmc.markers, cluster == "5")

FeaturePlot(Myeloid.clean, "CD28", min.cutoff = 0, order = T)

Myeloid.clean$cluster <- "S100A9_Mono"
Myeloid.clean$cluster[Myeloid.clean$integrated_snn_res.0.3 == 1] <- "HLA_Mono"
Myeloid.clean$cluster[Myeloid.clean$integrated_snn_res.0.3 == 2] <- "HyperInflam_Mono"
Myeloid.clean$cluster[Myeloid.clean$integrated_snn_res.0.3 == 3] <- "NAMPT_Mono"
Myeloid.clean$cluster[Myeloid.clean$integrated_snn_res.0.3 == 4] <- "CD16_Mono"
Myeloid.clean$cluster[Myeloid.clean$integrated_snn_res.0.3 == 5] <- "X5_Mono"
Myeloid.clean$cluster[Myeloid.clean$integrated_snn_res.0.3 == 6] <- "X6_Mono"
Myeloid.clean$cluster[Myeloid.clean$integrated_snn_res.0.3 == 7] <- "cDC"
Myeloid.clean$cluster[Myeloid.clean$integrated_snn_res.0.3 == 8] <- "pDC"

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
p_umap2

#ggsave(filename = "./Paper_Figures/Myeloid/Myeloid_annotated.pdf",p_umap2, width = 5, height = 4, dpi = 300)

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

count2.2$conditon <- factor(count2.2$conditon, levels = c("HD", "PTLDN", "RTH", "PTLD"))

comparasions <- list(c("HD", "PTLDN"), c("RTH","PTLD"))

p2.2 <- ggplot(count2.2, aes(x=conditon, y=Percent, color=conditon)) + geom_boxplot(alpha=0.7, outlier.alpha = 0)
p2.2 <- p2.2 + facet_wrap(.~Var2, scales = "free") + theme_cowplot()
p2.2 <- p2.2 + geom_jitter(color="black", size=0.4, alpha=0.9)
p2.2 <- p2.2 + stat_compare_means(method="t.test", comparisons = comparasions) + scale_color_manual(values=condition_color) + XZ_flip_x_label()
p2.2
#ggsave("../lyme_disease/Manuscript/Figures/Figure3/bypercent.pdf", p2.2, width = 8, height = 8, dpi = 300)

# Setting a seed for reproducibility
set.seed(123) 
Myeloid.clean$condition_level <- factor(Myeloid.clean$condition,levels = c("HD", "PTLDN", "RTH","PTLD"))

sampled_cells <- list()
for (i in levels(Myeloid.clean$condition_level)) {
  print(i)
  # Identify cells in the current condition
  cells_in_condition <- WhichCells(Myeloid.clean, expression = condition_level == i)
  sampled_cells[[i]] <- sample(cells_in_condition, size = 10000, replace = FALSE)
}

all_sampled_cells <- unlist(sampled_cells)
Myeloid.clean.sampled10000 <- subset(Myeloid.clean, cells = all_sampled_cells)
table(Myeloid.clean.sampled10000$condition_level)

p2.5 <- DimPlot(Myeloid.clean.sampled10000, group.by = "cluster", split.by = "condition_level") + scale_color_manual(values = Mono_color) + coord_fixed(ratio=1)
p2.5

#ggsave(filename = "./Paper_Figures/Myeloid/Myeloid_annotated_bycondition.pdf", p2.5, width = 12, height = 8, dpi = 300)


#https://ucdavis-bioinformatics-training.github.io/2021-August-Advanced-Topics-in-Single-Cell-RNA-Seq-Trajectory-and-Velocity/data_analysis/monocle_fixed
tol_high_contrast_palette <- c("#DDAA33", "#BB5566", "#004488")
tol_vibrant_palette <- c("#0077BB", "#33BBEE", "#009988",
                         "#EE7733", "#CC3311", "#EE3377",
                         "#BBBBBB")
tol_muted_palette <- c("#332288", "#88CCEE", "#44AA99",
                       "#117733", "#999933", "#DDCC77",
                       "#CC6677", "#882255", "#AA4499")

DefaultAssay(Myeloid.clean.sampled10000) <- "RNA"
cds <- as.cell_data_set(Myeloid.clean.sampled10000)

preprocess_cds(
  cds,
  method = c("PCA"),
  num_dim = 50)

cds <- cluster_cells(cds, reduction_method = "UMAP")

p1 <- plot_cells(cds, color_cells_by = "cluster", show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
wrap_plots(p1, p2)

integrated.sub <- subset(as.Seurat(cds, assay = NULL), monocle3_partitions == 1)
cds <- as.cell_data_set(integrated.sub)

#this step takes about 10 minutes
cds <- learn_graph(cds, use_partition = TRUE, verbose = F)

plot_cells(cds,
           color_cells_by = "cluster",
           label_cell_groups = TRUE,
           label_groups_by_cluster=T,
           label_leaves=F,
           label_branch_points=F)

cds <- order_cells(cds, root_cells = colnames(cds[,clusters(cds) == 5]))

p <- plot_cells(cds,
           color_cells_by = "pseudotime",
           group_cells_by = "cluster",
           label_cell_groups = FALSE,
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           trajectory_graph_color = "grey60")
p <- p +  facet_wrap(~ condition_level, nrow = 2) 
p

#ggsave(filename = "Paper_Figures/Monocle3_on_moncytes.pdf", p, height = 8, width = 10)

cds_graph_test_results <- graph_test(cds,
                                     neighbor_graph = "principal_graph",
                                     cores = 8)

rowData(cds)$gene_short_name <- row.names(rowData(cds))

head(cds_graph_test_results, error=FALSE, message=FALSE, warning=FALSE)

deg_ids <- rownames(subset(cds_graph_test_results[order(cds_graph_test_results$morans_I, decreasing = TRUE),], q_value < 0.05))

p <- plot_cells(cds,
           genes=head(deg_ids),
           show_trajectory_graph = FALSE,
           label_cell_groups = FALSE,
           label_leaves = FALSE)

my_genes <- row.names(subset(fData(cds), gene_short_name %in% c("CXCL1"))) 
cds_subset <- cds[my_genes,]

plot_genes_in_pseudotime(cds_subset) + facet_wrap(~ condition_level, nrow = 4)

gene_module_df <- find_gene_modules(cds[deg_ids,], resolution=c(10^seq(-6,-1)))





