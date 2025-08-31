
#Cytokine analysis for temp and treatment 
#Author: Xin Chen, Ph.D. 
#Date: 2025-05-24
#Last update: 2025-05-24

library(stringr)
library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(corrr)
library(factoextra)
library(FactoMineR)
library(patchwork)
library(reshape2)
library(ggrepel)
library(pheatmap)
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
#packageVersion("ggplot2")

#install.packages("ggplot2")


#mac
setwd("~/Library/CloudStorage/Dropbox/lyme_disease/NULISAseq/XCE2/")
getwd()

data.df <- read.csv("./rawdata/rawdata_count.csv", header = T, row.names = 1) %>% t() %>% data.frame()
data.df$SampleID <- rownames(data.df)
#orgnize metadata
metadata <- read.csv("./rawdata/metadata for Nulisa-seq.csv", header = T)

#metadata_temp <- filter(metadata_separated, separated_4 != "experiments2") %>% filter(separated_4 != "controls")
metadata_Media <- filter(metadata, Condition == "Meida")
metadata_Control <- filter(metadata, Condition == "Control")
metadata_treat <- filter(metadata, Condition != "Meida" & Condition != "Control")
metadata_treat

#colnames(metadata_temp) <- c("SampleID", "separated_1", "donorID", "treatment", "temperature")
colnames(metadata_treat) <- c("SampleID",  "treatment", "donorID", "Dilution")

data.df$SampleID <- sub("([A-Z])(\\d+)\\.(\\d+)", "\\1\\2-\\3", data.df$SampleID)
data.df$SampleID <- sub("\\.(DX\\d+)", "-\\1", data.df$SampleID)

data.treat <- left_join(metadata_treat, data.df, by=c("SampleID"="SampleID"))
data.treat

library(dplyr)
library(tidyr)
library(DESeq2)
library(ggplot2)
library(ggrepel)

### Filter for dilution 10
data.dil10 <- data.treat %>%
  filter(Dilution == 10)

### Create count matrix (remove metadata, transpose)
count_matrix <- data.dil10 %>%
  dplyr::select(-SampleID, -treatment, -donorID, -Dilution) %>%
  t()
colnames(count_matrix) <- data.dil10$SampleID  # Set sample IDs as colnames

### Create sample metadata
sample_metadata <- data.dil10 %>%
  dplyr::select(SampleID, treatment, donorID) %>%
  distinct() %>%
  arrange(SampleID)

### Match count matrix column order with sample metadata
count_matrix <- count_matrix[, sample_metadata$SampleID]

### Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = round(count_matrix),
  colData = sample_metadata,
  design = ~ donorID + treatment
)

### Run DESeq2
dds <- DESeq(dds)


Bb_W_res <- results(dds, contrast = c("treatment", "Bb_W", "None"))
FN_W_res <- results(dds, contrast = c("treatment", "FN_W", "None"))
SO_W_res <- results(dds, contrast = c("treatment", "SO_W", "None"))
LPS_W_res <- results(dds, contrast = c("treatment", "LPS", "None"))


library(dplyr)

# Convert to data frame and keep log2FC & padj
Bb_df <- as.data.frame(Bb_W_res) %>%
  rownames_to_column("gene") %>%
  dplyr::select(gene, log2FoldChange, padj) %>%
  rename(Bb_log2FC = log2FoldChange, Bb_padj = padj)

FN_df <- as.data.frame(FN_W_res) %>%
  rownames_to_column("gene") %>%
  dplyr::select(gene, log2FoldChange, padj) %>%
  rename(FN_log2FC = log2FoldChange, FN_padj = padj)

SO_df <- as.data.frame(SO_W_res) %>%
  rownames_to_column("gene") %>%
  dplyr::select(gene, log2FoldChange, padj) %>%
  rename(SO_log2FC = log2FoldChange, SO_padj = padj)

LPS_df <- as.data.frame(LPS_W_res) %>%
  rownames_to_column("gene") %>%
  dplyr::select(gene, log2FoldChange, padj) %>%
  rename(LPS_log2FC = log2FoldChange, LPS_padj = padj)

merged_df <- Reduce(function(x,y) full_join(x,y, by="gene"), list(LPS_df, Bb_df, FN_df, SO_df))
head(merged_df)


sig_genes <- merged_df %>%
  filter(
    (Bb_padj < 0.05) |
      (FN_padj < 0.05) |
      (SO_padj < 0.05)|
      (LPS_padj < 0.05)
  ) %>%
  dplyr::select(gene, LPS_log2FC, Bb_log2FC, FN_log2FC, SO_log2FC)

# Convert to matrix for heatmap
heatmap_mat <- sig_genes %>%
  column_to_rownames("gene") %>%
  as.matrix()


library(pheatmap)
heatmap_mat_ordered <- heatmap_mat[, c("Bb_log2FC","SO_log2FC","FN_log2FC","LPS_log2FC")]
colnames(heatmap_mat_ordered) <- c("Bb","SO","FN","LPS")
# Optional: scale by row to visualize relative changes
p1 <- pheatmap(
  heatmap_mat[, c("Bb_log2FC","SO_log2FC","FN_log2FC","LPS_log2FC")], # fixed order
  cluster_rows = TRUE,          # or FALSE if you want fixed row order
  cluster_cols = FALSE,         # no column clustering
  show_rownames = FALSE,        # hide gene names
  clustering_distance_rows = "euclidean",
  clustering_method = "complete",
  color = colorRampPalette(c("blue","white","red"))(50),
  border_color = NA,            # <--- remove the grid lines
  main = "Significant Cytokines - Log2FC"
)

ggsave("./log2DC from PG over PBS_DESeq_for paper.pdf", plot = p1, width = 7, height = 10)



library(pheatmap)
heatmap_mat_ordered <- heatmap_mat[, c("Bb_log2FC","SO_log2FC","FN_log2FC","LPS_log2FC")]
colnames(heatmap_mat_ordered) <- c("Bb","SO","FN","LPS")
# Optional: scale by row to visualize relative changes
p1 <- pheatmap(
  heatmap_mat[, c("Bb_log2FC","SO_log2FC","FN_log2FC","LPS_log2FC")], # fixed order
  cluster_rows = TRUE,          # or FALSE if you want fixed row order
  cluster_cols = FALSE,         # no column clustering
  show_rownames = TRUE,        # hide gene names
  clustering_distance_rows = "euclidean",
  clustering_method = "complete",
  color = colorRampPalette(c("blue","white","red"))(50),
  border_color = NA,            # <--- remove the grid lines
  main = "Significant Cytokines - Log2FC"
)

ggsave("./log2DC from PG over PBS_DESeq.pdf", plot = p1, width = 7, height = 30)




##########
### Run DESeq2
dds <- DESeq(dds)


Bb_D_res <- results(dds, contrast = c("treatment", "Bb_D", "None"))
FN_D_res <- results(dds, contrast = c("treatment", "FN_D", "None"))
SO_D_res <- results(dds, contrast = c("treatment", "SO_D", "None"))
LPS_D_res <- results(dds, contrast = c("treatment", "LPS", "None"))


library(dplyr)

# Convert to data frame and keep log2FC & padj
Bb_df <- as.data.frame(Bb_D_res) %>%
  rownames_to_column("gene") %>%
  dplyr::select(gene, log2FoldChange, padj) %>%
  rename(Bb_log2FC = log2FoldChange, Bb_padj = padj)

FN_df <- as.data.frame(FN_D_res) %>%
  rownames_to_column("gene") %>%
  dplyr::select(gene, log2FoldChange, padj) %>%
  rename(FN_log2FC = log2FoldChange, FN_padj = padj)

SO_df <- as.data.frame(SO_D_res) %>%
  rownames_to_column("gene") %>%
  dplyr::select(gene, log2FoldChange, padj) %>%
  rename(SO_log2FC = log2FoldChange, SO_padj = padj)

LPS_df <- as.data.frame(LPS_D_res) %>%
  rownames_to_column("gene") %>%
  dplyr::select(gene, log2FoldChange, padj) %>%
  rename(LPS_log2FC = log2FoldChange, LPS_padj = padj)

merged_df <- Reduce(function(x,y) full_join(x,y, by="gene"), list(LPS_df, Bb_df, FN_df, SO_df))
head(merged_df)


sig_genes <- merged_df %>%
  filter(
    (Bb_padj < 0.05) |
      (FN_padj < 0.05) |
      (SO_padj < 0.05)|
      (LPS_padj < 0.05)
  ) %>%
  dplyr::select(gene, LPS_log2FC, Bb_log2FC, FN_log2FC, SO_log2FC)

# Convert to matrix for heatmap
heatmap_mat <- sig_genes %>%
  column_to_rownames("gene") %>%
  as.matrix()


library(pheatmap)

# Optional: scale by row to visualize relative changes
p1=pheatmap(heatmap_mat,
            #scale = "row",
            clustering_distance_rows = "euclidean",
            clustering_distance_cols = "euclidean",
            clustering_method = "complete",
            color = colorRampPalette(c("blue","white","red"))(50),
            main = "Significant Cytokines - Log2FC")


ggsave("./log2DC from PG over PBS_DESeq_D.pdf", plot = p1, width = 7, height = 30)

ggsave("./log2DC from PG over PBS_DESeq_for paper_D.pdf", plot = p1, width = 7, height = 10)


#####

target_genes <- c(
  "IL17A.IL17F","IFNL1",
  "CCL24", "IL12p70", "IL18R1","CCL23","IFNG","CHI3L1", "SPP1",
  "IL22", "CCL2", "HLA.DRA", "IL17A","CCL22", "CCL17", "TNFRSF18", "TNFSF13"
  )

# Reorder your columns for desired order
heatmap_mat <- heatmap_mat[, c("Bb_log2FC","SO_log2FC","FN_log2FC","LPS_log2FC")]

# Keep only target genes in the same order as the list
genes_in_matrix <- target_genes[target_genes %in% rownames(heatmap_mat)]
heatmap_subset <- heatmap_mat[genes_in_matrix, ]

p2=pheatmap(heatmap_subset,
         cluster_rows = FALSE,   # keep the same order as target_genes
         cluster_cols = FALSE,   # keep column order Bb → SO → FN → LPS
         show_rownames = TRUE,   # change to FALSE if you want to hide labels
         color = colorRampPalette(c("blue","white","red"))(50),
         main = "Target Cytokines: Log2FC vs PBS"
)
p2


ggsave("./log2DC from PG over PBS_DESeq_distinct for paper.pdf", plot = p2, width = 3, height = 5)

###########
target_genes <- c(
  "CCL22", "IL17A", "HLA.DRA", "CCL2", "IL22", "SPP1",
  "CHI3L1", "IL18R1", "CCL24", "IFNG", "CCL23",
  "IL12p70",  "IFNL1", "IL17A.IL17F", "IL12p70", "CHI3L1", "SPP1",
  "VEGFA", "IL13", "CCL17", "TNFRSF18", "TNFSF13"
)