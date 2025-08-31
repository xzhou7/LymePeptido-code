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
data.long.treat <- data.treat %>% pivot_longer(cols = AGER:WNT7A,
                                               names_to = "cytokine", values_to = "value")

table(data.long.treat$treatment)


data.long.treat$treatment <- factor(data.long.treat$treatment, levels = c("None", "LPS","Bb_D","Bb_W","FN_D","FN_W",
                                                                          "SO_D", "SO_W"))
data.long.treat$donorID <- factor(data.long.treat$donorID)


#Filter to Dilution 10
#df_dil10 <- data.long.treat %>%
 # filter(Dilution == 10)
#Summarize by treatment × cytokine
#cytokine_summary <- df_dil10 %>%
#  group_by(treatment, cytokine) %>%
#  summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop")

#Pivot to wide format for heatmap
#cytokine_wide <- cytokine_summary %>%
#  pivot_wider(names_from = cytokine, values_from = mean_value) %>%
#  column_to_rownames("treatment") %>%
#  as.data.frame()
#cytokine_matrix <- as.matrix(cytokine_wide)

#optimal: Z-score normalize across treatments
#cytokine_matrix_scaled <- scale(cytokine_matrix)

#1) remove rolw with NA
#cytokine_matrix_clean <- cytokine_matrix_scaled[complete.cases(cytokine_matrix_scaled), ]
#cytokine_matrix_clean <- cytokine_matrix_scaled[, colSums(is.na(cytokine_matrix_scaled)) == 0]


#cytokine_matrix_clean$treatment <- factor(cytokine_matrix_clean$treatment, levels = c("None", "LPS","Bb_D","Bb_W","FN_D","FN_W",
                                                                          "SO_D", "SO_W", "SP_W","SP_D","VP_W","VP_D"))
#Heatmap to compare cytokine profiles across treatments
#pheatmap(cytokine_matrix_clean,
#         cluster_rows = TRUE,
#         cluster_cols = TRUE,
#         main = "Cytokine Induction Across Treatments (Dilution 10)")

#############if select the cytokine first###########
# Subset the dataset to only include selected treatments
data.filtered <- subset(data.long.treat, treatment %in% c("None", "LPS","Bb_D","Bb_W","FN_D","FN_W","SO_D", "SO_W"))

# Double-check the result
table(data.filtered$treatment)

df_dil10 <- data.filtered %>%
  filter(Dilution == 10)
#Summarize by treatment × cytokine
cytokine_summary <- df_dil10 %>%
  group_by(treatment, cytokine) %>%
  summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop")

#Pivot to wide format for heatmap
cytokine_wide <- cytokine_summary %>%
  pivot_wider(names_from = cytokine, values_from = mean_value) %>%
  column_to_rownames("treatment") %>%
  as.data.frame()
cytokine_matrix <- as.matrix(cytokine_wide)

#optimal: Z-score normalize across treatments
cytokine_matrix_scaled <- scale(cytokine_matrix)

#1) remove rolw with NA
cytokine_matrix_clean <- cytokine_matrix_scaled[complete.cases(cytokine_matrix_scaled), ]
cytokine_matrix_clean <- cytokine_matrix_scaled[, colSums(is.na(cytokine_matrix_scaled)) == 0]

# Step 1: Identify cytokines where the 'None' row is > 0
cytokines_to_remove <- colnames(cytokine_matrix_clean)[cytokine_matrix_clean["None", ] > 0]

# Step 2: Drop those cytokines (columns) from the matrix
cytokine_matrix_filtered <- cytokine_matrix_clean[, !(colnames(cytokine_matrix_clean) %in% cytokines_to_remove)]

# Step 3 (Optional): Check remaining cytokines
colnames(cytokine_matrix_filtered)

#identify cytokine greater than 1 in Bb_W
Bb_W_cytokines <- colnames(cytokine_matrix_clean)[cytokine_matrix_clean["Bb_W", ] >1]
#identify cytokine greater than 1 in Bb_D
Bb_D_cytokines <- colnames(cytokine_matrix_clean)[cytokine_matrix_clean["Bb_D", ] >1]
#identify cytokine greater than 1 in LPS
LPS_cytokines <- colnames(cytokine_matrix_clean)[cytokine_matrix_clean["LPS", ] >1]

#identify cytokine greater than 1 in SO_W
SOW_cytokines <- colnames(cytokine_matrix_clean)[cytokine_matrix_clean["LPS", ] >1]

#library(clusterProfiler)
#library(org.Hs.eg.db)

# Convert Bb_W cytokines
bb_w_entrez <- bitr(Bb_W_cytokines, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

# Convert Bb_D cytokines
bb_d_entrez <- bitr(Bb_D_cytokines, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

# Convert SO_D cytokines
so_w_entrez <- bitr(SOW_cytokines, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

# GO Biological Process
go_bb_w <- enrichGO(gene = bb_w_entrez$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", readable = TRUE)
go_bb_d <- enrichGO(gene = bb_d_entrez$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", readable = TRUE)
go_so_w <- enrichGO(gene = bb_d_entrez$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", readable = TRUE)
# KEGG Pathway (optional)
kegg_bb_w <- enrichKEGG(gene = bb_w_entrez$ENTREZID, organism = "hsa")
kegg_bb_d <- enrichKEGG(gene = bb_d_entrez$ENTREZID, organism = "hsa")


dotplot(go_bb_w, showCategory = 20, title = "GO Enrichment for Bb_W")
dotplot(go_bb_d, showCategory = 20, title = "GO Enrichment for Bb_D")
dotplot(go_so_w, showCategory = 20, title = "GO Enrichment for Bb_D")


pheatmap(cytokine_matrix_filtered,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Autoimmune-Related Monocyte Cytokines (Dilution 10)")


#highlight autoimmune cytokines:
autoimmune_monocyte_cytokines <- c("IL1B", "IL6", "TNF", "IL12B", "IL23","CXCL1", "CXCL2", "CXCL3", "CXCL4", "CXCL5", "IL19", "CCL1", "IL17A.IL17F", "CCL24", 
                                   "CXCL10", "CCL2", "CCL3", "CCL4", "TNFSF13B", "CCL7", "CCL8", "IFNB1", "TNFSF18", "IFNA1", "IFNA3", "TNFSF10", "CXCL11",
                                   "IL10", "CCL24", "TREM2", "CXCL16", "IL23", "MMP1", "CHI3L1", "CSF1", "IL12RB1", "CCL22", "MMP9", "NAMPT", "CSF2", "CSF3", "IL22", 
                                   "IL18BP")

# Only keep those cytokines that exist in your matrix
available_cytokines <- intersect(autoimmune_monocyte_cytokines, colnames(cytokine_matrix_clean))

cytokine_autoimmune_subset <- cytokine_matrix_clean[, available_cytokines]

#highly in a heatmap
pheatmap(cytokine_autoimmune_subset,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Autoimmune-Related Monocyte Cytokines (Dilution 10)")

################
#Top 10 Cytokines per Treatment
# Convert to data frame for easier handling
scaled_df <- as.data.frame(cytokine_matrix_scaled)

# Add treatment as a column
scaled_df$treatment <- rownames(scaled_df)

# Reshape to long format
library(tidyr)
long_scaled <- scaled_df %>%
  pivot_longer(-treatment, names_to = "cytokine", values_to = "z_score")

# For each treatment, get top 10 cytokines by z-score
top10_cytokines <- long_scaled %>%
  group_by(treatment) %>%
  arrange(desc(z_score)) %>%
  slice_head(n = 10) %>%
  ungroup()
print(top10_cytokines)
#write.csv(top10_cytokines, file = "./result/top10_cytokines_per_treatment.csv", row.names = FALSE)


#Plot Heatmap of Top 10 Cytokines per Treatment
top_cytokines <- unique(top10_cytokines$cytokine)
cytokine_top_matrix <- cytokine_matrix_scaled[, top_cytokines]
library(pheatmap)

pheatmap(cytokine_top_matrix,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "none",  # already z-scored
         main = "Top 10 Cytokines per Treatment",
         fontsize_row = 8,
         fontsize_col = 10)


top_cytokines <- unique(top10_cytokines$cytokine)
cytokine_top_matrix <- cytokine_matrix[, top_cytokines]
cytokine_matrix_scaled <- scale(cytokine_top_matrix)
library(pheatmap)

pheatmap(cytokine_matrix_scaled,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "none",  # already z-scored
         main = "Top 10 Cytokines per Treatment",
         fontsize_row = 8,
         fontsize_col = 10)

#################
# Filter to Dilution 10 and selected treatments
df_dil10 <- data.long.treat %>%
  filter(Dilution == 10, treatment %in% c("None", "Bb_D", "Bb_W", "LPS"))

# Summarize by treatment × cytokine
cytokine_summary <- df_dil10 %>%
  group_by(treatment, cytokine) %>%
  summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop")

# Pivot to wide format for heatmap
cytokine_wide <- cytokine_summary %>%
  pivot_wider(names_from = cytokine, values_from = mean_value) %>%
  column_to_rownames("treatment") %>%
  as.data.frame()
# Convert to matrix
cytokine_matrix <- as.matrix(cytokine_wide)

cytokine_matrix_clean <- cytokine_matrix_scaled[complete.cases(cytokine_matrix_scaled), ]
cytokine_matrix_clean <- cytokine_matrix_scaled[, colSums(is.na(cytokine_matrix_scaled)) == 0]

#Heatmap to compare cytokine profiles across treatments
pheatmap(cytokine_matrix_clean,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Cytokine Induction Across Treatments (Dilution 10)")

#Plot Heatmap of Top 10 Cytokines per Treatment
top_cytokines <- unique(top10_cytokines$cytokine)
cytokine_top_matrix <- cytokine_matrix_scaled[, top_cytokines]
library(pheatmap)

pheatmap(cytokine_top_matrix,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "none",  # already z-scored
         main = "Top 10 Cytokines per Treatment",
         fontsize_row = 8,
         fontsize_col = 10)


top_cytokines <- unique(top10_cytokines$cytokine)
cytokine_top_matrix <- cytokine_matrix[, top_cytokines]
cytokine_matrix_scaled <- scale(cytokine_top_matrix)
library(pheatmap)

pheatmap(cytokine_matrix_scaled,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "none",  # already z-scored
         main = "Top 10 Cytokines per Treatment",
         fontsize_row = 8,
         fontsize_col = 10)


############DESeq2 for peptidoglycan#################XC
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

### Get DE results for Bb_W vs Bb_D
res <- results(dds, contrast = c("treatment", "Bb_W", "SO_W"))

### Convert to data frame and add gene names
res_df <- as.data.frame(res)
res_df$cytokine <- rownames(res_df)

### Identify significant genes
res_df <- res_df %>%
  mutate(
    sig = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, "Significant", "Not Significant"),
    label = ifelse(sig == "Significant", cytokine, NA)
  )

### OPTIONAL: Add custom color groups (can adjust logic here)
res_df <- res_df %>%
  mutate(color_group = ifelse(sig == "Significant", "red", "black"))

### Label logic for volcano plot
res_df <- res_df %>%
  mutate(display_label = ifelse(color_group == "red" | log2FoldChange > 5 | log2FoldChange < -5, label, NA))

### Volcano Plot
p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = color_group)) +
  geom_point(alpha = 0.8) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(aes(label = display_label), size = 4, max.overlaps = 30) +
  labs(
    title = "Volcano Plot: DESeq2 (Bb_W vs Bb_D)",
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted p-value"
  ) +
  theme_minimal() +
  scale_color_manual(values = c("black" = "black", "red" = "red")) +
  ylim(0, 38)
p
### Save figure
#ggsave(filename = "../../Manuscript/Figures/Figure 6-other peptidoglycans/vocano_SPW_VPW.pdf",
#       plot = p, width = 5, height = 5, dpi = 300)


###########################

####monocyte_patient DEG gene#######
PTLD_mono <- read.csv(file = "~/Library/CloudStorage/Dropbox/lyme_disease/PTLD_HD_Monocytes/HD_PTLD.csv", header = TRUE, stringsAsFactors = FALSE)

sig_genes <- PTLD_mono %>%
  dplyr::filter(padj < 0.0005 & abs(log2FoldChange) > 1)

sig_gene_list <- sig_genes$X

###sign_gene from peptidoglycan Bb_W
BbW_sig_genes <- res %>%
  as.data.frame() %>%
  dplyr::filter(padj < 0.0005 & abs(log2FoldChange) > 1)


BbW_sig_genes_list <-rownames(BbW_sig_genes)
overlap_genes <- intersect(sig_gene_list, BbW_sig_genes_list)


library(dplyr)
library(ggrepel)
# Add color group: red for overlap genes, black otherwise
res_df <- res_df %>%
  mutate(
    color_group = ifelse(cytokine %in% overlap_genes, "red", "black"),
    display_label = ifelse(cytokine %in% overlap_genes | log2FoldChange > 5 | log2FoldChange < -5, cytokine, NA)
  )

# Create a new column for labels based on your conditions
res_df <- res_df %>%
  mutate(display_label = ifelse(color_group == "red" | (color_group == "black" & log2FoldChange < -5)| (color_group == "black" & log2FoldChange > 5), label, NA))

# Plot with updated label column
p11 = ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = color_group)) +
  geom_point(alpha = 0.8) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_text_repel(aes(label = display_label), size = 5, max.overlaps = 30) +
  labs(title = "Volcano Plot: DESeq2 (Bb_W vs Bb_D)",
       x = "Log2 Fold Change", y = "-Log10 Adjusted p-value") +
  theme_minimal() +
  scale_color_manual(values = c("black" = "black", "red" = "red")) +
  ylim(0, 40)
p11
ggsave(filename = "../../Manuscript/Figures/Figure 6-other peptidoglycans/vocano_BbW_BbD.pdf",p11, width = 5, height = 6, dpi = 300)

########XC


######TRY TO VITUALIZE THE overlap DEG cytokines from each stimulus compared to control using Upset and others
library(dplyr)

# Then separate the steps clearly
counts <- data.treat %>%
  filter(Dilution == 10) %>%
  dplyr::select(-SampleID, -treatment, -donorID, -Dilution)
# Now use as.matrix safely
count_matrix <- as.matrix(data.dil10)
sample_info <- data.treat %>%
  filter(Dilution == 10) %>%
  dplyr::select(SampleID, treatment, donorID, Dilution)

########
library(dplyr)

data.dil10 <- data.treat %>%
  filter(Dilution == 10)
library(tidyr)

# Remove metadata columns
count_matrix <- data.dil10 %>%
  dplyr::select(-SampleID, -treatment, -donorID, -Dilution) %>%
  t()

# Set sample IDs as column names
colnames(count_matrix) <- data.dil10$SampleID

sample_metadata <- data.dil10 %>%
  dplyr::select(SampleID, treatment, donorID) %>%
  distinct() %>%
  arrange(SampleID)

# Match order
count_matrix <- count_matrix[, sample_metadata$SampleID]

dds <- DESeqDataSetFromMatrix(
  countData = round(count_matrix),
  colData = sample_metadata,
  design = ~ donorID + treatment
)
dds <- DESeq(dds)
res_Bb_W <- results(dds, contrast = c("treatment", "Bb_W", "None"))
res_Bb_D <- results(dds, contrast = c("treatment", "Bb_D", "None"))
res_FN_W <- results(dds, contrast = c("treatment", "FN_W", "None"))
res_FN_D <- results(dds, contrast = c("treatment", "FN_D", "None"))
res_SO_W <- results(dds, contrast = c("treatment", "SO_W", "None"))
res_SO_D <- results(dds, contrast = c("treatment", "SO_D", "None"))#"None", "LPS","Bb_D","Bb_W","FN_D","FN_W",
res_SP_W <- results(dds, contrast = c("treatment", "SP_W", "None"))#"SO_D", "SO_W", "SP_W","SP_D","VP_W","VP_D"
res_SP_D <- results(dds, contrast = c("treatment", "SP_D", "None"))
res_VP_W <- results(dds, contrast = c("treatment", "VP_W", "None"))
res_VP_D <- results(dds, contrast = c("treatment", "VP_D", "None"))

library(dplyr)

deg_Bb_W <- as.data.frame(res_Bb_W) %>%
  rownames_to_column("cytokine") %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1) %>%
  pull(cytokine)

deg_Bb_D <- as.data.frame(res_Bb_D) %>%
  rownames_to_column("cytokine") %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1) %>%
  pull(cytokine)

deg_SO_W <- as.data.frame(res_SO_W) %>%
  rownames_to_column("cytokine") %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1) %>%
  pull(cytokine)

deg_SO_D <- as.data.frame(res_SO_D) %>%
  rownames_to_column("cytokine") %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1) %>%
  pull(cytokine)

deg_FN_W <- as.data.frame(res_FN_W) %>%
  rownames_to_column("cytokine") %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1) %>%
  pull(cytokine)

deg_FN_D <- as.data.frame(res_FN_D) %>%
  rownames_to_column("cytokine") %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1) %>%
  pull(cytokine)


deg_VP_W <- as.data.frame(res_VP_W) %>%
  rownames_to_column("cytokine") %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1) %>%
  pull(cytokine)

deg_VP_D <- as.data.frame(res_VP_D) %>%
  rownames_to_column("cytokine") %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1) %>%
  pull(cytokine)

deg_SP_W <- as.data.frame(res_SP_W) %>%
  rownames_to_column("cytokine") %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1) %>%
  pull(cytokine)

deg_SP_D <- as.data.frame(res_SP_D) %>%
  rownames_to_column("cytokine") %>%
  filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1) %>%
  pull(cytokine)

deg_list <- list(
  Bb_W = deg_Bb_W,
  Bb_D = deg_Bb_D,
  SO_W = deg_SO_W,
  SO_D = deg_SO_D,
  FN_W = deg_FN_W,
  FN_D = deg_FN_D,
  VP_W = deg_VP_W,
  VP_D = deg_VP_D,
  SP_W = deg_SP_W,
  SP_D = deg_SP_D
)

library(UpSetR)

# Convert list to binary matrix
deg_sets_matrix <- fromList(deg_list)


# Convert to UpSet input
deg_sets <- fromList(deg_list)
upset(deg_sets, order.by = "freq", main.bar.color = "steelblue")

upset(deg_sets_matrix,
      sets = c("Bb_W", "SO_W", "FN_W", "VP_W", "SP_W"),
      order.by = "freq",
      main.bar.color = "steelblue",
      sets.bar.color = "gray40")

upset(deg_sets_matrix,
      sets = c("SO_W", "SO_D"),
      order.by = "freq",
      main.bar.color = "steelblue",
      sets.bar.color = "gray40")

upset(deg_sets_matrix,
      sets = c("FN_W", "FN_D"),
      order.by = "freq",
      main.bar.color = "steelblue",
      sets.bar.color = "gray40")

upset(deg_sets_matrix,
      sets = c("SP_W", "SP_D"),
      order.by = "freq",
      main.bar.color = "steelblue",
      sets.bar.color = "gray40")

upset(deg_sets_matrix,
      sets = c("Bb_W", "Bb_D"),
      order.by = "freq",
      main.bar.color = "steelblue",
      sets.bar.color = "gray40")


upset(deg_sets_matrix,
      sets = c("VP_W", "VP_D"),
      order.by = "freq",
      main.bar.color = "steelblue",
      sets.bar.color = "gray40")

#####Venn dia
library(VennDiagram)

venn.plot <- venn.diagram(
  x = list(
    Bb_W = deg_Bb_W,
    Bb_D = deg_Bb_D
  ),
  filename = NULL,
  fill = c("red", "blue"),
  alpha = 0.5,
  cex = 1.2,
  cat.cex = 1.2,
  cat.pos = 0
)

grid::grid.draw(venn.plot)

venn.plot <- venn.diagram(
  x = list(
    VP_W = deg_VP_W,
    VP_D = deg_VP_D
  ),
  filename = NULL,
  fill = c("red", "blue"),
  alpha = 0.5,
  cex = 1.2,
  cat.cex = 1.2,
  cat.pos = 0
)

#try heatmap
deg_list <- list(
  Bb_W = deg_Bb_W,
  Bb_D = deg_Bb_D,
  SO_W = deg_SO_W,
  SO_D = deg_SO_D,
  FN_W = deg_FN_W,
  FN_D = deg_FN_D,
  VP_W = deg_VP_W,
  VP_D = deg_VP_D,
  SP_W = deg_SP_W,
  SP_D = deg_SP_D
)

# Get all unique cytokines across all DEG lists
all_cytokines <- unique(unlist(deg_list))

# Build a binary matrix
deg_matrix <- sapply(deg_list, function(set) {
  as.integer(all_cytokines %in% set)
})

rownames(deg_matrix) <- all_cytokines
library(pheatmap)

pheatmap(
  deg_matrix,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = c("white", "black"),
  display_numbers = TRUE,
  main = "DEG Cytokine Overlap Across Conditions",
  fontsize_row = 8,
  fontsize_col = 10
)

###log2FC for heatmap
library(dplyr)

# Function to extract log2FC from a DESeq2 result
extract_log2fc <- function(res_object, condition_name) {
  as.data.frame(res_object) %>%
    rownames_to_column("cytokine") %>%
    filter(!is.na(padj), padj < 0.05, abs(log2FoldChange) > 1) %>%
    dplyr::select(cytokine, log2FoldChange) %>%
    rename(!!condition_name := log2FoldChange)
}

fc_list <- list(
  Bb_W = extract_log2fc(res_Bb_W, "Bb_W"),
  Bb_D = extract_log2fc(res_Bb_D, "Bb_D"),
  SO_W = extract_log2fc(res_SO_W, "SO_W"),
  SO_D = extract_log2fc(res_SO_D, "SO_D"),
  FN_W = extract_log2fc(res_FN_W, "FN_W"),
  FN_D = extract_log2fc(res_FN_D, "FN_D"),
  VP_W = extract_log2fc(res_VP_W, "VP_W"),
  VP_D = extract_log2fc(res_VP_D, "VP_D"),
  SP_W = extract_log2fc(res_SP_W, "SP_W"),
  SP_D = extract_log2fc(res_SP_D, "SP_D")
)

# Reduce to one wide matrix
library(purrr)

log2fc_matrix <- purrr::reduce(fc_list, full_join, by = "cytokine") %>%
  column_to_rownames("cytokine")

# 1️⃣ Fill NAs with 0
log2fc_matrix[is.na(log2fc_matrix)] <- 0

# 2️⃣ Filter cytokines present in >=5 conditions (as you originally did)
log2fc_filtered <- log2fc_matrix[rowSums(log2fc_matrix != 0) >= 5, ]

# 3️⃣ Cytokines you want to FORCE INCLUDE
extra_cytokines <- c("IL1B", "IL6", "TNF", "IL12B", "IL23","CXCL1", "CXCL2", "CXCL3", "CXCL4", "CXCL5", "IL19", "CCL1", "IL17A.IL17F", "CCL24", 
                     "CXCL10", "CCL2", "CCL3", "CCL4", "TNFSF13B", "CCL7", "CCL8", "IFNB1", "TNFSF18", "IFNA1", "IFNA3", "TNFSF10", "CXCL11",
                     "IL10", "CCL24", "TREM2", "CXCL16", "IL23", "MMP1", "CHI3L1", "CSF1", "IL12RB1", "CCL22", "MMP9", "NAMPT", "CSF2", "CSF3", "IL22")  # <-- add any others

# 4️⃣ Merge existing and extra cytokines
existing_cytokines <- rownames(log2fc_filtered)
final_cytokines <- union(existing_cytokines, extra_cytokines)

# 5️⃣ Subset from full matrix to include both DEGs and forced-in cytokines
log2fc_final <- log2fc_matrix[rownames(log2fc_matrix) %in% final_cytokines, ]

# 6️⃣ Cytokines you want to REMOVE
remove_cytokines <- c("FGF21", "LILRB2", "SIRPA", "LGALS9", "MIF", "CD4", "S100A12", "IL2", "IL18BP", "VSTM1", "BST2", "IL7R", "OSM", "CXCL13", "CCL19", "VEGFA", "CSF3R")  # <-- add any you want to exclude

# 7️⃣ Apply exclusion
log2fc_final <- log2fc_final[!rownames(log2fc_final) %in% remove_cytokines, ]

# 8️⃣ Sort rows alphabetically (optional)
log2fc_final <- log2fc_final[order(rownames(log2fc_final)), ]

# 9️⃣ Apply your desired column order
desired_order <- c("Bb_D", "Bb_W", "FN_W", "SP_W", "VP_W", "SO_W", "FN_D", "SP_D", "VP_D", "SO_D")
desired_order <- desired_order[desired_order %in% colnames(log2fc_final)]  # make sure order matches columns
log2fc_final <- log2fc_final[, desired_order]

# 🔟 Finally: plot heatmap
pheatmap(
  log2fc_final,
  cluster_rows = TRUE,
  cluster_cols = FALSE,  # fix column order
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "Log2 Fold Change of Selected Cytokines",
  fontsize_row = 8,
  fontsize_col = 10
)

#ggsave(filename = "../../Manuscript/Figures/Figure 6-other peptidoglycans/heatmap_treatment.pdf",p2, width = 4, height = 12, dpi = 300)


#####
library(dplyr)
library(ggplot2)
library(ggrepel)
library(DESeq2)
library(readr)

# Treatments list
treatments <- c("LPS", "Bb_D", "Bb_W", "FN_D", "SP_D", "VP_D", "SO_D")

# Load monocyte DEG
PTLD_mono <- read.csv("~/Library/CloudStorage/Dropbox/lyme_disease/PTLD_HD_Monocytes/HD_PTLD.csv", header = TRUE)
sig_genes <- PTLD_mono %>% filter(padj < 0.0005 & abs(log2FoldChange) > 6)
sig_gene_list <- sig_genes$X

# Loop through all treatments
for (treat in treatments) {
  
  # Get DE result for each treatment
  res <- results(dds, contrast = c("treatment", treat, "None"))
  res_df <- as.data.frame(res)
  res_df$cytokine <- rownames(res_df)
  
  # Find treatment-specific overlap genes
  sig_genes_treat <- res_df %>% filter(padj < 0.0005 & abs(log2FoldChange) > 1)
  sig_genes_treat_list <- rownames(sig_genes_treat)
  overlap_genes <- intersect(sig_gene_list, sig_genes_treat_list)
  
  # Label significant genes
  res_df <- res_df %>%
    mutate(
      sig = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, "Significant", "Not Significant"),
      label = ifelse(sig == "Significant", cytokine, NA),
      color_group = ifelse(cytokine %in% overlap_genes, "red", "black"),
      display_label = ifelse(cytokine %in% overlap_genes | log2FoldChange > 5 | log2FoldChange < -5, cytokine, NA),
      display_label = ifelse(color_group == "red" | 
                               (color_group == "black" & log2FoldChange < -5) | 
                               (color_group == "black" & log2FoldChange > 5),
                             label, NA)
    )
  
  # Volcano plot for this treatment
  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = color_group)) +
    geom_point(alpha = 0.8) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_text_repel(aes(label = display_label), size = 4, max.overlaps = 30) +
    labs(title = paste("Volcano Plot:", treat, "vs None"),
         x = "Log2 Fold Change", y = "-Log10 Adjusted p-value") +
    theme_minimal() +
    scale_color_manual(values = c("black" = "black", "red" = "red")) +
    ylim(0, 38)
  
  # Show the plot
  print(p)
  
  # Save each plot if you want:
  ggsave(filename = paste0("./Volcano_", treat, ".pdf"), plot = p, width = 5, height = 5, dpi = 300)
}




#####################################

library(dplyr)
library(ggplot2)
library(ggrepel)
library(DESeq2)
library(patchwork)
library(readr)

# ========================
# 1. Load monocyte PTLD DEGs (highly significant)
# ========================

PTLD_mono <- read.csv("~/Library/CloudStorage/Dropbox/lyme_disease/PTLD_HD_Monocytes/HD_PTLD.csv", header = TRUE)

# Highly significant PTLD genes (very strict)
monocyte_threshold_padj <- 0.0005
monocyte_threshold_log2FC <- 6

sig_genes <- PTLD_mono %>%
  filter(padj < monocyte_threshold_padj & abs(log2FoldChange) > monocyte_threshold_log2FC)

sig_gene_list <- sig_genes$X

# ========================
# 2. Set treatment groups
# ========================

# For peptidoglycan W group
treatments_W <- c("Bb_W", "FN_W", "SP_W", "VP_W", "SO_W")

# For peptidoglycan D group
treatments_D <- c("Bb_D", "FN_D", "SP_D", "VP_D", "SO_D")

# Color palettes
mycol_W <- c("Bb_W" = "#fc3c46", "FN_W" = "#b12d30", "SP_W" = "#43b5e6", "SO_W" = "#58ac41", "VP_W" = "#f0bbcb")
mycol_D <- c("Bb_D" = "#fc3c46", "FN_D" = "#b12d30", "SP_D" = "#43b5e6", "SO_D" = "#58ac41", "VP_D" = "#f0bbcb")

# ========================
# 3. Function for volcano data generation
# ========================

generate_deg_list <- function(treatments) {
  
  deg_list <- list()
  
  for (treat in treatments) {
    
    # DESeq2 results
    res <- results(dds, contrast = c("treatment", treat, "None"))
    res_df <- as.data.frame(res)
    res_df$cytokine <- rownames(res_df)
    res_df$Group <- treat
    
    # Find overlap genes
    sig_res <- res_df %>%
      filter(padj < 0.0005 & abs(log2FoldChange) > 1)
    
    sig_genes_treat_list <- rownames(sig_res)
    overlap_genes <- intersect(sig_gene_list, sig_genes_treat_list)
    
    # Annotate significance & color
    res_df <- res_df %>%
      mutate(
        sig = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, "Significant", "Not Significant"),
        label = ifelse(sig == "Significant", cytokine, NA),
        color_group = ifelse(cytokine %in% overlap_genes, "red", "black"),
        display_label = ifelse(cytokine %in% overlap_genes, cytokine, NA)
      )
    
    deg_list[[treat]] <- res_df
  }
  return(bind_rows(deg_list))
}

# ========================
# 4. Generate volcano data for W and D groups
# ========================

deg_all_W <- generate_deg_list(treatments_W)
deg_all_D <- generate_deg_list(treatments_D)

# ========================
# 5. Volcano plot function
# ========================

make_volcano_facet <- function(deg_all, mycol) {
  
  highlight_df <- deg_all %>% filter(!is.na(display_label))
  
  p <- ggplot() +
    geom_point(data = deg_all, aes(x = log2FoldChange, y = -log10(padj)), size = 0.8, color = "grey") +
    geom_point(data = highlight_df, aes(x = log2FoldChange, y = -log10(padj), color = Group), size = 1.5) +
    geom_text_repel(data = highlight_df, aes(x = log2FoldChange, y = -log10(padj), label = display_label),
                    size = 3.5, max.overlaps = 30) +
    geom_vline(xintercept = c(-1, 1), color = "grey50", linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), color = "grey50", linetype = "dashed") +
    coord_flip(xlim = c(-12, 12)) +  # unified x-axis
    facet_grid(. ~ Group, scales = "fixed") +
    scale_color_manual(values = mycol) +
    labs(x = "log2 Fold Change", y = "-log10 Adjusted p-value") +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      axis.text = element_text(size = 10),
      axis.text.x = element_text(angle = 45, vjust = 0.8),
      strip.text.x = element_text(size = 10, face = "bold")
    )
  return(p)
}

# ========================
# 6. Generate facet volcano plots
# ========================

p_W <- make_volcano_facet(deg_all_W, mycol_W)
p_D <- make_volcano_facet(deg_all_D, mycol_D)

# ========================
# 7. Combine volcano panels (vertically)
# ========================

library(patchwork)

combined_plot <- p_D / p_W  # vertical stacking

# Show
combined_plot

# Save figure
#ggsave("./PTLD_Volcano_Facets.pdf", plot = combined_plot, width = 8, height = 8, dpi = 300)

# ========================
# 8. Generate overlap counts for summary
# ========================

# Overlap counts:
get_overlap_count <- function(treatments, deg_all) {
  sapply(treatments, function(treat) {
    deg_all %>% filter(Group == treat, color_group == "red") %>% nrow()
  })
}

overlap_counts_W <- get_overlap_count(treatments_W, deg_all_W)
overlap_counts_D <- get_overlap_count(treatments_D, deg_all_D)

# Barplot summary
overlap_df <- data.frame(
  Treatment = c(treatments_W, treatments_D),
  Overlap = c(overlap_counts_W, overlap_counts_D)
)

ggplot(overlap_df, aes(x = Treatment, y = Overlap, fill = Treatment)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(title = "Overlap with PTLD Monocytes", y = "Overlap Genes", x = "Treatment") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c(mycol_W, mycol_D))





###################
library(dplyr)
library(ggplot2)
library(ggrepel)
library(patchwork)

# ========================
# 1. Treatment groups (use same treatments)
# ========================
treatments_W <- c("Bb_W", "FN_W", "SP_W", "VP_W", "SO_W")
treatments_D <- c("Bb_D", "FN_D", "SP_D", "VP_D", "SO_D")

# Color palettes (same as before)
mycol_W <- c("Bb_W" = "#b12d30", "FN_W" = "#D2691E", "SP_W" = "#43b5e6", "SO_W" = "#58ac41", "VP_W" = "#800080")
mycol_D <- c("Bb_D" = "#b12d30", "FN_D" = "#D2691E", "SP_D" = "#43b5e6", "SO_D" = "#58ac41", "VP_D" = "#800080")

# ========================
# 2. Generate DEGs for volcano (no overlap, all DEGs)
# ========================
generate_all_deg <- function(treatments) {
  
  deg_list <- list()
  
  for (treat in treatments) {
    res <- results(dds, contrast = c("treatment", treat, "None"))
    res_df <- as.data.frame(res)
    res_df$cytokine <- rownames(res_df)
    res_df$Group <- treat
    
    res_df <- res_df %>%
      mutate(
        sig = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, "Significant", "Not Significant"),
        display_label = ifelse(abs(log2FoldChange) > 6, cytokine, NA)
      )
    
    deg_list[[treat]] <- res_df
  }
  
  return(bind_rows(deg_list))
}

# Generate data
deg_all_W <- generate_all_deg(treatments_W)
deg_all_D <- generate_all_deg(treatments_D)

# ========================
# 3. Volcano plot function
# ========================

make_simple_volcano <- function(deg_all, mycol) {
  highlight_df <- deg_all %>% filter(!is.na(display_label))
  
  p <- ggplot() +
    geom_point(data = deg_all, aes(x = log2FoldChange, y = -log10(padj)), size = 0.8, color = "grey") +
    geom_text_repel(data = highlight_df, aes(x = log2FoldChange, y = -log10(padj), label = display_label),
                    size = 3.5, max.overlaps = 30) +
    geom_vline(xintercept = c(-1, 1), color = "grey50", linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), color = "grey50", linetype = "dashed") +
    coord_flip(xlim = c(-12, 12)) +
    facet_grid(. ~ Group, scales = "fixed") +
    scale_color_manual(values = mycol) +
    labs(x = "log2 Fold Change", y = "-log10 Adjusted p-value") +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      axis.text = element_text(size = 10),
      axis.text.x = element_text(angle = 45, vjust = 0.8),
      strip.text.x = element_text(size = 10, face = "bold")
    )
  return(p)
}

# ========================
# 4. Create volcano plots
# ========================
p_W_all <- make_simple_volcano(deg_all_W, mycol_W)
p_D_all <- make_simple_volcano(deg_all_D, mycol_D)

# Combine
combined_plot_all <- p_D_all / p_W_all

# Show plot
combined_plot_all

# Save
ggsave("./Full_DEG_Volcano_Facets.pdf", plot = combined_plot_all, width = 16, height = 8, dpi = 300)




make_simple_volcano <- function(deg_all, mycol) {
  highlight_df <- deg_all %>% filter(!is.na(display_label))
  
  p <- ggplot() +
    geom_point(data = deg_all, aes(x = log2FoldChange, y = -log10(padj)), size = 0.8, color = "grey") +
    geom_text_repel(
      data = highlight_df,
      aes(x = log2FoldChange, y = -log10(padj), label = display_label, color = Group),
      size = 3, max.overlaps = 30, show.legend = FALSE
    ) +
    geom_vline(xintercept = c(-1, 1), color = "grey50", linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), color = "grey50", linetype = "dashed") +
    coord_flip(xlim = c(-12, 12)) +
    facet_grid(. ~ Group, scales = "fixed") +
    scale_color_manual(values = mycol) +
    labs(x = "log2 Fold Change", y = "-log10 Adjusted p-value") +
    theme_bw() +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      axis.text = element_text(size = 10),
      axis.text.x = element_text(angle = 45, vjust = 0.8),
      strip.text.x = element_text(size = 10, face = "bold")
    )
  return(p)
}
######
p_W_all <- make_simple_volcano(deg_all_W, mycol_W)
p_D_all <- make_simple_volcano(deg_all_D, mycol_D)

combined_plot_all <- p_W_all / p_D_all

combined_plot_all




library(patchwork)

# Remove x-axis label for p_D_all:
p_D_all_clean <- p_D_all + 
  theme(
    axis.title.x = element_blank(),  # remove x label
    plot.margin = margin(t = 5, b = 5) # small margin if needed
  )

# Combine with small space
combined_plot_all <- p_D_all_clean / p_W_all + 
  plot_layout(heights = c(1, 1), guides = "collect") & 
  theme(plot.margin = margin(2, 2, 2, 2))  # global margin

# Show
combined_plot_all
#ggsave("./Full_DEG_Volcano_Facets.pdf", plot = combined_plot_all, width = 10, height = 8, dpi = 300)





#################cytoknie bargraph#####################
library(ggplot2)
library(dplyr)
library(ggpubr)
library(multcompView)
library(tidyr)

#Cytokine.list <- c("CCL22", "CHI3L1", "CCL8", "IFNB1", "CCL7", "CCL8", "IL6", "CXCL2", "CXCL1", "CXCL2", "CXCL3", "IL1B")
#Cytokine.list <- c("CCL7", "CCL8")

Cytokine.list <- c("TNF","IL1B","CCL8", "CXCL13", "IL17A.IL17F", "IL12B", "IL6", "CSF2")
# The actual treatment names in your desired order
treatment_order <- c(
  "None", "Bb_W", "FN_W", "SP_W", "VP_W", "SO_W",
  "DN_W", "Bb_D", "FN_D", "SP_D", "VP_D", "SO_D", "LPS"
)

# The numeric labels you want to show
numeric_labels <- as.character(1:13)

data_sub <- data.long.treat %>%
  filter(Dilution == 10, cytokine %in% Cytokine.list)

data_sub <- data_sub %>%
  mutate(treatment = factor(treatment, levels = treatment_order))

# Define your color palette
mycolors <- c(
  "None"  = "#999999",  # grey
  "Bb_W"  = "#E69F00",  # orange
  "FN_W"  = "#56B4E9",  # light blue
  "SP_W"  = "#009E73",  # green
  "VP_W"  = "#F0E442",  # yellow
  "SO_W"  = "#D55E00",  # dark orange
  "DN_W"  = "#CC79A7",  # pink
  "Bb_D"  = "#B2182B",  # dark red
  "FN_D"  = "#2166AC",  # dark blue
  "SP_D"  = "#4DAF4A",  # medium green
  "VP_D"  = "#FFD92F",  # bright yellow
  "SO_D"  = "#984EA3",  # purple
  "LPS" = "#8DD3C7"   # teal (your 13th condition)
)

# Create a reusable plotting function
cytokine_anova_cld_plot <- function(data, cytokine_name, dilution_filter = "10") {
  
  # Subset data
  data_sub <- data %>% filter(Dilution == dilution_filter, cytokine == cytokine_name)
  
  # One-way ANOVA
  aov.test <- aov(value ~ treatment, data = data_sub)
  summary.aov <- summary(aov.test)
  pval <- summary.aov[[1]][["Pr(>F)"]][1]
  
  # Tukey post-hoc
  tukey <- TukeyHSD(aov.test)
  
  # Compact letter display
  tukey_letters <- multcompLetters4(aov.test, tukey)
  cld <- as.data.frame.list(tukey_letters$treatment)
  cld$treatment <- rownames(cld)
  rownames(cld) <- NULL
  colnames(cld)[1] <- "Letter"
  
  # Merge CLD for plotting
  data_for_plot <- data_sub %>%
    left_join(cld, by = "treatment")
  
  # Calculate y positions for letters
  letter_positions <- data_for_plot %>%
    group_by(treatment, Letter) %>%
    summarise(y.position = max(value) * 1.05, .groups = 'drop')
  
  # Create plot
  p <- ggplot(data_sub, aes(x = treatment, y = value)) + 
    geom_jitter(color = "black", width = 0.15, size = 1.5) + 
    geom_boxplot(aes(fill = treatment), alpha = 0.6, outlier.alpha = 0, size = 0.2) +
    geom_text(data = letter_positions, 
              aes(x = treatment, y = y.position, label = Letter),
              size = 4, vjust = 0) +  
    theme_minimal() + 
    labs(title = cytokine_name,
         subtitle = paste0("ANOVA p = ", signif(pval, 3))) + 
    theme(legend.position = "none", text = element_text(size = 14)) +
    scale_fill_manual(values = mycolors) +
    scale_x_discrete(labels = numeric_labels)
  
  return(p)
}

# Now loop through your cytokines
for (Cytokine in Cytokine.list){
  print(Cytokine)
  
  pmono <- cytokine_anova_cld_plot(data = data.long.treat, cytokine_name = Cytokine)
  
  print(pmono)
  
  path.mono <- paste0("./result/CLD/Dilute10_", Cytokine, "_ANOVA_CLD.pdf")
  ggsave(filename = path.mono, plot = pmono, width = 3, height = 5, dpi = 300)
}





################################after DESeq, run PCA for all samples################################
library(dplyr)
library(FactoMineR)
library(factoextra)
library(ggplot2)

# Define your full color map
mycolors <- c(
  "None"  = "#999999",   # grey (keep)
  "Bb_W"  = "#F8766D",   # bright coral
  "FN_W"  = "#00BFC4",   # cyan
  "SP_W"  = "#7CAE00",   # lime green
  "VP_W"  = "#C77CFF",   # purple
  "SO_W"  = "#FF61CC",   # pink
  "DN_W"  = "#00BA38",   # vivid green
  "Bb_D"  = "#D55E00",   # dark orange
  "FN_D"  = "#1F78B4",   # strong blue
  "SP_D"  = "#B2182B",   # dark red
  "VP_D"  = "#F0E442",   # yellow (keep)
  "SO_D"  = "#984EA3",   # deep purple (keep)
  "LPS"   = "#8DD3C7"    # teal (keep)
)


# Subset treatments you are analyzing here
selected_treatments <- c("None", "Bb_D", "Bb_W", "LPS", "FN_W","FN_D","SP_W","SP_D","VP_W", "VP_D","SO_W","SO_D" )


filtered_data <- data.treat %>%
  filter(treatment %in% selected_treatments, Dilution == 10)

cytokine_data_dplyr <- filtered_data %>%
  dplyr::select(SampleID, AGER:WNT7A) %>%
  dplyr::select(-c(CCL21, IL5RA))

# Set SampleID as rownames
rownames(cytokine_data_dplyr) <- cytokine_data_dplyr$SampleID
cytokine_data_dplyr <- dplyr::select(cytokine_data_dplyr, -SampleID)

# Remove zero-variance columns
cytokine_data_dplyr <- cytokine_data_dplyr[, apply(cytokine_data_dplyr, 2, function(x) sd(x, na.rm = TRUE) > 0)]

# PCA
cytokine_normalized <- scale(cytokine_data_dplyr)
spleen.pca <- PCA(cytokine_normalized, graph = FALSE)

sample_names <- rownames(get_pca_ind(spleen.pca)$coord)

meta_matched <- filtered_data %>%
  filter(SampleID %in% sample_names) %>%
  arrange(match(SampleID, sample_names))

# PCA plot with customized colors
p3 <- fviz_pca_ind(
  spleen.pca,
  col.ind = meta_matched$treatment,
  label = "none",
  repel = TRUE,
  addEllipses = TRUE,
  ellipse.level = 0.8,
  geom = "point",      # use only points
  shape = 16           # solid circles
) +
  scale_color_manual(values = mycolors) +
  ggtitle("PCA: Cytokine Response (Dilution = 10)")

print(p3)

spleen.pca$ind$group <- meta_matched$treatment
p3 <- fviz_pca_ind(
  spleen.pca,
  habillage = "group",  # now using internal column
  geom = c("point"),
  pointshape = 16,
  addEllipses = TRUE,
  ellipse.level = 0.8,
  label = "none",
  repel = TRUE
) +
  scale_color_manual(values = mycolors) +
  ggtitle("PCA: Cytokine Response (Dilution = 10)")


# First, make sure meta_matched matches your PCA data
sample_names <- rownames(pca_res$x)
meta_matched <- sample_metadata %>%
  filter(SampleID %in% sample_names) %>%
  arrange(match(SampleID, sample_names))

# Define your treatment order
my_levels <- c("None", "LPS", "Bb_D", "Bb_W", 
               "FN_D", "FN_W", "SO_D", "SO_W", 
               "SP_D", "SP_W", "VP_W", "VP_D")

# Apply this order to your meta_matched data frame
meta_matched$treatment <- factor(meta_matched$treatment, levels = my_levels)

# Now draw PCA with ellipses:
p3 <- fviz_pca_ind(
  pca_res,
  habillage = meta_matched$treatment,  # group for coloring and ellipse
  addEllipses = TRUE,
  ellipse.level = 0.8,  # your confidence interval
  repel = TRUE,
  label = "none"
) + 
  scale_color_manual(values = mycolors) + 
  ggtitle("PCA: Cytokine Response (6 Treatments, Dilution = 10)")

print(p3)

#ggsave(filename = "../../Manuscript/Figures/Figure 6-other peptidoglycans/PCA_rawcount_all.pdf",p3, width =4, height = 4, dpi = 300)



#######try to break out None, Bb_D, Bb_W, LPS############




#######try to break out "None", "FN_W",  "SO_W", "SP_W", "VP_W",############
selected_treatments <- c("None","Bb_W", "FN_W",  "SO_W", "SP_W", "VP_W")

mycolors <- c(
  "None"  = "#999999",   # grey (keep)
  "Bb_W"  = "#F8766D",   # bright coral
  "FN_W"  = "#00BFC4",   # cyan
  "SP_W"  = "#7CAE00",   # lime green
  "VP_W"  = "#C77CFF",   # purple
  "SO_W"  = "#FF61CC"   # pink
)


filtered_data <- data.treat %>%
  filter(treatment %in% selected_treatments, Dilution == 10)

cytokine_data_dplyr <- filtered_data %>%
  dplyr::select(SampleID, AGER:WNT7A) %>%
  dplyr::select(-c(CCL21, IL5RA))

# Set SampleID as rownames
rownames(cytokine_data_dplyr) <- cytokine_data_dplyr$SampleID
cytokine_data_dplyr <- dplyr::select(cytokine_data_dplyr, -SampleID)

# Remove zero-variance columns
cytokine_data_dplyr <- cytokine_data_dplyr[, apply(cytokine_data_dplyr, 2, function(x) sd(x, na.rm = TRUE) > 0)]

# PCA
cytokine_normalized <- scale(cytokine_data_dplyr)
spleen.pca <- PCA(cytokine_normalized, graph = FALSE)

sample_names <- rownames(get_pca_ind(spleen.pca)$coord)

meta_matched <- filtered_data %>%
  filter(SampleID %in% sample_names) %>%
  arrange(match(SampleID, sample_names))

# PCA plot with customized colors
p3 <- fviz_pca_ind(
  spleen.pca,
  col.ind = meta_matched$treatment,
  label = "none",
  repel = TRUE,
  addEllipses = TRUE,
  ellipse.level = 0.8,
  geom = "point",      # use only points
) +
  scale_color_manual(values = mycolors) +
  ggtitle("PCA: Cytokine Response (Dilution = 10)")

print(p3)

spleen.pca$ind$group <- meta_matched$treatment
p3 <- fviz_pca_ind(
  spleen.pca,
  habillage = "group",  # now using internal column
  geom = c("point"),
  pointshape = 16,
  addEllipses = TRUE,
  ellipse.level = 0.8,
  label = "none",
  repel = TRUE
) +
  scale_color_manual(values = mycolors) +
  ggtitle("PCA: Cytokine Response (Dilution = 10)")


# First, make sure meta_matched matches your PCA data
sample_names <- rownames(pca_res$x)
meta_matched <- sample_metadata %>%
  filter(SampleID %in% sample_names) %>%
  arrange(match(SampleID, sample_names))

# Define your treatment order
my_levels <- c("None", "FN_W","SO_W", "SP_W", "VP_W")

# Apply this order to your meta_matched data frame
meta_matched$treatment <- factor(meta_matched$treatment, levels = my_levels)

# Now draw PCA with ellipses:
p3 <- fviz_pca_ind(
  pca_res,
  habillage = meta_matched$treatment,  # group for coloring and ellipse
  addEllipses = TRUE,
  ellipse.level = 0.8,  # your confidence interval
  repel = TRUE,
  label = "none"
) + 
  scale_color_manual(values = mycolors) + 
  ggtitle("PCA: Cytokine Response (6 Treatments, Dilution = 10)")

print(p3)

#ggsave(filename = "../../Manuscript/Figures/Figure 6-other peptidoglycans/PCA_rawcount_W.pdf",p3, width =3, height = 4, dpi = 300)



#######try to break out None, "FN_D",  "SO_D", "SP_D", "VP_D"############

selected_treatments <- c("None","Bb_D", "FN_D",  "SO_D", "SP_D", "VP_D")

mycolors <- c(
  "None"  = "#999999",   # grey (keep)
  "Bb_D"  = "#D55E00",   # dark orange
  "FN_D"  = "#1F78B4",   # strong blue
  "SP_D"  = "#B2182B",   # dark red
  "VP_D"  = "#F0E442",   # yellow (keep)
  "SO_D"  = "#984EA3"  # deep purple (keep)
)

filtered_data <- data.treat %>%
  filter(treatment %in% selected_treatments, Dilution == 10)

cytokine_data_dplyr <- filtered_data %>%
  dplyr::select(SampleID, AGER:WNT7A) %>%
  dplyr::select(-c(CCL21, IL5RA))

# Set SampleID as rownames
rownames(cytokine_data_dplyr) <- cytokine_data_dplyr$SampleID
cytokine_data_dplyr <- dplyr::select(cytokine_data_dplyr, -SampleID)

# Remove zero-variance columns
cytokine_data_dplyr <- cytokine_data_dplyr[, apply(cytokine_data_dplyr, 2, function(x) sd(x, na.rm = TRUE) > 0)]

# PCA
cytokine_normalized <- scale(cytokine_data_dplyr)
spleen.pca <- PCA(cytokine_normalized, graph = FALSE)

sample_names <- rownames(get_pca_ind(spleen.pca)$coord)

meta_matched <- filtered_data %>%
  filter(SampleID %in% sample_names) %>%
  arrange(match(SampleID, sample_names))

# PCA plot with customized colors
p3 <- fviz_pca_ind(
  spleen.pca,
  col.ind = meta_matched$treatment,
  label = "none",
  repel = TRUE,
  addEllipses = TRUE,
  ellipse.level = 0.8,
  geom = "point",      # use only points
) +
  scale_color_manual(values = mycolors) +
  ggtitle("PCA: Cytokine Response (Dilution = 10)")

print(p3)

#ggsave(filename = "../../Manuscript/Figures/Figure 6-other peptidoglycans/PCA_rawcount_D.pdf",p3, width =3, height = 4, dpi = 300)
# First, make sure meta_matched matches your PCA data
sample_names <- rownames(pca_res$x)
meta_matched <- sample_metadata %>%
  filter(SampleID %in% sample_names) %>%
  arrange(match(SampleID, sample_names))

# Define your treatment order
my_levels <- c("None", "FN_D","SO_D", "SP_D", "VP_D", "Bb_D")

# Apply this order to your meta_matched data frame
meta_matched$treatment <- factor(meta_matched$treatment, levels = my_levels)

# Now draw PCA with ellipses:
p3 <- fviz_pca_ind(
  pca_res,
  habillage = meta_matched$treatment,  # group for coloring and ellipse
  addEllipses = TRUE,
  ellipse.level = 0.8,  # your confidence interval
  repel = TRUE,
  label = "none"
) + 
  scale_color_manual(values = mycolors) + 
  ggtitle("PCA: Cytokine Response (6 Treatments, Dilution = 10)")

print(p3)

#ggsave(filename = "../../Manuscript/Figures/Figure 6-other peptidoglycans/PCA_rawcount_D.pdf",p3, width =3, height = 4, dpi = 300)





#######try to break out None, "FN_D",  "SO_D", "SP_D", "VP_D"############

selected_treatments <- c("None","Bb_D", "Bb_W",  "LPS")
mycolors <- c(
  "None"  = "#999999",   # grey (keep)
  "Bb_W"  = "#F8766D",   # bright coral
  "Bb_D"  = "#D55E00",   # dark orangE
  "LPS"   = "#8DD3C7"    # teal (keep)
)

filtered_data <- data.treat %>%
  filter(treatment %in% selected_treatments, Dilution == 10)

cytokine_data_dplyr <- filtered_data %>%
  dplyr::select(SampleID, AGER:WNT7A) %>%
  dplyr::select(-c(CCL21, IL5RA))

# Set SampleID as rownames
rownames(cytokine_data_dplyr) <- cytokine_data_dplyr$SampleID
cytokine_data_dplyr <- dplyr::select(cytokine_data_dplyr, -SampleID)

# Remove zero-variance columns
cytokine_data_dplyr <- cytokine_data_dplyr[, apply(cytokine_data_dplyr, 2, function(x) sd(x, na.rm = TRUE) > 0)]

# PCA
cytokine_normalized <- scale(cytokine_data_dplyr)
spleen.pca <- PCA(cytokine_normalized, graph = FALSE)

sample_names <- rownames(get_pca_ind(spleen.pca)$coord)

meta_matched <- filtered_data %>%
  filter(SampleID %in% sample_names) %>%
  arrange(match(SampleID, sample_names))

# PCA plot with customized colors
p3 <- fviz_pca_ind(
  spleen.pca,
  col.ind = meta_matched$treatment,
  label = "none",
  repel = TRUE,
  addEllipses = TRUE,
  ellipse.level = 0.8,
  geom = "point",      # use only points
) +
  scale_color_manual(values = mycolors) +
  ggtitle("PCA: Cytokine Response (Dilution = 10)")

print(p3)

#ggsave(filename = "../../Manuscript/Figures/Figure 6-other peptidoglycans/PCA_rawcount_D.pdf",p3, width =3, height = 4, dpi = 300)
# First, make sure meta_matched matches your PCA data
sample_names <- rownames(pca_res$x)
meta_matched <- sample_metadata %>%
  filter(SampleID %in% sample_names) %>%
  arrange(match(SampleID, sample_names))

# Define your treatment order
my_levels <- c("None", "Bb_D","Bb_W", "LPS")

# Apply this order to your meta_matched data frame
meta_matched$treatment <- factor(meta_matched$treatment, levels = my_levels)

# Now draw PCA with ellipses:
p3 <- fviz_pca_ind(
  pca_res,
  habillage = meta_matched$treatment,  # group for coloring and ellipse
  addEllipses = TRUE,
  ellipse.level = 0.8,  # your confidence interval
  repel = TRUE,
  label = "none"
) + 
  scale_color_manual(values = mycolors) + 
  ggtitle("PCA: Cytokine Response (6 Treatments, Dilution = 10)")

print(p3)

#ggsave(filename = "../../Manuscript/Figures/Figure 6-other peptidoglycans/PCA_rawcount_D.pdf",p3, width =3, height = 4, dpi = 300)

levels(meta_matched$treatment)

