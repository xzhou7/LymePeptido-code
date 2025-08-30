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
library(ggunchull)
library(ggrepel)
library(tidydr)
library(ggsci)
library(Cairo)
library(DoubletFinder)
library(harmony)
library(ggpubr)
library(scales)
library(clusterProfiler)
library(org.Hs.eg.db)
library(EnhancedVolcano)

library(clusterProfiler)
#choose the label genes
# Set working directory
setwd("~/Library/CloudStorage/Dropbox/lyme_disease/Pepti_Spleen_singlecell/")

# Load DEG data
spleen_DEG <- read.csv("../Manuscript/Figures/Figure 5-immune organoids and autoreactive B cells/Spleen_DEG.csv", 
                       header = TRUE, row.names = 1)

# Separate SYMBOL and ENSEMBL genes
Peptido_DEG <- spleen_DEG

# Convert SYMBOL genes
converted_genes <- bitr(Peptido_DEG$gene, 
                        fromType = "SYMBOL", 
                        toType = c("ENSEMBL", "ENTREZID"), 
                        OrgDb = org.Hs.eg.db)

Peptido_DEG_annotated <- merge(Peptido_DEG, converted_genes, 
                               by.x = "gene", by.y = "SYMBOL", 
                               all.x = TRUE)
Peptido_DEG_annotated$SYMBOL <- Peptido_DEG_annotated$gene
Peptido_DEG_annotated <- Peptido_DEG_annotated[!grepl("^ENSG", Peptido_DEG_annotated$SYMBOL), ]

# Convert ENSEMBL-only genes
Peptido_DEG_ensembl <- Peptido_DEG[grepl("^ENSG", Peptido_DEG$gene), ]

covert_gene1 <- bitr(Peptido_DEG_ensembl$gene, 
                     fromType = "ENSEMBL", 
                     toType = c("SYMBOL", "ENTREZID"), 
                     OrgDb = org.Hs.eg.db)

Peptido_DEG_ensembl_ano <- merge(Peptido_DEG_ensembl, covert_gene1, 
                                 by.x = "gene", by.y = "ENSEMBL", 
                                 all.x = TRUE)

# Align columns
common_cols <- intersect(colnames(Peptido_DEG_annotated), colnames(Peptido_DEG_ensembl_ano))
Peptido_DEG_ensembl_ano <- Peptido_DEG_ensembl_ano[, common_cols]
Peptido_DEG_annotated <- Peptido_DEG_annotated[, common_cols]

# Combine and filter
Peptido_DEG_PATH <- rbind(Peptido_DEG_annotated, Peptido_DEG_ensembl_ano) %>%
  filter(p_val_adj <= 0.05)

# Extract DEGs by cell type
deg_list <- list(
  Bcell = filter(Peptido_DEG_PATH, celltype == "B_Cell"),
  CD4_T = filter(Peptido_DEG_PATH, celltype == "CD4T_Cell"),
  CD8_T = filter(Peptido_DEG_PATH, celltype == "CD8T_Cell"),
  Macrophage_DC = filter(Peptido_DEG_PATH, celltype == "Macrophage_DC"),
  NK = filter(Peptido_DEG_PATH, celltype == "NK_Cell")
)

# Significant up/downregulated genes
#sig_up <- lapply(deg_list, function(df) filter(df, avg_log2FC > 0.25) %>% pull(gene))
#sig_down <- lapply(deg_list, function(df) filter(df, avg_log2FC < -0.25) %>% pull(gene))

# Count and filter shared genes
#threshold_shared <- 5
#genes_up_shared <- names(table(unlist(sig_up))[table(unlist(sig_up)) >= threshold_shared])
#genes_down_shared <- names(table(unlist(sig_down))[table(unlist(sig_down)) >= threshold_shared])



#######check the top 20 genes (logFC)in each cel tyeps

 #Remove duplicates by gene and celltype before ranking
Peptido_DEG_PATH <- Peptido_DEG_PATH %>%
  arrange(gene, celltype, desc(abs(avg_log2FC))) %>%
  distinct(gene, celltype, .keep_all = TRUE)

 #Get top N upregulated genes per cell type (e.g., top 10)
top_genes_per_celltype <- Peptido_DEG_PATH %>%
  filter(avg_log2FC > 0.25) %>%
  group_by(celltype) %>%
  top_n(20, wt = avg_log2FC) %>%
  arrange(celltype, desc(avg_log2FC))



# View results
top_genes_per_celltype
#write.csv(file = "../Pepti_Spleen_singlecell/Result/top_genes_per_celltype.csv", top_genes_per_celltype )


# Get top N downregulated genes per cell type (e.g., top 10)
top_genes_per_celltype_down <- Peptido_DEG_PATH %>%
  filter(avg_log2FC < -0.25) %>%
  group_by(celltype) %>%
  top_n(20, wt = avg_log2FC) %>%
  arrange(celltype, desc(avg_log2FC))

# View results
top_genes_per_celltype_down
#write.csv(file = "../Pepti_Spleen_singlecell/Result/top_genes_per_celltype_down.csv", top_genes_per_celltype_down )

# Combine and clean gene list
#genes_combined <- unique(c(genes_up_shared, genes_down_shared))
#genes_combined <- genes_combined[!grepl("^ENSG", genes_combined)]

# Add manually selected genes
#new_genes <- c("CCL22", "CCL3", "CXCL8", "EBI3", "CXCL5", "CD83", "CXCL1", "CXCL2", "CXCL3", "CCL8", "IL10", "IL1A")

#genes selected from top 20 up
top20_up <- c("CXCL8", "IL1B", "CCL17", "CCL22", "EBI3", "CXCL8", "IL1B", "DUSP4", "GZMB", "HLA-DQA1", "HLA-DRA", "PFKFB3", 
              "CREM", "HLA-DPB1", "BHLHE40", "HLA-DPA1", "FCRL5", "CXCL5", "CXCL1", "CXCL3", "TNFAIP6", "CCL18", "CD300E", 
              "C3", "PDPN", "IL22")
TOP20_down <- c("ARRB2", "FAM210A", "SDHD", "GFM1", "DLAT", "TXNIP", "NKG7", "GBP5", "IFI44", "GBP1", "STAT1", "GBP2")

new_genes <- c("KIR2DL4", "TNF", "CXCL2", "")
#genes_combined <- Reduce(union, list(genes_combined, new_genes, top20_up, TOP20_down))
genes_combined <- Reduce(union, list(top20_up, TOP20_down, new_genes))

# Remove unwanted genes
genes_to_remove <- c("RPP21", "TSEN34", "RGL2", "PSMB8-AS1", "LENG1", "GPANK1", "CNOT3", 
                     "SAMHD1", "HLA-B", "LSM2", "RDH13", "SNHG32", "MDC1", "PRR3", "OSCAR", "FCGR3A")
genes_combined <- setdiff(genes_combined, genes_to_remove)

# Match and deduplicate genes
matched_peptido_genes <- Peptido_DEG_PATH %>%
  filter(gene %in% genes_combined) %>%
  distinct(gene, celltype, .keep_all = TRUE)

#write.csv(file = "../Pepti_Spleen_singlecell/Result/matched_peptido_genes.csv", matched_peptido_genes )



# Dynamically define gene sets based on celltype in Peptido_DEG_PATH
CD4T_genes      <- filter(matched_peptido_genes, celltype == "CD4T_Cell")$gene
CD8T_genes      <- filter(matched_peptido_genes, celltype == "CD8T_Cell")$gene
NK_genes        <- filter(matched_peptido_genes, celltype == "NK_Cell")$gene
B_genes         <- filter(matched_peptido_genes, celltype == "B_Cell")$gene
Macro_genes     <- filter(matched_peptido_genes, celltype == "Macrophage_DC")$gene

# Create the plot
p <- Peptido_DEG_PATH %>%
  filter(p_val_adj < 0.05) %>%
  ggplot(aes(x = celltype, y = avg_log2FC, label = gene)) +
  
  # Color genes from matched celltypes in black, others in light gray
  geom_jitter(aes(color = (gene %in% CD4T_genes     & celltype == "CD4T_Cell"& p_val_adj < 0.05) |
                    (gene %in% CD8T_genes     & celltype == "CD8T_Cell"& p_val_adj < 0.05) |
                    (gene %in% NK_genes       & celltype == "NK_Cell"& p_val_adj < 0.05) |
                    (gene %in% B_genes        & celltype == "B_Cell"& p_val_adj < 0.05) |
                    (gene %in% Macro_genes    & celltype == "Macrophage_DC")& p_val_adj < 0.05)) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "#DCDCDC")) +
  # Label genes from each cell type
  geom_text_repel(data = subset(Peptido_DEG_PATH, gene %in% CD4T_genes & celltype == "CD4T_Cell" & p_val_adj < 0.05), 
                  vjust = -1, hjust = 0.5, color = "black", size = 5) +
  geom_text_repel(data = subset(Peptido_DEG_PATH, gene %in% CD8T_genes & celltype == "CD8T_Cell" & p_val_adj < 0.05), 
                  vjust = -1, hjust = 0.5, color = "black", size = 5) +
  geom_text_repel(data = subset(Peptido_DEG_PATH, gene %in% NK_genes & celltype == "NK_Cell" & p_val_adj < 0.05), 
                  vjust = -1, hjust = 0.5, color = "black", size = 5) +
  geom_text_repel(data = subset(Peptido_DEG_PATH, gene %in% B_genes & celltype == "B_Cell" & p_val_adj < 0.05), 
                  vjust = -1, hjust = 0.5, color = "black", size = 5) +
  geom_text_repel(data = subset(Peptido_DEG_PATH, gene %in% Macro_genes & celltype == "Macrophage_DC" & p_val_adj < 0.05), 
                  vjust = -1, hjust = 0.5, color = "black", size = 5) +
  
  # Theme and axis
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x  = element_text(size = 14, angle = 0, hjust = 0.5),
    axis.text.y  = element_text(size = 14),
    plot.title   = element_text(size = 18, face = "bold")
  ) +
  labs(title = "Differentially Expressed Genes by Cell Types (Peptidoglycan Stimulation)",
       x = "Cell Type", y = "Log2 Fold Change") +
  #coord_cartesian(ylim = c(-2, 8)) +
  guides(color = FALSE)
p

# Save plot
ggsave(filename = "../Manuscript/Figures/Figure 5-immune organoids and autoreactive B cells/Peptido_DEG_Plot1.pdf", 
       plot = p, width = 10, height = 13, dpi = 300)

#library(magick)

# Read PNG
#img <- image_read("../Manuscript/Figures/Figure 5-immune organoids and autoreactive B cells/Peptido_DEG_Plot1.png")

# Save directly as PDF
#image_write(img, path = "../Manuscript/Figures/Figure 5-immune organoids and autoreactive B cells/Peptido_DEG_Plotpng.pdf", format = "pdf")


############Final code for paper


# Load necessary libraries
library(dplyr)
library(ggplot2)
library(ggrepel)

# Step 1: Prepare data
Peptido_DEG_PATH <- Peptido_DEG_PATH %>%
  arrange(gene, celltype, desc(abs(avg_log2FC))) %>%
  distinct(gene, celltype, .keep_all = TRUE)

# Step 2: Define gene sets
#new_genes <- c("CCL22", "CCL3", "CXCL8", "EBI3", "CXCL5", "CD83", "CXCL1", "CXCL2", "CXCL3", "CCL8", "IL10", "IL1A")

#top20_up1 <- c("CXCL8", "IL1B", "CCL17", "CCL22", "EBI3", "CXCL8", "IL1B", "DUSP4", "GZMB", 
#              "HLA-DQA1", "HLA-DRA", "PFKFB3", "CREM", "HLA-DPB1", "BHLHE40", "HLA-DPA1", 
#              "FCRL5", "CXCL5", "CXCL1", "CXCL3", "TNFAIP6", "CCL18", "CD300E", "IL22")

top20_up <- c( "IL1B", "EBI3",  "IL1B", "DUSP4","GZMB",
              "PFKFB3", "CREM", "CXCL1", "CXCL3", "TNFAIP6")
TOP20_down <- c( "FAM210A", "SDHD", "GFM1", "DLAT", "TXNIP")
new_genes <- c("KIR2DL4", "TNF", "CCL4", "CCL3", "CD83", "IRG1", "IFB1", "EDN1", "DUSP8", "CCL5", 
                "CD80", "SLC37A2", "HMOX2", "COQ9", "CXCL8", "IL22", "CXCL5")
genes_to_label <- unique(c(top20_up, TOP20_down, new_genes))

# Step 3: Filter and prepare labeled data
matched_peptido_genes <- Peptido_DEG_PATH %>%
  filter(gene %in% genes_to_label) %>%
  distinct(gene, celltype, .keep_all = TRUE)

# Step 4: Plot
p <- Peptido_DEG_PATH %>%
  filter(p_val_adj < 0.05) %>%
  ggplot(aes(x = celltype, y = avg_log2FC, label = gene)) +
  
  # Highlight only the selected genes
  geom_jitter(aes(color = gene %in% genes_to_label), size = 1.5, width = 0.25, alpha = 0.8) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "#DCDCDC")) +
  
  #Label selected genes
  geom_text_repel(data = subset(Peptido_DEG_PATH, gene %in% genes_to_label & p_val_adj < 0.05), 
                  aes(x = celltype, y = avg_log2FC, label = gene),
                  color = "black", size = 5, max.overlaps = Inf, box.padding = 0.5) +
  
  # Theme and axes
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x  = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 14),
    plot.title   = element_text(size = 18, face = "bold")
  ) +
  #labs(title = "Differentially Expressed Genes by Cell Type (Peptidoglycan Stimulation)",
  #     x = "Cell Type", y = "Log2 Fold Change") +
  scale_y_continuous(limits = NULL) +
  guides(color = FALSE)
p

# Step 5: Save
ggsave(filename = "../Manuscript/Figures/Figure 5-immune organoids and autoreactive B cells/Peptido_DEG_Plot.pdf", 
       plot = p, width = 8, height = 8, dpi = 300)


##########reduce the size of the paper, so no lable#####


# Load necessary libraries
library(dplyr)
library(ggplot2)
library(ggrepel)

# Step 1: Prepare data
Peptido_DEG_PATH <- Peptido_DEG_PATH %>%
  arrange(gene, celltype, desc(abs(avg_log2FC))) %>%
  distinct(gene, celltype, .keep_all = TRUE)

top20_up <- c( "IL1B", "EBI3",  "IL1B", "DUSP4","GZMB",
               "PFKFB3", "CREM", "CXCL1", "CXCL3", "TNFAIP6")
TOP20_down <- c( "FAM210A", "SDHD", "GFM1", "DLAT", "TXNIP")
new_genes <- c("KIR2DL4", "TNF", "CCL4", "CCL3", "CD83", "IRG1", "IFB1", "EDN1", "DUSP8", "CCL5", 
               "CD80", "SLC37A2", "HMOX2", "COQ9", "CXCL8", "IL22", "CXCL5")
genes_to_label <- unique(c(top20_up, TOP20_down, new_genes))

# Step 3: Filter and prepare labeled data
matched_peptido_genes <- Peptido_DEG_PATH %>%
  filter(gene %in% genes_to_label) %>%
  distinct(gene, celltype, .keep_all = TRUE)

# Step 4: Plot

p <- Peptido_DEG_PATH %>%
  filter(p_val_adj < 0.05) %>%
  ggplot(aes(x = celltype, y = avg_log2FC, label = gene)) +
  
  geom_jitter(aes(color = gene %in% genes_to_label), 
              size = 1.5, width = 0.25, alpha = 0.8) +
  scale_color_manual(values = c("TRUE" = "black", "FALSE" = "#DCDCDC")) +
  
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),   # remove x-axis title
    axis.text.x  = element_blank(),   # remove x-axis tick text
    axis.ticks.x = element_blank(),   # (optional) remove tick marks
    axis.title.y = element_blank(),
    axis.text.y  = element_text(size = 14),
    plot.title   = element_text(size = 18, face = "bold")
  ) +
  scale_y_continuous(limits = NULL) +
  guides(color = FALSE)

p


#####
p <- Peptido_DEG_PATH %>%
  filter(p_val_adj < 0.05) %>%
  ggplot(aes(x = celltype, y = avg_log2FC)) +
  
  # First plot grey background points
  geom_jitter(
    data = subset(Peptido_DEG_PATH, !(gene %in% genes_to_label) & p_val_adj < 0.05),
    color = "#DCDCDC", size = 1.5, width = 0.25, alpha = 0.8
  ) +
  
  # Then plot black highlighted points on top
  geom_jitter(
    data = subset(Peptido_DEG_PATH, gene %in% genes_to_label & p_val_adj < 0.05),
    color = "black", size = 1.5, width = 0.25, alpha = 0.9
  ) +
  
  theme_minimal() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y  = element_text(size = 20),
    plot.title   = element_text(size = 18, face = "bold")
  ) +
  scale_y_continuous(limits = NULL)

p

######
# Step 5: Save
ggsave(filename = "../Manuscript/Figures/Figure 5-immune organoids and autoreactive B cells/Peptido_DEG_Plot.png", 
       plot = p, width = 8, height = 8, dpi = 300)






############

p <- Peptido_DEG_PATH %>%
  filter(p_val_adj < 0.05) %>%
  ggplot(aes(x = celltype, y = avg_log2FC)) +
  
  # First plot grey background points
  geom_jitter(
    data = subset(Peptido_DEG_PATH, !(gene %in% genes_to_label) & p_val_adj < 0.05),
    color = "#DCDCDC", size = 1.5, width = 0.25, alpha = 0.8
  ) +
  
  #Label selected genes
  geom_text_repel(data = subset(Peptido_DEG_PATH, gene %in% genes_to_label & p_val_adj < 0.05), 
                  aes(x = celltype, y = avg_log2FC, label = gene),
                  color = "black", size = 5, max.overlaps = Inf, box.padding = 0.5) +
  
  # Then plot black highlighted points on top
  geom_jitter(
    data = subset(Peptido_DEG_PATH, gene %in% genes_to_label & p_val_adj < 0.05),
    color = "black", size = 1.5, width = 0.25, alpha = 0.9
  ) +
  
  theme_minimal() +
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x  = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y  = element_text(size = 14),
    plot.title   = element_text(size = 18, face = "bold")
  ) +
  scale_y_continuous(limits = NULL)

p





  







