# Step 3.0 Pathway Analysis: Clean Version

# Load Libraries
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
library(ggsci)
library(Cairo)
library(DoubletFinder)
library(harmony)
library(ggpubr)
library(scales)
library(clusterProfiler)
library(org.Hs.eg.db)
library(EnhancedVolcano)
library(pheatmap)

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
sig_up <- lapply(deg_list, function(df) filter(df, avg_log2FC > 0.25) %>% pull(gene))
sig_down <- lapply(deg_list, function(df) filter(df, avg_log2FC < -0.25) %>% pull(gene))

# Count and filter shared genes
threshold_shared <- 5
genes_up_shared <- names(table(unlist(sig_up))[table(unlist(sig_up)) >= threshold_shared])
genes_down_shared <- names(table(unlist(sig_down))[table(unlist(sig_down)) >= threshold_shared])

# Combine and clean gene list
genes_combined <- unique(c(genes_up_shared, genes_down_shared))
genes_combined <- genes_combined[!grepl("^ENSG", genes_combined)]

# Add manually selected genes
new_genes <- c("CCL22", "CCL3", "CXCL8", "EBI3", "CXCL5", "CD83", "CXCL1", "CXCL2", "CXCL3", "CCL8", "IL10", "IL1A")
genes_combined <- union(genes_combined, new_genes)

# Remove unwanted genes
genes_to_remove <- c("RPP21", "TSEN34", "RGL2", "PSMB8-AS1", "LENG1", "GPANK1", "CNOT3", 
                     "SAMHD1", "HLA-B", "LSM2", "RDH13", "SNHG32", "MDC1", "PRR3", "OSCAR", "FCGR3A")
genes_combined <- setdiff(genes_combined, genes_to_remove)

# Match and deduplicate genes
matched_peptido_genes <- Peptido_DEG_PATH %>%
  filter(gene %in% genes_combined) %>%
  distinct(gene, celltype, .keep_all = TRUE)

write.csv(file = "../Pepti_Spleen_singlecell/Result/matched_peptido_genes.csv", matched_peptido_genes )
# Pivot to heatmap format
heatmap_df <- matched_peptido_genes %>%
  dplyr::select(gene, celltype, avg_log2FC) %>%
  pivot_wider(names_from = celltype, values_from = avg_log2FC)

# Convert to matrix
heatmap_mat <- as.data.frame(heatmap_df)
rownames(heatmap_mat) <- heatmap_mat$gene
heatmap_mat$gene <- NULL
heatmap_mat <- as.matrix(heatmap_mat)
heatmap_mat[!is.finite(heatmap_mat)] <- 0

# Plot and save heatmap
output_file <- "~/Library/CloudStorage/Dropbox/lyme_disease/Manuscript/Figures/Figure 5-immune organoids and autoreactive B cells/up_down_gene_in_all_subsets_New.pdf"

pheatmap(heatmap_mat,
         scale = "none",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         fontsize_row = 12,
         fontsize_col = 12,
         main = "avg_log2FC Heatmap of DEGs Across Cell Types",
         filename = output_file,
         width = 4,
         height = 10)


###########
# Custom breakpoints
breaks_blue <- seq(-1.5, 0, length.out = 50)
breaks_red  <- seq(0, 7, length.out = 51)  # +1 to include 0 in both
breaks <- c(breaks_blue, breaks_red[-1])  # remove duplicate 0

# Custom color palette: blue to white, then white to red
colors <- c(
  colorRampPalette(c("blue", "white"))(49),  # -1.5 to 0
  colorRampPalette(c("white", "red"))(50)    # 0 to 7
)

# Generate the heatmap
pheatmap(heatmap_mat,
         scale = "none",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colors,
         breaks = breaks,
         fontsize_row = 12,
         fontsize_col = 12,
         main = "avg_log2FC Heatmap of DEGs Across Cell Types",
         filename = output_file,
         width = 4,
         height = 10)







