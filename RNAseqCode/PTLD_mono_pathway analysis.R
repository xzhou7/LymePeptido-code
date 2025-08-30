## Healthy Donor and PTLD Donor Monocytes RNAseq
# Updated: Mar 11, 2025
# Author: Xin Zhou, Ph.D.

# Load libraries
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(ReactomePA)
library(EnhancedVolcano)
library(enrichplot)
library(ggplot2)
library(ggtext)
library(dplyr)
library(stringr)
library(tidyverse)

# Set working directory (Mac)
setwd("~/Library/CloudStorage/Dropbox/lyme_disease/PTLD_HD_Monocytes/")

#-------------------------------
# 1️⃣ Load and preprocess data
#-------------------------------
gene.count <- read.csv(
  "./Data/Bulk_seq_momo_Count/aligned_to_hg19_-562735560/RNA_Express_AppResult-ds.72563998561b473690d3de501c4f12b1/differential/global/gene.counts.csv",
  header = TRUE, row.names = 1
)
meta.table <- read.csv("./Data/metadata.csv", header = TRUE, row.names = 1)

# Remove genes with 0 counts
gene.count.clean <- gene.count[rowSums(gene.count) > 0, ]

# Filter CD14 monocytes (Rest)
CD14_meta <- meta.table %>% filter(celltype == "CD14", status == "Rest")
rownames(CD14_meta) <- paste0("X", str_replace_all(rownames(CD14_meta), "-", "."))
CD14_cts <- gene.count[, rownames(CD14_meta)]
CD14_cts_clean <- CD14_cts[rowSums(CD14_cts) > 0, ]

#-------------------------------
# 2️⃣ DESeq2 Analysis
#-------------------------------
dds_CD14_rest <- DESeqDataSetFromMatrix(
  countData = CD14_cts_clean,
  colData = CD14_meta,
  design = ~ dex
)
dds_CD14_rest <- DESeq(dds_CD14_rest)
res_CD14_rest <- results(dds_CD14_rest, name = "dex_PTLD_vs_HD")

# Save DE results as dataframe
cd14_result.df <- as.data.frame(res_CD14_rest)
cd14_result.df$gene <- rownames(cd14_result.df)
cd14_result.df <- cd14_result.df %>% arrange(padj)

# Volcano plot (example)
EnhancedVolcano(cd14_result.df,
                lab = cd14_result.df$gene,
                x = 'log2FoldChange', y = 'pvalue') +
  ggtitle("CD14 Monocytes PTLD vs HD")

#-------------------------------
# 3️⃣ Prepare gene lists
#-------------------------------
universe.gene.list <- bitr(rownames(gene.count), 
                           fromType = "SYMBOL",
                           toType = "ENTREZID",
                           OrgDb = org.Hs.eg.db)

cd14.rest.up <- cd14_result.df %>% filter(padj < 0.05, log2FoldChange > 0)
cd14.rest.down <- cd14_result.df %>% filter(padj < 0.05, log2FoldChange < 0)

cd14.rest.up.df <- bitr(cd14.rest.up$gene, fromType = "SYMBOL",
                        toType = "ENTREZID", OrgDb = org.Hs.eg.db)
cd14.rest.down.df <- bitr(cd14.rest.down$gene, fromType = "SYMBOL",
                          toType = "ENTREZID", OrgDb = org.Hs.eg.db)

#-------------------------------
# 4️⃣ Enrichment Analysis
#-------------------------------

# GO (ALL ontologies)
ego.1.up <- enrichGO(gene = cd14.rest.up.df$ENTREZID,
                     universe = universe.gene.list$ENTREZID,
                     OrgDb = org.Hs.eg.db,
                     ont = "ALL",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.01,
                     qvalueCutoff = 0.05) %>%
  setReadable(OrgDb = org.Hs.eg.db)

# KEGG
kegg_up <- enrichKEGG(gene = cd14.rest.up.df$ENTREZID,
                      organism = 'hsa', pvalueCutoff = 0.05)
kegg_up_df <- as.data.frame(kegg_up)

# DO
dose_up <- enrichDO(gene = cd14.rest.up.df$ENTREZID,
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH",
                    minGSSize = 5, readable = TRUE)
dose_up_df <- as.data.frame(dose_up)

#-------------------------------
# 5️⃣ Define pathway categories
#-------------------------------
autoimmune_disease <- c(
  "systemic lupus erythematosus", "sjogren's syndrome",
  "systemic scleroderma", "rheumatoid arthritis",
  "arthritis", "vasculitis", "multiple sclerosis",
  "type 1 diabetes mellitus", "autoimmune disease",
  "inflammatory bowel disease"
)

signaling_pathway <- c(
  "TNF signaling pathway", "NF-kappa B signaling pathway",
  "C-type lectin receptor signaling pathway",
  "Cytokine-cytokine receptor interaction",
  "NOD-like receptor signaling pathway", "IL-17 signaling pathway",
  "Toll-like receptor signaling pathway",
  "Epstein-Barr virus infection", "Herpes simplex virus 1 infection"
)

viral_infections <- c(
  "Influenza A", "Epstein-Barr virus infection",
  "Herpes simplex virus 1 infection", "Human cytomegalovirus infection",
  "Coronavirus disease - COVID-19",
  "Human papillomavirus infection",
  "Human immunodeficiency virus 1 infection"
)

bacterial_infections <- c(
  "Tuberculosis", "Yersinia infection", "Salmonella infection",
  "Shigellosis", "Pathogenic E. coli infection",
  "Helicobacter pylori infection"
)

Biological_process <- c(
  "response to molecule of bacterial origin",
  "response to lipopolysaccharide",
  "positive regulation of cytokine production",
  "lymphocyte differentiation"
)

#-------------------------------
# 6️⃣ Filter pathways and combine
#-------------------------------
ego_Biological_process <- as.data.frame(ego.1.up) %>%
  filter(Description %in% Biological_process) %>%
  mutate(Source = "Biological Process")

dose_autoimmune <- dose_up_df %>% 
  filter(Description %in% autoimmune_disease) %>% 
  mutate(Source = "Autoimmune")

kegg_signaling <- kegg_up_df %>% 
  filter(Description %in% signaling_pathway) %>% 
  mutate(Source = "Signaling")

kegg_viral <- kegg_up_df %>% 
  filter(Description %in% viral_infections) %>% 
  mutate(Source = "Viral")

kegg_bacterial <- kegg_up_df %>% 
  filter(Description %in% bacterial_infections) %>% 
  mutate(Source = "Bacterial")

combined_df <- bind_rows(
  ego_Biological_process,
  dose_autoimmune,
  kegg_signaling,
  kegg_viral,
  kegg_bacterial
)

# Assign colors
combined_df$LabelColor <- case_when(
  combined_df$Source == "Autoimmune"         ~ "#B2182B",
  combined_df$Source == "Signaling"          ~ "#D95F0E",
  combined_df$Source == "Viral"              ~ "#4575B4",
  combined_df$Source == "Bacterial"          ~ "#1A9850",
  combined_df$Source == "Biological Process" ~ "#6A3D9A",
  TRUE ~ "black"
)

# Compute numeric GeneRatio
combined_df <- combined_df %>%
  mutate(GeneRatioNum = as.numeric(str_split_fixed(GeneRatio, "/", 2)[,1]) /
           as.numeric(str_split_fixed(GeneRatio, "/", 2)[,2])) %>%
  group_by(Source) %>%
  arrange(desc(GeneRatioNum), .by_group = TRUE) %>%
  ungroup()

# Colored HTML labels
combined_df$Description_colored <- paste0(
  "<span style='color:", combined_df$LabelColor, "'>",
  combined_df$Description, "</span>"
)
combined_df$Description_colored <- factor(
  combined_df$Description_colored,
  levels = rev(unique(combined_df$Description_colored))
)

#-------------------------------
# 7️⃣ Final dotplot
#-------------------------------
p_enrichment <- ggplot(combined_df, aes(
  x = GeneRatioNum,
  y = Description_colored,
  size = Count,
  color = p.adjust
)) +
  geom_point() +
  theme_bw() +
  labs(
    x = "Gene Ratio",
    y = NULL,
    size = "Gene Count",
    color = "Adjusted p-value",
    title = "Pathway Enrichment by Category (Largest GeneRatio on Top)"
  ) +
  scale_color_gradient(low = "darkred", high = "pink", trans = "log10") +
  theme(
    axis.text.y = element_markdown(size = 14),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 13),
    plot.title = element_text(size = 14, face = "bold")
  )

ggsave("PTLD_Mono_Pathway_Enrichment.pdf",
       plot = p_enrichment, width = 8, height = 9, dpi = 300)

ggsave("../Manuscript/Figures/Figure3-autoimmune monocytes/PTLD_Mono_Pathway_Enrichment.pdf",
       plot = p_enrichment, width = 8, height = 9, dpi = 300)



#########downregulation pathway###############

metabolic_diseases <- c("combined oxidative phosphorylation deficiency",
)

metabolic_pathways <- c("mitochondrial metabolism disease",
)


ego.1.down <- enrichGO(gene        = cd14.rest.down.df$ENTREZID,
                       universe      = universe.gene.list$ENTREZID,
                       OrgDb         = org.Hs.eg.db,
                       ont           = "ALL",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05)

summary(ego.1.down)
ego.1.down <- setReadable(ego.1.down, OrgDb = org.Hs.eg.db)
ego.1.down

ego.1.down2 <- pairwise_termsim(ego.1.down)
ego.1.down2.p1 <- treeplot(ego.1.down2)
ego.1.down2.p1
ego.1.down2.p2 <- treeplot(ego.1.down2, hclust_method = "average")
ego.1.down2.p2

dose_down <- enrichDO(gene = cd14.rest.down.df$ENTREZID, pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 5, readable = TRUE)

#write.csv(ego.1.down, file = "ego_PTLD_Down.csv", row.names = FALSE)

mito_terms <- c(
  "mitochondrial gene expression", "mitochondrial translation",
  "mitochondrial protein-containing complex", "mitochondrial ribosome",
  "mitochondrial matrix", "mitochondrial inner membrane",
  "mitochondrial small ribosomal subunit", "combined oxidative phosphorylation deficiency", 
  "mitochondrial metabolism disease"
)

ribosome_terms <- c(
  "ribonucleoprotein complex biogenesis", "ribosome biogenesis",
  "tRNA metabolic process", "organellar ribosome",
  "organellar small ribosomal subunit"
)

nuclear_terms <- c(
  "fibrillar center", "nuclear speck",
  "histone acetyltransferase complex", "ubiquitin ligase complex",
  "organelle inner membrane"
)

ego.1.down_df <- as.data.frame(ego.1.down)
dose_down_df <- as.data.frame(dose_down)


# Combine into a list for easier processing
all_terms <- c(mito_terms, ribosome_terms, nuclear_terms)

ego_down_filtered <- ego.1.down_df %>%
  filter(Description %in% all_terms) %>%
  mutate(Source = case_when(
    Description %in% mito_terms    ~ "Mitochondrial Function & Translation",
    Description %in% ribosome_terms ~ "Ribosome & RNA Biogenesis",
    Description %in% nuclear_terms  ~ "Nuclear / Regulatory Complexes",
    TRUE ~ "Other"
  ))

dose_down_filtered <- dose_down_df %>%
  filter(Description %in% all_terms) %>%
  mutate(Source = case_when(
    Description %in% mito_terms    ~ "Mitochondrial Function & Translation",
    Description %in% ribosome_terms ~ "Ribosome & RNA Biogenesis",
    Description %in% nuclear_terms  ~ "Nuclear / Regulatory Complexes",
    TRUE ~ "Other"
  ))

combined_down_df <- bind_rows(ego_down_filtered, dose_down_filtered)

# Compute numeric GeneRatio
combined_down_df <- combined_down_df %>%
  mutate(GeneRatioNum = as.numeric(str_split_fixed(GeneRatio, "/", 2)[,1]) /
           as.numeric(str_split_fixed(GeneRatio, "/", 2)[,2])) %>%
  group_by(Source) %>%
  arrange(desc(GeneRatioNum), .by_group = TRUE) %>%
  ungroup()

combined_down_df$LabelColor <- case_when(
  combined_down_df$Source == "Mitochondrial Function & Translation" ~ "#D95F0E", # orange
  combined_down_df$Source == "Ribosome & RNA Biogenesis"            ~ "#4575B4", # blue
  combined_down_df$Source == "Nuclear / Regulatory Complexes"       ~ "#1A9850", # green
  TRUE ~ "black"
)

combined_down_df$Description_colored <- paste0(
  "<span style='color:", combined_down_df$LabelColor, "'>",
  combined_down_df$Description, "</span>"
)

combined_down_df$Description_colored <- factor(
  combined_down_df$Description_colored,
  levels = rev(unique(combined_down_df$Description_colored))
)

library(ggplot2)
library(ggtext)

p_down <- ggplot(combined_down_df, aes(
  x = GeneRatioNum,
  y = Description_colored,
  size = Count,
  color = p.adjust
)) +
  geom_point() +
  theme_bw() +
  labs(
    x = "Gene Ratio",
    y = NULL,
    size = "Gene Count",
    color = "Adjusted p-value",
    title = "Downregulated Pathways (Grouped by Functional Category)"
  ) +
  scale_color_gradient(low = "darkred", high = "pink", trans = "log10") +
  theme(
    axis.text.y = element_markdown(size = 11),
    axis.text.x = element_text(size = 12),
    axis.title.x = element_text(size = 13),
    plot.title = element_text(size = 14, face = "bold")
  )

# Optional: save figure
ggsave("PTLD_Pathway_Downregulated_FunctionalGroups.pdf", 
       plot = p_down, width = 7, height = 5, dpi = 300)

ggsave("../Manuscript/Figures/Figure3-autoimmune monocytes/PTLD_Pathway_Downregulated_FunctionalGroups.pdf", 
       plot = p_down, width = 7, height = 5, dpi = 300)















