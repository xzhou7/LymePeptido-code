#Step 3.0 Pathway analysis  

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
library(harmony)
library(ggpubr)
library(scales)
library(clusterProfiler)
library(org.Hs.eg.db)
library(EnhancedVolcano)

setwd("~/Library/CloudStorage/Dropbox/lyme_disease/Pepti_Spleen_singlecell/")
getwd()

Peptido_DEG <- read.csv("../Manuscript/Figures/Figure 5-immune organoids and autoreactive B cells/Spleen_DEG.csv", header = T, row.names = 1)


# Subset genes that look like ENSEMBL IDs
Peptido_DEG_ensembl <- Peptido_DEG[grepl("^ENSG", Peptido_DEG$gene), ]

covert_gene1 <- bitr(Peptido_DEG_ensembl$gene, 
                     fromType = "ENSEMBL", 
                     toType = c("SYMBOL", "ENTREZID"), 
                     OrgDb = org.Hs.eg.db)

Peptido_DEG_ensembl_ano <-  merge(Peptido_DEG_ensembl, covert_gene1, 
                                  by.x = "gene", by.y = "ENSEMBL", 
                                  all.x = TRUE)
dim(Peptido_DEG_ensembl_ano)
Peptido_DEG_ensembl_ano$ENSEMBL <- Peptido_DEG_ensembl_ano$gene
colnames(Peptido_DEG_ensembl_ano)
Peptido_DEG_ensembl_ano[grepl("^ENSG", Peptido_DEG_ensembl_ano$SYMBOL), ]

#calcuate all other genes
converted_genes <- bitr(Peptido_DEG$gene, 
                        fromType = "SYMBOL", 
                        toType = c("ENSEMBL", "ENTREZID"), 
                        OrgDb = org.Hs.eg.db)

Peptido_DEG_annotated <- merge(Peptido_DEG, converted_genes, 
                               by.x = "gene", by.y = "SYMBOL", 
                               all.x = TRUE)
Peptido_DEG_annotated$SYMBOL <- Peptido_DEG_annotated$gene

Peptido_DEG_annotated[grepl("^ENSG", Peptido_DEG_annotated$SYMBOL), ]
Peptido_DEG_annotated<- Peptido_DEG_annotated[!grepl("^ENSG", Peptido_DEG_annotated$SYMBOL), ]

colnames(Peptido_DEG_annotated)

Peptido_DEG_ensembl_ano <- Peptido_DEG_ensembl_ano[, colnames(Peptido_DEG_annotated)]

Peptido_DEG_PATH <- rbind(Peptido_DEG_ensembl_ano,Peptido_DEG_annotated) %>%  filter(p_val_adj  <= 0.05)
Peptido_DEG_PATH$SYMBOL
Peptido_DEG_PATH[grepl("^ENSG", Peptido_DEG_PATH$SYMBOL), ]
#there should be no value from SYMBOL that start with ENSG

universe.gene <- unique(na.omit(Peptido_DEG_PATH$ENTREZID))

#ggsave(filename = "./DEGandGO/CD8Tcell_DEG.pdf",P.CD8.deg, width = 20, height = 15, dpi = 300)

CD8T_up <- filter(Peptido_DEG_PATH, celltype == "CD8T_Cell" & avg_log2FC > 0 & p_val_adj < 0.05)
CD8T_down <- filter(Peptido_DEG_PATH, celltype == "CD8T_Cell" & avg_log2FC < 0 & p_val_adj < 0.05)


#cnet_down_CD8

#ggsave(filename = "./DEGandGO/CD8_T_GO_Down_CNET.pdf", cnet_down_CD8, width = 20, height = 15, dpi = 300)

#Kegg
CD8T_up.kegg <- enrichKEGG(
  gene         = unique(CD8T_up$ENTREZID),
  organism     = "hsa",
  universe     = universe.gene,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05)

CD8T_down.kegg <- enrichKEGG(
  gene         = unique(CD8T_down$ENTREZID),
  organism     = "hsa",
  universe     = universe.gene,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05)

CD8T_up.kegg.result <- CD8T_up.kegg@result
CD8T_up.kegg.result$Class <- "UpRegulated"
CD8T_down.kegg.result <- CD8T_down.kegg@result
CD8T_down.kegg.result$Class <- "DownRegulated"
CD8T_kegg.result <- rbind(CD8T_up.kegg.result, CD8T_down.kegg.result)
CD8T_kegg.result$CellType <- "CD8T_Cell"

#ggsave(filename = "DEGandGO/CD8_KEGG_UpDown.pdf",CD8T_kegg_dotplot, width = 10, height = 10, dpi = 300)

CD8T_up.kegg.plot <- dotplot(CD8T_up.kegg, title = "CD8T Upregulated KEGG", showCategory = 50)
CD8T_down.kegg.plot <- dotplot(CD8T_down.kegg, title = "CD8T Downregulated KEGG",showCategory = 50)
CD8T_kegg_dotplot <- CD8T_up.kegg.plot + CD8T_down.kegg.plot
CD8T_kegg_dotplot

#######################################################################################XC
# Convert the enrichment result to a data frame
CD8_up.kegg <- setReadable(CD8T_up.kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
CD8T_down.kegg <- setReadable(CD8T_down.kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
CD8_up_kegg_df <- as.data.frame(CD8_up.kegg)

# Apply setReadable to both up and down objects
CD8T_up.kegg <- setReadable(CD8T_up.kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")


# Now extract readable results
CD8T_up.kegg.result <- as.data.frame(CD8T_up.kegg)
CD8T_up.kegg.result$Class <- "UpRegulated"

CD8T_down.kegg.result <- as.data.frame(CD8T_down.kegg)
CD8T_down.kegg.result$Class <- "DownRegulated"

# Combine for dotplot
CD8T_kegg.result <- rbind(CD8T_up.kegg.result, CD8T_down.kegg.result)


# Generate dotplot using readable results
CD8T_kegg_dotplot <- dotplot(enrichResult = CD8T_kegg.result, showCategory = 30)

#GET RESULT:
head(CD8_up.kegg@result)
# Convert the enrichment result to a data frame
CD8_up_kegg_df <- as.data.frame(CD8_up.kegg)


# Write to CSV
#write.csv(
#  CD8_up_kegg_df,
#  file = "~/Library/CloudStorage/Dropbox/lyme_disease/Manuscript/Figures/Figure 5-immune organoids and autoreactive B cells/CD8_up_KEGG_results.csv",
#  row.names = FALSE
#)

#######################################################################################XC
CD4_up.kegg <- setReadable(CD4_up.kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
CD4_down.kegg <- setReadable(CD4_down.kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

CD4_up_kegg_plot <- dotplot(CD4_up.kegg, title = "CD4 T Cell KEGG - Upregulated", showCategory = 40)
CD4_down_kegg_plot <- dotplot(CD4_down.kegg, title = "CD4 T Cell KEGG - Downregulated")

CD4_kegg_dotplot <- CD4_up_kegg_plot + CD4_down_kegg_plot
CD4_kegg_dotplot

#GET RESULT:
head(CD4_up.kegg@result)
# Convert the enrichment result to a data frame
CD4_up_kegg_df <- as.data.frame(CD4_up.kegg)

# Write to CSV

#write.csv(
#  CD4_up_kegg_df,
#  file = "~/Library/CloudStorage/Dropbox/lyme_disease/Manuscript/Figures/Figure 5-immune organoids and autoreactive B cells/CD4_up_KEGG_results.csv",
#  row.names = FALSE
#)






#######


#DOSE
CD8T_up.DOSE <- enrichDO(
  gene          = unique(CD8T_up$ENTREZID),
  ont           = "HDO",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  universe      = universe.gene,
  minGSSize     = 5,
  maxGSSize     = 500,
  qvalueCutoff  = 0.05,
  readable      = TRUE)

CD8T_down.DOSE <- enrichDO(
  gene          = unique(CD8T_down$ENTREZID),
  ont           = "HDO",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  universe      = universe.gene,
  minGSSize     = 5,
  maxGSSize     = 500,
  qvalueCutoff  = 0.05,
  readable      = TRUE)

head(CD8T_up.DOSE)
head(CD8T_down.DOSE)

CD8T_up_DOSE_plot <- dotplot(CD8T_up.DOSE, title = "CD8 T Cell DOSE - Upregulated", showCategory = 25)
CD8T_down_DOSE_plot <- dotplot(CD8T_down.DOSE, title = "CD8 T Cell DOSE - Downregulated", showCategory = 25)

CD8_DOSE <- CD8T_up_DOSE_plot + CD8T_down_DOSE_plot
CD8_DOSE

CD8_union_ENTREZID <- unique(c(CD8T_up$ENTREZID, CD8T_down$ENTREZID)) %>% na.omit()
CD8_all.DOSE <- enrichDO(
  gene          = CD8_union_ENTREZID,
  ont           = "HDO",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  universe      = universe.gene,
  minGSSize     = 5,
  maxGSSize     = 500,
  qvalueCutoff  = 0.05,
  readable      = FALSE)

CD8_all_DOSE_plot <- dotplot(CD8_all.DOSE, title = "CD8 T Cell DOSE - all", showCategory = 50)
CD8_all_DOSE_plot

#CD8_dose.ego_cnet <- cnetplot(
#  CD8_all.DOSE,
#  categorySize = "pvalue",
#  showCategory = 50,
#  foldChange = NULL,
 # circular = FALSE,
#  colorEdge = TRUE,
#  node_label = "category") + ggtitle("CD8 T Cell DOSE")

#CD8_dose.ego_cnet

#ggsave(filename = "./DEGandGO/DOSE_CD8_All.pdf", CD8_dose.ego_cnet, width = 20, height = 15, dpi = 300)

# Subset up- and down-regulated genes for Macrophage_DC
MacDC_deg_results <-filter(Peptido_DEG_PATH, celltype == "Macrophage_DC")

MacDC_deg_results_unique <- MacDC_deg_results %>%
  filter(!is.na(gene) & !is.na(avg_log2FC) & !is.na(p_val_adj)) %>%
  distinct(gene, .keep_all = TRUE)

P.MacDC.deg <- EnhancedVolcano(MacDC_deg_results_unique,
                               lab = MacDC_deg_results_unique$gene,
                               #selectLab = c("CXCL5", "CCL22", "EBI3","CCL4", "CXCL8"),
                               FCcutoff=0.15,
                               # xlim = c(-5.5, 5.5),
                               x = 'avg_log2FC',
                               y = 'p_val_adj', 
                               title = "Macrophage_DC")
P.MacDC.deg

#ggsave(filename = "./DEGandGO/MacDC_DEG.pdf",P.MacDC.deg, width = 20, height = 15, dpi = 300)

MacDC_up <- filter(Peptido_DEG_PATH, celltype == "Macrophage_DC" & avg_log2FC > 0 & p_val_adj < 0.05)
MacDC_down <- filter(Peptido_DEG_PATH, celltype == "Macrophage_DC" & avg_log2FC < 0 & p_val_adj < 0.05)

MacDC_up.ego <- enrichGO(
  gene          = unique(MacDC_up$ENTREZID),
  OrgDb         = org.Hs.eg.db,
  universe      = universe.gene,
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05,
  readable      = TRUE)

MacDC_down.ego <- enrichGO(
  gene          = unique(MacDC_down$ENTREZID),
  OrgDb         = org.Hs.eg.db,
  universe      = universe.gene,
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05,
  readable      = TRUE)

MacDC_up.goresult <- MacDC_up.ego@result
MacDC_up.goresult$Class <- "UpRegulated"
MacDC_down.goresult <- MacDC_down.ego@result
MacDC_down.goresult$Class <- "DownRegulated"
MacDC_goresult <- rbind(MacDC_up.goresult, MacDC_down.goresult)
MacDC_goresult$CellType <- "Macrophage_DC"

MacDC_up.dotplot <- dotplot(MacDC_up.ego, title = "Macrophage_DC Upregulated GO", showCategory = 25)
MacDC_down.dotplot <- dotplot(MacDC_down.ego, title = "Macrophage_DC Downregulated GO", showCategory = 25)
MacDC_dotplot <- MacDC_up.dotplot + MacDC_down.dotplot
MacDC_dotplot

#cnet_up_MacDC <- cnetplot(
#  MacDC_up.ego,
#  categorySize = "pvalue",
#  showCategory = 50,
#  foldChange = NULL,
#  circular = FALSE,
#  colorEdge = TRUE,
#  node_label = "category"
#) + ggtitle("Macrophage/DC Upregulated GO (BP)")

#cnet_up_MacDC

#ggsave(filename = "./DEGandGO/MacDC_GO_Up_CNET.pdf", cnet_up_MacDC, width = 20, height = 15, dpi = 300)

#cnet_down_MacDC <- cnetplot(
#  MacDC_down.ego,
#  categorySize = "pvalue",
 # showCategory = 50,
#  foldChange = NULL,
#  circular = FALSE,
#  colorEdge = TRUE,
#  node_label = "category"
#) + ggtitle("Macrophage/DC Downregulated GO (BP)")

#cnet_down_MacDC

#ggsave(filename = "./DEGandGO/MacDC_GO_Down_CNET.pdf", cnet_down_MacDC, width = 20, height = 15, dpi = 300)

#Kegg
MacDC_up.kegg <- enrichKEGG(
  gene         = unique(MacDC_up$ENTREZID),
  organism     = "hsa",
  universe     = universe.gene,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05)

MacDC_down.kegg <- enrichKEGG(
  gene         = unique(MacDC_down$ENTREZID),
  organism     = "hsa",
  universe     = universe.gene,
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05)

MacDC_up.kegg.result <- MacDC_up.kegg@result
MacDC_up.kegg.result$Class <- "UpRegulated"
MacDC_down.kegg.result <- MacDC_down.kegg@result
MacDC_down.kegg.result$Class <- "DownRegulated"
MacDC_kegg.result <- rbind(MacDC_up.kegg.result, MacDC_down.kegg.result)
MacDC_kegg.result$CellType <- "Macrophage_DC"

MacDC_up.kegg.plot <- dotplot(MacDC_up.kegg, title = "Macrophage_DC Upregulated KEGG", showCategory = 50)
MacDC_down.kegg.plot <- dotplot(MacDC_down.kegg, title = "Macrophage_DC Downregulated KEGG", showCategory = 50)

MacDC_kegg_dotplot <- MacDC_up.kegg.plot + MacDC_down.kegg.plot
MacDC_kegg_dotplot

#ggsave(filename = "./DEGandGO/MacDC_KEGG_DOT.pdf", MacDC_kegg_dotplot, width = 10, height = 7, dpi = 300)

# Filter to BP-only terms from the result
MacDC_up_BP <- MacDC_up.ego
MacDC_up_BP@result <- dplyr::filter(MacDC_up.ego@result, ONTOLOGY == "BP")

MacDC_down_BP <- MacDC_down.ego
MacDC_down_BP@result <- dplyr::filter(MacDC_down.ego@result, ONTOLOGY == "BP")

MacDC_up_dotplot_BP <- dotplot(MacDC_up_BP, title = "Macrophage_DC Upregulated GO (BP only)")
MacDC_down_dotplot_BP <- dotplot(MacDC_down_BP, title = "Macrophage_DC Downregulated GO (BP only)")

MacDC_dotplot_BP <- MacDC_up_dotplot_BP + MacDC_down_dotplot_BP

MacDC_dotplot_BP

#cnet_up <- cnetplot(
#  MacDC_up_BP,
#  categorySize = "pvalue",
#  showCategory = 50,
#  foldChange = NULL,
#  circular = FALSE,
#  colorEdge = TRUE,
#  node_label = "category"  # Only show pathway labels
#) + ggtitle("Macrophage_DC Upregulated GO (BP)")

#cnet_up

#cnet_down <- cnetplot(
#  MacDC_down_BP,
#  categorySize = "pvalue",
##  showCategory = 50,
#  foldChange = NULL,
#  circular = FALSE,
#  colorEdge = TRUE,
#  node_label = "category"  # Only show pathway labels
#) + ggtitle("Macrophage_DC Downregulated GO (BP)")

#cnet_down

#DOSE
MacDC_up.DOSE <- enrichDO(
  gene          = unique(MacDC_up$ENTREZID),
  ont           = "HDO",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  universe      = universe.gene,
  minGSSize     = 5,
  maxGSSize     = 500,
  qvalueCutoff  = 0.05,
  readable      = TRUE)

MacDC_down.DOSE <- enrichDO(
  gene          = unique(MacDC_down$ENTREZID),
  ont           = "HDO",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  universe      = universe.gene,
  minGSSize     = 5,
  maxGSSize     = 500,
  qvalueCutoff  = 0.05,
  readable      = TRUE)

MacDC_up.dose.result <- MacDC_up.DOSE@result
MacDC_up.dose.result$Class <- "UpRegulated"
MacDC_down.dose.result <- MacDC_down.DOSE@result
MacDC_down.dose.result$Class <- "DownRegulated"
MacDC_dose.result <- rbind(MacDC_up.dose.result, MacDC_down.dose.result)
MacDC_dose.result$CellType <- "Macrophage_DC"

MacDC_up_dose_plot <- dotplot(MacDC_up.DOSE, title = "Macrophage_DC DOSE - Upregulated", showCategory = 25)
MacDC_down_dose_plot <- dotplot(MacDC_down.DOSE, title = "Macrophage_DC DOSE - Downregulated", showCategory = 25)

MacDC_dose_dotplot <- MacDC_up_dose_plot + MacDC_down_dose_plot
MacDC_dose_dotplot

#ggsave(filename = "./DEGandGO/MacDC_DOSE.pdf", MacDC_dose_dotplot, width = 10, height = 7, dpi = 300)

#all dose
MacDC_union_ENTREZID <- unique(c(MacDC_up$ENTREZID, MacDC_down$ENTREZID)) %>% na.omit()

MacDC_all.DOSE <- enrichDO(
  gene          = MacDC_union_ENTREZID %>% as.character(),
  ont           = "HDO",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  #universe      = universe.gene,
  minGSSize     = 5,
  maxGSSize     = 500,
  qvalueCutoff  = 0.05,
  readable      = FALSE)

MacDC_all.DOSE

MacDC_all_DOSE_plot <- dotplot(MacDC_all.DOSE, title = "Macrophage/DC DOSE - all", showCategory = 50)
MacDC_all_DOSE_plot


#######################################################################################XC

MacDC_all.DOSE <- setReadable(MacDC_all.DOSE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

# Convert the enrichment result to a data frame
MacDC_all_DOSE_df <- as.data.frame(MacDC_all.DOSE)

# Write to CSV

#write.csv(
#  MacDC_all_DOSE_df,
#  file = "~/Library/CloudStorage/Dropbox/lyme_disease/Manuscript/Figures/Figure 5-immune organoids and autoreactive B cells/MacDC_all_DOSE_results.csv",
##  row.names = FALSE
#)


#Define Your Autoimmune Disease Terms
autoimmune_terms <- c(
  "autoimmune disease", "lupus erythematosus", "systemic lupus erythematosus",
  "arthritis", "rheumatoid arthritis", "Sjogren's syndrome", "scleroderma",
  "systemic scleroderma", "multiple sclerosis", "autoimmune disease of blood",
  "hypersensitivity reaction type IV disease", "hypersensitivity reaction disease"
)

#Subset Your Pathway Results
# Assume your data frame is called pathway_df
autoimmune_pathways_MacDC <- MacDC_all_DOSE_df[MacDC_all_DOSE_df$Description %in% autoimmune_terms, ]

#Extract and Split Gene Lists
# Split geneIDs and unlist them
gene_lists <- strsplit(autoimmune_pathways_MacDC$geneID, split = "/")

# Flatten into one vector
all_genes <- unlist(gene_lists)

#Identify Repetitive Genes
# Count occurrences of each gene
gene_counts <- table(all_genes)

# Filter to get only repeated genes (i.e., appears in >1 pathway)
repeated_genes <- names(gene_counts[gene_counts > 4])


#filter(Peptido_DEG_PATH,SYMBOL == "FCGR3A")
#filter(Peptido_DEG_PATH,SYMBOL == "IL10")
#filter(Peptido_DEG_PATH,SYMBOL == "IL1A")


#######################################################################################XC







#MacDC_dose.ego_cnet <- cnetplot(
#  MacDC_all.DOSE,
#  categorySize = "pvalue",
#  showCategory = 50,
#  foldChange = NULL,
#  circular = FALSE,
#  colorEdge = TRUE,
#  node_label = "category") + ggtitle("Macrophage/DC DOSE")

#MacDC_dose.ego_cnet

#ggsave(filename = "./DEGandGO/DOSE_MacDC_All.pdf", MacDC_dose.ego_cnet, width = 20, height = 15, dpi = 300)

# Subset up- and down-regulated genes for NK_Cell
NK_deg_results <-filter(Peptido_DEG_PATH, celltype == "NK_Cell")

NK_deg_results_unique <- NK_deg_results %>%
  filter(!is.na(gene) & !is.na(avg_log2FC) & !is.na(p_val_adj)) %>%
  distinct(gene, .keep_all = TRUE)
P.NK.deg <- EnhancedVolcano(NK_deg_results_unique,
                            lab = NK_deg_results_unique$gene,
                            # selectLab = c("IL17F", "CCL22", "EBI3"),
                            FCcutoff=0.15,
                            xlim = c(-8, 8),
                            x = 'avg_log2FC',
                            y = 'p_val_adj', 
                            title = "NK Cell")
P.NK.deg

#ggsave(filename = "./DEGandGO/NK_DEG.pdf",P.NK.deg, width = 20, height = 15, dpi = 300)

NK_up <- filter(Peptido_DEG_PATH, celltype == "NK_Cell" & avg_log2FC > 0 & p_val_adj < 0.05)
NK_down <- filter(Peptido_DEG_PATH, celltype == "NK_Cell" & avg_log2FC < 0 & p_val_adj < 0.05)

NK_up.ego <- enrichGO(
  gene          = unique(NK_up$ENTREZID),
  OrgDb         = org.Hs.eg.db,
  universe      = universe.gene,
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05,
  readable      = TRUE)

NK_down.ego <- enrichGO(
  gene          = unique(NK_down$ENTREZID),
  OrgDb         = org.Hs.eg.db,
  universe      = universe.gene,
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05,
  readable      = TRUE)

NK_up.goresult <- NK_up.ego@result
NK_up.goresult$Class <- "UpRegulated"
NK_down.goresult <- NK_down.ego@result
NK_down.goresult$Class <- "DownRegulated"

NK_goresult <- rbind(NK_up.goresult, NK_down.goresult)
NK_goresult$CellType <- "NK_Cell"

NK_up.dotplot <- dotplot(NK_up.ego, title = "NK_Cell Upregulated GO")
NK_down.dotplot <- dotplot(NK_down.ego, title = "NK_Cell Downregulated GO")
NK_dotplot <- NK_up.dotplot + NK_down.dotplot
NK_dotplot

#cnet_up_NK <- cnetplot(
 # NK_up.ego,
 # categorySize = "pvalue",
 # showCategory = 50,
 # foldChange = NULL,
 # circular = FALSE,
#  colorEdge = TRUE,
 # node_label = "category"
#) + ggtitle("NK Cell Upregulated GO (BP)")

#cnet_up_NK

#ggsave(filename = "./DEGandGO/NK_GO_Up_CNET.pdf", cnet_up_NK, width = 20, height = 15, dpi = 300)

#cnet_down_NK <- cnetplot(
#  NK_down.ego,
#  categorySize = "pvalue",
#  showCategory = 50,
#  foldChange = NULL,
#  circular = FALSE,
#  colorEdge = TRUE,
#  node_label = "category"
#) + ggtitle("NK Cell Downregulated GO (BP)")

#cnet_down_NK

#ggsave(filename = "./DEGandGO/NK_GO_Down_CNET.pdf", cnet_down_NK, width = 20, height = 15, dpi = 300)

#Kegg
NK_up.kegg <- enrichKEGG(
  gene          = unique(NK_up$ENTREZID),
  organism      = "hsa",
  universe      = universe.gene,
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05)
NK_up.kegg

NK_down.kegg <- enrichKEGG(
  gene          = unique(NK_down$ENTREZID),
  organism      = "hsa",
  universe      = universe.gene,
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05)
NK_down.kegg

NK_up.kegg.result <- NK_up.kegg@result
NK_up.kegg.result$Class <- "UpRegulated"
NK_down.kegg.result <- NK_down.kegg@result
NK_down.kegg.result$Class <- "DownRegulated"
NK_kegg.result <- rbind(NK_up.kegg.result, NK_down.kegg.result)
NK_kegg.result$CellType <- "NK_Cell"

NK_up.kegg.plot <- dotplot(NK_up.kegg, title = "NK_Cell Upregulated KEGG", showCategory = 25)
NK_down.kegg.plot <- dotplot(NK_down.kegg, title = "NK_Cell Downregulated KEGG", showCategory = 25)

NK_kegg_dotplot <- NK_up.kegg.plot + NK_down.kegg.plot
NK_kegg_dotplot

#ggsave(filename = "DEGandGO/NK_KEGG_UpDown.pdf",NK_kegg_dotplot, width = 10, height = 10, dpi = 300)

#DOSE
NK_up.DOSE <- enrichDO(
  gene          = unique(NK_up$ENTREZID),
  #ont           = "HDO",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  universe      = universe.gene,
  minGSSize     = 5,
  maxGSSize     = 500,
  qvalueCutoff  = 0.05,
  readable      = TRUE)

NK_down.DOSE <- enrichDO(
  gene          = unique(NK_down$ENTREZID),
  #ont           = "HDO",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  universe      = universe.gene,
  minGSSize     = 5,
  maxGSSize     = 500,
  qvalueCutoff  = 0.05,
  readable      = TRUE)

NK_up.dose.result <- NK_up.DOSE@result
NK_up.dose.result$Class <- "UpRegulated"
NK_down.dose.result <- NK_down.DOSE@result
NK_down.dose.result$Class <- "DownRegulated"
NK_dose.result <- rbind(NK_up.dose.result, NK_down.dose.result)
NK_dose.result$CellType <- "NK_Cell"

head(NK_up.DOSE@result)
NK_up_dose_plot <- dotplot(NK_up.DOSE, title = "NK_Cell DOSE - Upregulated", showCategory = 25)
NK_up_dose_plot
NK_down_dose_plot <- dotplot(NK_down.DOSE, title = "NK_Cell DOSE - Downregulated", showCategory = 25)
NK_down_dose_plot

NK_dose_dotplot <- NK_up_dose_plot + NK_down_dose_plot
NK_dose_dotplot

#all
NK_union_ENTREZID <- unique(c(NK_up$ENTREZID, NK_down$ENTREZID)) %>% na.omit()
NK_all.DOSE <- enrichDO(
  gene          = as.character(NK_union_ENTREZID),
  ont           = "HDO",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  universe      = universe.gene,  # Uncomment if needed
  minGSSize     = 5,
  maxGSSize     = 500,
  qvalueCutoff  = 0.05,
  readable      = FALSE
)

NK_all_DOSE_plot <- dotplot(NK_all.DOSE, title = "NK Cell DOSE - all", showCategory = 50)
NK_all_DOSE_plot


###xc
######XC

NK_all.DOSE <- setReadable(NK_all.DOSE, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

# Convert the enrichment result to a data frame
NK_all.DOSE_df <- as.data.frame(NK_all.DOSE)

# Write to CSV

#write.csv(
#  NK_all.DOSE_df,
#  file = "~/Library/CloudStorage/Dropbox/lyme_disease/Manuscript/Figures/Figure 5-immune organoids and autoreactive B cells/NK_all.DOSE_results.csv",
#  row.names = FALSE
#)
#########

# Add cell type and source annotation
B_up.kegg_result <- B_up.kegg@result
B_up.kegg_result$CellType <- "B_cell"
B_up.kegg_result$Source <- "KEGG"
B_up.kegg_result$Class <- "UpRegulated"

B_down.kegg_result <- B_down.kegg@result 
B_down.kegg_result$CellType <- "B_cell"
B_down.kegg_result$Source <- "KEGG"
B_down.kegg_result$Class <- "DownRegulated"

CD4T_up.kegg_result <- CD4_up.kegg@result
CD4T_up.kegg_result$CellType <- "CD4_T"
CD4T_up.kegg_result$Source <- "KEGG"
CD4T_up.kegg_result$Class <- "UpRegulated"

CD4T_down.kegg_result <- CD4_down.kegg@result
CD4T_down.kegg_result$CellType <- "CD4_T"
CD4T_down.kegg_result$Source <- "KEGG"
CD4T_down.kegg_result$Class <- "DownRegulated"

CD8T_up.kegg_result <- CD8T_up.kegg@result
CD8T_up.kegg_result$CellType <- "CD8_T"
CD8T_up.kegg_result$Source <- "KEGG"
CD8T_up.kegg_result$Class <- "UpRegulated"

CD8T_down.kegg_result <- CD8T_down.kegg@result
CD8T_down.kegg_result$CellType <- "CD8_T"
CD8T_down.kegg_result$Source <- "KEGG"
CD8T_down.kegg_result$Class <- "DownRegulated"

MacDC_up.kegg_result <- MacDC_up.kegg@result
MacDC_up.kegg_result$CellType <- "Macrophage_DC"
MacDC_up.kegg_result$Source <- "KEGG"
MacDC_up.kegg_result$Class <- "UpRegulated"

MacDC_down.kegg_result <- MacDC_down.kegg@result
MacDC_down.kegg_result$CellType <- "Macrophage_DC"
MacDC_down.kegg_result$Source <- "KEGG"
MacDC_down.kegg_result$Class <- "DownRegulated"


NK_up.kegg_result <- NK_up.kegg@result
NK_up.kegg_result$CellType <- "NK"
NK_up.kegg_result$Source <- "KEGG"
NK_up.kegg_result$Class <- "UpRegulated"

NK_down.kegg_result <- NK_down.kegg@result
NK_down.kegg_result$CellType <- "NK"
NK_down.kegg_result$Source <- "KEGG"
NK_down.kegg_result$Class <- "DownRegulated"

B_all.DOSE_result <- B_all.DOSE@result
B_all.DOSE_result$CellType <- "B_cell"
B_all.DOSE_result$Source <- "DOSE"
B_all.DOSE_result$Class <-"Not_Applicable"

CD4T_all.DOSE_result <- CD4_all.DOSE@result
CD4T_all.DOSE_result$CellType <- "CD4_T"
CD4T_all.DOSE_result$Source <- "DOSE"
CD4T_all.DOSE_result$Class <- "Not_Applicable"

CD8T_all.DOSE_result <- CD8_all.DOSE@result
CD8T_all.DOSE_result$CellType <- "CD8_T"
CD8T_all.DOSE_result$Source <- "DOSE"
CD8T_all.DOSE_result$Class <- "Not_Applicable"

MacDC_all.DOSE_result <- MacDC_all.DOSE@result
MacDC_all.DOSE_result$CellType <- "Macrophage_DC"
MacDC_all.DOSE_result$Source <- "DOSE"
MacDC_all.DOSE_result$Class <- "Not_Applicable"

NK_all.DOSE_result <- NK_all.DOSE@result
NK_all.DOSE_result$CellType <- "NK"
NK_all.DOSE_result$Source <- "DOSE"
NK_all.DOSE_result$Class <- "Not_Applicable"

#write.csv(GO_all, file = "../Manuscript/Figures/Figure 5-immune organoids and autoreactive B cells/Combined_GO_Results.csv", row.names = FALSE)

#Enrichment_All <- bind_rows(
#  B_goresult, CD4T_goresult, CD8T_goresult, MacDC_goresult, NK_goresult,
#  B_up.kegg_result, B_down.kegg_result,
 # CD4T_up.kegg_result, CD4T_down.kegg_result,
#  CD8T_up.kegg_result, CD8T_down.kegg_result,
#  MacDC_up.kegg_result, MacDC_down.kegg_result,
#  NK_up.kegg_result, NK_down.kegg_result,
#  B_all.DOSE_result, CD4T_all.DOSE_result,
 # CD8T_all.DOSE_result, MacDC_all.DOSE_result, NK_all.DOSE_result)
#dplyr::count(Enrichment_All, Source, CellType, Class)
#write.csv(Enrichment_All, file = "../Manuscript/Figures/Figure 5-immune organoids and autoreactive B cells/Combined_Enrichment_Results.csv", row.names = FALSE)


#########
###to show that most cell types show autoimmuendisease pathway
#Macrophage:MacDC_all_DOSE_df
#CD8: CD8_up_kegg_df
#CD4: CD4_up_kegg_df
#NK: NK_all.DOSE_df
#B: B_up_kegg_df


library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)

# convert gene symbols to Entrez IDs
gene_list <- upregulated_genes$gene
entrez_ids <- bitr(gene_list, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

# Reactome enrichment
reactome_result <- enrichPathway(gene=entrez_ids$ENTREZID, pvalueCutoff=0.05, readable=TRUE)



#NK
autoimmune_terms_NK <- NK_all.DOSE_df[grep("autoimmune|immune system|COVID|virus|infectious", 
                                           NK_all.DOSE_df$Description, 
                                           ignore.case = TRUE), ]

#######################################################################################XC

#MacDC
#Define Your Autoimmune Disease Terms
autoimmune_terms <- c(
  "autoimmune disease", "lupus erythematosus", "systemic lupus erythematosus",
  "arthritis", "rheumatoid arthritis", "Sjogren's syndrome", "scleroderma",
  "systemic scleroderma", "multiple sclerosis", "autoimmune disease of blood",
  "hypersensitivity reaction type IV disease", "hypersensitivity reaction disease",
  "bacterial infectious disease","tuberculosis","hepatitis B","infectious"
)

#Subset Your Pathway Results
# Assume your data frame is called pathway_df
autoimmune_pathways_MacDC <- MacDC_all_DOSE_df[MacDC_all_DOSE_df$Description %in% autoimmune_terms, ]


#CD4T
# Your combined autoimmune-related and infection-triggered pathway terms
autoimmune_related_terms <- c(
  "Graft-versus-host disease",
  "Type I diabetes mellitus",
  "Systemic lupus erythematosus",
  "Autoimmune thyroid disease",
  "Rheumatoid arthritis",
  "Inflammatory bowel disease",
  "Th17 cell differentiation",
  "Th1 and Th2 cell differentiation",
  "IL-17 signaling pathway",
  "TNF signaling pathway",
  "Antigen processing and presentation",
  "Herpes simplex virus 1 infection",
  "Human T-cell leukemia virus 1 infection",
  "Influenza A",
  "Viral myocarditis",
  "Toxoplasmosis",
  "Staphylococcus aureus infection",
  "NOD-like receptor signaling pathway","Tuberculosis"
)

# Assuming your enrichment result is in a data frame called 'CD4_up_kegg_df'
# (e.g., from ReactomePA, clusterProfiler, etc.)
filtered_enrich_df <- CD4_up_kegg_df[CD4_up_kegg_df$Description %in% autoimmune_related_terms, ]

# View the filtered enrichment pathways
print(filtered_enrich_df)



#CD8T


# Step 1: Define the list of relevant disease terms
selected_pathways <- c(
  "Systemic lupus erythematosus",
  "Autoimmune thyroid disease",
  "Type I diabetes mellitus",
  "Rheumatoid arthritis",
  "Inflammatory bowel disease",
  "Epstein-Barr virus infection",
  "Viral myocarditis","Asthma","Th17 cell differentiation", "Th1 and Th2 cell differentiation"
)

# Step 3: Filter the enrichment results based on matching 'Description'
filtered_df_CD8 <- CD8_up_kegg_df[CD8_up_kegg_df$Description %in% selected_pathways, ]



#B cells

target_keywords <- c(
  "autoimmune disease",
  "systemic lupus erythematosus",
  "lupus erythematosus",
  "rheumatoid arthritis",
  "inflammatory bowel disease",
  "ulcerative colitis",
  "demyelinating disease",
  "viral infectious disease",
  "encephalitis",
  "disease by infectious agent"
)

filtered_df_B <- B_up_kegg_df[B_up_kegg_df$Description %in% target_keywords, ]


# add a column to identify the cell type
autoimmune_pathways_MacDC$CellType <- "MacDC"
filtered_enrich_df$CellType <- "CD4"
filtered_df_CD8$CellType <- "CD8"
autoimmune_terms_NK$CellType <- "NK"
filtered_df_B$CellType <- "B"

#make same columne table for each
MacDC_clean <- autoimmune_pathways_MacDC[, c("CellType", "Description", "Count","p.adjust")]
CD4_clean <- filtered_enrich_df[, c("CellType", "Description", "Count","p.adjust")]
CD8_clean <- filtered_df_CD8[, c("CellType", "Description", "Count","p.adjust")]
NK_clean <- autoimmune_terms_NK[, c("CellType", "Description", "Count","p.adjust")]
B_clean <- filtered_df_B[, c("CellType", "Description", "Count","p.adjust")]

#Normalize Description before rbind()
MacDC_clean$Description <- tolower(MacDC_clean$Description)
CD4_clean$Description   <- tolower(CD4_clean$Description)
CD8_clean$Description   <- tolower(CD8_clean$Description)
NK_clean$Description    <- tolower(NK_clean$Description)
B_clean$Description     <- tolower(B_clean$Description)

# Combine all into one dataframe
combined_df <- rbind(MacDC_clean, CD4_clean, CD8_clean, NK_clean, B_clean)

# seperate autoimmune disease, viral bacterial infection

# Extract autoimmune-related pathways
autoimmune_df <- combined_df[grep("ulcerative|type|systemic| syndrome|scleroderma|arthritis|sclerosis|lupus|inflammatory|hypersensitivity|graft-versus-host|encephalitis|autoimmune|arthritis", combined_df$Description, ignore.case = TRUE), ]
autoimmune_df$Pathway <- "Auto_disease"

# Extract viral infection–related pathways
viral_df <- combined_df[grep("viral|Epstein|herpes|HIV|hepatitis|papilloma|influenza|covid-19|virus|agent|coronavirus", combined_df$Description, ignore.case = TRUE), ]
viral_df$Pathway <- "viral"

# Extract bac infection–related pathways
bac_df <- combined_df[grep("tuberculosis|staphylococcus aureus infection|bacterial", combined_df$Description, ignore.case = TRUE), ]
bac_df$Pathway <- "bac"


#autoimmue related pathway
auto_process_df <- combined_df[grep("th17|th1|tnf|nod|signaling", combined_df$Description, ignore.case = TRUE), ]
auto_process_df$Pathway <- "auto_signaling"

#combined as row
library(dplyr)

combined_df <- dplyr::bind_rows(
  autoimmune_df %>% mutate(Source = "autoimmune_df"),
  auto_process_df %>% mutate(Source = "auto_process_df"),
  viral_df %>% mutate(Source = "viral_df"),
  bac_df %>% mutate(Source = "bac_df")
)

#Force Y-axis (Description) to follow its current row order
# Convert Description to factor with levels in the current row order
combined_df$Description <- factor(combined_df$Description, levels = rev(unique(combined_df$Description)))

# Then plot
plot <- ggplot(combined_df, aes(x = CellType, y = Description, size = Count, color = p.adjust)) +
  geom_point() +
  theme_bw() +
  labs(title = "Pathway Enrichment Across Cell Types") +
  theme(axis.text.y = element_text(size = 11))

# Print the plot
print(plot)


###
library(ggplot2)
library(ggtext)

# Assign colors to source categories
combined_df$LabelColor <- dplyr::case_when(
  combined_df$Source == "autoimmune_df"     ~ "#D73027",
  combined_df$Source == "auto_process_df"   ~ "#FC8D59",
  combined_df$Source == "viral_df"          ~ "#4575B4",
  combined_df$Source == "bac_df"            ~ "#1A9850",
  TRUE ~ "black"
)


# Assign colors to source categories
combined_df$LabelColor <- dplyr::case_when(
  combined_df$Source == "autoimmune_df"     ~ "#B2182B",
  combined_df$Source == "auto_process_df"   ~ "#D95F0E",
  combined_df$Source == "viral_df"          ~ "#4575B4",
  combined_df$Source == "bac_df"            ~ "#1A9850",
  TRUE ~ "black"
)

# Generate HTML-styled Description labels
combined_df$Description_colored <- paste0(
  "<span style='color:", combined_df$LabelColor, "'>",
  combined_df$Description, "</span>"
)

# Convert to factor to control order (if needed)
combined_df$Description_colored <- factor(
  combined_df$Description_colored,
  levels = rev(unique(combined_df$Description_colored))
)
plot <- ggplot(combined_df, aes(x = CellType, y = Description_colored, size = Count, color = p.adjust)) +
  geom_point() +
  theme_bw() +
  labs(title = "Pathway Enrichment Across Cell Types") +
  theme(
    axis.text.y = element_markdown(size = 11),  # <- this enables color
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size =14),
    axis.title.y = element_blank()
  ) +
  scale_color_gradient(low = "darkred", high = "pink", trans = "log10")


print(plot)

#install.packages("ggtext")  # Only run this once
library(ggtext)

####
ggsave(filename = "~/Library/CloudStorage/Dropbox/lyme_disease/Manuscript/Figures/Figure 5-immune organoids and autoreactive B cells/ Pathway Enrichment Across Cell Types.pdf", plot, width = 5, height = 7, dpi = 300)

