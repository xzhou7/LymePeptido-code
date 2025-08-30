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
library(DoubletFinder)
library(harmony)
library(ggpubr)
library(scales)
library(clusterProfiler)
library(org.Hs.eg.db)
library(EnhancedVolcano)

#mac
setwd("Library/CloudStorage/Dropbox/lyme_disease/Pepti_Spleen_singlecell/")

#windows
setwd("D:/Dropbox/lyme_disease/Pepti_Spleen_singlecell/")

getwd()

#Peptido_DEG <- read.csv("../Manuscript/Figures/Figure 5-immune organoids and autoreactive B cells/Spleen_DEG.csv", header = T, row.names = 1)
ALL_DEG <- read.csv("../Manuscript/Figures/Figure 5-immune organoids and autoreactive B cells/ALL_DEG.csv", header=T, row.names = 1)
ALL_DEG$gene <- row.names(ALL_DEG)

load("./RObject/03.DEGfinal.RData")

Patient_DEG <- read.csv("../PTLD_HD_Monocytes/HD_PTLD.csv", header = T)
Patient_DEG_sig <- filter(Patient_DEG, padj <= 0.05 )

Peptido_DEG <- spleen_DEG
#Is there a lot of overlapping genes between the patient/control and peptreat/control
macrophage_DEG <- filter(Peptido_DEG, celltype == "Macrophage_DC") %>% filter( p_val_adj < 0.0005) 
length(macrophage_DEG$gene)
intersect_genes <- intersect(Patient_DEG_sig$X, macrophage_DEG$gene)
length(intersect_genes)
only_in_patient <- setdiff(Patient_DEG_sig$X, macrophage_DEG$gene)
length(only_in_patient)
only_in_all <- setdiff(macrophage_DEG$gene, Patient_DEG_sig$X)
length(only_in_all)

# Make a named vector of colors based on overlap
macrophage_DEG$color <- ifelse(macrophage_DEG$gene %in% intersect_genes, "purple", "black")

p.macrophageDEG <- EnhancedVolcano(
  macrophage_DEG,
  lab = rep("", nrow(macrophage_DEG)),  # No text labels at all
  x = 'avg_log2FC',
  y = 'p_val_adj',
  FCcutoff = 0.15,
  xlim = c(-8, 8),
  title = "Macrophage DEG",
  colCustom = setNames(macrophage_DEG$color, macrophage_DEG$gene),
  legendPosition = "none" )

p.macrophageDEG

Patient_DEG_sig$color <- ifelse(Patient_DEG_sig$X %in% intersect_genes, "purple", "black")
p.PatientDEG <- EnhancedVolcano(
  Patient_DEG_sig,
  lab = rep("", nrow(Patient_DEG_sig)),  # suppress all labels
  x = 'log2FoldChange',
  y = 'padj',
  title = "Patient DEG",
  colCustom = setNames(Patient_DEG_sig$color, Patient_DEG_sig$X),
  legendPosition = "none"  # remove legend
)

p.PatientDEG

# Only label purple genes
Patient_DEG_sig_overlap <- filter(Patient_DEG_sig, color == "purple")

# Generate volcano plot
p.PatientDEG2 <- EnhancedVolcano(
  Patient_DEG_sig_overlap,
  lab = Patient_DEG_sig_overlap$label,
  x = 'log2FoldChange',
  y = 'padj',
  title = "Patient DEG",
  colCustom = setNames(Patient_DEG_sig_overlap$color, Patient_DEG_sig_overlap$X),
  legendPosition = "none"
)

p.PatientDEG2

DefaultAssay(spleen.combined.sct) <- "RNA"
universe.gene <- rownames(spleen.combined.sct)

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

#B cell
B_deg_results <-filter(Peptido_DEG_PATH, celltype == "B_Cell")

B_deg_results_unique <- B_deg_results %>%
  filter(!is.na(gene) & !is.na(avg_log2FC) & !is.na(p_val_adj)) %>%
  distinct(gene, .keep_all = TRUE)

P.B.deg <- EnhancedVolcano(B_deg_results_unique,
                           lab = B_deg_results_unique$gene,
                           FCcutoff=0.15,
                           xlim = c(-5.5, 5.5),
                           x = 'avg_log2FC',
                           y = 'p_val_adj', 
                           title = "B cell")
P.B.deg

#ggsave(filename = "./DEGandGO/Bcell_DEG.pdf",P.B.deg, width = 20, height = 15, dpi = 300)

#check a specific gene in all caterogys
filter(Peptido_DEG_PATH,SYMBOL == "CCL22")
filter(Peptido_DEG_PATH,SYMBOL == "IL17F")
filter(Peptido_DEG_PATH,SYMBOL == "IL1B")
filter(Peptido_DEG_PATH,SYMBOL == "IL6")
filter(Peptido_DEG_PATH,SYMBOL == "TNF")
filter(Peptido_DEG_PATH,SYMBOL == "CCL3")
filter(Peptido_DEG_PATH,SYMBOL == "CCL4")
filter(Peptido_DEG_PATH,SYMBOL == "CXCL8")
filter(Peptido_DEG_PATH,SYMBOL == "CXCL10")
filter(Peptido_DEG_PATH,SYMBOL == "CXCL5")
filter(Peptido_DEG_PATH,SYMBOL == "EBI3")

B_up <- filter(Peptido_DEG_PATH, celltype == "B_Cell" & avg_log2FC > 0 )#& p_val_adj < 0.005
B_down <- filter(Peptido_DEG_PATH, celltype == "B_Cell" & avg_log2FC < -0.5 )#& p_val_adj < 0.005

B_up.ego <- enrichGO(gene  = unique(B_up$ENTREZID),
                     OrgDb         = org.Hs.eg.db,
                     universe      = universe.gene,
                     ont           = "ALL",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05,
                     readable      = TRUE)

B_down.ego <- enrichGO(gene        = unique(B_down$ENTREZID),
                       OrgDb         = org.Hs.eg.db,
                       universe      = universe.gene,
                       ont           = "ALL",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05,
                       readable      = TRUE)

B_up.goresult.up <- B_up.ego@result
B_up.goresult.up$Class <- "UpRegulated"
Bdown.goresult.down <- B_down.ego@result
Bdown.goresult.down$Class <- "DownRegulated"

B_up.goresult.up$Description

B_goresult <- rbind(B_up.goresult.up, Bdown.goresult.down)
B_goresult$CellType <- "B_cell"

B_goresult

B_up.dotplot <- dotplot(B_up.ego, title = "B_up.goresult.up",showCategory = 40)
B_up.dotplot

B_down.dotplot <- dotplot(B_down.ego, title = "Bdown.goresult.down")

B_dotplot <- B_up.dotplot + B_down.dotplot
B_dotplot

B_up.ego_cnet_up <- cnetplot(
  B_up.ego,
  categorySize = "pvalue",
  showCategory = 50,
  foldChange = NULL,
  circular = FALSE,
  colorEdge = TRUE,
  node_label = "category" 
) + ggtitle("B Upregulated GO")

B_up.ego_cnet_up
#ggsave("./DEGandGO/BCell_GO_Up.pdf", B_up.ego_cnet_up, width = 15, height = 15, dpi = 300)

# KEGG enrichment for upregulated B cell genes
B_up.kegg <- enrichKEGG(
  gene          = unique(B_up$ENTREZID),
  organism      = "hsa", 
  universe      = universe.gene,
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2)
B_up.kegg

B_down.kegg <- enrichKEGG(
  gene          = unique(B_down$ENTREZID),
  organism      = "hsa",
  universe      = universe.gene,
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2)
B_down.kegg

B_up.kegg <- setReadable(B_up.kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
B_down.kegg <- setReadable(B_down.kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
B_up.kegg
B_down.kegg

B_up_kegg_plot <- dotplot(B_up.kegg, title = "B Cell KEGG - Upregulated", showCategory = 40)
B_up_kegg_plot

#ggsave("./DEGandGO/B_up_kegg_plot.pdf", B_up_kegg_plot,width = 5, height = 7, dpi = 300)

B_down_kegg_plot <- dotplot(B_down.kegg, title = "B Cell KEGG - Downregulated")
B_down_kegg_plot
B_kegg_dotplot <- B_up_kegg_plot + B_down_kegg_plot
B_kegg_dotplot

############################################
B_up.DOSE <- enrichDO(gene          = unique(B_up$ENTREZID)%>% na.omit(),
              ont           = "HDO",
              pvalueCutoff  = 0.05,
              pAdjustMethod = "BH",
              universe      = universe.gene,
              minGSSize     = 5,
              maxGSSize     = 500,
              qvalueCutoff  = 0.05,
              readable      = FALSE)
head(B_up.DOSE)

B_down.DOSE <- enrichDO(gene          = unique(B_down$ENTREZID)%>% na.omit(),
                      ont           = "HDO",
                      pvalueCutoff  = 0.05,
                      pAdjustMethod = "BH",
                      universe      = universe.gene,
                      minGSSize     = 5,
                      maxGSSize     = 500,
                      qvalueCutoff  = 0.05,
                      readable      = FALSE)


head(B_up.DOSE)
head(B_down.DOSE)
B_up_DOSE_plot <- dotplot(B_up.DOSE, title = "B Cell DOSE - Upregulated", showCategory = 40)
B_up_DOSE_plot

B_down_DOSE_plot <- dotplot(B_down.DOSE, title = "B Cell DOSE - Downregulated",showCategory = 25)
B_down_DOSE_plot

B_union_ENTREZID <- unique(c(B_up$ENTREZID, B_down$ENTREZID)) %>% na.omit()
B_all.DOSE <- enrichDO(gene          = B_union_ENTREZID,
                        ont           = "HDO",
                        pvalueCutoff  = 0.05,
                        pAdjustMethod = "BH",
                        universe      = universe.gene,
                        minGSSize     = 5,
                        maxGSSize     = 500,
                        qvalueCutoff  = 0.05,
                        readable      = FALSE)
B_all_DOSE_plot <- dotplot(B_all.DOSE, title = "B Cell DOSE - all", showCategory = 50)
B_all_DOSE_plot

B_dose.ego_cnet_up <- cnetplot(
  B_all.DOSE,
  categorySize = "pvalue",
  showCategory = 50,
  foldChange = NULL,
  circular = FALSE,
  colorEdge = TRUE,
  node_label = "category" 
) + ggtitle("B All DOSE")

B_dose.ego_cnet_up

#ggsave(filename = "./DEGandGO/DOSE_B_All.pdf", B_dose.ego_cnet_up, width = 20, height = 15, dpi=300)

# Subset up- and down-regulated genes for CD4T_Cell##################################################################################
CD4_deg_results <-filter(Peptido_DEG_PATH, celltype == "CD4T_Cell")

CD4_deg_results_unique <- CD4_deg_results %>%
  filter(!is.na(gene) & !is.na(avg_log2FC) & !is.na(p_val_adj)) %>%
  distinct(gene, .keep_all = TRUE)

P.CD4.deg <- EnhancedVolcano(CD4_deg_results_unique,
                           lab = CD4_deg_results_unique$gene,
                           FCcutoff=0.15,
                           xlim = c(-5.5, 5.5),
                           x = 'avg_log2FC',
                           y = 'p_val_adj', 
                           title = "CD4T_Cell")
P.CD4.deg

#ggsave(filename = "./DEGandGO/CD4Tcell_DEG.pdf",P.CD4.deg, width = 20, height = 15, dpi = 300)

CD4T_up <- filter(Peptido_DEG_PATH, celltype == "CD4T_Cell" & avg_log2FC > 0 & p_val_adj < 0.05)
CD4T_down <- filter(Peptido_DEG_PATH, celltype == "CD4T_Cell" & avg_log2FC < 0& p_val_adj < 0.05)

CD4T_up.ego <- enrichGO(
  gene          = unique(CD4T_up$ENTREZID),
  OrgDb         = org.Hs.eg.db,
  universe      = universe.gene,
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05,
  readable      = TRUE)

CD4T_down.ego <- enrichGO(
  gene          = unique(CD4T_down$ENTREZID),
  OrgDb         = org.Hs.eg.db,
  universe      = universe.gene,
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05,
  readable      = TRUE)

CD4T_up.goresult <- CD4T_up.ego@result
CD4T_up.goresult$Class <- "UpRegulated"
CD4T_down.goresult <- CD4T_down.ego@result
CD4T_down.goresult$Class <- "DownRegulated"

CD4T_goresult <- rbind(CD4T_up.goresult, CD4T_down.goresult)
CD4T_goresult$CellType <- "CD4T_Cell"

CD4T_up.dotplot <- dotplot(CD4T_up.ego, title = "CD4T Upregulated GO")
CD4T_down.dotplot <- dotplot(CD4T_down.ego, title = "CD4T Downregulated GO")
CD4T_dotplot <- CD4T_up.dotplot + CD4T_down.dotplot
CD4T_dotplot

cnet_up <- cnetplot(
  CD4T_up.ego,
  categorySize = "pvalue",
  showCategory = 50,
  foldChange = NULL,
  circular = FALSE,
  colorEdge = TRUE,
  node_label = "category" 
) + ggtitle("CD4 T Upregulated GO (BP)")

cnet_up

#ggsave(filename = "./DEGandGO/CD4_T_GO_Up_CNET.pdf", cnet_up, width = 20, height = 15, dpi = 300)

cnet_down <- cnetplot(
  CD4T_down.ego,
  categorySize = "pvalue",
  showCategory = 50,
  foldChange = NULL,
  circular = FALSE,
  colorEdge = TRUE,
  node_label = "category" 
) + ggtitle("CD4 T Downregulated GO (BP)")

cnet_down
#ggsave(filename = "./DEGandGO/CD4_T_GO_Down_CNET.pdf", cnet_down, width = 20, height = 15, dpi = 300)

CD4_up.kegg <- enrichKEGG(
  gene          = unique(CD4T_up$ENTREZID),
  organism      = "hsa",
  universe      = universe.gene,
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2)

CD4_down.kegg <- enrichKEGG(
  gene          = unique(CD4T_down$ENTREZID),
  organism      = "hsa",
  universe      = universe.gene,
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.2)

CD4_up.kegg <- setReadable(CD4_up.kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
CD4_down.kegg <- setReadable(CD4_down.kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

CD4_up_kegg_plot <- dotplot(CD4_up.kegg, title = "CD4 T Cell KEGG - Upregulated", showCategory = 40)
CD4_down_kegg_plot <- dotplot(CD4_down.kegg, title = "CD4 T Cell KEGG - Downregulated")

CD4_kegg_dotplot <- CD4_up_kegg_plot + CD4_down_kegg_plot
CD4_kegg_dotplot

#ggsave("./DEGandGO/CD4_kegg_dotplot.pdf", CD4_kegg_dotplot, width = 10, height = 7, dpi = 300)

#DOSE
CD4T_up.DOSE <- enrichDO(
  gene          = unique(CD4T_up$ENTREZID),
  ont           = "HDO",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  universe      = universe.gene,   # optional background
  minGSSize     = 5,
  maxGSSize     = 500,
  qvalueCutoff  = 0.05,
  readable      = T)

CD4T_down.DOSE <- enrichDO(
  gene          = unique(CD4T_down$ENTREZID),
  ont           = "HDO",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  universe      = universe.gene,
  minGSSize     = 5,
  maxGSSize     = 500,
  qvalueCutoff  = 0.05,
  readable      = T)

head(CD4T_up.DOSE)
head(CD4T_down.DOSE)

CD4T_up_DOSE_plot <- dotplot(CD4T_up.DOSE, title = "CD4 T Cell DOSE - Upregulated", showCategory = 25)
CD4T_up_DOSE_plot

CD4T_down_DOSE_plot <- dotplot(CD4T_down.DOSE, title = "CD4 T Cell DOSE - Downregulated",showCategory = 25)
CD4T_down_DOSE_plot

#unified DOSE
CD4_union_ENTREZID <- unique(c(CD4T_up$ENTREZID, CD4T_down$ENTREZID)) %>% na.omit()

CD4_all.DOSE <- enrichDO(
  gene          = CD4_union_ENTREZID,
  ont           = "HDO",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  universe      = universe.gene,
  minGSSize     = 5,
  maxGSSize     = 500,
  qvalueCutoff  = 0.05,
  readable      = FALSE)

CD4_all_DOSE_plot <- dotplot(CD4_all.DOSE, title = "CD4 T Cell DOSE - all", showCategory = 50)
CD4_all_DOSE_plot

CD4_dose.ego_cnet <- cnetplot(
  CD4_all.DOSE,
  categorySize = "pvalue",
  showCategory = 50,
  foldChange = NULL,
  circular = FALSE,
  colorEdge = TRUE,
  node_label = "category"
) + ggtitle("CD4 T Cell DOSE")

CD4_dose.ego_cnet

ggsave(filename = "./DEGandGO/DOSE_CD4_All.pdf", CD4_dose.ego_cnet, width = 20, height = 15, dpi = 300)

# Subset up- and down-regulated genes for CD8T_Cell
CD8_deg_results <-filter(Peptido_DEG_PATH, celltype == "CD8T_Cell")
CD8_deg_results_unique <- CD8_deg_results %>%
  filter(!is.na(gene) & !is.na(avg_log2FC) & !is.na(p_val_adj)) %>%
  distinct(gene, .keep_all = TRUE)

P.CD8.deg <- EnhancedVolcano(CD8_deg_results,
                             lab = CD8_deg_results$gene,
                             FCcutoff=0.15,
                             xlim = c(-5.5, 5.5),
                             x = 'avg_log2FC',
                             y = 'p_val_adj', 
                             title = "CD8T_Cell")
P.CD8.deg

#ggsave(filename = "./DEGandGO/CD8Tcell_DEG.pdf",P.CD8.deg, width = 20, height = 15, dpi = 300)

CD8T_up <- filter(Peptido_DEG_PATH, celltype == "CD8T_Cell" & avg_log2FC > 0 & p_val_adj < 0.05)
CD8T_down <- filter(Peptido_DEG_PATH, celltype == "CD8T_Cell" & avg_log2FC < 0 & p_val_adj < 0.05)

CD8T_up.ego <- enrichGO(
  gene          = unique(CD8T_up$ENTREZID),
  OrgDb         = org.Hs.eg.db,
  universe      = universe.gene,
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05,
  readable      = TRUE)

CD8T_down.ego <- enrichGO(
  gene          = CD8T_down$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  universe      = universe.gene,
  ont           = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05,
  readable      = TRUE)

CD8T_up.goresult <- CD8T_up.ego@result
CD8T_up.goresult$Class <- "UpRegulated"
CD8T_down.goresult <- CD8T_down.ego@result
CD8T_down.goresult$Class <- "DownRegulated"

CD8T_goresult <- rbind(CD8T_up.goresult, CD8T_down.goresult)
CD8T_goresult$CellType <- "CD8T_Cell"

CD8T_up.dotplot <- dotplot(CD8T_up.ego, title = "CD8T Upregulated GO")
CD8T_down.dotplot <- dotplot(CD8T_down.ego, title = "CD8T Downregulated GO")
CD8T_dotplot <- CD8T_up.dotplot + CD8T_down.dotplot
CD8T_dotplot

cnet_up_CD8 <- cnetplot(
  CD8T_up.ego,
  categorySize = "pvalue",
  showCategory = 50,
  foldChange = NULL,
  circular = FALSE,
  colorEdge = TRUE,
  node_label = "category"
) + ggtitle("CD8 T Upregulated GO (BP)")

cnet_up_CD8

ggsave(filename = "./DEGandGO/CD8_T_GO_Up_CNET.pdf", cnet_up_CD8, width = 20, height = 15, dpi = 300)

cnet_down_CD8 <- cnetplot(
  CD8T_down.ego,
  categorySize = "pvalue",
  showCategory = 50,
  foldChange = NULL,
  circular = FALSE,
  colorEdge = TRUE,
  node_label = "category"
) + ggtitle("CD8 T Downregulated GO (BP)")

cnet_down_CD8

ggsave(filename = "./DEGandGO/CD8_T_GO_Down_CNET.pdf", cnet_down_CD8, width = 20, height = 15, dpi = 300)

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

CD8T_up.kegg.plot <- dotplot(CD8T_up.kegg, title = "CD8T Upregulated KEGG", showCategory = 50)
CD8T_down.kegg.plot <- dotplot(CD8T_down.kegg, title = "CD8T Downregulated KEGG",showCategory = 50)

CD8T_kegg_dotplot <- CD8T_up.kegg.plot + CD8T_down.kegg.plot
CD8T_kegg_dotplot

#ggsave(filename = "DEGandGO/CD8_KEGG_UpDown.pdf",CD8T_kegg_dotplot, width = 10, height = 10, dpi = 300)

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

CD8_dose.ego_cnet <- cnetplot(
  CD8_all.DOSE,
  categorySize = "pvalue",
  showCategory = 50,
  foldChange = NULL,
  circular = FALSE,
  colorEdge = TRUE,
  node_label = "category") + ggtitle("CD8 T Cell DOSE")

CD8_dose.ego_cnet

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

cnet_up_MacDC <- cnetplot(
  MacDC_up.ego,
  categorySize = "pvalue",
  showCategory = 50,
  foldChange = NULL,
  circular = FALSE,
  colorEdge = TRUE,
  node_label = "category"
) + ggtitle("Macrophage/DC Upregulated GO (BP)")

cnet_up_MacDC

ggsave(filename = "./DEGandGO/MacDC_GO_Up_CNET.pdf", cnet_up_MacDC, width = 20, height = 15, dpi = 300)

cnet_down_MacDC <- cnetplot(
  MacDC_down.ego,
  categorySize = "pvalue",
  showCategory = 50,
  foldChange = NULL,
  circular = FALSE,
  colorEdge = TRUE,
  node_label = "category"
) + ggtitle("Macrophage/DC Downregulated GO (BP)")

cnet_down_MacDC

ggsave(filename = "./DEGandGO/MacDC_GO_Down_CNET.pdf", cnet_down_MacDC, width = 20, height = 15, dpi = 300)

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

cnet_up <- cnetplot(
  MacDC_up_BP,
  categorySize = "pvalue",
  showCategory = 50,
  foldChange = NULL,
  circular = FALSE,
  colorEdge = TRUE,
  node_label = "category"  # Only show pathway labels
) + ggtitle("Macrophage_DC Upregulated GO (BP)")

cnet_up

cnet_down <- cnetplot(
  MacDC_down_BP,
  categorySize = "pvalue",
  showCategory = 50,
  foldChange = NULL,
  circular = FALSE,
  colorEdge = TRUE,
  node_label = "category"  # Only show pathway labels
) + ggtitle("Macrophage_DC Downregulated GO (BP)")

cnet_down

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

MacDC_dose.ego_cnet <- cnetplot(
  MacDC_all.DOSE,
  categorySize = "pvalue",
  showCategory = 50,
  foldChange = NULL,
  circular = FALSE,
  colorEdge = TRUE,
  node_label = "category") + ggtitle("Macrophage/DC DOSE")

MacDC_dose.ego_cnet

ggsave(filename = "./DEGandGO/DOSE_MacDC_All.pdf", MacDC_dose.ego_cnet, width = 20, height = 15, dpi = 300)

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

cnet_up_NK <- cnetplot(
  NK_up.ego,
  categorySize = "pvalue",
  showCategory = 50,
  foldChange = NULL,
  circular = FALSE,
  colorEdge = TRUE,
  node_label = "category"
) + ggtitle("NK Cell Upregulated GO (BP)")

cnet_up_NK

ggsave(filename = "./DEGandGO/NK_GO_Up_CNET.pdf", cnet_up_NK, width = 20, height = 15, dpi = 300)

cnet_down_NK <- cnetplot(
  NK_down.ego,
  categorySize = "pvalue",
  showCategory = 50,
  foldChange = NULL,
  circular = FALSE,
  colorEdge = TRUE,
  node_label = "category"
) + ggtitle("NK Cell Downregulated GO (BP)")

cnet_down_NK

ggsave(filename = "./DEGandGO/NK_GO_Down_CNET.pdf", cnet_down_NK, width = 20, height = 15, dpi = 300)

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
  ont           = "HDO",
  pvalueCutoff  = 0.05,
  pAdjustMethod = "BH",
  universe      = universe.gene,
  minGSSize     = 5,
  maxGSSize     = 500,
  qvalueCutoff  = 0.05,
  readable      = TRUE)

NK_down.DOSE <- enrichDO(
  gene          = unique(NK_down$ENTREZID),
  ont           = "HDO",
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

NK_dose.ego_cnet <- cnetplot(
  NK_all.DOSE,
  categorySize = "pvalue",
  showCategory = 50,
  foldChange = NULL,
  circular = FALSE,
  colorEdge = TRUE,
  node_label = "category"
) + ggtitle("NK Cell DOSE")

NK_dose.ego_cnet

ggsave(filename = "./DEGandGO/DOSE_NK_All.pdf", NK_dose.ego_cnet, width = 20, height = 15, dpi = 300)


# Add cell type and source annotation
B_goresult$CellType <- "B_cell"        
B_goresult$Source <- "GO"
CD4T_goresult$CellType <- "CD4_T"      
CD4T_goresult$Source <- "GO"
CD8T_goresult$CellType <- "CD8_T"
CD8T_goresult$Source <- "GO"
MacDC_goresult$CellType <- "Macrophage_DC"
MacDC_goresult$Source <- "GO"
NK_goresult$CellType <- "NK"
NK_goresult$Source <- "GO"

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

GO_all <- bind_rows(B_goresult, CD4T_goresult, CD8T_goresult,MacDC_goresult,NK_goresult)
GO_all
write.csv(GO_all, file = "../Manuscript/Figures/Figure 5-immune organoids and autoreactive B cells/Combined_GO_Results.csv", row.names = FALSE)

Enrichment_All <- bind_rows(
  B_goresult, CD4T_goresult, CD8T_goresult, MacDC_goresult, NK_goresult,
  B_up.kegg_result, B_down.kegg_result,
  CD4T_up.kegg_result, CD4T_down.kegg_result,
  CD8T_up.kegg_result, CD8T_down.kegg_result,
  MacDC_up.kegg_result, MacDC_down.kegg_result,
  NK_up.kegg_result, NK_down.kegg_result,
  B_all.DOSE_result, CD4T_all.DOSE_result,
  CD8T_all.DOSE_result, MacDC_all.DOSE_result, NK_all.DOSE_result)
dplyr::count(Enrichment_All, Source, CellType, Class)
write.csv(Enrichment_All, file = "../Manuscript/Figures/Figure 5-immune organoids and autoreactive B cells/Combined_Enrichment_Results.csv", row.names = FALSE)

#########################
#for hervK
DefaultAssay(spleen.combined.sct) <- "RNA"
expr_counts <- FetchData(spleen.combined.sct, vars = "ERVK3-1")
sum(expr_counts$ERVK3-1 > 0)

expr_df <- FetchData(spleen.combined.sct, vars = c("ERVK3-1", "celltype2"))
expr_df
expr_summary <- expr_df %>%
  group_by(celltype2) %>%
  summarise(
    total_cells = n(),
    expressing_cells = sum(`ERVK3-1` > 0),
    percent_expressing = round(100 * expressing_cells / total_cells, 2)
  )

expr_summary

# Fetch required data
ervk3_data <- FetchData(spleen.combined.sct, vars = c("ERVK3-1", "celltype2", "orig.ident"))

# Filter for expressing cells
ervk3_positive_counts <- ervk3_data %>%
  filter(`ERVK3-1` > 0) %>%
  group_by(celltype2, orig.ident) %>%
  summarise(expressing_cell_count = n(), .groups = "drop")

ervk3_positive_counts

# Fetch required data
ervk3_data <- FetchData(spleen.combined.sct, vars = c("ERVK3-1", "celltype2", "orig.ident"))

# Fetch required data
ervk3_data <- FetchData(spleen.combined.sct, vars = c("ERVK3-1", "celltype2", "orig.ident"))

# Calculate total and positive counts
ervk3_percentage <- ervk3_data %>%
  group_by(celltype2, orig.ident) %>%
  summarise(
    total_cells = n(),
    expressing_cells = sum(`ERVK3-1` > 0),
    percent_expressing = 100 * expressing_cells / total_cells,
    .groups = "drop")

ervk3_percentage


feature_to_plot <- "your_numeric_feature"

df.percent <- FetchData(spleen.combined.sct, vars = c(feature_to_plot, "celltype2")) %>%
  filter(!is.na(.data[[feature_to_plot]]), !is.na(celltype2)) %>%
  group_by(celltype2) %>%
  mutate(median_val = median(.data[[feature_to_plot]])) %>%
  ungroup() %>%
  mutate(celltype2 = fct_reorder(celltype2, median_val, .desc = TRUE))


ggplot(df.percent, aes(x = celltype2, y = .data[[feature_to_plot]])) +
  geom_boxplot(fill = "lightblue", outlier.size = 0.5) +
  theme_minimal() +
  labs(x = "Cell Type", y = feature_to_plot, title = paste(feature_to_plot, "Across Cell Types")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(0, NA)
