##Healthy Donor and PTLD Donor Monocytes RNAseq
#Updated Mar 11, 2025
#Author: Xin Zhou, Ph.D. 

library(DESeq2)
library(ggplot2)
library(clusterProfiler)
library(dplyr)
library(factoextra)
library(stringr)
library(ReactomePA)
library(reactome.db)
library(VennDiagram)
library(tidyverse)
library(ggsci)
library(org.Hs.eg.db)
library(DOSE)
library(EnhancedVolcano)
library(enrichplot)
library(apeglm)

#Mac  
#setwd("~/Library/CloudStorage/Box-Box/Xin.Chen.Shareable/R5_Lyme_Paper/Figure.2_Monocytes/Data/")

#Windows
#setwd("C:/Users/zhoux/Box/Xin.Chen.Shareable/PTLD_HD_Monocytes/Data/")
setwd("D:/Dropbox/lyme_disease/PTLD_HD_Monocytes/Data/")

gene.count <- read.csv("./Bulk_seq_momo_Count/aligned_to_hg19_-562735560/RNA_Express_AppResult-ds.72563998561b473690d3de501c4f12b1/differential/global/gene.counts.csv",header = T, row.names = 1)
meta.table <- read.csv("./metadata.csv", header = T, row.names = 1)
head(gene.count)

meta.table

#check sequencing depth
colSums(gene.count)
#check if any gene is empty
table(rowSums(gene.count) == 0)
gene.count.clean <- gene.count[!rowSums(gene.count)==0,]

#PCA plot
gene.count.clean.t <- as.data.frame(t(gene.count.clean))
gene.count.clean.t[1:5, 1:5]

res.pca <- prcomp(gene.count.clean.t, scale = TRUE)
fviz_eig(res.pca)

####if want label back, delete (label = "none")
p1 <- fviz_pca_ind(res.pca, habillage=meta.table$dex, label="none")
p1 <- p1 + theme_minimal() + scale_shape_manual(values=c(15,16,17,18)) + scale_color_aaas()+ 
  theme(legend.position = "bottom")
p1

p2 <- fviz_pca_ind(res.pca, habillage=meta.table$celltype,label = "none")
p2 <- p2 + theme_minimal() + scale_shape_manual(values=c(15,16,17,18)) + scale_color_d3()+ 
  theme(legend.position = "bottom")
p2

p3 <- fviz_pca_ind(res.pca, habillage=meta.table$status,label = "none")
p3 <- p3 + theme_minimal() + scale_shape_manual(values=c(15,16,17,18)) + scale_color_nejm()+ 
  theme(legend.position = "bottom")
p3 

p4 <- fviz_pca_ind(res.pca, habillage=meta.table$condition,label = "none")
p4 <- p4 + theme_minimal() + scale_shape_manual(values=c(15,16,17,18,19,20)) + scale_color_jama() + 
  theme(legend.position = "bottom")
p4 

p.pca <- p1+p2+p3+p4
p.pca

#######################
#seperate CD14 and CD16
CD14_meta <- filter(meta.table, celltype == "CD14" & status == "Rest")
rownames(CD14_meta) <- paste("X",str_replace(row.names(CD14_meta), "-", ".") %>% str_replace("-", ".")  %>% str_replace("-", "."), sep="")
CD14_cts <- dplyr::select(gene.count, row.names(CD14_meta))
CD14_cts_clean <- CD14_cts[!rowSums(CD14_cts)==0,]
dim(CD14_cts_clean)
dim(CD14_meta)

CD14_act <- filter(meta.table, celltype == "CD14" & status == "Act")
rownames(CD14_act) <- paste("X",str_replace(row.names(CD14_act), "-", ".") %>% str_replace("-", ".")  %>% str_replace("-", "."), sep="")
CD14_act_cts <- dplyr::select(gene.count, row.names(CD14_act))

CD14_act_cts_clean <- CD14_act_cts[!rowSums(CD14_act_cts)==0,]
dim(CD14_act_cts_clean)
dim(CD14_act)
CD14_act$treatment <- c("NLPS","NLPS","LPS", "NLPS","LPS", "LPS")

cd16_meta <- filter(meta.table, celltype == "CD16")
rownames(cd16_meta) <- paste("X",str_replace(row.names(cd16_meta), "-", ".") %>% str_replace("-", ".")  %>% str_replace("-", "."), sep="")
cd16_cts <- dplyr::select(gene.count, row.names(cd16_meta))
cd16_cts_clean <- cd16_cts[!rowSums(cd16_cts)==0,]
dim(cd16_cts_clean)
dim(cd16_meta)

#####################################################################################
#set up deseq2 object
colnames(gene.count)
dds_CD14_rest <- DESeqDataSetFromMatrix(countData = CD14_cts_clean,
                              colData = CD14_meta,
                              design= ~ dex)

dds_CD14_rest <- DESeq(dds_CD14_rest)
resultsNames(dds_CD14_rest)
res_CD14_rest <- results(dds_CD14_rest, name="dex_PTLD_vs_HD")
res_CD14_rest2 <- lfcShrink(dds_CD14_rest, coef="dex_PTLD_vs_HD", type="apeglm")

plotMA(res_CD14_rest, ylim=c(-2,2))
plotMA(res_CD14_rest2, ylim=c(-2,2))

res_CD14_rest

p.CD14 <-EnhancedVolcano(res_CD14_rest,
                lab = rownames(res_CD14_rest),
                x = 'log2FoldChange',
                y = 'pvalue')

p.CD14 <- p.CD14 + ggtitle(" CD14 monocytes Control/Lyme")
p.CD14

#####Figure 3D
p.CD14 <-EnhancedVolcano(res_CD14_rest,
                         lab = rownames(res_CD14_rest),
                         x = 'log2FoldChange',
                         selectLab = c('TMEM176B','ADH1A'),
                         y = 'pvalue',
                         pCutoff = 10e-8,
                         FCcutoff = 1.5,
                         col=c('black', 'blue', 'blue', 'red3'),
                         drawConnectors = TRUE,
                         widthConnectors = 0.75)

p.CD14 <- p.CD14 + ggtitle(" CD14 monocytes Control/Lyme")
p.CD14
#ggsave(filename = "../../lyme_disease/Manuscript/Figures/Figure3/monocytesDEG.pdf", p.CD14, width = 6, height = 7, dpi = 300)

res_CD14_rest[rownames(res_CD14_rest) == "CD14",]

cd14_result.df <- res_CD14_rest@listData %>% data.frame()
rownames(cd14_result.df) <- res_CD14_rest@rownames
cd14_result.df$gene <- res_CD14_rest@rownames
#reorder table based on p value
cd14_result.df <- cd14_result.df %>% arrange(padj)
#write.csv(file = "./Bulk_seq_momo_Count/DESEQ/CD14rest.DESEQ2.csv",cd14_result.df)
#write.csv(file="../../lyme_disease/Manuscript/Figures/Figure3/monocytesDEG_HD_PTLD.csv",cd14_result.df)
#pathway analysis

genelist <- rownames(gene.count)
universe.gene.list <- bitr(genelist, fromType = "SYMBOL",
                             toType = "ENTREZID",
                             OrgDb = org.Hs.eg.db)

cd14.rest.up <- filter(cd14_result.df, padj < 0.05 & log2FoldChange > 0)
cd14.rest.down <- filter(cd14_result.df, padj < 0.05 & log2FoldChange < 0)

cd14.rest.up.df <- bitr(cd14.rest.up$gene, fromType = "SYMBOL",
                toType = "ENTREZID",
                OrgDb = org.Hs.eg.db)

cd14.rest.down.df <- bitr(cd14.rest.down$gene, fromType = "SYMBOL",
                        toType = "ENTREZID",
                        OrgDb = org.Hs.eg.db)

#https://supfam.mrc-lmb.cam.ac.uk/SUPERFAMILY/cgi-bin/go.cgi
ego.1.up <- enrichGO(gene       = cd14.rest.up.df$ENTREZID,
                universe      = universe.gene.list$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05)

summary(ego.1.up)
ego.1.up <- setReadable(ego.1.up, OrgDb = org.Hs.eg.db)
#write.csv(file = "./Bulk_seq_momo_Count/DESEQ/GO_CD14rest.up.csv",ego.1.up@result)

dotplot(ego.1.up, showCategory=20)
p.GO <- enrichplot::goplot(ego.1.up,showCategory=20)
p.GO
#ggsave(filename = "./Bulk_seq_momo_Count/CD14Lyme_GO.pdf", p.GO, width = 15, height = 10, dpi=300)

ego.1.up2 <- pairwise_termsim(ego.1.up)
ego.1.up2.p1 <- treeplot(ego.1.up2)
ego.1.up2.p1
ego.1.up2.p2 <- treeplot(ego.1.up2, hclust_method = "average")
ego.1.up2.p2 <- ego.1.up2.p2 + scale_color_gradient(low = "pink", high = "darkred")
ego.1.up2.p2

ggsave(filename = "../../Manuscript/Figures/Figure3-autoimmune monocytes/UpregulatedGO.pdf",ego.1.up2.p2, width = 7, height = 5, dpi=300 )
#write.csv(file = "../../lyme_disease/Manuscript/Figures/Figure3/UpregulatedGO.csv",ego.1.up)

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
ego.1.down2.p2 <- ego.1.down2.p2 + scale_color_gradient(low = "pink", high = "darkred")
ego.1.down2.p2

ggsave(filename = "../../Manuscript/Figures/Figure3-autoimmune monocytes/DownregulatedGO.pdf",ego.1.down2.p2, width = 7, height = 5, dpi=300 )
write.csv(file = "../../lyme_disease/Manuscript/Figures/Figure3/DownregulatedGO.csv",ego.1.down)

#write.csv(file = "./Bulk_seq_momo_Count/DESEQ/GO_CD14rest.down.csv",ego.1.down@result)

dotplot(ego.1.down, showCategory=40)
p.GO2 <- goplot(ego.1.down,showCategory=20)
p.GO2
#####################################################################################
CD14_act$treatment <- factor(CD14_act$treatment, levels = c("NLPS", "LPS"))

dds_CD14_LPS <- DESeqDataSetFromMatrix(countData = CD14_act_cts_clean,
                                        colData = CD14_act,
                                        design= ~ treatment)

dds_CD14_LPS <- DESeq(dds_CD14_LPS)
resultsNames(dds_CD14_LPS)
CD14_LPS <- results(dds_CD14_LPS, name="treatment_LPS_vs_NLPS")
CD14_LPS2 <- lfcShrink(dds_CD14_LPS, coef="treatment_LPS_vs_NLPS", type="apeglm")

plotMA(CD14_LPS, ylim=c(-2,2))
plotMA(CD14_LPS2, ylim=c(-2,2))

lps.CD14 <-EnhancedVolcano(CD14_LPS,
                         lab = rownames(CD14_LPS),
                         x = 'log2FoldChange',
                         y = 'pvalue')

lps.CD14 <- lps.CD14 + ggtitle("CD14 monocytes Control/LPS treatment")
lps.CD14

CD14_LPS

CD14_LPS[rownames(CD14_LPS) == "CCL4",]
CD14_LPS[rownames(CD14_LPS) == "CD14",]


lps.cd14_result.df <- CD14_LPS@listData %>% data.frame()
dim(lps.cd14_result.df)
dim(CD14_LPS)

rownames(lps.cd14_result.df) <- CD14_LPS@rownames
lps.cd14_result.df$gene <- CD14_LPS@rownames
#reorder table based on p value
lps.cd14_result.df <- lps.cd14_result.df %>% arrange(padj)
#write.csv(file = "./Bulk_seq_momo_Count/DESEQ/CD14_LPS_NLPS.DESEQ2.csv",lps.cd14_result.df)

#pathway analysis
lps.cd14.rest.up <- filter(lps.cd14_result.df, padj < 0.05 & log2FoldChange > 0)
lps.cd14.rest.down <- filter(lps.cd14_result.df, padj < 0.05 & log2FoldChange < 0)

lps.cd14.up.df <- bitr(lps.cd14.rest.up$gene, fromType = "SYMBOL",
                        toType = "ENTREZID",
                        OrgDb = org.Hs.eg.db)

lps.cd14.down.df <- bitr(lps.cd14.rest.down$gene, fromType = "SYMBOL",
                          toType = "ENTREZID",
                          OrgDb = org.Hs.eg.db)

#https://supfam.mrc-lmb.cam.ac.uk/SUPERFAMILY/cgi-bin/go.cgi
lps.ego.1.up <- enrichGO(gene      = lps.cd14.up.df$ENTREZID,
                     universe      = universe.gene.list$ENTREZID,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05)

(summary(lps.ego.1.up))
lps.ego.1.up <- setReadable(lps.ego.1.up, OrgDb = org.Hs.eg.db)
lps.ego.1.up@result
#write.csv(file = "./Bulk_seq_momo_Count/DESEQ/GO_CD14_LPS.up.csv",lps.ego.1.up@result)

dotplot(lps.ego.1.up, showCategory=20)
enrichMap(lps.ego.1.up, vertex.label.cex=1.2, layout=igraph::layout.kamada.kawai)

lps.ego.1.down <- enrichGO(gene      = lps.cd14.down.df$ENTREZID,
                       universe      = universe.gene.list$ENTREZID,
                       OrgDb         = org.Hs.eg.db,
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05)

summary(lps.ego.1.down)
lps.ego.1.down <- setReadable(lps.ego.1.down, OrgDb = org.Hs.eg.db)
#write.csv(file = "./Bulk_seq_momo_Count/DESEQ/GO_CD14_LPS.down.csv",lps.ego.1.down@result)

dotplot(ego.1.down, showCategory=30)

#####################################################################################

dds_CD16_rest <- DESeqDataSetFromMatrix(countData = cd16_cts_clean,
                                        colData = cd16_meta,
                                        design= ~ dex)

dds_CD16_rest <- DESeq(dds_CD16_rest)
resultsNames(dds_CD16_rest)
res_CD16_rest <- results(dds_CD16_rest, name="dex_PTLD_vs_HD")
res_CD16_rest2 <- lfcShrink(dds_CD16_rest, coef="dex_PTLD_vs_HD", type="apeglm")

plotMA(res_CD16_rest, ylim=c(-2,2))
plotMA(res_CD16_rest2, ylim=c(-2,2))

res_CD16_rest

p.CD16 <-EnhancedVolcano(res_CD16_rest,
                         lab = rownames(res_CD16_rest),
                         x = 'log2FoldChange',
                         y = 'pvalue')

p.CD16 <- p.CD16 + ggtitle(" CD16 monocytes Control/Lyme")
p.CD16

CD16_result.df <- res_CD16_rest@listData %>% data.frame()
rownames(CD16_result.df) <- res_CD16_rest@rownames
CD16_result.df$gene <- res_CD16_rest@rownames
#reorder table based on p value
CD16_result.df <- CD16_result.df %>% arrange(padj)
CD16_result.df
#write.csv(file = "./Bulk_seq_momo_Count/DESEQ/CD16rest.DESEQ2.csv",CD16_result.df)

#pathway analysis
CD16.rest.up <- filter(CD16_result.df, padj < 0.05 & log2FoldChange > 0)
CD16.rest.down <- filter(CD16_result.df, padj < 0.05 & log2FoldChange < 0)

CD16.rest.up.df <- bitr(CD16.rest.up$gene, fromType = "SYMBOL",
                        toType = "ENTREZID",
                        OrgDb = org.Hs.eg.db)

CD16.rest.down.df <- bitr(CD16.rest.down$gene, fromType = "SYMBOL",
                          toType = "ENTREZID",
                          OrgDb = org.Hs.eg.db)

#https://supfam.mrc-lmb.cam.ac.uk/SUPERFAMILY/cgi-bin/go.cgi
CD16_ego.1.up <- enrichGO(gene       = CD16.rest.up.df$ENTREZID,
                     universe      = universe.gene.list$ENTREZID,
                     OrgDb         = org.Hs.eg.db,
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.01,
                     qvalueCutoff  = 0.05)

(summary(CD16_ego.1.up))
CD16_ego.1.up <- setReadable(CD16_ego.1.up, OrgDb = org.Hs.eg.db)
#write.csv(file = "./Bulk_seq_momo_Count/DESEQ/GO_CD16rest.up.csv",CD16_ego.1.up@result)

dotplot(CD16_ego.1.up, showCategory=30)
enrichMap(CD16_ego.1.up, vertex.label.cex=1.2, layout=igraph::layout.kamada.kawai)

CD16_ego.1.down <- enrichGO(gene       = CD16.rest.down.df$ENTREZID,
                       universe      = universe.gene.list$ENTREZID,
                       OrgDb         = org.Hs.eg.db,
                       ont           = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05)

summary(CD16_ego.1.down)
CD16_ego.1.down <- setReadable(CD16_ego.1.down, OrgDb = org.Hs.eg.db)
#write.csv(file = "./Bulk_seq_momo_Count/DESEQ/GO_CD16rest.down.csv",CD16_ego.1.down@result)

dotplot(ego.1.down, showCategory=30)

goplot(ego.1.down)

###########################################################################################
#this is to check how many genes/pathways are over lap between conditions
CD14_rest_LPS_up_overlap <- intersect(filter(cd14.rest.up, log2FoldChange > 1.5)$gene, filter(lps.cd14.rest.up,log2FoldChange > 1.5)$gene)
CD14_rest_LPS_up_overlap

CD14_rest_LPS_down_overlap <- intersect(filter(cd14.rest.down, log2FoldChange <  -1.5)$gene, filter(lps.cd14.rest.down,log2FoldChange < -1.5)$gene)
CD14_rest_LPS_down_overlap

dim(filter(lps.cd14.rest.up,log2FoldChange > 1.5))
dim(filter(cd14.rest.down, log2FoldChange <  -1.5))

intersect(ego.1.up@result$Description[1:20], lps.ego.1.up@result$Description[1:20])
intersect(ego.1.down@result$Description[1:20], lps.ego.1.down@result$Description[1:20])

dim(lps.cd14.rest.up)
dim(lps.cd14.rest.down)

lps.ego.1.up
lps.ego.1.down

dim(CD16.rest.up)
dim(CD16.rest.down)

CD16_rest_LPS_up_overlap <- intersect(filter(CD16.rest.up, log2FoldChange > 1.5)$gene, filter(lps.cd14.rest.up,log2FoldChange > 1.5)$gene)
CD16_rest_LPS_up_overlap

CD16_rest_LPS_down_overlap <- intersect(filter(CD16.rest.down, log2FoldChange <  -1.5)$gene, filter(lps.cd14.rest.down,log2FoldChange < -1.5)$gene)
CD16_rest_LPS_down_overlap

intersect(CD16_ego.1.up@result$Description[1:20], lps.ego.1.up@result$Description[1:20])
intersect(CD16_ego.1.down@result$Description[1:20], lps.ego.1.down@result$Description[1:20])

CD14_rest_LPS_up_overlap
##############################################################################################################
#Disease Onotology 
CD14_rest_LPS_up_overlap.df <- bitr(CD14_rest_LPS_up_overlap, fromType = "SYMBOL",
                        toType = "ENTREZID",
                        OrgDb = org.Hs.eg.db)

Overlap_LPS_CD14rest_GO_up <- enrichDO(gene = CD14_rest_LPS_up_overlap.df$ENTREZID,
              ont           = "DO",
              pvalueCutoff  = 0.05,
              pAdjustMethod = "BH",
              universe      = universe.gene.list$ENTREZID,
              minGSSize     = 5,
              maxGSSize     = 500,
              qvalueCutoff  = 0.05,
              readable      = FALSE)

Overlap_LPS_CD14rest_GO_up <- setReadable(Overlap_LPS_CD14rest_GO_up, OrgDb = org.Hs.eg.db)
summary(Overlap_LPS_CD14rest_GO_up)
p.over.DO <- dotplot(Overlap_LPS_CD14rest_GO_up, showCategory=15,orderBy = "GeneRatio")
p.over.DO
#ggsave(filename = "./Bulk_seq_momo_Count/Overlapping.DO.pdf", p.over.DO, height = 8, width = 10, dpi = 300)

PatientCD14_disease_up <-  enrichDO(gene = cd14.rest.up.df$ENTREZID,
                                       ont           = "DO",
                                       pvalueCutoff  = 0.05,
                                       pAdjustMethod = "BH",
                                       universe      = universe.gene.list$ENTREZID,
                                       minGSSize     = 5,
                                       maxGSSize     = 500,
                                       qvalueCutoff  = 0.05,
                                       readable      = FALSE)

PatientCD14_disease_up <- setReadable(PatientCD14_disease_up,OrgDb = org.Hs.eg.db)
summary(PatientCD14_disease_up)
p.CD14rest.DO <- dotplot(PatientCD14_disease_up, showCategory=15,orderBy = "GeneRatio")
p.CD14rest.DO
#ggsave(filename = "./Bulk_seq_momo_Count/CD14Patient.DO.pdf", p.CD14rest.DO, height = 8, width = 10, dpi = 300)








