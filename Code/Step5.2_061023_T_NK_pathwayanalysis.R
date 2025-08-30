#Pathway Analysis
#05-08-2024
#https://yulab-smu.top/biomedical-knowledge-mining-book/index.html
#https://yulab-smu.top/biomedical-knowledge-mining-book/wikipathways-analysis.html

library(Seurat)
library(MAST)
library(clusterProfiler)
library(enrichplot)
library(pathview)
library(dplyr)
library(stringr)
library(Scillus)
library(ggplot2)
library(ComplexHeatmap)
library(ggrepel)
library(ggpubr)
library(DESeq2)
library(EnhancedVolcano)
library(cowplot)
library(org.Hs.eg.db)
library(reshape2)

#windows
setwd("C:/Users/zhoux/Box/Xin.Chen.Shareable/R7_NR/")
getwd()

source("../083123-FOXP3 GZMB KO analysis/Integreted/Code/xztools.R")
source("./Code/00_colorKey.R")
patient_data <- read.csv("./Results/patient_data_addPTLDN.csv", header = T)
load("./Data/Step5.1_032723_T_NK_clean.RData")
PTLDN_HD_ALL <- read.csv("../lyme_disease/Manuscript/Figures/Figure2/All_lympho_PTLDN_HD_DEG.csv")

######################################################################################################
universe.gene <- bitr(rownames(T_NK.clean), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

#Treg
Treg_up <- filter(PTLDN_HD_ALL, celltype == "Treg" & avg_log2FC > 0)
Treg_down <- filter(PTLDN_HD_ALL, celltype == "Treg" & avg_log2FC < 0)

Treg_up$gene

Treg_up.df <- bitr(Treg_up$gene, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)

Treg_down.df <- bitr(Treg_down$gene, fromType = "SYMBOL",
                   toType = c("ENSEMBL", "ENTREZID"),
                   OrgDb = org.Hs.eg.db)

treg.up.ego <- enrichGO(gene  = Treg_up.df$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                universe      = universe.gene,
                ont           = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.01,
                qvalueCutoff  = 0.05,
                readable      = TRUE)

treg.down.ego <- enrichGO(gene        = Treg_down.df$ENTREZID,
                        OrgDb         = org.Hs.eg.db,
                        universe      = universe.gene,
                        ont           = "ALL",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.01,
                        qvalueCutoff  = 0.05,
                        readable      = TRUE)

head(treg.up.ego)
head(treg.down.ego)

treg.up.goresult.up <- treg.up.ego@result
treg.up.goresult.up$Class <- "UpRegulated"
treg.up.goresult.down <- treg.down.ego@result
treg.up.goresult.down$Class <- "DownRegulated"

treg.up.goresult <- rbind(treg.up.goresult.up,treg.up.goresult.down)
treg.up.goresult$CellType <- "Treg"

treg.up.dotplot <- dotplot(treg.up.ego, title = "Treg.up")
treg.down.dotplot <- dotplot(treg.down.ego, title = "Treg.down")

treg.dotplot <- treg.up.dotplot + treg.down.dotplot
treg.dotplot

# 
# treg.up.kegg <- enrichKEGG(gene  = Treg_up.df$ENTREZID,
#                            organism     = 'hsa',
#                            qvalueCutoff  = 0.05)
# 
# treg.down.kegg <- enrichKEGG(gene  = Treg_down.df$ENTREZID,
#                            organism     = 'hsa',
#                            qvalueCutoff  = 0.05)
# 
# head(treg.up.kegg)
# head(treg.down.kegg)
# 
# treg.up.wiki <- enrichWP(Treg_up.df$ENTREZID, organism = "Homo sapiens") 
# treg.down.wiki <- enrichWP(Treg_down.df$ENTREZID, organism = "Homo sapiens") 
# 
# head(treg.up.wiki)
# head(treg.down.wiki)

# GammaDelta_T
GammaDelta_T_up <- filter(PTLDN_HD_ALL, celltype == "GammaDelta_T" & avg_log2FC > 0)
GammaDelta_T_down <- filter(PTLDN_HD_ALL, celltype == "GammaDelta_T" & avg_log2FC < 0)

GammaDelta_T_up.df <- bitr(GammaDelta_T_up$gene, fromType = "SYMBOL",
                           toType = c("ENSEMBL", "ENTREZID"),
                           OrgDb = org.Hs.eg.db)

GammaDelta_T_down.df <- bitr(GammaDelta_T_down$gene, fromType = "SYMBOL",
                             toType = c("ENSEMBL", "ENTREZID"),
                             OrgDb = org.Hs.eg.db)

gammaDelta_T_up.ego <- enrichGO(gene  = GammaDelta_T_up.df$ENTREZID,
                                OrgDb         = org.Hs.eg.db,
                                universe      = universe.gene,
                                ont           = "ALL",
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.01,
                                qvalueCutoff  = 0.05,
                                readable      = TRUE)

gammaDelta_T_down.ego <- enrichGO(gene        = GammaDelta_T_down.df$ENTREZID,
                                  OrgDb         = org.Hs.eg.db,
                                  universe      = universe.gene,
                                  ont           = "ALL",
                                  pAdjustMethod = "BH",
                                  pvalueCutoff  = 0.01,
                                  qvalueCutoff  = 0.05,
                                  readable      = TRUE)

gammaDelta_T_up.goresult.up <- gammaDelta_T_up.ego@result
gammaDelta_T_up.goresult.up$Class <- "UpRegulated"
gammaDelta_T_up.goresult.down <- gammaDelta_T_down.ego@result
gammaDelta_T_up.goresult.down$Class <- "DownRegulated"

gammaDelta_T_up.goresult <- rbind(gammaDelta_T_up.goresult.up, gammaDelta_T_up.goresult.down)
gammaDelta_T_up.goresult$CellType <- "GammaDelta_T"

gammaDelta_T_up.dotplot <- dotplot(gammaDelta_T_up.ego, title = "GammaDelta_T.up")
gammaDelta_T_down.dotplot <- dotplot(gammaDelta_T_down.ego, title = "GammaDelta_T.down")

gammaDelta_T_dotplot <- gammaDelta_T_up.dotplot + gammaDelta_T_down.dotplot
gammaDelta_T_dotplot

# MAIT
MAIT_up <- filter(PTLDN_HD_ALL, celltype == "MAIT" & avg_log2FC > 0)
MAIT_down <- filter(PTLDN_HD_ALL, celltype == "MAIT" & avg_log2FC < 0)

MAIT_up.df <- bitr(MAIT_up$gene, fromType = "SYMBOL",
                   toType = c("ENSEMBL", "ENTREZID"),
                   OrgDb = org.Hs.eg.db)

MAIT_down.df <- bitr(MAIT_down$gene, fromType = "SYMBOL",
                     toType = c("ENSEMBL", "ENTREZID"),
                     OrgDb = org.Hs.eg.db)

mait_up.ego <- enrichGO(gene  = MAIT_up.df$ENTREZID,
                        OrgDb         = org.Hs.eg.db,
                        universe      = universe.gene,
                        ont           = "ALL",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.01,
                        qvalueCutoff  = 0.05,
                        readable      = TRUE)

mait_down.ego <- enrichGO(gene        = MAIT_down.df$ENTREZID,
                          OrgDb         = org.Hs.eg.db,
                          universe      = universe.gene,
                          ont           = "ALL",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.01,
                          qvalueCutoff  = 0.05,
                          readable      = TRUE)

mait_up.goresult.up <- mait_up.ego@result
mait_up.goresult.up$Class <- "UpRegulated"
mait_down.goresult.down <- mait_down.ego@result
mait_down.goresult.down$Class <- "DownRegulated"

mait_up.goresult <- rbind(mait_up.goresult.up, mait_down.goresult.down)
mait_up.goresult$CellType <- "MAIT"

mait_up.dotplot <- dotplot(mait_up.ego, title = "MAIT.up")
mait_down.dotplot <- dotplot(mait_down.ego, title = "MAIT.down")

mait_dotplot <- mait_up.dotplot + mait_down.dotplot
mait_dotplot

# Memory_CD4
Memory_CD4_up <- filter(PTLDN_HD_ALL, celltype == "Memory_CD4" & avg_log2FC > 0)
Memory_CD4_down <- filter(PTLDN_HD_ALL, celltype == "Memory_CD4" & avg_log2FC < 0)

Memory_CD4_up.df <- bitr(Memory_CD4_up$gene, fromType = "SYMBOL",
                         toType = c("ENSEMBL", "ENTREZID"),
                         OrgDb = org.Hs.eg.db)

Memory_CD4_down.df <- bitr(Memory_CD4_down$gene, fromType = "SYMBOL",
                           toType = c("ENSEMBL", "ENTREZID"),
                           OrgDb = org.Hs.eg.db)

memory_CD4_up.ego <- enrichGO(gene  = Memory_CD4_up.df$ENTREZID,
                              OrgDb         = org.Hs.eg.db,
                              universe      = universe.gene,
                              ont           = "ALL",
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.01,
                              qvalueCutoff  = 0.05,
                              readable      = TRUE)

memory_CD4_down.ego <- enrichGO(gene        = Memory_CD4_down.df$ENTREZID,
                                OrgDb         = org.Hs.eg.db,
                                universe      = universe.gene,
                                ont           = "ALL",
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.01,
                                qvalueCutoff  = 0.05,
                                readable      = TRUE)

memory_CD4_up.goresult.up <- memory_CD4_up.ego@result
memory_CD4_up.goresult.up$Class <- "UpRegulated"
memory_CD4_down.goresult.down <- memory_CD4_down.ego@result
memory_CD4_down.goresult.down$Class <- "DownRegulated"

memory_CD4_up.goresult <- rbind(memory_CD4_up.goresult.up, memory_CD4_down.goresult.down)
memory_CD4_up.goresult$CellType <- "Memory_CD4"

memory_CD4_up.dotplot <- dotplot(memory_CD4_up.ego, title = "Memory_CD4.up")
memory_CD4_down.dotplot <- dotplot(memory_CD4_down.ego, title = "Memory_CD4.down")

memory_CD4_dotplot <- memory_CD4_up.dotplot + memory_CD4_down.dotplot
memory_CD4_dotplot

# Memory_CD8
Memory_CD8_up <- filter(PTLDN_HD_ALL, celltype == "Memory_CD8" & avg_log2FC > 0)
Memory_CD8_down <- filter(PTLDN_HD_ALL, celltype == "Memory_CD8" & avg_log2FC < 0)

Memory_CD8_up.df <- bitr(Memory_CD8_up$gene, fromType = "SYMBOL",
                         toType = c("ENSEMBL", "ENTREZID"),
                         OrgDb = org.Hs.eg.db)

Memory_CD8_down.df <- bitr(Memory_CD8_down$gene, fromType = "SYMBOL",
                           toType = c("ENSEMBL", "ENTREZID"),
                           OrgDb = org.Hs.eg.db)

memory_CD8_up.ego <- enrichGO(gene  = Memory_CD8_up.df$ENTREZID,
                              OrgDb         = org.Hs.eg.db,
                              universe      = universe.gene,
                              ont           = "ALL",
                              pAdjustMethod = "BH",
                              pvalueCutoff  = 0.01,
                              qvalueCutoff  = 0.05,
                              readable      = TRUE)

memory_CD8_down.ego <- enrichGO(gene        = Memory_CD8_down.df$ENTREZID,
                                OrgDb         = org.Hs.eg.db,
                                universe      = universe.gene,
                                ont           = "ALL",
                                pAdjustMethod = "BH",
                                pvalueCutoff  = 0.01,
                                qvalueCutoff  = 0.05,
                                readable      = TRUE)

memory_CD8_up.goresult.up <- memory_CD8_up.ego@result
memory_CD8_up.goresult.up$Class <- "UpRegulated"
memory_CD8_down.goresult.down <- memory_CD8_down.ego@result
memory_CD8_down.goresult.down$Class <- "DownRegulated"

memory_CD8_up.goresult <- rbind(memory_CD8_up.goresult.up, memory_CD8_down.goresult.down)
memory_CD8_up.goresult$CellType <- "Memory_CD8"

memory_CD8_up.dotplot <- dotplot(memory_CD8_up.ego, title = "Memory_CD8.up")
memory_CD8_down.dotplot <- dotplot(memory_CD8_down.ego, title = "Memory_CD8.down")

memory_CD8_dotplot <- memory_CD8_up.dotplot + memory_CD8_down.dotplot
memory_CD8_dotplot

# NK
NK_up <- filter(PTLDN_HD_ALL, celltype == "NK" & avg_log2FC > 0)
NK_down <- filter(PTLDN_HD_ALL, celltype == "NK" & avg_log2FC < 0)

NK_up.df <- bitr(NK_up$gene, fromType = "SYMBOL",
                 toType = c("ENSEMBL", "ENTREZID"),
                 OrgDb = org.Hs.eg.db)

NK_down.df <- bitr(NK_down$gene, fromType = "SYMBOL",
                   toType = c("ENSEMBL", "ENTREZID"),
                   OrgDb = org.Hs.eg.db)

nk_up.ego <- enrichGO(gene  = NK_up.df$ENTREZID,
                      OrgDb         = org.Hs.eg.db,
                      universe      = universe.gene,
                      ont           = "ALL",
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)

nk_down.ego <- enrichGO(gene        = NK_down.df$ENTREZID,
                        OrgDb         = org.Hs.eg.db,
                        universe      = universe.gene,
                        ont           = "ALL",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.01,
                        qvalueCutoff  = 0.05,
                        readable      = TRUE)

nk_up.goresult.up <- nk_up.ego@result
nk_up.goresult.up$Class <- "UpRegulated"
nk_down.goresult.down <- nk_down.ego@result
nk_down.goresult.down$Class <- "DownRegulated"

nk_up.goresult <- rbind(nk_up.goresult.up, nk_down.goresult.down)
nk_up.goresult$CellType <- "NK"

nk_up.dotplot <- dotplot(nk_up.ego, title = "NK.up")
nk_down.dotplot <- dotplot(nk_down.ego, title = "NK.down")

nk_dotplot <- nk_up.dotplot + nk_down.dotplot
nk_dotplot

# TFH
TFH_up <- filter(PTLDN_HD_ALL, celltype == "TFH" & avg_log2FC > 0)
TFH_down <- filter(PTLDN_HD_ALL, celltype == "TFH" & avg_log2FC < 0)

TFH_up.df <- bitr(TFH_up$gene, fromType = "SYMBOL",
                  toType = c("ENSEMBL", "ENTREZID"),
                  OrgDb = org.Hs.eg.db)

TFH_down.df <- bitr(TFH_down$gene, fromType = "SYMBOL",
                    toType = c("ENSEMBL", "ENTREZID"),
                    OrgDb = org.Hs.eg.db)

tfh_up.ego <- enrichGO(gene  = TFH_up.df$ENTREZID,
                       OrgDb         = org.Hs.eg.db,
                       universe      = universe.gene,
                       ont           = "ALL",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05,
                       readable      = TRUE)

tfh_down.ego <- enrichGO(gene        = TFH_down.df$ENTREZID,
                         OrgDb         = org.Hs.eg.db,
                         universe      = universe.gene,
                         ont           = "ALL",
                         pAdjustMethod = "BH",
                         pvalueCutoff  = 0.01,
                         qvalueCutoff  = 0.05,
                         readable      = TRUE)

tfh_up.goresult.up <- tfh_up.ego@result
tfh_up.goresult.up$Class <- "UpRegulated"
tfh_down.goresult.down <- tfh_down.ego@result
tfh_down.goresult.down$Class <- "DownRegulated"

tfh_up.goresult <- rbind(tfh_up.goresult.up, tfh_down.goresult.down)
tfh_up.goresult$CellType <- "TFH"

tfh_up.dotplot <- dotplot(tfh_up.ego, title = "TFH.up")
tfh_down.dotplot <- dotplot(tfh_down.ego, title = "TFH.down")

tfh_dotplot <- tfh_up.dotplot + tfh_down.dotplot
tfh_dotplot

pathway_result <- rbind(treg.up.goresult, gammaDelta_T_up.goresult, mait_up.goresult,memory_CD4_up.goresult,
                        memory_CD8_up.goresult, nk_up.goresult, tfh_up.goresult)
#write.csv(file = "./Results/pathway.analysis.result.csv", pathway_result)

pathwayresult.analysis <- table(pathway_result$Description, pathway_result$Class) %>% data.frame() %>% 
  dcast(Var1 ~ Var2, value.var = "Freq")



write.csv(file = "./pathway.result.csv",pathwayresult.analysis)


