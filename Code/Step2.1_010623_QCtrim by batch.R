#Single cell Lyme disease Batch setup
#Author: Xin Chen, Ph.D.
#Date Created: 010523
# Assign batch and lane

#set up working directory
setwd("/Users/xinchen/Library/CloudStorage/Box-Box/Xin.Chen.Shareable/R7_NR")
getwd()
load("./Data/Step1.2_010523_assign_batch_lane.RData")

#load necessary package
library(data.table)
library(Seurat)
library(ggplot2)
library(stringr)
library(patchwork)

# QC for batch 1
#set a good estimate
pbmc.batch1 <- subset(pbmc,batch =="batch1")
plot1 <- FeatureScatter(pbmc.batch1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 <- plot1 + facet_wrap(.~colors, scales = "free")
plot1 <- plot1 + geom_vline(xintercept = 100) + geom_hline(yintercept = 20) +geom_vline(xintercept = 2500)
plot1 <- plot1 + geom_abline(intercept = 10, slope = 1/30)
plot1
png(filename ="./Files/1-Batch1_before.png")
plot1
dev.off()

# trim data
pbmc.batch1 <- subset(pbmc, subset = nFeature_RNA > 20 & nCount_RNA < 2500)
pbmc.batch1 <- subset(pbmc.batch1, subset = nFeature_RNA > 1/30 * nCount_RNA)
plot2 <- FeatureScatter(pbmc.batch1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2 <- plot2 + facet_wrap(.~colors, scales = "free")
plot2
png(filename ="./Files/2-Batch1_after.png")
plot2
dev.off()


# QC for batch 2
#set a good estimate
pbmc.batch2 <- subset(pbmc,batch =="batch2")
plot1 <- FeatureScatter(pbmc.batch2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 <- plot1 + facet_wrap(.~colors, scales = "free")
plot1 <- plot1 + geom_vline(xintercept = 100) + geom_hline(yintercept = 20) +geom_vline(xintercept = 3000)
plot1 <- plot1 + geom_abline(intercept = 30, slope = 1/35)
plot1

png(filename ="./Files/3-Batch2_before.png")
plot1
dev.off()

# trim data
pbmc.batch2 <- subset(pbmc.batch2, subset = nFeature_RNA > 30 & nCount_RNA < 3000)
pbmc.batch2 <- subset(pbmc.batch2, subset = nFeature_RNA > 1/35 * nCount_RNA)
plot2 <- FeatureScatter(pbmc.batch2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2 <- plot2 + facet_wrap(.~colors, scales = "free")
plot2
png(filename ="./Files/4-Batch2_after.png")
plot2
dev.off()
b2 <- plot2

# QC for batch 3
#set a good estimate
pbmc.batch3 <- subset(pbmc,batch =="batch3")
plot1 <- FeatureScatter(pbmc.batch3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 <- plot1 + facet_wrap(.~colors, scales = "free")
plot1 <- plot1 + geom_vline(xintercept = 50) + geom_hline(yintercept = 20) +geom_vline(xintercept = 1500)
plot1 <- plot1 + geom_abline(intercept = 10, slope = 1/20)
plot1

png(filename ="./Files/5-Batch3_before.png")
plot1
dev.off()

# trim data
pbmc.batch3 <- subset(pbmc.batch3, subset = nFeature_RNA > 10 & nCount_RNA < 1500)
pbmc.batch3 <- subset(pbmc.batch3, subset = nFeature_RNA > 0.05 * nCount_RNA)
plot2 <- FeatureScatter(pbmc.batch3, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2 <- plot2 + facet_wrap(.~colors, scales = "free")
plot2
png(filename ="./Files/6-Batch3_after.png")
plot2
dev.off()
b3 <-plot2

# QC for batch 4
#set a good estimate
pbmc.batch4 <- subset(pbmc,batch =="batch4")
plot1 <- FeatureScatter(pbmc.batch4, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 <- plot1 + facet_wrap(.~colors, scales = "free")
plot1 <- plot1 + geom_vline(xintercept = 50) + geom_hline(yintercept = 20) +geom_vline(xintercept = 1100)
plot1 <- plot1 + geom_abline(intercept = 10, slope = 1/15)
plot1

png(filename ="./Files/7-Batch4_before.png")
plot1
dev.off()

# trim data
pbmc.batch4 <- subset(pbmc.batch4, subset = nFeature_RNA > 20 & nCount_RNA < 1100)
pbmc.batch4 <- subset(pbmc.batch4, subset = nFeature_RNA > 1/15 * nCount_RNA)
plot2 <- FeatureScatter(pbmc.batch4, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2 <- plot2 + facet_wrap(.~colors, scales = "free")
plot2
png(filename ="./Files/8-Batch4_after.png")
plot2
dev.off()

b4 <-plot2
combine_trim_batch <- b1+b2+b3+b4

png(filename ="./Files/8-combine_Batch_after.png", width = 2000, height = 800)
combine_trim_batch
dev.off()


# QC for batch 5
#set a good estimate
pbmc.batch5 <- subset(pbmc,batch =="batch5")
plot1 <- FeatureScatter(pbmc.batch5, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 <- plot1 + facet_wrap(.~colors, scales = "free")
plot1 <- plot1 + geom_vline(xintercept = 100) + geom_hline(yintercept = 20) +geom_vline(xintercept = 3000)
plot1 <- plot1 + geom_abline(intercept = 30, slope = 1/35)
plot1

png(filename ="./Files/10-Batch5_before.png")
plot1
dev.off()

# trim data
pbmc.batch5 <- subset(pbmc.batch5, subset = nFeature_RNA > 30 & nCount_RNA < 3000)
pbmc.batch5 <- subset(pbmc.batch5, subset = nFeature_RNA > 1/35 * nCount_RNA)
plot2 <- FeatureScatter(pbmc.batch5, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2 <- plot2 + facet_wrap(.~colors, scales = "free")
plot2
png(filename ="./Files/11-Batch5_after.png")
plot2
dev.off()
b2 <- plot2

#merge pbmc with trimed batches
pbmc <- merge(pbmc.batch1, y = c(pbmc.batch2, pbmc.batch3, pbmc.batch4,pbmc.batch5))
pbmc[[]]
#Check batch after trim
table(pbmc$subject[pbmc$batch == "batch1"])
table(pbmc$subject[pbmc$batch == "batch2"])
table(pbmc$subject[pbmc$batch == "batch3"])
table(pbmc$subject[pbmc$batch == "batch4"])
table(pbmc$subject[pbmc$batch == "batch5"])

#check QC for combined data
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "batch", raster = FALSE)
plot1

png(filename ="./Files/9-after_trim_combine.png", width = 2000, height = 800)
plot1
dev.off()

save(pbmc, file = "./Data/Step2.1_010623_QC_trim_by_batch.RData")






