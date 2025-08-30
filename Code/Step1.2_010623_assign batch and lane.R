#Single cell Lyme disease Batch setup
#Author: Xin Chen, Ph.D.
#Date Created: 010523
# Assign batch and lane

#set up working directory
setwd("/Users/xinchen/Library/CloudStorage/Box-Box/Xin.Chen.Shareable/R7_NR")
getwd()
load("./Data/Step1.1_010523_load_data.RData")

#load necessary package
library(data.table)
library(Seurat)
library(ggplot2)
library(stringr)
library(patchwork)

table(pbmc$subject)

# set up lane
pbmc[["lane"]] <-rownames(pbmc[[]])
subject <-pbmc$subject
pbmc$lane[which(str_detect(subject, "HD1|RTH1V7b1|PTLD1V5"))] <- "124a"
pbmc$lane[which(str_detect(subject, "HD2|PTLD2V1|RTH2V5b1"))] <- "124b"
pbmc$lane[which(str_detect(subject, "HD8|PTLD1V3|RTH6V3"))] <- "124c"
pbmc$lane[which(str_detect(subject, "HD7|PTLD4V3|RTH5V3b2"))] <- "124d"
pbmc$lane[which(str_detect(subject, "HD5|PTLD2V5|RTH4V3|RTH4V7"))] <- "124f"
pbmc$lane[which(str_detect(subject, "HD6|PTLD3V3|RTH1V3|RTH1V5"))] <- "124h"
pbmc$lane[which(str_detect(subject, "PTLDN1|PTLDN2|PTLDN3"))] <- "124i"
pbmc$lane[which(str_detect(subject, "PTLDN4|RTH3V3"))] <- "124j"
pbmc$lane[which(str_detect(subject, "HD3|PTLDN5|PTLDN6"))] <- "124k"
pbmc$lane[which(str_detect(subject, "HD4|RTH3V5|PTLDN7|PTLDN8"))] <- "124L"
pbmc$lane[which(str_detect(subject, "HD9|PTLD1V2|RTH2V2|RTH2V5b5|RTH2V7"))] <- "124m"
pbmc$lane[which(str_detect(subject, "HD10|PTLD2V2|PTLD2V3|RTH6V3|RTH6V5"))] <- "124n"
pbmc$lane[which(str_detect(subject, "HD11|PTLD3V2|PTLD3V5|RTH5V2|RTH5V3b5"))] <- "124O"
pbmc$lane[which(str_detect(subject, "HD12|PTLD4V7|RTH4V2|RTH4V5|RTH1V2"))] <- "124P"
pbmc$lane[which(str_detect(subject, "PTLD2V7|PTLD3V7|RTH6V7|RTH5V7|RTH1V7b5"))] <- "124Q"
table(pbmc$lane)

# set up batch
pbmc[["batch"]] <-rownames(pbmc[[]])
pbmc$batch[which(str_detect(subject, "HD1|RTH1V7b1|PTLD1V5"))] <- "batch1"
pbmc$batch[which(str_detect(subject, "HD2|PTLD2V1|RTH2V5b1"))] <- "batch1"
pbmc$batch[which(str_detect(subject, "HD8|PTLD1V3|RTH6V3"))] <- "batch2"
pbmc$batch[which(str_detect(subject, "HD7|PTLD4V3|RTH5V3b2"))] <- "batch2"
pbmc$batch[which(str_detect(subject, "HD5|PTLD2V5|RTH4V3|RTH4V7"))] <- "batch3"
pbmc$batch[which(str_detect(subject, "HD6|PTLD3V3|RTH1V3|RTH1V5"))] <- "batch3"
pbmc$batch[which(str_detect(subject, "PTLDN1|PTLDN2|PTLDN3"))] <- "batch4"
pbmc$batch[which(str_detect(subject, "PTLDN4|RTH3V3"))] <- "batch4"
pbmc$batch[which(str_detect(subject, "HD3|PTLDN5|PTLDN6"))] <- "batch4"
pbmc$batch[which(str_detect(subject, "HD4|RTH3V5|PTLDN7|PTLDN8"))] <- "batch4"
pbmc$batch[which(str_detect(subject, "HD9|PTLD1V2|RTH2V2|RTH2V5b5|RTH2V7"))] <- "batch5"
pbmc$batch[which(str_detect(subject, "HD10|PTLD2V2|PTLD2V3|RTH6V3|RTH6V5"))] <- "batch5"
pbmc$batch[which(str_detect(subject, "HD11|PTLD3V2|PTLD3V5|RTH5V2|RTH5V3b5"))] <- "batch5"
pbmc$batch[which(str_detect(subject, "HD12|PTLD4V7|RTH4V2|RTH4V5|RTH1V2"))] <- "batch5"
pbmc$batch[which(str_detect(subject, "PTLD2V7|PTLD3V7|RTH6V7|RTH5V7|RTH1V7b5"))] <- "batch5"
table(pbmc$subject)
table(pbmc$batch)
pbmc[[]]

# to check if object assigned to right batch
table(pbmc$subject[pbmc$batch == "batch1"])
table(pbmc$subject[pbmc$batch == "batch2"])
table(pbmc$subject[pbmc$batch == "batch3"])
table(pbmc$subject[pbmc$batch == "batch4"])
table(pbmc$subject[pbmc$batch == "batch5"])

#save
save(pbmc, file = "./Data/Step1.2_010523_assign_batch_lane.RData")
