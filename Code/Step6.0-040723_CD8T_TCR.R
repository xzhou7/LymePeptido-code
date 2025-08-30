#Single cell Lyme disease post merge
#Author: Xin Chen, Ph.D.
#Date Created: 032723
#T cell analysis

#set up working directory
setwd(dir = "~/Library/CloudStorage/Box-Box/Xin.Chen.Shareable/R7_NR/CD8TCR")
getwd()

#load necessary package
library(immunarch)
library(dplyr)
library(reshape2)
library(stringr)
library(tibble)
library(factoextra)
library(immunarch)
immdata <- reLoad('example')

load("~/Library/CloudStorage/Box-Box/Xin.Chen.Shareable/R7_NR/Data/Step5.1_032923_CD8T_clean.RData")

CD8T.clean$cluster <- "Naive"
CD8T.clean$cluster[CD8T.clean$integrated_snn_res.1 == 2] <- "GZMB+KIR+"
CD8T.clean$cluster[CD8T.clean$integrated_snn_res.1 == 6] <- "GZMK+"
CD8T.clean$cluster[CD8T.clean$integrated_snn_res.1 == 7] <- "CCR6+CD161+"
CD8T.clean$cluster[CD8T.clean$integrated_snn_res.1 == 12] <- "TCRgd+"
CD8T.clean$cluster[CD8T.clean$integrated_snn_res.1 == 15] <- "GNLY+"
meta.seurat <- CD8T.clean[[]]
meta.seurat
meta.seurat.TCR <- select(meta.seurat,condition, subject, cluster)
meta.seurat.TCR$cell.id <- rownames(meta.seurat.TCR)
meta.seurat.TCR$cell.id <- gsub("....$", "", meta.seurat.TCR$cell.id)


table(meta.seurat.TCR$cell.id)
write.csv(file = "./meta.seurat.tcr.csv", meta.seurat.TCR)

meta.seurat.TCR <- read.csv("./meta.seurat.tcr.csv", header = T)
depth <- table(meta.seurat.TCR$subject) %>% as.data.frame()

depth


#load 124a
VDJ.124a <- read.csv("/Users/xinchen/Desktop/Lyme/data/VDJ/124a_1_Seq1merged124aTG_VDJ_perCell.csv", header = T)
ST.124a <- read.csv("/Users/xinchen/Desktop/Lyme/data/VDJ/124a_1_Seq1merged124aTG_Sample_Tag_Calls.csv", header = T)

ST.124a <- filter(ST.124a, Sample_Tag != "Multiplet") %>% filter(Sample_Tag != "Undetermined")
ST.124a$Subject <- "unknown"
ST.124a$Subject[(str_detect(ST.124a$Sample_Tag, "SampleTag01"))] <- "HD1"
ST.124a$Subject[(str_detect(ST.124a$Sample_Tag, "SampleTag02"))] <- "RTH1V7"
ST.124a$Subject[(str_detect(ST.124a$Sample_Tag, "SampleTag03"))] <- "PTLD1V5"
ST.124a$cell.id <- paste(ST.124a$Subject, ST.124a$Cell_Index, sep="_")

ST.124a.meta <- merge(ST.124a,meta.seurat.TCR, by='cell.id')
VDJ.124a <- filter(VDJ.124a, BCR_Paired_Chains == "TRUE" | TCR_Paired_Chains == "TRUE")
VDJ.124a.meta <- merge(VDJ.124a,ST.124a.meta, by='Cell_Index')
table(VDJ.124a.meta$subject)
#write.csv(file = "./VDJ.124a.meta.csv",VDJ.124a.meta)

#Creat a table for GLIPH


############################################################################################################################################
#load 124b
ST.124b <- read.csv("/Users/xinchen/Desktop/Lyme/data/VDJ/124b_1_Seq1newmerge124b_Sample_Tag_Calls.csv", header = T)
VDJ.124b <- read.csv("/Users/xinchen/Desktop/Lyme/data/VDJ/124b_1_Seq1newmerge124b_VDJ_perCell.csv", header = T)

ST.124b <- filter(ST.124b, Sample_Tag != "Multiplet") %>% filter(Sample_Tag != "Undetermined")
ST.124b$Subject <- "unknown"
ST.124b$Subject[(str_detect(ST.124b$Sample_Tag, "SampleTag04"))] <- "HD2"
ST.124b$Subject[(str_detect(ST.124b$Sample_Tag, "SampleTag05"))] <- "RTH2V5"
ST.124b$Subject[(str_detect(ST.124b$Sample_Tag, "SampleTag06"))] <- "PTLD2V1"
ST.124b$cell.id <- paste(ST.124b$Subject, ST.124b$Cell_Index, sep="_")
ST.124b.meta <- merge(ST.124b,meta.seurat.TCR, by='cell.id')
VDJ.124b <- filter(VDJ.124b, BCR_Paired_Chains == "TRUE" | TCR_Paired_Chains == "TRUE")
VDJ.124b.meta <- merge(VDJ.124b,ST.124b.meta, by="Cell_Index")
table(VDJ.124b.meta$subject)
#################################################################################################################################################
#load 124c
ST.124c <- read.csv("/Users/xinchen/Desktop/Lyme/data/VDJ/124c_seq1-124c_Sample_Tag_Calls.csv", header = T)
VDJ.124c <- read.csv("/Users/xinchen/Desktop/Lyme/data/VDJ/124c_seq1-124c_VDJ_perCell.csv", header = T)

ST.124c <- filter(ST.124c, Sample_Tag != "Multiplet") %>% filter(Sample_Tag != "Undetermined")
ST.124c$Subject <- "unknown"
ST.124c$Subject[(str_detect(ST.124c$Sample_Tag, "SampleTag07"))] <- "HD8"
ST.124c$Subject[(str_detect(ST.124c$Sample_Tag, "SampleTag08"))] <- "RTH6V3"
ST.124c$Subject[(str_detect(ST.124c$Sample_Tag, "SampleTag09"))] <- "PTLD1V3"
ST.124c$cell.id <- paste(ST.124c$Subject, ST.124c$Cell_Index, sep="_")
ST.124c.meta <- merge(ST.124c, meta.seurat.TCR, by='cell.id')
VDJ.124c <- filter(VDJ.124c, BCR_Paired_Chains == "TRUE" | TCR_Paired_Chains == "TRUE")
VDJ.124c.meta <- merge(VDJ.124c,ST.124c.meta, by="Cell_Index")
table(meta.seurat.TCR$cell.id)
table(ST.124c$cell.id)
table(meta.seurat.TCR$subject)

"HD8_584768" %in% meta.seurat.TCR$cell.id

#################################################################################################################################################
#load 124d
ST.124d <- read.csv("/Users/xinchen/Desktop/Lyme/data/VDJ/124d_1_seq1-124d_Sample_Tag_Calls.csv", header = T)
VDJ.124d <- read.csv("/Users/xinchen/Desktop/Lyme/data/VDJ/124d_1_seq1-124d_VDJ_perCell.csv", header = T)

ST.124d <- filter(ST.124d, Sample_Tag != "Multiplet") %>% filter(Sample_Tag != "Undetermined")
ST.124d$Subject <- "unknown"
ST.124d$Subject[(str_detect(ST.124d$Sample_Tag, "SampleTag10"))] <- "HD7"
ST.124d$Subject[(str_detect(ST.124d$Sample_Tag, "SampleTag11"))] <- "RTH5V3"
ST.124d$Subject[(str_detect(ST.124d$Sample_Tag, "SampleTag12"))] <- "PTLD4V3"
ST.124d$cell.id <- paste(ST.124d$Subject, ST.124d$Cell_Index, sep="_")
ST.124d.meta <- merge(ST.124d,meta.seurat.TCR, by='cell.id')
VDJ.124d <- filter(VDJ.124d, BCR_Paired_Chains == "TRUE" | TCR_Paired_Chains == "TRUE")
VDJ.124d.meta <- merge(VDJ.124d,ST.124d.meta, by="Cell_Index")
table(ST.124d$Subject)

#################################################################################################################################################
#load 124h
ST.124h <- read.csv("/Users/xinchen/Desktop/Lyme/data/VDJ/124h_21540Jix_seq2-30_Sample_Tag_Calls.csv", header = T)
VDJ.124h <- read.csv("/Users/xinchen/Desktop/Lyme/data/VDJ/124h_21540Jix_seq2-30_VDJ_perCell.csv", header = T)

ST.124h <- filter(ST.124h, Sample_Tag != "Multiplet") %>% filter(Sample_Tag != "Undetermined")
ST.124h$Subject <- "unknown"
ST.124h$Subject[(str_detect(ST.124h$Sample_Tag, "SampleTag08"))] <- "HD6"
ST.124h$Subject[(str_detect(ST.124h$Sample_Tag, "SampleTag09"))] <- "RTH1V3"
ST.124h$Subject[(str_detect(ST.124h$Sample_Tag, "SampleTag10"))] <- "RTH1V5"
ST.124h$Subject[(str_detect(ST.124h$Sample_Tag, "SampleTag11"))] <- "PTLD3V3"
ST.124h$cell.id <- paste(ST.124h$Subject, ST.124h$Cell_Index, sep="_")
ST.124h.meta <- merge(ST.124h,meta.seurat.TCR, by='cell.id')
VDJ.124h <- filter(VDJ.124h, BCR_Paired_Chains == "TRUE" | TCR_Paired_Chains == "TRUE")
VDJ.124h.meta <- merge(VDJ.124h,ST.124h.meta, by="Cell_Index")


#################################################################################################################################################
#load 124i
ST.124i <- read.csv("/Users/xinchen/Desktop/Lyme/data/VDJ/124i_21540Jix_seq2-31_Sample_Tag_Calls.csv", header = T)
VDJ.124i <- read.csv("/Users/xinchen/Desktop/Lyme/data/VDJ/124i_21540Jix_seq2-31_VDJ_perCell.csv", header = T)

ST.124i <- filter(ST.124i, Sample_Tag != "Multiplet") %>% filter(Sample_Tag != "Undetermined")
ST.124i$Subject <- "unknown"
ST.124i$Subject[(str_detect(ST.124i$Sample_Tag, "SampleTag01"))] <- "PTLDN1"
ST.124i$Subject[(str_detect(ST.124i$Sample_Tag, "SampleTag02"))] <- "PTLDN2"
ST.124i$Subject[(str_detect(ST.124i$Sample_Tag, "SampleTag03"))] <- "PTLDN3"
ST.124i$cell.id <- paste(ST.124i$Subject, ST.124i$Cell_Index, sep="_")
ST.124i.meta <- merge(ST.124i,meta.seurat.TCR, by='cell.id')
VDJ.124i <- filter(VDJ.124i, BCR_Paired_Chains == "TRUE" | TCR_Paired_Chains == "TRUE")
VDJ.124i.meta <- merge(VDJ.124i,ST.124i.meta, by="Cell_Index")


#################################################################################################################################################
#load 124j
ST.124j <- read.csv("/Users/xinchen/Desktop/Lyme/data/VDJ/124j_21540Jix_seq2-32_Sample_Tag_Calls.csv", header = T)
VDJ.124j <- read.csv("/Users/xinchen/Desktop/Lyme/data/VDJ/124j_21540Jix_seq2-32_VDJ_perCell.csv", header = T)

ST.124j <- filter(ST.124j, Sample_Tag != "Multiplet") %>% filter(Sample_Tag != "Undetermined") %>% filter(Sample_Tag != "SampleTag07_hs") %>% filter(Sample_Tag != "SampleTag06_hs")
ST.124j$Subject <- "unknown"
ST.124j$Subject[(str_detect(ST.124j$Sample_Tag, "SampleTag04"))] <- "PTLDN4"
ST.124j$Subject[(str_detect(ST.124j$Sample_Tag, "SampleTag05"))] <- "RTH3V3"

ST.124j$cell.id <- paste(ST.124j$Subject, ST.124j$Cell_Index, sep="_")
ST.124j.meta <- merge(ST.124j,meta.seurat.TCR, by='cell.id')
VDJ.124j <- filter(VDJ.124j, BCR_Paired_Chains == "TRUE" | TCR_Paired_Chains == "TRUE")
VDJ.124j.meta <- merge(VDJ.124j,ST.124j.meta, by="Cell_Index")


#################################################################################################################################################
#load 124k
ST.124k <- read.csv("/Users/xinchen/Desktop/Lyme/data/VDJ/124k_21540Jix_seq2-33_Sample_Tag_Calls.csv", header = T)
VDJ.124k <- read.csv("/Users/xinchen/Desktop/Lyme/data/VDJ/124k_21540Jix_seq2-33_VDJ_perCell.csv", header = T)

ST.124k <- filter(ST.124k, Sample_Tag != "Multiplet") %>% filter(Sample_Tag != "Undetermined")
ST.124k$Subject <- "unknown"
ST.124k$Subject[(str_detect(ST.124k$Sample_Tag, "SampleTag01"))] <- "HD3"
ST.124k$Subject[(str_detect(ST.124k$Sample_Tag, "SampleTag03"))] <- "PTLDN5"
ST.124k$Subject[(str_detect(ST.124k$Sample_Tag, "SampleTag04"))] <- "PTLDN6"
ST.124k$cell.id <- paste(ST.124k$Subject, ST.124k$Cell_Index, sep="_")
ST.124k.meta <- merge(ST.124k,meta.seurat.TCR, by='cell.id')
VDJ.124k <- filter(VDJ.124k, BCR_Paired_Chains == "TRUE" | TCR_Paired_Chains == "TRUE")
VDJ.124k.meta <- merge(VDJ.124k,ST.124k.meta, by="Cell_Index")
table(VDJ.124k.meta$subject)

#################################################################################################################################################
#load 124L
ST.124L <- read.csv("/Users/xinchen/Desktop/Lyme/data/VDJ/124L_21540Jix_seq2-34_Sample_Tag_Calls.csv", header = T)
VDJ.124L <- read.csv("/Users/xinchen/Desktop/Lyme/data/VDJ/124L_21540Jix_seq2-34_VDJ_perCell.csv", header = T)

ST.124L <- filter(ST.124L, Sample_Tag != "Multiplet") %>% filter(Sample_Tag != "Undetermined")
ST.124L$Subject <- "unknown"
ST.124L$Subject[(str_detect(ST.124L$Sample_Tag, "SampleTag08"))] <- "HD4"
ST.124L$Subject[(str_detect(ST.124L$Sample_Tag, "SampleTag09"))] <- "RTH3V5"
ST.124L$Subject[(str_detect(ST.124L$Sample_Tag, "SampleTag10"))] <- "PTLDN7"
ST.124L$Subject[(str_detect(ST.124L$Sample_Tag, "SampleTag11"))] <- "PTLDN8"

ST.124L$cell.id <- paste(ST.124L$Subject, ST.124L$Cell_Index, sep="_")
ST.124L.meta <- merge(ST.124L,meta.seurat.TCR, by='cell.id')
VDJ.124L <- filter(VDJ.124L, BCR_Paired_Chains == "TRUE" | TCR_Paired_Chains == "TRUE")
VDJ.124L.meta <- merge(VDJ.124L,ST.124L.meta, by="Cell_Index")
table(VDJ.124L.meta$subject)

#################################################################################################################################################
#load 124F
ST.124F <- read.csv("/Users/xinchen/Desktop/Lyme/data/VDJ/124F_21540Jix_seq2-29_Sample_Tag_Calls.csv", header = T)
VDJ.124F <- read.csv("/Users/xinchen/Desktop/Lyme/data/VDJ/124F_21540Jix_seq2-29_VDJ_perCell copy.csv", header = T)

ST.124F <- filter(ST.124F, Sample_Tag != "Multiplet") %>% filter(Sample_Tag != "Undetermined")
ST.124F$Subject <- "unknown"
ST.124F$Subject[(str_detect(ST.124F$Sample_Tag, "SampleTag06"))] <- "HD5"
ST.124F$Subject[(str_detect(ST.124F$Sample_Tag, "SampleTag07"))] <- "RTH4V3"
ST.124F$Subject[(str_detect(ST.124F$Sample_Tag, "SampleTag08"))] <- "RTH4V7"
ST.124F$Subject[(str_detect(ST.124F$Sample_Tag, "SampleTag09"))] <- "PTLD2V5"
ST.124F$cell.id <- paste(ST.124F$Subject, ST.124F$Cell_Index, sep="_")
ST.124F.meta <- merge(ST.124F,meta.seurat.TCR, by='cell.id')
VDJ.124F <- filter(VDJ.124F, BCR_Paired_Chains == "TRUE" | TCR_Paired_Chains == "TRUE")
VDJ.124F.meta <- merge(VDJ.124F,ST.124F.meta, by="Cell_Index")
table(VDJ.124F.meta$subject)

####################################################################################################

#load 124m
ST.124M <- read.csv("/Users/xinchen/Desktop/Lyme/data/VDJ/124m_1_seq1-22211MMD-124m_Sample_Tag_Calls.csv", header = T)
VDJ.124M <- read.csv("/Users/xinchen/Desktop/Lyme/data/VDJ/124m_1_seq1-22211MMD-124m_VDJ_perCell.csv", header = T)

ST.124M <- filter(ST.124M, Sample_Tag != "Multiplet") %>% filter(Sample_Tag != "Undetermined")
ST.124M$Subject <- "unknown"
ST.124M$Subject[(str_detect(ST.124M$Sample_Tag, "SampleTag01"))] <- "HD9"
ST.124M$Subject[(str_detect(ST.124M$Sample_Tag, "SampleTag02"))] <- "PTLD1V2"
ST.124M$Subject[(str_detect(ST.124M$Sample_Tag, "SampleTag03"))] <- "RTH2V2"
ST.124M$Subject[(str_detect(ST.124M$Sample_Tag, "SampleTag04"))] <- "RTH2V5_b5"
ST.124M$Subject[(str_detect(ST.124M$Sample_Tag, "SampleTag05"))] <- "RTH2V7"
ST.124M$cell.id <- paste(ST.124M$Subject, ST.124M$Cell_Index, sep="_")
ST.124M.meta <- merge(ST.124M,meta.seurat.TCR, by='cell.id')
VDJ.124M <- filter(VDJ.124M, BCR_Paired_Chains == "TRUE" | TCR_Paired_Chains == "TRUE")
VDJ.124M.meta <- merge(VDJ.124M,ST.124M.meta, by="Cell_Index")
table(VDJ.124M.meta$subject)

####################################################################################################

#load 124n
ST.124n <- read.csv("/Users/xinchen/Desktop/Lyme/data/VDJ/124n_1_seq1-22211MMD-124n_Sample_Tag_Calls.csv", header = T)
VDJ.124n <- read.csv("/Users/xinchen/Desktop/Lyme/data/VDJ/124n_1_seq1-22211MMD-124n_VDJ_perCell.csv", header = T)

ST.124n <- filter(ST.124n, Sample_Tag != "Multiplet") %>% filter(Sample_Tag != "Undetermined")
ST.124n$Subject <- "unknown"
ST.124n$Subject[(str_detect(ST.124n$Sample_Tag, "SampleTag06"))] <- "HD10"
ST.124n$Subject[(str_detect(ST.124n$Sample_Tag, "SampleTag07"))] <- "PTLD2V2"
ST.124n$Subject[(str_detect(ST.124n$Sample_Tag, "SampleTag08"))] <- "PTLD2V3"
ST.124n$Subject[(str_detect(ST.124n$Sample_Tag, "SampleTag09"))] <- "RTH6V3"
ST.124n$Subject[(str_detect(ST.124n$Sample_Tag, "SampleTag10"))] <- "RTH6V5"
ST.124n$cell.id <- paste(ST.124n$Subject, ST.124n$Cell_Index, sep="_")
ST.124n.meta <- merge(ST.124n,meta.seurat.TCR, by='cell.id')
VDJ.124n <- filter(VDJ.124n, BCR_Paired_Chains == "TRUE" | TCR_Paired_Chains == "TRUE")
VDJ.124n.meta <- merge(VDJ.124n,ST.124n.meta, by="Cell_Index")
table(VDJ.124n.meta$subject)

########################################################################################################

#load 124o
ST.124o <- read.csv("/Users/xinchen/Desktop/Lyme/data/VDJ/124o_2_seq1-22211MMD-124o_Sample_Tag_Calls.csv", header = T)
VDJ.124o <- read.csv("/Users/xinchen/Desktop/Lyme/data/VDJ/124o_1_seq1-22211MMD-124o_VDJ_perCell.csv", header = T)

ST.124o <- filter(ST.124o, Sample_Tag != "Multiplet") %>% filter(Sample_Tag != "Undetermined")
ST.124o$Subject <- "unknown"
ST.124o$Subject[(str_detect(ST.124o$Sample_Tag, "SampleTag10"))] <- "HD11"
ST.124o$Subject[(str_detect(ST.124o$Sample_Tag, "SampleTag11"))] <- "PTLD3V2"
ST.124o$Subject[(str_detect(ST.124o$Sample_Tag, "SampleTag12"))] <- "PTLD3V5"
ST.124o$Subject[(str_detect(ST.124o$Sample_Tag, "SampleTag01"))] <- "RTH5V2"
ST.124o$Subject[(str_detect(ST.124o$Sample_Tag, "SampleTag02"))] <- "RTH5V3_b5"
ST.124o$cell.id <- paste(ST.124o$Subject, ST.124o$Cell_Index, sep="_")
ST.124o.meta <- merge(ST.124o,meta.seurat.TCR, by='cell.id')
VDJ.124o <- filter(VDJ.124o, BCR_Paired_Chains == "TRUE" | TCR_Paired_Chains == "TRUE")
VDJ.124o.meta <- merge(VDJ.124o,ST.124o.meta, by="Cell_Index")
table(VDJ.124o.meta$subject)

##############################################################################################################

#load 124P
ST.124P <- read.csv("/Users/xinchen/Desktop/Lyme/data/VDJ/124p_1_seq1-22211MMD-124p_Sample_Tag_Calls.csv", header = T)
VDJ.124P <- read.csv("/Users/xinchen/Desktop/Lyme/data/VDJ/124p_1_seq1-22211MMD-124p_VDJ_perCell.csv", header = T)

ST.124P <- filter(ST.124o, Sample_Tag != "Multiplet") %>% filter(Sample_Tag != "Undetermined")
ST.124P$Subject <- "unknown"
ST.124P$Subject[(str_detect(ST.124P$Sample_Tag, "SampleTag10"))] <- "HD11"
ST.124P$Subject[(str_detect(ST.124P$Sample_Tag, "SampleTag11"))] <- "PTLD3V2"
ST.124P$Subject[(str_detect(ST.124P$Sample_Tag, "SampleTag12"))] <- "PTLD3V5"
ST.124P$Subject[(str_detect(ST.124P$Sample_Tag, "SampleTag01"))] <- "RTH5V2"
ST.124P$Subject[(str_detect(ST.124P$Sample_Tag, "SampleTag02"))] <- "RTH5V3_b5"
ST.124P$cell.id <- paste(ST.124P$Subject, ST.124P$Cell_Index, sep="_")
ST.124P.meta <- merge(ST.124P,meta.seurat.TCR, by='cell.id')
VDJ.124P <- filter(VDJ.124P, BCR_Paired_Chains == "TRUE" | TCR_Paired_Chains == "TRUE")
VDJ.124P.meta <- merge(VDJ.124P,ST.124P.meta, by="Cell_Index")
table(VDJ.124P.meta$subject)

################################################################################################################
#load 124Q
ST.124Q <- read.csv("/Users/xinchen/Desktop/Lyme/data/VDJ/124q_3_seq1-22211MMD-124q_Sample_Tag_Calls.csv", header = T)
VDJ.124Q <- read.csv("/Users/xinchen/Desktop/Lyme/data/VDJ/124q_3_seq1-22211MMD-124q_VDJ_perCell.csv", header = T)

ST.124Q <- filter(ST.124o, Sample_Tag != "Multiplet") %>% filter(Sample_Tag != "Undetermined")
ST.124Q$Subject <- "unknown"
ST.124Q$Subject[(str_detect(ST.124Q$Sample_Tag, "SampleTag10"))] <- "HD11"
ST.124Q$Subject[(str_detect(ST.124Q$Sample_Tag, "SampleTag11"))] <- "PTLD3V2"
ST.124Q$Subject[(str_detect(ST.124Q$Sample_Tag, "SampleTag12"))] <- "PTLD3V5"
ST.124Q$Subject[(str_detect(ST.124Q$Sample_Tag, "SampleTag01"))] <- "RTH5V2"
ST.124Q$Subject[(str_detect(ST.124Q$Sample_Tag, "SampleTag02"))] <- "RTH5V3_b5"
ST.124Q$cell.id <- paste(ST.124Q$Subject, ST.124Q$Cell_Index, sep="_")
ST.124Q.meta <- merge(ST.124Q,meta.seurat.TCR, by='cell.id')
VDJ.124Q <- filter(VDJ.124Q, BCR_Paired_Chains == "TRUE" | TCR_Paired_Chains == "TRUE")
VDJ.124Q.meta <- merge(VDJ.124Q,ST.124Q.meta, by="Cell_Index")
table(VDJ.124Q.meta$subject)

############################################################################################################################
VDJ.124a.meta

#make metadata
#check cell number for each project
list.124a <- table(VDJ.124a.meta$subject) %>% as.data.frame() %>% pull(Var1) %>% as.character()
list.124b <- table(VDJ.124b.meta$subject) %>% as.data.frame() %>% pull(Var1) %>% as.character()
list.124c <- table(VDJ.124c.meta$subject) %>% as.data.frame() %>% pull(Var1) %>% as.character()
list.124d <- table(VDJ.124d.meta$subject) %>% as.data.frame() %>% pull(Var1) %>% as.character()
list.124i <- table(VDJ.124i.meta$subject) %>% as.data.frame() %>% pull(Var1) %>% as.character()
list.124j <- table(VDJ.124j.meta$subject) %>% as.data.frame() %>% pull(Var1) %>% as.character()
list.124l <- table(VDJ.124L.meta$subject) %>% as.data.frame() %>% pull(Var1) %>% as.character()
list.124k <- table(VDJ.124k.meta$subject) %>% as.data.frame() %>% pull(Var1) %>% as.character()
list.124F <- table(VDJ.124F.meta$subject) %>% as.data.frame() %>% pull(Var1) %>% as.character()
list.124M <- table(VDJ.124M.meta$subject) %>% as.data.frame() %>% pull(Var1) %>% as.character()
list.124n <- table(VDJ.124n.meta$subject) %>% as.data.frame() %>% pull(Var1) %>% as.character()
list.124o <- table(VDJ.124o.meta$subject) %>% as.data.frame() %>% pull(Var1) %>% as.character()
list.124P <- table(VDJ.124P.meta$subject) %>% as.data.frame() %>% pull(Var1) %>% as.character()
list.124Q <- table(VDJ.124Q.meta$subject) %>% as.data.frame() %>% pull(Var1) %>% as.character()

immunarch_meta <- c(list.124a,list.124b,list.124c,list.124d,list.124i,list.124j,list.124l,list.124k,list.124F,list.124M,list.124n,list.124o, list.124P, list.124Q) %>% data.frame()
colnames(immunarch_meta) <- "Sample"
immunarch_meta$condition  <- "HD"
immunarch_meta$condition[str_detect(immunarch_meta$Sample, "PTLDN")] <- "PTLDN"
immunarch_meta$condition[str_detect(immunarch_meta$Sample, "PTLD")] <- "PTLD"
immunarch_meta$condition[str_detect(immunarch_meta$Sample, "RTH")] <- "RTH"
immunarch_meta
write.table(immunarch_meta, file = "./immunarch_meta.txt", append = FALSE, sep = "\t", dec = ".", row.names = F, col.names = T)

identical(colnames(VDJ.124a.meta),colnames(VDJ.124k.meta))

VDJall <- rbind(VDJ.124a.meta,VDJ.124b.meta,VDJ.124c.meta,VDJ.124d.meta,VDJ.124i.meta,VDJ.124j.meta,VDJ.124L.meta,VDJ.124k.meta, VDJ.124F.meta, VDJ.124M.meta,VDJ.124n.meta, VDJ.124o.meta, VDJ.124P.meta,VDJ.124Q.meta)

write.table(VDJall, file = "./VDJall.txt", append = FALSE, sep = "\t", dec = ".", row.names = F, col.names = T)

#filter(temp.1, CDR3.nt=="GCTCTGAAACCTCAGGGCGGATCTGAAAAGCTGGTC;GCCAGCAGCCAAGGATTCGGGAGCTCCTACAATGAGCAGTTC")

VDJall <- VDJall %>% mutate(CDR3.nt = paste(TCR_Alpha_Gamma_CDR3_Nucleotide_Dominant,TCR_Beta_Delta_CDR3_Nucleotide_Dominant, sep = ";"))  %>% 
  mutate(CDR3.aa = paste(TCR_Alpha_Gamma_CDR3_Translation_Dominant,TCR_Beta_Delta_CDR3_Translation_Dominant, sep = ";")) %>% 
  mutate(V.name = paste(TCR_Alpha_Gamma_V_gene_Dominant,TCR_Beta_Delta_V_gene_Dominant, sep = ";")) %>% 
  mutate(D.name = TCR_Beta_Delta_D_gene_Dominant) %>% 
  mutate(J.name = paste(TCR_Alpha_Gamma_J_gene_Dominant,TCR_Beta_Delta_J_gene_Dominant, sep = ";")) %>% 
  mutate(chain = paste(TCR_Alpha_Gamma_C_gene_Dominant,TCR_Beta_Delta_C_gene_Dominant, sep = ";"))

for (i in 1:24){
  print(immunarch_meta$Sample[i])
  j <- immunarch_meta$Sample[i]
  seqdepth <- depth$Freq[depth$Var1==j]
  print(seqdepth)
  temp.vdj <- filter(VDJall, subject == j)
  temp.2 <- temp.vdj %>% filter(CDR3.nt != ";") %>% 
    group_by(CDR3.nt) %>% 
    mutate(elementid = as.character(X)) %>%
    mutate(Clones = length(elementid)) %>%
    select(Clones, CDR3.nt, CDR3.aa,V.name,D.name,J.name)
  
  temp.1 <-temp.vdj %>% filter(CDR3.nt != ";") %>% 
    group_by(CDR3.nt) %>% 
    summarise(Barcode = str_c(X, collapse = ","))
  
  temp <- merge(unique(temp.2),temp.1, by="CDR3.nt", all.y=T) %>% 
    mutate(Proportion = round(Clones/seqdepth, 5)) %>% 
    select(Clones,Proportion, CDR3.aa,Barcode)
  temp <- temp[order(-temp$Clones),]
  
  path=paste0("./T_Immunarch/",j,".txt")
  write("# Exported from immunarch 0.5.6 https://immunarch.com",file = path, append = F)
  write.table(temp, file = path, sep = "\t", row.names = F, col.names = T,append = T)
}



SAMPLE.DATA <- repLoad("./T_Immunarch/")
SAMPLE.DATA

repExplore(SAMPLE.DATA$data, "lens") %>% vis()
repClonality(SAMPLE.DATA$data, "rare") %>% vis()
repOverlap(SAMPLE.DATA$data) %>% vis()
geneUsage(SAMPLE.DATA$data[[8]]) %>% vis()
repDiversity(SAMPLE.DATA$data,.method = "Chao1") %>% vis(.by = "condition", .meta = SAMPLE.DATA$meta)
repOverlapAnalysis(SAMPLE.DATA$data, "mds") %>% vis()


imm_ov1 <- repOverlap(SAMPLE.DATA$data, .method = "public", .verbose = F)
imm_ov2 <- repOverlap(SAMPLE.DATA$data, .method = "morisita", .verbose = F)
p1 <- vis(imm_ov1)
p2 <- vis(imm_ov2, .text.size = 2)
p1 + p2


vis(imm_ov1, "heatmap2")

temp.vdj <- filter(VDJall, subject == "HD1")
temp.2 <- temp.vdj %>% filter(CDR3.nt != ";") %>% 
  group_by(CDR3.nt) %>% 
  mutate(elementid = as.character(X)) %>%
  mutate(Clones = length(elementid)) %>%
  select(Clones, CDR3.nt, CDR3.aa,V.name,D.name,J.name,chain)

temp.1 <-temp.vdj %>% filter(CDR3.nt != ";") %>% 
  group_by(CDR3.nt) %>% 
  summarise(Barcode = str_c(X, collapse = ","))

temp <- merge(unique(temp.2),temp.1, by="CDR3.nt", all.y=T) %>% mutate(Proportion = round(Clones / sum(Clones), 5))%>% select(Clones, CDR3.nt, CDR3.aa,V.name,D.name,J.name,chain,Barcode)
temp <- temp[order(-temp$Clones),]


path=paste0("./T_Immunarch/","HD1",".txt")
path

write("# Exported from immunarch 0.5.6 https://immunarch.com",file = path, append = F)
write.table(temp, file = path, sep = "\t", row.names = F, col.names = T,append = T)



str_replace(temp.2$X, '" "', '","')
str_detect(temp.2$elementid, '" "')

count(temp.vdj,CDR3.nt)

table(table(VDJall.c$V.name))

VDJ.124a.meta.c
table(VDJ.124a.meta$Clone)

table(VDJ.124L.meta$immunarch_T_V.name, VDJ.124L.meta$subject)



VDJ.124a.meta.c$immunarch_T_V.name














#split subjects from project

HD1 <- filter(VDJ.124a.meta, subject == "HD1")
PTLDO1V5 <- filter(VDJ.124a.meta, subject == "PTLDO1V5")
RTH1V7 <- filter(VDJ.124a.meta, subject == "RTH1V7")

table(VDJ.124a.meta$subject)
write.csv(file = "./VDJ.HD1.csv",HD1)
write.csv(file = "./VDJ.PTLDO1V5.csv",PTLDO1V5)
write.csv(file = "./VDJ.RTH1V7.csv",RTH1V7)
SL.124a <- c("HD1","RTH1V7","PTLDO1V5")








#select HD1
VDJ.HD1 <- filter(VDJ.124a.meta, subject=="HD1")

#select HD1 T
VDJ.HD1.T <- filter(VDJ.HD1,TCR_Paired_Chains ==T)

#HD1 T cell Alpha Chain
sample.HD1.A <-VDJ.HD1.T %>% select(TCR_Alpha_Gamma_Read_Count,TCR_Alpha_Gamma_CDR3_Nucleotide_Dominant,TCR_Alpha_Gamma_CDR3_Translation_Dominant,
                                    TCR_Alpha_Gamma_V_gene_Dominant, TCR_Alpha_Gamma_J_gene_Dominant, TCR_Alpha_Gamma_C_gene_Dominant)
colnames(sample.HD1.A) <- c("Clones","CDR3.nt","CDR3.aa","V.name",	"J.name", "C.name")
add_column(sample.HD1.A, D.name = "NA", .after ="V.name")

#HD1 T cell Beta Chain
sample.HD1.B <-VDJ.HD1.T %>% select(TCR_Beta_Delta_Read_Count,TCR_Beta_Delta_CDR3_Nucleotide_Dominant,TCR_Beta_Delta_CDR3_Translation_Dominant,
                                    TCR_Beta_Delta_V_gene_Dominant,TCR_Beta_Delta_D_gene_Dominant, TCR_Beta_Delta_J_gene_Dominant, TCR_Beta_Delta_C_gene_Dominant)

colnames(sample.HD1.B) <- c("Clones","CDR3.nt","CDR3.aa","V.name", "D.name", "J.name", "C.name")

#select HD1 B
#HD1 B cell Alpha Chain
#HD1 B cell Beta Chain


VDJ.124a.meta
#select RTH1V7
#select RTH1V7 T cell
#RTH1V7 T cell Alpha Chain
#RTH1V7 T cell Beta Chain
#select RTH1V7 B
#RTH1V7 B cell Alpha Chain
#RTH1V7 B cell Beta Chain


VDJ.124.combined <- rbind(VDJ.124a.meta, VDJ.124b.meta, VDJ.124c.meta, VDJ.124d.meta, VDJ.124F.meta, 
                          VDJ.124h.meta, VDJ.124i.meta, VDJ.124j.meta, VDJ.124k.meta, VDJ.124L.meta)

GLIPHinput <- select(VDJ.124.combined, TCR_Beta_Delta_CDR3_Translation_Dominant,TCR_Beta_Delta_V_gene_Dominant,TCR_Beta_Delta_J_gene_Dominant,
                     TCR_Alpha_Gamma_CDR3_Translation_Dominant,subject,condition,TCR_Alpha_Gamma_Read_Count,TCR_Beta_Delta_Read_Count, CellType.1)

GLIPHinput$count <- GLIPHinput$TCR_Alpha_Gamma_Read_Count + GLIPHinput$TCR_Beta_Delta_Read_Count
GLIPHinput <- GLIPHinput %>% mutate(subject_condition = paste(subject, condition, sep="_")) %>% select(-subject, -condition, -TCR_Alpha_Gamma_Read_Count,-TCR_Beta_Delta_Read_Count)


write.csv(file = "./GLIPHinput1.csv", GLIPHinput)

CD4T <- c("X00_Naive.CD4.T", "X01_Memory.CD4.T", "X26_Naive.CD4.T", "X31_Naive.CD4.T")
GLIPHinput_CD4 <- GLIPHinput %>% filter(CellType.1 %in% CD4T)

CD8T <- unique(GLIPHinput$CellType.1[str_detect(GLIPHinput$CellType.1, "CD8")]) 
GLIPHinput_CD8 <- GLIPHinput %>% filter(CellType.1 %in% CD8T)

write.table(file = "./GLIPH_Input_CD4_raw.txt", GLIPHinput_CD4, sep="\t")
write.table(file = "./GLIPH_Input_CD8_raw.txt", GLIPHinput_CD8, sep="\t")

table(GLIPHinput$CellType.1)






