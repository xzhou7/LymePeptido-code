# Lyme Disease Single Cell 
# Author: Xin Chen
# Date: 2023-01-05
# To load data
#reference website
#https://stackoverflow.com/questions/13043928/selecting-data-frame-rows-based-on-partial-string-match-in-a-column
#https://stackoverflow.com/questions/38291794/extract-string-before/38295003

#Load necessary library
library(data.table)
library(Seurat)
library(ggplot2)
library(stringr)
library(patchwork)

#setting up working directory
setwd("/Users/xinchen/Library/CloudStorage/Box-Box/Xin.Chen.Shareable/R7_NR")
getwd()

# did not include HD1, PTLD1, RTH1 from 124d
#load data
#HD
data.df.HD1 <- read.csv("~/Desktop/Lyme/data/HD/HD1_Seq1merged124aTG_SampleTag01_hs_DBEC_MolsPerCell_HD0313.csv", header = T, row.names = 1)
#head(data.df.HD1)
seurat.csv.HD1 <- as.data.frame(t(data.df.HD1))

data.df.HD2 <- read.csv("~/Desktop/Lyme/data/HD/HD2-Seq1newmerge124b_SampleTag04_hs_DBEC_MolsPerCell_HD_2219.csv", header = T, row.names = 1)
seurat.csv.HD2 <- as.data.frame(t(data.df.HD2))

data.df.HD3 <- read.csv("~/Desktop/Lyme/data/HD/HD3-seq1ok124k_SampleTag01_hs_DBEC_MolsPerCell_HD_5247.csv", header = T, row.names = 1)
#head(data.df.HD3)
seurat.csv.HD3 <- as.data.frame(t(data.df.HD3))

data.df.HD4 <- read.csv("~/Desktop/Lyme/data/HD/HD4-seq1ok124L_SampleTag08_hs_DBEC_MolsPerCell_HD5249.csv", header = T, row.names = 1)
#head(data.df.HD4)
seurat.csv.HD4 <- as.data.frame(t(data.df.HD4))

data.df.HD5 <- read.csv("~/Desktop/Lyme/data/HD/HD5-seq1xc124f_SampleTag06_hs_DBEC_MolsPerCell_HD1088.csv", header = T, row.names = 1)
#head(data.df.HD5)
seurat.csv.HD5 <- as.data.frame(t(data.df.HD5))

data.df.HD6 <- read.csv("~/Desktop/Lyme/data/HD/HD6-seq1xc124h_SampleTag08_hs_DBEC_MolsPerCell_HD1093.csv", header = T, row.names = 1)
seurat.csv.HD6 <- as.data.frame(t(data.df.HD6))

data.df.HD7 <- read.csv("~/Desktop/Lyme/data/HD/HD7_seq1-124d_SampleTag10_hs_DBEC_MolsPerCell_HD1172.csv", header = T, row.names = 1)
seurat.csv.HD7 <- as.data.frame(t(data.df.HD7))

data.df.HD8 <- read.csv("~/Desktop/Lyme/data/HD/HD8_seq1-124c_SampleTag07_hs_DBEC_MolsPerCell_HD2220.csv", header = T, row.names = 1)
seurat.csv.HD8 <- as.data.frame(t(data.df.HD8))

data.df.HD9 <- read.csv("~/Desktop/Lyme/data/HD/NR_HD9-seq1-22211MMD-124m_SampleTag01_hs_DBEC_MolsPerCell.csv", header = T, row.names = 1)
seurat.csv.HD9 <- as.data.frame(t(data.df.HD9))

data.df.HD10 <- read.csv("~/Desktop/Lyme/data/HD/NR_HD10_seq1-22211MMD-124n_SampleTag06_hs_DBEC_MolsPerCell.csv", header = T, row.names = 1)
seurat.csv.HD10 <- as.data.frame(t(data.df.HD10))

data.df.HD11 <- read.csv("~/Desktop/Lyme/data/HD/NR_HD11_seq1-22211MMD-124o_SampleTag10_hs_DBEC_MolsPerCell.csv", header = T, row.names = 1)
seurat.csv.HD11 <- as.data.frame(t(data.df.HD11))

data.df.HD12 <- read.csv("~/Desktop/Lyme/data/HD/NR_HD12_seq1-22211MMD-124p_SampleTag03_hs_DBEC_MolsPerCell.csv", header = T, row.names = 1)
seurat.csv.HD12 <- as.data.frame(t(data.df.HD12))

#RTH
data.df.RTH1V2 <- read.csv("~/Desktop/Lyme/data/RTH/NR_RTH1V2_seq1-22211MMD-124p_SampleTag07_hs_DBEC_MolsPerCell.csv", header = T, row.names = 1)
seurat.csv.RTH1V2 <- as.data.frame(t(data.df.RTH1V2))

data.df.RTH1V3 <- read.csv("~/Desktop/Lyme/data/RTH/RTH1V3-seq1xc124h_SampleTag09_hs_DBEC_MolsPerCell_RTH_1-071-3.csv", header = T, row.names = 1)
seurat.csv.RTH1V3 <- as.data.frame(t(data.df.RTH1V3))

data.df.RTH1V5 <- read.csv("~/Desktop/Lyme/data/RTH/RTH1V5-seq1xc124h_SampleTag10_hs_DBEC_MolsPerCell_RTH_1-071-5.csv", header = T, row.names = 1)
seurat.csv.RTH1V5 <- as.data.frame(t(data.df.RTH1V5))

data.df.RTH1V7b1 <- read.csv("~/Desktop/Lyme/data/RTH/RTH1V7_Seq1merged124aTG_SampleTag02_hs_DBEC_MolsPerCell_RTH_01_071_7.csv", header = T, row.names = 1)
seurat.csv.RTH1V7b1 <- as.data.frame(t(data.df.RTH1V7b1))

data.df.RTH1V7b5 <- read.csv("~/Desktop/Lyme/data/RTH/NR_RTH1V7_b5_seq1-22211MMD-124q_SampleTag12_hs_DBEC_MolsPerCell.csv", header = T, row.names = 1)
seurat.csv.RTH1V7b5 <- as.data.frame(t(data.df.RTH1V7b5))

data.df.RTH2V2 <- read.csv("~/Desktop/Lyme/data/RTH/NR_RTH2V2_seq1-22211MMD-124m_SampleTag03_hs_DBEC_MolsPerCell.csv", header = T, row.names = 1)
seurat.csv.RTH2V2 <- as.data.frame(t(data.df.RTH2V2))

data.df.RTH2V5b1 <- read.csv("~/Desktop/Lyme/data/RTH/RTH2V5_b1-Seq1newmerge124b_SampleTag05_hs_DBEC_MolsPerCell_RTH_01-082-5.csv", header = T, row.names = 1)
seurat.csv.RTH2V5b1 <- as.data.frame(t(data.df.RTH2V5b1))

data.df.RTH2V5b5 <- read.csv("~/Desktop/Lyme/data/RTH/NR_RTH1V7_b5_seq1-22211MMD-124q_SampleTag12_hs_DBEC_MolsPerCell.csv", header = T, row.names = 1)
seurat.csv.RTH2V5b1 <- as.data.frame(t(data.df.RTH2V5b1))

data.df.RTH2V7 <- read.csv("~/Desktop/Lyme/data/RTH/NR_RTH2V7_seq1-22211MMD-124m_SampleTag05_hs_DBEC_MolsPerCell.csv", header = T, row.names = 1)
seurat.csv.RTH2V7 <- as.data.frame(t(data.df.RTH2V7))

data.df.RTH2V5b5 <- read.csv("~/Desktop/Lyme/data/RTH/NR_RTH2V5_b5_seq1-22211MMD-124m_SampleTag04_hs_DBEC_MolsPerCell.csv", header = T, row.names = 1)
seurat.csv.RTH2V5b5 <- as.data.frame(t(data.df.RTH2V5b5))

data.df.RTH3V3 <- read.csv("~/Desktop/Lyme/data/RTH/RTH3V3-seq1ok124j_SampleTag05_hs_DBEC_MolsPerCell_RTH_01-081-V3.csv", header = T, row.names = 1)
seurat.csv.RTH3V3 <- as.data.frame(t(data.df.RTH3V3))

data.df.RTH3V5 <- read.csv("~/Desktop/Lyme/data/RTH/RTH3V5-seq1ok124L_SampleTag09_hs_DBEC_MolsPerCell_RTH_01-081-V5.csv", header = T, row.names = 1)
seurat.csv.RTH3V5 <- as.data.frame(t(data.df.RTH3V5))

data.df.RTH4V2 <- read.csv("~/Desktop/Lyme/data/RTH/NR_RTH4V2_seq1-22211MMD-124p_SampleTag05_hs_DBEC_MolsPerCell.csv", header = T, row.names = 1)
seurat.csv.RTH4V2 <- as.data.frame(t(data.df.RTH4V2))

data.df.RTH4V3 <- read.csv("~/Desktop/Lyme/data/RTH/RTH4V3-seq1xc124f_SampleTag07_hs_DBEC_MolsPerCell_063_3_visit_RTH_1-063-3.csv", header = T, row.names = 1)
seurat.csv.RTH4V3 <- as.data.frame(t(data.df.RTH4V3))

data.df.RTH4V5 <- read.csv("~/Desktop/Lyme/data/RTH/NR_RTH4V5_seq1-22211MMD-124p_SampleTag06_hs_DBEC_MolsPerCell.csv", header = T, row.names = 1)
seurat.csv.RTH4V5 <- as.data.frame(t(data.df.RTH4V5))

data.df.RTH4V7 <- read.csv("~/Desktop/Lyme/data/RTH/RTH4V7-seq1xc124f_SampleTag08_hs_DBEC_MolsPerCell_063_7_visit_RTH_1-063-7.csv", header = T, row.names = 1)
seurat.csv.RTH4V7 <- as.data.frame(t(data.df.RTH4V7))

data.df.RTH5V2 <- read.csv("~/Desktop/Lyme/data/RTH/NR_RTH5V2_seq1-22211MMD-124o_SampleTag01_hs_DBEC_MolsPerCell.csv", header = T, row.names = 1)
seurat.csv.RTH5V2 <- as.data.frame(t(data.df.RTH5V2))

data.df.RTH5V3b2 <- read.csv("~/Desktop/Lyme/data/RTH/RTH5V3_seq1-124d_SampleTag11_hs_DBEC_MolsPerCell_RTH1-025-3.csv", header = T, row.names = 1)
seurat.csv.RTH5V3b2 <- as.data.frame(t(data.df.RTH5V3b2))

data.df.RTH5V3b5 <- read.csv("~/Desktop/Lyme/data/RTH/NR_RTH5V3_b5_seq1-22211MMD-124o_SampleTag02_hs_DBEC_MolsPerCell.csv", header = T, row.names = 1)
seurat.csv.RTH5V3b5 <- as.data.frame(t(data.df.RTH5V3b5))

data.df.RTH5V7 <- read.csv("~/Desktop/Lyme/data/RTH/NR_RTH5V7_seq1-22211MMD-124q_SampleTag11_hs_DBEC_MolsPerCell.csv", header = T, row.names = 1)
seurat.csv.RTH5V7 <- as.data.frame(t(data.df.RTH5V7))

data.df.RTH6V3 <- read.csv("~/Desktop/Lyme/data/RTH/NR_RTH6V3_seq1-22211MMD-124n_SampleTag09_hs_DBEC_MolsPerCell.csv", header = T, row.names = 1)
seurat.csv.RTH6V3 <- as.data.frame(t(data.df.RTH6V3))

data.df.RTH6V5 <- read.csv("~/Desktop/Lyme/data/RTH/NR_RTH6V5_seq1-22211MMD-124n_SampleTag10_hs_DBEC_MolsPerCell.csv", header = T, row.names = 1)
seurat.csv.RTH6V5 <- as.data.frame(t(data.df.RTH6V5))

data.df.RTH6V7 <- read.csv("~/Desktop/Lyme/data/RTH/NR_RTH6V7_seq1-22211MMD-124q_SampleTag09_hs_DBEC_MolsPerCell.csv", header = T, row.names = 1)
seurat.csv.RTH6V7 <- as.data.frame(t(data.df.RTH6V7))

#PTLDOs
data.df.PTLD1V2 <- read.csv("~/Desktop/Lyme/data/PTLDS/NR_PTLD1V2_seq1-22211MMD-124m_SampleTag02_hs_DBEC_MolsPerCell.csv", header = T, row.names = 1)
seurat.csv.PTLD1V2 <- as.data.frame(t(data.df.PTLD1V2))

data.df.PTLD1V3 <- read.csv("~/Desktop/Lyme/data/PTLDS/PTLD1V3_seq1-124c_SampleTag09_hs_DBEC_MolsPerCell_PTLD_1-023-3.csv", header = T, row.names = 1)
seurat.csv.PTLD1V3 <- as.data.frame(t(data.df.PTLD1V3))

data.df.PTLD1V5 <- read.csv("~/Desktop/Lyme/data/PTLDS/PTLD1V5_Seq1merged124aTG_SampleTag03_hs_DBEC_MolsPerCell_PTLDS_01-023-5.csv", header = T, row.names = 1)
seurat.csv.PTLD1V5 <- as.data.frame(t(data.df.PTLD1V5))

data.df.PTLD2V1 <- read.csv("~/Desktop/Lyme/data/PTLDS/PTLD2V1-Seq1newmerge124b_SampleTag06_hs_DBEC_MolsPerCell_PTLD_01-009-7.csv", header = T, row.names = 1)
seurat.csv.PTLD2V1 <- as.data.frame(t(data.df.PTLD2V1))

data.df.PTLD2V2 <- read.csv("~/Desktop/Lyme/data/PTLDS/NR_PTLD2V2_seq1-22211MMD-124n_SampleTag07_hs_DBEC_MolsPerCell.csv", header = T, row.names = 1)
seurat.csv.PTLD2V2 <- as.data.frame(t(data.df.PTLD2V2))

data.df.PTLD2V3 <- read.csv("~/Desktop/Lyme/data/PTLDS/NR_PTLD2V3_seq1-22211MMD-124n_SampleTag08_hs_DBEC_MolsPerCell.csv", header = T, row.names = 1)
seurat.csv.PTLD2V3 <- as.data.frame(t(data.df.PTLD2V3))

data.df.PTLD2V5 <- read.csv("~/Desktop/Lyme/data/PTLDS/PTLD2V5-seq1xc124f_SampleTag09_hs_DBEC_MolsPerCell_009_PTLD_1-009-5.csv", header = T, row.names = 1)
seurat.csv.PTLD2V5 <- as.data.frame(t(data.df.PTLD2V5))

data.df.PTLD2V7 <- read.csv("~/Desktop/Lyme/data/PTLDS/NR_PTLDO2V7_seq1-22211MMD-124q_SampleTag01_hs_DBEC_MolsPerCell.csv", header = T, row.names = 1)
seurat.csv.PTLD2V7 <- as.data.frame(t(data.df.PTLD2V7))

data.df.PTLD3V2 <- read.csv("~/Desktop/Lyme/data/PTLDS/NR_PTLDO3V2_seq1-22211MMD-124o_SampleTag11_hs_DBEC_MolsPerCell.csv", header = T, row.names = 1)
seurat.csv.PTLD3V2 <- as.data.frame(t(data.df.PTLD3V2))

data.df.PTLD3V3 <- read.csv("~/Desktop/Lyme/data/PTLDS/PTLD3V3-seq1xc124h_SampleTag11_hs_DBEC_MolsPerCell_PTLDS_1-045-3.csv", header = T, row.names = 1)
seurat.csv.PTLD3V3 <- as.data.frame(t(data.df.PTLD3V3))

data.df.PTLD3V5 <- read.csv("~/Desktop/Lyme/data/PTLDS/NR_PTLDO3V5_seq1-22211MMD-124o_SampleTag12_hs_DBEC_MolsPerCell.csv", header = T, row.names = 1)
seurat.csv.PTLD3V5 <- as.data.frame(t(data.df.PTLD3V5))

data.df.PTLD3V7 <- read.csv("~/Desktop/Lyme/data/PTLDS/NR_PTLDO3V7_seq1-22211MMD-124q_SampleTag08_hs_DBEC_MolsPerCell.csv", header = T, row.names = 1)
seurat.csv.PTLD3V7 <- as.data.frame(t(data.df.PTLD3V7))

data.df.PTLD4V3 <- read.csv("~/Desktop/Lyme/data/PTLDS/PTLD4V3_seq1-124d_SampleTag12_hs_DBEC_MolsPerCell_PTLD-1-025-3.csv", header = T, row.names = 1)
seurat.csv.PTLD4V3 <- as.data.frame(t(data.df.PTLD4V3))

data.df.PTLD4V7 <- read.csv("~/Desktop/Lyme/data/PTLDS/NR_PTLDO4V7_seq1-22211MMD-124p_SampleTag04_hs_DBEC_MolsPerCell.csv", header = T, row.names = 1)
seurat.csv.PTLD4V7 <- as.data.frame(t(data.df.PTLD4V7))

#Newbatch PTLD

data.df.PTLDN1 <- read.csv("~/Desktop/Lyme/data/NewPTLDS/PTLDN1-seq1ok124i_SampleTag01_hs_DBEC_MolsPerCell_07-483.csv", header = T, row.names = 1)
seurat.csv.PTLDN1 <- as.data.frame(t(data.df.PTLDN1))

data.df.PTLDN2 <- read.csv("~/Desktop/Lyme/data/NewPTLDS/PTLDN2-seq1ok124i_SampleTag02_hs_DBEC_MolsPerCell_07-484.csv", header = T, row.names = 1)
seurat.csv.PTLDN2 <- as.data.frame(t(data.df.PTLDN2))

data.df.PTLDN3 <- read.csv("~/Desktop/Lyme/data/NewPTLDS/PTLDN3-seq1ok124i_SampleTag03_hs_DBEC_MolsPerCell_07-485.csv", header = T, row.names = 1)
seurat.csv.PTLDN3 <- as.data.frame(t(data.df.PTLDN3))

data.df.PTLDN4 <- read.csv("~/Desktop/Lyme/data/NewPTLDS/PTLDN4-seq1ok124j_SampleTag04_hs_DBEC_MolsPerCell_07-478.csv", header = T, row.names = 1)
seurat.csv.PTLDN4 <- as.data.frame(t(data.df.PTLDN4))

data.df.PTLDN5 <- read.csv("~/Desktop/Lyme/data/NewPTLDS/PTLDN5-seq1ok124k_SampleTag03_hs_DBEC_MolsPerCell_07-472.csv", header = T, row.names = 1)
seurat.csv.PTLDN5 <- as.data.frame(t(data.df.PTLDN5))

data.df.PTLDN6 <- read.csv("~/Desktop/Lyme/data/NewPTLDS/PTLDN6-seq1ok124k_SampleTag04_hs_DBEC_MolsPerCell_07_473.csv", header = T, row.names = 1)
seurat.csv.PTLDN6 <- as.data.frame(t(data.df.PTLDN6))

data.df.PTLDN7 <- read.csv("~/Desktop/Lyme/data/NewPTLDS/PTLDN7-seq1ok124L_SampleTag10_hs_DBEC_MolsPerCell_07-459.csv", header = T, row.names = 1)
seurat.csv.PTLDN7 <- as.data.frame(t(data.df.PTLDN7))

data.df.PTLDN8 <- read.csv("~/Desktop/Lyme/data/NewPTLDS/PTLDN8-seq1ok124L_SampleTag11_hs_DBEC_MolsPerCell-07-479.csv", header = T, row.names = 1)
seurat.csv.PTLDN8 <- as.data.frame(t(data.df.PTLDN8))

#STP
#data.df.STP1 <- read.csv("~/Desktop/Lyme/data/STP/STP1-seq1ok124j_SampleTag06_hs_DBEC_MolsPerCell_STP43.csv", header = T, row.names = 1)
#seurat.csv.STP1 <- as.data.frame(t(data.df.STP1))

#data.df.STP2 <- read.csv("~/Desktop/Lyme/data/STP/STP2-seq1ok124j_SampleTag07_hs_DBEC_MolsPerCell_STP390.csv", header = T, row.names = 1)
#seurat.csv.STP2 <- as.data.frame(t(data.df.STP2))

#data.df.STP3 <- read.csv("~/Desktop/Lyme/data/STP/STP3-seq1xc124f_SampleTag10_hs_DBEC_MolsPerCell_STP277.csv", header = T, row.names = 1)
#seurat.csv.STP3 <- as.data.frame(t(data.df.STP3))

#data.df.STP4 <- read.csv("~/Desktop/Lyme/data/STP/STP4-seq1xc124f_SampleTag11_hs_DBEC_MolsPerCell_STP529.csv", header = T, row.names = 1)
#head(data.df.STP4)
#seurat.csv.STP4 <- as.data.frame(t(data.df.STP4))

#data.df.STP5 <- read.csv("~/Desktop/Lyme/data/STP/STP5-seq1xc124f_SampleTag12_hs_DBEC_MolsPerCell_STP633.csv", header = T, row.names = 1)
#head(data.df.STP5)
#seurat.csv.STP5 <- as.data.frame(t(data.df.STP5))

#data.df.STP6 <- read.csv("~/Desktop/Lyme/data/STP/STP6-seq1xc124h_SampleTag12_hs_DBEC_MolsPerCell_STP386.csv", header = T, row.names = 1)
#seurat.csv.STP6 <- as.data.frame(t(data.df.STP6))

#seperate rna and protein assay HD
#dim(data.df)

seurat.ab.HD1 <- seurat.csv.HD1[rownames(seurat.csv.HD1) %like% "pAbO",]
seurat.rna.HD1 <- seurat.csv.HD1[! rownames(seurat.csv.HD1) %like% "pAbO",]

seurat.ab.HD2 <- seurat.csv.HD2[rownames(seurat.csv.HD2) %like% "pAbO",]
seurat.rna.HD2 <- seurat.csv.HD2[! rownames(seurat.csv.HD2) %like% "pAbO",]

seurat.ab.HD3 <- seurat.csv.HD3[rownames(seurat.csv.HD3) %like% "pAbO",]
seurat.rna.HD3 <- seurat.csv.HD3[! rownames(seurat.csv.HD3) %like% "pAbO",]

seurat.ab.HD4 <- seurat.csv.HD4[rownames(seurat.csv.HD4) %like% "pAbO",]
seurat.rna.HD4 <- seurat.csv.HD4[! rownames(seurat.csv.HD4) %like% "pAbO",]

seurat.ab.HD5 <- seurat.csv.HD5[rownames(seurat.csv.HD5) %like% "pAbO",]
seurat.rna.HD5 <- seurat.csv.HD5[! rownames(seurat.csv.HD5) %like% "pAbO",]

seurat.ab.HD6 <- seurat.csv.HD6[rownames(seurat.csv.HD6) %like% "pAbO",]
seurat.rna.HD6 <- seurat.csv.HD6[! rownames(seurat.csv.HD6) %like% "pAbO",]

seurat.ab.HD7 <- seurat.csv.HD7[rownames(seurat.csv.HD7) %like% "pAbO",]
seurat.rna.HD7 <- seurat.csv.HD7[! rownames(seurat.csv.HD7) %like% "pAbO",]

seurat.ab.HD8 <- seurat.csv.HD8[rownames(seurat.csv.HD8) %like% "pAbO",]
seurat.rna.HD8 <- seurat.csv.HD8[! rownames(seurat.csv.HD8) %like% "pAbO",]

seurat.ab.HD9 <- seurat.csv.HD9[rownames(seurat.csv.HD9) %like% "pAbO",]
seurat.rna.HD9 <- seurat.csv.HD9[! rownames(seurat.csv.HD9) %like% "pAbO",]

seurat.ab.HD10 <- seurat.csv.HD10[rownames(seurat.csv.HD10) %like% "pAbO",]
seurat.rna.HD10 <- seurat.csv.HD10[! rownames(seurat.csv.HD10) %like% "pAbO",]

seurat.ab.HD11 <- seurat.csv.HD11[rownames(seurat.csv.HD11) %like% "pAbO",]
seurat.rna.HD11 <- seurat.csv.HD11[! rownames(seurat.csv.HD11) %like% "pAbO",]

seurat.ab.HD12 <- seurat.csv.HD12[rownames(seurat.csv.HD12) %like% "pAbO",]
seurat.rna.HD12 <- seurat.csv.HD12[! rownames(seurat.csv.HD12) %like% "pAbO",]

#RTH

seurat.ab.RTH1V2 <- seurat.csv.RTH1V2[rownames(seurat.csv.RTH1V2) %like% "pAbO",]
seurat.rna.RTH1V2 <- seurat.csv.RTH1V2[! rownames(seurat.csv.RTH1V2) %like% "pAbO",]

seurat.ab.RTH1V3 <- seurat.csv.RTH1V3[rownames(seurat.csv.RTH1V3) %like% "pAbO",]
seurat.rna.RTH1V3 <- seurat.csv.RTH1V3[! rownames(seurat.csv.RTH1V3) %like% "pAbO",]

seurat.ab.RTH1V5 <- seurat.csv.RTH1V5[rownames(seurat.csv.RTH1V5) %like% "pAbO",]
seurat.rna.RTH1V5 <- seurat.csv.RTH1V5[! rownames(seurat.csv.RTH1V5) %like% "pAbO",]

seurat.ab.RTH1V7b1 <- seurat.csv.RTH1V7b1[rownames(seurat.csv.RTH1V7b1) %like% "pAbO",]
seurat.rna.RTH1V7b1 <- seurat.csv.RTH1V7b1[! rownames(seurat.csv.RTH1V7b1) %like% "pAbO",]

seurat.ab.RTH1V7b5 <- seurat.csv.RTH1V7b5[rownames(seurat.csv.RTH1V7b5) %like% "pAbO",]
seurat.rna.RTH1V7b5 <- seurat.csv.RTH1V7b5[! rownames(seurat.csv.RTH1V7b5) %like% "pAbO",]

seurat.ab.RTH2V2 <- seurat.csv.RTH2V2[rownames(seurat.csv.RTH2V2) %like% "pAbO",]
seurat.rna.RTH2V2 <- seurat.csv.RTH2V2[! rownames(seurat.csv.RTH2V2) %like% "pAbO",]

seurat.ab.RTH2V5b1 <- seurat.csv.RTH2V5b1[rownames(seurat.csv.RTH2V5b1) %like% "pAbO",]
seurat.rna.RTH2V5b1 <- seurat.csv.RTH2V5b1[! rownames(seurat.csv.RTH2V5b1) %like% "pAbO",]

seurat.ab.RTH2V5b5 <- seurat.csv.RTH2V5b5[rownames(seurat.csv.RTH2V5b5) %like% "pAbO",]
seurat.rna.RTH2V5b5 <- seurat.csv.RTH2V5b5[! rownames(seurat.csv.RTH2V5b5) %like% "pAbO",]

seurat.ab.RTH2V7 <- seurat.csv.RTH2V7[rownames(seurat.csv.RTH2V7) %like% "pAbO",]
seurat.rna.RTH2V7 <- seurat.csv.RTH2V7[! rownames(seurat.csv.RTH2V7) %like% "pAbO",]

seurat.ab.RTH3V3 <- seurat.csv.RTH3V3[rownames(seurat.csv.RTH3V3) %like% "pAbO",]
seurat.rna.RTH3V3 <- seurat.csv.RTH3V3[! rownames(seurat.csv.RTH3V3) %like% "pAbO",]

seurat.ab.RTH3V5 <- seurat.csv.RTH3V5[rownames(seurat.csv.RTH3V5) %like% "pAbO",]
seurat.rna.RTH3V5 <- seurat.csv.RTH3V5[! rownames(seurat.csv.RTH3V5) %like% "pAbO",]

seurat.ab.RTH4V2 <- seurat.csv.RTH4V2[rownames(seurat.csv.RTH4V2) %like% "pAbO",]
seurat.rna.RTH4V2 <- seurat.csv.RTH4V2[! rownames(seurat.csv.RTH4V2) %like% "pAbO",]

seurat.ab.RTH4V3 <- seurat.csv.RTH4V3[rownames(seurat.csv.RTH4V3) %like% "pAbO",]
seurat.rna.RTH4V3 <- seurat.csv.RTH4V3[! rownames(seurat.csv.RTH4V3) %like% "pAbO",]

seurat.ab.RTH4V5 <- seurat.csv.RTH4V5[rownames(seurat.csv.RTH4V5) %like% "pAbO",]
seurat.rna.RTH4V5 <- seurat.csv.RTH4V5[! rownames(seurat.csv.RTH4V5) %like% "pAbO",]

seurat.ab.RTH4V7 <- seurat.csv.RTH4V7[rownames(seurat.csv.RTH4V7) %like% "pAbO",]
seurat.rna.RTH4V7 <- seurat.csv.RTH4V7[! rownames(seurat.csv.RTH4V7) %like% "pAbO",]

seurat.ab.RTH5V2 <- seurat.csv.RTH5V2[rownames(seurat.csv.RTH5V2) %like% "pAbO",]
seurat.rna.RTH5V2 <- seurat.csv.RTH5V2[! rownames(seurat.csv.RTH5V2) %like% "pAbO",]

seurat.ab.RTH5V3b2 <- seurat.csv.RTH5V3b2[rownames(seurat.csv.RTH5V3b2) %like% "pAbO",]
seurat.rna.RTH5V3b2 <- seurat.csv.RTH5V3b2[! rownames(seurat.csv.RTH5V3b2) %like% "pAbO",]

seurat.ab.RTH5V3b5 <- seurat.csv.RTH5V3b5[rownames(seurat.csv.RTH5V3b5) %like% "pAbO",]
seurat.rna.RTH5V3b5 <- seurat.csv.RTH5V3b5[! rownames(seurat.csv.RTH5V3b5) %like% "pAbO",]

seurat.ab.RTH5V7 <- seurat.csv.RTH5V7[rownames(seurat.csv.RTH5V7) %like% "pAbO",]
seurat.rna.RTH5V7 <- seurat.csv.RTH5V7[! rownames(seurat.csv.RTH5V7) %like% "pAbO",]

seurat.ab.RTH6V3 <- seurat.csv.RTH6V3[rownames(seurat.csv.RTH6V3) %like% "pAbO",]
seurat.rna.RTH6V3 <- seurat.csv.RTH6V3[! rownames(seurat.csv.RTH6V3) %like% "pAbO",]

seurat.ab.RTH6V5 <- seurat.csv.RTH6V5[rownames(seurat.csv.RTH6V5) %like% "pAbO",]
seurat.rna.RTH6V5 <- seurat.csv.RTH6V5[! rownames(seurat.csv.RTH6V5) %like% "pAbO",]

seurat.ab.RTH6V7 <- seurat.csv.RTH6V7[rownames(seurat.csv.RTH6V7) %like% "pAbO",]
seurat.rna.RTH6V7 <- seurat.csv.RTH6V7[! rownames(seurat.csv.RTH6V7) %like% "pAbO",]

#PTLDs
seurat.ab.PTLD1V2 <- seurat.csv.PTLD1V2[rownames(seurat.csv.PTLD1V2) %like% "pAbO",]
seurat.rna.PTLD1V2 <- seurat.csv.PTLD1V2[! rownames(seurat.csv.PTLD1V2) %like% "pAbO",]

seurat.ab.PTLD1V3 <- seurat.csv.PTLD1V3[rownames(seurat.csv.PTLD1V3) %like% "pAbO",]
seurat.rna.PTLD1V3 <- seurat.csv.PTLD1V3[! rownames(seurat.csv.PTLD1V3) %like% "pAbO",]

seurat.ab.PTLD1V5 <- seurat.csv.PTLD1V5[rownames(seurat.csv.PTLD1V5) %like% "pAbO",]
seurat.rna.PTLD1V5 <- seurat.csv.PTLD1V5[! rownames(seurat.csv.PTLD1V5) %like% "pAbO",]

seurat.ab.PTLD2V1 <- seurat.csv.PTLD2V1[rownames(seurat.csv.PTLD2V1) %like% "pAbO",]
seurat.rna.PTLD2V1 <- seurat.csv.PTLD2V1[! rownames(seurat.csv.PTLD2V1) %like% "pAbO",]

seurat.ab.PTLD2V2 <- seurat.csv.PTLD2V2[rownames(seurat.csv.PTLD2V2) %like% "pAbO",]
seurat.rna.PTLD2V2 <- seurat.csv.PTLD2V2[! rownames(seurat.csv.PTLD2V2) %like% "pAbO",]

seurat.ab.PTLD2V3 <- seurat.csv.PTLD2V3[rownames(seurat.csv.PTLD2V3) %like% "pAbO",]
seurat.rna.PTLD2V3 <- seurat.csv.PTLD2V3[! rownames(seurat.csv.PTLD2V3) %like% "pAbO",]

seurat.ab.PTLD2V5 <- seurat.csv.PTLD2V5[rownames(seurat.csv.PTLD2V5) %like% "pAbO",]
seurat.rna.PTLD2V5 <- seurat.csv.PTLD2V5[! rownames(seurat.csv.PTLD2V5) %like% "pAbO",]

seurat.ab.PTLD2V7 <- seurat.csv.PTLD2V7[rownames(seurat.csv.PTLD2V7) %like% "pAbO",]
seurat.rna.PTLD2V7 <- seurat.csv.PTLD2V7[! rownames(seurat.csv.PTLD2V7) %like% "pAbO",]

seurat.ab.PTLD3V2 <- seurat.csv.PTLD3V2[rownames(seurat.csv.PTLD3V2) %like% "pAbO",]
seurat.rna.PTLD3V2 <- seurat.csv.PTLD3V2[! rownames(seurat.csv.PTLD3V2) %like% "pAbO",]

seurat.ab.PTLD3V3 <- seurat.csv.PTLD3V3[rownames(seurat.csv.PTLD3V3) %like% "pAbO",]
seurat.rna.PTLD3V3 <- seurat.csv.PTLD3V3[! rownames(seurat.csv.PTLD3V3) %like% "pAbO",]

seurat.ab.PTLD3V5 <- seurat.csv.PTLD3V5[rownames(seurat.csv.PTLD3V5) %like% "pAbO",]
seurat.rna.PTLD3V5 <- seurat.csv.PTLD3V5[! rownames(seurat.csv.PTLD3V5) %like% "pAbO",]

seurat.ab.PTLD3V7 <- seurat.csv.PTLD3V7[rownames(seurat.csv.PTLD3V7) %like% "pAbO",]
seurat.rna.PTLD3V7 <- seurat.csv.PTLD3V7[! rownames(seurat.csv.PTLD3V7) %like% "pAbO",]

seurat.ab.PTLD4V3 <- seurat.csv.PTLD4V3[rownames(seurat.csv.PTLD4V3) %like% "pAbO",]
seurat.rna.PTLD4V3 <- seurat.csv.PTLD4V3[! rownames(seurat.csv.PTLD4V3) %like% "pAbO",]

seurat.ab.PTLD4V7 <- seurat.csv.PTLD4V7[rownames(seurat.csv.PTLD4V7) %like% "pAbO",]
seurat.rna.PTLD4V7 <- seurat.csv.PTLD4V7[! rownames(seurat.csv.PTLD4V7) %like% "pAbO",]

#PTLDN
seurat.ab.PTLDN1 <- seurat.csv.PTLDN1[rownames(seurat.csv.PTLDN1) %like% "pAbO",]
seurat.rna.PTLDN1 <- seurat.csv.PTLDN1[! rownames(seurat.csv.PTLDN1) %like% "pAbO",]

seurat.ab.PTLDN2 <- seurat.csv.PTLDN2[rownames(seurat.csv.PTLDN2) %like% "pAbO",]
seurat.rna.PTLDN2 <- seurat.csv.PTLDN2[! rownames(seurat.csv.PTLDN2) %like% "pAbO",]

seurat.ab.PTLDN3 <- seurat.csv.PTLDN3[rownames(seurat.csv.PTLDN3) %like% "pAbO",]
seurat.rna.PTLDN3 <- seurat.csv.PTLDN3[! rownames(seurat.csv.PTLDN3) %like% "pAbO",]

seurat.ab.PTLDN4 <- seurat.csv.PTLDN4[rownames(seurat.csv.PTLDN4) %like% "pAbO",]
seurat.rna.PTLDN4 <- seurat.csv.PTLDN4[! rownames(seurat.csv.PTLDN4) %like% "pAbO",]

seurat.ab.PTLDN5 <- seurat.csv.PTLDN5[rownames(seurat.csv.PTLDN5) %like% "pAbO",]
seurat.rna.PTLDN5 <- seurat.csv.PTLDN5[! rownames(seurat.csv.PTLDN5) %like% "pAbO",]

seurat.ab.PTLDN6 <- seurat.csv.PTLDN6[rownames(seurat.csv.PTLDN6) %like% "pAbO",]
seurat.rna.PTLDN6 <- seurat.csv.PTLDN6[! rownames(seurat.csv.PTLDN6) %like% "pAbO",]

seurat.ab.PTLDN7 <- seurat.csv.PTLDN7[rownames(seurat.csv.PTLDN7) %like% "pAbO",]
seurat.rna.PTLDN7 <- seurat.csv.PTLDN7[! rownames(seurat.csv.PTLDN7) %like% "pAbO",]

seurat.ab.PTLDN8 <- seurat.csv.PTLDN8[rownames(seurat.csv.PTLDN8) %like% "pAbO",]
seurat.rna.PTLDN8 <- seurat.csv.PTLDN8[! rownames(seurat.csv.PTLDN8) %like% "pAbO",]

#STP
#seurat.ab.STP1 <- seurat.csv.STP1[rownames(seurat.csv.STP1) %like% "pAbO",]
#seurat.rna.STP1 <- seurat.csv.STP1[! rownames(seurat.csv.STP1) %like% "pAbO",]

#seurat.ab.STP2 <- seurat.csv.STP2[rownames(seurat.csv.STP2) %like% "pAbO",]
#seurat.rna.STP2 <- seurat.csv.STP2[! rownames(seurat.csv.STP2) %like% "pAbO",]

#seurat.ab.STP3 <- seurat.csv.STP3[rownames(seurat.csv.STP3) %like% "pAbO",]
#seurat.rna.STP3 <- seurat.csv.STP3[! rownames(seurat.csv.STP3) %like% "pAbO",]

#seurat.ab.STP4 <- seurat.csv.STP4[rownames(seurat.csv.STP4) %like% "pAbO",]
#seurat.rna.STP4 <- seurat.csv.STP4[! rownames(seurat.csv.STP4) %like% "pAbO",]

#seurat.ab.STP5 <- seurat.csv.STP5[rownames(seurat.csv.STP5) %like% "pAbO",]
#seurat.rna.STP5 <- seurat.csv.STP5[! rownames(seurat.csv.STP5) %like% "pAbO",]

#seurat.ab.STP6 <- seurat.csv.STP6[rownames(seurat.csv.STP6) %like% "pAbO",]
#seurat.rna.STP6 <- seurat.csv.STP6[! rownames(seurat.csv.STP6) %like% "pAbO",]

#create seurat object using RNA assay and ab assays
#seurat object for DH
pbmc.HD1 <- CreateSeuratObject(counts = seurat.rna.HD1)
pbmc.HD1[["antibody"]] <- CreateAssayObject(counts = seurat.ab.HD1)

pbmc.HD2 <- CreateSeuratObject(counts = seurat.rna.HD2)
pbmc.HD2[["antibody"]] <- CreateAssayObject(counts = seurat.ab.HD2)

pbmc.HD3 <- CreateSeuratObject(counts = seurat.rna.HD3)
pbmc.HD3[["antibody"]] <- CreateAssayObject(counts = seurat.ab.HD3)

pbmc.HD4 <- CreateSeuratObject(counts = seurat.rna.HD4)
pbmc.HD4[["antibody"]] <- CreateAssayObject(counts = seurat.ab.HD4)

pbmc.HD5 <- CreateSeuratObject(counts = seurat.rna.HD5)
pbmc.HD5[["antibody"]] <- CreateAssayObject(counts = seurat.ab.HD5)

pbmc.HD6 <- CreateSeuratObject(counts = seurat.rna.HD6)
pbmc.HD6[["antibody"]] <- CreateAssayObject(counts = seurat.ab.HD6)

pbmc.HD7 <- CreateSeuratObject(counts = seurat.rna.HD7)
pbmc.HD7[["antibody"]] <- CreateAssayObject(counts = seurat.ab.HD7)

pbmc.HD8 <- CreateSeuratObject(counts = seurat.rna.HD8)
pbmc.HD8[["antibody"]] <- CreateAssayObject(counts = seurat.ab.HD8)

pbmc.HD9 <- CreateSeuratObject(counts = seurat.rna.HD9)
pbmc.HD9[["antibody"]] <- CreateAssayObject(counts = seurat.ab.HD9)

pbmc.HD10 <- CreateSeuratObject(counts = seurat.rna.HD10)
pbmc.HD10[["antibody"]] <- CreateAssayObject(counts = seurat.ab.HD10)

pbmc.HD11 <- CreateSeuratObject(counts = seurat.rna.HD11)
pbmc.HD11[["antibody"]] <- CreateAssayObject(counts = seurat.ab.HD11)

pbmc.HD12 <- CreateSeuratObject(counts = seurat.rna.HD12)
pbmc.HD12[["antibody"]] <- CreateAssayObject(counts = seurat.ab.HD12)

#seurat object for PTLD
pbmc.PTLD1V2 <- CreateSeuratObject(counts = seurat.rna.PTLD1V2)
pbmc.PTLD1V2[["antibody"]] <- CreateAssayObject(counts = seurat.ab.PTLD1V2)

pbmc.PTLD1V3 <- CreateSeuratObject(counts = seurat.rna.PTLD1V3)
pbmc.PTLD1V3[["antibody"]] <- CreateAssayObject(counts = seurat.ab.PTLD1V3)

pbmc.PTLD1V5 <- CreateSeuratObject(counts = seurat.rna.PTLD1V5)
pbmc.PTLD1V5[["antibody"]] <- CreateAssayObject(counts = seurat.ab.PTLD1V5)

pbmc.PTLD2V1 <- CreateSeuratObject(counts = seurat.rna.PTLD2V1)
pbmc.PTLD2V1[["antibody"]] <- CreateAssayObject(counts = seurat.ab.PTLD2V1)

pbmc.PTLD2V2 <- CreateSeuratObject(counts = seurat.rna.PTLD2V2)
pbmc.PTLD2V2[["antibody"]] <- CreateAssayObject(counts = seurat.ab.PTLD2V2)

pbmc.PTLD2V3 <- CreateSeuratObject(counts = seurat.rna.PTLD2V3)
pbmc.PTLD2V3[["antibody"]] <- CreateAssayObject(counts = seurat.ab.PTLD2V3)

pbmc.PTLD2V5 <- CreateSeuratObject(counts = seurat.rna.PTLD2V5)
pbmc.PTLD2V5[["antibody"]] <- CreateAssayObject(counts = seurat.ab.PTLD2V5)

pbmc.PTLD2V7 <- CreateSeuratObject(counts = seurat.rna.PTLD2V7)
pbmc.PTLD2V7[["antibody"]] <- CreateAssayObject(counts = seurat.ab.PTLD2V7)

pbmc.PTLD3V2 <- CreateSeuratObject(counts = seurat.rna.PTLD3V2)
pbmc.PTLD3V2[["antibody"]] <- CreateAssayObject(counts = seurat.ab.PTLD3V2)

pbmc.PTLD3V3 <- CreateSeuratObject(counts = seurat.rna.PTLD3V3)
pbmc.PTLD3V3[["antibody"]] <- CreateAssayObject(counts = seurat.ab.PTLD3V3)

pbmc.PTLD3V5 <- CreateSeuratObject(counts = seurat.rna.PTLD3V5)
pbmc.PTLD3V5[["antibody"]] <- CreateAssayObject(counts = seurat.ab.PTLD3V5)

pbmc.PTLD3V7 <- CreateSeuratObject(counts = seurat.rna.PTLD3V7)
pbmc.PTLD3V7[["antibody"]] <- CreateAssayObject(counts = seurat.ab.PTLD3V7)

pbmc.PTLD4V3 <- CreateSeuratObject(counts = seurat.rna.PTLD4V3)
pbmc.PTLD4V3[["antibody"]] <- CreateAssayObject(counts = seurat.ab.PTLD4V3)

pbmc.PTLD4V7 <- CreateSeuratObject(counts = seurat.rna.PTLD4V7)
pbmc.PTLD4V7[["antibody"]] <- CreateAssayObject(counts = seurat.ab.PTLD4V7)

#seurat object for newPTLD
pbmc.PTLDN1 <- CreateSeuratObject(counts = seurat.rna.PTLDN1)
pbmc.PTLDN1[["antibody"]] <- CreateAssayObject(counts = seurat.ab.PTLDN1)

pbmc.PTLDN2 <- CreateSeuratObject(counts = seurat.rna.PTLDN2)
pbmc.PTLDN2[["antibody"]] <- CreateAssayObject(counts = seurat.ab.PTLDN2)

pbmc.PTLDN3 <- CreateSeuratObject(counts = seurat.rna.PTLDN3)
pbmc.PTLDN3[["antibody"]] <- CreateAssayObject(counts = seurat.ab.PTLDN3)

pbmc.PTLDN4 <- CreateSeuratObject(counts = seurat.rna.PTLDN4)
pbmc.PTLDN4[["antibody"]] <- CreateAssayObject(counts = seurat.ab.PTLDN4)

pbmc.PTLDN5 <- CreateSeuratObject(counts = seurat.rna.PTLDN5)
pbmc.PTLDN5[["antibody"]] <- CreateAssayObject(counts = seurat.ab.PTLDN5)

pbmc.PTLDN6 <- CreateSeuratObject(counts = seurat.rna.PTLDN6)
pbmc.PTLDN6[["antibody"]] <- CreateAssayObject(counts = seurat.ab.PTLDN6)

pbmc.PTLDN7 <- CreateSeuratObject(counts = seurat.rna.PTLDN7)
pbmc.PTLDN7[["antibody"]] <- CreateAssayObject(counts = seurat.ab.PTLDN7)

pbmc.PTLDN8 <- CreateSeuratObject(counts = seurat.rna.PTLDN8)
pbmc.PTLDN8[["antibody"]] <- CreateAssayObject(counts = seurat.ab.PTLDN8)

#seurat object for RTH
pbmc.RTH1V2 <- CreateSeuratObject(counts = seurat.rna.RTH1V2)
pbmc.RTH1V2[["antibody"]] <- CreateAssayObject(counts = seurat.ab.RTH1V2)

pbmc.RTH1V3 <- CreateSeuratObject(counts = seurat.rna.RTH1V3)
pbmc.RTH1V3[["antibody"]] <- CreateAssayObject(counts = seurat.ab.RTH1V3)

pbmc.RTH1V5 <- CreateSeuratObject(counts = seurat.rna.RTH1V5)
pbmc.RTH1V5[["antibody"]] <- CreateAssayObject(counts = seurat.ab.RTH1V5)

pbmc.RTH1V7b1 <- CreateSeuratObject(counts = seurat.rna.RTH1V7b1)
pbmc.RTH1V7b1[["antibody"]] <- CreateAssayObject(counts = seurat.ab.RTH1V7b1)

pbmc.RTH1V7b5 <- CreateSeuratObject(counts = seurat.rna.RTH1V7b5)
pbmc.RTH1V7b5[["antibody"]] <- CreateAssayObject(counts = seurat.ab.RTH1V7b5)

pbmc.RTH2V2 <- CreateSeuratObject(counts = seurat.rna.RTH2V2)
pbmc.RTH2V2[["antibody"]] <- CreateAssayObject(counts = seurat.ab.RTH2V2)

pbmc.RTH2V5b1 <- CreateSeuratObject(counts = seurat.rna.RTH2V5b1)
pbmc.RTH2V5b1[["antibody"]] <- CreateAssayObject(counts = seurat.ab.RTH2V5b1)

pbmc.RTH2V5b5 <- CreateSeuratObject(counts = seurat.rna.RTH2V5b5)
pbmc.RTH2V5b5[["antibody"]] <- CreateAssayObject(counts = seurat.ab.RTH2V5b5)

pbmc.RTH2V7 <- CreateSeuratObject(counts = seurat.rna.RTH2V7)
pbmc.RTH2V7[["antibody"]] <- CreateAssayObject(counts = seurat.ab.RTH2V7)

pbmc.RTH3V3 <- CreateSeuratObject(counts = seurat.rna.RTH3V3)
pbmc.RTH3V3[["antibody"]] <- CreateAssayObject(counts = seurat.ab.RTH3V3)

pbmc.RTH3V5 <- CreateSeuratObject(counts = seurat.rna.RTH3V5)
pbmc.RTH3V5[["antibody"]] <- CreateAssayObject(counts = seurat.ab.RTH3V5)

pbmc.RTH4V2 <- CreateSeuratObject(counts = seurat.rna.RTH4V2)
pbmc.RTH4V2[["antibody"]] <- CreateAssayObject(counts = seurat.ab.RTH4V2)

pbmc.RTH4V3 <- CreateSeuratObject(counts = seurat.rna.RTH4V3)
pbmc.RTH4V3[["antibody"]] <- CreateAssayObject(counts = seurat.ab.RTH4V3)

pbmc.RTH4V5 <- CreateSeuratObject(counts = seurat.rna.RTH4V5)
pbmc.RTH4V5[["antibody"]] <- CreateAssayObject(counts = seurat.ab.RTH4V5)

pbmc.RTH4V7 <- CreateSeuratObject(counts = seurat.rna.RTH4V7)
pbmc.RTH4V7[["antibody"]] <- CreateAssayObject(counts = seurat.ab.RTH4V7)

pbmc.RTH5V2 <- CreateSeuratObject(counts = seurat.rna.RTH5V2)
pbmc.RTH5V2[["antibody"]] <- CreateAssayObject(counts = seurat.ab.RTH5V2)

pbmc.RTH5V3b2 <- CreateSeuratObject(counts = seurat.rna.RTH5V3b2)
pbmc.RTH5V3b2[["antibody"]] <- CreateAssayObject(counts = seurat.ab.RTH5V3b2)

pbmc.RTH5V3b5 <- CreateSeuratObject(counts = seurat.rna.RTH5V3b5)
pbmc.RTH5V3b5[["antibody"]] <- CreateAssayObject(counts = seurat.ab.RTH5V3b5)

pbmc.RTH5V7 <- CreateSeuratObject(counts = seurat.rna.RTH5V7)
pbmc.RTH5V7[["antibody"]] <- CreateAssayObject(counts = seurat.ab.RTH5V7)

pbmc.RTH6V3 <- CreateSeuratObject(counts = seurat.rna.RTH6V3)
pbmc.RTH6V3[["antibody"]] <- CreateAssayObject(counts = seurat.ab.RTH6V3)

pbmc.RTH6V5 <- CreateSeuratObject(counts = seurat.rna.RTH6V5)
pbmc.RTH6V5[["antibody"]] <- CreateAssayObject(counts = seurat.ab.RTH6V5)

pbmc.RTH6V7 <- CreateSeuratObject(counts = seurat.rna.RTH6V7)
pbmc.RTH6V7[["antibody"]] <- CreateAssayObject(counts = seurat.ab.RTH6V7)

#seurat object for STP

#pbmc.STP1 <- CreateSeuratObject(counts = seurat.rna.STP1)
#pbmc.STP1[["antibody"]] <- CreateAssayObject(counts = seurat.ab.STP1)

#pbmc.STP2 <- CreateSeuratObject(counts = seurat.rna.STP2)
#pbmc.STP2[["antibody"]] <- CreateAssayObject(counts = seurat.ab.STP2)

#pbmc.STP3 <- CreateSeuratObject(counts = seurat.rna.STP3)
#pbmc.STP3[["antibody"]] <- CreateAssayObject(counts = seurat.ab.STP3)

#pbmc.STP4 <- CreateSeuratObject(counts = seurat.rna.STP4)
#pbmc.STP4[["antibody"]] <- CreateAssayObject(counts = seurat.ab.STP4)

#pbmc.STP5 <- CreateSeuratObject(counts = seurat.rna.STP5)
#pbmc.STP5[["antibody"]] <- CreateAssayObject(counts = seurat.ab.STP5)

#pbmc.STP6 <- CreateSeuratObject(counts = seurat.rna.STP6)
#pbmc.STP6[["antibody"]] <- CreateAssayObject(counts = seurat.ab.STP6)


#merge seuratobject by batch
pbmc.project1 <- merge(pbmc.HD1, y = c(pbmc.RTH1V7b1, pbmc.PTLD1V5, pbmc.HD2,pbmc.RTH2V5b1, pbmc.PTLD2V1), add.cell.ids = c("HD1", "RTH1V7b1","PTLD1V5", "HD2", "RTH2V5b1", "PTLD2V1"), project = "project1")
pbmc.project2 <-merge(pbmc.HD7, y = c(pbmc.RTH5V3b2, pbmc.PTLD4V3, pbmc.HD8,pbmc.RTH6V3, pbmc.PTLD1V3 ), add.cell.ids = c("HD7", "RTH5V3b2", "PTLD4V3", "HD8", "RTH6V3", "PTLD1V3"), project = "project2")
pbmc.project3 <-merge(pbmc.HD5, y = c(pbmc.RTH4V3, pbmc.RTH4V7,pbmc.PTLD2V5), add.cell.ids = c("HD5", "RTH4V3", "RTH4V7", "PTLD2V5"), project = "project3")
pbmc.project4 <-merge(pbmc.HD6, y = c(pbmc.RTH1V3, pbmc.RTH1V5,pbmc.PTLD3V3), add.cell.ids = c("HD6", "RTH1V3", "RTH1V5", "PTLD3V3"), project = "project4")
pbmc.project5 <-merge(pbmc.PTLDN1, y = c(pbmc.PTLDN2,pbmc.PTLDN3,pbmc.PTLDN4, pbmc.RTH3V3), add.cell.ids = c("PTLDN1", "PTLDN2", "PTLDN3","PTLDN4", "RTH3V3"), project = "project5")
pbmc.project6 <-merge(pbmc.HD3, y = c(pbmc.PTLDN5, pbmc.PTLDN6,pbmc.HD4,pbmc.RTH3V5, pbmc.PTLDN7, pbmc.PTLDN8), add.cell.ids = c("HD3", "PTLDN5", "PTLDN6", "HD4","RTH3V5", "PTLDN7", "PTLDN8"), project = "project6")
pbmc.project7 <-merge(pbmc.HD9, y = c(pbmc.PTLD1V2, pbmc.RTH2V2, pbmc.RTH2V5b5, pbmc.RTH2V7, pbmc.HD10, pbmc.PTLD2V2, pbmc.PTLD2V3, pbmc.RTH6V3, pbmc.RTH6V5), add.cell.ids = c("HD9", "PTLD1V2","RTH2V2", "RTH2V5b5","RTH2V7","HD10", "PTLD2V2","PTLD2V3","RTH6V3", "RTH6V5"), project = "project7")
pbmc.project8 <-merge(pbmc.HD11, y = c(pbmc.PTLD3V2, pbmc.PTLD3V5, pbmc.RTH5V2, pbmc.RTH5V3b5, pbmc.HD12, pbmc.PTLD4V7, pbmc.RTH4V2, pbmc.RTH4V5, pbmc.RTH1V2), add.cell.ids = c("HD11", "PTLD3V2", "PTLD3V5", "RTH5V2","RTH5V3b5", "HD12", "PTLD4V7", "RTH4V2", "RTH4V5", "RTH1V2"), project = "project8")
pbmc.project9 <-merge(pbmc.PTLD2V7, y = c(pbmc.PTLD3V7, pbmc.RTH6V7, pbmc.RTH5V7, pbmc.RTH1V7b5), add.cell.ids = c("PTLD2V7", "PTLD3V7", "RTH6V7", "RTH5V7", "RTH1V7b5"), project = "project9")

pbmc <- merge(pbmc.project1, y = c(pbmc.project2,pbmc.project3,pbmc.project4,pbmc.project5, pbmc.project6, pbmc.project7, pbmc.project8, pbmc.project9))
pbmc[[]]


# set condition
pbmc[["condition"]] <-rownames(pbmc[[]])
pbmc$condition[which(str_detect(rownames(pbmc[[]]), "HD"))] <- "HD"
pbmc$condition[which(str_detect(rownames(pbmc[[]]), "RTH"))] <- "RTH"
pbmc$condition[which(str_detect(rownames(pbmc[[]]), "PTLD"))] <- "PTLD"
pbmc$condition[which(str_detect(rownames(pbmc[[]]), "PTLDN"))] <- "PTLDN"

table(pbmc$condition)

#pbmc$condition[which(str_detect(rownames(pbmc[[]]), "STP"))] <- "STP"

# set subject
pbmc$subject <- rownames(pbmc[[]])
pbmc$subject <- sub("\\_.*", "", pbmc[[]]$subject)
pbmc[[]]
table(pbmc$subject)


# to check if object assigned to right batch
table(pbmc$subject[pbmc$condition == "HD"])
table(pbmc$subject[pbmc$condition == "RTH"])
table(pbmc$subject[pbmc$condition == "PTLD"])
table(pbmc$subject[pbmc$condition == "PTLDN"])
#table(pbmc$subject[pbmc$condition == "STP"])

save(pbmc, file = "./Data/Step1.1_010523_load_data.RData")
