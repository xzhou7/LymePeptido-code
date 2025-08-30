
setwd(dir = "~/Library/CloudStorage/Box-Box/Xin.Chen.Shareable/R7_NR/TCR")
getwd()


#Install the Immunarch package
#install.packages("immunarch")          

#Load the library
library(immunarch)

#124a
#Separate the VDJ_Dominant_Contigs_AIRR file by sample
total <- read.table(file = "/Users/xinchen/Desktop/Lyme/Raw file/Updated VDJ/124a/124a-final5_VDJ_Dominant_Contigs_AIRR.tsv", sep = "\t", header = T)
demo_stcalls <- read.table("/Users/xinchen/Desktop/Lyme/Raw file/Updated VDJ/124a/124a-final5_Sample_Tag_Calls.csv", header = T, sep = ",")

##Identify which cell indices are assigned to which sample tags
HD1 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag01_hs")]
RTH1V7_b1 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag02_hs")]
PTLDO1V5 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag03_hs")]

#Subset the full AIRR file into sample-specific files
HD1 <- total[which(total$cell_id %in% HD1),]
RTH1V7_b1 <- total[which(total$cell_id %in% RTH1V7_b1),]
PTLDO1V5 <- total[which(total$cell_id %in% PTLDO1V5),]

#Write separated sample AIRR formatted files to a folder (load_immunarch_data)
write.table(x = HD1, file = "./load_immunarch_data/HD1.tsv", sep = "\t", quote = F, row.names = F)
write.table(x = RTH1V7_b1, file = "./load_immunarch_data/RTH1V7_b1.tsv", sep = "\t", quote = F, row.names = F)
write.table(x = PTLDO1V5, file = "./load_immunarch_data/PTLDO1V5.tsv", sep = "\t", quote = F, row.names = F)

#124b
#Separate the VDJ_Dominant_Contigs_AIRR file by sample
total <- read.table(file = "/Users/xinchen/Desktop/Lyme/Raw file/Updated VDJ/124b/124b-final5_VDJ_Dominant_Contigs_AIRR.tsv", sep = "\t", header = T)
demo_stcalls <- read.table("/Users/xinchen/Desktop/Lyme/Raw file/Updated VDJ/124b/124b-final5_Sample_Tag_Calls.csv", header = T, sep = ",")

##Identify which cell indices are assigned to which sample tags
sample1 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag04_hs")]
sample2 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag05_hs")]
sample3 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag06_hs")]

#Subset the full AIRR file into sample-specific files
sample1 <- total[which(total$cell_id %in% sample1),]
sample2 <- total[which(total$cell_id %in% sample2),]
sample3 <- total[which(total$cell_id %in% sample3),]

#Write separated sample AIRR formatted files to a folder (load_immunarch_data)
write.table(x = sample1, file = "./load_immunarch_data/HD2.tsv", sep = "\t", quote = F, row.names = F)
write.table(x = sample2, file = "./load_immunarch_data/RTH2V5_b1.tsv", sep = "\t", quote = F, row.names = F)
write.table(x = sample3, file = "./load_immunarch_data/PTLDO2V1.tsv", sep = "\t", quote = F, row.names = F)

#124C
#Separate the VDJ_Dominant_Contigs_AIRR file by sample
total <- read.table(file = "/Users/xinchen/Desktop/Lyme/Raw file/Updated VDJ/124c/Largest-Input-124c-final6_VDJ_Dominant_Contigs_AIRR.tsv.gz", sep = "\t", header = T)
demo_stcalls <- read.table("/Users/xinchen/Desktop/Lyme/Raw file/Updated VDJ/124c/Largest-Input-124c-final6_Sample_Tag_Calls.csv", header = T, sep = ",")

##Identify which cell indices are assigned to which sample tags
sample1 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag07_hs")]
sample2 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag08_hs")]
sample3 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag09_hs")]

#Subset the full AIRR file into sample-specific files
sample1 <- total[which(total$cell_id %in% sample1),]
sample2 <- total[which(total$cell_id %in% sample2),]
sample3 <- total[which(total$cell_id %in% sample3),]

#Write separated sample AIRR formatted files to a folder (load_immunarch_data)
write.table(x = sample1, file = "./load_immunarch_data/HD8.tsv", sep = "\t", quote = F, row.names = F)
write.table(x = sample2, file = "./load_immunarch_data/RTH6V3.tsv", sep = "\t", quote = F, row.names = F)
write.table(x = sample3, file = "./load_immunarch_data/PTLDO1V3.tsv", sep = "\t", quote = F, row.names = F)

#124D
#Separate the VDJ_Dominant_Contigs_AIRR file by sample
total <- read.table(file = "/Users/xinchen/Desktop/Lyme/Raw file/Updated VDJ/124d/124d-final5_VDJ_Dominant_Contigs_AIRR.tsv", sep = "\t", header = T)
demo_stcalls <- read.table("/Users/xinchen/Desktop/Lyme/Raw file/Updated VDJ/124d/124d-final5_Sample_Tag_Calls.csv", header = T, sep = ",")

##Identify which cell indices are assigned to which sample tags
sample1 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag10_hs")]
sample2 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag11_hs")]
sample3 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag12_hs")]

#Subset the full AIRR file into sample-specific files
sample1 <- total[which(total$cell_id %in% sample1),]
sample2 <- total[which(total$cell_id %in% sample2),]
sample3 <- total[which(total$cell_id %in% sample3),]

#Write separated sample AIRR formatted files to a folder (load_immunarch_data)
write.table(x = sample1, file = "./load_immunarch_data/HD7.tsv", sep = "\t", quote = F, row.names = F)
write.table(x = sample2, file = "./load_immunarch_data/RTH5V3_b2.tsv", sep = "\t", quote = F, row.names = F)
write.table(x = sample3, file = "./load_immunarch_data/PTLDO4V3.tsv", sep = "\t", quote = F, row.names = F)

#124h
#Separate the VDJ_Dominant_Contigs_AIRR file by sample
total <- read.table(file = "/Users/xinchen/Desktop/Lyme/Raw file/Updated VDJ/124h/124h-final5_VDJ_Dominant_Contigs_AIRR.tsv", sep = "\t", header = T)
demo_stcalls <- read.table("/Users/xinchen/Desktop/Lyme/Raw file/Updated VDJ/124h/124h-final5_Sample_Tag_Calls.csv", header = T, sep = ",")

##Identify which cell indices are assigned to which sample tags
sample1 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag08_hs")]
sample2 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag09_hs")]
sample3 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag10_hs")]
sample4 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag11_hs")]
sample5 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag12_hs")]

#Subset the full AIRR file into sample-specific files
sample1 <- total[which(total$cell_id %in% sample1),]
sample2 <- total[which(total$cell_id %in% sample2),]
sample3 <- total[which(total$cell_id %in% sample3),]
sample4 <- total[which(total$cell_id %in% sample4),]
sample5 <- total[which(total$cell_id %in% sample5),]

#Write separated sample AIRR formatted files to a folder (load_immunarch_data)
write.table(x = sample1, file = "./load_immunarch_data/HD6.tsv", sep = "\t", quote = F, row.names = F)
write.table(x = sample2, file = "./load_immunarch_data/RTH1V3.tsv", sep = "\t", quote = F, row.names = F)
write.table(x = sample3, file = "./load_immunarch_data/RTH1V5.tsv", sep = "\t", quote = F, row.names = F)
write.table(x = sample4, file = "./load_immunarch_data/PTLDO3V3.tsv", sep = "\t", quote = F, row.names = F)
write.table(x = sample5, file = "./load_immunarch_data/STP6.tsv", sep = "\t", quote = F, row.names = F)

#124i
#Separate the VDJ_Dominant_Contigs_AIRR file by sample
total <- read.table(file = "/Users/xinchen/Desktop/Lyme/Raw file/Updated VDJ/124i/124i-final5_VDJ_Dominant_Contigs_AIRR.tsv", sep = "\t", header = T)
demo_stcalls <- read.table("/Users/xinchen/Desktop/Lyme/Raw file/Updated VDJ/124i/124i-final5_Sample_Tag_Calls.csv", header = T, sep = ",")

##Identify which cell indices are assigned to which sample tags
sample1 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag01_hs")]
sample2 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag02_hs")]
sample3 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag03_hs")]

#Subset the full AIRR file into sample-specific files
sample1 <- total[which(total$cell_id %in% sample1),]
sample2 <- total[which(total$cell_id %in% sample2),]
sample3 <- total[which(total$cell_id %in% sample3),]

#Write separated sample AIRR formatted files to a folder (load_immunarch_data)
write.table(x = sample1, file = "./load_immunarch_data/PTLDN1.tsv", sep = "\t", quote = F, row.names = F)
write.table(x = sample2, file = "./load_immunarch_data/PTLDN2.tsv", sep = "\t", quote = F, row.names = F)
write.table(x = sample3, file = "./load_immunarch_data/PTLDN3.tsv", sep = "\t", quote = F, row.names = F)

#124j
#Separate the VDJ_Dominant_Contigs_AIRR file by sample
total <- read.table(file = "/Users/xinchen/Desktop/Lyme/Raw file/Updated VDJ/124j/124j-final5_VDJ_Dominant_Contigs_AIRR.tsv", sep = "\t", header = T)
demo_stcalls <- read.table("/Users/xinchen/Desktop/Lyme/Raw file/Updated VDJ/124j/124j-final5_Sample_Tag_Calls.csv", header = T, sep = ",")

##Identify which cell indices are assigned to which sample tags
sample1 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag04_hs")]
sample2 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag05_hs")]
sample3 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag06_hs")]
sample4 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag07_hs")]

#Subset the full AIRR file into sample-specific files
sample1 <- total[which(total$cell_id %in% sample1),]
sample2 <- total[which(total$cell_id %in% sample2),]
sample3 <- total[which(total$cell_id %in% sample3),]
sample4 <- total[which(total$cell_id %in% sample4),]

#Write separated sample AIRR formatted files to a folder (load_immunarch_data)
write.table(x = sample1, file = "./load_immunarch_data/PTLDN4.tsv", sep = "\t", quote = F, row.names = F)
write.table(x = sample2, file = "./load_immunarch_data/RTH3V3.tsv", sep = "\t", quote = F, row.names = F)
write.table(x = sample3, file = "./load_immunarch_data/STP1.tsv", sep = "\t", quote = F, row.names = F)
write.table(x = sample4, file = "./load_immunarch_data/STP2.tsv", sep = "\t", quote = F, row.names = F)

#124k
#Separate the VDJ_Dominant_Contigs_AIRR file by sample
total <- read.table(file = "/Users/xinchen/Desktop/Lyme/Raw file/Updated VDJ/124k/124k-final5_VDJ_Dominant_Contigs_AIRR.tsv", sep = "\t", header = T)
demo_stcalls <- read.table("/Users/xinchen/Desktop/Lyme/Raw file/Updated VDJ/124k/124k-final5_Sample_Tag_Calls.csv", header = T, sep = ",")

##Identify which cell indices are assigned to which sample tags
sample1 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag01_hs")]
sample2 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag03_hs")]
sample3 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag04_hs")]


#Subset the full AIRR file into sample-specific files
sample1 <- total[which(total$cell_id %in% sample1),]
sample2 <- total[which(total$cell_id %in% sample2),]
sample3 <- total[which(total$cell_id %in% sample3),]

#Write separated sample AIRR formatted files to a folder (load_immunarch_data)
write.table(x = sample1, file = "./load_immunarch_data/HD3.tsv", sep = "\t", quote = F, row.names = F)
write.table(x = sample2, file = "./load_immunarch_data/PTLDN5.tsv", sep = "\t", quote = F, row.names = F)
write.table(x = sample3, file = "./load_immunarch_data/PTLDN6.tsv", sep = "\t", quote = F, row.names = F)

#124L
#Separate the VDJ_Dominant_Contigs_AIRR file by sample
total <- read.table(file = "/Users/xinchen/Desktop/Lyme/Raw file/Updated VDJ/124L/124L-final5_VDJ_Dominant_Contigs_AIRR.tsv", sep = "\t", header = T)
demo_stcalls <- read.table("/Users/xinchen/Desktop/Lyme/Raw file/Updated VDJ/124L/124L-final5_Sample_Tag_Calls.csv", header = T, sep = ",")

##Identify which cell indices are assigned to which sample tags
sample1 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag08_hs")]
sample2 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag09_hs")]
sample3 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag10_hs")]
sample4 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag11_hs")]

#Subset the full AIRR file into sample-specific files
sample1 <- total[which(total$cell_id %in% sample1),]
sample2 <- total[which(total$cell_id %in% sample2),]
sample3 <- total[which(total$cell_id %in% sample3),]
sample4 <- total[which(total$cell_id %in% sample4),]

#Write separated sample AIRR formatted files to a folder (load_immunarch_data)
write.table(x = sample1, file = "./load_immunarch_data/HD4.tsv", sep = "\t", quote = F, row.names = F)
write.table(x = sample2, file = "./load_immunarch_data/RTH3V5.tsv", sep = "\t", quote = F, row.names = F)
write.table(x = sample3, file = "./load_immunarch_data/PTLDN7.tsv", sep = "\t", quote = F, row.names = F)
write.table(x = sample3, file = "./load_immunarch_data/PTLDN8.tsv", sep = "\t", quote = F, row.names = F)

#124m
#Separate the VDJ_Dominant_Contigs_AIRR file by sample
total <- read.table(file = "/Users/xinchen/Desktop/Lyme/Raw file/Updated VDJ/124m/xc124m-mRNA-final4_VDJ_Dominant_Contigs_AIRR.tsv", sep = "\t", header = T)
demo_stcalls <- read.table("/Users/xinchen/Desktop/Lyme/Raw file/Updated VDJ/124m/xc124m-mRNA-final4_Sample_Tag_Calls.csv", header = T, sep = ",")

##Identify which cell indices are assigned to which sample tags
sample1 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag01_hs")]
sample2 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag02_hs")]
sample3 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag03_hs")]
sample4 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag04_hs")]
sample5 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag05_hs")]

#Subset the full AIRR file into sample-specific files
sample1 <- total[which(total$cell_id %in% sample1),]
sample2 <- total[which(total$cell_id %in% sample2),]
sample3 <- total[which(total$cell_id %in% sample3),]
sample4 <- total[which(total$cell_id %in% sample4),]
sample5 <- total[which(total$cell_id %in% sample5),]

#Write separated sample AIRR formatted files to a folder (load_immunarch_data)
write.table(x = sample1, file = "./load_immunarch_data/HD9.tsv", sep = "\t", quote = F, row.names = F)
write.table(x = sample2, file = "./load_immunarch_data/PTLD1V2.tsv", sep = "\t", quote = F, row.names = F)
write.table(x = sample3, file = "./load_immunarch_data/RTH2V2.tsv", sep = "\t", quote = F, row.names = F)
write.table(x = sample4, file = "./load_immunarch_data/RTH2V5_b5.tsv", sep = "\t", quote = F, row.names = F)
write.table(x = sample5, file = "./load_immunarch_data/RTH2V7.tsv", sep = "\t", quote = F, row.names = F)


#124n
#Separate the VDJ_Dominant_Contigs_AIRR file by sample
total <- read.table(file = "/Users/xinchen/Desktop/Lyme/Raw file/Updated VDJ/124N/xc124n_VDJ_Dominant_Contigs_AIRR.tsv", sep = "\t", header = T)
demo_stcalls <- read.table("/Users/xinchen/Desktop/Lyme/Raw file/Updated VDJ/124N/xc124n_Sample_Tag_Calls.csv", header = T, sep = ",")

##Identify which cell indices are assigned to which sample tags
sample1 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag06_hs")]
sample2 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag07_hs")]
sample3 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag08_hs")]
sample4 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag09_hs")]
sample5 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag10_hs")]

#Subset the full AIRR file into sample-specific files
sample1 <- total[which(total$cell_id %in% sample1),]
sample2 <- total[which(total$cell_id %in% sample2),]
sample3 <- total[which(total$cell_id %in% sample3),]
sample4 <- total[which(total$cell_id %in% sample4),]
sample5 <- total[which(total$cell_id %in% sample5),]

#Write separated sample AIRR formatted files to a folder (load_immunarch_data)
write.table(x = sample1, file = "./load_immunarch_data/HD10.tsv", sep = "\t", quote = F, row.names = F)
write.table(x = sample2, file = "./load_immunarch_data/PTLD2V2.tsv", sep = "\t", quote = F, row.names = F)
write.table(x = sample3, file = "./load_immunarch_data/PTLD2V3.tsv", sep = "\t", quote = F, row.names = F)
write.table(x = sample4, file = "./load_immunarch_data/RTH6V3.tsv", sep = "\t", quote = F, row.names = F)
write.table(x = sample5, file = "./load_immunarch_data/RTH6V5.tsv", sep = "\t", quote = F, row.names = F)


#124O
#Separate the VDJ_Dominant_Contigs_AIRR file by sample
total <- read.table(file = "/Users/xinchen/Desktop/Lyme/Raw file/Updated VDJ/124O/124o-mRNA-final4_VDJ_Dominant_Contigs_AIRR.tsv", sep = "\t", header = T)
demo_stcalls <- read.table("/Users/xinchen/Desktop/Lyme/Raw file/Updated VDJ/124O/124o-mRNA-final4_Sample_Tag_Calls.csv", header = T, sep = ",")

##Identify which cell indices are assigned to which sample tags
sample1 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag10_hs")]
sample2 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag11_hs")]
sample3 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag12_hs")]
sample4 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag01_hs")]
sample5 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag02_hs")]

#Subset the full AIRR file into sample-specific files
sample1 <- total[which(total$cell_id %in% sample1),]
sample2 <- total[which(total$cell_id %in% sample2),]
sample3 <- total[which(total$cell_id %in% sample3),]
sample4 <- total[which(total$cell_id %in% sample4),]
sample5 <- total[which(total$cell_id %in% sample5),]

#Write separated sample AIRR formatted files to a folder (load_immunarch_data)
write.table(x = sample1, file = "./load_immunarch_data/HD11.tsv", sep = "\t", quote = F, row.names = F)
write.table(x = sample2, file = "./load_immunarch_data/PTLDO3V2.tsv", sep = "\t", quote = F, row.names = F)
write.table(x = sample3, file = "./load_immunarch_data/PTLDO3V5.tsv", sep = "\t", quote = F, row.names = F)
write.table(x = sample4, file = "./load_immunarch_data/RTH5V2.tsv", sep = "\t", quote = F, row.names = F)
write.table(x = sample5, file = "./load_immunarch_data/RTH5V3_b5.tsv", sep = "\t", quote = F, row.names = F)


#124P
#Separate the VDJ_Dominant_Contigs_AIRR file by sample
total <- read.table(file = "/Users/xinchen/Desktop/Lyme/Raw file/Updated VDJ/124P/xc124p_VDJ_Dominant_Contigs_AIRR.tsv", sep = "\t", header = T)
demo_stcalls <- read.table("/Users/xinchen/Desktop/Lyme/Raw file/Updated VDJ/124P/xc124p_Sample_Tag_Calls.csv", header = T, sep = ",")

##Identify which cell indices are assigned to which sample tags
sample1 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag03_hs")]
sample2 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag04_hs")]
sample3 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag05_hs")]
sample4 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag06_hs")]
sample5 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag07_hs")]

#Subset the full AIRR file into sample-specific files
sample1 <- total[which(total$cell_id %in% sample1),]
sample2 <- total[which(total$cell_id %in% sample2),]
sample3 <- total[which(total$cell_id %in% sample3),]
sample4 <- total[which(total$cell_id %in% sample4),]
sample5 <- total[which(total$cell_id %in% sample5),]

#Write separated sample AIRR formatted files to a folder (load_immunarch_data)
write.table(x = sample1, file = "./load_immunarch_data/HD12.tsv", sep = "\t", quote = F, row.names = F)
write.table(x = sample2, file = "./load_immunarch_data/PTLDO4V7.tsv", sep = "\t", quote = F, row.names = F)
write.table(x = sample3, file = "./load_immunarch_data/RTH4V2.tsv", sep = "\t", quote = F, row.names = F)
write.table(x = sample4, file = "./load_immunarch_data/RTH4V5.tsv", sep = "\t", quote = F, row.names = F)
write.table(x = sample5, file = "./load_immunarch_data/RTH1V2.tsv", sep = "\t", quote = F, row.names = F)


#124Q
#Separate the VDJ_Dominant_Contigs_AIRR file by sample
total <- read.table(file = "/Users/xinchen/Desktop/Lyme/Raw file/Updated VDJ/124Q/xc124q-mRNA-final4_VDJ_Dominant_Contigs_AIRR.tsv", sep = "\t", header = T)
demo_stcalls <- read.table("/Users/xinchen/Desktop/Lyme/Raw file/Updated VDJ/124Q/xc124q-mRNA-final4_Sample_Tag_Calls.csv", header = T, sep = ",")

##Identify which cell indices are assigned to which sample tags
sample1 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag01_hs")]
sample2 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag08_hs")]
sample3 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag09_hs")]
sample4 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag11_hs")]
sample5 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag12_hs")]

#Subset the full AIRR file into sample-specific files
sample1 <- total[which(total$cell_id %in% sample1),]
sample2 <- total[which(total$cell_id %in% sample2),]
sample3 <- total[which(total$cell_id %in% sample3),]
sample4 <- total[which(total$cell_id %in% sample4),]
sample5 <- total[which(total$cell_id %in% sample5),]

#Write separated sample AIRR formatted files to a folder (load_immunarch_data)
write.table(x = sample1, file = "./load_immunarch_data/PTLDO2V7.tsv", sep = "\t", quote = F, row.names = F)
write.table(x = sample2, file = "./load_immunarch_data/PTLDO3V7.tsv", sep = "\t", quote = F, row.names = F)
write.table(x = sample3, file = "./load_immunarch_data/RTH6V7.tsv", sep = "\t", quote = F, row.names = F)
write.table(x = sample4, file = "./load_immunarch_data/RTH5V7.tsv", sep = "\t", quote = F, row.names = F)
write.table(x = sample5, file = "./load_immunarch_data/RTH1V7_b5.tsv", sep = "\t", quote = F, row.names = F)

#124f
#Separate the VDJ_Dominant_Contigs_AIRR file by sample
total <- read.table(file = "/Users/xinchen/Desktop/Lyme/Raw file/Updated VDJ/124f/124f-final5_VDJ_Dominant_Contigs_AIRR.tsv", sep = "\t", header = T)
demo_stcalls <- read.table("/Users/xinchen/Desktop/Lyme/Raw file/Updated VDJ/124f/124f-final5_Sample_Tag_Calls.csv", header = T, sep = ",")

##Identify which cell indices are assigned to which sample tags
sample1 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag06_hs")]
sample2 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag07_hs")]
sample3 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag08_hs")]
sample4 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag09_hs")]
sample5 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag10_hs")]
sample6 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag11_hs")]
sample7 <- demo_stcalls$Cell_Index[which(demo_stcalls$Sample_Name == "SampleTag12_hs")]

#Subset the full AIRR file into sample-specific files
sample1 <- total[which(total$cell_id %in% sample1),]
sample2 <- total[which(total$cell_id %in% sample2),]
sample3 <- total[which(total$cell_id %in% sample3),]
sample4 <- total[which(total$cell_id %in% sample4),]
sample5 <- total[which(total$cell_id %in% sample5),]
sample6 <- total[which(total$cell_id %in% sample6),]
sample7 <- total[which(total$cell_id %in% sample7),]

#Write separated sample AIRR formatted files to a folder (load_immunarch_data)
write.table(x = sample1, file = "./load_immunarch_data/HD5.tsv", sep = "\t", quote = F, row.names = F)
write.table(x = sample2, file = "./load_immunarch_data/RTH4V3.tsv", sep = "\t", quote = F, row.names = F)
write.table(x = sample3, file = "./load_immunarch_data/RTH4V7.tsv", sep = "\t", quote = F, row.names = F)
write.table(x = sample4, file = "./load_immunarch_data/PTLDO2V5.tsv", sep = "\t", quote = F, row.names = F)
write.table(x = sample5, file = "./load_immunarch_data/STP3.tsv", sep = "\t", quote = F, row.names = F)
write.table(x = sample6, file = "./load_immunarch_data/STP4.tsv", sep = "\t", quote = F, row.names = F)
write.table(x = sample7, file = "./load_immunarch_data/STP5.tsv", sep = "\t", quote = F, row.names = F)

#Make metadata file for loading in data
meta <- data.frame(Sample = c("HD1", "RTH1V7_b1", "PTLDO1V5","HD2", "RTH2V5_b1", "PTLDO2V1", "HD8", "RTH6V3", 
                              "PTLDO1V3","HD7", "RTH5V3_b2", "PTLDO4V3", "HD6", "RTH1V3", "RTH1V5", "PTLDO3V3", 
                              "STP6", "PTLDN1", "PTLDN2", "PTLDN3", "PTLDN4", "RTH3V3", "STP1", "STP2", "HD5", 
                              "PTLDN5", "PTLDN6", "HD4", "RTH3V5", "PTLDN7", "PTLDN8", "HD9", "PTLD1V2", "RTH2V2", 
                              "RTH2V5_b5", "RTH2V7", "HD10", "PTLD2V2", "RTH6V3", "RTH6V5", "HD11", "PTLDO3V2", 
                              "PTLDO3V5", "RTH5V2", "RTH5V3_b5", "HD11", "PTLDO3V2", "PTLDO3V5", "RTH5V2", "RTH5V3_b5",
                              "PTLDO2V7", "PTLDO3V7", "RTH6V7", "RTH5V7", "RTH1V7_b5", "HD5","RTH4V3","RTH4V7","PTLDO2V5",
                              "STP3", "STP4", "STP5"), stringsAsFactors = F)

meta$Treatment <- c("HD", "RTH","PTLD","HD", "RTH","PTLD","HD", "RTH",
                    "PTLD","HD", "RTH","PTLD", "HD", "RTH", "RTH", "PTLD",
                    "STP", "PTLDN", "PTLDN", "PTLDN", "PTLDN", "RTH", "STP", "STP", "HD", 
                    "PTLDN", "PTLDN", "HD", "RTH", "PTLDN", "PTLDN", "HD", "PTLD", "RTH", 
                    "RTH","RTH", "HD","PTLD", "PTLD", "RTH", "HD","PTLD",
                     "PTLD", "RTH", "RTH", "HD", "PTLD", "PTLD", "RTH", "RTH", 
                    "PTLD", "PTLD", "RTH", "RTH", "RTH", "HD", "RTH", "RTH", "PTLD",
                    "STP", "STP", "STP")

meta <- as_tibble(meta)
meta <- write.table(meta, file = "./load_immunarch_data/metadata.txt", row.names = F, quote = F, col.names = T, sep = "\t")

#Load the data
immdata <- repLoad("./load_immunarch_data")

#Explore the data
repExplore(immdata$data, "lens") %>% vis()  # Visualise the length distribution of CDR3
repClonality(immdata$data, "homeo") %>% vis()  # Visualise the relative abundance of clonotypes
repDiversity(immdata$data) %>% vis(.meta = immdata$meta)  # Visualise the Chao1 diversity of repertoires
geneUsage(immdata$data) %>% vis()

################################################################################################################################################
#Vignettes for further analysis:

#Repertoire Overlap: https://immunarch.com/articles/web_only/v4_overlap.html
#Gene Usage: https://immunarch.com/articles/web_only/v5_gene_usage.html
#Diversity Estimation: https://immunarch.com/articles/web_only/v6_diversity.html
#Tracking clonotypes across time points: https://immunarch.com/articles/web_only/v8_tracking.html
#Kmer and sequence motif analysis and visualization: https://immunarch.com/articles/web_only/v9_kmers.html

