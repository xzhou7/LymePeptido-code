#Single cell Lyme disease Batch setup
#Author: Xin Chen, Ph.D.
#Date Created: Oct.31.2021
#Date Updated: Dec.8.2021
#analyze CD14 monocyte and predition


#load necessary package
library(data.table)
library(dplyr)
library(Seurat)
library(ggplot2)
library(ggpubr)
library(stringr)
library(patchwork)
library(metap)
library(multtest)
library(Rcpp)
library(ggsci)
library(cowplot)
library(WGCNA)
#set up working directory
setwd("/Users/xinchen/Library/CloudStorage/Box-Box/Xin.Chen.Shareable/R5/")
getwd()

load("/Users/xinchen/Library/CloudStorage/Box-Box/Xin.Chen.Shareable/R5/Data/0131_Step3.5_Cell.Type.2_Combined.RData")



scientific_theme1 <- theme_bw() + theme(plot.title = element_text(face = "bold", size = 12),
                                        legend.background = element_rect(fill = "white", size = 4, colour = "white"),
                                        legend.justification = c(0, 1),
                                        axis.ticks = element_line(colour = "grey70", size = 0.2),
                                        panel.grid.major = element_line(colour = "grey70", size = 0.2),
                                        panel.grid.minor = element_blank())

table(pbmc.clean$condition2)

Mono.list <-c("X14_CD16.Mono","X04_CD14_Mono","X05_CD14_Mono","X06_HLA_Mono","X07_Unk_Mono","X10_Unk_Mono") 

CD14Mono.list <- c("X04_CD14_Mono","X05_CD14_Mono","X06_HLA_Mono","X07_Unk_Mono","X10_Unk_Mono")

Idents(pbmc.clean) <- "condition2"
All.Mono <- subset(pbmc.clean, CellType.1 %in% Mono.list)

Idents(All.Mono) <- "condition"
All.Mono.clean <- subset(All.Mono, condition != "PTLDN")
All.Mono.clean$condition <- factor(All.Mono.clean$condition,levels=c("HD", "RTH", "PTLDO"))
DimPlot(All.Mono.clean, split.by = "condition", group.by = "CellType.1")

#check the percentage of X07 in all subject
bar.data.mono <- data.frame(table(All.Mono.clean$CellType.1, All.Mono.clean$subject))

table(bar.data.mono$condition)
bar.data.clean <- filter(bar.data.mono, Var1 %in% Mono.list)
bar.data.clean$condition <- factor(bar.data.clean$condition,levels=c("HD", "RTH", "PTLD"))
p1.1 <- ggplot(bar.data.clean, aes(fill=Var1, y=Freq, x=Var2)) +  geom_bar(position="fill", stat="identity") # + scale_fill_aaas()
p1.1 <- p1.1 + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylab("Percentage") + xlab("Subject")
p1.1 <- p1.1 + scale_fill_d3()
p1.1

#check the percentage of X07 in condition
bar.data.mono <- data.frame(table(All.Mono.clean$CellType.1, All.Mono.clean$subject))
bar.data.mono$condition <- "normal"
bar.data.mono$condition[str_detect(bar.data.mono$Var2, "HD")] <- "HD"
bar.data.mono$condition[str_detect(bar.data.mono$Var2, "PTLD")] <- "PTLD"
bar.data.mono$condition[str_detect(bar.data.mono$Var2, "RTH")] <- "RTH"

table(bar.data.mono$condition)

bar.data.clean <- filter(bar.data.mono, Var1 %in% Mono.list)
p1.1 <- ggplot(bar.data.clean, aes(fill=Var1, y=Freq, x=condition)) +  geom_bar(position="fill", stat="identity") # + scale_fill_aaas()
p1.1 <- p1.1 + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylab("Percentage") + xlab("Subject")
p1.1 <- p1.1 + scale_fill_d3()
p1.1


#Calculate monocyte in group
mono_percentage <- group_by(bar.data.clean, Var2) %>% mutate(percent = Freq/sum(Freq))
data <- data.frame(group = rep(LETTERS[1:3], each = 4),  # Create example data
                   subgroup = letters[1:4],
                   value = 1:12)
data 

#data_new1 <- transform(data,
#                       perc = ave(value,
#                                  group,
#                                  FUN = prop.table))
data_new1 <- transform(bar.data.clean,
                       perc = ave(Freq,
                                  Var2,
                                  FUN = prop.table))

data_new <- filter(data_new1, Var1 =="X07_Unk_Mono")

counts <- table(mtcars$gear)
barplot(data_new$perc, main="",xlab="condition")

?barplot

p <- filter(data_new1, Var1 %in% c("X07_Unk_Mono")) %>%
  ggplot(aes(x=condition,y=perc, color=condition)) + geom_point() 
p <- p + scientific_theme1 + ylab("% of cluster 7 in CD14 monocytes") +xlab("Condition")
p


p <- filter(data_new1, Var1 %in% c("X07_Unk_Mono")) %>%
  ggplot(aes(x=condition,y=perc, color=condition)) + geom_boxplot()
p <- p + scientific_theme1 + ylab("% of cluster 7 in CD14 monocytes") +xlab("Condition")
p











