#TCR analysis
#created 05-15-2024
#updated 06-17-2025
#https://vdjdb.cdr3.net/search
#install.packages("immunarch")

library(immunarch)
library(tidyverse)
library(ggplot2)
library(stringr)
library(ggpubr)
library(gridExtra)
library(grid)
library(purrr)
library(dplyr)
source("../Code/00_colorKey.R")
TCR_group_colors <- c("1" = "lightgrey",
                      "2-10" = "yellow",
                      "11-30" = "orange",
                      "31-50" = "darkred",
                      "More than 50" = "purple")

TCR_group_colors2 <- c("1" = "lightgrey",
                       "2-20" = "yellow",
                       "More than 20" = "orange")

setwd("/Users/xinchen/Library/CloudStorage/Dropbox/lyme_disease/R7_NR/TCR")
getwd()


immdata <- repLoad("./load_immunarch_data", .mode = "paired")
vdjdb <- read_tsv("./vdjdb/SearchTable-2025-06-17 10_03_31.947.tsv")

names(immdata)
# Typically includes: "data", "meta"
head(immdata$data[[1]])

result_df <- data.frame()
data_list <- list()

for (i in 1:1:length(immdata$data)) {
  # Assign the list element to a new variable
  name <- names(immdata$data)[i]
  assign(name, immdata$data[[i]])
  
  # Add a new column to this variable
  tmp <- get(name)
  tmp$subject <- name
  assign(name, tmp)
  
  cdr3_counts <- table(immdata$data[[i]]$CDR3.aa)
  
  # Categorize the counts into groups
  counts <- as.numeric(cdr3_counts)
  groups <- cut(counts, breaks = c(-Inf, 1, 20, Inf), labels = c("1", "2-20", "More than 20"))
  
  # Count the number of occurrences in each group
  group_counts <- table(groups)
  
  # Create a data frame with these counts
  tmp_df <- as.data.frame(t(group_counts))
  tmp_df$Percentage <- tmp_df$Freq / sum(tmp_df$Freq) * 100
  tmp_df$SubjectID <- name
  
  tmp_df <- tmp_df %>% arrange(desc(groups)) %>% mutate(lab.ypos = cumsum(Percentage) - 0.5 * Percentage)
  
  p.tep <- ggplot(tmp_df, aes(x="", y=Percentage, fill=groups)) +
    geom_bar(stat="identity", width=1,size = 0.1, color = "white") +
    coord_polar(theta = "y") +
    coord_polar("y", start=0) + scale_fill_manual(values = TCR_group_colors2) +
    theme_void() +
    labs(title=name) 
  
  path <- paste0("./pie_chart2/", name, ".pdf")
  ggsave(filename = path,  p.tep, width = 5, height = 4, dpi = 300)
  # Add this data frame to the result data frame
  result_df <- rbind(result_df, tmp_df)
  
  data_list[[i]] <- get(name)
}

combined_immdata <- dplyr::bind_rows(data_list)
combined_colonality <- result_df

combined_colonality$group <- "02.PTLDN"
combined_colonality$group[str_detect(combined_colonality$SubjectID, "HD")] <- "01.HD" 
combined_colonality$group[str_detect(combined_colonality$SubjectID, "PTLDO")] <- "04.PTLD" 
combined_colonality$group[str_detect(combined_colonality$SubjectID, "RTH")] <- "03.RTH" 

combined_colonality$cohort <- "SLICE_I"
combined_colonality$cohort[combined_colonality$group == "04.PTLD"] <- "SLICE_III"
combined_colonality$cohort[combined_colonality$group == "03.RTH"]<- "SLICE_III"

comparisons <- list(c("01.HD", "02.PTLDN"), c("03.RTH", "04.PTLD"))

super_expanded_clone <- combined_colonality %>% filter(groups == "More than 20") %>% 
  ggplot(aes(x=group, y=Percentage, color = group)) + geom_jitter() + geom_boxplot(alpha = 0.7, outlier.alpha = 0) + 
  facet_wrap(.~ cohort, scales = "free_x") + scale_color_manual(values = condition_color_number) + theme_minimal() 



super_expanded_clone <- combined_colonality %>%
  filter(groups == "More than 20") %>%
  ggplot(aes(x = group, y = Percentage, color = group)) +
  geom_jitter(width = 0.2, size = 2) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0, width = 0.6) +
  facet_wrap(. ~ cohort, scales = "free_x") +
  scale_color_manual(values = condition_color_number) +
  stat_compare_means(comparisons = comparisons,
                     method = "wilcox.test",
                     label = "p.signif",
                     size = 5,
                     label.y.npc = "top") +
  theme_minimal() +
  theme(strip.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
super_expanded_clone

#ggsave(filename = "../../Manuscript/Figures/Figure2/Clone_Percent.pdf", super_expanded_clone, height = 3.5, width = 3, dpi = 300)

###########statistics test#############

df <- result_df %>%
  filter(groups == "More than 20") %>%
  select(SubjectID, Percentage)

df$group <- "02.PTLDN"
df$group[str_detect(df$SubjectID, "HD")] <- "01.HD"
df$group[str_detect(df$SubjectID, "PTLDO")] <- "04.PTLD"
df$group[str_detect(df$SubjectID, "RTH")] <- "03.RTH"


df$cohort <- "SLICE_I"
df$cohort[df$group == "04.PTLD"] <- "SLICE_III"
df$cohort[df$group == "03.RTH"] <- "SLICE_III"

wilcox.test(Percentage ~ group, data = df %>% filter(cohort == "SLICE_I" & group %in% c("01.HD", "02.PTLDN")), exact = TRUE)

wilcox.test(Percentage ~ group, data = df %>% filter(cohort == "SLICE_III" & group %in% c("03.RTH", "04.PTLD")), exact = TRUE)

p1 <- wilcox.test(Percentage ~ group, data = df %>% filter(cohort == "SLICE_I" & group %in% c("01.HD", "02.PTLDN")), exact = TRUE)$p.value
p2 <- wilcox.test(Percentage ~ group, data = df %>% filter(cohort == "SLICE_III" & group %in% c("03.RTH", "04.PTLD")), exact = TRUE)$p.value

data.frame(
  Cohort = c("SLICE_I", "SLICE_III"),
  Comparison = c("HD vs PTLDN", "RTH vs PTLD"),
  P_value = c(p1, p2)
)


############


# Create annotation dataframe
# Round to 3 decimal places and format for display
pval_df <- data.frame(
  cohort = c("SLICE_I", "SLICE_III"),
  group1 = c("01.HD", "03.RTH"),
  group2 = c("02.PTLDN", "04.PTLD"),  # will compute dynamically
  p.value = c(p1, p2)
) %>%
  mutate(p.label = sprintf("p = %.3f", p.value))


# Get max y position per group per cohort
top_y <- combined_colonality %>%
  filter(groups == "More than 20") %>%
  group_by(cohort, group) %>%
  summarise(y = max(Percentage, na.rm = TRUE), .groups = "drop")

# Merge to get the top of each group being compared
top_pos <- top_y %>%
  group_by(cohort) %>%
  summarise(y.position = max(y, na.rm = TRUE) * 1.05, .groups = "drop")

# Merge with your p-value annotation table
pval_df <- left_join(pval_df, top_pos, by = "cohort")

super_expanded_clone <- combined_colonality %>%
  filter(groups == "More than 20") %>%
  ggplot(aes(x = group, y = Percentage, color = group)) +
  geom_jitter(width = 0.2, size = 2) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0, width = 0.6) +
  facet_wrap(. ~ cohort, scales = "free_x") +
  scale_color_manual(values = condition_color_number) +
  scale_x_discrete(labels = c(
    "01.HD" = "HD",
    "02.PTLDN" = "PTLDN",
    "03.RTH" = "RTH",
    "04.PTLD" = "PTLD"
  )) +
  stat_pvalue_manual(
    pval_df,
    label = "p.label",
    y.position = "y.position",
    tip.length = 0.01,
    size = 4
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14),
    axis.text.y = element_text(angle = 0, hjust = 0.5, size = 14),
    axis.title.y = element_text(size = 16),
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", size = 0.3),
    axis.ticks = element_line(color = "black", size = 0.3),
    axis.ticks.length = unit(0.2, "cm")
  )

super_expanded_clone

ggsave(filename = "../../Manuscript/Figures/Figure2/Clone_Percent_withPvalTop.pdf",
       plot = super_expanded_clone, height = 4, width = 3.5, dpi = 300)


# TCR analysis
# created 05-15-2024
# updated 07-02-2025
# https://vdjdb.cdr3.net/search

# Load libraries
library(immunarch)
library(tidyverse)
library(ggpubr)
library(gridExtra)
library(grid)

# Source color schemes
source("../Code/00_colorKey.R")

# Define custom TCR group colors
TCR_group_colors2 <- c("1" = "lightgrey", "2-20" = "yellow", "More than 20" = "orange")

# Set working directory
setwd("/Users/xinchen/Library/CloudStorage/Dropbox/lyme_disease/R7_NR/TCR")

# Load immunarch and VDJdb data
immdata <- repLoad("./load_immunarch_data", .mode = "paired")
vdjdb <- read_tsv("./vdjdb/SearchTable-2025-06-17 10_03_31.947.tsv")

# Initialize data containers
result_df <- data.frame()
data_list <- list()

# Loop through samples
for (i in 1:length(immdata$data)) {
  name <- names(immdata$data)[i]
  tmp <- immdata$data[[i]]
  tmp$subject <- name
  
  # Count CDR3s and categorize
  cdr3_counts <- table(tmp$CDR3.aa)
  counts <- as.numeric(cdr3_counts)
  groups <- cut(counts, breaks = c(-Inf, 1, 20, Inf), labels = c("1", "2-20", "More than 20"))
  group_counts <- table(groups)
  
  # Prepare pie chart data
  tmp_df <- as.data.frame(t(group_counts))
  tmp_df$Percentage <- tmp_df$Freq / sum(tmp_df$Freq) * 100
  tmp_df$SubjectID <- name
  tmp_df <- tmp_df %>% arrange(desc(groups)) %>% mutate(lab.ypos = cumsum(Percentage) - 0.5 * Percentage)
  
  # Save pie chart
  p.tep <- ggplot(tmp_df, aes(x="", y=Percentage, fill=groups)) +
    geom_bar(stat="identity", width=1, size=0.1, color="white") +
    coord_polar(theta = "y", start = 0) +
    scale_fill_manual(values = TCR_group_colors2) +
    theme_void() +
    labs(title = name)
  ggsave(filename = paste0("./pie_chart2/", name, ".pdf"), p.tep, width = 5, height = 4, dpi = 300)
  
  result_df <- rbind(result_df, tmp_df)
  data_list[[i]] <- tmp
}

# Combine clonality and assign groups
combined_colonality <- result_df
combined_colonality$group <- "02.PTLDN"
combined_colonality$group[str_detect(combined_colonality$SubjectID, "HD")] <- "01.HD"
combined_colonality$group[str_detect(combined_colonality$SubjectID, "PTLDO")] <- "04.PTLD"
combined_colonality$group[str_detect(combined_colonality$SubjectID, "RTH")] <- "03.RTH"

combined_colonality$cohort <- "SLICE_I"
combined_colonality$cohort[combined_colonality$group == "04.PTLD"] <- "SLICE_III"
combined_colonality$cohort[combined_colonality$group == "03.RTH"] <- "SLICE_III"

# Prepare for Wilcoxon test (exact p-values)
df <- combined_colonality %>%
  filter(groups == "More than 20") %>%
  select(SubjectID, Percentage)

df$group <- "02.PTLDN"
df$group[str_detect(df$SubjectID, "HD")] <- "01.HD"
df$group[str_detect(df$SubjectID, "PTLDO")] <- "04.PTLD"
df$group[str_detect(df$SubjectID, "RTH")] <- "03.RTH"

df$cohort <- "SLICE_I"
df$cohort[df$group == "04.PTLD"] <- "SLICE_III"
df$cohort[df$group == "03.RTH"] <- "SLICE_III"

# Exact Wilcoxon tests
p1 <- wilcox.test(Percentage ~ group, data = df %>% filter(cohort == "SLICE_I" & group %in% c("01.HD", "02.PTLDN")), exact = TRUE)$p.value
p2 <- wilcox.test(Percentage ~ group, data = df %>% filter(cohort == "SLICE_III" & group %in% c("03.RTH", "04.PTLD")), exact = TRUE)$p.value

# Format p-value annotation data
pval_df <- data.frame(
  cohort = c("SLICE_I", "SLICE_III"),
  group1 = c("01.HD", "03.RTH"),
  group2 = c("02.PTLDN", "04.PTLD"),
  p.value = c(p1, p2)
) %>%
  mutate(p.label = sprintf("p = %.3f", p.value))

# Dynamically calculate y.position
top_y <- combined_colonality %>%
  filter(groups == "More than 20") %>%
  group_by(cohort) %>%
  summarise(y.position = max(Percentage, na.rm = TRUE) * 1.05, .groups = "drop")

pval_df <- left_join(pval_df, top_y, by = "cohort")

# Plot
super_expanded_clone <- combined_colonality %>%
  filter(groups == "More than 20") %>%
  ggplot(aes(x = group, y = Percentage, color = group)) +
  geom_jitter(width = 0.2, size = 2) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0, width = 0.6) +
  facet_wrap(. ~ cohort, scales = "free_x") +
  scale_color_manual(values = condition_color_number) +
  stat_pvalue_manual(
    pval_df,
    label = "p.label",
    y.position = "y.position",
    tip.length = 0.01,
    size = 4
  ) +
  theme_minimal() +
  theme(strip.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")
super_expanded_clone

#####
super_expanded_clone <- combined_colonality %>%
  filter(groups == "More than 20") %>%
  ggplot(aes(x = group, y = Percentage, color = group)) +
  geom_jitter(width = 0.2, size = 2) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0, width = 0.6) +
  facet_wrap(. ~ cohort, scales = "free_x") +
  scale_color_manual(values = condition_color_number) +
  stat_pvalue_manual(
    pval_df,
    label = "p.label",
    y.position = "y.position",
    tip.length = 0.01,
    size = 4
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    panel.grid.major = element_blank(),  # remove major grid
    panel.grid.minor = element_blank(),  # remove minor grid
    axis.line = element_line(color = "black", size = 0.3)  # add axis lines
  )

super_expanded_clone



#####
# Save
ggsave(filename = "../../Manuscript/Figures/Figure2/Clone_Percent_withPvalTop.pdf",
       plot = super_expanded_clone, height = 4.5, width = 3.5, dpi = 300)











