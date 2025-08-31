# Load libraries
library(ggplot2)
library(dplyr)
library(ggpubr)
library(multcompView)
library(tidyr)
# mac
setwd("~/Library/CloudStorage/Dropbox/lyme_disease/NULISAseq/XCE2/")
getwd()

data.df <- read.csv("./rawdata/05-13-25 NULISAseq_Inflammation_XinChen_Report.csv", header = TRUE, row.names = 1) %>%
  t() %>% data.frame()
data.df$SampleID <- rownames(data.df)

# Organize metadata
metadata <- read.csv("./rawdata/metadata for Nulisa-seq.csv", header = TRUE)
metadata_Media <- filter(metadata, Condition == "Meida")
metadata_Control <- filter(metadata, Condition == "Control")
metadata_treat <- filter(metadata, Condition != "Meida" & Condition != "Control")
colnames(metadata_treat) <- c("SampleID", "treatment", "donorID", "Dilution")

data.df$SampleID <- sub("([A-Z])(\\d+)\\.(\\d+)", "\\1\\2-\\3", data.df$SampleID)
data.df$SampleID <- sub("\\.(DX\\d+)", "-\\1", data.df$SampleID)

data.treat <- left_join(metadata_treat, data.df, by = "SampleID")
data.long.treat <- data.treat %>% pivot_longer(cols = AGER:WNT7A,
                                               names_to = "cytokine", values_to = "value")

table(data.long.treat$treatment)

# Set working directory for plots
setwd("~/Library/CloudStorage/Dropbox/lyme_disease/Manuscript/Figures/Figure 6-other peptidoglycans/CLD_2")
getwd()

# Load libraries
library(ggplot2)
library(dplyr)
library(ggpubr)
library(multcompView)
library(tidyr)

# Cytokine targets
Cytokine.list <- c("TNF", "IL1B", "CCL8", "CXCL13", "IL17A.IL17F", "IL12B", "IL6", "CSF2", "CCL22", "CHI3L1", "CCL7", "IFNB1")

# Define allowed treatments
treatment_order <- c("None", "Bb_W", "FN_W", "SO_W", "Bb_D", "FN_D", "SO_D", "LPS")

# Clean and restrict treatment values
data.long.treat <- data.long.treat %>%
  mutate(treatment = trimws(treatment)) %>%  # remove spaces
  filter(treatment %in% treatment_order) %>% # keep only specified
  mutate(treatment = factor(treatment, levels = treatment_order))
table(data.long.treat$treatment)


# Number labels only
numeric_labels <- as.character(1:length(treatment_order))
names(numeric_labels) <- treatment_order

# Define color palette
mycolors <- c(
  "None"  = "#999999", "Bb_W"  = "#E69F00", "FN_W"= "#56B4E9", "SO_W"  = "#4DAF4A",
  "Bb_D"  = "#B2182B", "FN_D"  = "#FFD92F", "SO_D"  = "#984EA3","LPS"  = "#CC79A7"
)
cytokine_anova_cld_plot <- function(data, cytokine_name, dilution_filter = "10") {
  data_sub <- data %>%
    filter(Dilution == dilution_filter, cytokine == cytokine_name) %>%
    mutate(treatment = factor(treatment, levels = treatment_order))
  
  aov.test <- aov(value ~ treatment, data = data_sub)
  pval <- summary(aov.test)[[1]][["Pr(>F)"]][1]
  tukey <- TukeyHSD(aov.test)
  tukey_letters <- multcompLetters4(aov.test, tukey)
  
  cld <- as.data.frame.list(tukey_letters$treatment)
  cld$treatment <- rownames(cld)
  rownames(cld) <- NULL
  colnames(cld)[1] <- "Letter"
  
  data_for_plot <- left_join(data_sub, cld, by = "treatment")
  
  letter_positions <- data_for_plot %>%
    group_by(treatment, Letter) %>%
    summarise(y.position = max(value) * 1.05, .groups = 'drop') %>%
    mutate(treatment = factor(treatment, levels = treatment_order))
  
  p <- ggplot(data_sub, aes(x = as.numeric(treatment), y = value)) +  # convert x to numeric
    geom_hline(yintercept = 0, color = "black", size = 0.3) +         # horizontal line at 0
    geom_vline(xintercept = 0, color = "black", size = 0.3) +         # vertical line at 0
    geom_jitter(color = "black", width = 0.15, size = 1.5) +
    geom_boxplot(aes(fill = treatment, group = treatment), alpha = 0.6, 
                 outlier.alpha = 0, size = 0.2) +
    geom_text(data = letter_positions,
              aes(x = as.numeric(treatment), y = y.position, label = Letter),
              size = 4, vjust = 0) +
    labs(title = cytokine_name,
         subtitle = paste0("ANOVA p = ", signif(pval, 3)),
         fill = "Treatment", x = NULL, y = NULL) +
    scale_fill_manual(values = mycolors) +
    scale_x_continuous(breaks = 1:length(treatment_order),
                       labels = numeric_labels,
                       expand = c(0.05,0.05)) +  # avoid box overlapping with y-axis
    theme_minimal() +
    theme(
      legend.position = "none",
      text = element_text(size = 14),
      axis.text.x = element_text(angle = 0, hjust = 1),
      panel.grid = element_blank()
    )
  
  return(p)
}

####
cytokine_anova_cld_plot <- function(data, cytokine_name, dilution_filter = "10") {
  data_sub <- data %>%
    filter(Dilution == dilution_filter, cytokine == cytokine_name) %>%
    mutate(treatment = factor(treatment, levels = treatment_order))
  
  aov.test <- aov(value ~ treatment, data = data_sub)
  pval <- summary(aov.test)[[1]][["Pr(>F)"]][1]
  tukey <- TukeyHSD(aov.test)
  tukey_letters <- multcompLetters4(aov.test, tukey)
  
  cld <- as.data.frame.list(tukey_letters$treatment)
  cld$treatment <- rownames(cld)
  rownames(cld) <- NULL
  colnames(cld)[1] <- "Letter"
  
  data_for_plot <- left_join(data_sub, cld, by = "treatment")
  
  letter_positions <- data_for_plot %>%
    group_by(treatment, Letter) %>%
    summarise(y.position = max(value) * 1.05, .groups = 'drop') %>%
    mutate(treatment = factor(treatment, levels = treatment_order))
  
  p <- ggplot(data_sub, aes(x = treatment, y = value)) +
    geom_boxplot(aes(fill = treatment), alpha = 0.6, outlier.alpha = 0, size = 0.2) +
    geom_jitter(color = "black", width = 0.15, size = 1.5) +
    geom_text(data = letter_positions,
              aes(x = treatment, y = y.position, label = Letter),
              size = 4, vjust = 0) +
    labs(title = cytokine_name,
         subtitle = paste0("ANOVA p = ", signif(pval, 3)),
         fill = "Treatment", x = NULL, y = NULL) +
    scale_fill_manual(values = mycolors) +
    scale_x_discrete(labels = numeric_labels, drop = FALSE) +
    
    # Use classic theme to get full axis with ticks
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      axis.ticks.length = unit(0.15, "cm"),      # control tick length
      axis.line = element_line(size = 0.4),       # axis line thickness
      axis.ticks = element_line(size = 0.4),      # tick thickness
      axis.text.x = element_text(angle = 0, hjust = 1)
    )
  
  return(p)
}


p_test <- cytokine_anova_cld_plot(data.long.treat, "TNF")
print(p_test)


for (Cytokine in Cytokine.list) {
  print(Cytokine)
  pmono <- cytokine_anova_cld_plot(data = data.long.treat, cytokine_name = Cytokine)
  ggsave(filename = paste0("0618_Dilute10_", Cytokine, "_ANOVA_ClassicTicks.pdf"),
         plot = pmono, width = 3, height = 5, dpi = 300)
}

####
# Loop and save plots
for (Cytokine in Cytokine.list) {
  print(Cytokine)
  pmono <- cytokine_anova_cld_plot(data = data.long.treat, cytokine_name = Cytokine)
  print(pmono)
  ggsave(filename = paste0("./update_for_paper/0618_Dilute10_", Cytokine, "_ANOVA_CLDv2.pdf"),
         plot = pmono, width = 3, height = 5, dpi = 300)
}

table(data.long.treat$treatment)


########Yaxis start from 0###

cytokine_anova_cld_plot <- function(data, cytokine_name, dilution_filter = "10") {
  data_sub <- data %>%
    filter(Dilution == dilution_filter, cytokine == cytokine_name) %>%
    mutate(treatment = factor(treatment, levels = treatment_order))
  
  aov.test <- aov(value ~ treatment, data = data_sub)
  pval <- summary(aov.test)[[1]][["Pr(>F)"]][1]
  tukey <- TukeyHSD(aov.test)
  tukey_letters <- multcompLetters4(aov.test, tukey)
  
  cld <- as.data.frame.list(tukey_letters$treatment)
  cld$treatment <- rownames(cld)
  rownames(cld) <- NULL
  colnames(cld)[1] <- "Letter"
  
  data_for_plot <- left_join(data_sub, cld, by = "treatment")
  
  letter_positions <- data_for_plot %>%
    group_by(treatment, Letter) %>%
    summarise(y.position = max(value) * 1.05, .groups = 'drop') %>%
    mutate(treatment = factor(treatment, levels = treatment_order))
  
  p <- ggplot(data_sub, aes(x = treatment, y = value)) +
    geom_boxplot(aes(fill = treatment), alpha = 0.6, outlier.alpha = 0, size = 0.2) +
    geom_jitter(color = "black", width = 0.15, size = 1.5) +
    geom_text(data = letter_positions,
              aes(x = treatment, y = y.position, label = Letter),
              size = 4, vjust = 0) +
    labs(title = cytokine_name,
         subtitle = paste0("ANOVA p = ", signif(pval, 3)),
         fill = "Treatment", x = NULL, y = NULL) +
    scale_fill_manual(values = mycolors) +
    scale_x_discrete(labels = numeric_labels, drop = FALSE) +
    expand_limits(y = 0) +  # <--- Force y-axis to start at 0
    theme_classic(base_size = 14) +
    theme(
      legend.position = "none",
      axis.ticks.length = unit(0.15, "cm"),  # control tick length
      axis.line = element_line(size = 0.4),  # axis line thickness
      axis.ticks = element_line(size = 0.4), # tick thickness
      axis.text.x = element_text(angle = 0, hjust = 1)
    )
  
  return(p)
}

# Loop and save plots
for (Cytokine in Cytokine.list) {
  print(Cytokine)
  pmono <- cytokine_anova_cld_plot(data = data.long.treat, cytokine_name = Cytokine)
  print(pmono)
  ggsave(filename = paste0("./with-Y_zero/0618_Dilute10_", Cytokine, "_ANOVA_CLDv2.pdf"),
         plot = pmono, width = 3, height = 5, dpi = 300)
}





