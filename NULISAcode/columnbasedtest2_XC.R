#Cytokine analysis for temp and treatment 
#Author: Xin Chen, Ph.D. 
#Date: 2025-05-24
#Last update: 2025-05-24

library(stringr)
library(tidyverse)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(corrr)
library(factoextra)
library(FactoMineR)
library(patchwork)
library(reshape2)
library(ggrepel)

#mac
setwd("~/Library/CloudStorage/Dropbox/lyme_disease/NULISAseq/XCE2/")
getwd()


###############################Bargraph betwen ND, Bb, and LPS, CHECK IF THERE ARE MORE CYTOKINES ARE DIFERENTIAL BETWEEN THEM
data.df <- read.csv("./rawdata/05-13-25 NULISAseq_Inflammation_XinChen_Report.csv", header = T, row.names = 1) %>% t() %>% data.frame()
data.df$SampleID <- rownames(data.df)

#orgnize metadata
metadata <- read.csv("./rawdata/metadata for Nulisa-seq.csv", header = T)

#metadata_temp <- filter(metadata_separated, separated_4 != "experiments2") %>% filter(separated_4 != "controls")
metadata_Media <- filter(metadata, Condition == "Meida")
metadata_Control <- filter(metadata, Condition == "Control")
metadata_treat <- filter(metadata, Condition != "Meida" & Condition != "Control")
metadata_treat

#colnames(metadata_temp) <- c("SampleID", "separated_1", "donorID", "treatment", "temperature")
colnames(metadata_treat) <- c("SampleID",  "treatment", "donorID", "Dilution")

data.df$SampleID <- sub("([A-Z])(\\d+)\\.(\\d+)", "\\1\\2-\\3", data.df$SampleID)
data.df$SampleID <- sub("\\.(DX\\d+)", "-\\1", data.df$SampleID)

data.treat <- left_join(metadata_treat, data.df, by=c("SampleID"="SampleID"))
data.treat
data.long.treat <- data.treat %>% pivot_longer(cols = AGER:WNT7A,
    names_to = "cytokine", values_to = "value")

table(data.long.treat$treatment)


data.long.treat$treatment <- factor(data.long.treat$treatment, levels = c("None", "LPS","Bb_D","Bb_W","FN_D","FN_W",
                                                                          "SO_D", "SO_W", "SP_W","SP_D","VP_W","VP_D"))
data.long.treat$donorID <- factor(data.long.treat$donorID)
  
#Filter to Dilution 10
df_dil10 <- data.long.treat %>%
  filter(Dilution == 10)

########################################################################################XC
library(dplyr)

df_filtered <- df_dil10 %>%
  filter(treatment %in% c("None", "Bb_W", "LPS"))

library(purrr)

# Function to run test for one cytokine
test_cytokine <- function(data, test_group) {
  subset <- data %>% filter(treatment %in% c("None", test_group))
  
  results <- subset %>%
    group_by(cytokine) %>%
    summarise(
      p_value = tryCatch(wilcox.test(value ~ treatment)$p.value, error = function(e) NA),
      .groups = "drop"
    ) %>%
    mutate(treatment = test_group)
  
  return(results)
}

# Run tests for both Bb_W and LPS
res_bbw <- test_cytokine(df_filtered, "Bb_W")
res_lps <- test_cytokine(df_filtered, "LPS")

# Combine results
de_results <- bind_rows(res_bbw, res_lps) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  filter(!is.na(p_value))

sig_cytokines <- de_results %>%
  filter(p_adj < 0.05) %>%
  pull(cytokine) %>%
  unique()

library(ggplot2)

df_plot <- df_filtered %>%
  filter(cytokine %in% sig_cytokines)

ggplot(df_plot, aes(x = treatment, y = value, fill = treatment)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.2, size = 1, alpha = 0.6) +
  facet_wrap(~cytokine, scales = "free_y") +
  theme_minimal() +
  labs(title = "Differentially Expressed Cytokines (p.adj < 0.05)",
       y = "Expression Value") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

########################################################################################XC

#Summarize by treatment × cytokine
cytokine_summary <- df_dil10 %>%
  group_by(treatment, cytokine) %>%
  summarise(mean_value = mean(value, na.rm = TRUE), .groups = "drop")

#Pivot to wide format for heatmap
cytokine_wide <- cytokine_summary %>%
  pivot_wider(names_from = cytokine, values_from = mean_value) %>%
  column_to_rownames("treatment") %>%
  as.data.frame()
cytokine_matrix <- as.matrix(cytokine_wide)

#select for None/Bb/LPS



#optimal: Z-score normalize across treatments
cytokine_matrix_scaled <- scale(cytokine_matrix)

#1) remove rolw with NA
cytokine_matrix_clean <- cytokine_matrix_scaled[complete.cases(cytokine_matrix_scaled), ]
cytokine_matrix_clean <- cytokine_matrix_scaled[, colSums(is.na(cytokine_matrix_scaled)) == 0]

#Option 2: Replace NAs with 0 or a small value
cytokine_matrix_clean <- cytokine_matrix_scaled
cytokine_matrix_clean[is.na(cytokine_matrix_clean)] <- 0
# Remove rows where all values are zero
cytokine_matrix_nozero <- cytokine_matrix[rowSums(cytokine_matrix == 0) < ncol(cytokine_matrix), ]
cytokine_matrix_scaled <- scale(cytokine_matrix)
#Heatmap to compare cytokine profiles across treatments
pheatmap(cytokine_matrix_nozero,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Cytokine Induction Across Treatments (Dilution 10)")



cytokine_matrix_nozero_cols <- cytokine_matrix[, colSums(cytokine_matrix == 0) < nrow(cytokine_matrix)]
# Remove all-zero rows
cytokine_matrix <- cytokine_matrix[rowSums(cytokine_matrix == 0) < ncol(cytokine_matrix), ]

# Remove all-zero columns
cytokine_matrix <- cytokine_matrix[, colSums(cytokine_matrix == 0) < nrow(cytokine_matrix)]
cytokine_matrix_scaled <- scale(cytokine_matrix)
pheatmap(cytokine_matrix_scaled,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Cytokine Induction Across Treatments (Dilution 10)")


#highlight autoimmune cytokines:
autoimmune_monocyte_cytokines <- c("IL1B", "IL6", "TNF", "IL12B", "IL23A", 
                                   "CXCL10", "CCL2", "CCL3", "CCL4", "TNFSF13B")

# Only keep those cytokines that exist in your matrix
available_cytokines <- intersect(autoimmune_monocyte_cytokines, colnames(cytokine_matrix_clean))

cytokine_autoimmune_subset <- cytokine_matrix_clean[, available_cytokines]

#highly in a heatmap
pheatmap(cytokine_autoimmune_subset,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Autoimmune-Related Monocyte Cytokines (Dilution 10)")




########################################################################################XC
library(ggplot2)
library(dplyr)
library(ggpubr)
Cytokine.list <- c("CCL22", "CHI3L1", "CCL8", "IFNB1")

for (Cytokine in Cytokine.list){
  print(Cytokine)
  
  # Subset your data for each cytokine
  data_sub <- data.long.treat %>% filter(Dilution == "10", cytokine == Cytokine)
  
  # One-way ANOVA
  aov.test <- aov(value ~ treatment, data = data_sub)
  summary.aov <- summary(aov.test)
  pval <- summary.aov[[1]][["Pr(>F)"]][1]
  
  # Tukey post-hoc test
  tukey <- TukeyHSD(aov.test)
  tukey_df <- as.data.frame(tukey$treatment)
  tukey_df$Comparison <- rownames(tukey_df)
  rownames(tukey_df) <- NULL
  
  # Split the comparison string to get group1 and group2
  tukey_df <- tukey_df %>%
    separate(Comparison, into = c("group1", "group2"), sep = "-")
  
  # Filter significant comparisons only
  tukey_sig <- tukey_df %>% filter(`p adj` < 0.05)
  
  # If no significant pairs, skip p-value annotation
  if(nrow(tukey_sig) > 0){
    # Create y.position to avoid label overlap
    y_max <- max(data_sub$value, na.rm = TRUE)
    tukey_sig <- tukey_sig %>%
      arrange(`p adj`) %>%
      mutate(y.position = y_max * 1.05 + (row_number() - 1) * (y_max * 0.1),
             p.format = paste0("p=", signif(`p adj`, 3)))
    
    # Plot with significant pairs
    pmono <- data_sub %>% 
      ggplot(aes(x = treatment, y = value)) + 
      geom_jitter(aes(color = donorID)) + 
      geom_boxplot(alpha = 0, outlier.alpha = 0) +
      stat_pvalue_manual(tukey_sig, 
                         label = "p.format", 
                         y.position = "y.position",
                         xmin = "group1", xmax = "group2",
                         color = "red", size = 4) +
      theme_minimal() + 
      ggtitle(paste0(Cytokine, "\nANOVA p = ", signif(pval, 3))) +
      theme(legend.position = "none", text = element_text(size = 14))
    
  } else {
    # If no significant post-hoc pairs, plot without p-values
    pmono <- data_sub %>% 
      ggplot(aes(x = treatment, y = value)) + 
      geom_jitter(aes(color = donorID)) + 
      geom_boxplot(alpha = 0, outlier.alpha = 0) +
      theme_minimal() + 
      ggtitle(paste0(Cytokine, "\nANOVA p = ", signif(pval, 3))) +
      theme(legend.position = "none", text = element_text(size = 14))
  }
  
  print(pmono)
  
  path.mono <- paste0("./result/Dilute10_", Cytokine, "_ANOVA_Tukey.pdf")
  ggsave(filename = path.mono, plot = pmono, width = 5, height = 5, dpi = 300)
}



library(ggplot2)
library(dplyr)
library(ggpubr)
library(multcompView)

for (Cytokine in Cytokine.list){
  print(Cytokine)
  
  data_sub <- data.long.treat %>% filter(Dilution == "10", cytokine == Cytokine)
  
  # One-way ANOVA
  aov.test <- aov(value ~ treatment, data = data_sub)
  summary.aov <- summary(aov.test)
  pval <- summary.aov[[1]][["Pr(>F)"]][1]
  
  # Tukey post-hoc test
  tukey <- TukeyHSD(aov.test)
  
  # Generate Compact Letter Display
  tukey_letters <- multcompLetters4(aov.test, tukey)
  
  # Extract the letters for each group
  cld <- as.data.frame.list(tukey_letters$treatment)
  cld$treatment <- rownames(cld)
  rownames(cld) <- NULL
  colnames(cld)[1] <- "Letter"
  
  # Merge CLD back to data for plotting
  data_for_plot <- data_sub %>%
    left_join(cld, by = "treatment")
  
  # Calculate position for letter placement
  letter_positions <- data_for_plot %>%
    group_by(treatment, Letter) %>%
    summarise(y.position = max(value) * 1.05, .groups = 'drop')
  
  # Main plot
  pmono <- ggplot(data_sub, aes(x = treatment, y = value)) + 
    geom_jitter(aes(color = donorID)) + 
    geom_boxplot(alpha = 0, outlier.alpha = 0) +
    geom_text(data = letter_positions, 
              aes(x = treatment, y = y.position, label = Letter),
              size = 5, vjust = 0) +  # Adjust text position here
    theme_minimal() + 
    labs(title = Cytokine,
         subtitle = paste0("ANOVA p = ", signif(pval, 3))) + 
    theme(legend.position = "none", text = element_text(size = 14))
  
  print(pmono)
  
  path.mono <- paste0("./result/Dilute10_", Cytokine, "_CLD.pdf")
  ggsave(filename = path.mono, plot = pmono, width = 3, height = 5, dpi = 300)
}

#Each treatment group gets a letter (A, B, C…) based on Tukey’s results.

#Groups that share the same letter are not significantly different.






table(data.long.treat$cytokine)
Cytokine <- "CCL22"
Cytokine <- "CCL22"
p1 <- data.long.treat %>%filter(cytokine == Cytokine) %>% ggplot(aes(x=treatment, y= value)) + 
  geom_jitter(aes(color=donorID)) + geom_boxplot(alpha=0, outlier.alpha = 0)
p1 <- p1 + facet_wrap(.~Dilution, scales = "free_x") + theme_minimal() + ggtitle(Cytokine) + theme(legend.position = "none")
p1

#Shapiro-Wilk normality test
shapiro.test(filter(data.long.treat, cytokine == Cytokine)$value)

#loop results seperately for monocytes and spleen
Cytokine.list <- unique(data.long.treat$cytokine)

for (Cytokine in Cytokine.list){
  print(Cytokine)
  
  pmono <- data.long.treat %>% filter(Dilution == "10") %>% filter(cytokine == Cytokine) %>% 
    ggplot(aes(x=treatment, y= value)) + 
    geom_jitter(aes(color=donorID)) + geom_boxplot(alpha=0, outlier.alpha = 0)
  pmono <- pmono + stat_compare_means(method = "wilcox.test", ref.group = "None",  label = "p.signif", color = "red")+ 
    theme_minimal() + ggtitle(Cytokine) + theme(legend.position = "none")
  pmono
  
  pspleen <- data.long.treat %>% filter(cytokine == Cytokine) %>% 
    ggplot(aes(x=treatment, y= value)) + 
    geom_jitter(aes(color=donorID)) + geom_boxplot(alpha=0, outlier.alpha = 0)
  pspleen <- pspleen + facet_wrap(.~Dilution, scales = "free_x") + theme_minimal() + ggtitle(Cytokine) + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  pspleen
  
  path.mono <- paste0("./result/Dilute10_", Cytokine, ".pdf")
  ggsave(filename = path.mono, plot = pmono, width = 5, height = 5, dpi = 300)
  
  path.spleen <- paste0("./result/", Cytokine, ".pdf")
  ggsave(filename = path.spleen, plot = pspleen, width = 5, height = 5, dpi = 300)
  
}



########################################################################XC
#library(dplyr)
#library(ggplot2)
#library(ggpubr)

# Define list of treatments to keep
included_treatments <- c("None", "LPS", "Bb_D", "Bb_W")

# Get cytokine list
Cytokine.list <- unique(data.long.treat$cytokine)

for (Cytokine in Cytokine.list) {
  print(Cytokine)
  
  # ----------------- Monocyte Plot (Dilution 10 only) -----------------
  pmono <- data.long.treat %>%
    filter(Dilution == 10,
           treatment %in% included_treatments,
           cytokine == Cytokine) %>%
    ggplot(aes(x = treatment, y = value)) +
    geom_jitter(aes(color = donorID), width = 0.2, size = 1.5, alpha = 0.7) +
    geom_boxplot(alpha = 0, outlier.shape = NA) +
    stat_compare_means(method = "wilcox.test", ref.group = "None",
                       label = "p.signif", color = "red") +
    theme_minimal() +
    ggtitle(Cytokine) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 30, hjust = 1))
  
  # Save monocyte plot
  path.mono <- paste0("./result/XC_Dilute10_", Cytokine, ".pdf")
  ggsave(filename = path.mono, plot = pmono, width = 5, height = 5, dpi = 300)
}
########################################################################XC
####################################PCA
data.treat$treatment
cytokine_data_dplyr <- data.treat %>% dplyr::select(AGER:WNT7A) %>% dplyr::select(-c(CCL21, IL5RA))

cytokine_normalized <- scale(cytokine_data_dplyr)
head(cytokine_normalized)

spleen.pca <- PCA(cytokine_normalized,  graph = FALSE)
get_eig(spleen.pca)

ind.spleen <- get_pca_ind(spleen.pca)
ind.spleen

#cytokines contributed to PC1
pfv1 <- fviz_contrib(spleen.pca, choice = "var", axes = 1, top = 30)
pfv1 <- pfv1 + ggtitle("cytokines contributed to PC1 (LPS vs Bori)")
pfv1

#cytokines contributed to PC2
pfvz2 <- fviz_contrib(spleen.pca, choice = "var", axes = 2, top = 30)
pfvz2 <- pfvz2 + ggtitle("cytokines contributed to PC2 (none vs stimu)")
pfvz2

pfv1 / pfvz2

p3 <- fviz_pca_ind(spleen.pca,  col.ind = str_extract(data.treat$treatment, "^[^_]+"), label = "none",  repel = T)
p3 <- p3 + ggtitle("Difference between Stimuli")
p3

p3 <- fviz_pca_ind(spleen.pca,  col.ind = data.treat$treatment, label = "none",  repel = T)
p3 <- p3 + ggtitle("Difference between Stimuli")
p3


#####################
library(dplyr)

# Step 1: Filter Dilution = 10 and select cytokines
cytokine_data_dplyr <- data.treat %>%
  filter(Dilution == 10) %>%
  select(SampleID, AGER:WNT7A) %>%
  select(-c(CCL21, IL5RA))  # remove unwanted cytokines

# Step 2: Set SampleID as rownames
rownames(cytokine_data_dplyr) <- cytokine_data_dplyr$SampleID

# Step 3: Remove SampleID column (now it's rownames)
cytokine_data_dplyr <- select(cytokine_data_dplyr, -SampleID)

# Step 4: Scale and PCA
cytokine_normalized <- scale(cytokine_data_dplyr)
spleen.pca <- PCA(cytokine_normalized, graph = FALSE)

sample_names <- rownames(get_pca_ind(spleen.pca)$coord)

meta_matched <- data.treat %>%
  filter(SampleID %in% sample_names, Dilution == 10) %>%
  arrange(match(SampleID, sample_names))

library(factoextra)

p3 <- fviz_pca_ind(
  spleen.pca,
  col.ind = meta_matched$treatment,
  label = "none",
  repel = TRUE
) +
  ggtitle("PCA of Cytokine Response (Dilution = 10) by Treatment")

p3

##########

p3 <- fviz_pca_ind(
  spleen.pca,
  axes = c(1, 2),
  col.ind = meta_matched$treatment,  # coloring by treatment only
  label = "none",
  repel = TRUE,
  geom = "point"  # force points only, no shape mapping
) +
  ggtitle("PCA of Cytokine Response (Dilution = 10) by Treatment")

p3 + scale_shape_manual(values = rep(16, 12))  # use filled circles only

print(p3)

get_eig(spleen.pca)

#############PCA for none, Bb and LPS
library(dplyr)
selected_treatments <- c("None", "Bb_D", "Bb_W", "LPS")
filtered_data <- data.treat #%>%
  #filter(treatment %in% selected_treatments, Dilution == 10)

cytokine_data_dplyr <- filtered_data %>%
  dplyr::select(SampleID, AGER:WNT7A) %>%
  dplyr::select(-c(CCL21, IL5RA))

# Set SampleID as rownames
rownames(cytokine_data_dplyr) <- cytokine_data_dplyr$SampleID
cytokine_data_dplyr <- dplyr::select(cytokine_data_dplyr, -SampleID)

# Remove zero-variance columns (optional but recommended)
cytokine_data_dplyr <- cytokine_data_dplyr[, apply(cytokine_data_dplyr, 2, function(x) sd(x, na.rm = TRUE) > 0)]
library(FactoMineR)

cytokine_normalized <- scale(cytokine_data_dplyr)
spleen.pca <- PCA(cytokine_normalized, graph = FALSE)

sample_names <- rownames(get_pca_ind(spleen.pca)$coord)

meta_matched <- filtered_data %>%
  filter(SampleID %in% sample_names) %>%
  arrange(match(SampleID, sample_names))

library(factoextra)

p3 <- fviz_pca_ind(
  spleen.pca,
  col.ind = meta_matched$treatment,
  label = "none",
  repel = TRUE,
  addEllipses = TRUE,
  ellipse.level = 0.8  # default is 0.95
) +
  ggtitle("PCA: Cytokine Response (6 Treatments, Dilution = 10)")

print(p3)

ggsave(filename = "../../Manuscript/Figures/Figure 6-other peptidoglycans/PCA_Bb and LPS.pdf",p3, width =4, height = 3, dpi = 300)
#ggsave(filename = "../../Manuscript/Figures/Figure 6-other peptidoglycans/PCA_Bb and LPS2.pdf",p3, width =4, height = 3, dpi = 300)


# Extract variable contributions (loadings)
var <- get_pca_var(spleen.pca)

# View contributions to PC1 and PC2
contrib_df <- data.frame(
  Cytokine = rownames(var$contrib),
  PC1 = var$contrib[, "Dim.1"],
  PC2 = var$contrib[, "Dim.2"]
)

# Sort top contributing cytokines to PC1
top_PC1 <- contrib_df %>%
  arrange(desc(PC1)) %>%
  head(10)

# Sort top contributing cytokines to PC2
top_PC2 <- contrib_df %>%
  arrange(desc(PC2)) %>%
  head(10)

# Display top contributors
print(top_PC1)
print(top_PC2)




###########
#############PCA for digest
library(dplyr)
selected_treatments <- c("Bb_D","SO_D","VP_D","FN_D","SP_D", "LPS")


#("None", "LPS","Bb_D","Bb_W","FN_D","FN_W",
#  "SO_D", "SO_W", "SP_W","SP_D","VP_W","VP_D"
  
filtered_data <- data.treat %>%
  filter(treatment %in% selected_treatments, Dilution == 10)

cytokine_data_dplyr <- filtered_data %>%
  dplyr::select(SampleID, AGER:WNT7A) %>%
  dplyr::select(-c(CCL21, IL5RA))

# Set SampleID as rownames
rownames(cytokine_data_dplyr) <- cytokine_data_dplyr$SampleID
cytokine_data_dplyr <- dplyr::select(cytokine_data_dplyr, -SampleID)

# Remove zero-variance columns (optional but recommended)
cytokine_data_dplyr <- cytokine_data_dplyr[, apply(cytokine_data_dplyr, 2, function(x) sd(x, na.rm = TRUE) > 0)]
library(FactoMineR)

cytokine_normalized <- scale(cytokine_data_dplyr)
spleen.pca <- PCA(cytokine_normalized, graph = FALSE)

sample_names <- rownames(get_pca_ind(spleen.pca)$coord)

meta_matched <- filtered_data %>%
  filter(SampleID %in% sample_names) %>%
  arrange(match(SampleID, sample_names))

library(factoextra)

p3 <- fviz_pca_ind(
  spleen.pca,
  col.ind = meta_matched$treatment,
  label = "none",
  repel = TRUE,
  addEllipses = TRUE,
  ellipse.level = 0.5  # default is 0.95
) +
  ggtitle("PCA: Cytokine Response (6 Treatments, Dilution = 10)")

print(p3)

ggsave(filename = "../../Manuscript/Figures/Figure 6-other peptidoglycans/PCA_Digested form and LPS.pdf",p3, width = 4, height = 3, dpi = 300)

##########
selected_treatments <- c("Bb_W","SO_W","VP_W","FN_W","SP_W", "LPS")


#("None", "LPS","Bb_D","Bb_W","FN_D","FN_W",
#  "SO_D", "SO_W", "SP_W","SP_D","VP_W","VP_D"

filtered_data <- data.treat %>%
  filter(treatment %in% selected_treatments, Dilution == 10)

cytokine_data_dplyr <- filtered_data %>%
  select(SampleID, AGER:WNT7A) %>%
  select(-c(CCL21, IL5RA))

# Set SampleID as rownames
rownames(cytokine_data_dplyr) <- cytokine_data_dplyr$SampleID
cytokine_data_dplyr <- select(cytokine_data_dplyr, -SampleID)

# Remove zero-variance columns (optional but recommended)
cytokine_data_dplyr <- cytokine_data_dplyr[, apply(cytokine_data_dplyr, 2, function(x) sd(x, na.rm = TRUE) > 0)]
library(FactoMineR)

cytokine_normalized <- scale(cytokine_data_dplyr)
spleen.pca <- PCA(cytokine_normalized, graph = FALSE)

sample_names <- rownames(get_pca_ind(spleen.pca)$coord)

meta_matched <- filtered_data %>%
  filter(SampleID %in% sample_names) %>%
  arrange(match(SampleID, sample_names))

library(factoextra)

p3 <- fviz_pca_ind(
  spleen.pca,
  col.ind = meta_matched$treatment,
  label = "none",
  repel = TRUE,
  addEllipses = TRUE,
  ellipse.level = 0.3  # default is 0.95
) +
  ggtitle("PCA: Cytokine Response (6 Treatments, Dilution = 10)")

print(p3)

ggsave(filename = "../../Manuscript/Figures/Figure 6-other peptidoglycans/PCA_Intact form and LPS.pdf",p3, width = 4, height = 3, dpi = 300)

# Extract variable contributions (loadings)
var <- get_pca_var(spleen.pca)

# View contributions to PC1 and PC2
contrib_df <- data.frame(
  Cytokine = rownames(var$contrib),
  PC1 = var$contrib[, "Dim.1"],
  PC2 = var$contrib[, "Dim.2"]
)

# Sort top contributing cytokines to PC1
top_PC1 <- contrib_df %>%
  arrange(desc(PC1)) %>%
  head(10)

# Sort top contributing cytokines to PC2
top_PC2 <- contrib_df %>%
  arrange(desc(PC2)) %>%
  head(10)

# Display top contributors
print(top_PC1)
print(top_PC2)






#ggsave(filename = "./PCA.Stain.Differ.pdf", p3, width = 6, height = 6, dpi = 300)

data.treat$treatment_suffix <- ifelse(
  str_detect(data.treat$treatment, "_"),
  str_extract(data.treat$treatment, "(?<=_).*"),
  data.treat$treatment
)

p3.1 <- fviz_pca_ind(spleen.pca,  col.ind = data.treat$treatment_suffix, label = "none",  repel = T)
p3.1 <- p3.1 + ggtitle("Difference between Digest (D) and Whole (W) PG_bp")
p3.1
#ggsave(filename = "./PCA.WholeandDigest.pdf", p3.1, width = 6, height = 6, dpi = 300)

#Radar plot
library(fmsb)
#install.packages("fmsb")
library(purrr)
data.treat_10only <- filter(data.treat, Dilution == "10")
data.treat_10AGG <- aggregate(dplyr::select(data.treat_10only, AGER:WNT7A), by = list(data.treat_10only$treatment), FUN = mean )
dim(data.treat_10AGG)

radar_vars <- c("OSM", "CD83", "IL15", "VSTM1", "WNT16", 
                "CCL24", "CCL11", "FGF23", "TDFB3")

radar_vars <- c("IL1B", "CCL22", "IL6", "CCL3", "IL17A.IL17F", 
                "IL17B", "IL17C", "IL17F", "IL17RA")

radar_vars <- c("CSF3","IL5","VEGFA","CCL4","CXCL13","IL10",
               "CCL8","CCL2","CCL7","CCL11","IL17RA","CCL24",
               "CD40","CXCL16","CSF2","CCL17","CCL3","S100A12",
               "S100A9","IL18BP","IFNB1","CCL1","IL18","CXCL10",
               "CXCL5","CXCL3","CXCL1","IL6")

radar_data <- data.treat_10AGG[, c("Group.1", radar_vars)]

radar_scaled <- radar_data
radar_scaled[,-1] <- as.data.frame(scale(radar_data[,-1], center = FALSE, 
                                         scale = apply(radar_data[,-1], 2, max)))

par(mfrow = c(3, 4))  # Set up a 3x4 plot grid
for (i in 1:nrow(radar_scaled)) {
  row_data <- radar_scaled[i, ]
  
  # Add max and min rows required for fmsb radar plot
  radar_plot_data <- rbind(rep(1, length(radar_vars)),
                           rep(0, length(radar_vars)),
                           row_data[,-1])
  colnames(radar_plot_data) <- radar_vars
  
  radarchart(radar_plot_data,
             axistype = 1,
             pcol = "red",
             pfcol = rgb(1, 0, 0, 0.4),
             plwd = 2,
             cglcol = "grey",
             cglty = 1,
             axislabcol = "grey",
             caxislabels = seq(0, 1, 0.2),
             title = row_data$Group.1)
}

###########\
key_cytokines <- c("OSM", "CD83", "IL15", "VSTM1", "WNT16","CCL24", "CCL11", "FGF23", "TDFB3")
library(dplyr)
library(tidyr)
library(fmsb)

# Step 1: Select treatments and cytokines
selected_treatments <- c("Bb_W","SO_W","VP_W","FN_W","SP_W")
key_cytokines <- c("OSM", "CD83", "IL15", "VSTM1", "WNT16","CCL24", "CCL11", "FGF23")

# Step 2: Filter for selected treatments
filtered_data <- data.treat_10only %>%
  filter(treatment %in% selected_treatments)

# Step 3: Calculate mean expression for selected cytokines
cytokine_means <- filtered_data %>%
  select(treatment, all_of(key_cytokines)) %>%
  group_by(treatment) %>%
  summarize(across(everything(), ~mean(.x, na.rm = TRUE)))

# Step 4: Reshape and prepare for radar
radar_matrix <- as.data.frame(cytokine_means)
rownames(radar_matrix) <- radar_matrix$treatment
radar_matrix <- radar_matrix[, -1]  # remove treatment column
radar_data <- as.data.frame(t(radar_matrix))  # transpose

# Step 5: Add max/min rows
max_row <- apply(radar_data, 2, max, na.rm = TRUE)
min_row <- apply(radar_data, 2, min, na.rm = TRUE)
radar_ready <- rbind(max_row, min_row, radar_data)

# Step 6: Plot radar chart
radarchart(radar_ready,
           axistype = 1,
           pcol = c("blue", "red", "green", "purple", "orange"),
           plwd = 2,
           plty = 1,
           cglcol = "grey", cglty = 1, caxislabels = NULL,
           cglwd = 0.8, title = "Key Cytokine Expression Across Treatments")

legend("topright", legend = colnames(radar_ready), 
       col = c("blue", "red", "green", "purple"), 
       lty = 1, lwd = 2)

library(fmsb)

# Example: Create the radar plot with custom axis labels
radarchart(radar_ready,
           axistype = 1,
           pcol = c("blue", "red", "green", "purple", "orange"),
           plwd = 2,
           plty = 1,
           cglcol = "grey", cglty = 1,
           axislabcol = "black",
           caxislabels = seq(5, 25, 5),  # adjust based on your min/max values
           cglwd = 0.8,
           title = "Cytokine Expression Across PG Treatments")

legend("topright",
       legend = colnames(radar_ready),
       col = c("blue", "red", "green", "purple", "orange"),
       lty = 1, lwd = 2)

library(fmsb)

data <- as.data.frame(rbind(
  c(10, 10, 10, 10, 10),   # max
  c(0, 0, 0, 0, 0),       # min
  c(7, 5, 6, 3, 4),       # sample 1
  c(6, 7, 3, 5, 6)        # sample 2
))
colnames(data) <- c("A", "B", "C", "D", "E")

radarchart(data,
           pcol = c("red", "blue"),
           plwd = 2)

# Convert your first column to row names
radar_fixed <- radar_ready
rownames(radar_fixed) <- radar_fixed[[1]]
radar_fixed <- radar_fixed[, -1]  # remove the first column (cytokine names)

# Make sure all values are numeric
radar_fixed[] <- lapply(radar_fixed, as.numeric)

