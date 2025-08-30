# Luminex Data Analysis
# Third Experiment (LPS/Pepti mono stimulation)
# Author: Xin Chen, Ph.D. 
# Date: 2024-12-30
# Last Updated: 2025-06-22
#Cleanup for 2nd luminex


library(ggplot2)
library(dplyr)
library(stringr)
library(FactoMineR)
library(corrr)
library(factoextra)
library(patchwork)
library(vegan)
library(ggpubr)

setwd("~/Library/CloudStorage/Dropbox/lyme_disease/peptiglygen_stimulation/Mono_LPS_Pepti_Sti/2nd luminex experiment/")
getwd()

# Load data
raw.data.df <- read.csv(file = "./2nd_raw_count.csv", header = TRUE, row.names = 1)

####
# Make rownames unique by appending suffixes
raw.data.df$Name_unique <- make.unique(as.character(raw.data.df$Name))
rownames(raw.data.df) <- raw.data.df$Name_unique
meta.cols <- c("Well", "Name", "Type", "Name_unique")
cytokine.df <- raw.data.df[, !(colnames(raw.data.df) %in% meta.cols)]
cytokine.df <- log2(cytokine.df + 1)

# Load required libraries
library(dplyr)
library(tibble)

# Step 1: Group and average duplicated samples by Name, remove control samples
cytokine.df <- raw.data.df %>%
  group_by(Name) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  filter(!Name %in% c("S1", "S2", "S3", "S4", "S5", "S6", "S7", "SB"))

# Step 2: Add treatment Group column based on Name
cytokine.df$Group <- case_when(
  grepl("none",     cytokine.df$Name, ignore.case = TRUE) ~ "none",
  grepl("Di_pepti", cytokine.df$Name, ignore.case = TRUE) ~ "Di_pepti",
  grepl("LPS",      cytokine.df$Name, ignore.case = TRUE) ~ "LPS",
  grepl("W_pepti",  cytokine.df$Name, ignore.case = TRUE) ~ "W_pepti",
  grepl("L18MDP",   cytokine.df$Name, ignore.case = TRUE) ~ "L18MDP",
  TRUE ~ "unknown"
)

# Step 3: Extract Donor ID from sample Name (e.g., XC-A1-1362_LPS → 1362)
cytokine.df$Donor <- sub(".*?-([0-9]+)_.*", "\\1", cytokine.df$Name)

# Step 4: Add Source column (PBMC vs Spleen)
cytokine.df$Source <- case_when(
  grepl("Sleen|Spleen", cytokine.df$Name, ignore.case = TRUE) ~ "Spleen",
  TRUE ~ "PBMC"
)

# Step 5: Set Name as rownames
cytokine.df <- cytokine.df %>% column_to_rownames("Name")

# Step 6: Log2-transform cytokine values (exclude metadata columns)
meta_cols <- c("Group", "Donor", "Source")
cytokine_values <- cytokine.df[, !colnames(cytokine.df) %in% meta_cols]
cytokine_log <- log2(cytokine_values + 1)

# Step 7: Combine log-transformed values with metadata
cytokine.df <- cbind(cytokine_log, cytokine.df[, meta_cols])


head(cytokine.df[, 1:5])

# List all sample names labeled as "Spleen"
#spleen_labeled <- rownames(cytokine.df)[cytokine.df$Source == "Spleen"]
# Print them
#print(spleen_labeled)

##############PCA#####################NO helpful

cytokine_only <- cytokine.df[, !(colnames(cytokine.df) %in% c("Group", "Donor", "Source"))]

pca_res <- prcomp(cytokine_only, scale. = TRUE)
pca_df <- as.data.frame(pca_res$x)
pca_df$Group <- cytokine.df$Group
pca_df$Source <- cytokine.df$Source

ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, shape = Source)) +
  geom_point(size = 3) +
  labs(title = "PCA of Cytokine Profiles") +
  theme_minimal()


############SUBSET OUT none, LPS and two pepti###########

# Subset cytokine.df to keep only specified groups
cytokine.subset <- cytokine.df %>%
  filter(Group %in% c("none", "W_pepti", "Di_pepti", "LPS"))

#PCA

# Step 1: Extract cytokine matrix (exclude metadata)
cytokine_matrix <- cytokine.subset[, !(colnames(cytokine.subset) %in% c("Group", "Donor", "Source"))]

# Step 2: Run PCA
pca_res <- prcomp(cytokine_matrix, scale. = TRUE)

# Step 3: Make PCA results into a data frame for plotting
pca_df <- as.data.frame(pca_res$x)
pca_df$Group <- cytokine.subset$Group
pca_df$Source <- cytokine.subset$Source
pca_df$Donor <- cytokine.subset$Donor

# Step 4: Plot
library(ggplot2)

ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, shape = Source)) +
  geom_point(size = 3) +
  labs(
    title = "PCA of Cytokine Profiles",
    x = paste0("PC1 (", round(100 * summary(pca_res)$importance[2, 1], 1), "% variance)"),
    y = paste0("PC2 (", round(100 * summary(pca_res)$importance[2, 2], 1), "% variance)")
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "right")











