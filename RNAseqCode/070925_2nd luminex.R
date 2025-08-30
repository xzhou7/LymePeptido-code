##################2nd luminex data#########
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(stringr)
library(FactoMineR)
library(corrr)
library(factoextra)
library(patchwork)
library(vegan)
library(ggpubr)

# Set working directory
setwd("~/Library/CloudStorage/Dropbox/lyme_disease/peptiglygen_stimulation/Mono_LPS_Pepti_Sti/2nd luminex experiment/")
getwd()


# Read and annotate raw data
raw.data.df <- read.csv(file = "./2nd_raw_count.csv", header = TRUE, row.names = 1)

# Subset rows where Name contains "none", "Di_pepti", or "LPS"
raw.data.df <- raw.data.df[grepl("none|Di_pepti|LPS", raw.data.df$Name), ]


# Step 2: Add treatment column based on Name content
raw.data.df <- raw.data.df %>%
  mutate(
    treatment = case_when(
      str_detect(Name, "none") ~ "none",
      str_detect(Name, "Di_pepti") ~ "peptidoglycan",
      str_detect(Name, "LPS") ~ "LPS",
      TRUE ~ NA_character_
    ),
    treatment = factor(treatment, levels = c("none", "peptidoglycan", "LPS"))
  )

# Optional: Check result
table(raw.data.df$treatment)
head(raw.data.df[, c("Name", "treatment")])


# Add donor column by extracting the number after the second dash
raw.data.df <- raw.data.df %>%
  mutate(
    donor = str_extract(Name, "(?<=-)[0-9]{4}(?=_)")
  )

# Optional: Convert to numeric if needed
raw.data.df$donor <- as.numeric(raw.data.df$donor)

# Check results
head(raw.data.df[, c("Name", "donor", "treatment")])

# Remove rows with NA in donor column
raw.data.df <- raw.data.df %>%
  filter(!is.na(donor))

##################bargraph cytokines##########

library(ggplot2)
library(ggpubr)
library(patchwork)

# Cytokines of interest
selected_cytokines <- c("TNF.a", "IL.6", "IL.1b", "MIP.1a.CCL3", "MIP.1b.CCL4", "MDC.CCL22")

# Set treatment color palette (case matches actual treatment levels)
treatment_colors <- c("none" = "#3d51a3",      # Blue
                      "LPS" = "#ee2926",       # Red
                      "peptidoglycan" = "#3eb54a")  # Green

# Create and store plots in a list
plot_list <- list()
for (cytokine.i in selected_cytokines) {
  ymax <- max(raw.data.df[[cytokine.i]], na.rm = TRUE)
  
  p <- ggplot(raw.data.df, aes_string(x = "treatment", y = cytokine.i, fill = "treatment")) +
    geom_boxplot(outlier.alpha = 0, alpha = 0.4) +
    geom_jitter(color = "black", size = 0.6, width = 0.2) +
    stat_compare_means(comparisons = list(c("none", "peptidoglycan")), 
                       method = "wilcox.test", label.y = ymax * 1.1, color = "black") +
    stat_compare_means(comparisons = list(c("LPS", "none")), 
                       method = "wilcox.test", label.y = ymax * 1.0, color = "black") +
    stat_compare_means(comparisons = list(c("LPS", "peptidoglycan")), 
                       method = "wilcox.test", label.y = ymax * 0.9, color = "black") +
    scale_fill_manual(values = treatment_colors) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 9),
      axis.title.y = element_text(size = 10)
    ) +
    labs(title = cytokine.i, y = "log2 MFI")
  
  plot_list[[cytokine.i]] <- p
}

# Combine all plots into a 2x3 layout
combined_plot <- wrap_plots(plot_list, nrow = 2, ncol = 3)

# Show the plot
print(combined_plot)

# Save to PDF
ggsave("cytokine_panel.pdf", combined_plot, width = 10, height = 6, dpi = 300)









###########################



# Add treatment from row names
raw.data.df$treatment <- tolower(sub(".*_", "", row.names(raw.data.df)))

raw.data.df$treatment <- sub(".*_", "", row.names(raw.data.df))
raw.data.df$treatment <- factor(raw.data.df$treatment, levels = c("none", "L18MDP", "Peptidoglycan", "LPS", "R848"))
raw.data.df$subjectid <- sub("_.*", "", row.names(raw.data.df))

# Subset data for desired comparisons
data.sub <- raw.data.df %>%
  filter(treatment %in% c("none", "LPS", "Peptidoglycan"))

#########Identify Cytokines Upregulated by Both LPS and Peptidoglycan for bargraph#################
library(tidyr)


# Subset to relevant treatments
data.sub <- raw.data.df %>%
  filter(treatment %in% c("none", "LPS", "Peptidoglycan"))

# Select numeric cytokine columns
cytokine_cols <- raw.data.df %>%
  select(where(is.numeric)) %>%
  colnames()


# Compute group means
cytokine_means <- data.sub %>%
  group_by(treatment) %>%
  summarise(across(all_of(cytokine_cols), ~mean(.x, na.rm = TRUE))) %>%
  pivot_longer(-treatment, names_to = "cytokine", values_to = "mean_expr") %>%
  pivot_wider(names_from = treatment, values_from = mean_expr)




# Calculate log2 fold-change
cytokine_means <- cytokine_means %>%
  mutate(log2FC_LPS = LPS - none,
         log2FC_Pep = Peptidoglycan - none,
         log2FC_max = pmax(log2FC_LPS, log2FC_Pep, na.rm = TRUE))

top20_upregulated <- cytokine_means %>%
  arrange(desc(log2FC_max)) %>%
  slice_head(n = 20) %>%
  pull(cytokine)


selected_cytokines <- top20_upregulated




# Define cytokines of interest
selected_cytokines <- c("TNF.a", "IL.6", "IL.1b", "MIP.1a.CCL3", "MIP.1b.CCL4", "MDC.CCL22")

# Define consistent treatment colors
treatment_colors <- c("none" = "#3d51a3",      # Blue
                      "LPS" = "#ee2926",       # Red
                      "Peptidoglycan" = "#3eb54a")  # Green

# Create and store plots in a list
plot_list <- list()
for (cytokine.i in selected_cytokines) {
  ymax <- max(data.sub[[cytokine.i]], na.rm = TRUE)
  
  p <- ggplot(data.sub, aes_string(x = "treatment", y = cytokine.i, fill = "treatment")) +
    geom_boxplot(outlier.alpha = 0, alpha = 0.4) +
    geom_jitter(color = "black", size = 0.6, width = 0.2) +
    stat_compare_means(comparisons = list(c("none", "Peptidoglycan")), 
                       method = "wilcox.test", label.y = ymax * 1.1, color = "black") +
    stat_compare_means(comparisons = list(c("LPS", "none")), 
                       method = "wilcox.test", label.y = ymax * 1.0, color = "black") +
    stat_compare_means(comparisons = list(c("LPS", "Peptidoglycan")), 
                       method = "wilcox.test", label.y = ymax * 0.9, color = "black") +
    scale_fill_manual(values = treatment_colors) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 9),
      axis.title.y = element_text(size = 10)
    ) +
    labs(title = cytokine.i, y = "log2 MFI")
  
  plot_list[[cytokine.i]] <- p
}

# Arrange all plots into a 2x3 grid
combined_plot <- wrap_plots(plot_list, nrow = 2, ncol = 3)

# Display in viewer
combined_plot

# Save to PDF
ggsave("cytokine_panel.pdf", combined_plot, width = 10, height = 6, dpi = 300)

