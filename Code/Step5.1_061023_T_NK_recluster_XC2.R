#Last Updated 04/09/2024

library(Seurat)
#library(MAST)
library(monocle3)
library(dplyr)
library(stringr)
library(Scillus)
library(ggplot2)
library(ComplexHeatmap)
library(ggrepel)
library(ggpubr)
library(DESeq2)
library(EnhancedVolcano)
library(cowplot)
library(reshape2)
#install.packages("devtools")
library(devtools)
#install_github("RGLab/MAST")


#mac
setwd("/Users/xinchen/Library/CloudStorage/Dropbox/lyme_disease/R7_NR/")

source("./Code/xztools.R")
source("./Code/00_colorKey.R")
patient_data <- read.csv("./Results/patient_data_addPTLDN.csv", header = T)


load("./Data/Step5_070523_T.NK.clean_Annotated.RData")
T.NK.clean_Annotated$celltype[T.NK.clean_Annotated$integrated_snn_res.0.7 == "8"] <- "GammaDelta_T"
T.NK.clean_Annotated$celltype[T.NK.clean_Annotated$integrated_snn_res.0.7 == "10"] <- "TFH_CD4"
T.NK.clean_Annotated$celltype[T.NK.clean_Annotated$integrated_snn_res.0.7 == "13"] <- "Proliferating_Cell"

table(T.NK.clean_Annotated$celltype)
table(T.NK.clean_Annotated$subject)

library(dplyr)
library(stringr)
library(stringr)

T.NK.clean_Annotated$group <- case_when(
  str_detect(T.NK.clean_Annotated$subject, "^HD") ~ "HD",
  str_detect(T.NK.clean_Annotated$subject, "^RTH") ~ "RTH",
  str_detect(T.NK.clean_Annotated$subject, "^PTLDN") ~ "PTLDN",
  str_detect(T.NK.clean_Annotated$subject, "^PTLD") ~ "PTLD",
  TRUE ~ "Other"
)
table(T.NK.clean_Annotated$subject)

library(stringr)

# Subset PTLD and RTH samples
#colnames(T.NK.clean_Annotated@meta.data)
library(stringr)

# Define subjects to exclude
remove_subjects <- c("RTH2V5b1", "RTH5V3b5", "PTLD2V1")

# Select cell barcodes that match PTLD or RTH and are not in the removal list
ptld_rth_cells <- colnames(T.NK.clean_Annotated)[
  str_detect(T.NK.clean_Annotated$subject, "^PTLD|^RTH") &
    !(T.NK.clean_Annotated$subject %in% remove_subjects)
]

# Subset Seurat object
ptld_rth <- subset(T.NK.clean_Annotated, cells = ptld_rth_cells)

ptld_rth$group <- case_when(
  str_detect(ptld_rth$subject, "^PTLD") ~ "PTLD",
  str_detect(ptld_rth$subject, "^RTH") ~ "RTH"
)

ptld_rth$visit <- case_when(
  str_detect(ptld_rth$subject, "V2") ~ "v2",
  str_detect(ptld_rth$subject, "V3") ~ "v3",
  str_detect(ptld_rth$subject, "V5") ~ "v5",
  str_detect(ptld_rth$subject, "V7") ~ "v7",
  TRUE ~ NA_character_
)

# Optional: check distribution
table(ptld_rth$visit)


library(dplyr)

meta <- ptld_rth@meta.data

# Count cells per group × visit × celltype
cell_counts <- meta %>%
  filter(!is.na(visit)) %>%
  group_by(group, visit, celltype) %>%
  summarise(count = n(), .groups = "drop")

# Total cells per group × visit
totals <- meta %>%
  filter(!is.na(visit)) %>%
  group_by(group, visit) %>%
  summarise(total = n(), .groups = "drop")

# Join and calculate percentages
cell_percentages <- left_join(cell_counts, totals, by = c("group", "visit")) %>%
  mutate(percent = 100 * count / total)

library(ggplot2)

# Ensure correct visit order
cell_percentages$visit <- factor(cell_percentages$visit, levels = c("v2", "v3", "v5", "v7"))

# Plot
ggplot(cell_percentages, aes(x = visit, y = percent, color = group, group = group)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  facet_wrap(~ celltype, scales = "free_y") +
  ylab("Cell Percentage (%)") +
  xlab("Visit") +
  theme_minimal() +
  theme(legend.title = element_blank())


library(dplyr)

meta <- ptld_rth@meta.data

subject_percentages <- meta %>%
  filter(!is.na(visit)) %>%
  group_by(subject, group, visit, celltype) %>%
  summarise(count = n(), .groups = "drop") %>%
  left_join(
    meta %>%
      group_by(subject, visit) %>%
      summarise(total = n(), .groups = "drop"),
    by = c("subject", "visit")
  ) %>%
  mutate(percent = 100 * count / total)



library(ggpubr)
library(dplyr)
library(ggplot2)

plots <- list()

# Define color manually: RTH = blue, PTLD = default (e.g., red)
manual_colors <- c("RTH" = "#0072B5FF", "PTLD" = "#BC3C29FF")  # Blue and reddish-orange

# Define selected cell types
selected_celltypes <- c("Activated_CD4", "GZMB_CD8", "GZMK_CD8","CD56bright_NK", "CD56dim_NK", 
                         "MAIT_T", "TFH_CD4")
#selected_celltypes <- c("MAIT_T", "TFH_CD4")

for (ct in selected_celltypes) {
  df_plot <- subject_percentages %>% filter(celltype == ct)
  
  p <- ggboxplot(df_plot, x = "visit", y = "percent", color = "group", add = "jitter") +
    stat_compare_means(aes(group = group), method = "t.test", label = "p.format") +  # show numeric p
    scale_color_manual(values = manual_colors) +
    ggtitle(ct) +
    ylab("Percentage") +
    theme_minimal() +
    theme(legend.title = element_blank())
  
  plots[[ct]] <- p
}

# Combine all plots into one object
library(patchwork)
p <- wrap_plots(plots, ncol = 5)
p
# Save combined plot
ggsave("./Results/061925/CD4_CD8_PTLD_vs_RTH_pvalue.pdf", plot = p, width = 10, height = 5, dpi=300)



##########

manual_colors <- c("RTH" = "#0072B5FF", "PTLD" = "#BC3C29FF")
selected_celltypes <- c("MAIT_T", "TFH_CD4")
plots <- list()

for (ct in selected_celltypes) {
  df_plot <- subject_percentages %>% filter(celltype == ct)
  df_plot$group <- factor(df_plot$group, levels = c("RTH", "PTLD"))  # enforce order
  
  p <- ggboxplot(df_plot, x = "visit", y = "percent", color = "group", add = "jitter", 
                 palette = manual_colors) +
    stat_compare_means(aes(group = group), method = "t.test", label = "p.format") +
    ggtitle(ct) +
    ylab("Percentage") +
    theme_minimal() +
    theme(legend.title = element_blank())
  
  plots[[ct]] <- p
}

library(patchwork)
p <- wrap_plots(plots)
ggsave("./Results/061925/Tfh_MAIT_PTLD_vs_RTH_pvalue.pdf", plot = p, width = 10, height = 3.5, dpi = 300)


############

manual_colors <- c("RTH" = "#0072B5FF", "PTLD" = "#BC3C29FF")
selected_celltypes <- c("MAIT_T", "TFH_CD4")

for (ct in selected_celltypes) {
  df_plot <- subject_percentages %>% filter(celltype == ct)
  df_plot$group <- factor(df_plot$group, levels = c("RTH", "PTLD"))
  
  p <- ggboxplot(df_plot, x = "visit", y = "percent", color = "group", add = "jitter", 
                 palette = manual_colors) +
    stat_compare_means(aes(group = group), method = "t.test", label = "p.format") +
    ggtitle(ct) +
    ylab("Percentage") +
    theme_minimal() +
    theme(legend.title = element_blank())
  
  # Save each plot separately
  ggsave(
    filename = paste0("./Results/061925/", ct, "_PTLD_vs_RTH_pvalue.pdf"),
    plot = p,
    width = 5, height = 4, dpi = 300
  )
}



############# I want to check the same for PTLD and HD##########not significance#############


library(stringr)

# Subset cell barcodes for HD and PTLDN only
hd_ptldn_cells <- colnames(T.NK.clean_Annotated)[
  str_detect(T.NK.clean_Annotated$subject, "^HD|^PTLDN")
]

# Subset Seurat object
hd_ptldn <- subset(T.NK.clean_Annotated, cells = hd_ptldn_cells)

# Assign group labels
hd_ptldn$group <- case_when(
  str_detect(hd_ptldn$subject, "^HD") ~ "HD",
  str_detect(hd_ptldn$subject, "^PTLDN") ~ "PTLDN"
)

meta <- hd_ptldn@meta.data

subject_percentages <- meta %>%
  group_by(subject, group, celltype) %>%
  summarise(count = n(), .groups = "drop") %>%
  left_join(
    meta %>%
      group_by(subject) %>%
      summarise(total = n(), .groups = "drop"),
    by = "subject"
  ) %>%
  mutate(percent = 100 * count / total)

library(ggpubr)
library(patchwork)

manual_colors <- c("HD" = "#009E73", "PTLDN" = "#D55E00")  # Green and red-orange
plots <- list()

# Optional: restrict to specific cell subsets
selected_celltypes <- c("Activated_CD4", "CD56bright_NK", "CD56dim_NK", 
                        "GZMB_CD8", "GZMK_CD8", "MAIT_T", "TFH_CD4")

for (ct in selected_celltypes) {
  df_plot <- subject_percentages %>% filter(celltype == ct)
  
  p <- ggboxplot(df_plot, x = "group", y = "percent", color = "group", add = "jitter") +
    stat_compare_means(method = "t.test", label = "p.format") +
    scale_color_manual(values = manual_colors) +
    ggtitle(ct) +
    ylab("Percentage") +
    theme_minimal() +
    theme(legend.position = "none")
  
  plots[[ct]] <- p
}

# Combine and save
p <- wrap_plots(plots, ncol = 3)
ggsave("./Results/HD_vs_PTLDN_Selected_Celltypes.pdf", plot = p, width = 12, height = 10, dpi = 300 )



library(ggpubr)
library(patchwork)

manual_colors <- c("HD" = "#009E73", "PTLDN" = "#D55E00")  # Green and red-orange
plots <- list()

# Loop over all available cell types
for (ct in unique(subject_percentages$celltype)) {
  df_plot <- subject_percentages %>% filter(celltype == ct)
  
  p <- ggboxplot(df_plot, x = "group", y = "percent", color = "group", add = "jitter") +
    stat_compare_means(method = "t.test", label = "p.format") +
    scale_color_manual(values = manual_colors) +
    ggtitle(ct) +
    ylab("Percentage") +
    theme_minimal() +
    theme(legend.position = "none")
  
  plots[[ct]] <- p
}

# Combine and save
p <- wrap_plots(plots, ncol = 3)
ggsave("./Results/061925/HD_vs_PTLDN_All_Celltypes.pdf", plot = p, width = 14, height = 20)

###############To identify differentially expressed genes (DEGs) between PTLD and RTH at any visit (v2, v3, v5, or v7)#######

ptld_rth_cells <- colnames(T.NK.clean_Annotated)[
  str_detect(T.NK.clean_Annotated$subject, "^PTLD|^RTH") &
    !(T.NK.clean_Annotated$subject %in% c("RTH2V5b1", "RTH5V3b5", "PTLD2V1"))
]
ptld_rth <- subset(T.NK.clean_Annotated, cells = ptld_rth_cells)

ptld_rth$group <- case_when(
  str_detect(ptld_rth$subject, "^PTLD") ~ "PTLD",
  str_detect(ptld_rth$subject, "^RTH") ~ "RTH"
)

ptld_rth$visit <- case_when(
  str_detect(ptld_rth$subject, "V2") ~ "v2",
  str_detect(ptld_rth$subject, "V3") ~ "v3",
  str_detect(ptld_rth$subject, "V5") ~ "v5",
  str_detect(ptld_rth$subject, "V7") ~ "v7",
  TRUE ~ NA_character_
)


library(Seurat)
deg_list <- list()

for (v in c("v2", "v3", "v5", "v7")) {
  seurat_sub <- subset(ptld_rth, subset = visit == v)
  
  # Make sure group has both PTLD and RTH
  if (length(unique(seurat_sub$group)) == 2) {
    de_result <- FindMarkers(seurat_sub, 
                             ident.1 = "PTLD", ident.2 = "RTH", 
                             group.by = "group", 
                             logfc.threshold = 0.25, 
                             min.pct = 0.1, 
                             test.use = "wilcox")
    de_result$gene <- rownames(de_result)
    de_result$visit <- v
    deg_list[[v]] <- de_result
  }
}

# Combine results
deg_combined <- do.call(rbind, deg_list)


deg_signif <- deg_combined %>%
  filter(p_val_adj < 0.05)

# List of unique DE genes across visits
unique_deg_genes <- unique(deg_signif$gene)
length(unique_deg_genes)

# Top 10 DEGs by lowest adj p-value
head(deg_signif[order(deg_signif$p_val_adj), ], 10)

# Save to CSV
write.csv(deg_signif, "./Results/PTLD_vs_RTH_DEGs_by_visit.csv", row.names = FALSE)


install.packages("VennDiagram")  # if not already installed
library(VennDiagram)







