setwd("~/Library/CloudStorage/Dropbox/lyme_disease/R7_NR/")
#load necessary package
library(data.table)
library(dplyr)
library(Seurat)
library(ggplot2)
library(stringr)
library(patchwork)
library(metap)
library(multtest)
library(Rcpp)
library(cowplot)

load("./Data/Step4_022723_Clean_cluster.RData")

#identical(pbmc.clean$seurat_clusters, pbmc.clean$integrated_snn_res.1.2)
Idents(pbmc.clean) <- pbmc.clean$integrated_snn_res.1.2
pbmc.clean <- RenameIdents(object = pbmc.clean, '0' = "X00_NavieCD4T")
pbmc.clean <- RenameIdents(object = pbmc.clean, '1' = "X01_NaiveCD4T")
pbmc.clean <- RenameIdents(object = pbmc.clean, '2' = "X02_CD14mono")
pbmc.clean <- RenameIdents(object = pbmc.clean, '3' = "X03_Treg")
pbmc.clean <- RenameIdents(object = pbmc.clean, '4' = "X04_NK")
pbmc.clean <- RenameIdents(object = pbmc.clean, '5' = "X05_B")
pbmc.clean <- RenameIdents(object = pbmc.clean, '6' = "X06_CD14mono")
pbmc.clean <- RenameIdents(object = pbmc.clean, '7' = "X07_GZMB_CD8T")
pbmc.clean <- RenameIdents(object = pbmc.clean, '8' = "X08_CD4T")
pbmc.clean <- RenameIdents(object = pbmc.clean, '9' = "X09_CD14mono")
pbmc.clean <- RenameIdents(object = pbmc.clean, '10' = "X10_TH17_CD4T")
pbmc.clean <- RenameIdents(object = pbmc.clean, '11' = "X11_GZMK_CD8T")
pbmc.clean <- RenameIdents(object = pbmc.clean, '12' = "X12_NaiveCD8T")
pbmc.clean <- RenameIdents(object = pbmc.clean, '13' = "X13_IgD_B")
pbmc.clean <- RenameIdents(object = pbmc.clean, '14' = "X14_CD14mono")
pbmc.clean <- RenameIdents(object = pbmc.clean, '15' = "X15_CD4T")
pbmc.clean <- RenameIdents(object = pbmc.clean, '16' = "X16_CD14mono")
pbmc.clean <- RenameIdents(object = pbmc.clean, '17' = "X17_NK")
pbmc.clean <- RenameIdents(object = pbmc.clean, '18' = "X18_CD16mono")
pbmc.clean <- RenameIdents(object = pbmc.clean, '19' = "X19_myeloid")
pbmc.clean <- RenameIdents(object = pbmc.clean, '20' = "X20_CD4T")
pbmc.clean <- RenameIdents(object = pbmc.clean, '21' = "X21_DC")
pbmc.clean <- RenameIdents(object = pbmc.clean, '22' = "X22_IntermediateMono")
pbmc.clean <- RenameIdents(object = pbmc.clean, '23' = "X23_T")
pbmc.clean <- RenameIdents(object = pbmc.clean, '24' = "X24_CD14mono")
pbmc.clean <- RenameIdents(object = pbmc.clean, '25' = "X25_B")
pbmc.clean <- RenameIdents(object = pbmc.clean, '26' = "X26_NKT")
pbmc.clean <- RenameIdents(object = pbmc.clean, '27' = "X27_DR+GZMB+DC")
pbmc.clean <- RenameIdents(object = pbmc.clean, '28' = "X28_CD14mono")
pbmc.clean <- RenameIdents(object = pbmc.clean, '29' = "X29_NKT")


pbmc.clean$CellType.1 <- Idents(pbmc.clean)
pbmc.clean$CellType.1 <- factor(pbmc.clean$CellType.1, levels = sort(unique(pbmc.clean$CellType.1),decreasing=T))
Idents(pbmc.clean) <- pbmc.clean$CellType.1
p.label <- DimPlot(pbmc.clean,label = T,raster=FALSE) + coord_equal()
p.label

subjectID <- c("PTLD1V2","PTLD1V3","PTLD1V5","PTLD2V2","PTLD2V3","PTLD2V5","PTLD2V7","PTLD3V2",
               "PTLD3V3","PTLD3V5","PTLD3V7","PTLD4V3","RTH1V2","RTH1V3",  
               "RTH1V5","RTH1V7b1","RTH2V2","RTH2V5b1","RTH2V7","RTH3V3","RTH3V5",  
               "RTH4V2","RTH4V3","RTH4V5","RTH4V7","RTH5V2","RTH5V3b5","RTH5V7","RTH6V3",  
               "RTH6V5","RTH6V7")

PTLD_RTH <- subset(pbmc.clean, subject %in% subjectID)
table(PTLD_RTH$subject)
#DimPlot(PTLD_RTH)
library(stringr)

# Create a new object (or modify in place)
PTLD_RTH$donor <- str_extract(PTLD_RTH$subject, "^[A-Z]+\\d+")
PTLD_RTH$visit <- str_extract(PTLD_RTH$subject, "V\\d+[a-z0-9]*")


PTLD_RTH$visit_group <- str_replace(PTLD_RTH$visit, "b[0-9]+", "")  # V3b2 → V3
PTLD_RTH$sample_id <- paste(PTLD_RTH$donor, PTLD_RTH$visit_group, sep = "_")
table(PTLD_RTH$donor)
table(PTLD_RTH$visit_group)



#######Longitudinal Analysis Methods
#1. Pseudobulk + Linear Mixed Models
#Recommended: aggregates at subject-timepoint level and models repeated measures
#a. Create Pseudobulk Matrix by Subject + Timepoint
PTLD_RTH$sample_id <- paste(PTLD_RTH$donor, PTLD_RTH$visit_group, sep = "_")
pb_matrix <- AggregateExpression(PTLD_RTH, group.by = "sample_id", assays = "RNA", slot = "counts")$RNA

#Step 2: Convert to logCPM (counts per million)

#PCA
library(edgeR)

# Create DGEList
dge <- DGEList(counts = pb_matrix)

# Filter lowly expressed genes (optional but recommended)
keep_genes <- filterByExpr(dge)
dge <- dge[keep_genes, , keep.lib.sizes = FALSE]

# Normalize and calculate logCPM
dge <- calcNormFactors(dge)
logCPM <- cpm(dge, log = TRUE)
#write.csv(logCPM, "./Results/061725/061725_logCPM.csv", row.names = FALSE)

#Step 3: Create Sample Metadata Table
meta <- data.frame(
  sample_id = colnames(logCPM),
  donor = str_extract(colnames(logCPM), "^[A-Z]+\\d+"),
  visit_group = str_extract(colnames(logCPM), "V\\d+"),
  condition = ifelse(str_detect(colnames(logCPM), "^PTLD"), "PTLD", "RTH")
)

write.csv(meta, "./Results/061725/061725_meta.csv", row.names = FALSE)






#Step 4: Run Linear Mixed Models (per gene)

library(lme4)
library(broom.mixed)

results_list <- list()

for (gene in rownames(logCPM)) {
  df <- data.frame(
    expr = logCPM[gene, ],
    donor = meta$donor,
    visit_group = meta$visit_group,
    condition = meta$condition
  )
  
  model <- tryCatch(
    lmer(expr ~ visit_group * condition + (1 | donor), data = df),
    error = function(e) NULL
  )
  
  if (!is.null(model)) {
    tidy_model <- broom.mixed::tidy(model)
    tidy_model$gene <- gene
    results_list[[gene]] <- tidy_model
  }
}

model_results <- do.call(rbind, results_list)

#Step 5: Extract Significant Terms (e.g., interactions)

library(dplyr)

top_hits <- model_results %>%
  filter(term %in% c("conditionPTLD",
                     "visit_groupV3:conditionPTLD",
                     "visit_groupV5:conditionPTLD",
                     "visit_groupV7:conditionPTLD")) %>%
  arrange(p.value)

head(top_hits, 10)

write.csv(model_results, "./061725_pseudobulk_LMM_results.csv", row.names = FALSE)






# PCA
pca_res <- prcomp(t(logCPM), scale. = TRUE)
pca_df <- as.data.frame(pca_res$x)
pca_df$sample_id <- rownames(pca_df)
pca_df$donor <- str_extract(pca_df$sample_id, "^[A-Z]+\\d+")
pca_df$visit <- str_extract(pca_df$sample_id, "V\\d+")

ggplot(pca_df, aes(x = PC1, y = PC2, color = visit, label = donor)) +
  geom_point(size = 3) +
  geom_text_repel() +
  theme_cowplot()


library(ggpubr)
comparisons <- list(c("PTLD", "RTH"))
library(ggplot2)
library(ggpubr)
library(cowplot)
library(ggsci)

count2.3 <- count2.3 %>%
  rename(condition = conditon)

library(ggpubr)

# Define your comparison groups
comparisons <- list(c("PTLD", "RTH"))

# Boxplot with numeric p-values
p2.3 <- ggplot(count2.3, aes(x = condition, y = Percent, color = condition)) +
  geom_boxplot(alpha = 0.9, outlier.shape = NA) +
  geom_jitter(color = "black", size = 0.4, alpha = 0.9, width = 0.2) +
  facet_wrap(~ Var2, scales = "free") +
  stat_compare_means(comparisons = comparisons, label = "p.format", method = "wilcox.test") +
  scale_color_aaas() +
  theme_cowplot()

p2.3

ggsave(filename = "./Results/061725percentage_PTLD_RTH.pdf", p2.3, width=12, height = 20, dpi = 300)

pbmc.clean <- PTLD_RTH
count2.3 <- data.frame(table(pbmc.clean$subject, pbmc.clean$CellType.1))
total.bysubject <- data.frame(table(pbmc.clean$subject))
colnames(total.bysubject)[2] <- "Total"
count2.3$conditon <- "HD"
count2.3$conditon[str_detect(count2.3$Var1, "PTLD")] <- "PTLD"
count2.3$conditon[str_detect(count2.3$Var1, "PTLDN")] <- "PTLDN"
count2.3$conditon[str_detect(count2.3$Var1, "RTH")] <- "RTH"
count2.3 <- merge(count2.3, total.bysubject, by="Var1")
count2.3$Percent <- count2.3$Freq/count2.3$Total

count2.3$SLICE <- "SLICE_I"
count2.3$SLICE[str_detect(count2.3$conditon, "PTLDN")] <- "SLICE_III"
count2.3$SLICE[str_detect(count2.3$conditon, "HD")] <- "SLICE_III"

count2.3$conditon <- factor(count2.3$conditon, levels = c("HD", "PTLDN", "RTH", "PTLD"))

#count2.3 <- filter(count2.3, ! Var1 %in% c("PTLD2V1","RTH1V7b5","RTH2V5b5","RTH5V3b2"))     



comparasions <- list("RTH","PTLD")
p2.3 <- ggplot(count2.3, aes(x=conditon, y=Percent, color=conditon)) + geom_boxplot(alpha=0.9)
p2.3 <- p2.3 + facet_wrap(.~Var2, scales = "free") + theme_cowplot() + scale_color_aaas()
p2 <- p2 + geom_jitter(color="black", size=0.4, alpha=0.9)
p2.3 <- p2.3 + stat_compare_means(comparisons = comparasions)
p2.3

library(dplyr)
library(stringr)

count2.3 <- count2.3 %>%
  rename(condition = conditon) %>%  # fix the column name
  mutate(
    donor = str_extract(Var1, "^[A-Z]+\\d+"),
    visit = str_extract(Var1, "V\\d+[a-z0-9]*"),
    visit_group = str_replace(visit, "b[0-9]+", "")  # V3b2 → V3
  ) %>%
  filter(condition %in% c("PTLD", "RTH"))  # include only relevant groups

library(lme4)
library(broom.mixed)
#install.packages("broom.mixed")

celltypes <- unique(count2.3$Var2)

results <- lapply(celltypes, function(ct) {
  df <- filter(count2.3, Var2 == ct)
  
  if (n_distinct(df$donor) >= 3 && n_distinct(df$condition) > 1) {
    model <- tryCatch(
      lmer(Percent ~ condition * visit_group + (1 | donor), data = df),
      error = function(e) NULL
    )
    
    if (!is.null(model)) {
      out <- tidy(model) %>%
        filter(effect == "fixed") %>%
        mutate(celltype = ct)
      return(out)
    }
  }
  return(NULL)
})

model_results_condition <- do.call(rbind, results)

library(tidyr)

model_results_condition %>%
  filter(term %in% c("conditionPTLD",
                     "conditionPTLD:visit_groupV3",
                     "conditionPTLD:visit_groupV5",
                     "conditionPTLD:visit_groupV7")) %>%
  arrange(p.value) %>%
  select(celltype, term, estimate, std.error, p.value)

colnames(model_results_condition)
install.packages("lmerTest")
library(lmerTest)
model <- tryCatch(
  lmerTest::lmer(Percent ~ condition * visit_group + (1 | donor), data = df),
  error = function(e) NULL
)

library(lmerTest)  # instead of lme4
library(broom)

results <- lapply(celltypes, function(ct) {
  df <- filter(count2.3, Var2 == ct)
  
  if (n_distinct(df$donor) >= 3 && n_distinct(df$condition) > 1) {
    model <- tryCatch(
      lmer(Percent ~ condition * visit_group + (1 | donor), data = df),
      error = function(e) NULL
    )
    
    if (!is.null(model)) {
      out <- broom::tidy(model) %>%
        filter(effect == "fixed") %>%
        mutate(celltype = ct)
      return(out)
    }
  }
  return(NULL)
})

model_results_condition <- do.call(rbind, results)

model_results_condition %>%
  filter(term %in% c("conditionPTLD",
                     "conditionPTLD:visit_groupV3",
                     "conditionPTLD:visit_groupV5",
                     "conditionPTLD:visit_groupV7")) %>%
  arrange(p.value) %>%
  select(celltype, term, estimate, std.error, p.value)

write.csv(model_results_condition, "./Results/061725model_results_condition.csv", row.names = FALSE)

library(ggplot2)
library(dplyr)

# Filter interaction terms only
interaction_terms <- model_results_condition %>%
  filter(grepl("conditionPTLD:visit_group", term)) %>%
  filter(!is.na(p.value)) %>%
  arrange(p.value) %>%
  mutate(significant = p.value < 0.05)

# Plot top 20
top_n_interactions <- interaction_terms %>%
  slice_min(p.value, n = 20)

ggplot(top_n_interactions, aes(x = reorder(celltype, estimate), y = estimate, color = significant)) +
  geom_point(size = 3) +
  coord_flip() +
  facet_wrap(~term, scales = "free_y") +
  labs(title = "Top 20 Significant Interaction Effects: Condition × Visit",
       x = "Cell Type", y = "Effect Size (Estimate)") +
  scale_color_manual(values = c("TRUE" = "#8C1515", "FALSE" = "gray")) +
  theme_minimal()


##########ploe longirudinal expression for top genes:
library(tidyr)
library(dplyr)
library(ggplot2)
library(cowplot)

# Convert logCPM matrix to long format
logCPM_df <- as.data.frame(logCPM) %>%
  tibble::rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample_id", values_to = "logCPM") %>%
  left_join(meta, by = "sample_id")

# Select top genes from model results (or use any list you want)
top_genes <- model_results %>%
  filter(term == "conditionPTLD") %>%
  arrange(p.value) %>%
  distinct(gene) %>%
  slice_head(n = 4) %>%
  pull(gene)

# Plot expression over visit for top genes
for (g in top_genes) {
  p <- logCPM_df %>%
    filter(gene == g) %>%
    ggplot(aes(x = visit_group, y = logCPM, group = donor, color = condition)) +
    geom_point() +
    geom_line() +
    labs(title = paste("Longitudinal Expression:", g),
         x = "Visit", y = "logCPM") +
    facet_wrap(~gene) +
    theme_cowplot()
  
  print(p)
}

model_results %>%
  filter(term == "conditionPTLD") %>%
  arrange(p.value) %>%
  select(gene, estimate, p.value) %>%
  head(10)
any(!is.na(model_results$p.value))
unique(model_results$term)

library(dplyr)
library(ggplot2)

# Count samples per donor × visit × condition
sample_counts <- meta %>%
  count(donor, visit_group, condition)

# Bubble plot of sample counts
p.count <- ggplot(sample_counts, aes(x = visit_group, y = donor, size = n, color = condition)) +
  geom_point(alpha = 0.8) +
  scale_size_continuous(range = c(3, 10)) +
  theme_minimal(base_size = 14) +
  labs(title = "Sample Count per Donor × Visit × Condition",
       x = "Visit", y = "Donor", size = "n cells")

# Save the plot
ggsave("./Results/061725/sample_count_per_donor.pdf", p.count, width = 7, height = 5, dpi = 300)

##########Find genes where expression differs over time between PTLD and RTH, stratified by cell type.

expr ~ visit_group * condition + (1 | donor)

###Step 1: Get top differential genes per cell type
library(dplyr)

# Filter meaningful model terms
top_genes_all <- model_results %>%
  filter(term %in% c("conditionRTH", 
                     "visit_groupV3:conditionRTH",
                     "visit_groupV5:conditionRTH",
                     "visit_groupV7:conditionRTH")) %>%
  filter(p.value < 0.05) %>%
  arrange(p.value)

# View top
top_genes_all %>% group_by(celltype) %>% slice_min(p.value, n = 5)

top_genes_all$gene
# [1] "KLRC1" "KLRC3" "SLAMF1" "CD68" "S100A10" "CXCR4" "PRDM1" "LY86" "CXCL1"

library(ggplot2)
library(cowplot)

# We'll assume meta has: sample_id, donor, condition, visit_group, celltype
# and logCPM has genes as rows, sample IDs as columns

# Filter genes present in logCPM
genes_to_plot <- intersect(top_genes_all$gene, rownames(logCPM))

# Output directory
dir.create("./Results/061725/gene_plots", showWarnings = FALSE)

for (gene in genes_to_plot) {
  
  # Match expression to metadata
  expr_vals <- logCPM[gene, match(meta$sample_id, colnames(logCPM))]
  
  plot_df <- meta %>%
    mutate(expr = expr_vals,
           gene = gene)
  
  # Plot longitudinal expression
  p <- ggplot(plot_df, aes(x = visit_group, y = expr, group = donor, color = condition)) +
    geom_line(alpha = 0.6) +
    geom_point(size = 2) +
    theme_cowplot() +
    labs(title = paste("Longitudinal Expression:", gene),
         x = "Visit", y = "logCPM") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save to PDF
  ggsave(filename = paste0("./Results/061725/gene_plots/", gene, "_trajectory.pdf"),
         plot = p, width = 6, height = 4)
}

#######################################Δ-expression (RTH − PTLD) over time for the top variable genes
meta <- meta[!duplicated(meta$sample_id), ]
logCPM <- logCPM[, !duplicated(colnames(logCPM))]

library(matrixStats)
var_genes <- rowVars(as.matrix(logCPM))
top_genes <- names(sort(var_genes, decreasing = TRUE))[1:10]

library(reshape2)
logCPM_long <- melt(logCPM[top_genes, ])
colnames(logCPM_long) <- c("gene", "sample_id", "expr")

library(dplyr)
df <- left_join(logCPM_long, meta, by = "sample_id")

delta_expr <- df %>%
  group_by(gene, visit_group, condition) %>%
  summarise(mean_expr = mean(expr, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = condition, values_from = mean_expr) %>%
  mutate(delta = RTH - PTLD)

library(ggplot2)

ggplot(delta_expr, aes(x = visit_group, y = delta, color = gene, group = gene)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  labs(title = "Δ-Expression (RTH − PTLD) Over Time",
       x = "Visit",
       y = "Δ logCPM") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

