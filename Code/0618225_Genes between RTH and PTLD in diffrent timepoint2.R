####check which gene expression differentally expressin between RTH and TPLD on V7


setwd("~/Library/CloudStorage/Dropbox/lyme_disease/R7_NR/")
getwd()

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


load("./Data/Step4.1_022823_rename_clean_cluster.RData")

selected_celltypes <- c(
  "X02_CD14mono", "X06_CD14mono", "X14_CD14mono", "X18_CD16mono", 
  "X22_IntermediateMono", "X24_CD14mono", "X28_CD14mono", 
  "X17_NK", "X16_CD14mono", "X15_CD4T", "X11_GZMK_CD8T", "X03_Treg"
)
selected_celltypes <- c("X16_CD14mono","X03_Treg")
library(Seurat)
library(edgeR)
library(dplyr)
library(stringr)
library(reshape2)

results_list <- list()

for (ct in selected_celltypes) {
  subset_obj <- subset(pbmc.clean, idents = ct)
  
  # Assign metadata
  subset_obj$donor <- str_extract(subset_obj$subject, "^[A-Z]+\\d+")
  subset_obj$visit <- str_extract(subset_obj$subject, "V\\d+[a-z0-9]*")
  subset_obj$visit_group <- str_replace(subset_obj$visit, "b[0-9]+", "")
  subset_obj$sample_id <- paste(subset_obj$donor, subset_obj$visit_group, sep = "_")
  
  # Aggregate counts per sample
  pb_mat <- AggregateExpression(subset_obj, group.by = "sample_id", assays = "RNA", slot = "counts")$RNA
  
  # Convert to logCPM
  dge <- DGEList(counts = pb_mat)
  keep_genes <- filterByExpr(dge)
  dge <- dge[keep_genes,, keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  logCPM <- cpm(dge, log = TRUE)
  
  # Metadata
  meta <- data.frame(
    sample_id = colnames(logCPM),
    donor = str_extract(colnames(logCPM), "^[A-Z]+\\d+"),
    visit_group = str_extract(colnames(logCPM), "V\\d+"),
    condition = ifelse(str_detect(colnames(logCPM), "^PTLD"), "PTLD", "RTH")
  )
  
  # Reshape for t-tests
  logCPM_df <- as.data.frame(logCPM)
  logCPM_df$gene <- rownames(logCPM_df)
  long <- melt(logCPM_df, id.vars = "gene", variable.name = "sample_id", value.name = "expr")
  df <- left_join(long, meta, by = "sample_id")
  
  # Run t-test for each gene
  stat_result <- df %>%
    group_by(gene) %>%
    summarise(
      p_val = tryCatch(t.test(expr ~ condition)$p.value, error = function(e) NA),
      delta = mean(expr[condition == "RTH"]) - mean(expr[condition == "PTLD"]),
      .groups = "drop"
    ) %>%
    mutate(p_adj = p.adjust(p_val, method = "fdr"),
           celltype = ct)
  
  results_list[[ct]] <- stat_result
}

all_results <- bind_rows(results_list)
top_hits <- all_results %>% filter(p_adj < 0.05) %>% group_by(celltype) %>% slice_max(abs(delta), n = 1)
print(top_hits)

library(ggplot2)

for (i in 1:nrow(top_hits)) {
  gene <- top_hits$gene[i]
  ct <- top_hits$celltype[i]
  
  subset_obj <- subset(pbmc.clean, idents = ct)
  subset_obj$donor <- str_extract(subset_obj$subject, "^[A-Z]+\\d+")
  subset_obj$visit <- str_extract(subset_obj$subject, "V\\d+[a-z0-9]*")
  subset_obj$visit_group <- str_replace(subset_obj$visit, "b[0-9]+", "")
  subset_obj$sample_id <- paste(subset_obj$donor, subset_obj$visit_group, sep = "_")
  
  pb_mat <- AggregateExpression(subset_obj, group.by = "sample_id", assays = "RNA", slot = "counts")$RNA
  dge <- DGEList(counts = pb_mat)
  keep_genes <- filterByExpr(dge)
  dge <- dge[keep_genes,, keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  logCPM <- cpm(dge, log = TRUE)
  
  expr_df <- data.frame(expr = logCPM[gene, ])
  expr_df$sample_id <- rownames(expr_df)
  expr_df <- left_join(expr_df, meta, by = "sample_id")
  
  p <- ggplot(expr_df, aes(x = visit_group, y = expr, color = condition)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.8) +
    theme_bw() +
    labs(title = paste0(ct, ": ", gene), y = "logCPM", x = "Visit")
  
  ggsave(paste0("./Results/061725/", ct, "_", gene, "_visit_boxplot.pdf"), plot = p, width = 6, height = 4)
}
