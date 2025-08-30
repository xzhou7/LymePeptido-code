#xztools


theme_scientific <- function(base_size = 14, base_family = "serif") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      plot.background = element_blank(),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      legend.key = element_blank(),
      axis.line.y.left  = element_line(color = "black"),
      axis.line.x.bottom = element_line(color = "black"),
      axis.text.y.left = element_text(size = base_size),
      axis.text.x.bottom = element_text(size = base_size),
      axis.title = element_text(size = base_size + 2),
      legend.position = "none", # Hide legend
      plot.title = element_text(hjust = 0.5, size = base_size + 4, face = "bold"), # Center title, make it bigger and bold
      axis.title.y.right = element_blank(), # Remove right y-axis title
      axis.ticks.y.right = element_blank(), # Remove right y-axis ticks
      axis.text.y.right = element_blank(), # Remove right y-axis text
    )
}


theme_umap_clean <- function(base_size = 14, base_family = "serif") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      plot.background = element_blank(),
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_blank(),
      legend.key = element_blank(),
      axis.line = element_blank(), # Remove axis lines
      axis.text = element_blank(), # Remove axis text (ticks)
      axis.ticks = element_blank(), # Remove axis ticks
      axis.title = element_blank(), # Remove axis titles
      legend.position = "none", # Hide legend
      plot.title = element_text(hjust = 0.5, size = base_size + 4, face = "bold") # Center title, make it bigger and bold
    )
}


T.NK.celltype_colors <- c(
  "Naive_CD4" = "#BCBD22FF",
  "CD56dim_NK" = "#FF7F0EFF",
  "Activated_CD4" = "#2CA02CFF",
  "Activated_CD4_2" =  "#a6cee3",
  "GZMK_CD8" = "#6a3d9a",
  "Naive_CD8" = "#FFBB78FF",
  "MAIT_T" = "#E377C2FF",
  "GammaDelta_T" = "#fb8072",
  "GZMB_CD8" = "#8C564BFF",
  "CD56bright_NK" = "#1F77B4FF",
  "TFH_CD4" = "#17BECFFF",
  "Treg" = "#D62728FF",
  "LEF1_NK" = "#98DF8AFF",
  "Proliferating_Cell" = "#7F7F7FFF",
  "IFN_NK" = "darkred",
  "Other" = "#C5B0D5FF"
)
# 
# Tonsil_celltype_colors <- c(
#   "Activated B Cell" = "#D62728FF",
#   "Cycling T Cell" = "#FFBB78FF",
#   "FCRL4_MBC" = "#2CA02CFF",
#   "GC Dark Zone" =  "#800080FF",
#   "GC Light Zone" ="#6a3d9a",
#   "GD T Cell" = "#1F77B4FF",
#   "IFN B Cell" = "darkred",
#   "MBC" = "#fb8072",
#   "Memory CD4 T Cell" = "#FF7F0EFF",
#   "Memory CD8 T Cell" = "#8C564BFF",
#   "Myeloid Cells" = "#17BECFFF",
#   "Naive B Cell" = "#a6cee3",
#   "Naive CD4 T Cell" = "#9467BDFF",
#   "Naive CD8 T Cell" = "#7F7F7FFF",
#   "Plasmablast_kappa" = "#BCBD22FF",
#   "Plasmablast_lambda" = "#E377C2FF",
#   "pre-BCR B Cell" = "#98DF8AFF"
# )
# 
# Tonsil_celltype_colors <- c(
#   "Activated B Cell" = "#D62728FF",
#   "Cycling T Cell" = "#FF7F00FF",  # darker, more saturated orange
#   "FCRL4_MBC" = "#2CA02CFF",
#   "GC Dark Zone" =  "#800080FF",  # already adjusted
#   "GC Light Zone" ="#483D8BFF",  # darker blue (DarkSlateBlue)
#   "GD T Cell" = "#1F77B4FF",
#   "IFN B Cell" = "#e31a1c",
#   "MBC" = "#8B0000FF",  # darker red (DarkRed)
#   "Memory CD4 T Cell" = "#FF7F0EFF",
#   "Memory CD8 T Cell" = "#8C564BFF",
#   "Myeloid Cells" = "#104E8BFF",  # darker blue (DodgerBlue4)
#   "Naive B Cell" = "#006400FF",  # darker green (DarkGreen)
#   "Naive CD4 T Cell" = "#8A2BE2FF",  # more saturated purple (BlueViolet)
#   "Naive CD8 T Cell" = "#7F7F7FFF",
#   "Plasmablast_kappa" = "#228B22FF",  # more saturated green (ForestGreen)
#   "Plasmablast_lambda" = "#E377C2FF",
#   "pre-BCR B Cell" = "#8B4513FF"  # darker brown (SaddleBrown)
# )

Tonsil_celltype_colors <- c(
  "Activated B Cell" = "#1b9e77",
  "Cycling T Cell" = "#66a61e",
  "FCRL4_MBC" = "#7570b3",
  "GC Dark Zone" = "#e7298a",
  "GC Light Zone" = "#f781bf",
  "GD T Cell" = "#e6ab02",
  "IFN B Cell" =  "#4daf4a",
  "MBC" ="#ff7f00",
  "Memory CD4 T Cell" =  "#e31a1c",
  "Memory CD8 T Cell" = "#17BECFFF",
  "Myeloid Cells" ="#a6761d",
  "Naive B Cell" = "#999999",
  "Naive CD4 T Cell" =  "#104E8BFF",
  "Naive CD8 T Cell" = "#800080FF",
  "Plasmablast_kappa" = "#a65628",
  "Plasmablast_lambda" = "#d95f02",
  "pre-BCR B Cell" = "#984ea3"
)
# 
# Tonsil_celltype_colors <- c(
#   "Activated B Cell" = "#D62728FF",
#   "Cycling T Cell" = "#d95f02",
#   "FCRL4_MBC" = "#2CA02CFF",
#   "GC Dark Zone" =  "#800080FF",
#   "GC Light Zone" = "#66a61e",
#   "GD T Cell" = "#1F77B4FF",
#   "IFN B Cell" = "darkred",
#   "MBC" = "#8B0000FF",
#   "Memory CD4 T Cell" = "#e31a1c",
#   "Memory CD8 T Cell" = "#8C564BFF",
#   "Myeloid Cells" = "#4daf4a",
#   "Naive B Cell" = "#006400FF",
#   "Naive CD4 T Cell" = "#8A2BE2FF",
#   "Naive CD8 T Cell" = "#7F7F7FFF",
#   "Plasmablast_kappa" = "#a65628",
#   "Plasmablast_lambda" = "#E377C2FF",
#   "pre-BCR B Cell" = "#8B4513FF"
# )

tonsil_origin_colors <- c(
  "Tonsil" = "#66c2a5",   # green from d3 color palette
  "Tonsil_Organoid" = "#377eb8",  # blue from d3 color palette
  "Organoid_LAIV" = "#d7191c"  # red from d3 color palette
)





