###############################
#Color code

condition_color <- c("PTLD" = "#BC3C29FF",
                     "RTH" = "#0072B5FF",
                     "PTLDN" = "#E18727FF",
                     "HD" = "#20854EFF")

condition_color_spleen <- c("Peptidoglycan" = "#BC3C29FF",
                     "Control" = "#0072B5FF")

condition_color_number <- c("04.PTLD" = "#BC3C29FF",
                     "03.RTH" = "#0072B5FF",
                     "02.PTLDN" = "#E18727FF",
                     "01.HD" = "#20854EFF")

Mono_color1 = c(
  "S100A9_Mono" = ggsci::pal_d3()(n=9)[9], 
  "HLA_Mono" = ggsci::pal_d3()(n=9)[3], 
  "HyperInflam_Mono" = ggsci::pal_d3()(n=9)[4],
  "NAMPT_Mono" = ggsci::pal_d3()(n=9)[2],
  "CD16_Mono" = ggsci::pal_d3()(n=9)[6],
  "X5_Mono" = ggsci::pal_d3()(n=9)[5],
  "X6_Mono" = ggsci::pal_d3()(n=9)[8],
  "cDC" = ggsci::pal_d3()(n=9)[7],
  "pDC" = ggsci::pal_d3()(n=9)[1]
)

Mono_color <- c(
  "S100A9_Mono"      = "#BCBD22",  # 9
  "HLA_Mono"         = "#2CA02C",  # 3
  "HyperInflam_Mono"= "#D62728",  # 4
  "NAMPT_Mono"       = "#FF7F0E",  # 2
  "CD16_Mono"        = "#8C564B",  # 6
  "X5_Mono"          = "#9467BD",  # 5
  "X6_Mono"          = "#7F7F7F",  # 8
  "cDC"              = "#E377C2",  # 7
  "pDC"              = "#1F77B4"   # 1
)


stanford_red <- c("#8C1515", "#D3D3D3")

my_palette <- c("#3B99D4", 
"#8ED14B", 
"#F06B49",
"#ECC2F1", 
"#82C7C3", 
"#19413E",
"#1776EB", 
"#F5B2AC", 
"#533085",
"#89363A",
"#D92B45", 
"#60C9FF", 
"#1B9F2E", 
"#BA217D", 
"#635019", 
"#E3698A", 
"#076B82", 
"#A86A16")

my_palette_alpha_0.7 <- sapply(my_palette, function(color) {
  adjustcolor(color, alpha.f = 0.7)
}) %>% as.character()

cluster_color = c(
  "S100A9_Mono" = ggsci::pal_d3()(n=9)[9], 
  "HLA_Mono" = ggsci::pal_d3()(n=9)[3], 
  "HyperInflam_Mono" = ggsci::pal_d3()(n=9)[4],
  "NAMPT_Mono" = ggsci::pal_d3()(n=9)[2],
  "CD16_Mono" = ggsci::pal_d3()(n=9)[6],
  "X5_Mono" = ggsci::pal_d3()(n=9)[5],
  "X6_Mono" = ggsci::pal_d3()(n=9)[8],
  "cDC" = ggsci::pal_d3()(n=9)[7],
  "pDC" = ggsci::pal_d3()(n=9)[1]
)

XZ_flip_x_label <- function() {
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
}

pep_cluster2_color <- c("Mem_B" = "#1B9F2E",
                        "Naive_B" = "#82C7C3",
                        "GerminalCenter_B" = "#8ED14B",
                        "Plasmablast" = "#19413E",
                        
                        "Activated_NK" = "#BA217D",
                        "Rest_NK" = "#ECC2F1",
                        "Myeloid_cells" = "#A86A16",
                        
                        "GZMK_CD8T" = "#533085",
                        "GZMB_CD8T" = "#F06B49",
                        "Gamma_Delta_T" = "#E3698A",
                        "Naive_CD8T" = "#D92B45",
                        "Mem_CD8T" = "#89363A",
                        
                        "Naive_CD4T" = "#60C9FF",
                        "Mem_CD4" = "#076B82",
                        "Treg" = "#1776EB",
                        
                        "Undetermined" = "#D3D3D3"
                        )

Spleen_Org_Color <- c(
  "Spleen_Control" = "#1F77B4FF", # ggsci::pal_d3()(n=9)[1]
  "Spleen_Peptidoglycan" = "#FF7F0EFF",#ggsci::pal_d3()(n=9)[2]
  "LPS" = "#D62728FF") #ggsci::pal_d3()(n=9)[4])

spleen_myeloid_color <- c(
  "memory B cells"                        = "#BCBD22",  # 9
  "cytotoxic T cells"                     = "#2CA02C",  # 3
  "monocyte-derived inflam macrophages"  = "#D62728",  # 4
  "inflam tissue macrophages"            = "#FF7F0E",  # 2
  "anti-inflam tissue macrophages"       = "#8C564B",  # 6
  "differentiating monocytes"            = "#9467BD",  # 5
  "dendritic cells"                      = "#7F7F7F",  # 8
  "hematopoietic stem and progenitor cells" = "#E377C2",  # 7
  "proliferating myeloid cells"  = "#1F77B4"   # 1
)

