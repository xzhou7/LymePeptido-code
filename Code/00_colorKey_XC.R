###############################
#Color code

condition_color <- c("PTLD" = "#BC3C29FF",
                     "RTH" = "#0072B5FF",
                     "PTLDN" = "#E18727FF",
                     "HD" = "#20854EFF")


condition_color_number <- c("04.PTLD" = "#BC3C29FF",
                     "03.RTH" = "#0072B5FF",
                     "02.PTLDN" = "#E18727FF",
                     "01.HD" = "#20854EFF")

Mono_color = c(
  "S100A9_Mono" = ggsci::pal_d3()(n=9)[9], 
  "HLA_Mono" = ggsci::pal_d3()(n=9)[3], 
  "HyperInflam_Mono" = ggsci::pal_d3()(n=9)[4],
  "NAMPT_Mono" = ggsci::pal_d3()(n=9)[2],
  "CD16_Mono" = ggsci::pal_d3()(n=9)[6],
  "LEF1_Mono" = ggsci::pal_d3()(n=9)[5],
  "IL32_Mono" = ggsci::pal_d3()(n=9)[8],
  "cDC" = ggsci::pal_d3()(n=9)[7],
  "pDC" = ggsci::pal_d3()(n=9)[1]
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
"#19413E", 
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
  "cDC_Mono" = ggsci::pal_d3()(n=9)[7],
  "pDC_Mono" = ggsci::pal_d3()(n=9)[1]
)

XZ_flip_x_label <- function() {
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
}


