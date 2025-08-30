#script for patient symptome data
library("RColorBrewer")

data_PTLDN <- data.frame(
  PID = c("07-472", "07-473", "07-477", "07-478", "07-483", "07-484", "07-485"),
  gendr_dg = c(1, 1, 1, 1, 0, 1, 1),
  ageyr_dg = c(68, 69, 22, 44, 48, 52, 66),
  time_en = c(308, 1140, 262, 2469, 5096, 281, 791),
  serogrp = c(2, 2, 2, 0, 2, 0, 2),
  sero_mr = c(1, 1, 1, 1, 1, 1, 1),
  ax_cur = c(1, 1, 0, 0, 0, 0, 0),
  rash_mr = c(0, 0, 0, 0, 0, 0, 0),
  neuro_mr = c(1, 0, 0, 0, 0, 0, 0),
  arthritis_mr = c(0, 0, 0, 0, 0, 0, 0),
  carditis_mr = c(0, 0, 0, 0, 0, 0, 0),
  flulike_mr = c(1, 1, 1, 1, 1, 1, 1),
  tot_plqs = c(9, 6, 4, 4, 10, 13, 3),
  neuro_plqs = c(3, 3, 3, 3, 3, 6, 0),
  pcs = c(31.27, 42.95, 47.46, 43.95, 20.41, 30.25, 42.93),
  mcs = c(47.47, 61.3, 50.86, 40.58, 36.91, 42.04, 55.58)
)

data_PTLD <- data.frame(
  pid = c("01-009", "01-009", "01-009", "01-009", "01-009", "01-023", "01-023", "01-023", "01-025", "01-025", "01-045", "01-045", "01-045", "01-045"),
  visit = c(1, 2, 3, 5, 7, 2, 3, 5, 3, 7, 2, 3, 5, 7),
  age = c(64.7, 64.7, 64.8, 65.3, 66.7, 38.5, 38.6, 39, 54, 55.9, 42.1, 42.1, 42.6, 44),
  gendr_dg = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1),
  durat_days = c(3, 3, 3, 3, 3, 4, 4, 4, 42, 42, 3, 3, 3, 3),
  dissem = c(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0),
  serogrp = c(1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1),
  numb_sx = c(12, 17, 14, 9, 11, 24, 19, 23, 10, 9, 18, 9, 5, 7),
  pcs = c(49.18, 43.42, 37.54, 44.58, 48.17, 11.17, 18.33, 21.45, 43.21, 49.88, 32.99, 51.27, 60.85, 58.08),
  mcs = c(60.63, 39.18, 52.79, 38.02, 56.16, 50.84, 33, 36.14, 36, 26.16, 43.55, 41.66, 35.45, 45.7),
  time_en = c(0, 21, 30, 182, 720, 21, 30, 182, 30, 720, 21, 30, 182, 720),
  time = c(3, 24, 33, 185, 723, 25, 34, 186, 72, 762, 24, 33, 185, 723),
  group = rep("PTLD", 14)
)

rth <- data.frame(
  pid = c("01-037", "01-037", "01-037", "01-053", "01-053", "01-053", "01-063", "01-063", "01-063", "01-063", "01-071", "01-071", "01-071", "01-071", "01-082", "01-082", "01-082"),
  visit = c(3, 5, 7, 2, 3, 7, 2, 3, 5, 7, 2, 3, 5, 7, 2, 5, 7),
  age = c(28.6, 29, 30.6, 59, 59.1, 61.1, 51.1, 51.2, 51.6, 53.1, 55.8, 55.8, 56.3, 57.8, 44.9, 45.3, 46.9),
  gendr_dg = c(1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1),
  durat_days = c(3, 3, 3, 4, 4, 4, 13, 13, 13, 13, 8, 8, 8, 8, 4, 4, 4),
  dissem = c(0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
  serogrp = c(1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1),
  numb_sx = c(2, 2, 0, 11, 4, 4, 10, 0, 0, 0, 12, 0, 0, 1, 16, 4, 4),
  pcs = c(54.82, 57.94, 52.84, 42.08, 38.53, 44.25, 55.91, 60.67, 59.92, 60.55, 44.08, 56.86, 57.16, 50.46, 44.63, 56.8, 55.13),
  mcs = c(60.32, 58.21, 60.39, 45.37, 55.32, 49.6, 50.16, 53.58, 55.27, 53.9, 55.62, 57.22, 56.77, 60.77, 48.77, 57.32, 57.01),
  time_en = c(30, 182, 720, 21, 30, 720, 21, 30, 182, 720, 21, 30, 182, 720, 21, 182, 720),
  time = c(33, 185, 723, 25, 34, 724, 34, 43, 195, 733, 29, 38, 190, 728, 25, 186, 724),
  group = c(rep("RTH", 17))
)

rth_data_ID <- data.frame(
  RTH_Patient = c("RTH1", "RTH2", "RTH3", "RTH4", "RTH5", "RTH6", "PTLD1", "PTLD2", "PTLD3", "PTLD4"),
  PID = c("01-071", "01-082", "01-081", "01-063", "01-053", "01-037","01-023", "01-009", "01-045", "01-025")
)

ptldn_data <- data.frame(
  PTLDN_Patient = c("PTLDN1", "PTLDN2", "PTLDN3", "PTLDN4", "PTLDN5", "PTLDN6", "PTLDN7", "PTLDN8"),
  PID = c("07-483", "07-484", "07-485", "07-478", "07-472", "07-473", "07-459", "07-479")
)


data_long <- rbind(data_PTLD,rth)

data_long$time_en <- factor(data_long$time_en, levels = c("0", "21", "30", "182", "720"))

data_long <- left_join(data_long, rth_data_ID, by = c("pid"="PID"))

last_point <- data_long %>%
  group_by(RTH_Patient) %>% 
  slice(which.max(time_en))

blue_palette <- brewer.pal(9, "Blues")[4:9] # select from 4th to allow for difference in shades
red_palette <- brewer.pal(9, "Reds")[4:8] # select from 4th to allow for difference in shades

patient_colors <- c(blue_palette, red_palette)
names(patient_colors) <- c(paste0("RTH", 1:6), paste0("PTLD", 1:4))

rth.ptld.plot<- filter(data_long, time_en != 0) %>% ggplot(aes(x=time_en, y=numb_sx, group=RTH_Patient, color=RTH_Patient)) + geom_point() +
  geom_line(linetype = "dashed") + geom_text_repel(data=last_point, aes(label=RTH_Patient), size=2) + theme_classic() + scale_color_manual(values = patient_colors)

rth.ptld.plot <- rth.ptld.plot + ggtitle("Slice 1 Study") + xlab("Days Post Infection") + ylab("Number of Symptoms")
rth.ptld.plot

#ggsave("./Results/Lyme.Symp.pdf",rth.ptld.plot, width = 4, height = 3, dpi = 300)
#ggsave("../lyme_disease/Manuscript/Figures/Figure1/Sympt.Pa.pdf",rth.ptld.plot, width = 4, height = 3, dpi = 300)

data_long %>%
  group_by(RTH_Patient) %>%
  summarise(
    min_numb_sx = min(numb_sx),
    max_numb_sx = max(numb_sx),
    diff_numb_sx = max_numb_sx - min_numb_sx
  )
