## Figure 1

## figure1.R: generate figure 1 of SCFA paper
## Author: Andrew Oliver
## Date: Feb 22, 2024
## to run: 

## docker run --rm -it \
## -v ~/Downloads/SCFA-Analysis/figure_scripts:/home/scripts \
## -v ~/Downloads/SCFA-Analysis-DATA/data/:/home/data \
## -w /home/docker \
## scfa_analysis:rstudio bash -c "Rscript Figure1.R"

## load libraries ==============================================================
library(dplyr)
library(ggpubr)
library(ggplot2)
library(wesanderson)
library(corrr)
library(patchwork)
library(cowplot)

## set working directory =======================================================
setwd("/home")

## source data =================================================================
source("/home/scripts/pre_process_raw_scfas.R")
anthropometrics <- read.csv("/home/data/FL100_age_sex_bmi.csv")

## wrangle data or source helper functions =====================================
source("/home/scripts/partial_regression.R")

## make plot 1A ================================================================
plot1a <- fecal_scfas %>%
  dplyr::select(., subject_id, acetate, propionate, butyrate) %>%
  dplyr::mutate(., new_total = (acetate + propionate + butyrate)) %>%
  reshape2::melt(., id.vars = c("subject_id", "new_total")) %>%
  ggplot() + aes(x = reorder(subject_id, -new_total) , weight = value) + 
  geom_bar(aes(fill = variable)) +
  scale_fill_manual(values = wesanderson::wes_palette("AsteroidCity1")) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(),
        legend.position = "none",
        panel.background = element_blank(),
        axis.text = element_text(colour = "black")) +
  labs(x = "", y = "Fecal\nSCFA Abundance (nmol/mg)") 

## make plot 1B ================================================================
plot1b_main <- scfa_plasma_dedup %>%
  dplyr::select(., subject_id, p_acetic_acid_nmol, p_propionic_acid_nmol, p_butyric_acid_nmol, p_scfa_nmol_total) %>%
  reshape2::melt(., id.vars = c("subject_id", "p_scfa_nmol_total")) %>%
  ggplot() + aes(x = reorder(subject_id, -p_scfa_nmol_total), weight = value) + 
  geom_bar(aes(fill = variable)) +
  scale_fill_manual(values = wesanderson::wes_palette("AsteroidCity1")) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(),
        legend.position = "none",
        panel.background = element_blank(),
        axis.text = element_text(colour = "black")) +
  labs(x = "", y = "Plasma\nSCFA Abundance (nmol/ul)") 

## make plot 1B (inset) ========================================================
plot1b_inset <- scfa_plasma_dedup %>%
  dplyr::select(., subject_id, p_propionic_acid_nmol, p_butyric_acid_nmol, p_scfa_nmol_total) %>%
  reshape2::melt(., id.vars = c("subject_id", "p_scfa_nmol_total")) %>%
  ggplot() + aes(x = reorder(subject_id, -p_scfa_nmol_total), weight = value) + 
  geom_bar(aes(fill = variable)) +
  scale_fill_manual(values = wesanderson::wes_palette("AsteroidCity1")[2-3]) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(),
        legend.position = "none",
        panel.background = element_blank(),
        axis.text = element_text(colour = "black")) +
  labs(x = "", y = " SCFA Abundance (nmol/ul)") 

plot1b <-
  ggdraw() +
  draw_plot(plot1b_main) +
  draw_plot(plot1b_inset, x = 0.68, y = .65, width = .3, height = .3)

## make plot 1D ================================================================
plot1d <- merge(fecal_scfas, scfa_plasma_dedup, by = "subject_id") %>% 
  dplyr::select(., acetate, propionate, butyrate, 
                p_acetic_acid_nmol, p_propionic_acid_nmol, 
                p_butyric_acid_nmol, total_scfa,p_scfa_nmol_total) %>% 
  corrr::correlate(method = "pearson") %>% corrr::autoplot()

## make plot 1E ================================================================
## fecal ph vs fecal scfas

fecal_ph_fecal_scfa <- merge(readr::read_delim("/home/data/FL100_stool_variables.txt"), 
                             fecal_scfas, by = "subject_id") %>% 
  merge(., anthropometrics, by = "subject_id")


set.seed(123)
fecal_ph_fecal_scfa$normalized <- bestNormalize::bestNormalize(fecal_ph_fecal_scfa$total_scfa, allow_orderNorm = F, k = 10, r = 10)$x.t
print(bestNormalize::bestNormalize(fecal_ph_fecal_scfa$total_scfa, allow_orderNorm = F, k = 10, r = 10)$chosen_transform)

model <- lm(normalized ~ as.numeric(fecal_ph) + as.numeric(bmi) + as.numeric(age) + as.factor(sex.factor) + as.numeric(bristol_num) + as.numeric(st_wt), data = fecal_ph_fecal_scfa)
partial_regression <- car::avPlots(model, main = paste("Total SCFAs vs. Fecal pH"))

tmp_fecal_ph <- cor.test(partial_regression$`as.numeric(fecal_ph)`[,1], partial_regression$`as.numeric(fecal_ph)`[,2])

plot1e <- ggplot(as.data.frame(partial_regression$`as.numeric(fecal_ph)`)) + 
  aes(x = `as.numeric(fecal_ph)`, y = normalized) + 
  geom_point() + 
  geom_smooth(method = "lm", color = "purple", se = TRUE, linetype = "dashed") + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), axis.text = element_text(colour = "black")) +
  labs(x = "Fecal pH (normalized) | covariates", y = "Total SCFA (normalized) | covariates")

## make plot 1C ================================================================
fecal_scfas_anthro <- merge(anthropometrics, fecal_scfas, by = "subject_id")
tmp_inner <- data.frame()
tmp_outer <- data.frame()

for (scfa in c("acetate","propionate", "butyrate")) {
  print(scfa)
  
  fecal_scfas_anthro$normalized <- bestNormalize::bestNormalize(fecal_scfas_anthro[[scfa]], allow_orderNorm = F, k = 10, r = 10)$x.t
  print(bestNormalize::bestNormalize(fecal_scfas_anthro[[scfa]], allow_orderNorm = F, k = 10, r = 10)$chosen_transform)
  
  model <- lm(normalized ~ as.numeric(bmi) + as.numeric(age) + as.factor(sex.factor), data = fecal_scfas_anthro)
  partial_regression <- car::avPlots(model, main = paste(scfa, "vs. Anthro"))
  
  tmp_bmi <- cor.test(partial_regression$`as.numeric(bmi)`[,1], partial_regression$`as.numeric(bmi)`[,2])
  tmp_inner_bmi <- data.frame(normalized_factor = partial_regression$`as.numeric(bmi)`[,1],
                              normalized_scfa = partial_regression$`as.numeric(bmi)`[,2],
                              scfa_name = rep(scfa, length(NROW(partial_regression$`as.numeric(bmi)`))),
                              factor_name = rep("bmi", length(NROW(partial_regression$`as.numeric(bmi)`))))
  
  tmp_age <- cor.test(partial_regression$`as.numeric(age)`[,1], partial_regression$`as.numeric(age)`[,2])
  tmp_inner_age <- data.frame(normalized_factor = partial_regression$`as.numeric(age)`[,1],
                              normalized_scfa = partial_regression$`as.numeric(age)`[,2],
                              scfa_name = rep(scfa, length(NROW(partial_regression$`as.numeric(age)`))),
                              factor_name = rep("age", length(NROW(partial_regression$`as.numeric(age)`))))
  
  tmp_sex <- cor.test(partial_regression$`as.factor(sex.factor)Male`[,1], partial_regression$`as.factor(sex.factor)Male`[,2])
  tmp_inner_sex <- data.frame(normalized_factor = partial_regression$`as.factor(sex.factor)Male`[,1],
                              normalized_scfa = partial_regression$`as.factor(sex.factor)Male`[,2],
                              scfa_name = rep(scfa, length(NROW(partial_regression$`as.factor(sex.factor)Male`))),
                              factor_name = rep("sex", length(NROW(partial_regression$`as.factor(sex.factor)Male`))))
  df_anthro <- data.frame(factor = c("BMI", "Age", "Sex"), 
                          corr_estimate = c(tmp_bmi$estimate, tmp_age$estimate, tmp_sex$estimate), 
                          p_value = c(tmp_bmi$p.value, tmp_age$p.value, tmp_sex$p.value))
  df_anthro$p_adjust <- p.adjust(df_anthro$p_value)
  tmp_outer <- rbind(tmp_outer, tmp_inner_bmi, tmp_inner_age, tmp_inner_sex)
  
  print(df_anthro)
}

tmp_outer$scfa_name <- factor(x = tmp_outer$scfa_name, levels = c("acetate", "propionate", "butyrate"), ordered = T)
plot1c_1 <- ggplot(subset(tmp_outer, tmp_outer$factor_name == "bmi")) + 
  aes(x = normalized_factor, y = normalized_scfa) + 
  geom_point(aes(color = scfa_name), alpha = 0.1) + 
  geom_smooth(method = "lm", aes(color = scfa_name), se = FALSE) + 
  scale_color_manual(values = wesanderson::wes_palette("AsteroidCity1")) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), axis.text = element_text(colour = "black")) +
  labs(x = "BMI | age,sex", y = "Fecal SCFA (normalized)")
plot1c_2 <- ggplot(subset(tmp_outer, tmp_outer$factor_name == "age")) + 
  aes(x = normalized_factor, y = normalized_scfa) + 
  geom_point(aes(color = scfa_name), alpha = 0.1, position = position_jitter()) + 
  geom_smooth(method = "lm", aes(color = scfa_name), se = FALSE) + 
  scale_color_manual(values = wesanderson::wes_palette("AsteroidCity1")) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), axis.text = element_text(colour = "black")) +
  labs(x = "AGE | bmi,sex", y = "Fecal SCFA (normalized)")
plot1c_3 <- tmp_outer %>% dplyr::filter(., factor_name == "sex") %>%
  dplyr::mutate(., sex_factor = ifelse(as.numeric(normalized_factor) > 0, "male", "female")) %>%
  ggplot() + 
  aes(x = sex_factor, y = normalized_scfa) + 
  geom_boxplot(aes(fill = scfa_name)) + 
  scale_fill_manual(values = wesanderson::wes_palette("AsteroidCity1")) + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), axis.text = element_text(colour = "black")) +
  labs(x = "SEX | bmi,age", y = "Fecal SCFA (normalized)")

design <- "
AAEFGH
BBEFGI
"

plot1a + plot1b + plot1c_1 + plot1c_2 + plot1c_3 + plot1d + plot1e + patchwork::plot_layout(design = design) + patchwork::plot_annotation(tag_levels = c("A", "B", "C", "C", "C", "D", "E"))
dir.create(path = "/home/scripts/output_figures", showWarnings = TRUE)
ggsave(filename = "/home/scripts/output_figures/Figure1.png", 
       device = "png",
       width = 20, 
       height = 8,
       units = "in",
       dpi = 400)
dev.off()
