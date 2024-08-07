## Supplemental Figure 4

## supplemental_figure4.R: generate supp figure 4 of SCFA paper
## Author: Andrew Oliver
## Date: Jul 16, 2024
## docker run --rm -it \
## -v ~/Downloads/SCFA-Analysis/figure_scripts:/home/scripts \
## -v ~/Downloads/SCFA-Analysis-DATA/data/:/home/data \
## -w /home/docker \
## scfa_analysis:rstudio bash -c "Rscript supplemental_figure4.R"

## load libraries ==============================================================
library(dplyr)
library(ggpubr)
library(ggplot2)
library(ggh4x)
library(wesanderson)
library(patchwork)

## set working directory =======================================================
setwd("/home/")

## source SCFA data ============================================================
set.seed(123)
source("/home/scripts/pre_process_raw_scfas.R")

## wrangle data or source helper functions =====================================
inflammation_markers <- merge(readr::read_csv("/home/data/CTSC24532USDAWHNRCNu-GIMarkers7Oct2021_DATA_2021-10-07_1627.csv"), 
                              readr::read_csv("/home/data/CRP_WBC_9102021.csv"), 
                              by = "subject_id") %>%
  dplyr::select(., dplyr::any_of(c("subject_id", "fecal_calprotectin", "crp_bd1"))) %>%
  merge(., fecal_scfas, by = "subject_id", all = T) %>%
  merge(., scfa_plasma_dedup, by = "subject_id", all = T) %>%
  dplyr::mutate(., crp_bd1_cutoff = ifelse(crp_bd1 >= 10000, "high", "low")) %>%
  dplyr::mutate(., fecal_cal_cutoff = ifelse(fecal_calprotectin >= 100, "high", "low")) %>%
  reshape2::melt(., id.vars = c("subject_id", "fecal_cal_cutoff", "crp_bd1_cutoff")) %>%
  dplyr::filter(., variable %in% c("acetate","propionate", 
                                   "butyrate", 
                                   "total_scfa", "acetate_norm_ratio_dist", 
                                   "propionate_norm_ratio_dist", 
                                   "butyrate_norm_ratio_dist", "p_acetic_acid_nmol",
                                   "p_propionic_acid_nmol", "p_butyric_acid_nmol", 
                                   "p_scfa_nmol_total")) 

## clean up data, order factors
inflammation_markers$fecal_cal_cutoff <- factor(x = inflammation_markers$fecal_cal_cutoff, levels = c("low", "high"), ordered = T)
inflammation_markers$crp_bd1_cutoff <- factor(x = inflammation_markers$crp_bd1_cutoff, levels = c("low", "high"), ordered = T)
inflammation_markers$variable <- factor(x = inflammation_markers$variable, levels = c("acetate", 
                                                                                      "acetate_norm_ratio_dist", "propionate", 
                                                                                      "propionate_norm_ratio_dist", "butyrate", 
                                                                                      "butyrate_norm_ratio_dist", 
                                                                                      "total_scfa", "p_acetic_acid_nmol", "p_propionic_acid_nmol", "p_butyric_acid_nmol", "p_scfa_nmol_total"), ordered = T)
## create an overall SCFA factor for ratio and raw
inflammation_markers <- inflammation_markers %>%
  dplyr::mutate(overall_scfa = case_when(
    variable == "acetate" ~ 'fecal acetate',
    variable == "acetate_norm_ratio_dist" ~ 'fecal acetate',
    variable == "propionate" ~ 'fecal propionate',
    variable == "propionate_norm_ratio_dist" ~ 'fecal propionate',
    variable == "butyrate" ~ 'fecal butyrate',
    variable == "butyrate_norm_ratio_dist" ~ 'fecal butyrate',
    variable == "total_scfa" ~ 'fecal total SCFAs',
    variable == "p_acetic_acid_nmol" ~ 'plasma acetate',
    variable == "p_propionic_acid_nmol" ~ 'plasma propionate',
    variable == "p_butyric_acid_nmol" ~ 'plasma butyrate',
    variable == "p_scfa_nmol_total" ~ 'plasma total SCFAs',
    TRUE ~ 'what_is_this'      
  )
  )

## create a factor for color palletes
inflammation_markers <- inflammation_markers %>%
  dplyr::mutate(twin_scfas = case_when(
    variable == "acetate" ~ 'acetate',
    variable == "acetate_norm_ratio_dist" ~ 'acetate',
    variable == "propionate" ~ 'propionate',
    variable == "propionate_norm_ratio_dist" ~ 'propionate',
    variable == "butyrate" ~ 'butyrate',
    variable == "butyrate_norm_ratio_dist" ~ 'butyrate',
    variable == "total_scfa" ~ 'total_fecal',
    variable == "p_acetic_acid_nmol" ~ 'acetate',
    variable == "p_propionic_acid_nmol" ~ 'propionate',
    variable == "p_butyric_acid_nmol" ~ 'butyrate',
    variable == "p_scfa_nmol_total" ~ 'total_plasma',
    TRUE ~ 'what_is_this'      
  )
  )

## color pallete based on SCFA
inflammation_markers <- inflammation_markers %>%
  dplyr::mutate(color_fills = case_when(
    twin_scfas == "acetate" ~ '#0A9F9D',
    twin_scfas == "propionate" ~ '#CEB175',
    twin_scfas == "butyrate" ~ '#E54E21',
    twin_scfas == "total_fecal" ~ '#634F25',
    twin_scfas == "total_plasma" ~ '#D71515',
    TRUE ~ 'what_is_this'      
  )
  )

## create a ratio/raw factor
inflammation_markers <- inflammation_markers %>%
  dplyr::mutate(ratio_raw = case_when(
    variable == "acetate" ~ 'raw',
    variable == "acetate_norm_ratio_dist" ~ 'ratio',
    variable == "propionate" ~ 'raw',
    variable == "propionate_norm_ratio_dist" ~ 'ratio',
    variable == "butyrate" ~ 'raw',
    variable == "butyrate_norm_ratio_dist" ~ 'ratio',
    variable == "total_scfa" ~ 'raw',
    variable == "p_acetic_acid_nmol" ~ 'raw',
    variable == "p_propionic_acid_nmol" ~ 'raw',
    variable == "p_butyric_acid_nmol" ~ 'raw',
    variable == "p_scfa_nmol_total" ~ 'raw',
    TRUE ~ 'what_is_this'      
  )
  )

## order SCFAs
inflammation_markers$twin_scfas <- factor(x = inflammation_markers$twin_scfas, levels = c("acetate","propionate","butyrate", "total_fecal","total_plasma"), ordered = T)
## order ratios/raw
inflammation_markers$ratio_raw <- factor(x = inflammation_markers$ratio_raw, levels = c("raw", "ratio"), ordered = T)
## create DF for inflammation distribution
inflammation_distribution <- merge(readr::read_csv("/home/data/CTSC24532USDAWHNRCNu-GIMarkers7Oct2021_DATA_2021-10-07_1627.csv"), 
                              readr::read_csv("/home/data/CRP_WBC_9102021.csv"), 
                              by = "subject_id") %>%
  dplyr::select(., dplyr::any_of(c("subject_id", "fecal_calprotectin", "crp_bd1"))) %>%
  merge(., fecal_scfas, by = "subject_id", all = T) %>%
  merge(., scfa_plasma_dedup, by = "subject_id", all = T) %>%
  dplyr::filter(., subject_id %in% c(fecal_scfas$subject_id, scfa_plasma_dedup$subject_id)) %>%
  dplyr::select(., subject_id, crp_bd1, fecal_calprotectin)

## make supp figure 4A =========================================================

fecal_cal_density <- ggplot(inflammation_distribution, aes(x=fecal_calprotectin)) + 
  geom_density(alpha=.3, fill = "#DCAA79") +
  theme_bw() +
  geom_vline(xintercept = 100, linetype = "dashed", color = "red") +
  labs(x = "Fecal Calprotectin (ug/g)")

crp_density <- ggplot(inflammation_distribution, aes(x=crp_bd1)) + 
  geom_density(alpha=.3, fill = "#FFB3B3") +
  theme_bw() +
  geom_vline(xintercept = 10000, linetype = "dashed", color = "red") +
  labs(x = "C-reactive Protein (ug/L)")

## make supp figure 4B =========================================================
sup_figure4B <- inflammation_markers %>% 
  dplyr::filter(., crp_bd1_cutoff != "NA") %>%
  ggplot() + 
  aes(x = crp_bd1_cutoff, y = value) +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.4) + 
  geom_boxplot(outlier.alpha = 0, aes(alpha = crp_bd1_cutoff, fill = twin_scfas)) + 
  facet_nested(.~ overall_scfa + ratio_raw, scales = "free", independent = "all") +
  ggpubr::stat_pwc(hjust = 0.5, vjust = 1, p.adjust.by = "panel") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.background = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = wesanderson::wes_palette("AsteroidCity1")) +
  scale_alpha_manual(values = c(1, 0.4)) +
  labs(x = "C-reactive protein", y = "SCFA abundance\n(raw: nmol/mg or ul, ratio: Δ expected ratio)")
  

## make supp figure 4C =========================================================
sup_figure4C <- inflammation_markers %>% 
  dplyr::filter(., fecal_cal_cutoff != "NA") %>%
  ggplot() + 
  aes(x = fecal_cal_cutoff, y = value) +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.4) + 
  geom_boxplot(outlier.alpha = 0, aes(alpha = fecal_cal_cutoff, fill = twin_scfas)) + 
  facet_nested(.~ overall_scfa + ratio_raw, scales = "free", independent = "all") +
  ggpubr::stat_pwc(hjust = 0.5, vjust = 1) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.background = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = wesanderson::wes_palette("AsteroidCity1")) +
  scale_alpha_manual(values = c(1, 0.4)) +
  labs(x = "Fecal calprotectin", y = "SCFA abundance\n(raw: nmol/mg or ul, ratio: Δ expected ratio)")


design <- "A
B
C
D"

fecal_cal_density + crp_density + sup_figure4B + sup_figure4C + patchwork::plot_layout(design = design) + patchwork::plot_annotation(tag_levels = c("A"))

dir.create(path = "/home/scripts/output_figures", showWarnings = TRUE)
ggsave(filename = "/home/scripts/output_figures/supplemental_figure4.png", 
       device = "png",
       width = 12, 
       height = 12,
       units = "in",
       dpi = 400)
dev.off()
