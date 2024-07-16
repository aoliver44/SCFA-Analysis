## Supplemental Figure 4

## supplemental_figure4.R: generate supp figure 4 of SCFA paper
## Author: Andrew Oliver
## Date: Jul 16, 2024
## to run: docker run --rm -it -p 8787:8787 -e PASSWORD=yourpasswordhere -v `PWD`:/home/docker aoliver44/scfa_rstudio:1.0

## load libraries ==============================================================
library(dplyr)
library(ggpubr)
library(ggplot2)
library(ggh4x)
library(wesanderson)
library(patchwork)

## set working directory =======================================================
setwd("/home/docker")

## source SCFA data ============================================================
source("/home/docker/github/SCFA-Analysis/figure_scripts/pre_process_raw_scfas.R")

## wrangle data or source helper functions =====================================
inflammation_markers <- merge(readr::read_csv("/home/docker/data/CTSC24532USDAWHNRCNu-GIMarkers7Oct2021_DATA_2021-10-07_1627.csv"), 
                              readr::read_csv("/home/docker/data/CRP_WBC_9102021.csv"), 
                              by = "subject_id") %>%
  dplyr::select(., dplyr::any_of(c("subject_id", "fecal_calprotectin", "crp_bd1"))) %>%
  merge(., fecal_scfas, by = "subject_id", all = T) %>%
  merge(., scfa_plasma_dedup, by = "subject_id", all = T) %>%
  dplyr::mutate(., crp_bd1_cutoff = ifelse(crp_bd1 >= 10000, "high", "low")) %>%
  dplyr::mutate(., fecal_cal_cutoff = ifelse(fecal_calprotectin >= 100, "high", "low")) %>%
  reshape2::melt(., id.vars = c("subject_id", "fecal_cal_cutoff", "crp_bd1_cutoff")) %>%
  dplyr::filter(., variable %in% c("acetate","propionate", 
                                   "new_butyrate", "new_isobutyrate", 
                                   "total_scfa", "acetate_norm_ratio_dist", 
                                   "propionate_norm_ratio_dist", 
                                   "new_butyrate_norm_ratio_dist", "p_acetic_acid_nmol",
                                   "p_propionic_acid_nmol", "p_butyric_acid_nmol", 
                                   "p_scfa_nmol_total")) 

inflammation_markers$fecal_cal_cutoff <- factor(x = inflammation_markers$fecal_cal_cutoff, levels = c("low", "high"), ordered = T)
inflammation_markers$crp_bd1_cutoff <- factor(x = inflammation_markers$crp_bd1_cutoff, levels = c("low", "high"), ordered = T)
inflammation_markers$variable <- factor(x = inflammation_markers$variable, levels = c("acetate", 
                                                                                      "acetate_norm_ratio_dist", "propionate", 
                                                                                      "propionate_norm_ratio_dist", "new_butyrate", 
                                                                                      "new_butyrate_norm_ratio_dist", "new_isobutyrate", 
                                                                                      "total_scfa", "p_acetic_acid_nmol", "p_propionic_acid_nmol", "p_butyric_acid_nmol", "p_scfa_nmol_total"), ordered = T)
inflammation_markers <- inflammation_markers %>%
  dplyr::mutate(overall_scfa = case_when(
    variable == "acetate" ~ 'fecal acetate',
    variable == "acetate_norm_ratio_dist" ~ 'fecal acetate',
    variable == "propionate" ~ 'fecal propionate',
    variable == "propionate_norm_ratio_dist" ~ 'fecal propionate',
    variable == "new_butyrate" ~ 'fecal butyrate',
    variable == "new_butyrate_norm_ratio_dist" ~ 'fecal butyrate',
    variable == "new_isobutyrate" ~ 'fecal isobutyrate',
    variable == "total_scfa" ~ 'fecal total SCFAs',
    variable == "p_acetic_acid_nmol" ~ 'plasma acetate',
    variable == "p_propionic_acid_nmol" ~ 'plasma propionate',
    variable == "p_butyric_acid_nmol" ~ 'plasma butyrate',
    variable == "p_scfa_nmol_total" ~ 'plasma total SCFAs',
    TRUE ~ 'what_is_this'      
  )
  )
inflammation_markers <- inflammation_markers %>%
  dplyr::mutate(ratio_raw = case_when(
    variable == "acetate" ~ 'raw',
    variable == "acetate_norm_ratio_dist" ~ 'ratio',
    variable == "propionate" ~ 'raw',
    variable == "propionate_norm_ratio_dist" ~ 'ratio',
    variable == "new_butyrate" ~ 'raw',
    variable == "new_butyrate_norm_ratio_dist" ~ 'ratio',
    variable == "new_isobutyrate" ~ 'raw',
    variable == "total_scfa" ~ 'raw',
    variable == "p_acetic_acid_nmol" ~ 'raw',
    variable == "p_propionic_acid_nmol" ~ 'raw',
    variable == "p_butyric_acid_nmol" ~ 'raw',
    variable == "p_scfa_nmol_total" ~ 'raw',
    TRUE ~ 'what_is_this'      
  )
  )

inflammation_markers$ratio_raw <- factor(x = inflammation_markers$ratio_raw, levels = c("raw", "ratio"), ordered = T)

inflammation_distribution <- merge(readr::read_csv("/home/docker/data/CTSC24532USDAWHNRCNu-GIMarkers7Oct2021_DATA_2021-10-07_1627.csv"), 
                              readr::read_csv("/home/docker/data/CRP_WBC_9102021.csv"), 
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
  geom_boxplot(outlier.alpha = 0) + 
  facet_nested(.~ overall_scfa + ratio_raw, scales = "free", independent = "all") +
  ggpubr::stat_pwc(hjust = 0.5, vjust = 1, p.adjust.by = "panel") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.background = element_blank())

## make supp figure 4C =========================================================
sup_figure4C <- inflammation_markers %>% 
  dplyr::filter(., fecal_cal_cutoff != "NA") %>%
  ggplot() + 
  aes(x = fecal_cal_cutoff, y = value) +
  geom_point(position = position_jitter(width = 0.2), alpha = 0.4) + 
  geom_boxplot(outlier.alpha = 0) + 
  facet_nested(.~ overall_scfa + ratio_raw, scales = "free", independent = "all") +
  ggpubr::stat_pwc(hjust = 0.5, vjust = 1) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.background = element_blank())


design <- "A
B
C
D"

fecal_cal_density + crp_density + sup_figure4B + sup_figure4C + patchwork::plot_layout(design = design) + patchwork::plot_annotation(tag_levels = c("A"))
