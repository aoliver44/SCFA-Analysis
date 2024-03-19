## Figure 3

## figure3.R: generate figure 3 of SCFA paper
## Author: Andrew Oliver
## Date: March 18, 2024
## to run: docker run --rm -it -p 8787:8787 -e PASSWORD=yourpasswordhere -v `PWD`:/home/docker aoliver44/scfa_analysis:1.1

## load libraries ==============================================================
library(dplyr)
library(ggplot2)
library(wesanderson)
library(shapviz)

## set working directory =======================================================
setwd("/home/docker")

## source SCFA data ============================================================
source("/home/docker/scripts/polished_scripts/pre_process_raw_scfas.R")
stool_vars <- readr::read_delim("/home/docker/data/FL100_stool_variables.txt") %>%
  dplyr::select(., subject_id, st_wt, fecal_calprotectin, StoolConsistencyClass, bristol_num)

## wrangle data or source helper functions =====================================

## acetate
acetate_shap <- new.env()
load(file = "/home/docker/combined_ml_results/best_models/acetate_food_taxaHFE/ML_r_workspace.rds", envir = acetate_shap)

## propionate
propionate_shap <- new.env()
load(file = "/home/docker/combined_ml_results/best_models/propionate_food_taxaHFE/ML_r_workspace.rds", envir = propionate_shap)

## butyrate
butyrate_shap <- new.env()
load(file = "/home/docker/combined_ml_results/best_models/new_butyrate_food_taxaHFE/ML_r_workspace.rds", envir = butyrate_shap)

## make figure =================================================================

acetate_shap_plot <- shapviz::sv_importance(acetate_shap$sv_full, kind = "bee", show_numbers = TRUE, bee_width = 0.2, max_display = 5) + 
  labs(x = paste0("predicts low acetate", "< SHAP >","predicts high acetate")) + 
  theme_bw(base_size = 14) + theme(axis.text = element_text(colour = "black"), 
                                   panel.grid.minor = element_blank(), 
                                   panel.background = element_blank(), 
                                   panel.border = element_rect(linewidth = 2))

propionate_shap_plot <- shapviz::sv_importance(propionate_shap$sv_full, kind = "bee", show_numbers = TRUE, bee_width = 0.2, max_display = 5) + 
  labs(x = paste0("predicts low propionate", "< SHAP >","predicts high propionate")) + 
  theme_bw(base_size = 14) + theme(axis.text = element_text(colour = "black"), 
                                   panel.grid.minor = element_blank(), 
                                   panel.background = element_blank(), 
                                   panel.border = element_rect(linewidth = 2))

butyrate_shap_plot <- shapviz::sv_importance(butyrate_shap$sv_full, kind = "bee", show_numbers = TRUE, bee_width = 0.2, max_display = 5) + 
  labs(x = paste0("predicts low butyrate", "< SHAP >","predicts high butyrate")) + 
  theme_bw(base_size = 14) + theme(axis.text = element_text(colour = "black"), 
                                   panel.grid.minor = element_blank(), 
                                   panel.background = element_blank(), 
                                   panel.border = element_rect(linewidth = 2))


design <- "A
B
C"

acetate_shap_plot + propionate_shap_plot + butyrate_shap_plot + patchwork::plot_layout(design = design) + patchwork::plot_annotation(tag_levels = "A")
