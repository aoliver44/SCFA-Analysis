## Figure 2

## figure2.R: generate figure 2 of SCFA paper
## Author: Andrew Oliver
## Date: March 18, 2024
## to run: docker run --rm -it -p 8787:8787 -e PASSWORD=yourpasswordhere -v `PWD`:/home/docker aoliver44/scfa_analysis:1.1

## load libraries ==============================================================
library(dplyr)
library(readr)
library(ggplot2)
library(wesanderson)
library(patchwork)
library(ggh4x)

## set working directory =======================================================
setwd("/home/docker")

## source data =================================================================
path <- "/home/docker/combined_ml_results/combined_results.csv"
combined_results <- readr::read_csv(path) %>% janitor::clean_names()
combined_results$dataset <- gsub(pattern = "_new_", replacement = "_", combined_results$dataset)
combined_results$dataset <- gsub(pattern = "new_", replacement = "", combined_results$dataset)

## wrangle data or source helper functions =====================================
## scfas
fecal_scfa_list <- c("acetate","propionate", "butyrate", "total_scfa", "acetate_norm_ratio_dist", "propionate_norm_ratio_dist", 
                     "butyrate_norm_ratio_dist")
serum_scfa_list <- c("p_acetic_acid_nmol","p_propionic_acid_nmol", "p_butyric_acid_nmol", 
                     "p_scfa_nmol_total")

all_scfas <- c(fecal_scfa_list, serum_scfa_list)
all_scfas_pattern <- paste(all_scfas, collapse="|")

## calc MAE percent increase over null model
combined_results <- combined_results %>%
  dplyr::mutate(., percent_change = ((estimate - null_model_avg) / (null_model_avg + 0.0000000000000001))) 

combined_results$scfa <- combined_results$dataset
combined_results$scfa <- gsub(x = combined_results$scfa, pattern = c("_food_taxaHFE|_microbe_taxaHFE|asa_fecal_|asa_serum_|ffq_fecal_|ffq_serum_|humann_pwys_fecal_|humann_pwys_serum_|ifdp_fecal_|ifdp_serum_|monosacc_fecal_|monosacc_serum_|inflammation_"), replacement = "")
combined_results$data_group <- combined_results$dataset
combined_results$data_group <- gsub(x = combined_results$data_group, pattern = c(all_scfas_pattern), replacement = "")
combined_results$data_group <- gsub(x = combined_results$data_group, pattern = c("_fecal_|_serum_"), replacement = "")
combined_results$data_group <- gsub(x = combined_results$data_group, pattern = "_food", replacement = "food")
combined_results$data_group <- gsub(x = combined_results$data_group, pattern = "_microbe", replacement = "microbe")
combined_results$data_group <- gsub(x = combined_results$data_group, pattern = "inflammation_", replacement = "inflammation")

combined_results <- combined_results %>%
  dplyr::mutate(., scfa_type = ifelse(scfa %in% fecal_scfa_list, "fecal", "serum"))

combined_results <- combined_results %>%
  dplyr::mutate(., overall_type = ifelse(data_group %in% c("microbe_taxaHFE", "ifdp", "ifdp", "humann_pwys"), "microbial", 
                                         ifelse(data_group %in% c("asa", "ffq", "food_taxaHFE", "monosacc"), "diet", "inflammation")))

## calculate mean percent change for MAE
combined_results_summary_mae <- combined_results %>% 
  dplyr::group_by(., dataset, metric, scfa, data_group, scfa_type, overall_type) %>% 
  dplyr::summarise(., mean_percent_change = mean(percent_change), mean_original_score = mean(estimate)) %>%
  dplyr::filter(., metric == "mae")


## make plot 2 ================================================================

combined_results_summary_mae %>% ggplot() + 
  aes(x = scfa, weight = mean_percent_change) + 
  geom_bar() + 
  facet_nested(.~ overall_type + data_group + scfa_type, scales = "free_x", independent = "x") +
  #facet_grid(. ~ data_group + scfa_type, scales = "free_x") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1), 
        panel.background = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()) + 
  labs(x = "", y = "mean MAE percent change from Null Model\n(lower is better)")
