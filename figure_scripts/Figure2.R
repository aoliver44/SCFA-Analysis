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
path_lf <- "/home/docker/combined_ml_results/combined_results_lf.csv"
combined_results <- readr::read_csv(path_lf) %>% janitor::clean_names()
## make new column based on if taxaHFE, taxaHFE-ml or neither were run
combined_results <- combined_results %>%
  dplyr::mutate(., program = dplyr::case_when( 
    grepl("_lf_", dataset) ~ "taxahfe-ml",
    grepl("taxaHFE", dataset) ~ "taxahfe",
    .default = "non-taxahfe"))

## make an old name column before you make changes below, for the best seeds
combined_results$old_name <- combined_results$dataset

## taxahfe-ml dataset names need to match the taxahfe dataset names
combined_results$dataset <- gsub(pattern = "taxahfe_acetate_lf_food", replacement = "acetate_food_taxaHFE", x = combined_results$dataset)
combined_results$dataset <- gsub(pattern = "taxahfe_acetate_lf_microbe", replacement = "acetate_microbe_taxaHFE", x = combined_results$dataset)
combined_results$dataset <- gsub(pattern = "taxahfe_acetate_lf_norm_ratio_dist_food", replacement = "acetate_norm_ratio_dist_food_taxaHFE", x = combined_results$dataset)
combined_results$dataset <- gsub(pattern = "taxahfe_acetate_lf_norm_ratio_dist_microbe", replacement = "acetate_norm_ratio_dist_microbe_taxaHFE", x = combined_results$dataset)
combined_results$dataset <- gsub(pattern = "taxahfe_new_butyrate_lf_food", replacement = "new_butyrate_food_taxaHFE", x = combined_results$dataset)
combined_results$dataset <- gsub(pattern = "taxahfe_new_butyrate_lf_microbe", replacement = "new_butyrate_microbe_taxaHFE", x = combined_results$dataset)
combined_results$dataset <- gsub(pattern = "taxahfe_new_butyrate_lf_norm_ratio_dist_food", replacement = "new_butyrate_norm_ratio_dist_food_taxaHFE", x = combined_results$dataset)
combined_results$dataset <- gsub(pattern = "taxahfe_new_butyrate_lf_norm_ratio_dist_microbe", replacement = "new_butyrate_norm_ratio_dist_microbe_taxaHFE", x = combined_results$dataset)
combined_results$dataset <- gsub(pattern = "taxahfe_p_acetic_acid_lf_nmol_food", replacement = "p_acetic_acid_nmol_food_taxaHFE", x = combined_results$dataset)
combined_results$dataset <- gsub(pattern = "taxahfe_p_acetic_acid_lf_nmol_microbe", replacement = "p_acetic_acid_nmol_microbe_taxaHFE", x = combined_results$dataset)
combined_results$dataset <- gsub(pattern = "taxahfe_p_butyric_acid_lf_nmol_food", replacement = "p_butyric_acid_nmol_food_taxaHFE", x = combined_results$dataset)
combined_results$dataset <- gsub(pattern = "taxahfe_p_butyric_acid_lf_nmol_microbe", replacement = "p_butyric_acid_nmol_microbe_taxaHFE", x = combined_results$dataset)
combined_results$dataset <- gsub(pattern = "taxahfe_p_propionic_acid_lf_nmol_food", replacement = "p_propionic_acid_nmol_food_taxaHFE", x = combined_results$dataset)
combined_results$dataset <- gsub(pattern = "taxahfe_p_propionic_acid_lf_nmol_microbe", replacement = "p_propionic_acid_nmol_microbe_taxaHFE", x = combined_results$dataset)
combined_results$dataset <- gsub(pattern = "taxahfe_propionate_lf_food", replacement = "propionate_food_taxaHFE", x = combined_results$dataset)
combined_results$dataset <- gsub(pattern = "taxahfe_propionate_lf_microbe", replacement = "propionate_microbe_taxaHFE", x = combined_results$dataset)
combined_results$dataset <- gsub(pattern = "taxahfe_propionate_lf_norm_ratio_dist_food", replacement = "propionate_norm_ratio_dist_food_taxaHFE", x = combined_results$dataset)
combined_results$dataset <- gsub(pattern = "taxahfe_propionate_lf_norm_ratio_dist_microbe", replacement = "propionate_norm_ratio_dist_microbe_taxaHFE", x = combined_results$dataset)
combined_results$dataset <- gsub(pattern = "taxahfe_p_scfa_lf_nmol_total_food", replacement = "p_scfa_nmol_total_food_taxaHFE", x = combined_results$dataset)
combined_results$dataset <- gsub(pattern = "taxahfe_p_scfa_lf_nmol_total_microbe", replacement = "p_scfa_nmol_total_microbe_taxaHFE", x = combined_results$dataset)
combined_results$dataset <- gsub(pattern = "taxahfe_total_scfa_lf_food", replacement = "total_scfa_food_taxaHFE", x = combined_results$dataset)
combined_results$dataset <- gsub(pattern = "taxahfe_total_scfa_lf_microbe", replacement = "total_scfa_microbe_taxaHFE", x = combined_results$dataset)


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

## clean up names
combined_results$scfa_type <- gsub(pattern = "serum", replacement = "plasma", x = combined_results$scfa_type)
combined_results$data_group <- gsub(pattern = "^asa", replacement = "ASA24 Dietary Recalls", x = combined_results$data_group)
combined_results$data_group <- gsub(pattern = "^ffq", replacement = "Food Frequency Questionnaire", x = combined_results$data_group)
combined_results$data_group <- gsub(pattern = "food_taxaHFE", replacement = "Dietary Taxonomy", x = combined_results$data_group)
combined_results$data_group <- gsub(pattern = "inflammation", replacement = "Inflammation", x = combined_results$data_group)
combined_results$data_group <- gsub(pattern = "humann_pwys", replacement = "HUMAnN3 Pathways", x = combined_results$data_group)
combined_results$data_group <- gsub(pattern = "ifdp", replacement = "Inferred Fiber Degradation Profile", x = combined_results$data_group)
combined_results$data_group <- gsub(pattern = "microbe_taxaHFE", replacement = "Microbial Taxonomy", x = combined_results$data_group)
combined_results$data_group <- gsub(pattern = "monosacc", replacement = "Dietary Monosaccharides", x = combined_results$data_group)
combined_results$overall_type <- gsub(pattern = "diet", replacement = "Diet", x = combined_results$overall_type)
combined_results$overall_type <- gsub(pattern = "microbial", replacement = "Microbial", x = combined_results$overall_type)
combined_results$overall_type <- gsub(pattern = "inflammation", replacement = "Inflammation", x = combined_results$overall_type)
combined_results$scfa <- gsub(pattern = "acetate_norm_ratio_dist", replacement = "acetate-ratio", x = combined_results$scfa)
combined_results$scfa <- gsub(pattern = "butyrate_norm_ratio_dist", replacement = "butyrate-ratio", x = combined_results$scfa)
combined_results$scfa <- gsub(pattern = "propionate_norm_ratio_dist", replacement = "propionate-ratio", x = combined_results$scfa)
combined_results$scfa <- gsub(pattern = "total_scfa", replacement = "total SCFAs", x = combined_results$scfa)
combined_results$scfa <- gsub(pattern = "p_acetic_acid_nmol", replacement = "plasma acetate", x = combined_results$scfa)
combined_results$scfa <- gsub(pattern = "p_butyric_acid_nmol", replacement = "plasma butyrate", x = combined_results$scfa)
combined_results$scfa <- gsub(pattern = "p_propionic_acid_nmol", replacement = "plasma propionate", x = combined_results$scfa)
combined_results$scfa <- gsub(pattern = "p_scfa_nmol_total", replacement = "plasma total SCFAs", x = combined_results$scfa)

## calculate mean percent change for MAE
combined_results_summary_mae <- combined_results %>% 
  dplyr::group_by(., dataset, metric, scfa, data_group, scfa_type, overall_type, program) %>% 
  dplyr::summarise(., mean_percent_change = mean(percent_change), mean_original_score = mean(estimate)) %>%
  dplyr::filter(., metric == "mae") %>%
  dplyr::filter(., !is.na(program)) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(., data_group, scfa_type, program) %>%
  dplyr::mutate(., facet_mean = mean(mean_percent_change)) 


## make plot 2 ================================================================

## clean up order
combined_results_summary_mae$data_group <- factor(x = combined_results_summary_mae$data_group, 
                                                  levels = c("ASA24 Dietary Recalls", "Food Frequency Questionnaire", 
                                                             "Dietary Taxonomy", "Dietary Monosaccharides", "Inflammation",
                                                             "HUMAnN3 Pathways", "Inferred Fiber Degradation Profile", "Microbial Taxonomy"),
                                                  ordered = TRUE)

## add colors
combined_results_summary_mae <- combined_results_summary_mae %>%
  dplyr::mutate(., Color = dplyr::case_when( 
                   grepl("acetate", scfa) ~ "#0A9F9D",
                   grepl("butyrate", scfa) ~ "#E54E21",
                   grepl("propionate", scfa) ~ "#CEB175",
                   scfa == "total SCFAs" ~ "#634F25",
                   scfa == "plasma total SCFAs" ~ "#D71515",
                   .default = "other"))

## plot mae difference from null model!
ggplot() + 
  geom_bar(aes(fill = Color, x = scfa, weight = mean_percent_change), data = subset(combined_results_summary_mae, combined_results_summary_mae$program != "taxahfe-ml")) + 
  geom_errorbar(aes(y = facet_mean, ymin = facet_mean, ymax = facet_mean, x = scfa), data = combined_results_summary_mae, linetype = "solid") + 
  facet_nested(.~ overall_type + data_group + scfa_type, scales = "free_x", independent = "x") +
  #facet_grid(. ~ data_group + scfa_type, scales = "free_x") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1), 
      panel.background = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()) + 
  labs(x = "", y = "mean MAE percent change from Null Model\n(lower is better)") +
  scale_fill_identity() +
  scale_linetype_manual("solid")

## get a list of best seeds:
best_seeds <- combined_results %>% 
  dplyr::filter(., metric == "mae") %>%
  dplyr::group_by(., old_name) %>%
  dplyr::arrange(., estimate) %>%
  dplyr::slice(.,1) %>%
  dplyr::select(., seed, old_name)

## check and make sure you have 10 obs per dataset
combined_results %>% dplyr::filter(., metric == "mae") %>% group_by(., dataset, scfa, overall_type, program) %>% tally() %>% View()

## get coeffecient of variation for the MAE of the models
cv <- function(x) 100*( sd(x)/mean(x))
combined_results %>% filter(., metric == "mae", program != "taxahfe-ml") %>%
  group_by(., dataset)%>%
  summarise(cv= cv(estimate)) %>% summarize(., mean_cv = mean(cv), sd_cv =sd(cv))

## null models CV
combined_results %>% filter(., metric == "mae", program != "taxahfe-ml") %>%
  group_by(., dataset)%>%
  summarise(cv= cv(null_model_avg)) %>% summarize(., mean_cv = mean(cv), sd_cv =sd(cv))

## write raw ml results to file
write.csv(combined_results, file = "/home/docker/plots/ml_results_for_publishing.csv", quote = F, row.names = F)

