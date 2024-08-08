## Supplemental Figure 8

## supplemental_figure7.R: generate supplemental figure 8 of SCFA paper
## Author: Andrew Oliver
## Date: Feb 29, 2024
## docker run --rm -it \
## -v ~/Downloads/SCFA-Analysis/figure_scripts:/home/scripts \
## -v ~/Downloads/SCFA-Analysis-DATA/data/:/home/data \
## -w /home/docker \
## scfa_analysis:rstudio bash -c "Rscript supplemental_figure8.R"

## load libraries ==============================================================
library(dplyr)
library(ggplot2)
library(wesanderson)
library(readr)

## set working directory =======================================================
setwd("/home/")

## source SCFA data ============================================================
set.seed(123)
source("/home/scripts/pre_process_raw_scfas.R")
source("/home/scripts/partial_regression.R")

## load up other data ==========================================================
## load up l2_processed_cheese
acetate_food_taxaHFE_env <- new.env()
load("/home/data/acetate_food_taxaHFE/ML_r_workspace.rds", envir = acetate_food_taxaHFE_env)
taxaHFE_results <- acetate_food_taxaHFE_env$input %>% dplyr::select(., subject_id, l3_processed_cheeses_and_cheese_spreads)

## merge with SCFA data
taxaHFE_results_scfa <- merge(fecal_scfas, taxaHFE_results) %>%
  merge(read.csv("/home/data/FL100_age_sex_bmi.csv"), by = "subject_id")
stool_vars <- readr::read_delim("/home/data/FL100_stool_variables.txt") %>%
  dplyr::select(., subject_id, st_wt, fecal_calprotectin, StoolConsistencyClass, bristol_num)
taxaHFE_results_scfa <- merge(taxaHFE_results_scfa, stool_vars, by = "subject_id")


## plot figures ================================================================

for (scfa in c("acetate", "propionate",
               "butyrate")) {
  
  set.seed(123)
  taxaHFE_results_scfa$normalized <- bestNormalize::bestNormalize(taxaHFE_results_scfa[[scfa]], allow_orderNorm = F, k = 10, r = 10)$x.t
  print(bestNormalize::bestNormalize(taxaHFE_results_scfa[[scfa]], allow_orderNorm = F, k = 10, r = 10)$chosen_transform)
  
  model <- lm(normalized ~ as.numeric(l3_processed_cheeses_and_cheese_spreads) + as.numeric(bmi) + as.numeric(age) + as.factor(sex.factor) + as.numeric(bristol_num) + as.numeric(st_wt), data = taxaHFE_results_scfa)
  partial_regression <- car::avPlots(model)
  
  tmp_fecal_ph <- cor.test(partial_regression$`as.numeric(l3_processed_cheeses_and_cheese_spreads)`[,1], partial_regression$`as.numeric(l3_processed_cheeses_and_cheese_spreads)`[,2])
  tmp_annotation_data <- PartialCorrelationNew(scfas = scfa, independent = "l3_processed_cheeses_and_cheese_spreads", df = taxaHFE_results_scfa)
  
  plot <- ggplot(as.data.frame(partial_regression$`as.numeric(l3_processed_cheeses_and_cheese_spreads)`)) + 
    aes(x = `as.numeric(l3_processed_cheeses_and_cheese_spreads)`, y = normalized) + 
    geom_point(alpha = 0.1) + 
    geom_smooth(method = "lm", color = "red", se = F, linetype = "dashed") + 
    theme_bw() + 
    theme(panel.grid.minor = element_blank(), axis.text = element_text(colour = "black")) +
    annotate("text", x = 15, y = -2, label= paste0("r = ", round(tmp_annotation_data$cor_estimate, 3), "\np = ", round(tmp_annotation_data$cor_p_value,3))) +
    labs(x = "(L3) Processed cheeses & \ncheese spreads (normalized) | covariates", y = paste0(scfa, " (normalized)\n| covariates"))
  
  assign(x = paste0(scfa, "_cheese"), value = plot)
  count = count + 1
}

design = "ABC"

acetate_cheese + propionate_cheese + butyrate_cheese + patchwork::plot_layout(design = design) + patchwork::plot_annotation(tag_levels = "A")

dir.create(path = "/home/scripts/output_figures", showWarnings = TRUE)
ggsave(filename = "/home/scripts/output_figures/supplemental_figure8.png", 
       device = "png",
       width = 11, 
       height = 4,
       units = "in",
       dpi = 400)
dev.off()
