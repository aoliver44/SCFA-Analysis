## Supplemental Figure 7

## supplemental_figure7.R: generate supplemental figure 7 of SCFA paper
## Author: Andrew Oliver
## Date: Feb 29, 2024
## docker run --rm -it \
## -v ~/Downloads/SCFA-Analysis/figure_scripts:/home/scripts \
## -v ~/Downloads/SCFA-Analysis-DATA/data/:/home/data \
## -w /home/docker \
## scfa_analysis:rstudio bash -c "Rscript supplemental_figure7.R"

## load libraries ==============================================================
library(dplyr)
library(ggplot2)
library(vegan)

## set working directory =======================================================
setwd("/home/")

## source SCFA data ============================================================
set.seed(123)
source("/home/scripts/pre_process_raw_scfas.R")
source("/home/scripts/partial_regression.R")
anthropometrics <- readr::read_csv("/home/data/FL100_age_sex_bmi.csv")
stool_vars <- readr::read_delim("/home/data/FL100_stool_variables.txt")

## read in propionate permanova data ===========================================
## to get the model coeffiecents
df <- readr::read_csv(file = "/home/data/permanova_results/propionate_norm_ratio_dist/propionate_norm_ratio_dist_microbe_taxaHFE.csv") %>%
  tibble::column_to_rownames(., var = "subject_id")
metadata <- merge(stool_vars, anthropometrics, by = "subject_id", all = T) %>%
  dplyr::select(., dplyr::any_of(c("subject_id", "st_wt", "StoolConsistencyClass", "sex.factor", "bmi", "age")))

df <- merge(metadata, df, by.x = "subject_id", by.y = "row.names")
df <- df %>% tidyr::drop_na()

permanova_results <- vegan::adonis(formula = df[8:NCOL(df)] ~ as.numeric(bmi) + 
                                     as.numeric(age) + 
                                     as.factor(sex.factor) + 
                                     as.factor(StoolConsistencyClass) + 
                                     as.numeric(st_wt) +
                                     as.numeric(feature_of_interest),
                                   method = "bray", 
                                   permutations = 999, 
                                   data = df, by = "terms")

coef <- coefficients(permanova_results)["as.numeric(feature_of_interest)",]
top.coef <- as.data.frame(coef) %>% tibble::rownames_to_column(., var = "taxa") %>% dplyr::arrange(., desc(coef))

coef_plot <- ggplot(data = top.coef) + 
  aes(x = reorder(taxa, coef), weight = coef) +
  geom_bar() + 
  coord_flip() + 
  labs(x = "Taxa", y = "PERMANOVA Coefficents") + 
  theme_bw(base_size = 14)

## propionate
propionate_shap <- new.env()
load(file = "/home/data/propionate_norm_ratio_dist_microbe_taxaHFE/ML_r_workspace.rds", envir = propionate_shap)

propionate_shap_plot <- shapviz::sv_importance(propionate_shap$sv_full, kind = "bee", show_numbers = TRUE, bee_width = 0.2, max_display = 5) + 
  labs(x = paste0("predicts low propionate-ratio", "< SHAP >","predicts high propionate-ratio")) + 
  theme_bw(base_size = 14) + theme(axis.text = element_text(colour = "black"), 
                                   panel.grid.minor = element_blank(), 
                                   panel.background = element_blank(), 
                                   panel.border = element_rect(linewidth = 2))

design <- "AB"

coef_plot + propionate_shap_plot + patchwork::plot_layout(design = design) + patchwork::plot_annotation(tag_levels = "A")

dir.create(path = "/home/scripts/output_figures", showWarnings = TRUE)
ggsave(filename = "/home/scripts/output_figures/supplemental_figure7.png", 
       device = "png",
       width = 12, 
       height = 6,
       units = "in",
       dpi = 400)
dev.off()
