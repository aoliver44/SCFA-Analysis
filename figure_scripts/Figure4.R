## Figure 4

## figure4.R: generate figure 4 of SCFA paper
## Author: Andrew Oliver
## Date: March 18, 2024
## to run: docker run --rm -it -p 8787:8787 -e PASSWORD=yourpasswordhere -v `PWD`:/home/docker aoliver44/scfa_analysis:1.1

## load libraries ==============================================================
library(dplyr)
library(readr)
library(ggplot2)
library(metacoder)
library(car)
library(bestNormalize)

## set working directory =======================================================
setwd("/home/docker")

## source SCFA data ============================================================
source("/home/docker/github/SCFA-Analysis/figure_scripts/pre_process_raw_scfas.R")
stool_vars <- readr::read_delim("/home/docker/data/FL100_stool_variables.txt") %>%
  dplyr::select(., subject_id, st_wt, fecal_calprotectin, StoolConsistencyClass, bristol_num)

## wrangle data or source helper functions =====================================

# grab the samples you want based on the shap plot (highest 50 lowest 50 shap values)
# this is to compare the L2 potato in the metacoder plot
butyrate_ml_env <- new.env()
load(file = "/home/docker/combined_ml_results/best_models/new_butyrate_food_taxaHFE/ML_r_workspace.rds", envir = butyrate_ml_env)
subject_id_df <- as.data.frame(butyrate_ml_env$input$subject_id)
colnames(subject_id_df)[1] <- "subject_id"
shap_values <- cbind(subject_id_df, butyrate_ml_env$sv_full$S)
shap_values <- shap_values %>% select(., subject_id, l2_white_potatoes_and_puerto_rican_starchy_vegetables)
shap_values <- cbind(shap_values, select(butyrate_ml_env$sv_full$X, l2_white_potatoes_and_puerto_rican_starchy_vegetables) %>% rename(., "abundance" = 1))
high_shap <- shap_values %>% arrange(., desc(l2_white_potatoes_and_puerto_rican_starchy_vegetables)) %>% slice_head(n = 50) %>% pull(subject_id)
low_shap <- shap_values %>% arrange(., desc(l2_white_potatoes_and_puerto_rican_starchy_vegetables)) %>% slice_tail(n = 50) %>% pull(subject_id)
new_samples <- c(high_shap, low_shap)

## load in food tree and format for metacoder
## loading in food tree that is seen by taxaHFE - fills out tree
metaphlan <- read.delim("/home/docker/combined_ml_results/for_potatoes/new_butyrate_food_taxaHFE_raw_data.tsv", check.names = F)
features_tree <- "L2_White_potatoes_and_Puerto_Rican_starchy_vegetables"
metaphlan <- metaphlan %>% dplyr::filter(., grepl(pattern = features_tree, pathString))
metaphlan <- metaphlan %>% dplyr::filter(., depth == 7) %>% dplyr::select(., pathString, 11:305)
metaphlan$pathString <- gsub(pattern = "taxaTree/", replacement = "", x = metaphlan$pathString)
metaphlan$pathString <- gsub(pattern = "L4_/", replacement = "L4_unclassified/", x = metaphlan$pathString)
metaphlan$pathString <- gsub(pattern = "L5_/", replacement = "L5_unclassified/", x = metaphlan$pathString)

metaphlan <- metaphlan %>%
  tidyr::separate(., col = pathString, into = c("kingdom", "phylum", "class", "order", "family", "genus"),
                  sep = "\\/", extra = "merge", remove = T)

metaphlan$kingdom <- gsub(pattern = "L1_", replacement = "L1__", x = metaphlan$kingdom)
metaphlan$phylum <- gsub(pattern = "L2_", replacement = "L2__", x = metaphlan$phylum)
metaphlan$class <- gsub(pattern = "L3_", replacement = "L3__", x = metaphlan$class)
metaphlan$order <- gsub(pattern = "L4_", replacement = "L4__", x = metaphlan$order)
metaphlan$family <- gsub(pattern = "L5_", replacement = "L5__", x = metaphlan$family)
metaphlan$genus <- paste0("L6__", metaphlan$genus)

metaphlan$clade_name <- paste0(metaphlan$kingdom, "|", metaphlan$phylum, "|", metaphlan$class, "|", metaphlan$order, "|", metaphlan$family,"|", metaphlan$genus)

metaphlan <- metaphlan %>% dplyr::select(., clade_name, 7:301)

## load in butyrate data
fecal_scfas <- fecal_scfas %>% dplyr::filter(., subject_id %in% new_samples)
fecal_scfas <- fecal_scfas %>% dplyr::mutate(., tertile = ntile(butyrate, 3))
fecal_scfas <- fecal_scfas %>% dplyr::filter(., tertile %in% c(1,3)) %>%
  dplyr::mutate(., tertile = ifelse(tertile == 1, "low", "high"))
fecal_scfas$tertile <- factor(x = fecal_scfas$tertile, levels = c("low", "high"), ordered = T)
metaphlan <- metaphlan %>% dplyr::select(., clade_name, dplyr::all_of(as.factor(new_samples)))

## do the metacoder
obj <- metacoder::parse_tax_data(metaphlan,
                                 class_cols = "clade_name", # the column that contains taxonomic information
                                 class_sep = "|", # The character used to separate taxa in the classification
                                 class_regex = "^(.+)__(.+)$", # Regex identifying where the data for each taxon is
                                 class_key = c(tax_rank = "info", # A key describing each regex capture group
                                               tax_name = "taxon_name"))

obj$data$tax_abund <- metacoder::calc_taxon_abund(obj, "tax_data", cols = fecal_scfas$subject_id)
obj$data$tax_occ <- metacoder::calc_n_samples(obj, "tax_abund", cols = fecal_scfas$subject_id)
obj$data$type_abund <- metacoder::calc_group_mean(obj, "tax_abund",
                                                  cols = fecal_scfas$subject_id, 
                                                  groups = fecal_scfas$tertile)

obj$data$diff_table <- metacoder::compare_groups(obj,
                                                 data = "tax_abund",
                                                 cols = fecal_scfas$subject_id, # What columns of sample data to use
                                                 groups = fecal_scfas$tertile) # What category each sample is assigned to
set.seed(1) # This makes the plot appear the same each time it is run 
## Treatment 1 = TAN
## Treatment 2 = BLUE
## in the bottom graph, more blue = more associated with High butyrate tertile
plot4a <- obj %>% 
  metacoder::heat_tree(., 
                       node_label = taxon_names,
                       node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
                       node_color = mean_diff, # A column from `obj$data$diff_table`
                       #node_color_interval = c(0.8, 1.5), # The range of `log2_median_ratio` to display
                       node_color_range = c("cyan", "gray", "tan"), # The color palette used
                       node_size_axis_label = "OTU count",
                       node_label_size_range = c(0.01, 0.03),
                       node_color_axis_label = "Mean Abundance Difference\n(blue = more abundant in high tertile)",
                       layout = "davidson-harel", # The primary layout algorithm
                       initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations

#View(merge(select(obj$data$tax_data, taxon_id, clade_name), obj$data$diff_table, by = "taxon_id"))

## food tree ~ stacked barplot
metaphlan <- read.delim("/home/docker/combined_ml_results/for_potatoes/new_butyrate_food_taxaHFE_raw_data.tsv", check.names = F)
features_tree <- "L2_White_potatoes_and_Puerto_Rican_starchy_vegetables"
metaphlan <- metaphlan %>% dplyr::filter(., grepl(pattern = features_tree, pathString))
metaphlan <- metaphlan %>% dplyr::filter(., depth == 7) %>% dplyr::select(., pathString, 11:305)
metaphlan$pathString <- gsub(pattern = "taxaTree/", replacement = "", x = metaphlan$pathString)
metaphlan$pathString <- gsub(pattern = "L4_/", replacement = "L4_unclassified/", x = metaphlan$pathString)
metaphlan$pathString <- gsub(pattern = "L5_/", replacement = "L5_unclassified/", x = metaphlan$pathString)

metaphlan <- metaphlan %>%
  tidyr::separate(., col = pathString, into = c("kingdom", "phylum", "class", "order", "family", "genus"),
                  sep = "\\/", extra = "merge", remove = T)

metaphlan <- metaphlan%>% dplyr::select(., class, 7:301) %>% reshape2::melt(., id.vars = "class") %>% dplyr::rename("class" = 1, "subject_id" = 2)

plot_features <- metaphlan %>% 
  group_by(., class) %>% 
  summarise(., abund = mean(value)) %>% 
  arrange(desc(abund)) %>% 
  slice_head(., n = 9) %>% 
  pull(., class)

metaphlan <- metaphlan %>% 
  mutate(., plot_feature = ifelse(class %in% plot_features, class, "zlow_abundant"))

plot4b <- ggplot(metaphlan) + aes(x = "total_cohort", weight = value) + 
  geom_bar(aes(fill = plot_feature), position = position_fill()) + 
  scale_y_continuous(labels = scales::percent_format()) +
  ggsci::scale_fill_d3(palette = "category20") +
  labs(x = "", y = "Relative Abundance")

## butyrate vs l2 potatoes partial regression
butyrate_shap <- new.env()
load(file = "/home/docker/combined_ml_results/best_models/new_butyrate_food_taxaHFE/ML_r_workspace.rds", envir = butyrate_shap)

butyrate_shap$input <- merge(butyrate_shap$input, dplyr::select(stool_vars, subject_id, bristol_num), by = "subject_id")

set.seed(123)
butyrate_shap$input$normalized <- bestNormalize::bestNormalize(butyrate_shap$input$label, allow_orderNorm = T, k = 10, r = 20)$x.t
print(bestNormalize::bestNormalize(butyrate_shap$input$label, allow_orderNorm = T, k = 10, r = 20)$chosen_transform)

model <- lm(normalized ~ as.numeric(l2_white_potatoes_and_puerto_rican_starchy_vegetables) + as.numeric(bmi) + as.numeric(age) + as.factor(sex_factor) + as.numeric(bristol_num) + as.numeric(st_wt), data = butyrate_shap$input)
partial_regression <- car::avPlots(model)

tmp_fecal_ph <- cor.test(partial_regression$`as.numeric(l2_white_potatoes_and_puerto_rican_starchy_vegetables)`[,1], partial_regression$`as.numeric(l2_white_potatoes_and_puerto_rican_starchy_vegetables)`[,2])

plot4c <- ggplot(as.data.frame(partial_regression$`as.numeric(l2_white_potatoes_and_puerto_rican_starchy_vegetables)`)) + 
  aes(x = `as.numeric(l2_white_potatoes_and_puerto_rican_starchy_vegetables)`, y = normalized) + 
  geom_point(alpha = 0.1) + 
  geom_smooth(method = "lm", color = wesanderson::wes_palette("AsteroidCity1")[3], se = F, linetype = "dashed") + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), axis.text = element_text(colour = "black"), panel.background = element_blank(), 
        panel.border = element_rect(linewidth = 2)) +
  labs(x = "L2 White Potatoes (normalized) | covariates", y = paste0("Fecal butyrate", " (normalized) | covariates"))


design <- "
AAAB
AAAC"

plot4a + plot4b + plot4c + patchwork::plot_layout(design = design) + patchwork::plot_annotation(tag_levels = "A")
