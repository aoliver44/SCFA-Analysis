## Figure 5

## figure5.R: generate figure 5 of SCFA paper
## Author: Andrew Oliver
## Date: March 18, 2024
## to run: docker run --rm -it -p 8787:8787 -e PASSWORD=yourpasswordhere -v `PWD`:/home/docker aoliver44/scfa_analysis:1.1

## load libraries ==============================================================
library(dplyr)
library(shapviz)
library(ggplot2)
library(data.tree)
library(metacoder)


## set working directory =======================================================
setwd("/home/docker")

## source SCFA data ============================================================
source("/home/docker/scripts/polished_scripts/pre_process_raw_scfas.R")

## wrangle data or source helper functions =====================================

## make plots ==================================================================
##plot5a

butyrate_micro_ml_env <- new.env()
load(file = "/home/docker/combined_ml_results/best_models/new_butyrate_norm_ratio_dist_microbe_taxaHFE/ML_r_workspace.rds", envir = butyrate_micro_ml_env)

sv_full <- shapviz::shapviz(butyrate_micro_ml_env$shap_explainations_full, X = butyrate_micro_ml_env$shap_data_full)

shap_values <- as.data.frame(butyrate_micro_ml_env$sv_full$S)
raw_values <- as.data.frame(butyrate_micro_ml_env$sv_full$X)
colnames(raw_values) <- paste0(colnames(raw_values), "_raw")

features <- c("g_roseburia", "g_lachnospiraceae_unclassified", "f_eubacteriales_unclassified", "c_firmicutes_unclassified")

shap_values <- shap_values %>% dplyr::select(., dplyr::any_of(matches(features)))
raw_values <- raw_values %>% dplyr::select(., dplyr::any_of(matches(features)))

all_data <- cbind(shap_values, raw_values)
all_data$index <- seq(1:NROW(all_data))

all_data_melt <- reshape2::melt(all_data, id.vars = "index")
all_data_melt <- all_data_melt %>% 
  dplyr::mutate(., type = ifelse(grepl("raw", all_data_melt$variable), "abundance", "shap"))
all_data_melt$variable <- gsub("_raw", "", x = all_data_melt$variable)
all_data_melt <- all_data_melt %>% tidyr::spread(., key = "type", value = "value")

## make shap plot
plot5a <- ggplot(all_data_melt) + aes(x = abundance, y = shap) + 
  geom_point(aes(color = variable), alpha = 0.3) +
  geom_smooth(aes(color = variable)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw(base_size = 14) +
  theme(
    panel.background = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black")
  ) +
  ggsci::scale_color_aaas() +
  labs(x = "Relative Abundance", y = "SHAPley Value") 


## plot5b
metaphlan <- read.delim("/home/docker/data/merged_metaphlan_v4-0-6.txt", check.names = F)
features_tree <- "g__Roseburia|g__Lachnospiraceae_unclassified|f__Eubacteriales_unclassified|c__Firmicutes_unclassified"
metaphlan <- metaphlan %>% dplyr::filter(., grepl("s__", clade_name)) %>% dplyr::filter(., !grepl("t__", clade_name))
metaphlan <- metaphlan %>% dplyr::filter(., grepl(pattern = features_tree, clade_name))

setwd("/home/docker")
source("/home/docker/scripts/polished_scripts/pre_process_raw_scfas.R")
fecal_scfas <- fecal_scfas %>% dplyr::filter(., subject_id %in% colnames(metaphlan))
fecal_scfas <- fecal_scfas %>% dplyr::mutate(., tertile = ntile(butyrate_norm_ratio_dist, 3))
fecal_scfas <- fecal_scfas %>% dplyr::filter(., tertile %in% c(1,3)) %>%
  dplyr::mutate(., tertile = ifelse(tertile == 1, "low", "high"))
metaphlan <- metaphlan %>% dplyr::select(., clade_name, dplyr::any_of(fecal_scfas$subject_id))

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
plot5b <- obj %>% 
  metacoder::heat_tree(., 
                       node_label = taxon_names,
                       node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
                       node_color = log2_median_ratio, # A column from `obj$data$diff_table`
                       node_color_interval = c(-2, 2), # The range of `log2_median_ratio` to display
                       node_color_range = c("cyan", "gray", "tan"), # The color palette used
                       node_size_axis_label = "OTU count",
                       node_label_size_range = c(0.01, 0.03),
                       node_color_axis_label = "Log 2 ratio of median proportions",
                       layout = "davidson-harel", # The primary layout algorithm
                       initial_layout = "reingold-tilford") # The layout algorithm that initializes node locations

## plot5c
butyrate_humann_ml_env <- new.env()
load(file = "/home/docker/combined_ml_results/best_models/humann_pwys_fecal_new_butyrate_norm_ratio_dist/ML_r_workspace.rds", envir = butyrate_humann_ml_env)

sv_full <- shapviz::shapviz(butyrate_humann_ml_env$shap_explainations_full, X = butyrate_humann_ml_env$shap_data_full)

shap_values <- as.data.frame(butyrate_humann_ml_env$sv_full$S)
raw_values <- as.data.frame(butyrate_humann_ml_env$sv_full$X)
colnames(raw_values) <- paste0(colnames(raw_values), "_raw")

features <- c("thisynara_pwy_superpathway_of_thiamine_diphosphate_biosynthesis_iii_eukaryotes", 
              "pwy_7237_myo_chiro_and_scyllo_inositol_degradation", 
              "pwy_6527_stachyose_degradation")

shap_values <- shap_values %>% dplyr::select(., dplyr::any_of(matches(features)))
raw_values <- raw_values %>% dplyr::select(., dplyr::any_of(matches(features)))

all_data <- cbind(shap_values, raw_values)
all_data$index <- seq(1:NROW(all_data))

all_data_melt <- reshape2::melt(all_data, id.vars = "index")
all_data_melt <- all_data_melt %>% 
  dplyr::mutate(., type = ifelse(grepl("raw", all_data_melt$variable), "abundance", "shap"))
all_data_melt$variable <- gsub("_raw", "", x = all_data_melt$variable)
all_data_melt <- all_data_melt %>% tidyr::spread(., key = "type", value = "value")

## make shap plot
plot5c <- ggplot(all_data_melt) + aes(x = abundance, y = shap) + 
  geom_point(aes(color = variable), alpha = 0.3) +
  geom_smooth(aes(color = variable)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_bw(base_size = 14) +
  theme(
    panel.background = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(color = "black")
  ) +
  ggsci::scale_color_aaas() +
  labs(x = "Relative Abundance", y = "SHAPley Value") 

## make plot 5d
## load in humann data
humann <- readr::read_delim("/home/docker/data/merged_pathabundance-cpm.tsv") %>%
  dplyr::rename(., "pathway" = 1)
source("/home/docker/scripts/polished_scripts/pre_process_raw_scfas.R")
fecal_scfas <- fecal_scfas %>% dplyr::mutate(., tertile = ntile(butyrate_norm_ratio_dist, 3))
fecal_scfas <- fecal_scfas %>% dplyr::filter(., tertile %in% c(1,3)) %>%
  dplyr::mutate(., tertile = ifelse(tertile == 1, "low", "high"))

humann <- humann %>% 
  dplyr::filter(., grepl("THISYNARA-PWY", pathway)) %>%
  dplyr::filter(., grepl("\\|", pathway)) %>%
  tidyr::separate(., col = "pathway", into = c("pathway", "species"), sep = "\\.s__", remove = T)  %>%
  dplyr::select(., -pathway)
colnames(humann) <- gsub(pattern = ".extendedFrags_Abundance", replacement = "", x = colnames(humann))
humann$species <- tidyr::replace_na(humann$species, "unclassified")
humann_melt <- reshape2::melt(humann, id.vars = c("species"))

plot_features <- humann_melt %>% 
  group_by(., species) %>% 
  summarise(., abund = mean(value)) %>% 
  arrange(desc(abund)) %>% 
  slice_head(., n = 9) %>% 
  pull(., species)

humann_melt <- humann_melt %>% 
  mutate(., plot_feature = ifelse(species %in% plot_features, species, "zlow_abundant"))

humann_melt_meta <- merge(humann_melt, fecal_scfas, by.x = "variable", by.y = "subject_id")

## see what taxa are changing (percent wise) with the thiamine gene abundance
#humann_melt_meta %>% dplyr::group_by(., plot_feature, tertile) %>% 
#  dplyr::summarise(., mean_abund = mean(value)) %>%
#  tidyr::spread(., tertile, mean_abund) %>% 
#  dplyr::mutate(., pct_increase = ((high-low) / low )*100) 

humann_melt_meta$tertile <- factor(x = humann_melt_meta$tertile, levels = c("low", "high"), ordered = T)
plot5d <- ggplot(humann_melt_meta) + aes(x = tertile, y = value) +
  stat_summary(fun=mean, geom="bar", position = "stack", aes(fill = plot_feature)) +
  #geom_bar(aes(fill = plot_feature)) + 
  facet_wrap(.~ tertile, scales = "free_x") +
  ggsci::scale_fill_d3() +
  theme_bw() +
  theme(panel.background = element_blank(),
        axis.text = element_text(colour = "black"),
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()) +
  labs(x = "", y = "Copies per million")

design <- "
AB
CD"

plot5a + plot5b + plot5c + plot5d + patchwork::plot_layout(design = design) + patchwork::plot_annotation(tag_levels = "A")
