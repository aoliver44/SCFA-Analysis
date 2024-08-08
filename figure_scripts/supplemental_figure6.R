## Supplemental Figure 6

## supplemental_figure6.R: generate supplemental figure 6 of SCFA paper
## Author: Andrew Oliver
## Date: Feb 29, 2024
## docker run --rm -it \
## -v ~/Downloads/SCFA-Analysis/figure_scripts:/home/scripts \
## -v ~/Downloads/SCFA-Analysis-DATA/data/:/home/data \
## -w /home/docker \
## scfa_analysis:rstudio bash -c "Rscript supplemental_figure6.R"

## load libraries ==============================================================
library(dplyr)
library(picante)
library(ggplot2)
library(wesanderson)
library(vegan)

## set working directory =======================================================
setwd("/home/")

## source SCFA data ============================================================
set.seed(123)
source("/home/scripts/pre_process_raw_scfas.R")
source("/home/scripts/partial_regression.R")
anthropometrics <- readr::read_csv("/home/data/FL100_age_sex_bmi.csv")
stool_vars <- readr::read_delim("/home/data/FL100_stool_variables.txt")

## calculate phylogenetic diversity of microbiome ==============================

## read in tree
tree <- ape::read.tree("/home/data/mpa_vOct22_CHOCOPhlAnSGB_202212.nwk")
community <- readr::read_delim("/home/data/merged_metaphlan_v4-0-6.txt") %>% 
  dplyr::filter(., grepl("t__", clade_name)) %>%
  tidyr::separate(., col = "clade_name", into = c("L1","L2","L3","L4","L5","L6","L7","L8"), sep = "\\|", extra = "merge")
community$L8 <- gsub(pattern = "SGB", replacement = "", x = community$L8)
community$L8 <- gsub(pattern = "t__", replacement = "", x = community$L8)
community$L8 <- gsub(pattern = "_group", replacement = "", x = community$L8)

community <- community %>% select(., -dplyr::any_of(c("L1","L2","L3","L4","L5","L6","L7")))
community <- community %>% tibble::column_to_rownames(., var = "L8") %>% t() %>% as.data.frame()

metaphlan_faith_pd <- picante::pd(community, tree)

metaphlan_faith_pd <- metaphlan_faith_pd %>% tibble::rownames_to_column(., var = "subject_id")

## combine with SCFA data ======================================================

metaphlan_faith_pd_scfa <- merge(anthropometrics, stool_vars, by = "subject_id") %>%
  merge(., metaphlan_faith_pd, by = "subject_id") %>%
  merge(., fecal_scfas, by = "subject_id") %>%
  merge(., scfa_plasma_dedup, by = "subject_id")
## the above merge created an error of sorts. Samples were dropped due to the 
## merge with plasma SCFAs (we didnt use them so the merge was unnecessary). The
## values in the paper are slightly different than if you use all possible fecal
## samples (all lower by about 0.02-0.06). The significance does not change.
## conduct partial correlations ================================================

## fecal SCFAs
PartialCorrelationNew(scfas = "fecal", independent = "PD", df = metaphlan_faith_pd_scfa, remove_outliers = F)
PartialCorrelationNew(scfas = "serum", independent = "PD", df = metaphlan_faith_pd_scfa, remove_outliers = F)

## plot fecal SCFAs ~ PD =======================================================

for (scfa_tmp in c("acetate", "acetate_norm_ratio_dist", "propionate", "propionate_norm_ratio_dist",
               "butyrate", "butyrate_norm_ratio_dist", "total_scfa")) {
  
  set.seed(123)
  metaphlan_faith_pd_scfa$normalized <- bestNormalize::bestNormalize(metaphlan_faith_pd_scfa[[scfa_tmp]], allow_orderNorm = F, k = 10, r = 10)$x.t
  
  model <- lm(normalized ~ as.numeric(PD) + as.numeric(bmi) + as.numeric(age) + as.factor(sex.factor) + as.numeric(bristol_num) + as.numeric(st_wt), data = metaphlan_faith_pd_scfa)
  partial_regression <- car::avPlots(model)
  
  tmp_fecal_ph <- cor.test(partial_regression$`as.numeric(PD)`[,1], partial_regression$`as.numeric(PD)`[,2])
  
  tmp_annotation_data <- PartialCorrelationNew(scfa = scfa_tmp, independent = "PD", df = metaphlan_faith_pd_scfa)
  
  plot <- ggplot(as.data.frame(partial_regression$`as.numeric(PD)`)) + 
    aes(x = `as.numeric(PD)`, y = normalized) + 
    geom_point(alpha = 0.1) + 
    geom_smooth(method = "lm", color = "red", se = F, linetype = "dashed") + 
    theme_bw() + 
    theme(panel.grid.minor = element_blank(), axis.text = element_text(colour = "black")) +
    annotate("text", x = 10, y = -2, label= paste0("r = ", round(tmp_annotation_data$cor_estimate, 4), "\np = ", round(tmp_annotation_data$cor_p_value,4))) +
    labs(x = "Microbiome PD \n(normalized) | covariates", y = paste0(scfa_tmp, " (normalized)\n| covariates"))
  
  assign(x = paste0(scfa_tmp, "_pd"), value = plot)
  count = count + 1
}

## generate permanova figure
metadata <- merge(stool_vars, anthropometrics, by = "subject_id", all = T) %>%
  dplyr::select(., dplyr::any_of(c("subject_id", "st_wt", "StoolConsistencyClass", "sex.factor", "bmi", "age")))

## list of files:
list_of_files <- c(paste0("_microbe_taxaHFE_level_", seq(1:8), ".csv"),
                   "_microbe_taxaHFE.csv",
                   "_microbe_taxaHFE_no_sf.csv")

## make empty dataframe to store results
results_df <- data.frame(scfa = character(),
                         filename=character(), 
                         r_sq=numeric(),
                         p_value=numeric(),
                         n_features=numeric())

## loop over directories
for (directory in list.dirs("/home/data/permanova_results")[2:12]) {
  setwd(directory)
  scfa_name <- print(basename(directory))
  
  ## loop over levels
  for (file in list_of_files) {
    
    df <- readr::read_csv(file = paste0(scfa_name, file)) %>%
      tibble::column_to_rownames(., var = "subject_id")
    
    df <- merge(metadata, df, by.x = "subject_id", by.y = "row.names")
    df <- df %>% tidyr::drop_na()
    
    permanova_results <- vegan::adonis2(formula = df[8:NCOL(df)] ~ as.numeric(bmi) + 
                                          as.numeric(age) + 
                                          as.factor(sex.factor) + 
                                          as.factor(StoolConsistencyClass) + 
                                          as.numeric(st_wt) +
                                          as.numeric(feature_of_interest),
                                        method = "bray", 
                                        permutations = 999, 
                                        data = df, by = "terms")
    
    
    ## append the results to the results df
    results_df <- results_df %>% dplyr::add_row(
      scfa = scfa_name,
      filename = file,
      r_sq = permanova_results$R2[6],
      p_value = permanova_results$`Pr(>F)`[6],
      n_features = (NCOL(df) - 7)
    )
    
  }
}

## bars that are have low alpha are not significant PERMANOVAs
results_df$filename <- gsub(pattern = "_microbe_taxaHFE_level_", replacement = "level_", x = results_df$filename)
results_df$filename <- gsub(pattern = "_microbe_taxaHFE.csv", replacement = "taxaHFE", x = results_df$filename)
results_df$filename <- gsub(pattern = "_microbe_taxaHFE_no_sf.csv", replacement = "taxaHFE_no_sf", x = results_df$filename)
results_df$filename <- gsub(pattern = ".csv", replacement = "", x = results_df$filename)

permanova_r2 <- results_df %>%
  dplyr::mutate(alpha = ifelse(p_value > 0.05, "Not_significant", "Significant")) %>% 
  ggplot() + 
  aes(x = filename, weight = r_sq, alpha = alpha) + 
  geom_bar(aes(fill = filename)) + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(x = "", y = "PERMANOVA R2") +
  scale_alpha_discrete(range = c(0.15, 0.9)) +
  facet_grid(. ~ scfa) +
  ggsci::scale_fill_futurama()

design = "ABCD
EFGG
HHHH"

acetate_pd + acetate_norm_ratio_dist_pd + propionate_pd + propionate_norm_ratio_dist_pd + butyrate_pd + butyrate_norm_ratio_dist_pd + total_scfa_pd + permanova_r2 + patchwork::plot_layout(design = design) + patchwork::plot_annotation(tag_levels = "A")

dir.create(path = "/home/scripts/output_figures", showWarnings = TRUE)
ggsave(filename = "/home/scripts/output_figures/supplemental_figure6.png", 
       device = "png",
       width = 12, 
       height = 12,
       units = "in",
       dpi = 400)
dev.off()
