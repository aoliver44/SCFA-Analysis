## Table 1

## Table1.R: generate table 1 of SCFA paper
## Author: Andrew Oliver
## Date: July 16, 2024
## to run: docker run --rm -it -p 8787:8787 -e PASSWORD=yourpasswordhere -v `PWD`:/home/docker aoliver44/scfa_analysis:1.1

## load libraries ==============================================================
library(patchwork)
library(ggpubr)
library(dplyr)

## set working directory =======================================================
setwd("/home/docker")

## source data =================================================================
source("/home/docker/github/SCFA-Analysis/figure_scripts/pre_process_raw_scfas.R")

inflammation <- merge(readr::read_csv("/home/docker/data/CTSC24532USDAWHNRCNu-GIMarkers7Oct2021_DATA_2021-10-07_1627.csv"), 
                      readr::read_csv("/home/docker/data/CRP_WBC_9102021.csv"), 
                      by = "subject_id")

## source helper functions =====================================================
source("/home/docker/github/SCFA-Analysis/figure_scripts/partial_regression.R")

## wrangle data ================================================================

## add in covariates and scfas to inflammation variables
inflammation <- inflammation %>%
  merge(., fecal_scfas, by = "subject_id", all = T) %>%
  merge(., scfa_plasma_dedup, by = "subject_id", all = T) %>%
  merge(., anthropometrics, by = "subject_id", all = T) %>%
  merge(., fecal_vars %>% dplyr::select(., subject_id, st_wt, StoolConsistencyClass, bristol_num), by = "subject_id", all = T)

## create empty dataframe to store results
all_samples_inflammation <- data.frame(scfa=factor(), factor=factor(), estimate=numeric(),
                                       p_value=numeric(),p.adjust=numeric(), tobit_estimate=numeric(), 
                                       tobit_pvalue=numeric())

## Run partial correlations on inflammation variables
## only run subclincal (not frank inflammation for CRP or fecal calprotectin)
for (inflammation_var in c("fecal_calprotectin", "fecal_mpo", "fecal_neopterin", "plasma_lbp_bd1", "crp_bd1", "wbc_bd1")) {
  ## full dataset
  tmp_fecal <- PartialCorrelationNew(scfas = "fecal", independent = inflammation_var, df = subset(inflammation, inflammation$crp_bd1 < 10000 & inflammation$fecal_calprotectin < 100))
  tmp_serum <- PartialCorrelationNew(scfas = "serum", independent = inflammation_var, df = subset(inflammation, inflammation$crp_bd1 < 10000 & inflammation$fecal_calprotectin < 100))
  all_samples_inflammation <- rbind(all_samples_inflammation, tmp_fecal, tmp_serum)
}

## create varaibales for fecal and plasma SCFAs
all_samples_inflammation <- all_samples_inflammation %>% mutate(., scfa_type = ifelse(scfa %in% c("p_acetic_acid_nmol","p_propionic_acid_nmol", "p_butyric_acid_nmol", "p_scfa_nmol_total"), "plasma", "fecal"))

## p.adjust (not used because these are directed hypotheses)
all_samples_inflammation <- all_samples_inflammation  %>%
  dplyr::group_by(scfa) %>%
  dplyr::mutate(., regression_p.adjust = ifelse(scfa_type == "fecal" | scfa == "p_acetic_acid_nmol" | scfa == "p_scfa_nmol_total", p.adjust(lm_p_value, method='fdr'), p.adjust(tobit_pvalue, method='fdr'))) 

## Reduce to just what is significant & make final table
all_samples_inflammation_filtered <- all_samples_inflammation %>%
  dplyr::select(., scfa_type, scfa, factor, cor_estimate, cor_p_value, lm_p_value, n_ind, tobit_estimate, tobit_pvalue) %>%
  dplyr::filter(., cor_p_value < 0.05) %>%
  dplyr::mutate_if(., is.numeric, ~round(.x, digits = 3)) %>%
  dplyr::rename(., "Type" = "scfa_type",
                "SCFA" = "scfa",
                "Correlation Estimate" = "cor_estimate",
                "Correlation p-value" = "cor_p_value",
                "Regression p-value" = "lm_p_value",
                "N Individuals" = "n_ind",
                "Tobit Estimate" = "tobit_estimate", 
                "Tobit p-value" = "tobit_pvalue") %>% 
  dplyr::arrange(Type)

knitr::kable(all_samples_inflammation_filtered)
