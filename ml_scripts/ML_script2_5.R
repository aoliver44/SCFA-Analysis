#!/usr/bin/env Rscript

## ML_script2_5.R to prepare the data for the SCFA ML effort
## Author: Andrew Oliver
## Date: Dec 8, 2023
## assumptions: working in this dir on spitfire: /share/lemaylab/aoliver/SCFA
## to run: singularity run -W `pwd` --bind `pwd`:/home/docker /share/lemaylab/aoliver/software/nutrition_tools.sif bash -c 'Rscript /home/docker/ML_script2_5.R'

## load libs
library(dplyr)
library(readr)

## combine taxaHFE outputs with covariates

setwd("/home/docker/input_data/")
response_vars <- gsub(x=list.files(pattern = "_microbe_taxaHFE.csv"), pattern = "_microbe_taxaHFE.csv", replacement= "")

for (var in response_vars) {
  taxaHFE_file_food <- read.csv(file = paste0("/home/docker/input_data/", var, "_food_taxaHFE.csv"))
  taxaHFE_file_microbe <- read.csv(file = paste0("/home/docker/input_data/", var, "_microbe_taxaHFE.csv"))
  taxaHFE_covariates <- read.csv(file = paste0("/home/docker/raw_data/taxahfe_", var, ".csv"))
  
  if ("stool_consistency_class" %in% colnames(taxaHFE_covariates)) {
    taxaHFE_covariates <- taxaHFE_covariates %>%
      dplyr::mutate(stoolc_soft = ifelse(stool_consistency_class=='soft', 1, 0),
                    stoolc_normal = ifelse(stool_consistency_class=='normal', 1, 0),
                    stoolc_hard = ifelse(stool_consistency_class=='hard', 1, 0)) %>%
      dplyr::select(., -stool_consistency_class)
  }
  
  dietML_food <- merge(taxaHFE_covariates, taxaHFE_file_food, by = "subject_id") %>%
    dplyr::select(., -var)
  dietML_microbe <- merge(taxaHFE_covariates, taxaHFE_file_microbe, by = "subject_id") %>%
    dplyr::select(., -var)
  
  write.csv(file = paste0("/home/docker/input_data/", var, "_food_taxaHFE.csv"), x = dietML_food, quote = F, row.names = F)
  write.csv(file = paste0("/home/docker/input_data/", var, "_microbe_taxaHFE.csv"), x = dietML_microbe, quote = F, row.names = F)
  
}