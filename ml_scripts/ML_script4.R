#!/usr/bin/env Rscript

## ML_script4.R to collect the results of all the dietML runs
## Author: Andrew Oliver
## Date: Nov 15, 2023
## assumptions: working in this dir on spitfire: /share/lemaylab/aoliver/SCFA
## to run: singularity run -W `pwd` --bind `pwd`:/home/docker /share/lemaylab/aoliver/software/nutrition_tools.sif bash -c 'Rscript /home/docker/ML_script4.R'

## load libraries
library(dplyr)
library(readr)

## create directory for all
dir.create("/home/docker/combined_ml_results")

## copy file and rename if from where it came from
new_name <- list.dirs("/home/docker/dietML_results", full.names=F)
new_name <- new_name[2:length(new_name)]

for (ml_result in seq(1:length(list.files(path = "/home/docker/dietML_results/", recursive=T, pattern="ml_results.csv")))) {

  ## read in and add the column of where it came from
  tmp <- suppressMessages(readr::read_csv(file = paste0("/home/docker/dietML_results/", new_name[ml_result], "/ml_results.csv"), col_names = T)) %>%
    dplyr::mutate(., dataset = new_name[ml_result]) %>% 
    janitor::clean_names()
  readr::write_csv(x = tmp, file = paste0("/home/docker/combined_ml_results/", new_name[ml_result], ".csv"))
  
}

## append all results together
setwd("/home/docker/combined_ml_results/")
files <- list.files(pattern = "\\.csv$")
DF <-  read.csv(files[1])
for (f in files[-1]) DF <- rbind(DF, read.csv(f))   
write.csv(DF, "combined_results.csv", row.names=FALSE, quote=FALSE)



