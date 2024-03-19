#!/usr/bin/env Rscript

## ML_script1.R to prepare the data for the SCFA ML effort
## Author: Andrew Oliver
## Date: Nov 1, 2023
## assumptions: working in this dir on spitfire: /share/lemaylab/aoliver/SCFA
## to run: singularity run -W `pwd` --bind `pwd`:/home/docker /share/lemaylab/aoliver/software/nutrition_tools.sif bash -c 'Rscript /home/docker/ML_script1.R'

## load libs
library(dplyr)
library(readr)

## output dir
outpath <- "/home/docker/input_data/"

## read in data
path = "/home/docker/raw_data/"
#source("/home/docker/scripts/pre_process_raw_scfas.R") generates fecal_scfas.csv and scfa_plasma_dedup.csv
fecal_scfas <- readr::read_csv(paste0(path, "fecal_scfas.csv"))
scfa_plasma_dedup <- readr::read_csv(paste0(path, "scfa_plasma_dedup.csv"))
ffq <- readr::read_csv(paste0(path, "FL100_FFQ_cleaned_all_dt.csv"))
asa24 <- readr::read_csv(paste0(path, "FL100_ASA24_avgs_cleaned_all_dt.csv"))
pd_fat <- readr::read_delim(paste0(path, "qiime1_alphadiv_fat.txt"))
pd_fiber <- readr::read_delim(paste0(path, "qiime1_alphadiv_fiber.txt"))
pd_protein <- readr::read_delim(paste0(path, "qiime1_alphadiv_protein.txt"))
pd_carbs <- readr::read_delim(paste0(path, "qiime1_alphadiv_carb.txt"))
pd_all <- readr::read_delim(paste0(path, "qiime1_alphadiv_filt.txt"))
monosacc <- readr::read_csv(paste0(path, "diet_monosacc_fl100_111422_forAO.csv")) %>%
  dplyr::select(., 1:11)
hei_ffq <- readr::read_csv(paste0(path, "HEI FFQ_scores_12072021.csv")) %>%
  dplyr::select(., subject_id, hei_ffq_totalscore)
hei_asa <- readr::read_delim(paste0(path, "FL100_HEI_n378.txt")) %>%
  dplyr::select(., subject_id, hei_asa24_totalscore)
asa24_fiber <- readr::read_delim(paste0(path, "ASA24_average_fiber_summary_variables.txt"))
ifdp <- readr::read_csv(paste0(path, "IFDP_combined_results.csv")) %>%
  tibble::column_to_rownames(var = "fiber") %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(., var = "subject_id") %>%
  dplyr::mutate(., subject_id = gsub(x = subject_id, pattern = "X", replacement = ""))
humann_pwys <- readr::read_delim(paste0(path, "merged_pathabundance-cpm_nostrat.tsv"))
colnames(humann_pwys) <- gsub(pattern = '.extendedFrags_Abundance', replacement = '', x = colnames(humann_pwys))
humann_pwys <- humann_pwys %>%
  dplyr::rename(., "pwy" = 1) %>%
  tibble::column_to_rownames(., var = "pwy") %>%
  t() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(., var = "subject_id") %>%
  janitor::clean_names()
anthropometrics <- readr::read_csv(paste0(path, "FL100_age_sex_bmi.csv")) %>%
  dplyr::mutate(., sex.factor = ifelse(sex.factor == "Male", 1, 0))
stool_vars <- readr::read_delim(paste0(path, "FL100_stool_variables.txt")) %>%
  dplyr::select(., subject_id, st_wt, StoolConsistencyClass) %>%
  dplyr::mutate(stoolc_soft = ifelse(StoolConsistencyClass=='soft', 1, 0),
                stoolc_normal = ifelse(StoolConsistencyClass=='normal', 1, 0),
                stoolc_hard = ifelse(StoolConsistencyClass=='hard', 1, 0)) %>%
  dplyr::select(., -StoolConsistencyClass)
stool_vars_taxaHFE_cov <- readr::read_delim(paste0(path, "FL100_stool_variables.txt")) %>%
  dplyr::select(., subject_id, st_wt, StoolConsistencyClass)
inflammation_markers <- merge(readr::read_csv(paste0(path, "CTSC24532USDAWHNRCNu-GIMarkers7Oct2021_DATA_2021-10-07_1627.csv")), 
                              readr::read_csv(paste0(path, "CRP_WBC_9102021.csv")), 
                              by = "subject_id") %>%
  dplyr::select(., -dplyr::any_of(c("redcap_survey_identifier", "fecal_ph", "gi_markers_complete")))


## scfas
fecal_scfa_list <- c("acetate","propionate", "new_butyrate", "total_scfa", "acetate_norm_ratio_dist", "propionate_norm_ratio_dist",
                     "new_butyrate_norm_ratio_dist")
serum_scfa_list <- c("p_acetic_acid_nmol","p_propionic_acid_nmol", "p_butyric_acid_nmol",
                     "p_scfa_nmol_total")

## subset the scfa dataframes to the ones in the lists above
fecal_scfas <- fecal_scfas %>% dplyr::select(., subject_id, dplyr::any_of(fecal_scfa_list))
serum_scfas <- scfa_plasma_dedup %>% dplyr::select(., subject_id, dplyr::any_of(serum_scfa_list))

## combine FFQ data and ASA data and Monosacc data:

#######
## FFQ =========================================================================
#######

## FFQ final fecal (plus covariates)
ffq_fecal <- merge(fecal_scfas, serum_scfas, by = "subject_id", all = T) %>%
  merge(., anthropometrics, by = "subject_id", all = T) %>% # covariates
  merge(., stool_vars, by = "subject_id", all = T) %>% # covariates
  merge(., ffq, by = "subject_id", all = T) %>%
  merge(., hei_ffq, by = "subject_id", all = T)

for (scfa in fecal_scfa_list) {
  
  ## pull out SCFAs not being tested
  not_used_scfas <- fecal_scfa_list[fecal_scfa_list != scfa]
  ffq_final_tmp <- ffq_fecal %>%
    dplyr::select(., -dplyr::any_of(not_used_scfas), -dplyr::any_of(serum_scfa_list)) %>%
    tidyr::drop_na() %>%
    janitor::clean_names()
  
  ## write file
  write.csv(x = ffq_final_tmp, file = paste0(outpath, "ffq_fecal_", scfa, ".csv"), quote = F, row.names = F)
  
}

## FFQ final serum (plus covariates)
ffq_serum <- merge(fecal_scfas, serum_scfas, by = "subject_id", all = T) %>%
  merge(., anthropometrics, by = "subject_id", all = T) %>% # covariates
  merge(., ffq, by = "subject_id", all = T) %>%
  merge(., hei_ffq, by = "subject_id", all = T)

for (scfa in serum_scfa_list) {
  
  ## pull out SCFAs not being tested
  not_used_scfas <- serum_scfa_list[serum_scfa_list != scfa]
  ffq_serum_tmp <- ffq_serum %>%
    dplyr::select(., -dplyr::any_of(not_used_scfas), -dplyr::any_of(fecal_scfa_list)) %>%
    tidyr::drop_na()%>%
    janitor::clean_names()
  
  ## write file
  write.csv(x = ffq_serum_tmp, file = paste0(outpath, "ffq_serum_", scfa, ".csv"), quote = F, row.names = F)
  
}

#########
## ASA24 =======================================================================
#########

## ASA final fecal (plus covariates)
asa_fecal <- merge(fecal_scfas, serum_scfas, by = "subject_id", all = T) %>%
  merge(., anthropometrics, by = "subject_id", all = T) %>% # covariates
  merge(., stool_vars, by = "subject_id", all = T) %>% # covariates
  merge(., asa24, by = "subject_id", all = T) %>%
  merge(., pd_all, by = "subject_id", all = T) %>%
  merge(., pd_carbs, by = "subject_id", all = T) %>%
  merge(., pd_fat, by = "subject_id", all = T) %>%
  merge(., pd_fiber, by = "subject_id", all = T) %>%
  merge(., pd_protein, by = "subject_id", all = T) %>%
  merge(., hei_asa, by = "subject_id", all = T) %>%
  merge(., asa24_fiber, by = "subject_id", all = T)

for (scfa in fecal_scfa_list) {
  
  ## pull out SCFAs not being tested
  not_used_scfas <- fecal_scfa_list[fecal_scfa_list != scfa]
  asa_fecal_tmp <- asa_fecal %>%
    dplyr::select(., -dplyr::any_of(not_used_scfas), -dplyr::any_of(serum_scfa_list)) %>%
    tidyr::drop_na() %>%
    janitor::clean_names()
  
  ## write file
  write.csv(x = asa_fecal_tmp, file = paste0(outpath, "asa_fecal_", scfa, ".csv"), quote = F, row.names = F)
  
}

## ASA final serum (plus covariates)
asa_serum <- merge(fecal_scfas, serum_scfas, by = "subject_id", all = T) %>%
  merge(., anthropometrics, by = "subject_id", all = T) %>% # covariates
  merge(., asa24, by = "subject_id", all = T) %>%
  merge(., pd_all, by = "subject_id", all = T) %>%
  merge(., pd_carbs, by = "subject_id", all = T) %>%
  merge(., pd_fat, by = "subject_id", all = T) %>%
  merge(., pd_fiber, by = "subject_id", all = T) %>%
  merge(., pd_protein, by = "subject_id", all = T) %>%
  merge(., hei_asa, by = "subject_id", all = T) %>%
  merge(., asa24_fiber, by = "subject_id", all = T)

for (scfa in serum_scfa_list) {
  
  ## pull out SCFAs not being tested
  not_used_scfas <- serum_scfa_list[serum_scfa_list != scfa]
  asa_serum_tmp <- asa_serum %>%
    dplyr::select(., -dplyr::any_of(not_used_scfas), -dplyr::any_of(fecal_scfa_list)) %>%
    tidyr::drop_na() %>%
    janitor::clean_names()
  
  ## write file
  write.csv(x = asa_serum_tmp, file = paste0(outpath, "asa_serum_", scfa, ".csv"), quote = F, row.names = F)
  
}

##################
## Monosaccarhides =============================================================
##################

## Monosacc final fecal (plus covariates)
monosacc_fecal <- merge(fecal_scfas, serum_scfas, by = "subject_id", all = T) %>%
  merge(., anthropometrics, by = "subject_id", all = T) %>% # covariates
  merge(., stool_vars, by = "subject_id", all = T) %>% # covariates
  merge(., monosacc, by = "subject_id", all = T)

for (scfa in fecal_scfa_list) {
  
  ## pull out SCFAs not being tested
  not_used_scfas <- fecal_scfa_list[fecal_scfa_list != scfa]
  monosacc_fecal_tmp <- monosacc_fecal %>%
    dplyr::select(., -dplyr::any_of(not_used_scfas), -dplyr::any_of(serum_scfa_list)) %>%
    tidyr::drop_na() %>%
    janitor::clean_names()
  
  ## write file
  write.csv(x = monosacc_fecal_tmp, file = paste0(outpath, "monosacc_fecal_", scfa, ".csv"), quote = F, row.names = F)
  
}

## FFQ final serum (plus covariates)
monosacc_serum <- merge(fecal_scfas, serum_scfas, by = "subject_id", all = T) %>%
  merge(., anthropometrics, by = "subject_id", all = T) %>% # covariates
  merge(., monosacc, by = "subject_id", all = T)

for (scfa in serum_scfa_list) {
  
  ## pull out SCFAs not being tested
  not_used_scfas <- serum_scfa_list[serum_scfa_list != scfa]
  monosacc_serum_tmp <- monosacc_serum %>%
    dplyr::select(., -dplyr::any_of(not_used_scfas), -dplyr::any_of(fecal_scfa_list)) %>%
    tidyr::drop_na() %>%
    janitor::clean_names()
  
  ## write file
  write.csv(x = monosacc_serum_tmp, file = paste0(outpath, "monosacc_serum_", scfa, ".csv"), quote = F, row.names = F)
  
}

##################
## IFDP ========================================================================
##################

## IFDP final fecal (plus covariates)
ifdp_fecal <- merge(fecal_scfas, serum_scfas, by = "subject_id", all = T) %>%
  merge(., anthropometrics, by = "subject_id", all = T) %>% # covariates
  merge(., stool_vars, by = "subject_id", all = T) %>% # covariates
  merge(., ifdp, by = "subject_id", all = T)

for (scfa in fecal_scfa_list) {
  
  ## pull out SCFAs not being tested
  not_used_scfas <- fecal_scfa_list[fecal_scfa_list != scfa]
  ifdp_fecal_tmp <- ifdp_fecal %>%
    dplyr::select(., -dplyr::any_of(not_used_scfas), -dplyr::any_of(serum_scfa_list)) %>%
    tidyr::drop_na() %>%
    janitor::clean_names()
  
  ## write file
  write.csv(x = ifdp_fecal_tmp, file = paste0(outpath, "ifdp_fecal_", scfa, ".csv"), quote = F, row.names = F)
  
}

## FFQ final serum (plus covariates)
ifdp_serum <- merge(fecal_scfas, serum_scfas, by = "subject_id", all = T) %>%
  merge(., anthropometrics, by = "subject_id", all = T) %>% # covariates
  merge(., ifdp, by = "subject_id", all = T)

for (scfa in serum_scfa_list) {
  
  ## pull out SCFAs not being tested
  not_used_scfas <- serum_scfa_list[serum_scfa_list != scfa]
  ifdp_serum_tmp <- ifdp_serum %>%
    dplyr::select(., -dplyr::any_of(not_used_scfas), -dplyr::any_of(fecal_scfa_list)) %>%
    tidyr::drop_na() %>%
    janitor::clean_names()
  
  ## write file
  write.csv(x = ifdp_serum_tmp, file = paste0(outpath, "ifdp_serum_", scfa, ".csv"), quote = F, row.names = F)
  
}

##################
## humann pwys =================================================================
##################

## humann final fecal (plus covariates)
humann_fecal <- merge(fecal_scfas, serum_scfas, by = "subject_id", all = T) %>%
  merge(., anthropometrics, by = "subject_id", all = T) %>% # covariates
  merge(., stool_vars, by = "subject_id", all = T) %>% # covariates
  merge(., humann_pwys, by = "subject_id", all = T)

for (scfa in fecal_scfa_list) {
  
  ## pull out SCFAs not being tested
  not_used_scfas <- fecal_scfa_list[fecal_scfa_list != scfa]
  humann_fecal_tmp <- humann_fecal %>%
    dplyr::select(., -dplyr::any_of(not_used_scfas), -dplyr::any_of(serum_scfa_list)) %>%
    tidyr::drop_na() %>%
    janitor::clean_names()
  
  ## write file
  write.csv(x = humann_fecal_tmp, file = paste0(outpath, "humann_pwys_fecal_", scfa, ".csv"), quote = F, row.names = F)
  
}

## FFQ final serum (plus covariates)
humann_serum <- merge(fecal_scfas, serum_scfas, by = "subject_id", all = T) %>%
  merge(., anthropometrics, by = "subject_id", all = T) %>% # covariates
  merge(., humann_pwys, by = "subject_id", all = T)

for (scfa in serum_scfa_list) {
  
  ## pull out SCFAs not being tested
  not_used_scfas <- serum_scfa_list[serum_scfa_list != scfa]
  humann_serum_tmp <- humann_serum %>%
    dplyr::select(., -dplyr::any_of(not_used_scfas), -dplyr::any_of(fecal_scfa_list)) %>%
    tidyr::drop_na() %>%
    janitor::clean_names()
  
  ## write file
  write.csv(x = humann_serum_tmp, file = paste0(outpath, "humann_pwys_serum_", scfa, ".csv"), quote = F, row.names = F)
  
}

##################
## TaxaHFE =====================================================================
##################

## taxahfe final fecal (plus covariates)
taxaHFE_fecal <- merge(fecal_scfas, serum_scfas, by = "subject_id", all = T) %>%
  merge(., stool_vars_taxaHFE_cov, by = "subject_id", all = T) %>%
  merge(., anthropometrics, by = "subject_id", all = T) # covariates

for (scfa in fecal_scfa_list) {
  
  ## pull out SCFAs not being tested
  not_used_scfas <- fecal_scfa_list[fecal_scfa_list != scfa]
  taxaHFE_fecal_tmp <- taxaHFE_fecal %>%
    dplyr::select(., -dplyr::any_of(not_used_scfas), -dplyr::any_of(serum_scfa_list)) %>%
    tidyr::drop_na() %>%
    janitor::clean_names()
  
  ## write file
  write.csv(x = taxaHFE_fecal_tmp, file = paste0(path, "taxahfe_", scfa, ".csv"), quote = F, row.names = F)
  
}

## FFQ final serum (plus covariates)
taxaHFE_serum <- merge(fecal_scfas, serum_scfas, by = "subject_id", all = T) %>%
  merge(., anthropometrics, by = "subject_id", all = T) # covariates

for (scfa in serum_scfa_list) {
  
  ## pull out SCFAs not being tested
  not_used_scfas <- serum_scfa_list[serum_scfa_list != scfa]
  taxaHFE_serum_tmp <- taxaHFE_serum %>%
    dplyr::select(., -dplyr::any_of(not_used_scfas), -dplyr::any_of(fecal_scfa_list)) %>%
    tidyr::drop_na() %>%
    janitor::clean_names()
  
  ## write file
  write.csv(x = taxaHFE_serum_tmp, file = paste0(path, "taxahfe_", scfa, ".csv"), quote = F, row.names = F)
  
}

##################
## Inflammation markers ========================================================
##################

## taxahfe final fecal (plus covariates)
inflammation_fecal <-  merge(fecal_scfas, serum_scfas, by = "subject_id", all = T) %>%
  merge(., anthropometrics, by = "subject_id", all = T) %>% # covariates
  merge(., stool_vars, by = "subject_id", all = T) %>% # covariates
  merge(., inflammation_markers, by = "subject_id", all = T)

for (scfa in fecal_scfa_list) {
  
  ## pull out SCFAs not being tested
  not_used_scfas <- fecal_scfa_list[fecal_scfa_list != scfa]
  inflammation_fecal_tmp <- inflammation_fecal %>%
    dplyr::select(., -dplyr::any_of(not_used_scfas), -dplyr::any_of(serum_scfa_list)) %>%
    tidyr::drop_na() %>%
    janitor::clean_names()
  
  ## write file
  write.csv(x = inflammation_fecal_tmp, file = paste0(outpath, "inflammation_", scfa, ".csv"), quote = F, row.names = F)
  
}

## FFQ final serum (plus covariates)
inflammation_serum <- merge(fecal_scfas, serum_scfas, by = "subject_id", all = T) %>%
  merge(., anthropometrics, by = "subject_id", all = T) %>% # covariates
  merge(., inflammation_markers, by = "subject_id", all = T)

for (scfa in serum_scfa_list) {
  
  ## pull out SCFAs not being tested
  not_used_scfas <- serum_scfa_list[serum_scfa_list != scfa]
  inflammation_serum_tmp <- inflammation_serum %>%
    dplyr::select(., -dplyr::any_of(not_used_scfas), -dplyr::any_of(fecal_scfa_list)) %>%
    tidyr::drop_na() %>%
    janitor::clean_names()
  
  ## write file
  write.csv(x = inflammation_serum_tmp, file = paste0(outpath, "inflammation_", scfa, ".csv"), quote = F, row.names = F)
  
}

