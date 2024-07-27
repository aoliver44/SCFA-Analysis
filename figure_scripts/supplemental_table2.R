## Supplemental Table 2

## supplemental_table2.R: script to generate supp table 2 in SCFA paper
## Author: Andrew Oliver
## Date: Jul 16, 2024
## to run: docker run --rm -it -p 8787:8787 -e PASSWORD=yourpasswordhere -v `PWD`:/home/docker aoliver44/scfa_analysis:1.1

## set working directory =======================================================
setwd("/home/docker")

## source data =================================================================
source("/home/docker/github/SCFA-Analysis/figure_scripts/pre_process_raw_scfas.R")
source("/home/docker/github/SCFA-Analysis/figure_scripts/partial_regression.R")

hei_ffq <- readr::read_csv("/home/docker/data/HEI FFQ_scores_12072021.csv") %>%
  dplyr::select(., subject_id, hei_ffq_totalscore)
hei_asa <- readr::read_delim("/home/docker/data/FL100_HEI_n378.txt") %>%
  dplyr::select(., subject_id, hei_asa24_totalscore)
asa24_fiber <- readr::read_delim("/home/docker/data/ASA24_average_fiber_summary_variables.txt")
ffq_fiber_vars <- readr::read_csv("/home/docker/data/fibergroups_fl100cohort.csv") %>% 
  dplyr::select(., subject_id, fibe_per1000_ffq, dt_fibe, dt_fiber_sol)
pd_carbs <- readr::read_delim("/home/docker/data/qiime1_alphadiv_carb.txt") %>% 
  dplyr::select(., subject_id, PD_carbs)
pd_fiber <- readr::read_delim("/home/docker/data/qiime1_alphadiv_fiber.txt") %>% 
  dplyr::select(., subject_id, PD_fiber)
anthropometrics <- read.csv("/home/docker/data/FL100_age_sex_bmi.csv")
stool_vars <- readr::read_delim("/home/docker/data/FL100_stool_variables.txt") %>%
  dplyr::select(., subject_id, st_wt, fecal_calprotectin, StoolConsistencyClass, bristol_num)
ethnicity <- readr::read_csv("/home/docker/data/DEXA_ethnicities04272020.csv") #%>% dplyr::filter(., Ethnicity %in% c("White", "Hispanic")) %>% droplevels()

## Make Supp Table 2 ===========================================================

## Merge all data used in analysis
all_fiber_vars_scfa <- merge(fecal_scfas, hei_ffq, by = "subject_id", all = T) %>%
  merge(., hei_asa, by = "subject_id", all = T) %>%
  merge(., asa24_fiber, by = "subject_id", all = T) %>%
  merge(., ffq_fiber_vars, by = "subject_id", all = T) %>%
  merge(., pd_carbs, by = "subject_id", all = T) %>%
  merge(., pd_fiber, by = "subject_id", all = T) %>%
  merge(., scfa_plasma_dedup, by = "subject_id", all = T) %>%
  merge(., anthropometrics, by = "subject_id", all = T) %>%
  merge(., stool_vars, by = "subject_id", all = T) %>%
  merge(., ethnicity, by = "subject_id", all = T) %>%
  merge(., readr::read_csv("/home/docker/data/CRP_WBC_9102021.csv"), by = "subject_id", all = T)

## remove NA subject IDS 
all_fiber_vars_scfa <- all_fiber_vars_scfa %>%
  dplyr::filter(., subject_id != "NA")

## Create empty dataframe to store results
all_samples_fiber <- data.frame(scfa=factor(), factor=factor(), cor_estimate=numeric(),
                                p_value=numeric(),p.adjust=numeric(), tobit_estimate=numeric(), 
                                tobit_pvalue=numeric())

## run partial correlations of fiber variables on full dataset
for (fiber_var in c("hei_ffq_totalscore", "hei_asa24_totalscore", "avg_fibe_tnfs", "perKcal_fiber_tnfs", "fibe_per1000_ffq", "dt_fibe", "dt_fiber_sol", "PD_carbs", "PD_fiber")) {
  ## full dataset
  tmp_fecal <- PartialCorrelationNew(scfas = "fecal", independent = fiber_var, df = all_fiber_vars_scfa)
  tmp_serum <- PartialCorrelationNew(scfas = "serum", independent = fiber_var, df = all_fiber_vars_scfa)
  all_samples_fiber <- rbind(all_samples_fiber, tmp_fecal, tmp_serum)
}

## subset dataset to just individuals without frank inflammation
all_fiber_vars_scfa_noInflam <- all_fiber_vars_scfa %>%
  dplyr::filter(., subject_id != "NA") %>%
  dplyr::filter(., fecal_calprotectin < 100) %>%
  dplyr::filter(., crp_bd1 < 10000)

## create new empty dataframe to store data in
no_inflam_fiber <- data.frame(scfa=factor(), factor=factor(), cor_estimate=numeric(),p_value=numeric(),p.adjust=numeric())

## run partial correlations of fiber variables on no-inflammation dataset
for (fiber_var in c("hei_ffq_totalscore", "hei_asa24_totalscore", "avg_fibe_tnfs", "perKcal_fiber_tnfs", "fibe_per1000_ffq", "dt_fibe", "dt_fiber_sol", "PD_carbs", "PD_fiber")) {
  ## full dataset
  tmp_fecal <- PartialCorrelationNew(scfas = "fecal", independent = fiber_var, df = all_fiber_vars_scfa_noInflam)
  tmp_serum <- PartialCorrelationNew(scfas = "serum", independent = fiber_var, df = all_fiber_vars_scfa_noInflam)
  no_inflam_fiber <- rbind(no_inflam_fiber, tmp_fecal, tmp_serum)
}

## create variable to keep track of whether all samples or no-inflam samples were used
all_samples_fiber$inf_samples <- "kept"
no_inflam_fiber$inf_samples <- "dropped"

## combine datasets together
fiber_scfa <- rbind(all_samples_fiber, no_inflam_fiber)

## create a plasma or fecal SCFA factor
fiber_scfa <- fiber_scfa %>% mutate(., scfa_type = ifelse(scfa %in% c("p_acetic_acid_nmol","p_propionic_acid_nmol", "p_butyric_acid_nmol", "p_scfa_nmol_total"), "plasma", "fecal"))

## create a p-adjust pvalue by grouping by SCFA and inflammation group (full data vs no-inflamm)
## and correcting p-values within these groups
fiber_scfa <- fiber_scfa %>%
  dplyr::group_by(scfa, inf_samples) %>%
  dplyr::mutate(., regression_p.adjust = ifelse(scfa_type == "fecal" | scfa == "p_acetic_acid_nmol" | scfa == "p_scfa_nmol_total", p.adjust(lm_p_value, method='fdr'), p.adjust(tobit_pvalue, method='fdr'))) 

## write to file
write.csv(fiber_scfa, file = "/home/docker/plots/directed_hypotheses.csv", quote = F, row.names = F)
