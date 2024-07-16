## Supplemental Table 1

## supplemental_table1.R: generate supp table 1 of SCFA paper
## Author: Andrew Oliver
## Date: Jul 16, 2024
## to run: docker run --rm -it -p 8787:8787 -e PASSWORD=yourpasswordhere -v `PWD`:/home/docker aoliver44/scfa_rstudio:1.0

## load libraries ==============================================================
library(dplyr)

## set working directory =======================================================
setwd("/home/docker")

## source SCFA data ============================================================
source("/home/docker/github/SCFA-Analysis/figure_scripts/pre_process_raw_scfas.R")
anthropometrics <- read.csv("/home/docker/data/FL100_age_sex_bmi.csv")

## make supp table 1 ===========================================================

fecal_scfa_anthro <- merge(anthropometrics, fecal_scfas, by = "subject_id")
fecal_bin_n <- fecal_scfa_anthro %>% dplyr::mutate(., bin = dplyr::case_when(sex.factor == "Male" & between(age, 18,33) & between(bmi, 0,24.99) ~ "1",
                                   sex.factor == "Male" & between(age, 18,33) & between(bmi, 25,29.99) ~ "2",
                                   sex.factor == "Male" & between(age, 18,33) & between(bmi, 30,50) ~ "3",
                                   sex.factor == "Male" & between(age, 34,49) & between(bmi, 0,24.99) ~ "4",
                                   sex.factor == "Male" & between(age, 34,49) & between(bmi, 25,29.99) ~ "5",
                                   sex.factor == "Male" & between(age, 34,49) & between(bmi, 30,50) ~ "6",
                                   sex.factor == "Male" & between(age, 50,66) & between(bmi, 0,24.99) ~ "7",
                                   sex.factor == "Male" & between(age, 50,66) & between(bmi, 25,29.99) ~ "8",
                                   sex.factor == "Male" & between(age, 50,66) & between(bmi, 30,50) ~ "9",
                                   sex.factor == "Female" & between(age, 18,33) & between(bmi, 0,24.99) ~ "10",
                                   sex.factor == "Female" & between(age, 18,33) & between(bmi, 25,29.99) ~ "11",
                                   sex.factor == "Female" & between(age, 18,33) & between(bmi, 30,50) ~ "12",
                                   sex.factor == "Female" & between(age, 34,49) & between(bmi, 0,24.99) ~ "13",
                                   sex.factor == "Female" & between(age, 34,49) & between(bmi, 25,29.99) ~ "14",
                                   sex.factor == "Female" & between(age, 34,49) & between(bmi, 30,50) ~ "15",
                                   sex.factor == "Female" & between(age, 50,66) & between(bmi, 0,24.99) ~ "16",
                                   sex.factor == "Female" & between(age, 50,66) & between(bmi, 25,29.99) ~ "17",
                                   sex.factor == "Female" & between(age, 50,66) & between(bmi, 30,50) ~ "18",
                                   TRUE ~ NA)) %>% 
  dplyr::group_by(bin) %>% dplyr::tally() %>% dplyr::arrange((as.numeric(bin)))

plasma_scfa_anthro <- merge(anthropometrics, scfa_plasma_dedup, by = "subject_id")
plasma_bin_n <- plasma_scfa_anthro %>% dplyr::mutate(., bin = dplyr::case_when(sex.factor == "Male" & between(age, 18,33) & between(bmi, 0,24.99) ~ "1",
                                                       sex.factor == "Male" & between(age, 18,33) & between(bmi, 25,29.99) ~ "2",
                                                       sex.factor == "Male" & between(age, 18,33) & between(bmi, 30,50) ~ "3",
                                                       sex.factor == "Male" & between(age, 34,49) & between(bmi, 0,24.99) ~ "4",
                                                       sex.factor == "Male" & between(age, 34,49) & between(bmi, 25,29.99) ~ "5",
                                                       sex.factor == "Male" & between(age, 34,49) & between(bmi, 30,50) ~ "6",
                                                       sex.factor == "Male" & between(age, 50,66) & between(bmi, 0,24.99) ~ "7",
                                                       sex.factor == "Male" & between(age, 50,66) & between(bmi, 25,29.99) ~ "8",
                                                       sex.factor == "Male" & between(age, 50,66) & between(bmi, 30,50) ~ "9",
                                                       sex.factor == "Female" & between(age, 18,33) & between(bmi, 0,24.99) ~ "10",
                                                       sex.factor == "Female" & between(age, 18,33) & between(bmi, 25,29.99) ~ "11",
                                                       sex.factor == "Female" & between(age, 18,33) & between(bmi, 30,50) ~ "12",
                                                       sex.factor == "Female" & between(age, 34,49) & between(bmi, 0,24.99) ~ "13",
                                                       sex.factor == "Female" & between(age, 34,49) & between(bmi, 25,29.99) ~ "14",
                                                       sex.factor == "Female" & between(age, 34,49) & between(bmi, 30,50) ~ "15",
                                                       sex.factor == "Female" & between(age, 50,66) & between(bmi, 0,24.99) ~ "16",
                                                       sex.factor == "Female" & between(age, 50,66) & between(bmi, 25,29.99) ~ "17",
                                                       sex.factor == "Female" & between(age, 50,66) & between(bmi, 30,50) ~ "18",
                                                       TRUE ~ NA)) %>% 
  dplyr::group_by(bin) %>% dplyr::tally() %>% dplyr::arrange((as.numeric(bin)))


supp_table1 <- data.frame(sex=c(rep("MALE", 9),  rep("FEMALE", 9)), 
                          age=rep(c(rep("18-33",3), rep("34-49",3), rep("50-65",3)), 2), 
                          bmi=rep(c("0-24.99", "25-29.99", "30-50"), 6), 
                          fecal_n=fecal_bin_n$n, 
                          plasma_n=plasma_bin_n$n)





