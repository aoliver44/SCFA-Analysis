## Supplemental Figure3

## supplemental_figure3.R: generate supp figure 3 of SCFA paper
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

anthropometrics <- read.csv("/home/docker/data/FL100_age_sex_bmi.csv")
food_codes <- readr::read_delim(file = "/home/docker/data/food_tree_data/fl100_newick_taxonomy_nowater.txt", delim = "\t")
food_tree <- readr::read_delim(file = "/home/docker/data/food_tree_data/fl100_otu_abundance.txt", delim = "\t")
fndds <- readxl::read_xlsx(path = "/home/docker/data/food_tree_data/2015-2016 FNDDS At A Glance - FNDDS Nutrient Values.xlsx", skip = 1)
stool_vars <- readr::read_delim("/home/docker/data/FL100_stool_variables.txt")
ingredient_butyrate <- readr::read_csv("/home/docker/data/ingredient_butyrate_073123.csv")
asa_data_qc <- readr::read_delim("/home/docker/data/fl100_recalls_qcd.csv", delim = ",")

## wrangle data or source helper functions =====================================

## read in data

## change vars to dates
stool_vars$stool_collection <- as.Date(stool_vars$stool_collection, "%m/%d/%y")
stool_vars$visit2_date <- as.Date(stool_vars$visit2_date, "%m/%d/%y")

## put in stool filters
stool_vars <- stool_vars %>% dplyr::filter(., AfterV2 == 0) %>% dplyr::filter(., diff_time_hrs < 24)

## use un-averaged asa recalls to find record closest to poop
asa_data_qc <- asa_data_qc %>% tidyr::separate(., col = "Occ_Time", into = c("Occ_Time", "time_hrs"), extra = "merge", sep = " ")
asa_data_qc$Occ_Time <- as.Date(asa_data_qc$Occ_Time, "%m/%d/%y")
asa_data_qc$FoodCode <- paste0("fc_", asa_data_qc$FoodCode)

## merge ASA dietary records with stool records
stool_asa <- merge(asa_data_qc, stool_vars, by.x = "UserName", by.y = "subject_id", no.dups = F)

## calculate time difference between dietary record and stool collection
stool_asa$stool_collect_dietary_collect <- difftime(stool_asa$stool_collection, stool_asa$Occ_Time, units = "days")

## grab the individuals and the dietary records with the smallest positive difference (food record before stool sample)
## basically just removing anything with a negative day time stamp
best_records <- stool_asa %>% 
  dplyr::filter(stool_collect_dietary_collect > -1) %>% 
  dplyr::group_by(., UserName) %>% 
  dplyr::summarise(., min(stool_collect_dietary_collect)) %>%
  dplyr::rename(., "stool_collect_dietary_collect" = `min(stool_collect_dietary_collect)`)

stool_asa_filter <- merge(best_records, stool_asa, by = c("UserName", "stool_collect_dietary_collect"))

## summarize ASA butyrate per individual for that most close to poop day
stool_asa_filter_1day <- stool_asa_filter %>% 
  dplyr::select(., UserName, S040, RecallNo, stool_collect_dietary_collect) %>% 
  dplyr::group_by(., UserName, RecallNo, stool_collect_dietary_collect) %>% 
  dplyr::summarise(., sum_butyrate_bulk = sum(S040)) %>%
  rename(., "subject_id" = "UserName")

## summarize ingredientized ASA butyrate per individual for that most close to poop day
ingredient_butyrate <- ingredient_butyrate %>% 
  dplyr::select(., UserName, butyrate_consumed_g, RecallNo) %>%
  dplyr::rename(., "RecallNo_Ingred" = "RecallNo") %>%
  dplyr::rename(., "subject_id" = "UserName")

## merge bulk with ingredient
ingredient_butyrate_1day <- merge(stool_asa_filter_1day, ingredient_butyrate, by = "subject_id", no.dups = F)
ingredient_butyrate_1day <- ingredient_butyrate_1day %>% dplyr::filter(., RecallNo == RecallNo_Ingred)
ingredient_butyrate_1day <- ingredient_butyrate_1day %>% 
  dplyr::group_by(., subject_id) %>% 
  summarize(., sum_ingred = sum(butyrate_consumed_g))

close_to_poop <- merge(stool_asa_filter_1day, ingredient_butyrate_1day, by = "subject_id")

## remove samples who's nearest dietary record is more than 3 days before collection of
## poop. This is probably too far to relate diet record to poop.
close_to_poop <- close_to_poop %>% dplyr::filter(., stool_collect_dietary_collect < 3)
close_to_poop <- close_to_poop %>% tidyr::drop_na()

## test how close the one day is
cor.test(close_to_poop$sum_butyrate_bulk, close_to_poop$sum_ingred)

## merge with scfa data
butyrate_1_day <- merge(close_to_poop, fecal_scfas, by = "subject_id")
butyrate_1_day <- merge(butyrate_1_day, scfa_plasma_dedup, by = "subject_id") %>%
  dplyr::rename(., "Mixed_dish_butyrate" = "sum_butyrate_bulk",  "Ingredientized_butyrate" = "sum_ingred")

## melt into long form for ggplot
butyrate_melt <- butyrate_1_day %>% 
  dplyr::select(., subject_id, Mixed_dish_butyrate, Ingredientized_butyrate, butyrate, p_butyric_acid_nmol, butyrate_norm, p_butyric_acid_nmol_norm) %>% 
  dplyr::rename(., "plasma butyrate\nrelative abundance" = "p_butyric_acid_nmol_norm", "plasma butyrate" = "p_butyric_acid_nmol", "fecal butyrate" = "butyrate", "fecal butyrate\nrelative abundance" = "butyrate_norm") %>%
  reshape2::melt(., id.vars = c("subject_id", "Ingredientized_butyrate", "Mixed_dish_butyrate"))

## make supp figure 3 ==========================================================

ggplot(data = butyrate_melt) +
  aes(x = Mixed_dish_butyrate, y = value) + 
  geom_point() + 
  facet_wrap(variable ~ ., scales = "free", nrow = 1) + 
  theme_bw() + 
  ggpubr::stat_cor(size = 2.5) + 
  labs(y = "Butyrate (nmol / mg or ml)", x = "Dietary Butyrate\n(quantified from ASA24 dietary recals)")



