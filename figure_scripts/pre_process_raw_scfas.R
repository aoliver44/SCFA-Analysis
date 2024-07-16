## pre process raw SCFA data

## pre_process_raw_scfas.R: preprocess SCFA data
##  - take mean of duplicated plasma SCFA measurements
##  - remove fecal samples based on fecal filters
##  - if plasma reading is NA for propionate and butyrate, make it 10
##    - this is the limit of detection, which tobit model will handle
##  - create plasma variables in moles to match fecal samples
##  - calculate fecal ratio

## Author: Andrew Oliver
## Date: Feb 23, 2024
## to run: docker run --rm -it -p 8787:8787 -e PASSWORD=yourpasswordhere -v `PWD`:/home/docker aoliver44/scfa_rstudio:1.0

library(dplyr)

`%!in%` <- Negate(`%in%`)
## read in PLASMA SCFA data
scfa_plasma <- readr::read_delim("data/plasma_scfas.csv", delim = ",", col_types = "fdddf")

## any missing or NA in propionate or butyrate, make it 10
## this is because you supply 10 as the value that is left-censored in 
## censored models.
scfa_plasma[c("p_butyric_acid", "p_propionic_acid")][is.na(scfa_plasma[c("p_butyric_acid", "p_propionic_acid")])] <- 10

## create a list of plasma samples we have multiples of
scfa_duplicated <- scfa_plasma %>% filter(., duplicates != "") %>% pull(., subject_id) %>% droplevels()
scfas <- c("p_acetic_acid", "p_propionic_acid", "p_butyric_acid")

## create a empty DF for which to place the averaged out duplicated samples
deduplicated <- data.frame(subject_id=character(), 
                           p_acetic_acid=numeric(), 
                           p_propionic_acid=numeric(),
                           p_butyric_acid=numeric())

## loop through duplicated samples
for (subject in unique(scfa_duplicated)) {
  
  ## for each duplicated sample, pull it out and all the duplicated SCFA data
  scfa_duplicated_tmp <- scfa_plasma %>% filter(., subject_id == subject) %>% 
    select(., subject_id, p_acetic_acid, p_propionic_acid, p_butyric_acid)
  
  ## check and make sure you dont have a SCFA that is completely missing data
  ## you shouldnt with this data
  if (max((sum(is.na(scfa_duplicated_tmp$p_acetic_acid))), 
          (sum(is.na(scfa_duplicated_tmp$p_propionic_acid))), 
          (sum(is.na(scfa_duplicated_tmp$p_butyric_acid)))) > NROW(scfa_duplicated_tmp)) {
    stop("You have duplicated subject_id with completely missing SCFA data for at least one SCFA")
  }
  
  ## create a vector of the mean SCFA abundance of the values you have,
  ## ignoring any NAs
  missing_means <- colMeans(scfa_duplicated_tmp[2:4], na.rm = T)
  count = 1
  
  ## add this SCFA mean to the formerly empty, deduplicated DF
  deduplicated <- deduplicated %>% tibble::add_row(., subject_id=subject, 
                  p_acetic_acid=missing_means[1], 
                  p_propionic_acid=missing_means[2],
                  p_butyric_acid=missing_means[3])
}

## remove duplicated from raw data and add back in averaged values
scfa_plasma_dedup <- scfa_plasma %>% dplyr::select(., -duplicates) %>%
  dplyr::filter(., subject_id %!in% scfa_duplicated)
scfa_plasma_dedup <- rbind(scfa_plasma_dedup, deduplicated)

## normalize to nmol/ul to be like fecal SCFAs (nmol / mg)
scfa_plasma_dedup$p_butyric_acid_nmol <- scfa_plasma_dedup$p_butyric_acid / (1000 * 88.11)
scfa_plasma_dedup$p_acetic_acid_nmol <- scfa_plasma_dedup$p_acetic_acid / (1000 * 60.052)
scfa_plasma_dedup$p_propionic_acid_nmol <- scfa_plasma_dedup$p_propionic_acid / (1000 * 74.08)

## total here
scfa_plasma_dedup$p_scfa_nmol_total <- scfa_plasma_dedup$p_butyric_acid_nmol + 
  scfa_plasma_dedup$p_acetic_acid_nmol + 
  scfa_plasma_dedup$p_propionic_acid_nmol
  
## create relative abundance values
scfa_plasma_dedup$p_butyric_acid_nmol_norm <- scfa_plasma_dedup$p_butyric_acid_nmol / scfa_plasma_dedup$p_scfa_nmol_total
scfa_plasma_dedup$p_acetic_acid_nmol_norm <- scfa_plasma_dedup$p_acetic_acid_nmol / scfa_plasma_dedup$p_scfa_nmol_total
scfa_plasma_dedup$p_propionic_acid_nmol_norm <- scfa_plasma_dedup$p_propionic_acid_nmol / scfa_plasma_dedup$p_scfa_nmol_total

## write to file
readr::write_delim("data/plasma_scfas_normalized.csv", delim = ",", x = scfa_plasma_dedup)

## READ IN FECAL SCFA ==========================================================
fecal_scfas <- readr::read_delim("data/fecal_scfa_fl100.csv", delim = ",", col_types = "fdddd") %>% dplyr::select(., -butyrate)
new_fecal_butyrate <- readr::read_delim("data/butyrate_isob_new_integration_12-06-22.csv", delim = ",", col_types = "fdd")
fecal_scfas <- merge(new_fecal_butyrate, fecal_scfas, by = "subject_id")
fecal_scfas <- fecal_scfas %>% rename(., "butyrate" = "new_butyrate", "isobutyrate" = "new_isobutyrate")
fecal_vars <- read.delim("/home/docker/data/FL100_stool_variables.txt")

## get rid of fecal samples that are >24 hrs or after visit 1
fecal_vars <- fecal_vars %>%
  dplyr::filter(., diff_time_hrs < 24) %>%
  dplyr::filter(., AfterV2 == 0)

fecal_scfas <- fecal_scfas %>% dplyr::filter(., subject_id %in% fecal_vars$subject_id)

## take relative abundane of SCFA
fecal_scfas$acetate_norm <- fecal_scfas$acetate / fecal_scfas$total_scfa
fecal_scfas$butyrate_norm <- fecal_scfas$butyrate / fecal_scfas$total_scfa
fecal_scfas$propionate_norm <- fecal_scfas$propionate / fecal_scfas$total_scfa

## dist to norm (60:20:20)
fecal_scfas$acetate_norm_ratio_dist <- fecal_scfas$acetate_norm - 0.6
fecal_scfas$butyrate_norm_ratio_dist <- fecal_scfas$butyrate_norm - 0.2
fecal_scfas$propionate_norm_ratio_dist <- fecal_scfas$propionate_norm - 0.2

## get rid of isobutyrate - further evidence makes it seem like it is not isobutyrate
fecal_scfas <- fecal_scfas %>% dplyr::select(., -isobutyrate)

