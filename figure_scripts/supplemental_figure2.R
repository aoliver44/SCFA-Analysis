## Figure 2

## figure2.R: generate figure 2 of SCFA paper
## Author: Andrew Oliver
## Date: March 18, 2024
## to run: docker run --rm -it -p 8787:8787 -e PASSWORD=yourpasswordhere -v `PWD`:/home/docker aoliver44/scfa_analysis:1.1

## load libraries ==============================================================
library(ggplot2)
library(patchwork)

## set working directory =======================================================
setwd("/home/docker")

## source helper scripts  ======================================================
set.seed(123)
source("/home/docker/github/SCFA-Analysis/figure_scripts/pre_process_raw_scfas.R")
source("/home/docker/github/SCFA-Analysis/figure_scripts/partial_regression.R")

## Read in data and wrangle   ==================================================
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

all_fiber_vars_scfa_noInflam <- all_fiber_vars_scfa %>%
  dplyr::filter(., subject_id != "NA") %>%
  dplyr::filter(., fecal_calprotectin < 100) %>%
  dplyr::filter(., crp_bd1 < 10000)

all_fiber_vars_scfa$butyrate_norm_ratio_dist_normalized <- bestNormalize::bestNormalize(all_fiber_vars_scfa$butyrate_norm_ratio_dist, allow_orderNorm = F, k = 10, r = 10)$x.t

model <- lm(as.numeric(butyrate_norm_ratio_dist_normalized) ~ as.numeric(hei_ffq_totalscore) + as.numeric(bmi) + as.numeric(age) + as.factor(sex.factor) + as.numeric(bristol_num) + as.numeric(st_wt), data = all_fiber_vars_scfa)
partial_regression_tmp <- car::avPlots(model)

tmp_fecal_ph <- cor.test(partial_regression_tmp$`as.numeric(hei_ffq_totalscore)`[,1], partial_regression_tmp$`as.numeric(hei_ffq_totalscore)`[,2])

tmp_annotation_data <- PartialCorrelationNew(scfas = "butyrate_norm_ratio_dist", independent = "hei_ffq_totalscore", df = all_fiber_vars_scfa)

butyrate_hei_ffq <- ggplot(as.data.frame(partial_regression_tmp$`as.numeric(hei_ffq_totalscore)`)) + 
  aes(x = `as.numeric(hei_ffq_totalscore)`, y = `as.numeric(butyrate_norm_ratio_dist_normalized)`) + 
  geom_point(alpha = 0.1) + 
  geom_smooth(method = "lm", color = "red", se = F, linetype = "dashed") + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), axis.text = element_text(colour = "black")) +
  annotate("text", x = 15, y = -2, label= paste0("r = ", round(tmp_annotation_data$cor_estimate, 4), "\np = ", round(tmp_annotation_data$cor_p_value,4), "\np.adjust = 0.013")) + 
  labs(x = "HEI FFQ Totalscore (normalized) | covariates", y = paste0("Fecal butyrate-ratio (normalized)\n| covariates"))

set.seed(123)
all_fiber_vars_scfa_noInflam$butyrate_norm_ratio_dist_normalized <- bestNormalize::bestNormalize(all_fiber_vars_scfa_noInflam$butyrate_norm_ratio_dist, allow_orderNorm = F, k = 10, r = 10)$x.t
print(bestNormalize::bestNormalize(all_fiber_vars_scfa_noInflam$butyrate_norm_ratio_dist_normalized, allow_orderNorm = T, k = 10, r = 10)$chosen_transform)

model <- lm(as.numeric(butyrate_norm_ratio_dist_normalized) ~ as.numeric(hei_ffq_totalscore) + as.numeric(bmi) + as.numeric(age) + as.factor(sex.factor) + as.numeric(bristol_num) + as.numeric(st_wt), data = all_fiber_vars_scfa_noInflam)
partial_regression_tmp <- car::avPlots(model)

tmp_fecal_ph <- cor.test(partial_regression_tmp$`as.numeric(hei_ffq_totalscore)`[,1], partial_regression_tmp$`as.numeric(hei_ffq_totalscore)`[,2])

tmp_annotation_data <- PartialCorrelationNew(scfas = "butyrate_norm_ratio_dist", independent = "hei_ffq_totalscore", df = all_fiber_vars_scfa_noInflam)

butyrate_hei_ffq_no_inflam <- ggplot(as.data.frame(partial_regression_tmp$`as.numeric(hei_ffq_totalscore)`)) + 
  aes(x = `as.numeric(hei_ffq_totalscore)`, y = `as.numeric(butyrate_norm_ratio_dist_normalized)`) + 
  geom_point(alpha = 0.1) + 
  geom_smooth(method = "lm", color = "red", se = F, linetype = "dashed") + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), axis.text = element_text(colour = "black")) +
  annotate("text", x = 15, y = -3, label= paste0("r = ", round(tmp_annotation_data$cor_estimate, 4), "\np = ", round(tmp_annotation_data$cor_p_value,4), "\np.adjust = 0.064")) +
  labs(x = "HEI FFQ Totalscore (normalized) | covariates\nNo high inflammation individuals", y = paste0("Fecal butyrate-ratio (normalized)\n| covariates"))  

set.seed(123)
all_fiber_vars_scfa$butyrate_norm_ratio_dist_normalized <- bestNormalize::bestNormalize(all_fiber_vars_scfa$butyrate_norm_ratio_dist, allow_orderNorm = F, k = 10, r = 10)$x.t
print(bestNormalize::bestNormalize(all_fiber_vars_scfa$butyrate_norm_ratio_dist_normalized, allow_orderNorm = T, k = 10, r = 10)$chosen_transform)

model <- lm(as.numeric(butyrate_norm_ratio_dist_normalized) ~ as.numeric(fibe_per1000_ffq) + as.numeric(bmi) + as.numeric(age) + as.factor(sex.factor) + as.numeric(bristol_num) + as.numeric(st_wt), data = all_fiber_vars_scfa)
partial_regression_tmp <- car::avPlots(model)

tmp_fecal_ph <- cor.test(partial_regression_tmp$`as.numeric(fibe_per1000_ffq)`[,1], partial_regression_tmp$`as.numeric(fibe_per1000_ffq)`[,2])

tmp_annotation_data <- PartialCorrelationNew(scfas = "butyrate_norm_ratio_dist", independent = "fibe_per1000_ffq", df = all_fiber_vars_scfa)

butyrate_fiber_cal_corrected <- ggplot(as.data.frame(partial_regression_tmp$`as.numeric(fibe_per1000_ffq)`)) + 
  aes(x = `as.numeric(fibe_per1000_ffq)`, y = `as.numeric(butyrate_norm_ratio_dist_normalized)`) + 
  geom_point(alpha = 0.1) + 
  geom_smooth(method = "lm", color = "red", se = F, linetype = "dashed") + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), axis.text = element_text(colour = "black")) +
  annotate("text", x = 10, y = -2.5, label= paste0("r = ", round(tmp_annotation_data$cor_estimate, 4), "\np = ", round(tmp_annotation_data$cor_p_value,4), "\np.adjust > 0.05")) +
  labs(x = "Fiber per 1000 cal (FFQ) (normalized) | covariates", y = paste0("Fecal butyrate-ratio (normalized)\n| covariates")) 

set.seed(123)
all_fiber_vars_scfa$butyrate_norm_ratio_dist_normalized <- bestNormalize::bestNormalize(all_fiber_vars_scfa$butyrate_norm_ratio_dist, allow_orderNorm = F, k = 10, r = 10)$x.t
print(bestNormalize::bestNormalize(all_fiber_vars_scfa$butyrate_norm_ratio_dist_normalized, allow_orderNorm = T, k = 10, r = 10)$chosen_transform)

model <- lm(as.numeric(butyrate_norm_ratio_dist_normalized) ~ as.numeric(hei_asa24_totalscore) + as.numeric(bmi) + as.numeric(age) + as.factor(sex.factor) + as.numeric(bristol_num) + as.numeric(st_wt), data = all_fiber_vars_scfa)
partial_regression_tmp <- car::avPlots(model)

tmp_fecal_ph <- cor.test(partial_regression_tmp$`as.numeric(hei_asa24_totalscore)`[,1], partial_regression_tmp$`as.numeric(hei_asa24_totalscore)`[,2])

tmp_annotation_data <- PartialCorrelationNew(scfas = "butyrate_norm_ratio_dist", independent = "hei_asa24_totalscore", df = all_fiber_vars_scfa)

butyrate_fiber_hei_asa <- ggplot(as.data.frame(partial_regression_tmp$`as.numeric(hei_asa24_totalscore)`)) + 
  aes(x = `as.numeric(hei_asa24_totalscore)`, y = `as.numeric(butyrate_norm_ratio_dist_normalized)`) + 
  geom_point(alpha = 0.1) + 
  geom_smooth(method = "lm", color = "red", se = F, linetype = "dashed") + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), axis.text = element_text(colour = "black")) +
  annotate("text", x = 20, y = -2, label= paste0("r = ", round(tmp_annotation_data$cor_estimate, 4), "\np = ", round(tmp_annotation_data$cor_p_value,4), "\np.adjust > 0.05")) +
  labs(x = "HEI ASA24 Totalscore (normalized) | covariates", y = paste0("Fecal butyrate-ratio (normalized)\n| covariates")) 

set.seed(123)
all_fiber_vars_scfa_prop <- all_fiber_vars_scfa %>% dplyr::filter(., p_propionic_acid_nmol != "NA") %>% dplyr::select(., p_propionic_acid_nmol, dt_fiber_sol, age, sex.factor, bmi) %>% tidyr::drop_na()
all_fiber_vars_scfa_prop$p_propionic_acid_normalized <- bestNormalize::bestNormalize(all_fiber_vars_scfa_prop$p_propionic_acid_nmol, allow_orderNorm = T, k = 10, r = 10)$x.t
print(bestNormalize::bestNormalize(all_fiber_vars_scfa_prop$p_propionic_acid_normalized, allow_orderNorm = T, k = 10, r = 10)$chosen_transform)

tobit_model_minus_scfa <- VGAM::vglm(as.numeric(dt_fiber_sol) ~ as.numeric(bmi) + as.numeric(age) + as.factor(sex.factor), tobit(Lower = min(as.numeric(na.omit(all_fiber_vars_scfa_prop$p_propionic_acid_normalized)))), data = all_fiber_vars_scfa_prop)
tobit_model_minus_independent <- VGAM::vglm(as.numeric(p_propionic_acid_normalized) ~ as.numeric(bmi) + as.numeric(age) + as.factor(sex.factor), tobit(Lower = min(as.numeric(na.omit(all_fiber_vars_scfa_prop$p_propionic_acid_normalized)))), data = all_fiber_vars_scfa_prop)

set.seed(123)
tmp_fecal_ph <- cor.test(tobit_model_minus_scfa@residuals[,1], tobit_model_minus_independent@residuals[,1], use = "complete.obs")

p_propionate_sol_fiber <- ggplot() + 
  aes(x = tobit_model_minus_scfa@residuals[,1], y = tobit_model_minus_independent@residuals[,1]) + 
  geom_point(alpha = 0.1) + 
  geom_smooth(method = "lm", color = "red", se = F, linetype = "dashed") + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), axis.text = element_text(colour = "black")) +
  annotate("text", x = 10, y = -1.5, label= paste0("r = ", round(tmp_fecal_ph$estimate, 4), "\np = ", round(tmp_fecal_ph$p.value,4), "\np.adjust > 0.05")) +
  labs(x = "Soluble Fiber (FFQ) (normalized) | covariates", y = paste0("Plasma propionate (normalized)\n| covariates")) 

design = "ABC
DEF"

butyrate_hei_ffq + butyrate_hei_ffq_no_inflam + butyrate_fiber_cal_corrected + butyrate_fiber_hei_asa + p_propionate_sol_fiber + p_propionate_sol_fiber + patchwork::plot_layout(design = design) + patchwork::plot_annotation(tag_levels = "A")

