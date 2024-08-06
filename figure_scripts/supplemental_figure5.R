## Figure 2

## figure2.R: generate figure 2 of SCFA paper
## Author: Andrew Oliver
## Date: Feb 23, 2024
## docker run --rm -it \
## -v ~/Downloads/SCFA-Analysis/figure_scripts:/home/scripts \
## -v ~/Downloads/SCFA-Analysis-DATA/data/:/home/data \
## -w /home/docker \
## scfa_analysis:rstudio bash -c "Rscript supplemental_figure5.R"

## load libraries ==============================================================
library(dplyr)
library(ggpubr)
library(ggplot2)
library(ggh4x)
library(wesanderson)

## set working directory =======================================================
setwd("/home/")

## source SCFA data ============================================================
set.seed(123)
source("/home/scripts/pre_process_raw_scfas.R")
source("/home/scripts/partial_regression.R")


## wrangle data or source helper functions =====================================
anthropometrics <- read.csv("/home/data/FL100_age_sex_bmi.csv")
inflammation_markers <- merge(readr::read_csv("/home/data/CTSC24532USDAWHNRCNu-GIMarkers7Oct2021_DATA_2021-10-07_1627.csv"), 
                              readr::read_csv("/home/data/CRP_WBC_9102021.csv"), 
                              by = "subject_id") %>%
  dplyr::select(., dplyr::any_of(c("subject_id", "fecal_calprotectin", "crp_bd1"))) %>%
  merge(., fecal_scfas, by = "subject_id", all = T) %>%
  merge(., scfa_plasma_dedup, by = "subject_id", all = T) %>%
  dplyr::mutate(., crp_bd1_cutoff = ifelse(crp_bd1 >= 10000, "high", "low")) %>%
  dplyr::mutate(., fecal_cal_cutoff = ifelse(fecal_calprotectin >= 100, "high", "low")) %>%
  reshape2::melt(., id.vars = c("subject_id", "fecal_cal_cutoff", "crp_bd1_cutoff")) %>%
  dplyr::filter(., variable %in% c("acetate","propionate", 
                                   "butyrate", 
                                   "total_scfa", "acetate_norm_ratio_dist", 
                                   "propionate_norm_ratio_dist", 
                                   "butyrate_norm_ratio_dist", "p_acetic_acid_nmol",
                                   "p_propionic_acid_nmol", "p_butyric_acid_nmol", 
                                   "p_scfa_nmol_total")) 

## clean up data, order factors
inflammation_markers$fecal_cal_cutoff <- factor(x = inflammation_markers$fecal_cal_cutoff, levels = c("low", "high"), ordered = T)
inflammation_markers$crp_bd1_cutoff <- factor(x = inflammation_markers$crp_bd1_cutoff, levels = c("low", "high"), ordered = T)
inflammation_markers$variable <- factor(x = inflammation_markers$variable, levels = c("acetate", 
                                                                                      "acetate_norm_ratio_dist", "propionate", 
                                                                                      "propionate_norm_ratio_dist", "butyrate", 
                                                                                      "butyrate_norm_ratio_dist", 
                                                                                      "total_scfa", "p_acetic_acid_nmol", "p_propionic_acid_nmol", "p_butyric_acid_nmol", "p_scfa_nmol_total"), ordered = T)
## create an overall SCFA factor for ratio and raw
inflammation_markers <- inflammation_markers %>%
  dplyr::mutate(overall_scfa = case_when(
    variable == "acetate" ~ 'fecal acetate',
    variable == "acetate_norm_ratio_dist" ~ 'fecal acetate',
    variable == "propionate" ~ 'fecal propionate',
    variable == "propionate_norm_ratio_dist" ~ 'fecal propionate',
    variable == "butyrate" ~ 'fecal butyrate',
    variable == "butyrate_norm_ratio_dist" ~ 'fecal butyrate',
    variable == "total_scfa" ~ 'fecal total SCFAs',
    variable == "p_acetic_acid_nmol" ~ 'serum acetate',
    variable == "p_propionic_acid_nmol" ~ 'serum propionate',
    variable == "p_butyric_acid_nmol" ~ 'serum butyrate',
    variable == "p_scfa_nmol_total" ~ 'serum total SCFAs',
    TRUE ~ 'what_is_this'      
  )
  )

## create a ratio/raw factor
inflammation_markers <- inflammation_markers %>%
  dplyr::mutate(ratio_raw = case_when(
    variable == "acetate" ~ 'raw',
    variable == "acetate_norm_ratio_dist" ~ 'ratio',
    variable == "propionate" ~ 'raw',
    variable == "propionate_norm_ratio_dist" ~ 'ratio',
    variable == "butyrate" ~ 'raw',
    variable == "butyrate_norm_ratio_dist" ~ 'ratio',
    variable == "total_scfa" ~ 'raw',
    variable == "p_acetic_acid_nmol" ~ 'raw',
    variable == "p_propionic_acid_nmol" ~ 'raw',
    variable == "p_butyric_acid_nmol" ~ 'raw',
    variable == "p_scfa_nmol_total" ~ 'raw',
    TRUE ~ 'what_is_this'      
  )
  )
## order ratios/raw
inflammation_markers$ratio_raw <- factor(x = inflammation_markers$ratio_raw, levels = c("raw", "ratio"), ordered = T)
## create DF for inflammation distribution
inflammation <- merge(readr::read_csv("/home/data/CTSC24532USDAWHNRCNu-GIMarkers7Oct2021_DATA_2021-10-07_1627.csv"), 
                      readr::read_csv("/home/data/CRP_WBC_9102021.csv"), 
                      by = "subject_id") %>%
  merge(., fecal_scfas, by = "subject_id", all = T) %>%
  merge(., scfa_plasma_dedup, by = "subject_id", all = T) %>%
  merge(., anthropometrics, by = "subject_id", all = T) %>%
  merge(., fecal_vars %>% dplyr::select(., subject_id, st_wt, StoolConsistencyClass, bristol_num), by = "subject_id", all = T)

## make supp fig 5 =============================================================

inflammation_subset <- subset(inflammation, inflammation$fecal_calprotectin < 100 & inflammation$crp_bd1 < 10000)

set.seed(123)
inflammation_prop <- inflammation_subset %>% dplyr::filter(., p_propionic_acid_nmol != "NA") %>% dplyr::select(., p_propionic_acid_nmol, fecal_mpo, age, sex.factor, bmi) %>% tidyr::drop_na()
inflammation_prop$p_propionic_acid_normalized <- bestNormalize::bestNormalize(inflammation_prop$p_propionic_acid_nmol, allow_orderNorm = T, k = 10, r = 10)$x.t
print(bestNormalize::bestNormalize(inflammation_prop$p_propionic_acid_normalized, allow_orderNorm = T, k = 10, r = 10)$chosen_transform)

tobit_model_minus_scfa <- VGAM::vglm(as.numeric(fecal_mpo) ~ as.numeric(bmi) + as.numeric(age) + as.factor(sex.factor), tobit(Lower = min(as.numeric(na.omit(inflammation_prop$p_propionic_acid_normalized)))), data = inflammation_prop)
tobit_model_minus_independent <- VGAM::vglm(as.numeric(p_propionic_acid_normalized) ~ as.numeric(bmi) + as.numeric(age) + as.factor(sex.factor), tobit(Lower = min(as.numeric(na.omit(inflammation_prop$p_propionic_acid_normalized)))), data = inflammation_prop)

set.seed(123)
tmp_fecal_ph <- cor.test(tobit_model_minus_scfa@residuals[,1], tobit_model_minus_independent@residuals[,1], use = "complete.obs")

p_propionate_mpo <- ggplot() + 
  aes(x = tobit_model_minus_scfa@residuals[,1], y = tobit_model_minus_independent@residuals[,1]) + 
  geom_point(alpha = 0.1) + 
  geom_smooth(method = "lm", color = "red", se = F, linetype = "dashed") + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), axis.text = element_text(colour = "black")) +
  annotate("text", x = 1500, y = 1.5, label= paste0("r = ", round(tmp_fecal_ph$estimate, 4), "\np = ", round(tmp_fecal_ph$p.value,4), "\np.adjust = 0.0002")) +
  labs(x = "Fecal MPO (normalized) | covariates", y = paste0("Plasma propionate (normalized)\n| covariates")) 


count = 1
results_df_cor <- PartialCorrelationNew(scfas = "fecal", independent = "plasma_lbp_bd1", df = subset(inflammation, inflammation$fecal_calprotectin < 100 & inflammation$crp_bd1 < 10000))
for (scfa_tmp in c("acetate", "propionate", "butyrate", "total_scfa")) {
  
  set.seed(123)
  inflammation_subset$normalized <- bestNormalize::bestNormalize(inflammation_subset[[scfa_tmp]], allow_orderNorm = T, k = 10, r = 10)$x.t
  model <- lm(normalized ~ as.numeric(plasma_lbp_bd1) + as.numeric(bmi) + as.numeric(age) + as.factor(sex.factor) + as.numeric(bristol_num) + as.numeric(st_wt), data = inflammation_subset)
  partial_regression <- car::avPlots(model)
  set.seed(123)
  tmp_fecal_ph <- cor.test(partial_regression$`as.numeric(plasma_lbp_bd1)`[,1], partial_regression$`as.numeric(plasma_lbp_bd1)`[,2])
 plot <- ggplot(as.data.frame(partial_regression$`as.numeric(plasma_lbp_bd1)`)) + 
    aes(x = `as.numeric(plasma_lbp_bd1)`, y = normalized) + 
    geom_point(alpha = 0.1) + 
    geom_smooth(method = "lm", color = wesanderson::wes_palette("AsteroidCity1")[count], se = F, linetype = "dashed") + 
    theme_bw() + 
    theme(panel.grid.minor = element_blank(), axis.text = element_text(colour = "black")) +
    annotate("text", x = 15, y = 1.5, label= 
               paste0("r = ", round(subset(results_df_cor, results_df_cor$scfa == scfa_tmp)$cor_estimate, 4), 
                      "\np = ", round(subset(results_df_cor, results_df_cor$scfa == scfa_tmp)$cor_p_value,4), 
                      "\np.adjust > 0.05")) +
    labs(x = "Plasma LBP (normalized) | covariates", y = paste0(scfa_tmp, " (normalized) | covariates"))
 
 assign(x = paste0(scfa_tmp, "_plasma_lbp_bd1"), value = plot)
 count = count + 1
}

design = "ABC
DEF"

p_propionate_mpo + acetate_plasma_lbp_bd1 + propionate_plasma_lbp_bd1 + butyrate_plasma_lbp_bd1 + total_scfa_plasma_lbp_bd1 + p_propionate_mpo + patchwork::plot_layout(design = design) + patchwork::plot_annotation(tag_levels = "A")

dir.create(path = "/home/scripts/output_figures", showWarnings = TRUE)
ggsave(filename = "/home/scripts/output_figures/supplemental_figure5.png", 
       device = "png",
       width = 10, 
       height = 7,
       units = "in",
       dpi = 400)
dev.off()
