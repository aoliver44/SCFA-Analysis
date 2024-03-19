## Partial regression function

## partial_regression.R: function to perform partial regressions for SCFA paper
## Author: Andrew Oliver
## Date: Feb 22, 2024
## to run: singularity run -W `pwd` --bind `pwd`:/home/docker /share/lemaylab/aoliver/software/nutrition_tools.sif bash -c 'Rscript /home/docker/ML_script1.R'

## load libraries
library(dplyr)
library(VGAM)
library(bestNormalize)
library(car)

## partial regression function
PartialCorrelationNew <- function(scfas, independent, df, remove_outliers=FALSE) {
  
  if (scfas == "serum") {
    scfa_list <- c("p_acetic_acid_nmol", "p_propionic_acid_nmol", "p_butyric_acid_nmol", "p_scfa_nmol_total")
    formula <- paste0("normalized ~ ", independent, " + as.numeric(bmi) + as.numeric(age) + as.factor(sex.factor)")
    formula_independent <- paste0(independent, " ~ as.numeric(bmi) + as.numeric(age) + as.factor(sex.factor)")
    formula_scfa <- paste0("normalized ~ as.numeric(bmi) + as.numeric(age) + as.factor(sex.factor)")
    
    covariates <- c("age", "sex.factor", "bmi")
  } else if (scfas == "fecal") {
    scfa_list <- c("acetate","propionate", "butyrate", "isobutyrate", "total_scfa", "acetate_norm_ratio_dist", "propionate_norm_ratio_dist", "butyrate_norm_ratio_dist")
    formula <- paste0("normalized ~ ", independent, " + as.numeric(bmi) + as.numeric(age) + as.factor(sex.factor) + as.numeric(st_wt) + as.numeric(bristol_num)")
    covariates <- c("age", "sex.factor", "bmi", "st_wt", "bristol_num")
  } else {
    scfa_list <- c(scfas)
    formula <- paste0("normalized ~ ", independent, " + as.numeric(bmi) + as.numeric(age) + as.factor(sex.factor) + as.numeric(st_wt) + as.numeric(bristol_num)")
    covariates <- c("age", "sex.factor", "bmi", "st_wt", "bristol_num")
  }
  
  df_results <- data.frame(scfa=factor(), factor=factor(), cor_estimate=numeric(),
                           p_value=numeric(),p.adjust=numeric(), tobit_estimate=numeric(), 
                           tobit_pvalue=numeric())
  
  for (scfa in scfa_list) {
    
    print(paste0("############# ",scfa, " v. ", independent, " #############"))
    df <- df %>% dplyr::select(., dplyr::all_of(scfa_list), all_of(independent), all_of(covariates)) %>% tidyr::drop_na()
    df$normalized <- df[[scfa]]
    
    ## define outlier removal
    if (remove_outliers == TRUE) { df <- df %>% dplyr::filter(., df[[independent]] %!in% boxplot.stats(df[[independent]])$out) }
    
    model <- lm(formula, data = df)
    
    if (unlist(shapiro.test(resid(model))[2]) < 0.05) {
      
      df$normalized <- NULL
      set.seed(123)
      df$normalized <- bestNormalize::bestNormalize(df[[scfa]], allow_orderNorm = F, k = 10, r = 10)$x.t
      rm(model)
      model <- lm(formula, data = df)
      print(bestNormalize::bestNormalize(df[[scfa]], allow_orderNorm = F, k = 10, r = 10)$chosen_transform)
      print(shapiro.test(resid(model))[2])
      
      if (unlist(shapiro.test(resid(model))[2]) < 0.05) {
        print("Trying Order Norm...")
        df$normalized <- NULL
        set.seed(123)
        df$normalized <- bestNormalize::bestNormalize(df[[scfa]], allow_orderNorm = T, k = 10, r = 20)$x.t
        rm(model)
        model <- lm(formula, data = df)
        print(bestNormalize::bestNormalize(df[[scfa]], allow_orderNorm = T, k = 10, r = 20)$chosen_transform)
        print(shapiro.test(resid(model))[2])
        
        if (unlist(shapiro.test(resid(model))[2]) < 0.05) {
          print("Trying Lambert...")
          df$normalized <- NULL
          set.seed(123)
          df$normalized <- bestNormalize::bestNormalize(df[[scfa]], allow_orderNorm = T, allow_lambert_s = T, allow_lambert_h = T, k = 10, r = 20)$x.t
          rm(model)
          model <- lm(formula, data = df)
          print(bestNormalize::bestNormalize(df[[scfa]], allow_orderNorm = T, allow_lambert_s = T, allow_lambert_h = T, k = 10, r = 20)$chosen_transform)
          print(shapiro.test(resid(model))[2])
          }
        } 
    } else { 
      print("No normalization needed") 
      print(shapiro.test(resid(model))[2])
    }
    
    partial_regression <- car::avPlots(model, main = paste(scfa, "vs.", independent))
    assign(x = "partial_regression", value = partial_regression, envir = .GlobalEnv)
    cor_method = "pearson"
    
    tmp_factor <- cor.test(as.data.frame(partial_regression[1])[,1], as.data.frame(partial_regression[1])[,2], method = cor_method)
    tmp_bmi <- cor.test(partial_regression$`as.numeric(bmi)`[,1], partial_regression$`as.numeric(bmi)`[,2], method = cor_method)
    tmp_age <- cor.test(partial_regression$`as.numeric(age)`[,1], partial_regression$`as.numeric(age)`[,2], method = cor_method)
    tmp_sex <- cor.test(partial_regression$`as.factor(sex.factor)Male`[,1], partial_regression$`as.factor(sex.factor)Male`[,2], method = cor_method)
    
    if (scfa %in% c("p_acetic_acid_nmol", "p_scfa_nmol_total")) {
      df_anthro <- data.frame(scfa = c(rep(scfa, 4)), factor = c(independent, "BMI", "Age", "Sex"), 
                              cor_estimate = c(tmp_factor$estimate, tmp_bmi$estimate, tmp_age$estimate, tmp_sex$estimate), 
                              cor_p_value = c(tmp_factor$p.value, tmp_bmi$p.value, tmp_age$p.value, tmp_sex$p.value),
                              lm_p_value = c(summary(model)$coefficients[2,4], summary(model)$coefficients[3,4], summary(model)$coefficients[4,4], summary(model)$coefficients[5,4]),
                              n_ind = c(rep(nrow(df), 4)),
                              tobit_estimate=c(rep(NA, 4)), 
                              tobit_pvalue=c(rep(NA, 4)))
    }
    
    if (scfas != "serum") {
      tmp_stwt <- cor.test(partial_regression$`as.numeric(st_wt)`[,1], partial_regression$`as.numeric(st_wt)`[,2], method = cor_method)
      tmp_bristol <- cor.test(partial_regression$`as.numeric(bristol_num)`[,1], partial_regression$`as.numeric(bristol_num)`[,2], method = cor_method)
      df_anthro <- data.frame(scfa = c(rep(scfa, 6)), factor = c(independent, "BMI", "Age", "Sex", "stool_weight", "Bristol Score"), 
                              cor_estimate = c(tmp_factor$estimate, tmp_bmi$estimate, tmp_age$estimate, tmp_sex$estimate, tmp_stwt$estimate, tmp_bristol$estimate), 
                              cor_p_value = c(tmp_factor$p.value, tmp_bmi$p.value, tmp_age$p.value, tmp_sex$p.value, tmp_stwt$p.value, tmp_bristol$p.value),
                              lm_p_value = c(summary(model)$coefficients[2,4], summary(model)$coefficients[3,4], summary(model)$coefficients[4,4], summary(model)$coefficients[5,4], summary(model)$coefficients[6,4], summary(model)$coefficients[7,4]),
                              n_ind = c(rep(nrow(df), 6)),
                              tobit_estimate=c(rep(NA, 6)), 
                              tobit_pvalue=c(rep(NA, 6)))
      
    } else if (scfas == "serum" & scfa %in% c("p_propionic_acid_nmol", "p_butyric_acid_nmol")) {

      tobit_model <- VGAM::vglm(formula, tobit(Lower = min(as.numeric(na.omit(df$normalized)))), data = df)
      tobit_model_minus_scfa <- VGAM::vglm(formula_independent, tobit(Lower = min(as.numeric(na.omit(df$normalized)))), data = df)
      tobit_model_minus_independent<- VGAM::vglm(formula_scfa, tobit(Lower = min(as.numeric(na.omit(df$normalized)))), data = df)
      
      tmp_factor <- cor.test(tobit_model_minus_independent@residuals[,1], tobit_model_minus_scfa@residuals[,1])
      
      df_anthro <- data.frame(scfa = c(rep(scfa, 4)), factor = c(independent, "BMI", "Age", "Sex"), 
                              cor_estimate = c(tmp_factor$estimate, tmp_bmi$estimate, tmp_age$estimate, tmp_sex$estimate), 
                              cor_p_value = c(tmp_factor$p.value, tmp_bmi$p.value, tmp_age$p.value, tmp_sex$p.value),
                              lm_p_value = c(rep(NA, 4)),
                              n_ind = c(rep(nrow(df), 4)),
                              tobit_estimate=c(coef(summary(tobit_model))[,1][3], 
                                               coef(summary(tobit_model))[,1][4], 
                                               coef(summary(tobit_model))[,1][5], 
                                               coef(summary(tobit_model))[,1][6]), 
                              tobit_pvalue=c(coef(summary(tobit_model))[,4][3], 
                                             coef(summary(tobit_model))[,4][4], 
                                             coef(summary(tobit_model))[,4][5], 
                                             coef(summary(tobit_model))[,4][6]))
      
    }
    
    df_results <- rbind(df_results, df_anthro[1,1:8])
    
  }
  return(df_results)
}