#-----------------------------------------------------------------------------#
# File name: 0_functions.R
# Scope: functions for amputation and imputation
#-----------------------------------------------------------------------------#

# [0] Initialisation ----
library(tidyverse)
library(mice) # ampute and MI
library(missForest)

# Note: original data available upon request from the authors
load("~/ds_OHCA_final_20240610.RData")

outcome_var = "Bystander_CPR"
exposure_var = "Alert_issued"
other_covariate_var = "First_rhythm"
miss.vars = c("Age", "Arrest_witnessed_by", "Time_call")
covariates_dep_vars = c("Gender", "Location_of_arrest")
fullmodel_var = c(outcome_var, exposure_var, covariates_dep_vars, other_covariate_var, miss.vars)

# Labels
final_label_list <- list(
  Bystander_CPR = "Bystander CPR",
  Age = "Age",
  Gender = "Gender",
  Arrest_witnessed_by = "Witness type",
  Time_call = "Time of call",
  Location_of_arrest = "Arrest location",
  Alert_issued = "Alert issued",
  Alert_accepted = "Alert accepted",
  First_rhythm = "First rhythm"
)

ds_methods <- tibble(method = c("No data missing (reference)", "Complete case", "mean/mode", "M.ind (mean/mode)", "Missforest", "kNN"),
                     abbr = c("REF", "CC", "MM", "MI", "MF", "kNN"))

ds_coeff <- tibble(coeff = c("Alert_issuedYes", "Age", 
                             "Arrest_witnessed_byBystander - family",
                             "Arrest_witnessed_byBystander - healthcare provider", 
                             "Arrest_witnessed_byBystander - lay person", 
                             "Time_call06:00–18:59", "Time_call19:00–23:59",
                             "GenderMale", "Location_of_arrestPublic", 
                             "First_rhythmShockable", "First_rhythmUnknown"),
                   coeff_label = c("Alert Issued (Yes)", "Age", 
                                   "Witness Type (family)",
                                   "Witness Type (healthcare provider)", 
                                   "Witness Type (lay person)", 
                                   "Time call (06:00–18:59)", "Time call (19:00–23:59)",
                                   "Gender (Male)", "Arrest Location (Public)", 
                                   "First Rhythm (Shockable)", "First Rhythm (Unknown)"))

# [1] Generate data: create amputation ----
get.amputation <- function(amp.seed, pmiss, ds_ref, fullmodel_var, miss.vars, covariates_dep_vars){
  # X = matrix with all covariates (missing and complete)
  set.seed(amp.seed)
  ampdata <- ampute(data = ds_ref[, c(miss.vars, covariates_dep_vars)], 
                    prop = pmiss, 
                    patterns = matrix(c(0,1,1,1,1,
                                        1,0,1,1,1,
                                        1,1,0,1,1), ncol=length(c(miss.vars, covariates_dep_vars)), nrow=length(miss.vars), byrow=T),
                    freq = c(0.1, 0.45, 0.45), 
                    mech = "MAR", 
                    weights = matrix(c(0,0,0,1,0,
                                       0,0,0,0,1,
                                       0,0,0,0,1), ncol=length(c(miss.vars, covariates_dep_vars)), nrow=length(miss.vars), byrow=T)
  )
  M_matrix <- +is.na(ampdata$amp)
  M_missing <- M_matrix[, miss.vars]
  colnames(M_missing) <- paste0("m.ind_", colnames(M_missing))
  ds_amputed <- ds_ref %>% dplyr::select(all_of(fullmodel_var))
  ds_amputed[, miss.vars][M_missing == 1] <- NA
  
  return(list(ampdata=ampdata, ds_amputed=ds_amputed, M_missing=M_missing))
}


# [2] Impute and analyse data ----
# i.full data before deletetion
# ii. complete-case analysis
# iii. SI: ...

get.analysis <- function(run, ds_ref, fullmodel_var, miss.vars, ds_amputed, M_missing, ds_methods, ds_coeff, final_label_list, print.output = F) {

    ## [1] Imputation ----
    imp.seed <- run + 10
    set.seed(imp.seed)
    data_list <- list(ds_fulldata = ds_ref %>% dplyr::select(all_of(fullmodel_var)),
                      ds_CC       = na.omit(ds_amputed),
                      ds_imp1_mean_mode = ds_amputed %>% mutate(Age = Hmisc::impute(Age, fun = mean),
                                                                Arrest_witnessed_by = Hmisc::impute(Arrest_witnessed_by, fun = mode),
                                                                Time_call = Hmisc::impute(Time_call, fun = mode)),
                      ds_MI_mean_mode = NA,
                      ds_imp2_missForest = try(missForest::missForest(as.data.frame(ds_amputed))$ximp),
                      ds_imp3_kNN = try(VIM::kNN(ds_amputed, variable = miss.vars))
    )
    
    data_list[["ds_MI_mean_mode"]] = cbind(data_list$ds_imp1_mean_mode, M_missing)
    
    data_list <- lapply(data_list, function(ds) { 
      if (!inherits(ds, "try-error")) {
        ds = ds
      }else{
        ds = data_list$ds_fulldata %>% mutate(Bystander_CPR = NA)
      }
      
      for (var in names(ds_amputed)) {
        Hmisc::label(ds[[var]]) <- final_label_list[[var]]
      }
      ds$Arrest_witnessed_by <- relevel(ds$Arrest_witnessed_by, ref = "Not witnessed")
      return(ds)
    })
    
    ## [2] Analysis ----
    mod_fit_list <- lapply(data_list, function(ds) { 
      try(glm(Bystander_CPR ~ ., data=ds, family = "binomial"))
    })
    
    tidy_list <- lapply(mod_fit_list, function(model) {
      if (!inherits(model, "try-error")) {
        broom::tidy(model, exponentiate = F, conf.int = TRUE)
      }else{
        tibble(term = unique(ds_coeff$coeff),
               estimate = NA, std.error = NA, statistic = NA,
               p.value = NA, conf.low = NA, conf.high = NA)
      }
    }) 

    MSE_list <- lapply(mod_fit_list, function(mod){
      if (!inherits(mod, "try-error")) {
        mean(mod$residuals^2)
      }else{
        NA
      }
    })
    
    out_MSE <- as.data.frame(unlist(MSE_list)) %>% mutate(run = run, method = unique(ds_methods$method)) %>% rename("model_MSE" = "unlist(MSE_list)") 
    
    AIC_list <- lapply(mod_fit_list, function(mod){
      if (!inherits(mod, "try-error")) {
        AIC(mod)
      }else{
        NA
      }
    })
    
    out_AIC <- as.data.frame(unlist(AIC_list)) %>% mutate(run = run, method = unique(ds_methods$method)) %>% rename("model_AIC" = "unlist(AIC_list)") 
    
    out_table <- c()
    for (var_i in unique(ds_coeff$coeff)){
      for (i in seq_along(tidy_list)){
        values <- tidy_list[[i]] %>% filter(term == var_i) %>% 
          mutate(run = run,
                 method = unique(ds_methods$method)[i],
                 coeff = var_i) %>%
          dplyr::select(run, coeff, method, estimate, std.error, conf.low, conf.high)
        out_table <- rbind(out_table, values)
      }
    }
    
    out_table_final <- out_table %>% group_by(coeff) %>% 
      mutate(method = forcats::fct_relevel(method, rev(unique(ds_methods$method))),
             method_group = case_when(method %in% c("Missforest", "kNN") ~ "ML approach",
                                      grepl("M.ind", method) ~ "Missing-indicator",
                                      grepl("reference", method) ~ "Reference",
                                      TRUE ~ "Statistical approach"),
             true_value = estimate[method == "No data missing (reference)"],
             bias = estimate - true_value,
             estimate_SE = case_when(coeff == "Age" ~ paste0(sprintf("%.4f", estimate), " (", round(std.error, 3), ")"),
                                     T ~ paste0(sprintf("%.2f", estimate), " (", round(std.error, 3), ")")),
             coverage = ifelse(true_value >= conf.low & true_value <= conf.high, 1, 0))
    out_table_final <- merge(merge(out_table_final, out_MSE, by = c("run", "method")),
                             out_AIC, by = c("run", "method"))
    out_table_final$method_group <- factor(out_table_final$method_group, levels=c("Reference", "Statistical approach", "ML approach", "Missing-indicator"))
    return(out_table_final)
}

  
