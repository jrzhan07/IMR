#-----------------------------------------------------------------------------#
# File name: 1_CompleteSimulation.R
# Scope: Complete simulation study - need original data
#-----------------------------------------------------------------------------#

# [0] Initialisation ----
library(mice) # for ampute
library(foreach) # for using foreach loop
library(dplyr) 
library(rsimsum) # For computing performance measures with MCSE
source("0_functions.R")

set.seed(576819506)

# [1] run simulation study
results_0.1 <- foreach(i = 1:1000) %do% {    
  
  ## pmiss = 0.1
  ds_sim1 = get.amputation(amp.seed=i, pmiss=0.2, ds_ref=ds_OHCA_final, fullmodel_var, miss.vars, covariates_dep_vars)
  sim_res1 = get.analysis(run = i + 10, ds_ref = ds_OHCA_final, fullmodel_var, miss.vars, 
                          ds_amputed = ds_sim1$ds_amputed, M_missing = ds_sim1$M_missing,
                          ds_methods, ds_coeff, final_label_list)
  sim_res1$pmiss = 0.1
  
  sim_res <- sim_res1
  
  attr(sim_res, "seed")<-.Random.seed  
  cat(".")  
  if (i%%50==0) cat("\n")
  sim_res
}
save(results_0.1, file=paste0("7_pmiss0.1_simres_", today(), ".RData"))

results_0.2 <- foreach(i = 1:1000) %do% {    
  
  ## pmiss = 0.2
  ds_sim2 = get.amputation(amp.seed=i, pmiss=0.2, ds_ref=ds_OHCA_final, fullmodel_var, miss.vars, covariates_dep_vars)
  sim_res2 = get.analysis(run = i + 10, ds_ref = ds_OHCA_final, fullmodel_var, miss.vars, 
                          ds_amputed = ds_sim2$ds_amputed, M_missing = ds_sim2$M_missing,
                          ds_methods, ds_coeff, final_label_list)
  sim_res2$pmiss = 0.2

  sim_res <- sim_res2
  
  attr(sim_res, "seed")<-.Random.seed  
  cat(".")  
  if (i%%50==0) cat("\n")
  sim_res
}
save(results_0.2, file=paste0("7_pmiss0.2_simres_", today(), ".RData"))

results_0.3 <- foreach(i = 1:1000) %do% {    
  
  ## pmiss = 0.3
  ds_sim3 = get.amputation(amp.seed=i, pmiss=0.3, ds_ref=ds_OHCA_final, fullmodel_var, miss.vars, covariates_dep_vars)
  sim_res3 = get.analysis(run = i + 10, ds_ref = ds_OHCA_final, fullmodel_var, miss.vars, 
                          ds_amputed = ds_sim3$ds_amputed, M_missing = ds_sim3$M_missing,
                          ds_methods, ds_coeff, final_label_list)
  sim_res3$pmiss = 0.3
  
  sim_res <- sim_res3
  
  attr(sim_res, "seed")<-.Random.seed  
  cat(".")  
  if (i%%50==0) cat("\n")
  sim_res
  
}
save(results_0.3, file=paste0("7_pmiss0.3_simres_", today(), ".RData"))

results_0.4 <- foreach(i = 1:1000) %do% {    
  
  ## pmiss = 0.4
  ds_sim4 = get.amputation(amp.seed=i, pmiss=0.4, ds_ref=ds_OHCA_final, fullmodel_var, miss.vars, covariates_dep_vars)
  sim_res4 = get.analysis(run = i + 10, ds_ref = ds_OHCA_final, fullmodel_var, miss.vars, 
                          ds_amputed = ds_sim4$ds_amputed, M_missing = ds_sim4$M_missing,
                          ds_methods, ds_coeff, final_label_list)
  sim_res4$pmiss = 0.4
  
  sim_res <- sim_res4
  
  attr(sim_res, "seed")<-.Random.seed  
  cat(".")  
  if (i%%50==0) cat("\n")
  sim_res
  
}
save(results_0.4, file=paste0("7_pmiss0.4_simres_", today(), ".RData"))

results.df_0.1 = bind_rows(results_0.1, .id = "run")
results.df_0.2 = bind_rows(results_0.2, .id = "run")
results.df_0.3 = bind_rows(results_0.3, .id = "run")
results.df_0.4 = bind_rows(results_0.4, .id = "run")
results.df = rbind(results.df_0.1, results.df_0.2, results.df_0.3, results.df_0.4)
save(results.df, file=paste0("7_10to40_simres_", today(), ".RData"))
summary(results.df)
