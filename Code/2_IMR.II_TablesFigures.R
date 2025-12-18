#-----------------------------------------------------------------------------#
# File name: 2_IMR.II_TablesFigures.R
# Scope: Generate tables and figures for manuscript
#-----------------------------------------------------------------------------#

# [0] Initialisation ----
library(tidyr)
library(forcats)
library(ggplot2)
library(patchwork)
library(table1)
library(rsimsum) # For computing performance measures with MCSE

source("0_functions.R")
# source("0_checkfailure (Paper II).R")
## For Paper I ----
load("7_10to40_simres_2024-11-14.RData")
results.df_PaperI = results.df 

## For Paper II ----
load("7_20to40_simres_2025-12-13.RData")
results.df_PaperII = results.df %>% filter(method != "No data missing (reference)") %>%
                                     rename("model_MSE" = "MSE", "model_AIC"="AIC") %>%
                                     select(all_of(names(results.df_PaperI)))
results.df = bind_rows(results.df_PaperI, results.df_PaperII) %>% filter(pmiss != 0.1)
ds_methods_all = rbind(ds_methods, ds_methods_MI[-1,])


col1 <- c("darkgrey", rep("black", 5), rep("darkblue", 3)) 
color_map <- setNames(col1, ds_methods_all$abbr)
color_map_pmiss <- setNames(c("darkgrey", "#0072B2", "#E69F00", "#2ca02c", "#8E44AD"), seq(0,40, by=10))
color_map_m = c("5" = "firebrick", "10" = "dodgerblue", "15" = "forestgreen", "20" = "gold")
n.sim = 1000
missvar_pattern <- "Age|Arrest_witnessed_by|Time_call"
obsvar_pattern <- "Alert|Gender|Arrest_witnessedYes|Location_of_arrest|First_rhythm"


# [1] Performance measure ----
s <- simsum(data = results.df, estvarname = "estimate", true = "true_value", se = "std.error", 
            methodvar = "method", ref = "No data missing (reference)", by=c("pmiss", "coeff"))
summary_s <- s$summ


# [2] Outputs ----
results.df_plot <- merge(merge(results.df, ds_coeff, by = "coeff"),
                         ds_methods_all %>% rename("method_abbr" = "abbr"), by = "method") %>%
                   mutate(coeff_label = case_when(grepl("Time call", coeff_label) ~ sub("Time call", "Call Time", coeff_label),
                                                  T ~ coeff_label),
                          method_abbr = ifelse(method_abbr == "kNN", "KNN", method_abbr),
                          method_pmiss = paste0(method, pmiss),
                          miss.prop = ifelse(method_abbr == "REF", 0, pmiss*100))
results.df_plot$miss.prop <- factor(results.df_plot$miss.prop, levels = sort(unique(results.df_plot$miss.prop)))
results.df_plot$method_abbr <- factor(results.df_plot$method_abbr, levels = c("REF", "CC", "MM", "MI", "MF", "KNN", "PMM", "RF", "CART"),
                                      labels = c("REF", "CC", "MM", "MxI", "MF", "KNN", "PMM", "RF", "CART"))
results.df_plot_SI = results.df_plot %>% filter(!method_abbr %in% c("PMM", "RF", "CART" ))

results.df_plot_sum <- results.df_plot %>% group_by(miss.prop, pmiss, coeff, coeff_label, method, method_abbr, method_group, true_value) %>% 
                                            mutate(run = as.numeric(run)) %>%
                                            summarize(estimate_mean = mean(estimate), 
                                                       conf.low_mean = mean(conf.low), conf.high_mean = mean(conf.high),
                                                       SE_mean = mean(std.error),
                                                       cov = sum(coverage)/n.sim*100,
                                                       estimate_var = var(estimate)) %>% ungroup()
ds_plot_sum_final <- results.df_plot_sum %>% mutate(estimate_SE = case_when(coeff == "Age" ~ paste0(sprintf("%.4f", estimate_mean), " (", round(SE_mean, 3), ")"),
                                                                            grepl("ime_call", coeff) ~ paste0(sprintf("%.3f", estimate_mean), " (", round(SE_mean, 3), ")"),
                                                                            T ~ paste0(sprintf("%.2f", estimate_mean), " (", round(SE_mean, 3), ")")),
                                                   bias = estimate_mean - true_value, 
                                                   MSE = bias^2 + estimate_var) %>% 
                                             mutate(coeff_label = case_when(coeff_label == "Witness Type (family)" ~ "Witness Type (Family)",
                                                                            coeff_label == "Witness Type (healthcare provider)" ~ "Witness Type (Healthcare provider)",
                                                                            coeff_label == "Witness Type (lay person)" ~ "Witness Type (Lay person)",
                                                                            T ~ coeff_label))


## Figure 2: beta ----
p1a = ggplot(ds_plot_sum_final %>% filter(grepl(missvar_pattern, coeff)), 
            aes(x = method_abbr, y = estimate_mean, group = miss.prop)) + 
            facet_wrap(~ coeff_label, nrow = 2, scales = "free_y") +
            geom_hline(aes(yintercept = true_value), linetype = "dashed", col = "grey", linewidth=1) +
            geom_errorbar(aes(ymin = conf.low_mean, ymax = conf.high_mean, color = miss.prop), width = 0.08, position = position_dodge(width = 0.5)) +
            geom_point(aes(group = miss.prop, color = miss.prop), size = 3, shape = 19, position = position_dodge(width = 0.5)) + 
            scale_color_manual(values = color_map_pmiss)+
            labs(x = "Method", y = expression(beta), color=c("Missingness proportion (%)")) + 
            theme(axis.title.x = element_text(face = "bold", size = 14),
                  axis.title.y = element_text(face = "bold", size = 14),
                  axis.text.y = element_text(size = 12),
                  strip.placement = "outside",
                  panel.spacing = unit(1, "lines"),
                  axis.text.x = element_text(size = 12, colour = color_map)) +
            theme(legend.position = "none",
                  panel.grid = element_line(color = "gray", linetype = "dotted"),
                  panel.background = element_rect(fill = "white", color = "black"),  
                  panel.border = element_rect(color = "black", linewidth = 1, fill = NA),
                  strip.background = element_rect(fill = "grey", color = "black"),  
                  strip.text = element_text(color = "black", face = "bold", size = 14),
                  legend.key = element_rect(color = NA, fill = NA),
                  legend.title = element_text(size = 14, face = "bold"),   
                  legend.text = element_text(size = 12))  

p1b = ggplot(ds_plot_sum_final %>% filter(grepl(obsvar_pattern, coeff)), 
            aes(x = method_abbr, y = estimate_mean, group = miss.prop)) + 
            facet_wrap(~ coeff_label, nrow = 2, scales = "free_y") +
            geom_hline(aes(yintercept = true_value), linetype = "dashed", col = "grey", linewidth=1) +
            geom_errorbar(aes(ymin = conf.low_mean, ymax = conf.high_mean, color = miss.prop), width = 0.08, position = position_dodge(width = 0.5)) +
            geom_point(aes(group = miss.prop, color = miss.prop), size = 3, shape = 19, position = position_dodge(width = 0.5)) + 
            scale_color_manual(values = color_map_pmiss)+
            labs(x = "Method", y = expression(beta), color=c("Missingness proportion (%)")) + 
            theme(axis.title.x = element_text(face = "bold", size = 14),
                  axis.title.y = element_text(face = "bold", size = 14),
                  axis.text.y = element_text(size = 12),
                  strip.placement = "outside",
                  panel.spacing = unit(1, "lines"),
                  axis.text.x = element_text(size = 12, colour = color_map)) +
            theme(legend.position = "bottom",
                  panel.grid = element_line(color = "gray", linetype = "dotted"),
                  panel.background = element_rect(fill = "white", color = "black"),  
                  panel.border = element_rect(color = "black", linewidth = 1, fill = NA),
                  strip.background = element_rect(fill = "grey", color = "black"),  
                  strip.text = element_text(color = "black", face = "bold", size = 14),
                  legend.key = element_rect(color = NA, fill = NA),
                  legend.title = element_text(size = 14, face = "bold"),   
                  legend.text = element_text(size = 12))  

p1a / p1b + plot_annotation(tag_levels = 'a') 

## Figure S1: bias ----
p2a = ggplot(results.df_plot %>% filter(grepl(missvar_pattern, coeff) & method_abbr != "REF"), 
       aes(x = method_abbr, y = bias)) + 
  facet_wrap(~ coeff_label, nrow = 2, scales = "free_y") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth=1) +
  geom_boxplot(aes(group=method_pmiss, color=miss.prop)) +
  scale_color_manual(values = color_map_pmiss[3:5])+
  labs(x = "Method", y = "Bias", color=c("Missingness proportion (%)")) + 
  theme(axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.y = element_text(size = 12),
        strip.placement = "outside",
        panel.spacing = unit(1, "lines"),
        axis.text.x = element_text(size = 12, colour = color_map[-1])) +
  theme(legend.position = "none",
        panel.grid = element_line(color = "gray", linetype = "dotted"),
        panel.background = element_rect(fill = "white", color = "black"),  
        panel.border = element_rect(color = "black", linewidth = 1, fill = NA),
        strip.background = element_rect(fill = "grey", color = "black"),  
        strip.text = element_text(color = "black", face = "bold", size = 14),
        legend.key = element_rect(color = NA, fill = NA),
        legend.title = element_text(size = 14, face = "bold"),   
        legend.text = element_text(size = 12))  

p2b = ggplot(results.df_plot %>% filter(grepl(obsvar_pattern, coeff) & method_abbr != "REF"), 
       aes(x = method_abbr, y = bias)) + 
  facet_wrap(~ coeff_label, nrow = 2, scales = "free_y") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth=1) +
  geom_boxplot(aes(group=method_pmiss, color=miss.prop)) +
  scale_color_manual(values = color_map_pmiss)+
  labs(x = "Method", y = "Bias", color=c("Missingness proportion (%)")) + 
  theme(axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.y = element_text(size = 12),
        strip.placement = "outside",
        panel.spacing = unit(1, "lines"),
        axis.text.x = element_text(size = 12, colour = color_map[-1])) +
  theme(legend.position = "bottom",
        panel.grid = element_line(color = "gray", linetype = "dotted"),
        panel.background = element_rect(fill = "white", color = "black"),  
        panel.border = element_rect(color = "black", linewidth = 1, fill = NA),
        strip.background = element_rect(fill = "grey", color = "black"),  
        strip.text = element_text(color = "black", face = "bold", size = 14),
        legend.key = element_rect(color = NA, fill = NA),
        legend.title = element_text(size = 14, face = "bold"),   
        legend.text = element_text(size = 12))  

p2a / p2b + plot_annotation(tag_levels = 'a') 

## Figure S4: RMSE ----
ds_plot_RMSE <- results.df_plot %>% group_by(miss.prop, run, method) %>% 
  filter(method_abbr != "REF" & coeff == "Age" & miss.prop != 10) %>%
  mutate(model_RMSE = sqrt(model_MSE))
ds_plot_RMSE$method_abbr <- factor(ds_plot_RMSE$method_abbr, levels = c("CC", "MM", "MxI", "MF", "KNN", "PMM", "RF", "CART"),
                                   labels = c("CC", "MM", "MxI", "MF", "KNN", "PMM (m=5)", "RF (m=5)", "CART (m=5)"))

ref_RMSE <- unique(sqrt(results.df_plot %>% filter(method_abbr == "REF" & coeff == "Age") %>% pull(model_MSE)))

plot_RMSE = ggplot(ds_plot_RMSE, aes(x = method_abbr, y = model_RMSE, group=method_pmiss)) + 
  geom_hline(yintercept = ref_RMSE, linetype = "dashed", color = "grey", linewidth=1) +
  geom_boxplot(aes(color=miss.prop)) +
  scale_color_manual(values = color_map_pmiss)+
  labs(x = "Method", y = "RMSE", color=c("Missingness proportion (%)")) +
  coord_cartesian(ylim = c(min(ds_plot_RMSE$model_RMSE), max(ds_plot_RMSE$model_RMSE))) +
  theme(axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.y = element_text(size = 12),
        strip.placement = "outside",
        panel.spacing = unit(1, "lines"),
        axis.text.x = element_text(size = 12, colour = color_map[-1])) +
  theme(legend.position = "bottom",
        panel.grid = element_line(color = "gray", linetype = "dotted"),
        panel.background = element_rect(fill = "white", color = "black"),  
        panel.border = element_rect(color = "black", linewidth = 1, fill = NA),
        strip.background = element_rect(fill = "grey", color = "black"),  
        strip.text = element_text(color = "black", face = "bold", size = 14),
        legend.key = element_rect(color = NA, fill = NA),
        legend.title = element_text(size = 14, face = "bold"),   
        legend.text = element_text(size = 12))  


## Figure S5: AIC ----
ds_plot_AIC <- results.df_plot %>% group_by(miss.prop, run, method) %>% 
  filter(method_abbr != "REF" & coeff == "Age")
ds_plot_AIC$method_abbr <- factor(ds_plot_AIC$method_abbr, levels = c("CC", "MM", "MxI", "MF", "KNN", "PMM", "RF", "CART"),
                                  labels = c("CC", "MM", "MxI", "MF", "KNN", "PMM (m=5)", "RF (m=5)", "CART (m=5)"))

plot_AIC = ggplot(ds_plot_AIC, aes(x = method_abbr, y = model_AIC)) + 
  geom_hline(yintercept = 13344.06, linetype = "dashed", color = "grey", linewidth=1) +
  geom_boxplot(aes(group=method_pmiss, color=miss.prop)) +
  scale_color_manual(values = color_map_pmiss)+
  labs(x = "Method", y = "AIC", color=c("Missingness proportion (%)")) + 
  coord_cartesian(ylim = c(min(ds_plot_AIC$model_AIC), max(ds_plot_AIC$model_AIC))) +
  theme(axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.y = element_text(size = 12),
        strip.placement = "outside",
        panel.spacing = unit(1, "lines"),
        axis.text.x = element_text(size = 12, colour = color_map[-1])) +
  theme(legend.position = "bottom",
        panel.grid = element_line(color = "gray", linetype = "dotted"),
        panel.background = element_rect(fill = "white", color = "black"),  
        panel.border = element_rect(color = "black", linewidth = 1, fill = NA),
        strip.background = element_rect(fill = "grey", color = "black"),  
        strip.text = element_text(color = "black", face = "bold", size = 14),
        legend.key = element_rect(color = NA, fill = NA),
        legend.title = element_text(size = 14, face = "bold"),   
        legend.text = element_text(size = 12))  

## Figure 3: beta for MI (m = 5 to 20) ----
# p = 20%
load("7_pmiss0.2_simres_m5_2025-04-22.RData")
load("7_pmiss0.2_simres_m10_2025-04-22.RData")
load("7_pmiss0.2_simres_m15_2025-04-23.RData")
load("7_pmiss0.2_simres_m20_2025-04-23.RData")
results.df_m5 = bind_rows(results_0.2_no.MI_5, .id = "run")
results.df_m10 = bind_rows(results_0.2_no.MI10, .id = "run")
results.df_m15 = bind_rows(results_0.2_no.MI15, .id = "run")
results.df_m20 = bind_rows(results_0.2_no.MI20, .id = "run")

results.df_m = bind_rows(results.df_m5, results.df_m10, results.df_m15, results.df_m20)
results.df = results.df_m

results.df_plot <- merge(merge(results.df, ds_coeff, by = "coeff"),
                         ds_methods_all %>% rename("method_abbr" = "abbr"), by = "method") %>%
                   mutate(coeff_label = case_when(grepl("Time call", coeff_label) ~ sub("Time call", "Call Time", coeff_label),
                                                 T ~ coeff_label),
                          method_abbr = ifelse(method_abbr == "kNN", "KNN", method_abbr),
                          method_pmiss = paste0(method, pmiss),
                          miss.prop = ifelse(method_abbr == "REF", 0, pmiss*100))

results.df_plot$miss.prop <- factor(results.df_plot$miss.prop, levels = sort(unique(results.df_plot$miss.prop)))#c(0, 10, 20, 30, 40))
results.df_plot$method_abbr <- factor(results.df_plot$method_abbr, levels = c("REF", "CC", "MM", "MI", "MF", "KNN", "PMM", "RF", "CART"))
results.df_plot$m = factor(results.df_plot$m)

n.sim = 100
results.df_plot_sum <- results.df_plot %>% group_by(m, miss.prop, pmiss, coeff, coeff_label, method, method_abbr, method_group, true_value) %>% 
                                           mutate(run = as.numeric(run)) %>%
                                           summarize(estimate_mean = mean(estimate), 
                                                    conf.low_mean = mean(conf.low), conf.high_mean = mean(conf.high),
                                                    SE_mean = mean(std.error),
                                                    cov = sum(coverage)/n.sim*100,
                                                    estimate_var = var(estimate)) %>% ungroup()

ds_plot_sum_final <- results.df_plot_sum %>% mutate(estimate_SE = case_when(coeff == "Age" ~ paste0(sprintf("%.4f", estimate_mean), " (", round(SE_mean, 3), ")"),
                                                                            grepl("ime_call", coeff) ~ paste0(sprintf("%.3f", estimate_mean), " (", round(SE_mean, 3), ")"),
                                                                            T ~ paste0(sprintf("%.2f", estimate_mean), " (", round(SE_mean, 3), ")")),
                                                    bias = estimate_mean - true_value, 
                                                    MSE = bias^2 + estimate_var) %>%  
                                             mutate(coeff_label = case_when(coeff_label == "Witness Type (family)" ~ "Witness Type (Family)",
                                                   coeff_label == "Witness Type (healthcare provider)" ~ "Witness Type (Healthcare provider)",
                                                   coeff_label == "Witness Type (lay person)" ~ "Witness Type (Lay person)",
                                                   T ~ coeff_label)) %>%
                                            mutate(m = factor(m))



p1a = ggplot(ds_plot_sum_final %>% filter(grepl(missvar_pattern, coeff)), 
             aes(x = method_abbr, y = estimate_mean, group = m)) + 
  facet_wrap(~ coeff_label, nrow = 2, scales = "free_y") +
  geom_hline(aes(yintercept = true_value), linetype = "dashed", col = "grey", linewidth=1) +
  geom_errorbar(aes(ymin = conf.low_mean, ymax = conf.high_mean, color = m), width = 0.08, position = position_dodge(width = 0.5)) +
  geom_point(aes(group = m, color = m), size = 3, shape = 19, position = position_dodge(width = 0.5)) + 
  scale_color_manual(values = color_map_m)+
  labs(x = "Method", y = expression(beta), color=c("Number of multiple imputations (m)")) + 
  theme(axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.y = element_text(size = 12),
        strip.placement = "outside",
        panel.spacing = unit(1, "lines"),
        axis.text.x = element_text(size = 12, colour = "darkblue")) +
  theme(legend.position = "none",
        panel.grid = element_line(color = "gray", linetype = "dotted"),
        panel.background = element_rect(fill = "white", color = "black"),  
        panel.border = element_rect(color = "black", linewidth = 1, fill = NA),
        strip.background = element_rect(fill = "grey", color = "black"),  
        strip.text = element_text(color = "black", face = "bold", size = 14),
        legend.key = element_rect(color = NA, fill = NA),
        legend.title = element_text(size = 14, face = "bold"),   
        legend.text = element_text(size = 12))  

p1b = ggplot(ds_plot_sum_final %>% filter(grepl(obsvar_pattern, coeff)), 
             aes(x = method_abbr, y = estimate_mean, group = m)) + 
  facet_wrap(~ coeff_label, nrow = 2, scales = "free_y") +
  geom_hline(aes(yintercept = true_value), linetype = "dashed", col = "grey", linewidth=1) +
  geom_errorbar(aes(ymin = conf.low_mean, ymax = conf.high_mean, color = m), width = 0.08, position = position_dodge(width = 0.5)) +
  geom_point(aes(group = m, color = m), size = 3, shape = 19, position = position_dodge(width = 0.5)) + 
  scale_color_manual(values = color_map_m)+
  labs(x = "Method", y = expression(beta), color=c("Number of multiple imputations (m)")) + 
  theme(axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.y = element_text(size = 12),
        strip.placement = "outside",
        panel.spacing = unit(1, "lines"),
        axis.text.x = element_text(size = 12, colour = "darkblue")) +
  theme(legend.position = "bottom",
        panel.grid = element_line(color = "gray", linetype = "dotted"),
        panel.background = element_rect(fill = "white", color = "black"),  
        panel.border = element_rect(color = "black", linewidth = 1, fill = NA),
        strip.background = element_rect(fill = "grey", color = "black"),  
        strip.text = element_text(color = "black", face = "bold", size = 14),
        legend.key = element_rect(color = NA, fill = NA),
        legend.title = element_text(size = 14, face = "bold"),   
        legend.text = element_text(size = 12))  

p1a / p1b + 
  plot_annotation(tag_levels = 'a') 

## Figure S2: bias for MI (m = 5 to 20) ----
p2a = ggplot(results.df_plot %>% filter(pmiss == 0.2 & grepl(missvar_pattern, coeff) & method_abbr != "REF"), 
             aes(x = method_abbr, y = bias)) + 
  facet_wrap(~ coeff_label, nrow = 2, scales = "free_y") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth=1) +
  geom_boxplot(aes(group = interaction(method_pmiss, m), color=m)) +
  scale_color_manual(values = color_map_m)+
  labs(x = "Method", y = "Bias", color=c("Number of multiple imputations (m)")) + 
  theme(axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.y = element_text(size = 12),
        strip.placement = "outside",
        panel.spacing = unit(1, "lines"),
        axis.text.x = element_text(size = 12, colour = "darkblue")) +
  theme(legend.position = "none",
        panel.grid = element_line(color = "gray", linetype = "dotted"),
        panel.background = element_rect(fill = "white", color = "black"),  
        panel.border = element_rect(color = "black", linewidth = 1, fill = NA),
        strip.background = element_rect(fill = "grey", color = "black"),  
        strip.text = element_text(color = "black", face = "bold", size = 14),
        legend.key = element_rect(color = NA, fill = NA),
        legend.title = element_text(size = 14, face = "bold"),   
        legend.text = element_text(size = 12))  

p2b = ggplot(results.df_plot %>% filter(pmiss == 0.2) %>% filter(grepl(obsvar_pattern, coeff) & method_abbr != "REF"), 
             aes(x = method_abbr, y = bias)) + 
  facet_wrap(~ coeff_label, nrow = 2, scales = "free_y") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth=1) +
  geom_boxplot(aes(group=interaction(method_pmiss, m), color=m)) +
  scale_color_manual(values = color_map_m)+
  labs(x = "Method", y = "Bias", color=c("Number of multiple imputations (m)")) + 
  theme(axis.title.x = element_text(face = "bold", size = 14),
        axis.title.y = element_text(face = "bold", size = 14),
        axis.text.y = element_text(size = 12),
        strip.placement = "outside",
        panel.spacing = unit(1, "lines"),
        axis.text.x = element_text(size = 12, colour = "darkblue")) +
  theme(legend.position = "bottom",
        panel.grid = element_line(color = "gray", linetype = "dotted"),
        panel.background = element_rect(fill = "white", color = "black"),  
        panel.border = element_rect(color = "black", linewidth = 1, fill = NA),
        strip.background = element_rect(fill = "grey", color = "black"),  
        strip.text = element_text(color = "black", face = "bold", size = 14),
        legend.key = element_rect(color = NA, fill = NA),
        legend.title = element_text(size = 14, face = "bold"),   
        legend.text = element_text(size = 12))  

p2a / p2b + 
  plot_annotation(tag_levels = 'a') 






