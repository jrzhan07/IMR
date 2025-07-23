#-----------------------------------------------------------------------------#
# File name: 2_TablesFigures.R
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
load("7_pmiss0.1_simres_2024-11-14.RData")
load("7_pmiss0.2_simres_2024-11-14.RData")
load("7_pmiss0.3_simres_2024-11-14.RData")
load("7_pmiss0.4_simres_2024-11-14.RData")
load("7_10to40_simres_2024-11-14.RData")

my.render.cont <- c(.="Mean (SD)", .="Median [Min, Max]", .="Median (Q1 - Q3)")
my.render.cat <- "FREQ (PCTnoNA%)"

pvalue <- function(x, ...) {
  y <- unlist(x)
  g <- factor(rep(1:length(x), times=sapply(x, length)))
  if (is.numeric(y)) {
    # p <- t.test(y ~ g)$p.value
    # p <- t.test(y ~ g, paired=T)$p.value
    p <- wilcox.test(y ~ g)$p.value
    # p <- kruskal.test(y ~ g)$p.value
    
    c("", sub("<", "<", format.pval(p, digits=3, eps=0.001)))
  } else {
    expected <- chisq.test(table(y, g))$expected
    if(all(expected > 5)){
      p <- chisq.test(table(y, g))$p.value
      
      c("", sub("<", "<", format.pval(p, digits=3, eps=0.001)))
    }else{
      p <- fisher.test(table(y, g), workspace = 2e8)$p.value
      
      c("", paste0(sub("<", "<", format.pval(p, digits=3, eps=0.001)), "*"))
    }
  }
}

col1 <- c("darkgrey", rep("black", 5)) 
color_map <- setNames(col1, ds_methods$abbr)
color_map_pmiss <- setNames(c("darkgrey", "#0072B2", "#E69F00", "#2ca02c", "#8E44AD"), seq(0,40, by=10))

# [1] Performance measure ----
s <- simsum(data = results.df, estvarname = "estimate", true = "true_value", se = "std.error", 
            methodvar = "method", ref = "No data missing (reference)", by=c("pmiss", "coeff"))
summary_s <- s$summ


# [2] Outputs ----

## [2.1] Table1 ----
# Note: Tables 1 require original data 
# table1(~ . | Bystander_CPR + Alert_issued, data=ds_OHCA_final[,-c(1,2)], render.continuous = my.render.cont, render.categorical=my.render.cat)
# table1(~ . | Alert_issued + Bystander_CPR, data=ds_OHCA_final[,-c(1,2)], render.continuous = my.render.cont, render.categorical=my.render.cat)

## [2.2] Table2 ----
n.sim = 1000
results.df_plot <- merge(merge(results.df, ds_coeff, by = "coeff"),
                         ds_methods %>% rename("method_abbr" = "abbr"), by = "method") %>%
                   mutate(coeff_label = case_when(grepl("Time call", coeff_label) ~ sub("Time call", "Call Time", coeff_label),
                                                  T ~ coeff_label),
                          method_abbr = ifelse(method_abbr == "kNN", "KNN", method_abbr),
                          method_abbr = ifelse(method_abbr == "MI", "MxI", method_abbr),
                          method_pmiss = paste0(method, pmiss),
                          miss.prop = ifelse(method_abbr == "REF", 0, pmiss*100))
results.df_plot$miss.prop <- factor(results.df_plot$miss.prop, levels = c(0, 10, 20, 30, 40))
results.df_plot$method_abbr <- factor(results.df_plot$method_abbr, levels = c("REF", "CC", "MM", "MxI", "MF", "KNN"))

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
                                                   MSE = bias^2 + estimate_var) 

ds_table2 <- ds_plot_sum_final %>% filter(pmiss == 0.2) %>% select(coeff_label, method, estimate_SE) %>%
                                    pivot_wider(id_cols = coeff_label, names_from = "method", values_from = "estimate_SE") %>%
                                    select(coeff_label, `No data missing (reference)`, `Complete case`, `mean/mode`, `M.ind (mean/mode)`, `Missforest`, `kNN`) %>%
                                    rename("Complete case (CC)" = "Complete case",
                                           "mean/mode (MM)" = "mean/mode",
                                           "Missing-indicator (MI)" = "M.ind (mean/mode)",
                                           "missForest (MF)" = "Missforest",
                                           "K-Nearest Neighbour (KNN)" = "kNN") %>%
                                    mutate(coeff_label = case_when(coeff_label == "Gender (Male)" ~ "Male",
                                                                   coeff_label == "Alert Issued (Yes)" ~ "Alert issued",
                                                                   coeff_label == "Witness Type (family)" ~ "Bystander - family",
                                                                   coeff_label == "Witness Type (healthcare provider)" ~ "Bystander - healthcare provider",
                                                                   coeff_label == "Witness Type (lay person)" ~ "Bystander - lay person",
                                                                   coeff_label == "Arrest Location (Public)" ~ "Public arrest location",
                                                                   coeff_label == "First Rhythm (Shockable)" ~ "Shockable",
                                                                   coeff_label == "First Rhythm (Unknown)" ~ "Unknown",
                                                                   coeff_label == "Call Time (06:00–18:59)" ~ "06:00–18:59",
                                                                   coeff_label == "Call Time (19:00–23:59)" ~ "19:00–23:59",
                                                                   T ~ coeff_label)) %>%
                                    mutate(coeff_label = factor(coeff_label, levels = c("Alert issued", "Male", "Age", 
                                                                                        "Bystander - family", "Bystander - healthcare provider", "Bystander - lay person",
                                                                                        "06:00–18:59", "19:00–23:59", "Public arrest location", "Shockable", "Unknown"))) %>%
                                    arrange(coeff_label)

## [2.3] Figures ----
## Figure 1: beta ----
ds_plot_sum_final <- ds_plot_sum_final %>% filter(pmiss != 0.1) %>% 
                                           mutate(coeff_label = case_when(coeff_label == "Witness Type (family)" ~ "Witness Type (Family)",
                                                                          coeff_label == "Witness Type (healthcare provider)" ~ "Witness Type (Healthcare provider)",
                                                                          coeff_label == "Witness Type (lay person)" ~ "Witness Type (Lay person)",
                                                                          T ~ coeff_label))

missvar_pattern <- "Age|Arrest_witnessed_by|Time_call"
obsvar_pattern <- "Alert|Gender|Arrest_witnessedYes|Location_of_arrest|First_rhythm"


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

p1 = p1a / p1b + 
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(size = 16, face = 'bold'))


## Figure 2: bias ----
p2a = ggplot(results.df_plot %>% filter(pmiss != 0.1 & grepl(missvar_pattern, coeff) & method_abbr != "REF"), 
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

p2b = ggplot(results.df_plot %>% filter(pmiss != 0.1) %>% filter(grepl(obsvar_pattern, coeff) & method_abbr != "REF"), 
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

p2 = p2a / p2b + 
  plot_annotation(tag_levels = 'a') &
  theme(plot.tag = element_text(size = 16, face = 'bold'))

## Figure S1: RMSE ----
ds_plot_RMSE <- results.df_plot %>% group_by(run, method) %>% 
                filter(method_abbr != "REF" & coeff == "Age") %>%
                mutate(model_RMSE = sqrt(model_MSE))
ref_RMSE <- unique(sqrt(results.df_plot %>% filter(method_abbr == "REF" & coeff == "Age") %>% pull(model_MSE)))

pS1 = ggplot(ds_plot_RMSE, aes(x = method_abbr, y = model_RMSE)) + 
  geom_hline(yintercept = ref_RMSE, linetype = "dashed", color = "grey", linewidth=1) +
  geom_boxplot(aes(group=method_pmiss, color=miss.prop)) +
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

## Figure S2: AIC ----
ds_plot_AIC <- results.df_plot %>% group_by(run, method) %>% filter(method_abbr != "REF" & coeff == "Age")

pS2 = ggplot(ds_plot_AIC, aes(x = method_abbr, y = model_AIC)) + 
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
