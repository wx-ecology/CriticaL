# visualize modeling results 

library(ggplot2)
library(tidyverse)
library(cowplot)

## ----- FIG 1 : percolation distance vs displacement --------- ## 
move.sum.10.hi <- read_rds(paste0("./data/movement/ready_data/move_10d_hi.rds"))
move.sum.10.me <- read_rds(paste0("./data/movement/ready_data/move_10d_me.rds"))

p1 <- 
  move.sum.10.me %>% ggplot( aes(x = log_pd, y = log_Displacement_km)) +
  geom_point(color = "#B66DDD", alpha = 0.3) +
  geom_smooth(color = "#8B5BA5") +
  theme_minimal() +
  ylim(-5, 5) +
  xlim(-3.5,3.5) +
  xlab ("percolation distance (km, log scale)") +
  ylab ("10-day median displacement (km, log scale)")

p2 <- 
  move.sum.10.hi %>% ggplot( aes(x = log_pd, y = log_Displacement_km)) +
  geom_point(color = "orange", alpha = 0.3) +
  geom_smooth(color = "#ed883b") +
  theme_minimal() +
  ylim(-5, 5) +
  xlim(-3.5,3.5) +
  xlab ("percolation distance (km, log scale)") +
  ylab ("10-day long-distance displacement (km, log scale)")

plot_grid(p1, p2, labels = c('A', 'B'), label_size = 12)

move.sum.1.hi <- read_rds(paste0("./data/movement/ready_data/move_1d_hi.rds"))
move.sum.1.me <- read_rds(paste0("./data/movement/ready_data/move_1d_me.rds"))

ggsave("./visualization/FIGURES/FIG_10d_displacementVSpercolation.png")

p3 <- 
  move.sum.1.me %>% ggplot( aes(x = log_pd, y = log_Displacement_km)) +
  geom_point(color = "#B66DDD", alpha = 0.3) +
  geom_smooth(color = "#8B5BA5") +
  theme_minimal() +
  xlab ("percolation distance (km, log scale)") +
  ylab ("1-day median displacement (km, log scale)")

p4 <- 
  move.sum.1.hi %>% ggplot( aes(x = log_pd, y = log_Displacement_km)) +
  geom_point(color = "orange", alpha = 0.3) +
  geom_smooth(color = "#ed883b") +
  theme_minimal() +
  xlab ("percolation distance (km, log scale)") +
  ylab ("1-day long-distance displacement (km, log scale)")

plot_grid(p3, p4, labels = c('A', 'B'), label_size = 12)
ggsave("./visualization/FIGURES/FIG_1d_displacementVSpercolation.png")


# ------ TABLE: model comparisons --------- ## #figure does not make sense here
library(stringr)
mod_comp <- read_csv("./results/ModResults_allSpp.csv")

# # Bar plot
# dat.HFI <- 
#   mod_comp %>% filter(`disturbance amount variable` == "log_HFI") %>%
#   mutate(group = paste0(`temporal scale`, "d_", str_sub(`spatial scale`, 6,7))) %>%
#   dplyr::select(group, `configuration variable`, AIC) %>%
#   pivot_wider(names_from = `configuration variable`, values_from = AIC) %>%
#   rowwise() %>% # Apply operations row-wise
#   mutate(max_AIC = max(c(null, log_pd)),
#          across(null:log_pd, ~ . - max_AIC)) %>%
#   ungroup() %>% # Remove rowwise grouping
#   dplyr::select(-max_AIC) %>%
#   pivot_longer(cols = 2:3, names_to = "configuration measures", values_to = "delta_AIC") %>%
#   mutate(`configuration measures` = case_when( `configuration measures` == "log_pd" ~ "percolation distance",
#                                                .default = "null"))
# dat.HFI$`configuration measures` = factor(dat.HFI$`configuration measures`, 
#                                           levels = c("percolation distance", "null"))
# 
# ggplot(dat.HFI, aes(x = group, y = delta_AIC, fill = `configuration measures`)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   scale_y_reverse() + # Reverse for lower AIC = better
#   labs(
#     x = " ",
#     y = "delta AIC (lower is better)"
#   ) +
#   theme_minimal() +
#   scale_fill_manual(values=c("#ED883B", "#C6C6C6")) +
#   theme(legend.position = "bottom")
# 
# dat.bd <- 
#   mod_comp %>% filter(`disturbance amount variable` == "log_bd") %>%
#   mutate(group = paste0(`temporal scale`, "d_", str_sub(`spatial scale`, 6,7))) %>%
#   dplyr::select(group, `configuration variable`, AIC) %>%
#   pivot_wider(names_from = `configuration variable`, values_from = AIC) %>%
#   rowwise() %>% # Apply operations row-wise
#   mutate(max_AIC = max(c(null, log_dist_2_build, log_pd, `log_dist_2_build+log_pd`)),
#          across(null:`log_dist_2_build+log_pd`, ~ . - max_AIC)) %>%
#   ungroup() %>% # Remove rowwise grouping
#   dplyr::select(-max_AIC) %>%
#   pivot_longer(cols = 2:5, names_to = "configuration measures", values_to = "delta_AIC") %>%
#   dplyr::select(group, `configuration measures`, delta_AIC) %>%
#   mutate(`configuration measures` = case_when( `configuration measures` == "log_dist_2_build" ~ "distance to building",
#                                                `configuration measures` == "log_pd" ~ "percolation distance",
#                                                `configuration measures` == "log_dist_2_build+log_pd" ~ "distance to building \n  + \n percolation distance",
#                             .default = "null"))
# 
# dat.bd$`configuration measures` <- factor(dat.bd$`configuration measures`, levels = c("distance to building \n  + \n percolation distance",
#                                           "percolation distance","distance to building",  "null"))
# 
# ggplot(dat.bd, aes(x = group, y = delta_AIC, fill = `configuration measures`)) +
#   geom_bar(stat = "identity", position = "dodge") +
#   scale_y_reverse() + # Reverse for lower AIC = better
#   labs(
#     x = " ",
#     y = "delta AIC (lower is better)"
#   ) +
#   theme_minimal() +
#   scale_fill_manual(values=c("#ED883B", "#C6C6C6", "#E1E0E0", "#DADADA")) +
#   theme(legend.position = "bottom")
# 
# ggsave("./visualization/FIGURES/FIG_supp_model_comps_building_density.png", width = 8.5, height = 3.6)

## ----- FIG: model coefficients (all spp) --------- ## 

folder_path <- "./results/models"  
models <- list.files(folder_path, pattern = "allspp.*pd|pd.*allspp", full.names = TRUE)
for (model in models ){
  m = read_rds(model)
  pd = intervals(m, which = "fixed")$fixed["log_pd",]
  ndvi = intervals(m, which = "fixed")$fixed["scale_NDVI",]
  bodymass = intervals(m, which = "fixed")$fixed["log_BodyMass_kg",]
}

##########################################
# summarize spp vs ind results --------
##########################################
models <- list.files(folder_path, pattern = "sppVSind", full.names = TRUE)
fixed = tibble()
for (model in models ) {
  m <-  read_rds(model)
  name <-  str_extract(model, "(?<=sppVSind_).*?(?=.rds)")
  fixed.m <-  data.frame(intervals(m, which = "fixed")$fixed) %>% 
    bind_cols(name = name) %>%
    rownames_to_column(var = "variable") %>%
    filter(variable != "(Intercept)")
  fixed = rbind(fixed, fixed.m)
}

fixed <- fixed %>% mutate(
  temporal_scale = case_when(str_sub(name, 1,3) == "10d" ~ "10-day",
                             str_sub(name, 1,3) == "1d_" ~ "1-day"),
  spatial_scale = case_when(str_sub(name, 4,6) == "_hi" | str_sub(name, 4,6) == "hi_" ~ "long-distance",
                            str_sub(name, 4,6) == "_me" | str_sub(name, 4,6) == "me_" ~ "median"),
  amount_variable = case_when(str_sub(name, -3,-1) == "_bd" ~ "building density",
                              str_sub(name, -3,-1) == "HFI" ~ "HFI",
                              str_sub(name, -3,-1) == "HMI" ~ "HMI")
  ) %>%
  dplyr::select(-name) 

fixed %>% filter(amount_variable %in% c("HMI", "building density"),
                 variable %in% c("log_pd_spp", "log_pd_ind")) %>%
  mutate(
    amount_variable = case_when(amount_variable == "HMI" ~ "Human Modification Index", .default = "Building Density"),
    amount_variable = factor(amount_variable, levels = c("Human Modification Index", "Building Density")),
    color = factor(paste0(variable, spatial_scale), levels = 
                          c("log_pd_indlong-distance", "log_pd_indmedian",
                            "log_pd_spplong-distance", "log_pd_sppmedian"))) %>%
  ggplot(aes(x = temporal_scale, y = est., color = color)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 0.8) +
  geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.4)) +
  scale_color_manual(values = c("#b126fc", "#e9bfff", "#FFA500", "#ffe1a8")) +
  theme_minimal() +
  labs( x = "temporal scale",
        y = "Coefficient (95 % CI)") +
  facet_wrap(vars(amount_variable), nrow = 1)
ggsave("./visualization/FIGURES/FIG_sppVSind.png")

# make a confidence interval plot for all variables 

## ----- FIG: model coefficients (ind VS spp) --------- ## 
