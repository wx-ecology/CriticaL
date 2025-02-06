# visualize modeling results 

library(ggplot2)
library(tidyverse)
library(cowplot)

##############################################################
## ----- FIG 1 : percolation distance vs displacement --------- 
##############################################################
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

# ggsave("./visualization/FIGURES/FIG_10d_displacementVSpercolation.png")

move.sum.1.hi <- read_rds(paste0("./data/movement/ready_data/move_1d_hi.rds"))
move.sum.1.me <- read_rds(paste0("./data/movement/ready_data/move_1d_me.rds"))

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
# ggsave("./visualization/FIGURES/FIG_1d_displacementVSpercolation.png")

##############################################################
## ----- FIG : species adapt pd figure  --------------------- 
##############################################################
move.sum.10.hi <- read_rds(paste0("./data/movement/ready_data/move_10d_hi.rds"))
summary <- move.sum.10.hi %>% group_by(Binomial, Diet) %>%
  summarise(mean_disp = mean(Displacement_km, na.rm = T),
            max_disp = quantile(Displacement_km, 0.95, na.rm = T),
            min_disp = quantile(Displacement_km, 0.05, na.rm = T))

move.sum.10.hi$Binomial <- factor(move.sum.10.hi$Binomial, 
                                  levels = summary$Binomial[order(summary$mean_disp)])

disp <- ggplot(move.sum.10.hi, aes(x = Binomial, y = Displacement_km, color = as.factor(Diet))) +
  geom_pointrange(
    data = summary,
    aes(y = mean_disp, ymin = min_disp, ymax = max_disp),
    size = 1.2,  # Adjust line thickness
    position = position_dodge(width = 0.4)
  ) +
  geom_jitter(
    width = 0.2,  # Adjust jitter width
    size = 2,  # Point size
    alpha = 0.3  # Transparency
  ) +
  scale_color_manual(values = c("#ba8caa", "#93b7be", "#facb0f")) +
  coord_flip() +
  theme_minimal() +
  theme(
    legend.position = "none"
  )

move.sum.10.hi <- read_rds(paste0("./data/movement/ready_data/move_10d_hi.rds"))
summary <- move.sum.10.hi %>% group_by(Binomial, Diet) %>%
  summarise(mean_pd = mean(pd_adpt_km, na.rm = T),
            max_pd = quantile(pd_adpt_km, 0.95, na.rm = T),
            min_pd = quantile(pd_adpt_km, 0.05, na.rm = T))

move.sum.10.hi$Binomial <- factor(move.sum.10.hi$Binomial, 
                                 levels = summary$Binomial[order(summary$mean_pd)])

pd <- ggplot(move.sum.10.hi, aes(x = Binomial, y = pd_adpt_km, color = as.factor(Diet))) +
  geom_pointrange(
    data = summary,
    aes(y = mean_pd, ymin = min_pd, ymax = max_pd),
    size = 1.2,  # Adjust line thickness
    position = position_dodge(width = 0.4)
  ) +
  geom_jitter(
    width = 0.2,  # Adjust jitter width
    size = 2,  # Point size
    alpha = 0.3  # Transparency
  ) +
  scale_color_manual(values = c("#ba8caa", "#93b7be", "#facb0f")) +
  coord_flip() +
  theme_minimal() +
  theme(
    legend.position = "none"
  )

plot_grid(disp, pd, labels = c('A', 'B'), label_size = 12)
ggsave("./visualization/FIGURES/FIG_disp_pd_by_spp.png")

##############################################################
# ------ TABLE: model comparisons --------- 
##############################################################
## igure does not make sense here and the figure code below is just attemps
library(stringr)
mod_comp <- read_csv("./results/ModResults_allSpp/scale_zscored/ModResults_allSpp.csv")

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

##########################################
## ----- FIG: model coefficients (all spp) ---------
##########################################

folder_path <- "./results/models/zscored"  
models <- list.files(folder_path, pattern = "allspp.*pd|pd.*allspp", full.names = TRUE)
fixed = tibble()
for (model in models ){
  m = read_rds(model)
  name <-  str_extract(model, "(?<=allspp_).*?(?=.rds)")
  fixed.m <-  data.frame(intervals(m, which = "fixed")$fixed) %>% 
    bind_cols(name = name) %>%
    rownames_to_column(var = "variable") %>%
    filter(variable != "(Intercept)")
  fixed = rbind(fixed, fixed.m)
}

fixed <- fixed %>% filter(!str_detect(name, "slog_bd_w_pd")) %>%
  mutate(
    temporal_scale = case_when(str_sub(name, 1,3) == "10d" ~ "10-day",
                               str_sub(name, 1,3) == "1d_" ~ "1-day"),
    spatial_scale = case_when(str_sub(name, 4,6) == "_hi" | str_sub(name, 4,6) == "hi_" ~ "long-distance",
                              str_sub(name, 4,6) == "_me" | str_sub(name, 4,6) == "me_" ~ "median"),
    amount_variable = case_when(str_sub(name, 14,15) %in% c("d_", "_w") ~ "building density",
                                str_sub(name, 14,15) %in% c("HM", "MI") ~ "HMI",
                                str_sub(name, 14,15) %in% c("HF", "FI") ~ "HFI")  
    ) %>%
  dplyr::select(-name) 

p_HFI <- fixed %>% filter(amount_variable %in% c("HFI"),
                 !variable %in% c("DietHerbivore", "DietOmnivore")) %>%
  ggplot(aes(x = variable, y = est., color = spatial_scale)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 0.8) +
  geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.4)) +
  scale_color_manual(values = c("#B66DDD", "#FFA500")) +
  theme_minimal() +
  labs( x = " ",
        y = " ",
        title = "A. HFI model") +
  coord_flip() +
  facet_wrap(vars(temporal_scale), nrow = 1)
p_HMI <- fixed %>% filter(amount_variable %in% c("HMI"),
                          !variable %in% c("DietHerbivore", "DietOmnivore")) %>%
  ggplot(aes(x = variable, y = est., color = spatial_scale)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 0.8) +
  geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.4)) +
  scale_color_manual(values = c("#B66DDD", "#FFA500")) +
  theme_minimal() +
  labs( x = " ",
        y = " ",
        title = "B. HMI model") +
  coord_flip() +
  facet_wrap(vars(temporal_scale), nrow = 1)
p_bd <- fixed %>% filter(amount_variable %in% c("building density"),
                         !variable %in% c("DietHerbivore", "DietOmnivore")) %>%
  ggplot(aes(x = variable, y = est., color = spatial_scale)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey", linewidth = 0.8) +
  geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.4)) +
  scale_color_manual(values = c("#B66DDD", "#FFA500")) +
  theme_minimal() +
  labs( x = " ",
        y = "Coefficient (95 % CI)",
        title = "C. Building density model") +
  coord_flip() +
  facet_wrap(vars(temporal_scale), nrow = 1)

plot_grid(p_HFI, p_HMI, p_bd, label_size = 12, nrow = 3)
ggsave("./visualization/FIGURES/FIG_allspp_allmod_estimates_updated.png", width = 6.7, height = 10.4)
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

fixed %>% filter(amount_variable %in% c("HFI", "building density"),
                 variable %in% c("log_pd_spp", "log_pd_ind")) %>%
  mutate(
    amount_variable = case_when(amount_variable == "HFI" ~ "Human Footprint Index", .default = "Building Density"),
    amount_variable = factor(amount_variable, levels = c("Human Footprint Index", "Building Density")),
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
ggsave("./visualization/FIGURES/FIG_sppVSind_scaledHFI.png")

# make a confidence interval plot for all variables 

## ----- FIG: model coefficients (ind VS spp) --------- ## 
