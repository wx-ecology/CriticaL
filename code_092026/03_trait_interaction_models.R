# ============================================================================
# 03_trait_interaction_models.R
# ----------------------------------------------------------------------------
# Testing whether species traits (body mass and diet)
# modulate movement responses to (1) Settlement porosity (configuration) and
# (2) Human Footprint Index (composition).
#
# This addresses whether trait-based sensitivity is specific to landscape
# configuration or extends to composition as well, and helps identify which
# species traits predict greater vulnerability to settlement-driven barriers.
#
# Models follow the same framework as 01_main_models.R, with two
# preprocessing differences specific to the interaction analysis:
#   - log body mass is median-centered (rather than left uncentered) so the
#     main effects of porosity and HFI are interpreted at the median body
#     mass (~55 kg) of species in our dataset rather than at body mass = 1 kg.
#   - Diet is releveled with Herbivore as the reference category, so the
#     main effects of porosity and HFI represent the herbivore slope, which
#     is the most ecologically relevant reference for our analysis (large
#     herbivores show the strongest responses to settlement porosity).
#
# Four interaction models are fit per displacement scale (HFI is the only
# composition variable considered, for clarity):
#   (1) Diet      x Porosity  (configuration trait interaction)
#   (2) BodyMass  x Porosity  (configuration trait interaction)
#   (3) Diet      x HFI       (composition trait interaction)
#   (4) BodyMass  x HFI       (composition trait interaction)
#
# A follow-up herbivore-only model (BodyMass x Porosity) further tests
# whether body size mediates porosity sensitivity within the herbivore guild
# (the diet group showing the strongest configuration response). Order is
# dropped from the random effects since herbivores span fewer taxonomic
# orders.
#
# Diet-specific marginal slopes and pairwise contrasts are extracted from
# the Diet x Porosity and Diet x HFI models using `emtrends`.
#
# Output:
#   results/models/Mods_trait_interactions_<time>d_<scale>.rds
#       Tibble of fitted models (4 candidate structures per scale).
#   results/models/Mods_herbivore_bm_porosity_<time>d_<scale>.rds
#       Herbivore-only BodyMass x Porosity model per scale.
#   results/tables/trait_interaction_fixed_effects.csv
#   results/tables/herbivore_trait_fixed_effects.csv
#   results/tables/diet_marginal_slopes.csv
#   results/tables/diet_pairwise_comparisons.csv
# ============================================================================

source("./code/00_setup.R")

library(emmeans)   # emtrends() for diet-specific slopes


# ============================================================================
# Helper: preprocess data for trait interaction models
# ----------------------------------------------------------------------------
# Differs from `load_and_prep_data()` in two ways:
#   (i) log body mass is median-centered rather than left uncentered
#  (ii) Diet is releveled with Herbivore as the reference category
# ============================================================================

prep_trait_interaction_data <- function(dat) {
  dat %>%
    mutate(
      log_BodyMass_kg = log_BodyMass_kg -
        median(log_BodyMass_kg, na.rm = TRUE),
      Diet            = relevel(factor(Diet), ref = "Herbivore")
    )
}


# ============================================================================
# Helper: fit the four trait-interaction models
# ============================================================================

fit_trait_models <- function(dat) {
  
  ctrl <- list(opt = "optim", msMaxIter = 1000, msMaxEval = 1000)
  rand <- ~ 1 | Order / Family / Genus / Species
  corr <- corHaversine(form = ~ lon + lat, mimic = "corExp")
  
  message("    fitting Diet x Porosity...")
  m_diet_por <- lme(
    log_Displacement_km ~ log_BodyMass_kg + scale_NDVI +
      Diet * pd_adpt_km + HFI,
    correlation = corr, random = rand, control = ctrl,
    method = "ML", data = dat
  )
  
  message("    fitting BodyMass x Porosity...")
  m_bm_por <- lme(
    log_Displacement_km ~ log_BodyMass_kg * pd_adpt_km +
      scale_NDVI + Diet + HFI,
    correlation = corr, random = rand, control = ctrl,
    method = "ML", data = dat
  )
  
  message("    fitting Diet x HFI...")
  m_diet_hfi <- lme(
    log_Displacement_km ~ log_BodyMass_kg + scale_NDVI +
      Diet * HFI + pd_adpt_km,
    correlation = corr, random = rand, control = ctrl,
    method = "ML", data = dat
  )
  
  message("    fitting BodyMass x HFI...")
  m_bm_hfi <- lme(
    log_Displacement_km ~ log_BodyMass_kg * HFI +
      scale_NDVI + Diet + pd_adpt_km,
    correlation = corr, random = rand, control = ctrl,
    method = "ML", data = dat
  )
  
  tibble(
    model_label = c("Diet x Porosity", "BodyMass x Porosity",
                    "Diet x HFI",      "BodyMass x HFI"),
    model       = list(m_diet_por, m_bm_por, m_diet_hfi, m_bm_hfi)
  )
}


# ============================================================================
# Helper: fit herbivore-only BodyMass x Porosity model
# ----------------------------------------------------------------------------
# Order is dropped from random effects since herbivores span fewer
# taxonomic orders.
# ============================================================================

fit_herbivore_model <- function(dat) {
  
  dat_herb <- dat %>% filter(Diet == "Herbivore")
  
  lme(
    log_Displacement_km ~ log_BodyMass_kg * pd_adpt_km +
      scale_NDVI + HFI,
    correlation = corHaversine(form = ~ lon + lat, mimic = "corExp"),
    random      = ~ 1 | Family / Genus / Species,
    control     = list(opt = "optim", msMaxIter = 1000, msMaxEval = 1000),
    method      = "ML",
    data        = dat_herb
  )
}


# ============================================================================
# Run trait-interaction and herbivore models across all displacement scales
# ============================================================================

for (time_scale_days in c(10, 1)) {
  for (spatial_scale in c("hi", "me")) {
    
    scale_label <- paste0(time_scale_days, "d_", spatial_scale)
    message("\n=== Trait interaction models: ", scale_label, " ===")
    
    dat <- load_and_prep_data(time_scale_days, spatial_scale) %>%
      prep_trait_interaction_data()
    
    # Four main interaction models
    models <- fit_trait_models(dat)
    write_rds(models, file.path(
      models_out_dir,
      paste0("Mods_trait_interactions_", scale_label, ".rds")
    ))
    
    # Herbivore-only follow-up
    message("    fitting Herbivore-only BodyMass x Porosity...")
    m_herb <- fit_herbivore_model(dat)
    write_rds(m_herb, file.path(
      models_out_dir,
      paste0("Mods_herbivore_bm_porosity_", scale_label, ".rds")
    ))
    
    message("  Saved models for ", scale_label)
  }
}


# ============================================================================
# Result extraction: fixed-effect summaries
# ============================================================================

extract_trait_results <- function(time_scale_days, spatial_scale) {
  
  scale_label <- paste0(time_scale_days, "d_", spatial_scale)
  in_file     <- file.path(
    models_out_dir,
    paste0("Mods_trait_interactions_", scale_label, ".rds")
  )
  models <- read_rds(in_file)
  
  models %>%
    rowwise() %>%
    mutate(
      fixef = list(
        as.data.frame(summary(model)$tTable) %>%
          rownames_to_column("Parameter")
      )
    ) %>%
    ungroup() %>%
    dplyr::select(model_label, fixef) %>%
    unnest(fixef) %>%
    mutate(scale = scale_label) %>%
    dplyr::select(scale, model_label, Parameter,
                  Estimate = Value, SE = Std.Error,
                  DF, t_value = `t-value`, p_value = `p-value`)
}

extract_herbivore_results <- function(time_scale_days, spatial_scale) {
  
  scale_label <- paste0(time_scale_days, "d_", spatial_scale)
  in_file     <- file.path(
    models_out_dir,
    paste0("Mods_herbivore_bm_porosity_", scale_label, ".rds")
  )
  m_herb <- read_rds(in_file)
  
  as.data.frame(summary(m_herb)$tTable) %>%
    rownames_to_column("Parameter") %>%
    mutate(scale = scale_label) %>%
    dplyr::select(scale, Parameter,
                  Estimate = Value, SE = Std.Error,
                  DF, t_value = `t-value`, p_value = `p-value`)
}

fixef_all      <- list()
herb_fixef_all <- list()

for (time_scale_days in c(10, 1)) {
  for (spatial_scale in c("hi", "me")) {
    label <- paste0(time_scale_days, "d_", spatial_scale)
    fixef_all[[label]]      <- extract_trait_results(time_scale_days, spatial_scale)
    herb_fixef_all[[label]] <- extract_herbivore_results(time_scale_days, spatial_scale)
  }
}

write_csv(bind_rows(fixef_all),
          file.path(tables_out_dir, "trait_interaction_fixed_effects.csv"))
write_csv(bind_rows(herb_fixef_all),
          file.path(tables_out_dir, "herbivore_trait_fixed_effects.csv"))


# ============================================================================
# Diet-specific marginal slopes and pairwise contrasts (emtrends)
# ----------------------------------------------------------------------------
# Extracts per-diet slopes (and pairwise differences in slope) from the
# Diet x Porosity and Diet x HFI models. These quantify how strongly each
# diet group's movement responds to porosity / HFI.
# ============================================================================

extract_emtrends <- function(mod, predictor, model_label, scale_label) {
  
  slopes <- emtrends(mod, ~ Diet, var = predictor) %>%
    summary(infer = TRUE) %>%
    as.data.frame() %>%
    rename(
      Estimate = !!paste0(predictor, ".trend"),
      t_value  = t.ratio,
      p_value  = p.value
    ) %>%
    mutate(
      Predictor = predictor,
      Model     = model_label,
      Scale     = scale_label
    ) %>%
    dplyr::select(Scale, Model, Predictor, Diet,
                  Estimate, SE, DF, t_value, p_value)
  
  pairs <- emtrends(mod, pairwise ~ Diet, var = predictor)$contrasts %>%
    summary(infer = TRUE) %>%
    as.data.frame() %>%
    rename(
      Contrast = contrast,
      Estimate = estimate,
      t_value  = t.ratio,
      p_value  = p.value
    ) %>%
    mutate(
      Predictor = predictor,
      Model     = model_label,
      Scale     = scale_label
    ) %>%
    dplyr::select(Scale, Model, Predictor, Contrast,
                  Estimate, SE, DF, t_value, p_value)
  
  list(slopes = slopes, pairs = pairs)
}

diet_slopes_all <- list()
diet_pairs_all  <- list()

for (time_scale_days in c(10, 1)) {
  for (spatial_scale in c("hi", "me")) {
    
    scale_label <- paste0(time_scale_days, "d_", spatial_scale)
    
    models <- read_rds(file.path(
      models_out_dir,
      paste0("Mods_trait_interactions_", scale_label, ".rds")
    ))
    
    m_diet_por <- models$model[[which(models$model_label == "Diet x Porosity")]]
    m_diet_hfi <- models$model[[which(models$model_label == "Diet x HFI")]]
    
    em_por <- extract_emtrends(m_diet_por, "pd_adpt_km",
                               "Diet x Porosity", scale_label)
    em_hfi <- extract_emtrends(m_diet_hfi, "HFI",
                               "Diet x HFI", scale_label)
    
    diet_slopes_all[[scale_label]] <- bind_rows(em_por$slopes, em_hfi$slopes)
    diet_pairs_all[[scale_label]]  <- bind_rows(em_por$pairs,  em_hfi$pairs)
  }
}

write_csv(bind_rows(diet_slopes_all),
          file.path(tables_out_dir, "diet_marginal_slopes.csv"))
write_csv(bind_rows(diet_pairs_all),
          file.path(tables_out_dir, "diet_pairwise_comparisons.csv"))

message("\nTrait interaction outputs written to ", tables_out_dir)