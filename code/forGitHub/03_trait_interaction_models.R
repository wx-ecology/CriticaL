# ============================================================================
# 03_trait_interaction_models.R
# ----------------------------------------------------------------------------
# Supplementary analysis: testing whether species traits (body mass and diet)
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
# Output:
#   results/models/Mods_trait_interactions_<time>d_<scale>.rds
#       Tibble of fitted models (4 candidate structures per scale).
#   results/tables/trait_interaction_fixed_effects.csv
# ============================================================================

source("./code/00_setup.R")


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
# Run trait-interaction models across all displacement scales
# ============================================================================

for (time_scale_days in c(10, 1)) {
  for (spatial_scale in c("hi", "me")) {

    message("\n=== Trait interaction models: ", time_scale_days, "d_", spatial_scale, " ===")

    dat <- load_and_prep_data(time_scale_days, spatial_scale) %>%
      prep_trait_interaction_data()

    models <- fit_trait_models(dat)

    out_file <- file.path(
      models_out_dir,
      paste0("Mods_trait_interactions_", time_scale_days, "d_", spatial_scale, ".rds")
    )
    write_rds(models, out_file)
    message("  Saved: ", out_file)
  }
}


# ============================================================================
# Result extraction: AIC tables and fixed-effect summaries
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

fixef_all <- list()

for (time_scale_days in c(10, 1)) {
  for (spatial_scale in c("hi", "me")) {
    label <- paste0(time_scale_days, "d_", spatial_scale)
    fixef_all[[label]] <- extract_trait_results(time_scale_days, spatial_scale)
  }
}

write_csv(bind_rows(fixef_all),
          file.path(tables_out_dir, "trait_interaction_fixed_effects.csv"))

message("\nTrait interaction fixed-effect table written to ", tables_out_dir)
