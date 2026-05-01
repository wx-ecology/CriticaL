# ============================================================================
# 01_main_models.R
# ----------------------------------------------------------------------------
# Main analysis: testing whether settlement configuration (Settlement
# porosity) affects mammalian movement beyond what is explained by landscape
# composition.
#
# For each combination of:
#   - displacement time window: 1-day or 10-day
#   - displacement summary metric: median or 0.95 quantile (long-distance)
#   - composition variable: HFI, HMI, or settlement cover (log_bd)
#
# we fit three nested linear mixed-effects models:
#
#   (1) Baseline   : composition only (no configuration)
#   (2) Additive   : composition + Settlement porosity
#   (3) Interaction: composition * Settlement porosity
#
# All models share the same random-effect structure (nested taxonomic
# hierarchy) and exponential spatial-correlation structure (Haversine
# distance). Models are fit by maximum likelihood (ML) so that AIC-based
# comparisons are valid.
#
# Output:
#   results/models/Mods_globalbarrier_logdisp_pd_<time>d_<scale>.rds
#       List of fitted models for each (composition, configuration) combo.
#   results/tables/main_model_AIC_table.csv
#       AIC table comparing the three models per composition variable.
#   results/tables/main_model_fixed_effects.csv
#       Fixed-effect estimates from the best-supported model per combo.
#
# Run after: 00_setup.R (sourced automatically below)
# ============================================================================

source("./code/00_setup.R")


# ============================================================================
# Helper: fit the three nested models for a given composition variable
# ============================================================================

fit_composition_models <- function(dat, comp_var, optimizer = "optim") {

  # Build formula strings dynamically so the same code can handle HFI, HMI,
  # or log_bd as the composition variable.
  baseline_fml    <- as.formula(
    paste0("log_Displacement_km ~ log_BodyMass_kg + scale_NDVI + Diet + ",
           comp_var)
  )
  additive_fml    <- as.formula(
    paste0("log_Displacement_km ~ log_BodyMass_kg + scale_NDVI + Diet + ",
           comp_var, " + pd_adpt_km")
  )
  interaction_fml <- as.formula(
    paste0("log_Displacement_km ~ log_BodyMass_kg + scale_NDVI + Diet + ",
           comp_var, " * pd_adpt_km")
  )

  ctrl <- list(opt = optimizer, msMaxIter = 5000, msMaxEval = 5000,
               tolerance = 1e-6, niterEM = 50)
  rand <- ~ 1 | Order / Family / Genus / Species
  corr <- corHaversine(form = ~ lon + lat, mimic = "corExp")

  message("    fitting baseline (composition only)...")
  m_baseline <- lme(baseline_fml,    correlation = corr, random = rand,
                    control = ctrl, method = "ML", data = dat)

  message("    fitting additive (composition + porosity)...")
  m_additive <- lme(additive_fml,    correlation = corr, random = rand,
                    control = ctrl, method = "ML", data = dat)

  message("    fitting interaction (composition * porosity)...")
  m_inter    <- lme(interaction_fml, correlation = corr, random = rand,
                    control = ctrl, method = "ML", data = dat)

  tibble(
    comp_var = comp_var,
    conf_var = c(NA_character_, "pd_adpt_km", "pd_adpt_km*"),
    model    = list(m_baseline, m_additive, m_inter)
  )
}


# ============================================================================
# Run all model sets across displacement scales
# ============================================================================

for (time_scale_days in c(10, 1)) {
  for (spatial_scale in c("hi", "me")) {

    message("\n=== Fitting models: ", time_scale_days, "d_", spatial_scale, " ===")

    dat <- load_and_prep_data(time_scale_days, spatial_scale)

    # Some HMI / log_bd combinations did not converge with the default
    # `optim` optimizer and required `nlminb`. We use this empirically
    # determined optimizer choice for consistency across runs.
    models <- bind_rows(
      fit_composition_models(dat, "HFI",    optimizer = "optim"),
      fit_composition_models(dat, "HMI",    optimizer = "nlminb"),
      fit_composition_models(dat, "log_bd", optimizer = "nlminb")
    )

    out_file <- file.path(
      models_out_dir,
      paste0("Mods_globalbarrier_logdisp_pd_",
             time_scale_days, "d_", spatial_scale, ".rds")
    )
    saveRDS(models, out_file)
    message("  Saved: ", out_file)
  }
}


# ============================================================================
# Result extraction: AIC tables and fixed-effect summaries
# ============================================================================

extract_main_results <- function(time_scale_days, spatial_scale) {

  scale_label <- paste0(time_scale_days, "d_", spatial_scale)
  in_file     <- file.path(
    models_out_dir,
    paste0("Mods_globalbarrier_logdisp_pd_", scale_label, ".rds")
  )
  models <- readRDS(in_file)

  # AIC table: one row per fitted model
  aic_tbl <- models %>%
    rowwise() %>%
    mutate(
      AIC = AIC(model),
      df  = attr(logLik(model), "df")
    ) %>%
    ungroup() %>%
    group_by(comp_var) %>%
    mutate(deltaAIC = AIC - min(AIC)) %>%
    ungroup() %>%
    mutate(scale = scale_label) %>%
    dplyr::select(scale, comp_var, conf_var, df, AIC, deltaAIC)

  # Fixed-effect coefficients (all models, long format)
  fixef_tbl <- models %>%
    rowwise() %>%
    mutate(
      fixef = list(
        as.data.frame(summary(model)$tTable) %>%
          rownames_to_column("Parameter")
      )
    ) %>%
    ungroup() %>%
    dplyr::select(comp_var, conf_var, fixef) %>%
    unnest(fixef) %>%
    mutate(scale = scale_label) %>%
    dplyr::select(scale, comp_var, conf_var, Parameter,
                  Estimate = Value, SE = Std.Error,
                  DF, t_value = `t-value`, p_value = `p-value`)

  list(aic = aic_tbl, fixef = fixef_tbl)
}

aic_all   <- list()
fixef_all <- list()

for (time_scale_days in c(10, 1)) {
  for (spatial_scale in c("hi", "me")) {
    res <- extract_main_results(time_scale_days, spatial_scale)
    label <- paste0(time_scale_days, "d_", spatial_scale)
    aic_all[[label]]   <- res$aic
    fixef_all[[label]] <- res$fixef
  }
}

write_csv(bind_rows(aic_all),
          file.path(tables_out_dir, "main_model_AIC_table.csv"))
write_csv(bind_rows(fixef_all),
          file.path(tables_out_dir, "main_model_fixed_effects.csv"))

message("\nMain model AIC and fixed-effect tables written to ", tables_out_dir)
