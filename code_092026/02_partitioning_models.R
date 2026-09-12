# ============================================================================
# 02_partitioning_models.R
# ----------------------------------------------------------------------------
# Distinguishing behavioral plasticity from environmental filtering as
# alternative (but non-mutually-exclusive) mechanisms underlying the link
# between Settlement porosity and animal movement.
#
# Settlement porosity is partitioned into two within-species-centered
# components:
#   - pd_km_spp: species-mean porosity, capturing where each species occurs
#                (environmental filtering signal)
#   - pd_km_ind: individual deviation from the species mean, capturing
#                within-species movement responses to local porosity
#                conditions (behavioral plasticity signal)
#
# Both components are standardized (z-scored) so their effect sizes can be
# directly compared.
#
# Four candidate models are fit per displacement scale, all using HFI as the
# composition variable:
#   (1) none           : both components additive, no interaction
#   (2) inter-spp      : HFI * pd_km_spp + pd_km_ind
#   (3) inter-ind      : pd_km_spp + HFI * pd_km_ind
#   (4) inter-spp-ind  : HFI * pd_km_spp + HFI * pd_km_ind
#
# Output:
#   results/models/Mods_sppVSind_<time>d_<scale>.rds
#       Tibble of fitted models (4 candidate structures per scale).
#   results/tables/partitioning_AIC_table.csv
#   results/tables/partitioning_fixed_effects.csv
# ============================================================================

source("./code/00_setup.R")


# ============================================================================
# Helper: prepare data with within-species partitioning of porosity
# ============================================================================

prep_partitioned_data <- function(dat) {
  dat %>%
    group_by(Binomial) %>%
    mutate(
      pd_km_spp = mean(pd_adpt_km, na.rm = TRUE),
      pd_km_ind = pd_adpt_km - pd_km_spp
    ) %>%
    ungroup() %>%
    mutate(
      # Standardize so that effect sizes of pd_km_spp and pd_km_ind are
      # directly comparable in the model output.
      pd_km_spp = as.numeric(scale(pd_km_spp, center = TRUE, scale = TRUE)),
      pd_km_ind = as.numeric(scale(pd_km_ind, center = TRUE, scale = TRUE))
    )
}


# ============================================================================
# Helper: fit the four candidate partitioning models (HFI-based)
# ============================================================================

fit_partitioning_models <- function(dat) {

  # Model 1: additive (no interactions)
  fml_none <- log_Displacement_km ~ log_BodyMass_kg + scale_NDVI + Diet +
                                    HFI + pd_km_spp + pd_km_ind
  # Model 2: HFI x species-mean porosity
  fml_spp  <- log_Displacement_km ~ log_BodyMass_kg + scale_NDVI + Diet +
                                    HFI * pd_km_spp + pd_km_ind
  # Model 3: HFI x within-individual porosity
  fml_ind  <- log_Displacement_km ~ log_BodyMass_kg + scale_NDVI + Diet +
                                    pd_km_spp + HFI * pd_km_ind
  # Model 4: both interactions
  fml_both <- log_Displacement_km ~ log_BodyMass_kg + scale_NDVI + Diet +
                                    HFI * pd_km_spp + HFI * pd_km_ind

  ctrl <- list(opt = "optim", msMaxIter = 1000, msMaxEval = 1000)
  rand <- ~ 1 | Order / Family / Genus / Species
  corr <- corHaversine(form = ~ lon + lat, mimic = "corExp")

  message("    fitting [none]...")
  m1 <- lme(fml_none, correlation = corr, random = rand,
            control = ctrl, method = "ML", data = dat)

  message("    fitting [inter-spp]...")
  m2 <- lme(fml_spp,  correlation = corr, random = rand,
            control = ctrl, method = "ML", data = dat)

  message("    fitting [inter-ind]...")
  m3 <- lme(fml_ind,  correlation = corr, random = rand,
            control = ctrl, method = "ML", data = dat)

  message("    fitting [inter-spp-ind]...")
  m4 <- lme(fml_both, correlation = corr, random = rand,
            control = ctrl, method = "ML", data = dat)

  tibble(
    conf_var = c("pd+ind+spp", "pd*spp+ind", "pd*ind+spp", "pd*indspp"),
    model    = list(m1, m2, m3, m4)
  )
}


# ============================================================================
# Run all partitioning models across displacement scales
# ============================================================================

for (time_scale_days in c(10, 1)) {
  for (spatial_scale in c("hi", "me")) {

    message("\n=== Partitioning models: ", time_scale_days, "d_", spatial_scale, " ===")

    dat <- load_and_prep_data(time_scale_days, spatial_scale) %>%
      prep_partitioned_data()

    models <- fit_partitioning_models(dat)

    out_file <- file.path(
      models_out_dir,
      paste0("Mods_sppVSind_", time_scale_days, "d_", spatial_scale, ".rds")
    )
    write_rds(models, out_file)
    message("  Saved: ", out_file)
  }
}


# ============================================================================
# Result extraction: AIC tables and fixed-effect summaries
# ============================================================================

extract_partitioning_results <- function(time_scale_days, spatial_scale) {

  scale_label <- paste0(time_scale_days, "d_", spatial_scale)
  in_file     <- file.path(
    models_out_dir,
    paste0("Mods_sppVSind_", scale_label, ".rds")
  )
  models <- read_rds(in_file)

  aic_tbl <- models %>%
    rowwise() %>%
    mutate(
      AIC = AIC(model),
      df  = attr(logLik(model), "df")
    ) %>%
    ungroup() %>%
    mutate(
      deltaAIC = AIC - min(AIC),
      scale    = scale_label
    ) %>%
    dplyr::select(scale, conf_var, df, AIC, deltaAIC)

  fixef_tbl <- models %>%
    rowwise() %>%
    mutate(
      fixef = list(
        as.data.frame(summary(model)$tTable) %>%
          rownames_to_column("Parameter")
      )
    ) %>%
    ungroup() %>%
    dplyr::select(conf_var, fixef) %>%
    unnest(fixef) %>%
    mutate(scale = scale_label) %>%
    dplyr::select(scale, conf_var, Parameter,
                  Estimate = Value, SE = Std.Error,
                  DF, t_value = `t-value`, p_value = `p-value`)

  list(aic = aic_tbl, fixef = fixef_tbl)
}

aic_all   <- list()
fixef_all <- list()

for (time_scale_days in c(10, 1)) {
  for (spatial_scale in c("hi", "me")) {
    res <- extract_partitioning_results(time_scale_days, spatial_scale)
    label <- paste0(time_scale_days, "d_", spatial_scale)
    aic_all[[label]]   <- res$aic
    fixef_all[[label]] <- res$fixef
  }
}

write_csv(bind_rows(aic_all),
          file.path(tables_out_dir, "partitioning_AIC_table.csv"))
write_csv(bind_rows(fixef_all),
          file.path(tables_out_dir, "partitioning_fixed_effects.csv"))

message("\nPartitioning AIC and fixed-effect tables written to ", tables_out_dir)
