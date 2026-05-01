# ============================================================================
# 00_setup.R
# ----------------------------------------------------------------------------
# Shared dependencies and helper functions for all model scripts.
#
# This script is sourced by 01_main_models.R, 02_partitioning_models.R, and
# 03_trait_interaction_models.R. It loads required packages, defines the
# Haversine spatial-correlation function used throughout the analysis, and
# sets the input/output paths used by all downstream scripts.
#
# To reproduce the analysis, ensure the following structure is in place:
#
#   project_root/
#   |-- code/
#   |   |-- 00_setup.R                    (this file)
#   |   |-- 01_main_models.R
#   |   |-- 02_partitioning_models.R
#   |   |-- 03_trait_interaction_models.R
#   |   `-- HaversineLMEfunctions.R
#   |-- data/
#   |   |-- global_barrier_mod_input_10d_hi.rds
#   |   |-- global_barrier_mod_input_10d_me.rds
#   |   |-- global_barrier_mod_input_1d_hi.rds
#   |   `-- global_barrier_mod_input_1d_me.rds
#   `-- results/
#       |-- models/                       (created if missing)
#       `-- tables/                       (created if missing)
# ============================================================================

# ---- Required packages -----------------------------------------------------
# Models and statistics
library(nlme)         # linear mixed-effects models
library(AICcmodavg)   # AIC-based model comparison and predictions
library(MuMIn)        # marginal and conditional R^2

# Data wrangling
library(tidyverse)    # readr, dplyr, tidyr, etc.

# ---- Custom Haversine spatial-correlation function -------------------------
# Defines `corHaversine`, which computes great-circle distances between
# observation locations for use in the spatial-correlation structure of
# `lme()` calls. This is required because nlme's built-in correlation
# structures use Euclidean distances, which are inappropriate for
# global-scale lon/lat data.
source("./code/HaversineLMEfunctions.R")

# ---- Paths -----------------------------------------------------------------
# Paths are relative to the repository root. Adjust the values below if your
# local directory structure differs from the layout described above.
data_dir         <- "..."   # e.g. "./data"
results_dir      <- "..."   # e.g. "./results"
models_out_dir   <- file.path(results_dir, "models")
tables_out_dir   <- file.path(results_dir, "tables")

dir.create(models_out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tables_out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Helper: load and preprocess the input dataset -------------------------
# Reads the species-individual displacement dataset for a given time window
# and movement scale, applies the standard transformations used across all
# model scripts, and filters out species with fewer than 5 individuals.
#
# Arguments:
#   time_scale_days : numeric; either 1 or 10
#   spatial_scale   : character; either "hi" (0.95 quantile of displacements)
#                     or "me" (median of displacements)
# Returns:
#   tibble of preprocessed individual-level data ready for modeling.
load_and_prep_data <- function(time_scale_days, spatial_scale) {

  fname <- file.path(
    data_dir,
    paste0("global_barrier_mod_input_",
           time_scale_days, "d_", spatial_scale, ".rds")
  )

  dat <- read_rds(fname) %>%
    mutate(
      log_BodyMass_kg     = log(BodyMass_kg),
      scale_NDVI          = as.vector(NDVI * 0.0001),  # restore native NDVI scale
      log_Displacement_km = log(Displacement_km),
      # Mean-center continuous predictors so that interaction-term coefficients
      # are interpretable at the population mean of each predictor.
      log_bd     = scale(log(bd_adpt), scale = FALSE),
      HFI        = scale(HFI,          scale = FALSE),
      HMI        = scale(HMI,          scale = FALSE),
      pd_adpt_km = scale(pd_adpt_km,   scale = FALSE)
    ) %>%
    rename(lon = Longitude, lat = Latitude)

  # Restrict to species with at least 5 individuals to ensure model
  # convergence and adequate sample representativeness.
  spp_keep <- dat %>%
    group_by(Binomial) %>%
    summarize(n = length(unique(ID)), .groups = "drop") %>%
    filter(n >= 5) %>%
    pull(Binomial)

  dat %>% filter(Binomial %in% spp_keep)
}
