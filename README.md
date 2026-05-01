# Settlement porosity and global mammalian movement

This repository contains the code and data needed to reproduce the
statistical models and tables reported in:

> [Author list]. _[Paper title]_. [Journal], [Year]. DOI: [forthcoming]

A frozen archive of this repository is also deposited on Zenodo:
[Zenodo DOI: forthcoming]

[note: data not visible for now.]

## Overview

This project is a subproject of the Global Barrier Project.

We use a global telemetry dataset of terrestrial mammals to test whether
the spatial **configuration** of human-modified landscapes — quantified as
*Settlement porosity*, a percolation-based metric of how settlement
patches are arranged — affects mammalian movement beyond the effect of
landscape **composition** (the amount or intensity of human modification,
quantified as Human Footprint Index, Human Modification Index, or
settlement cover).

This repository implements the three sets of analyses reported in the
manuscript:

1. **Main models** — Three nested models (composition only; composition +
   porosity; composition × porosity) for each combination of displacement
   scale (1-day, 10-day) and summary metric (median, 0.95 quantile),
   testing whether configuration explains variance in mammalian movement
   beyond composition at different scales.

2. **Partitioning models** — Within-species centering of porosity into
   species-level (environmental filtering) and individual-level
   (behavioral plasticity) components, used to distinguish the two
   non-mutually-exclusive mechanisms underlying the porosity–movement
   relationship. Uses HFI as the composition variable.

3. **Trait-based interaction models** — Supplementary models testing
   whether body mass and diet modulate movement responses to configuration (i.e. porosity) and
   to composition (i.e. HFI).

## Repository structure

```
.
├── README.md                     This file
├── LICENSE                       Code license
├── code/forGitHub
│   ├── 00_setup.R                Shared dependencies and helper functions
│   ├── 01_main_models.R          Main composition vs. configuration models
│   ├── 02_partitioning_models.R  Behavioral plasticity vs. environmental
│   │                              filtering partitioning models
│   ├── 03_trait_interaction_models.R
│   │                             Supplementary trait × landscape models
│   └── HaversineLMEfunctions.R   Custom Haversine spatial-correlation
│                                  function for use with nlme::lme()
├── data/
│   ├── global_barrier_mod_input_10d_hi.rds
│   ├── global_barrier_mod_input_10d_me.rds
│   ├── global_barrier_mod_input_1d_hi.rds
│   └── global_barrier_mod_input_1d_me.rds
└── results/
    ├── models/                   Created on first run; stores fitted models
    └── tables/                   Created on first run; stores AIC and
                                   fixed-effect tables in CSV format
```

## How to reproduce the analysis

### 1. System requirements

- R version ≥ 4.3.2 (developed and tested on 4.3.2)
- A working internet connection is **not** required to reproduce the
  models — all input data are provided in `data/`.

### 2. Required R packages

The scripts depend on the following CRAN packages:

| Package      | Purpose                                              |
|--------------|------------------------------------------------------|
| `nlme`       | Linear mixed-effects models                          |
| `AICcmodavg` | AIC-based model comparison and predictions           |
| `MuMIn`      | Marginal and conditional R^2                         |
| `tidyverse`  | Data wrangling (`readr`, `dplyr`, `tidyr`, ...)      |

Install them with:

```r
install.packages(c("nlme", "AICcmodavg", "MuMIn", "tidyverse"))
```

The custom `corHaversine` correlation structure is provided as part of
this repository in `code/HaversineLMEfunctions.R` and is sourced
automatically by `00_setup.R`.

### 3. Run order

Each analysis script is self-contained and can be run independently. From
the repository root:

```r
# Main models (composition vs. configuration)
source("./code/01_main_models.R")

# Partitioning models (behavioral plasticity vs. environmental filtering)
source("./code/02_partitioning_models.R")

# Trait-based interaction models
source("./code/03_trait_interaction_models.R")
```

Each script will:
1. Load and preprocess the four input `.rds` files from `data/`.
2. Fit the relevant set of mixed-effects models (this is the slowest
   step; expect minutes to ~hour per script depending on hardware,
   primarily because of the spatial-correlation structure).
3. Save fitted models as `.rds` files in `results/models/`.
4. Extract AIC tables and fixed-effect summaries and write them as CSV
   files to `results/tables/`.

All scripts use **maximum likelihood** (`method = "ML"`) so that AIC
values are directly comparable across models.

## Data description

Each of the four input `.rds` files contains an individual-level dataset
of mammal movement, environmental covariates, and species traits. The
files differ in the displacement metric used as the response variable:

| File                                        | Displacement window | Summary metric          |
|---------------------------------------------|---------------------|-------------------------|
| `global_barrier_mod_input_1d_me.rds`        | 1 day               | Median                  |
| `global_barrier_mod_input_1d_hi.rds`        | 1 day               | 0.95 quantile (long-distance) |
| `global_barrier_mod_input_10d_me.rds`       | 10 days             | Median                  |
| `global_barrier_mod_input_10d_hi.rds`       | 10 days             | 0.95 quantile (long-distance) |

Each row represents one individual animal (i.e. data have been collapsed
to per-individual averages of all covariates prior to modeling).

### Data dictionary

| Column                | Type      | Description                                                       |
|-----------------------|-----------|-------------------------------------------------------------------|
| `ID`                  | character | Unique individual identifier                                      |
| `Binomial`            | character | Species scientific binomial name                                  |
| `Order`               | character | Taxonomic order                                                   |
| `Family`              | character | Taxonomic family                                                  |
| `Genus`               | character | Taxonomic genus                                                   |
| `Species`             | character | Taxonomic species (specific epithet)                              |
| `Diet`                | factor    | Dietary guild: `Carnivore`, `Herbivore`, or `Omnivore`            |
| `BodyMass_kg`         | numeric   | Species-average body mass (kg)                                    |
| `Displacement_km`     | numeric   | Individual-level displacement (median or 0.95 quantile, depending on file) |
| `Longitude`           | numeric   | Individual-mean longitude (decimal degrees)                       |
| `Latitude`            | numeric   | Individual-mean latitude (decimal degrees)                        |
| `pd_adpt_km`          | numeric   | Settlement porosity (km), derived from the Global Settlement Percolation dataset; matched to the species' mobility scale |
| `bd_adpt`             | numeric   | Settlement cover (%) at the same matched scale                    |
| `HFI`                 | numeric   | Human Footprint Index averaged over the displacement path          |
| `HMI`                 | numeric   | Human Modification Index averaged over the displacement path       |
| `NDVI`                | numeric   | NDVI averaged over the displacement path (scaled value, multiply by 0.0001 for original) |

Within each script, continuous predictors are mean-centered (and, in
`03_trait_interaction_models.R`, log body mass is additionally
median-centered) prior to fitting; the raw values in the `.rds` files are
provided unscaled so that downstream users can apply alternative
transformations if desired.

## Output files

Running all three scripts produces:

- `results/models/Mods_globalbarrier_logdisp_pd_<time>d_<scale>.rds`
  Main models (3 per composition variable × 3 composition variables = 9
  per scale).
- `results/models/Mods_sppVSind_<time>d_<scale>.rds`
  Partitioning models (4 candidate structures per scale).
- `results/models/Mods_trait_interactions_<time>d_<scale>.rds`
  Trait-interaction models (4 candidate structures per scale).
- `results/tables/main_model_AIC_table.csv`
- `results/tables/main_model_fixed_effects.csv`
- `results/tables/partitioning_AIC_table.csv`
- `results/tables/partitioning_fixed_effects.csv`
- `results/tables/trait_interaction_fixed_effects.csv`

Figure-generation code is not included in this repository, as figures are
purely visualizations of model output that can be regenerated from the
model objects and tables produced above.

## License

Code in this repository is released under the [MIT License](LICENSE).
The accompanying telemetry data are released under [data license:
forthcoming].

## Citation

If you use code or data from this repository, please cite both the paper
and the Zenodo archive:

> [Paper citation — to be added on acceptance]
>
> [Zenodo citation — to be added on archive creation]

## Contact

For questions about the analysis or reproducibility, contact:

- Wenjing Xu — wenjingxu [at] umass.edu — University of Massachusetts Amherst
