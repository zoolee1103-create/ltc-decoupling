# =============================================================================
# 00_setup.R  —  Package installation and global settings
# Why Care Systems Fail: Institutional Decoupling in Ageing Societies
# Tested on R 4.3.1 (2023-06-16)
# =============================================================================

required_packages <- c(
  # Project infrastructure
  "here",          # Portable file paths (here::here())
  "tictoc",        # Execution timing
  # Data wrangling
  "tidyverse",     # Core: dplyr, tidyr, purrr, ggplot2, readr, stringr
  "haven",         # Read Stata .dta files
  "readxl",        # Read Excel
  "lubridate",     # Date handling
  # Econometrics — core
  "fixest",        # TWFE and IV via feols(); wald() for LR tests
  "lmtest",        # Coefficient tests
  "sandwich",      # HC/clustered variance estimators
  # Wild Cluster Bootstrap (N=32 clusters; Cameron & Miller 2015)
  "fwildclusterboot",  # boottest() for WCB inference
  # Staggered DiD robust estimators
  "did",           # Callaway & Sant'Anna (2021) att_gt()
  "bacondecomp",   # Goodman-Bacon decomposition
  "HonestDiD",     # Rambachan & Roth (2023) sensitivity
  # Synthetic control
  "Synth",         # Abadie, Diamond & Hainmueller (2010)
  # Causal mediation
  "mediation",     # Imai, Keele & Yamamoto (2010)
  "EValue",        # VanderWeele & Ding (2017) E-values
  # Measurement models
  "lavaan",        # Confirmatory factor analysis
  "semTools",      # Multi-group CFA invariance tests
  # Entropy / IPW weighting
  "WeightIt",      # Entropy balancing and IPW
  "cobalt",        # Covariate balance after weighting
  # Survival / count outcomes
  "MASS",          # miscellaneous statistical utilities (e.g., MASS::kde2d, negative-binomial family)
  # Tables and figures
  "modelsummary",  # Regression tables
  "gt",            # Table formatting
  "kableExtra",    # LaTeX/HTML tables
  "ggplot2",       # Figures (loaded via tidyverse but explicit for clarity)
  "patchwork",     # Multi-panel figures
  "ggrepel",       # Non-overlapping labels
  "scales",        # Axis formatting
  "RColorBrewer",  # Colour palettes
  "broom",         # tidy() for model objects
  "Cairo"          # High-quality PDF output
)

new_pkgs <- required_packages[
  !(required_packages %in% installed.packages()[,"Package"])
]
if (length(new_pkgs) > 0) {
  message("Installing ", length(new_pkgs), " missing packages...")
  install.packages(new_pkgs, dependencies = TRUE,
                   repos = "https://cloud.r-project.org")
}
invisible(lapply(required_packages, library, character.only = TRUE,
                 warn.conflicts = FALSE, quietly = TRUE))

# --- Global options -----------------------------------------------------------
options(scipen = 999, digits = 4,
        dplyr.summarise.inform = FALSE,
        mc.cores = max(1L, parallel::detectCores() - 1L))

set.seed(20230901)  # Global seed; per-script seeds set locally for bootstrap

# --- Project paths ------------------------------------------------------------
proj_root  <- here::here()
data_dir   <- file.path(proj_root, "data")
output_dir <- file.path(proj_root, "output")
fig_dir    <- file.path(output_dir, "figures")
tab_dir    <- file.path(output_dir, "tables")
for (d in c(data_dir, fig_dir, tab_dir))
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)

# --- Regime colour palette ----------------------------------------------------
regime_colours <- c(
  "Beveridgean statist"      = "#2166AC",
  "Bismarckian corporatist"  = "#D6604D",
  "Residualist liberal"       = "#F4A582",
  "East Asian developmental" = "#4DAC26"
)

cat("\u2713 Setup complete | R", as.character(getRversion()), "\n")
