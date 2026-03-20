# =============================================================================
# 00_generate_synthetic_demo.R
# Synthetic demonstration datasets for code testing (no real data)
#
# Produces:
#   data/synthetic_demo.rds            — country-year panel  (N = 704)
#   data/synthetic_demo_person.rds     — person-level pooled (N ≈ 10,000)
#   data/analytic_panel.rds            — alias for downstream scripts
#   data/person_level_pooled.rds       — alias for downstream scripts
#   data/expert_panel_ratings.csv      — Delphi weight file (synthetic)
#
# Results from demo data will NOT match paper estimates.
# =============================================================================

source(here::here("R/00_setup.R"))
set.seed(42)

n_countries <- 32
years       <- 2000:2022

regimes <- rep(
  c("Beveridgean statist","Bismarckian corporatist",
    "Residualist liberal","East Asian developmental"),
  times = c(9, 8, 7, 8)
)
countries <- tibble(
  country_iso    = paste0("C", sprintf("%02d", seq_len(n_countries))),
  welfare_regime = regimes,
  oecd_member    = as.integer(seq_len(n_countries) <= 22)
)

# ----- Country-year panel ---------------------------------------------------
panel <- expand_grid(country_iso = countries$country_iso, year = years) %>%
  left_join(countries, by = "country_iso") %>%
  mutate(
    log_gdp_pc                 = rnorm(n(), 10.0, 0.8),
    old_age_dep_ratio          = rnorm(n(), 0.28, 0.07),
    ltc_exp_gdp_pct            = pmax(rnorm(n(), 1.5, 0.8), 0.1),
    govt_effectiveness         = rnorm(n(), 0.5,  0.7),
    icrg_quality               = rnorm(n(), 3.0,  0.5),
    # Fidelity weights (CCSD administrative indicators)
    legal_transposition_rate   = pmin(pmax(rnorm(n(), 0.80, 0.12), 0), 1),
    reg_authority_index        = pmin(pmax(rnorm(n(), 60,  15),   0), 100),
    commission_spec_completeness = pmin(pmax(rnorm(n(), 0.65, 0.15), 0), 1),
    needs_assessment_irr       = pmin(pmax(rnorm(n(), 0.70, 0.12), 0), 1),
    service_delivery_adequacy  = pmin(pmax(rnorm(n(), 0.72, 0.13), 0), 1),
    hhi_provider               = pmin(pmax(rnorm(n(), 0.20, 0.10), 0), 1),
    admin_delivery_rate        = pmin(pmax(rnorm(n(), 0.78, 0.12), 0), 1),
    entitlement_scope_index    = pmin(pmax(rnorm(n(), 55,  20),   0), 100),
    # Healthcare sector (IV exclusion restriction test)
    entitlement_scope_health   = pmin(pmax(rnorm(n(), 70,  15),   0), 100),
    GIS_health                 = pmin(pmax(rnorm(n(), 0.65, 0.12), 0), 1),
    # Outcome components
    adl_z    = rnorm(n(), 0, 1),
    hosp_z   = rnorm(n(), 0, 1),
    care_z   = rnorm(n(), 0, 1),
    social_z = rnorm(n(), 0, 1),   # social care withdrawal / informal-care exit rate
    # Reform variables
    reform_year = sample(c(2006, 2008, 2015, 2017, NA_integer_),
                          n(), replace = TRUE,
                          prob = c(.05, .05, .05, .05, .80)),
    first_treat_year = replace_na(reform_year, 0L),
    # Attrition / selection
    attrited  = rbinom(n(), 1, 0.05),
    observed  = rbinom(n(), 1, 0.95),
    ipaw      = 1.0,
    # Admin capacity controls (for robustness specs)
    admin_capacity_index = rnorm(n(), 0, 1),
    ltc_data_coverage    = pmin(pmax(rnorm(n(), 0.80, 0.15), 0), 1),
    # Falsification outcomes
    road_traffic_mortality = pmax(rnorm(n(), 6.0, 2.5), 0),
    cancer_incidence_std   = rnorm(n(), 0, 1),
    # Data quality weights
    dq_weight = case_when(
      country_iso %in% c("C28","C29") ~ 0.62,
      country_iso %in% c("C25","C26") ~ 0.79,
      TRUE ~ 1.0
    )
  ) %>%
  mutate(FEI = (adl_z + hosp_z + care_z + social_z) / 4)

saveRDS(panel, file.path(data_dir, "synthetic_demo.rds"))
saveRDS(panel, file.path(data_dir, "analytic_panel.rds"))
message("\u2713 Country-year panel: ", nrow(panel), " rows")

# ----- Expert panel ratings (Delphi; synthetic calibrated to published means) -
# Real file structure mirrors data/expert_panel_ratings.csv in CCSD
expert_ratings <- tibble(
  boundary    = c("v1v2", "v2v3", "v3v4", "v4v5"),
  mean_rating = c(3.48, 3.15, 2.82, 3.55),   # 1-5 Likert scale
  sd_rating   = c(0.58, 0.63, 0.71, 0.60),
  icc_lower   = c(0.74, 0.71, 0.69, 0.72),
  icc_upper   = c(0.87, 0.85, 0.84, 0.86)
)
write_csv(expert_ratings, file.path(data_dir, "expert_panel_ratings.csv"))
message("\u2713 Expert panel ratings saved.")

# ----- Person-level pooled dataset ------------------------------------------
n_persons <- 10000
cohorts   <- rep(c("SHARE","ELSA","CHARLS","KLoSA","LASI"), each = 2000)

person_data <- tibble(
  id              = seq_len(n_persons),
  cohort          = cohorts,
  country_iso     = sample(panel$country_iso, n_persons, replace = TRUE),
  year            = sample(2005:2020,          n_persons, replace = TRUE),
  age             = rnorm(n_persons,  72,  9),
  female          = rbinom(n_persons,  1, 0.54),
  edu_yrs         = rnorm(n_persons,  10,  4),
  adl_baseline    = rbinom(n_persons,  1, 0.22),
  adl_t1          = rbinom(n_persons,  1, 0.28),   # outcome BEFORE subsetting
  # P1 mediator: days to assessment (continuous; ~4 months mean)
  days_to_assessment     = pmax(rnorm(n_persons, 122, 65), 7),
  # P2 outcome
  informal_care_primary  = rbinom(n_persons, 1, 0.26),
  # P3 exposure and outcome (BOTH binary)
  provider_change_24m    = rbinom(n_persons, 1, 0.18),
  care_level_increase_24m = rbinom(n_persons, 1, 0.20),
  social_z = rnorm(n_persons, 0, 1),  # social care withdrawal z-score  # binary 0/1
  care_need_severity     = sample(1:5, n_persons, replace = TRUE),
  assessed_eligible      = rbinom(n_persons, 1, 0.65),
  log_gdp_pc             = rnorm(n_persons, 10.0, 0.8),
  old_age_dep_ratio      = rnorm(n_persons,  0.28, 0.07)
)

saveRDS(person_data, file.path(data_dir, "synthetic_demo_person.rds"))
saveRDS(person_data, file.path(data_dir, "person_level_pooled.rds"))
message("\u2713 Person-level synthetic data: ", nrow(person_data), " rows")
message("\u2713 All synthetic demo assets saved.")
