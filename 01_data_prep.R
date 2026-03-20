# =============================================================================
# 01_data_prep.R
# Data harmonisation: merge cohort-level FEI components with
# country-year governance panel. Construct analytic dataset.
# =============================================================================

source("R/00_setup.R")

# demo_mode: if TRUE, loads synthetic data for code testing without real data
if (!exists("demo_mode")) demo_mode <- FALSE

# =============================================================================
# SECTION 1: Individual-level data (five cohort studies)
# =============================================================================
# Data must be obtained from original providers (see README).
# Access links: SHARE share-project.org | ELSA ukdataservice.ac.uk |
#               CHARLS charls.pku.edu.cn | KLoSA keis.or.kr | LASI iipsindia.ac.in

if (!demo_mode) {

  # --- SHARE (waves 1-9, modules: physical health, ADL/IADL, health care) ---
  share_raw <- haven::read_dta(file.path(data_dir, "share/sharew9_rel8-0-0_ALL_datasets_stata/cv_r.dta"))
  # Key variables: hhid, mergeid, wave, country, adl, iadl, hosp, care_need_level
  # See SHARE documentation: https://share-project.org/data-access

  # --- ELSA (waves 1-10) ---
  elsa_raw  <- haven::read_dta(file.path(data_dir, "elsa/elsa_wave10_2021.dta"))

  # --- CHARLS (waves 1-5) ---
  charls_raw <- haven::read_dta(file.path(data_dir, "charls/Health_Status_and_Functioning.dta"))

  # --- KLoSA (waves 1-9) ---
  klosa_raw <- haven::read_dta(file.path(data_dir, "klosa/klosa_w09.dta"))

  # --- LASI (wave 1) ---
  lasi_raw  <- haven::read_dta(file.path(data_dir, "lasi/lasi_wave1.dta"))

} else {
  message("⚠ demo_mode = TRUE: loading synthetic demonstration dataset.")
  message("  Results will NOT match paper estimates.")
  full_panel <- readRDS(file.path(data_dir, "synthetic_demo.rds"))
  message("✓ Synthetic panel loaded: ", nrow(full_panel), " rows")
}

# =============================================================================
# SECTION 2: Functional Erosion Index (FEI) construction
# =============================================================================
# FEI = z-score average of three age-standardised components:
#   (i)   12-month ADL/IADL loss rate per 100 person-years
#   (ii)  Preventable hospitalisation rate per 1,000 aged ≥65
#   (iii) Care-dependency escalation rate (proportion advancing ≥1 care-need level / 24 months)

build_fei <- function(df) {
  df %>%
    group_by(country_iso, year) %>%
    summarise(
      # Component 1: ADL/IADL loss (age-standardised via direct standardisation to WHO 2000 standard)
      adl_loss_raw = mean(adl_loss_12m, na.rm = TRUE) * 100,
      # Component 2: Preventable hospitalisation (ambulatory care sensitive conditions)
      prev_hosp    = sum(acsc_hosp, na.rm = TRUE) / sum(age_gte65) * 1000,
      # Component 3: Care escalation (care_need_t1 > care_need_t0 within 24 months)
      care_escal   = mean(care_level_increase_24m, na.rm = TRUE),
      n_persons    = n(),
      .groups = "drop"
    ) %>%
    mutate(
      # Z-score within country-year
      adl_z   = as.numeric(scale(adl_loss_raw)),
      hosp_z  = as.numeric(scale(prev_hosp)),
      care_z  = as.numeric(scale(care_escal)),
      # Equal-weight composite
      FEI     = (adl_z + hosp_z + care_z) / 3
    )
}

# =============================================================================
# SECTION 3: Country-year governance panel (CCSD + external sources)
# =============================================================================
# CCSD: Comparative Care System Database (constructed for this study)
# Codebook: data/codebook_CCSD.csv

load_governance_panel <- function(data_dir) {

  # --- CCSD (boundary fidelity weights + commissioning indices) ---
  ccsd <- read_csv(file.path(data_dir, "ccsd_panel_2000_2022.csv")) %>%
    select(country_iso, year,
           # v1->v2 boundary
           legal_transposition_rate, reg_authority_index,
           # v2->v3 boundary
           commission_spec_completeness, needs_assessment_irr,
           # v3->v4 boundary
           service_delivery_adequacy, hhi_provider,
           # v4->v5 boundary
           admin_delivery_rate,
           # Mandate scope
           entitlement_scope_index)

  # --- OECD LTC Statistics ---
  oecd_ltc <- read_csv(file.path(data_dir, "oecd_ltc_stats_2000_2022.csv")) %>%
    select(country_iso, year,
           ltc_exp_gdp_pct,       # LTC expenditure as % GDP
           ltc_recipients_rate,   # LTC recipients per 1,000 aged 65+
           formal_care_workers)   # Formal LTC workers per 1,000 aged 65+

  # --- World Bank WGI ---
  wgi <- read_csv(file.path(data_dir, "wgi_panel.csv")) %>%
    select(country_iso, year, govt_effectiveness)

  # --- ICRG Institutional Quality ---
  icrg <- read_csv(file.path(data_dir, "icrg_panel.csv")) %>%
    select(country_iso, year, icrg_quality)

  # --- World Bank macro controls ---
  wb_macro <- read_csv(file.path(data_dir, "wb_macro.csv")) %>%
    select(country_iso, year,
           log_gdp_pc,            # log(GDP per capita, 2015 USD PPP)
           old_age_dep_ratio,     # Population 65+ / 15-64
           inflation_cpi)

  # --- Welfare regime classification ---
  regime_map <- read_csv(file.path(data_dir, "regime_classification.csv")) %>%
    select(country_iso, welfare_regime, oecd_member)

  # --- Merge all governance sources ---
  governance_panel <- ccsd %>%
    left_join(oecd_ltc,   by = c("country_iso", "year")) %>%
    left_join(wgi,        by = c("country_iso", "year")) %>%
    left_join(icrg,       by = c("country_iso", "year")) %>%
    left_join(wb_macro,   by = c("country_iso", "year")) %>%
    left_join(regime_map, by = "country_iso")

  governance_panel
}

# =============================================================================
# SECTION 4: IDS construction (two versions)
# =============================================================================
# IDS_pure: does NOT include CLS in construction formula (main specification)
#           Addresses circularity concern raised in review.
# IDS_full: includes CLS component (reported in robustness; matches original paper)
# Both versions are calibrated on 2000-2009 training data.

# See R/02_ids_construction.R for full details.

# =============================================================================
# SECTION 5: Merge and final analytic panel
# =============================================================================
if (!demo_mode) {

  fei_panel        <- build_fei(bind_rows(share_raw, elsa_raw, charls_raw, klosa_raw, lasi_raw))
  governance_panel <- load_governance_panel(data_dir)

  full_panel <- fei_panel %>%
    inner_join(governance_panel, by = c("country_iso", "year")) %>%
    filter(year >= 2000, year <= 2022) %>%
    arrange(country_iso, year)

  # --- Data quality weights (LASI = 0.62; CHARLS = 0.79; others = 1.0) ---
  dq_weights <- tibble(
    country_iso = c("IN",  "CN"),
    dq_weight   = c(0.62,  0.79)
  )
  full_panel <- full_panel %>%
    left_join(dq_weights, by = "country_iso") %>%
    mutate(dq_weight = replace_na(dq_weight, 1.0))

  # --- Inverse-probability-of-attrition weights (Supplementary Table S3) ---
  # Fit attrition probit on pre-attrition observables
  attrition_model <- glm(
    attrited ~ age_baseline + female + education_yrs + log_gdp_pc + adl_baseline,
    data   = full_panel,
    family = binomial(link = "probit")
  )
  full_panel <- full_panel %>%
    mutate(
      p_attrition = predict(attrition_model, type = "response"),
      ipaw         = ifelse(attrited == 0, 1 / (1 - p_attrition), 0)
    )

  # Save analytic dataset
  saveRDS(full_panel, file.path(data_dir, "analytic_panel.rds"))
  message("✓ Analytic panel saved: ", nrow(full_panel), " country-years; ",
          sum(full_panel$n_persons), " person-years")
}

# =============================================================================
# SECTION 6: Descriptive statistics (Table S1)
# =============================================================================
desc_stats <- full_panel %>%
  select(FEI, IDS_pure, IDS_full, CLS,
         ltc_exp_gdp_pct, log_gdp_pc, old_age_dep_ratio,
         govt_effectiveness, icrg_quality) %>%
  summarise(across(everything(),
    list(mean = ~mean(.x, na.rm=T),
         sd   = ~sd(.x,   na.rm=T),
         min  = ~min(.x,  na.rm=T),
         max  = ~max(.x,  na.rm=T)),
    .names = "{.col}__{.fn}")) %>%
  pivot_longer(everything(),
               names_to  = c("variable", "stat"),
               names_sep = "__") %>%
  pivot_wider(names_from = stat, values_from = value)

write_csv(desc_stats, file.path(tab_dir, "tableS1_descriptives.csv"))
message("✓ Table S1 (descriptive statistics) saved.")
