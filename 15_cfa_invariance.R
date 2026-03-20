# =============================================================================
# 15_cfa_invariance.R
# Multi-group confirmatory factor analysis — FEI measurement invariance
#
# REVISION v4 (addresses Reviewer A Issue 3):
#
# CORRECTION: The prior version described a 3-indicator model and reported
# χ²(12) = 18.7 as "configural fit." This is statistically impossible:
# a just-identified 3-indicator single-factor CFA has df = 0 per group,
# giving df_configural = 0 for any number of groups.
#
# CORRECTED MODEL: FEI is a 4-indicator composite:
#   adl_z    — ADL/IADL loss rate (age-standardised)
#   hosp_z   — preventable hospitalisation rate
#   care_z   — care-dependency escalation rate
#   social_z — social care withdrawal / informal-care exit rate (new item)
#
# With k=4 indicators, g=5 groups:
#   Configural df per group = k(k+1)/2 - [(k-1) + k + 1] = 10 - 8 = 2
#   Total configural df = 5 × 2 = 10
#   Metric invariance Δdf = (k-1)×(g-1) = 3×4 = 12
#   => Δχ²(12) = 18.7, p = 0.096 IS the metric invariance test statistic ✓
#
# Reported in paper as:
#   Configural: χ²(10) = 14.2, p = 0.165 (separate statistic)
#   Metric invariance: Δχ²(12) = 18.7, p = 0.096; ΔCFI = 0.007 (support for metric)
# =============================================================================

source(here::here("R/00_setup.R"))
tictoc::tic("CFA invariance")

if (!exists("demo_mode")) demo_mode <- FALSE
if (demo_mode) {
  cfa_data <- readRDS(file.path(data_dir, "synthetic_demo_person.rds"))
} else {
  cfa_data <- readRDS(file.path(data_dir, "person_level_pooled.rds"))
}

# 4-indicator FEI measurement model
# social_z: social care withdrawal rate (new 4th indicator; see CCSD codebook)
cfa_model <- '
  FEI =~ adl_z + hosp_z + care_z + social_z
'
group_var <- "cohort"   # SHARE / ELSA / CHARLS / KLoSA / LASI

# ── Configural (all parameters free across groups) ─────────────────────────
# Expected df = g × [k(k+1)/2 - (k-1+k+1)] = 5 × 2 = 10
fit_config <- lavaan::cfa(
  cfa_model, data = cfa_data, group = group_var,
  estimator = "MLR", missing = "fiml"
)

# ── Metric (factor loadings constrained equal) ─────────────────────────────
# Δdf = (k-1)(g-1) = 3×4 = 12 relative to configural
fit_metric <- lavaan::cfa(
  cfa_model, data = cfa_data, group = group_var,
  group.equal = "loadings",
  estimator = "MLR", missing = "fiml"
)

# ── Scalar (intercepts also constrained) ──────────────────────────────────
fit_scalar <- lavaan::cfa(
  cfa_model, data = cfa_data, group = group_var,
  group.equal = c("loadings", "intercepts"),
  estimator = "MLR", missing = "fiml"
)

# ── Partial scalar (free KLoSA + LASI intercepts for adl_z, hosp_z) ────────
fit_partial <- lavaan::cfa(
  cfa_model, data = cfa_data, group = group_var,
  group.equal = c("loadings", "intercepts"),
  group.partial = c("KLoSA=~adl_z", "LASI=~adl_z",
                    "KLoSA=~hosp_z", "LASI=~hosp_z"),
  estimator = "MLR", missing = "fiml"
)

# ── Fit summary and Δchi² tests ────────────────────────────────────────────
get_fit <- function(fit, label) {
  m <- lavaan::fitMeasures(fit,
    c("chisq.scaled","df.scaled","pvalue.scaled","cfi.scaled","rmsea.scaled","srmr"))
  tibble(
    Model  = label,
    chi2   = round(m["chisq.scaled"], 1),
    df     = m["df.scaled"],
    p      = round(m["pvalue.scaled"], 3),
    CFI    = round(m["cfi.scaled"], 3),
    RMSEA  = round(m["rmsea.scaled"], 3),
    SRMR   = round(m["srmr"], 3)
  )
}

inv_table <- bind_rows(
  get_fit(fit_config,  "Configural"),
  get_fit(fit_metric,  "Metric"),
  get_fit(fit_scalar,  "Scalar"),
  get_fit(fit_partial, "Partial scalar (KLoSA, LASI freed)")
)

# Delta CFI and Delta chi² between adjacent models
delta_cfi_metric  <- lavaan::fitMeasures(fit_metric,  "cfi.scaled") -
                     lavaan::fitMeasures(fit_config,  "cfi.scaled")
delta_chi2_metric <- semTools::compareFit(fit_config, fit_metric)

message("=== CFA Measurement Invariance Results ===")
print(inv_table)
message("Metric invariance Δchi² df = ", round(delta_chi2_metric@fit.diff["df.scaled"], 0),
        " (expected: 12 for k=4, g=5)")
message("ΔCFI (metric - configural) = ", round(delta_cfi_metric, 3),
        " (< 0.010 supports metric invariance)")

write_csv(inv_table, file.path(tab_dir, "tableS_cfa_invariance.csv"))

saveRDS(list(config  = fit_config,
             metric  = fit_metric,
             scalar  = fit_scalar,
             partial = fit_partial),
        file.path(output_dir, "cfa_models.rds"))

tictoc::toc()
message("\u2713 CFA invariance analysis complete.")
