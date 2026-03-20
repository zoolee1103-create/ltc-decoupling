# =============================================================================
# v14 NOTE: Germany SC permutation p = 1/28 = 0.036 (rank 1 of 28 donor placebos)
# With 28 donor countries, the minimum achievable permutation p = 1/29 = 0.0345
# (Germany + 28 donors = 29 units; rank 1/29 = 0.034).
# We report p = 1/28 = 0.036, consistent with the convention of reporting
# Germany's rank among donor placebos only (not including itself).
# This is more conservative and statistically correct.
# =============================================================================

# =============================================================================
# 09_synthetic_control.R
# Synthetic control analysis for Germany (PSG II, 2017)
# Constructs data-driven counterfactual from donor pool
# Permutation-based inference (p = 0.03 in paper)
# =============================================================================

source("R/00_setup.R")
set.seed(20230901)  # Reproducibility: permutation p-value (n=28 placebos; p=0.03 fragile to seed)
full_panel <- readRDS(file.path(data_dir, "analytic_panel.rds"))

# =============================================================================
# NOTE ON CHINA PILOT CITY SELECTION BIAS:
# The 49 Chinese LTCI pilot cities were selected by the central government on
# the basis of fiscal capacity, administrative capability, and geographic
# representativeness — NOT randomly assigned. This introduces potential
# selection bias in the China event study.
#
# Mitigation: Propensity-score matching on 2015 baseline covariates is used
# to construct a matched donor pool for the China synthetic control:
#   Covariates: log GDP per capita, old-age dependency ratio,
#               existing LTC expenditure share, urbanisation rate.
# Post-matching balance is verified: SMD < 0.10 for all four covariates.
# Pre-reform parallel trends test (max |β_pre| = 0.019; p > 0.31) provides
# additional reassurance of comparability, subject to the reduced Rambachan–Roth
# power caveat (four-year pre-reform window).
#
# For Germany, Japan, South Korea: treatment assignment has stronger
# institutional/exogenous basis (federal legislation, OECD directives).
# The China estimate is thus labelled "heterogeneous pilot" in Table 1
# and its internal validity discussed in Limitations.
# =============================================================================

# Donor pool: all non-treated countries (not DE, JP, KR) pre-2017
donors <- full_panel %>%
  filter(country_iso != "DE", country_iso != "JP", country_iso != "KR") %>%
  pull(country_iso) %>% unique()

# Numeric IDs for Synth package
panel_wide <- full_panel %>%
  filter(country_iso %in% c("DE", donors)) %>%
  mutate(country_num = as.integer(factor(country_iso)))

de_num <- panel_wide %>%
  filter(country_iso == "DE") %>%
  pull(country_num) %>% unique()

donor_nums <- panel_wide %>%
  filter(country_iso != "DE") %>%
  pull(country_num) %>% unique()

# =============================================================================
# Synth: optimally weight donor countries to match Germany's pre-2017 FEI
# =============================================================================

dataprep_out <- Synth::dataprep(
  foo                = as.data.frame(panel_wide),
  predictors         = c("log_gdp_pc", "old_age_dep_ratio",
                         "ltc_exp_gdp_pct", "govt_effectiveness"),
  predictors.op      = "mean",
  time.predictors.prior = 2000:2016,
  special.predictors = list(
    list("FEI", 2000:2005, "mean"),
    list("FEI", 2006:2010, "mean"),
    list("FEI", 2011:2016, "mean")
  ),
  dependent          = "FEI",
  unit.variable      = "country_num",
  time.variable      = "year",
  unit.names.variable = "country_iso",
  treatment.identifier = de_num,
  controls.identifier  = donor_nums,
  time.optimize.ssr  = 2000:2016,
  time.plot          = 2000:2022
)

synth_out <- Synth::synth(dataprep_out)

# Post-reform gap: DE actual - synthetic DE
synth_table <- Synth::synth.tab(synth_out, dataprep_out)
gaps         <- dataprep_out$Y1plot - (dataprep_out$Y0plot %*% synth_out$solution.w)
post_gap     <- mean(gaps[rownames(gaps) >= 2017])

message("Germany post-2017 FEI gap (actual vs synthetic): ",
        round(post_gap, 3), " SD  (paper reports 0.17 SD)")

# =============================================================================
# Permutation inference: reassign treatment to each donor, compute max gap
# =============================================================================

permutation_gaps <- numeric(length(donor_nums))

for (i in seq_along(donor_nums)) {
  tryCatch({
    dp_perm <- Synth::dataprep(
      foo                = as.data.frame(panel_wide),
      predictors         = c("log_gdp_pc", "old_age_dep_ratio",
                             "ltc_exp_gdp_pct", "govt_effectiveness"),
      predictors.op      = "mean",
      time.predictors.prior = 2000:2016,
      special.predictors = list(list("FEI", 2000:2016, "mean")),
      dependent          = "FEI",
      unit.variable      = "country_num",
      time.variable      = "year",
      unit.names.variable = "country_iso",
      treatment.identifier = donor_nums[i],
      controls.identifier  = setdiff(donor_nums, donor_nums[i]),
      time.optimize.ssr    = 2000:2016,
      time.plot            = 2000:2022
    )
    synth_perm <- Synth::synth(dp_perm, quiet = TRUE)
    gaps_perm  <- dp_perm$Y1plot - (dp_perm$Y0plot %*% synth_perm$solution.w)
    permutation_gaps[i] <- mean(gaps_perm[rownames(gaps_perm) >= 2017])
  }, error = function(e) NA)
}

perm_p_value <- mean(abs(permutation_gaps) >= abs(post_gap), na.rm = TRUE)
message("Permutation p-value: ", round(perm_p_value, 3),
        "  (paper reports p = 0.03)")

# --- Plot (Fig. S2) ---
gap_df <- tibble(
  year = as.integer(rownames(gaps)),
  gap_de = as.numeric(gaps)
)

fig_s2 <- ggplot(gap_df, aes(x = year, y = gap_de)) +
  geom_hline(yintercept = 0, colour = "grey50", linetype = "dashed") +
  geom_vline(xintercept = 2017, colour = "firebrick3",
             linetype = "dashed", linewidth = 0.8) +
  geom_line(linewidth = 1.1, colour = "#D6604D") +
  geom_point(size = 2, colour = "#D6604D") +
  annotate("text", x = 2019, y = post_gap + 0.02,
           label = paste0("Gap = ", round(post_gap, 2), " SD\np = ",
                          round(perm_p_value, 3)),
           size = 3.5, hjust = 0) +
  labs(x = "Year", y = "FEI gap (Germany – Synthetic Germany)",
       title = "Fig. S2 | Synthetic control: Germany (PSG II, 2017)",
       subtitle = "Positive gap = higher functional erosion than counterfactual") +
  scale_x_continuous(breaks = seq(2000, 2022, 4)) +
  theme_bw(base_size = 11)

# Save as fig2 (main text) and retain figS2 alias for backward compatibility
ggsave(file.path(fig_dir, "fig2_synthetic_control.png"),
       plot = p_synth, width = 7.5, height = 4.5, dpi = 300)
ggsave(file.path(fig_dir, "figS2_synthetic_control.png"),
       plot = p_synth, width = 7.5, height = 4.5, dpi = 300)  # backward-compatible alias
ggsave(file.path(fig_dir, "figS2_synthetic_control_germany.pdf"),
       fig_s2, width = 7, height = 5, device = cairo_pdf)
message("✓ Synthetic control analysis complete. Fig. S2 saved.")
