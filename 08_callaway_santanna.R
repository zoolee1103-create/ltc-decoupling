# =============================================================================
# 08_callaway_santanna.R
# Callaway & Sant'Anna (2021) group-time ATT estimator
# Robust to heterogeneous treatment effects in staggered DiD
# =============================================================================

source("R/00_setup.R")
full_panel <- readRDS(file.path(data_dir, "analytic_panel.rds"))

# Numeric country ID required by did package
full_panel <- full_panel %>%
  mutate(country_id = as.integer(factor(country_iso)))

# =============================================================================
# Callaway-Sant'Anna ATT(g,t)
# =============================================================================

cs_result <- did::att_gt(
  yname         = "FEI",
  tname         = "year",
  idname        = "country_id",
  gname         = "first_treat_year",  # year of first treatment (0 if never)
  xformla       = ~ log_gdp_pc + old_age_dep_ratio + govt_effectiveness,
  data          = full_panel,
  panel         = TRUE,
  control_group = "nevertreated",
  bstrap        = TRUE,
  biters        = 1000,
  clustervars   = "country_id"
)

# Aggregate to simple treatment effect
cs_agg <- did::aggte(cs_result, type = "simple")
message("Callaway-Sant'Anna GATT = ", round(cs_agg$overall.att, 3),
        "  SE = ", round(cs_agg$overall.se, 3))
message("(Paper reports GATT = 0.138, SE = 0.029)")

# Event-study aggregation
cs_es <- did::aggte(cs_result, type = "dynamic", min_e = -5, max_e = 5)
did::ggdid(cs_es) +
  labs(title = "CS Event-study (robustness)",
       x = "Periods relative to treatment",
       y = "ATT (SD units)") +
  theme_bw(base_size = 11)
ggsave(file.path(fig_dir, "figS_cs_eventstudy.pdf"), width = 7, height = 5)

saveRDS(cs_result, file.path(output_dir, "cs_result.rds"))
message("✓ Callaway-Sant'Anna estimation complete.")
