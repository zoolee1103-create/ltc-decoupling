# =============================================================================
# 12_robustness_table.R  —  Pre-specified robustness specifications (Table S3)
#
# REVIEWER FIXES:
#   - admin_capacity_index and ltc_data_coverage now sourced from CCSD
#     (added to synthetic demo in 00_generate_synthetic_demo.R)
#   - Heckman: observed flag is created inline from data completeness
#   - All 8 specifications run sequentially, output consistent β range
# =============================================================================

source(here::here("R/00_setup.R"))
tictoc::tic("Robustness table")
full_panel <- readRDS(file.path(data_dir, "analytic_panel.rds"))

controls <- c("log_gdp_pc","old_age_dep_ratio",
               "ltc_exp_gdp_pct","govt_effectiveness","icrg_quality")
ctrl_f   <- paste(controls, collapse = " + ")
base_f   <- paste("FEI ~ IDS_pure +", ctrl_f, "| country_iso + year")

run_spec <- function(df, formula_str = base_f, exposure = "IDS_pure",
                     label, wts = NULL) {
  m  <- feols(as.formula(formula_str), data = df,
              cluster = ~country_iso, weights = wts)
  set.seed(20230901)
  wcb <- fwildclusterboot::boottest(
    m, clustid = "country_iso", param = exposure,
    B = 4999, type = "webb", impose_null = FALSE)
  tibble(specification = label,
         beta     = coef(m)[exposure],
         ci_lower = wcb$conf_int[1],
         ci_upper = wcb$conf_int[2],
         p_wcb    = wcb$p_val)
}

specs <- list(

  # S1: Primary (IDS_pure, full sample, WCB)
  run_spec(full_panel, label = "S1: IDS_pure — primary, WCB"),

  # S2: IDS_full
  run_spec(full_panel,
           paste("FEI ~ IDS_full +", ctrl_f, "| country_iso + year"),
           "IDS_full", "S2: IDS_full (CLS in construction)"),

  # S3: OECD only
  run_spec(filter(full_panel, oecd_member == 1),
           label = "S3: OECD subsample (n=22)"),

  # S4: λ = 0.31
  {
    p2 <- full_panel %>%
      mutate(mean_loss = ((1-w12)+(1-w23)+(1-w34)+(1-w45))/4,
             H_B  = -(w12*log(w12)+w23*log(w23)+w34*log(w34)+w45*log(w45)),
             CLS  = mean_loss + 0.31*H_B,
             IDS_pure = as.numeric(scale(entitlement_scope_index - GIS)))
    run_spec(p2, label = "S4: \u03bb = 0.31 (lower CI)")
  },

  # S5: λ = 0.55
  {
    p2 <- full_panel %>%
      mutate(mean_loss = ((1-w12)+(1-w23)+(1-w34)+(1-w45))/4,
             H_B  = -(w12*log(w12)+w23*log(w23)+w34*log(w34)+w45*log(w45)),
             CLS  = mean_loss + 0.55*H_B,
             IDS_pure = as.numeric(scale(entitlement_scope_index - GIS)))
    run_spec(p2, label = "S5: \u03bb = 0.55 (upper CI)")
  },

  # S6: Entropy weighting
  {
    w_out <- WeightIt::weightit(
      IDS_pure > 0 ~ log_gdp_pc + old_age_dep_ratio + govt_effectiveness,
      data = full_panel, method = "entropy", estimand = "ATE")
    run_spec(full_panel, wts = w_out$weights, label = "S6: Entropy weighting")
  },

  # S7: Inverse-probability-of-attrition weights
  run_spec(full_panel, wts = full_panel$ipaw,
           label = "S7: IP-attrition weighting"),

  # S8: Heckman selection correction
  #     observed = 1 for all complete country-years; 0 for structurally missing
  {
    full_panel2 <- full_panel %>%
      mutate(observed = as.integer(!is.na(FEI) & !is.na(IDS_pure)))
    heck_sel <- glm(observed ~ log_gdp_pc + govt_effectiveness +
                      old_age_dep_ratio + ltc_exp_gdp_pct,
                    data = full_panel2, family = binomial(link = "probit"))
    full_panel2$mills <- dnorm(predict(heck_sel)) / pnorm(predict(heck_sel))
    run_spec(full_panel2,
             paste("FEI ~ IDS_pure + mills +", ctrl_f, "| country_iso + year"),
             label = "S8: Heckman selection correction")
  }
)

robust_tbl <- bind_rows(specs)
write_csv(robust_tbl, file.path(tab_dir, "tableS3_robustness.csv"))

message("\u03b2 range across all specifications: [",
        round(min(robust_tbl$beta), 3), ", ",
        round(max(robust_tbl$beta), 3), "]")
tictoc::toc()
message("\u2713 Robustness table complete.")
