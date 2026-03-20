# =============================================================================
# 10_regime_moderation.R
# Welfare regime moderation of the IDS–FEI gradient (Fig. 4a)
#
# Welfare regime × IDS_pure interaction tested via fixest::wald() [joint
# significance] and ΔAIC [model fit]; both preferred over anova() for feols
# objects (anova() is not supported for feols).
# =============================================================================

source(here::here("R/00_setup.R"))
tictoc::tic("Regime moderation")
full_panel <- readRDS(file.path(data_dir, "analytic_panel.rds"))

controls <- c("log_gdp_pc","old_age_dep_ratio",
               "ltc_exp_gdp_pct","govt_effectiveness","icrg_quality")
ctrl_f <- paste(controls, collapse = " + ")
regimes <- c("Beveridgean statist","Bismarckian corporatist",
             "Residualist liberal","East Asian developmental")

# =============================================================================
# 1. Regime-stratified TWFE + Wild Cluster Bootstrap per regime
# =============================================================================

fit_regime <- function(r) {
  df <- filter(full_panel, welfare_regime == r)
  m  <- feols(as.formula(paste("FEI ~ IDS_pure +", ctrl_f,
                                "| country_iso + year")),
              data = df, cluster = ~country_iso)
  # WCB where G >= 6; skip for very small subgroups
  if (n_distinct(df$country_iso) >= 6) {
    set.seed(20230901)
    wcb <- fwildclusterboot::boottest(
      m, clustid = "country_iso", param = "IDS_pure",
      B = 4999, type = "webb", impose_null = FALSE)
    ci <- wcb$conf_int
  } else {
    ci <- confint(m)["IDS_pure", ]
  }
  tibble(regime = r, beta = coef(m)["IDS_pure"],
         ci_lo = ci[1], ci_hi = ci[2],
         n = n_distinct(df$country_iso))
}

regime_coefs <- map_dfr(regimes, fit_regime)
print(regime_coefs)

delta_bism_east <- diff(regime_coefs$beta[c(2, 4)])
message("\u0394\u03b2 Bismarckian \u2013 East Asian = ", round(delta_bism_east, 3))

# =============================================================================
# 2. Pooled vs interacted model + wald() for interaction significance
#    CORRECTION: use fixest::wald() not anova() for feols objects
# =============================================================================

m_pooled <- feols(as.formula(paste("FEI ~ IDS_pure +", ctrl_f,
                                    "| country_iso + year")),
                  data = full_panel, cluster = ~country_iso)

m_interacted <- feols(
  as.formula(paste("FEI ~ i(welfare_regime, IDS_pure, ref='Beveridgean statist') +",
                   ctrl_f, "| country_iso + year")),
  data = full_panel, cluster = ~country_iso
)

# Joint Wald test for heterogeneity of IDS slopes across regimes
wald_het <- fixest::wald(m_interacted,
                          keep = "welfare_regime.*IDS_pure")
message("Wald test for regime moderation: F(", wald_het$df1, ",",
        wald_het$df2, ") = ", round(wald_het$stat, 2),
        "  p = ", round(wald_het$p, 4))

# AIC comparison
delta_aic <- AIC(m_interacted) - AIC(m_pooled)
message("\u0394AIC (interacted \u2013 pooled) = ", round(delta_aic, 1),
        " (negative = better fit with interactions)")

# =============================================================================
# 3. Boundary-specific fidelity loss by regime
# =============================================================================

boundary_by_regime <- full_panel %>%
  group_by(welfare_regime) %>%
  summarise(across(c(w12,w23,w34,w45),
                   list(loss = ~mean(1-.x, na.rm=TRUE))), .groups="drop")
write_csv(boundary_by_regime,
          file.path(tab_dir, "boundary_losses_by_regime.csv"))

# =============================================================================
# 4. Fig. 4a forest plot
# =============================================================================

pooled_beta <- coef(m_pooled)["IDS_pure"]
fig4a <- ggplot(regime_coefs, aes(x = beta, y = reorder(regime, beta),
                                   colour = regime)) +
  geom_vline(xintercept = pooled_beta, linetype = "dashed",
             colour = "grey30", linewidth = 0.7) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi),
                 height = 0.2, linewidth = 0.9) +
  geom_point(size = 4) +
  geom_text(aes(label = paste0("n=", n, "  \u03b2=", round(beta,3))),
            hjust = -0.1, size = 3.1) +
  scale_colour_manual(values = regime_colours, guide = "none") +
  scale_x_continuous(limits = c(0, 0.32), breaks = seq(0, 0.30, 0.05)) +
  labs(x = "IDS\u2013FEI coefficient (\u03b2, 95% WCB CI)", y = NULL,
       title = "a",
       subtitle = paste0("\u0394AIC = ", round(delta_aic,1),
                         "  Wald F = ", round(wald_het$stat,1),
                         "  p = ", round(wald_het$p,3))) +
  theme_bw(base_size = 11) +
  theme(axis.text.y = element_text(size = 10),
        panel.grid.minor = element_blank())

saveRDS(fig4a, file.path(output_dir, "fig4a_plot.rds"))
tictoc::toc()
message("\u2713 Regime moderation complete.")
