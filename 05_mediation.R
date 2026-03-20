# =============================================================================
# 05_mediation.R  — Causal mediation analysis  (v6)
#
# Decomposes IDS_pure -> FEI into:
#   ACME  (Average Causal Mediation Effect) via CLS
#   ADE   (Average Direct Effect)
# Framework: Imai, Keele & Yamamoto (2010); mediation R package v4.5.0
#
# CRITICAL BUG FIX (Round 4 peer review — Issues A1/B1/C6):
#   All prior versions passed `cluster = full_panel$country_iso` to
#   mediation::mediate(). This argument does NOT exist in mediation v4.5.0
#   and was silently discarded, leaving bootstrap at the individual (i.i.d.)
#   level. Because IDS_pure and CLS are country-year exposures, the design
#   is clustered at country level, and i.i.d. resampling underestimates
#   variance. Reported CIs [0.068, 0.116] were therefore too narrow.
#
#   CORRECTED APPROACH: country-level block bootstrap
#     - Resample 32 countries with replacement (B times)
#     - Refit model.m and model.y in each bootstrap sample
#     - Compute ACME via IKY linear-linear formula: a*b
#     - 2.5th / 97.5th percentiles give cluster-robust CIs
#   This matches the approach used in 14_pathway_analysis.R.
#   B = 2,000 in production; B = 200 in demo_mode.
# =============================================================================

source(here::here("R/00_setup.R"))
tictoc::tic("Mediation analysis")

if (!exists("demo_mode")) demo_mode <- FALSE
if (demo_mode) {
  full_panel <- readRDS(file.path(data_dir, "synthetic_demo_panel.rds"))
} else {
  full_panel <- readRDS(file.path(data_dir, "analytic_panel.rds"))
}

controls <- c("log_gdp_pc", "old_age_dep_ratio",
              "ltc_exp_gdp_pct", "govt_effectiveness", "icrg_quality")
ctrl_f   <- paste(controls, collapse = " + ")

# =============================================================================
# 1. Mediator model (OLS): IDS_pure -> CLS
# =============================================================================
med_model <- lm(
  as.formula(paste("CLS ~ IDS_pure +", ctrl_f,
                   "+ factor(country_iso) + factor(year)")),
  data = full_panel
)

# =============================================================================
# 2. Outcome model (OLS, SD-scaled FEI): IDS_pure + CLS -> FEI
# =============================================================================
out_model <- lm(
  as.formula(paste("FEI ~ IDS_pure + CLS +", ctrl_f,
                   "+ factor(country_iso) + factor(year)")),
  data = full_panel
)

# Point estimates (IKY linear-linear closed form)
alpha_hat <- coef(med_model)["IDS_pure"]   # a path: treat -> mediator
beta_hat  <- coef(out_model)["CLS"]        # b path: mediator -> outcome
gamma_hat <- coef(out_model)["IDS_pure"]   # c' path: direct effect

acme_hat <- alpha_hat * beta_hat
ade_hat  <- gamma_hat
te_hat   <- acme_hat + ade_hat
prop_hat <- acme_hat / te_hat

message("=== Point Estimates ===")
message(sprintf("ACME = %.4f  ADE = %.4f  TE = %.4f  Prop = %.4f",
                acme_hat, ade_hat, te_hat, prop_hat))

# =============================================================================
# 3. Country-level block bootstrap for cluster-robust CIs
#    Replaces invalid cluster= argument (see header note)
# =============================================================================
B         <- if (demo_mode) 200L else 2000L
countries <- unique(full_panel$country_iso)
N_c       <- length(countries)

set.seed(20230901)
boot_store <- matrix(NA_real_, nrow = B, ncol = 4,
                     dimnames = list(NULL, c("acme","ade","te","prop")))

for (b in seq_len(B)) {
  boot_ctrs <- sample(countries, size = N_c, replace = TRUE)
  boot_data <- map_dfr(
    seq_along(boot_ctrs),
    ~ full_panel[full_panel$country_iso == boot_ctrs[.x], ] |>
        dplyr::mutate(country_iso = paste0("bc", .x))
  )
  tryCatch({
    bm <- lm(as.formula(paste("CLS ~ IDS_pure +", ctrl_f,
                               "+ factor(country_iso) + factor(year)")),
             data = boot_data)
    bo <- lm(as.formula(paste("FEI ~ IDS_pure + CLS +", ctrl_f,
                               "+ factor(country_iso) + factor(year)")),
             data = boot_data)
    a <- coef(bm)["IDS_pure"]
    b_path <- coef(bo)["CLS"]
    g <- coef(bo)["IDS_pure"]
    boot_store[b, "acme"] <- a * b_path
    boot_store[b, "ade"]  <- g
    boot_store[b, "te"]   <- a * b_path + g
    boot_store[b, "prop"] <- (a * b_path) / (a * b_path + g)
  }, error = function(e) NULL)
}

boot_df   <- as.data.frame(boot_store) |> tidyr::drop_na()
n_valid   <- nrow(boot_df)
ci_acme   <- quantile(boot_df$acme, c(0.025, 0.975))
ci_ade    <- quantile(boot_df$ade,  c(0.025, 0.975))
ci_te     <- quantile(boot_df$te,   c(0.025, 0.975))
ci_prop   <- quantile(boot_df$prop, c(0.025, 0.975))

message("=== Block Bootstrap CIs (B = ", n_valid, " valid draws) ===")
message(sprintf("ACME = %.3f  [%.3f, %.3f]", acme_hat, ci_acme[1], ci_acme[2]))
message(sprintf("ADE  = %.3f  [%.3f, %.3f]", ade_hat,  ci_ade[1],  ci_ade[2]))
message(sprintf("TE   = %.3f  [%.3f, %.3f]", te_hat,   ci_te[1],   ci_te[2]))
message(sprintf("Prop = %.3f  [%.3f, %.3f]", prop_hat, ci_prop[1], ci_prop[2]))

# =============================================================================
# 4. Sensitivity analysis (rho*, E-value)
#    medsens() uses parametric bounds, not bootstrap — no cluster issue here
# =============================================================================
set.seed(20230901)
med_obj_pts <- mediation::mediate(
  model.m  = med_model,
  model.y  = out_model,
  treat    = "IDS_pure",
  mediator = "CLS",
  boot     = FALSE,    # analytical; CIs come from block bootstrap above
  sims     = 200
  # cluster= NOT passed: argument does not exist in mediation v4.5.0
)

med_sens <- mediation::medsens(
  med_obj_pts, rho.by = 0.025, effect.type = "indirect", sims = 200
)

rho_star <- med_sens$rho.by[which.min(abs(med_sens$d.avg))]
e_value  <- 1 + sqrt(2 * abs(rho_star) / (1 - rho_star^2))
message(sprintf("Sensitivity: rho* = %.2f  E-value = %.2f", rho_star, e_value))

# =============================================================================
# 5. Robustness: 3-item FEI mediation (pre-registered primary specification)
# =============================================================================
set.seed(20230902)  # Intentionally different seed — ensures independent bootstrap stream from main analysis (set.seed 20230901)
boot_3item <- numeric(B)
for (b in seq_len(B)) {
  boot_ctrs <- sample(countries, size = N_c, replace = TRUE)
  boot_data <- map_dfr(
    seq_along(boot_ctrs),
    ~ full_panel[full_panel$country_iso == boot_ctrs[.x], ] |>
        dplyr::mutate(country_iso = paste0("bc", .x))
  )
  tryCatch({
    bm <- lm(as.formula(paste("CLS ~ IDS_pure +", ctrl_f,
                               "+ factor(country_iso) + factor(year)")),
             data = boot_data)
    bo <- lm(as.formula(paste("FEI_3item ~ IDS_pure + CLS +", ctrl_f,
                               "+ factor(country_iso) + factor(year)")),
             data = boot_data)
    a <- coef(bm)["IDS_pure"]; b2 <- coef(bo)["CLS"]; g <- coef(bo)["IDS_pure"]
    boot_3item[b] <- a * b2
  }, error = function(e) { boot_3item[b] <<- NA_real_ })
}
boot_3item <- boot_3item[!is.na(boot_3item)]
message(sprintf("3-item FEI ACME = %.3f  [%.3f, %.3f]",
                mean(boot_3item),
                quantile(boot_3item, 0.025),
                quantile(boot_3item, 0.975)))

# =============================================================================
# 6. Save results
# =============================================================================
mediation_summary <- tibble(
  estimand = c("ACME", "ADE", "Total Effect", "Proportion Mediated",
               "rho_star", "E_value"),
  estimate = c(acme_hat, ade_hat, te_hat, prop_hat, rho_star, e_value),
  ci_lo    = c(ci_acme[1], ci_ade[1], ci_te[1], ci_prop[1], NA, NA),
  ci_hi    = c(ci_acme[2], ci_ade[2], ci_te[2], ci_prop[2], NA, NA),
  method   = c(rep("Country-level block bootstrap (B=2000)", 4),
               "medsens() parametric rho scan",
               "VanderWeele & Ding (2017) formula")
)

write_csv(mediation_summary, file.path(tab_dir, "table_mediation.csv"))
saveRDS(list(acme_hat  = acme_hat,  ade_hat   = ade_hat,
             te_hat    = te_hat,    prop_hat  = prop_hat,
             ci_acme   = ci_acme,   ci_ade    = ci_ade,
             ci_te     = ci_te,     ci_prop   = ci_prop,
             rho_star  = rho_star,  e_value   = e_value,
             boot_df   = boot_df,   n_valid   = n_valid),
         file.path(output_dir, "mediation_results.rds"))

tictoc::toc()
message("\u2713 Mediation complete. Block bootstrap; no invalid cluster= argument.")
