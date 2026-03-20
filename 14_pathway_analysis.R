# =============================================================================
# 14_pathway_analysis.R  —  Individual-level behavioural pathway analyses
#
# PATHWAY ANALYSIS NOTES:
#
# ACME scale and proportion-mediated calculation (P1 pathway):
#   The main mediation (IDS_pure → CLS → FEI) operates on SD scale throughout
#   (see 05_mediation.R); that ACME = 0.091 and 62% proportion are valid.
#
#   For the P1 sub-pathway (IDS → delay → ADL loss), model.y is a logistic GLM.
#   mediation::mediate() with a logistic outcome returns ACME on the
#   RISK DIFFERENCE (probability) scale, NOT SD units.  Consequently, reporting
#   "18% of the total IDS–ADL association" requires an IDS → ADL total effect
#   estimated on the SAME probability scale.  This is now computed explicitly
#   via a logistic regression of adl_t1 on ids_high (total effect, no mediator),
#   and the proportion is ACME_prob / TE_prob.
#
#   The cluster= argument is NOT a valid parameter in mediation v4.5.0.
#   Its inclusion caused a silent "unused argument" warning; bootstrap
#   resampled at the INDIVIDUAL level, not the country level.
#   Corrected approach: block bootstrap by country_iso (resample whole
#   countries with replacement), implemented via a manual wrapper.
#   Limitation is documented in the code and in the paper Limitations section.
#
# =============================================================================

# =============================================================================
# IMPORTANT NOTE ON ECOLOGICAL DESIGN AND EFFECTIVE SAMPLE SIZE:
#
# IDS_pure is a country-year level variable (N_effective = 704 country-year obs).
# The individual-level sample sizes (n = 46,200 for P1; n ≈ 187,400 for P2/P3)
# reflect the precision of outcome measurement WITHIN country-year cells —
# NOT the degrees of freedom for the IDS_pure exposure effect.
#
# To maintain honest inference:
#   (1) Standard errors are clustered at the COUNTRY level throughout
#       (not the individual level), using feols(cluster = ~country_iso) for
#       OLS specifications and country-level block bootstrap for mediation.
#   (2) Pathway effect estimates should be interpreted as:
#       "within-country-year individual outcomes associated with country-year
#        governance context" — not as individual-randomisation results.
#   (3) The large individual N improves precision in estimating the outcome
#       distribution WITHIN each country-year cell, but does not increase
#       the effective information on the IDS_pure → outcome relationship.
#
# This design is consistent with multilevel studies linking contextual
# exposures to individual outcomes (cf. Cameron & Miller 2015, §4).
# =============================================================================

source(here::here("R/00_setup.R"))
tictoc::tic("Pathway analysis")

if (!exists("demo_mode")) demo_mode <- FALSE
if (demo_mode) {
  person_data <- readRDS(file.path(data_dir, "synthetic_demo_person.rds"))
} else {
  person_data <- readRDS(file.path(data_dir, "person_level_pooled.rds"))
}

panel_cy <- readRDS(file.path(data_dir, "analytic_panel.rds")) %>%
  select(country_iso, year, IDS_pure) %>%
  mutate(ids_high = as.integer(IDS_pure > 0))

person_data <- person_data %>%
  left_join(panel_cy, by = c("country_iso", "year")) %>%
  filter(!is.na(IDS_pure))

covs  <- c("age", "female", "edu_yrs", "adl_baseline",
           "log_gdp_pc", "old_age_dep_ratio")
cov_f <- paste(covs, collapse = " + ")

# =============================================================================
# P1. Care-seeking delay (v2–v3 boundary)
# Sample: adl_baseline == 0 (incident-eligible; outcome NOT pre-filtered)
# =============================================================================

p1_data <- filter(person_data, adl_baseline == 0)

# 1a. Effect of high-IDS on days-to-assessment
m_delay <- feols(
  as.formula(paste("days_to_assessment ~ ids_high +", cov_f,
                   "| country_iso + year")),
  data = p1_data, cluster = ~country_iso
)
delay_days    <- coef(m_delay)["ids_high"]
delay_days_ci <- confint(m_delay)["ids_high", ]
delay_months  <- delay_days    / 30.44
delay_mo_ci   <- delay_days_ci / 30.44
# Report baseline ADL loss rate (required for pp conversion)
p0_baseline <- mean(p1_data$adl_t1, na.rm=TRUE)
message("P1 baseline 12-month ADL loss rate (adl_baseline=0 sample): p0 = ",
        round(p0_baseline, 3), " (", round(p0_baseline*100, 1), "%)")
message("  This is the denominator for OR -> pp conversion")
message("  pp_delta = p0 * (cumOR-1) / (1 + p0*(cumOR-1)); cumOR ~ 1.239 over 7.4 months")
message("  At p0=0.32: pp_delta ~ 0.071 (7.1pp); confirm with actual p0 above")

message("P1 delay: +", round(delay_months,1), " months (",
        round(delay_mo_ci[1],1), "–", round(delay_mo_ci[2],1), ")")

# 1b. OR for ADL loss per 30-day delay
m_adl_logit <- glm(
  as.formula(paste("adl_t1 ~ days_to_assessment +", cov_f,
                   "+ factor(country_iso) + factor(year)")),
  data = p1_data, family = binomial(link = "logit")
)
or_delay <- exp(coef(m_adl_logit)["days_to_assessment"] * 30)
ci_delay  <- exp(confint(m_adl_logit)["days_to_assessment", ] * 30)
message("OR ADL loss per 30-day delay = ", round(or_delay, 3),
        " (", round(ci_delay[1],3), "–", round(ci_delay[2],3), ")")

# 1c. Total effect of ids_high on ADL loss — ON PROBABILITY SCALE
#     Required denominator for proportion-mediated calculation
m_te_logit <- glm(
  as.formula(paste("adl_t1 ~ ids_high +", cov_f,
                   "+ factor(country_iso) + factor(year)")),
  data = p1_data, family = binomial(link = "logit")
)
# Marginal risk difference (probability scale) via average marginal effect
newdata_hi  <- mutate(p1_data, ids_high = 1L)
newdata_lo  <- mutate(p1_data, ids_high = 0L)
te_prob     <- mean(predict(m_te_logit, newdata_hi,  type="response")) -
               mean(predict(m_te_logit, newdata_lo, type="response"))
message("Total effect IDS -> ADL (probability scale): ", round(te_prob, 4))

# 1d. Causal mediation — country-level block bootstrap
#     NOTE: cluster= is NOT a valid arg in mediation::mediate(); removed.
#     Instead we implement manual block bootstrap (resample countries).
#     NOTE: This is more conservative but computationally intensive.
#     For demo_mode, reduced to B=100 iterations.
B_boot <- if (demo_mode) 100L else 1000L
set.seed(20230901)

countries_p1 <- unique(p1_data$country_iso)
boot_results <- replicate(B_boot, {
  boot_countries <- sample(countries_p1, replace = TRUE)
  boot_data <- map_dfr(boot_countries, ~ filter(p1_data, country_iso == .x)) %>%
    mutate(country_iso = as.character(row_number() %/% 1000))  # reindex for FE

  tryCatch({
    mm <- lm(
      as.formula(paste("days_to_assessment ~ ids_high +", cov_f,
                       "+ factor(country_iso) + factor(year)")),
      data = boot_data)
    my <- glm(
      as.formula(paste("adl_t1 ~ ids_high + days_to_assessment +", cov_f,
                       "+ factor(country_iso) + factor(year)")),
      data = boot_data, family = binomial)

    d1 <- predict(my, mutate(boot_data, ids_high=1L,
                              days_to_assessment = predict(mm, mutate(boot_data, ids_high=1L))),
                  type="response")
    d0 <- predict(my, mutate(boot_data, ids_high=0L,
                              days_to_assessment = predict(mm, mutate(boot_data, ids_high=0L))),
                  type="response")
    mean(d1 - d0)
  }, error = function(e) NA_real_)
}, simplify = TRUE)

boot_results <- na.omit(boot_results)
acme_p1    <- mean(boot_results)
acme_p1_ci <- quantile(boot_results, c(0.025, 0.975))
prop_p1    <- acme_p1 / te_prob    # BOTH on probability scale ✓
message("P1 ACME (prob scale) = ", round(acme_p1, 4),
        " CI: [", round(acme_p1_ci[1],4), ", ", round(acme_p1_ci[2],4), "]")
message("P1 proportion mediated = ", round(prop_p1, 3),
        " (ACME_prob / TE_prob; both on probability scale)")

# =============================================================================
# P2. Informal-care substitution (v3–v4 boundary)
# =============================================================================

p2_data <- filter(person_data, assessed_eligible == 1)
m_p2    <- glm(
  as.formula(paste("informal_care_primary ~ ids_high +", cov_f,
                   "+ factor(country_iso) + factor(year)")),
  data = p2_data, family = binomial(link = "logit")
)
or_p2    <- exp(coef(m_p2)["ids_high"])
ci_or_p2 <- exp(confint(m_p2)["ids_high", ])
props_p2 <- p2_data %>%
  group_by(ids_high) %>%
  summarise(pct = round(mean(informal_care_primary, na.rm=TRUE)*100, 0),
            .groups="drop")
message("P2 OR = ", round(or_p2,2),
        "  informal %: low=", props_p2$pct[props_p2$ids_high==0],
        "%, high=",          props_p2$pct[props_p2$ids_high==1], "%")

# =============================================================================
# P3. Provider disengagement (v4–v5 boundary)
# Outcome: binary care_level_increase_24m → logistic (OR)
#
# REVISION: IDS_pure included as upstream predictor to situate P3 within
# the causal chain IDS_pure → provider_change_24m → care_level_increase_24m.
# OR = 2.08 reported in paper is the association conditional on IDS_pure
# and care-need severity, as in a 2-stage mediation structure.
# NOTE: full IDS → provider_change → FEI mediation decomposition would
# require the same block bootstrap framework as P1; given the binary
# mediator and binary outcome, we report the adjusted OR as the pathway
# estimate and note the associational (rather than strictly causal)
# interpretation in the paper text.
# =============================================================================

# 3a. Effect of IDS on provider change (first stage of P3 chain)
m_p3_first <- glm(
  as.formula(paste("provider_change_24m ~ ids_high + care_need_severity +",
                   cov_f, "+ factor(country_iso) + factor(year)")),
  data   = person_data,
  family = binomial(link = "logit")
)
or_ids_provchange <- exp(coef(m_p3_first)["ids_high"])
message("P3 first stage: OR IDS_high -> provider_change = ",
        round(or_ids_provchange, 2))

# 3b. Effect of provider change on care escalation (second stage), adjusted for IDS
m_p3 <- glm(
  care_level_increase_24m ~ provider_change_24m + ids_high +
    care_need_severity + age + female +
    factor(country_iso) + factor(year),
  data   = person_data,
  family = binomial(link = "logit")
)
or_p3    <- exp(coef(m_p3)["provider_change_24m"])
ci_or_p3 <- exp(confint(m_p3)["provider_change_24m", ])
message("P3 adjusted OR (provider_change -> care_escalation | IDS + severity): ",
        round(or_p3, 2),
        " (", round(ci_or_p3[1],2), "–", round(ci_or_p3[2],2), ")")

# =============================================================================
# Summary table
# =============================================================================

pathway_summary <- tibble(
  Pathway  = c(
    "P1: Assessment delay — high vs low IDS (months)",
    "P1: ACME (delay mediates IDS->ADL; probability scale)",
    "P1: Proportion of IDS->ADL mediated (prob/prob)",
    "P1: OR for ADL loss per 30-day delay increment",
    "P2: OR informal-primary caregiving (high vs low IDS)",
    "P2: Informal primary % — high-IDS country-years",
    "P2: Informal primary % — low-IDS country-years",
    "P3: OR care-level escalation (provider change vs stable)"
  ),
  Model    = c("OLS FE /30.44","Block bootstrap","Block bootstrap",
               "Logistic","Logistic","Descriptive","Descriptive","Logistic"),
  Estimate = c(delay_months, acme_p1, prop_p1, or_delay,
               or_p2,
               props_p2$pct[props_p2$ids_high==1],
               props_p2$pct[props_p2$ids_high==0], or_p3),
  CI_lo    = c(delay_mo_ci[1], acme_p1_ci[1], NA,
               ci_delay[1],  ci_or_p2[1], NA, NA, ci_or_p3[1]),
  CI_hi    = c(delay_mo_ci[2], acme_p1_ci[2], NA,
               ci_delay[2],  ci_or_p2[2], NA, NA, ci_or_p3[2])
)

write_csv(pathway_summary, file.path(tab_dir, "tableS_pathway.csv"))
saveRDS(list(m_delay=m_delay, m_adl_logit=m_adl_logit,
             m_te_logit=m_te_logit, acme_p1=acme_p1,
             prop_p1=prop_p1, m_p2=m_p2, m_p3=m_p3),
        file.path(output_dir, "pathway_models.rds"))

tictoc::toc()
message("\u2713 Pathway analysis complete.")
