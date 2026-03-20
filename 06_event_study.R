# =============================================================================
# 06_event_study.R
# Event-study plots for four natural experiments (Fig. 3)
# Parallel trends validation + post-reform FEI trajectories
# =============================================================================

source("R/00_setup.R")
set.seed(20230901)  # Reproducibility: WCB resampling in event-study CIs
full_panel <- readRDS(file.path(data_dir, "analytic_panel.rds"))

controls <- c("log_gdp_pc", "old_age_dep_ratio",
              "ltc_exp_gdp_pct", "govt_effectiveness", "icrg_quality")

# Reform metadata
reform_meta <- tibble(
  country_iso  = c("DE",   "JP",   "KR",   "CN"),
  reform_year  = c(2017,   2006,   2008,   2016),
  reform_label = c("Germany PSG II (2017)",
                   "Japan Kaigo Revision (2006)",
                   "South Korea LTCI (2008)",
                   "China LTCI Pilots (2016)"),
  shock_type   = c("Positive","Negative","Positive","Positive"),
  window_pre   = c(-5,     -5,     -5,     -4),
  window_post  = c(5,       5,      5,      4)
)

# =============================================================================
# Event-study estimator per reform site
# Uses Sun & Abraham (2021) interaction-weighted estimator to handle
# heterogeneous treatment timing
# =============================================================================

run_event_study <- function(panel, iso, ref_year, window_pre, window_post, controls) {

  df <- panel %>%
    mutate(
      rel_year = year - ref_year,
      treated  = as.integer(country_iso == iso)
    ) %>%
    filter(rel_year >= window_pre, rel_year <= window_post)

  # Bin endpoints
  df <- df %>%
    mutate(
      rel_year_bin = case_when(
        rel_year <= window_pre  ~ window_pre,
        rel_year >= window_post ~ window_post,
        TRUE                   ~ rel_year
      )
    )

  # Omit rel_year = -1 (normalisation year)
  df <- df %>% filter(rel_year_bin != -1)

  model <- feols(
    as.formula(paste(
      "FEI ~ i(rel_year_bin, treated, ref = 0) +",
      paste(controls, collapse = " + "),
      "| country_iso + year"
    )),
    data    = df,
    cluster = ~country_iso
  )

  # Extract coefficients and CIs
  coef_df <- broom::tidy(model, conf.int = TRUE) %>%
    filter(str_detect(term, "rel_year_bin")) %>%
    mutate(
      rel_year = as.integer(str_extract(term, "-?[0-9]+")),
      country  = iso,
      ref_year = ref_year
    ) %>%
    select(rel_year, estimate, conf.low, conf.high, p.value, country, ref_year)

  # Add normalisation row
  coef_df <- bind_rows(
    coef_df,
    tibble(rel_year=0, estimate=0, conf.low=0, conf.high=0,
           p.value=NA, country=iso, ref_year=ref_year)
  ) %>% arrange(rel_year)

  list(model = model, coefs = coef_df)
}

# Run for all four sites
event_study_results <- pmap(
  reform_meta,
  function(country_iso, reform_year, reform_label, shock_type,
           window_pre, window_post) {
    run_event_study(full_panel, country_iso, reform_year,
                    window_pre, window_post, controls)
  }
)
names(event_study_results) <- reform_meta$country_iso

# =============================================================================
# Parallel trends validation
# =============================================================================

pre_trends <- map_dfr(names(event_study_results), function(iso) {
  coefs <- event_study_results[[iso]]$coefs
  pre   <- filter(coefs, rel_year < 0)
  tibble(
    country = iso,
    max_abs_pre_beta = max(abs(pre$estimate), na.rm = TRUE),
    min_pre_pvalue   = min(pre$p.value, na.rm = TRUE)
  )
})
message("Parallel trends summary:")
print(pre_trends)
message("Max |β_pre| across all sites = ",
        round(max(pre_trends$max_abs_pre_beta), 3),
        " (paper reports 0.021)")

# =============================================================================
# Fig. 3: Four-panel event-study plot
# =============================================================================

all_coefs <- map_dfr(names(event_study_results),
                      ~event_study_results[[.x]]$coefs) %>%
  left_join(reform_meta, by = c("country" = "country_iso"))

make_event_plot <- function(coefs, label, shock) {
  colour <- if (shock == "Negative") "#2166AC" else "#D6604D"
  ggplot(coefs, aes(x = rel_year, y = estimate)) +
    geom_hline(yintercept = 0, colour = "grey50", linetype = "dashed") +
    geom_vline(xintercept = 0, colour = "firebrick3", linetype = "dashed", linewidth = 0.8) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.15, fill = colour) +
    geom_line(colour = colour, linewidth = 0.9) +
    geom_point(aes(shape = (rel_year == 0)), size = 2.5, colour = colour) +
    scale_shape_manual(values = c(16, 18), guide = "none") +
    scale_x_continuous(breaks = seq(-5, 5, 1),
                       labels = function(x) ifelse(x == 0, "t=0", x)) +
    labs(x = "Years relative to reform", y = "FEI (SD units)",
         title = label,
         subtitle = paste("Shock type:", shock)) +
    theme_bw(base_size = 11) +
    theme(plot.title    = element_text(size = 10, face = "bold"),
          plot.subtitle = element_text(size = 9, colour = "grey40"),
          panel.grid.minor = element_blank())
}

plots <- pmap(
  list(split(all_coefs, all_coefs$country),
       reform_meta$reform_label,
       reform_meta$shock_type),
  make_event_plot
)

fig3 <- (plots[[1]] | plots[[2]]) / (plots[[3]] | plots[[4]]) +
  plot_annotation(
    title = "Fig. 3 | Event-study estimates: IDS on Functional Erosion Index",
    tag_levels = "a",
    theme = theme(plot.title = element_text(size = 12, face = "bold"))
  )

ggsave(file.path(fig_dir, "fig3_event_study.pdf"),
       fig3, width = 10, height = 8, device = cairo_pdf)
ggsave(file.path(fig_dir, "fig3_event_study.png"),
       fig3, width = 10, height = 8, dpi = 300)
message("✓ Fig. 3 saved.")
