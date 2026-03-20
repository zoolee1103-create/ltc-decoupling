# =============================================================================
# 07_goodman_bacon.R
# Goodman-Bacon (2021) decomposition
# Confirms all 2×2 comparison pairs carry positive weights
# =============================================================================

source("R/00_setup.R")
full_panel <- readRDS(file.path(data_dir, "analytic_panel.rds"))

# Goodman-Bacon requires staggered binary treatment indicator
# Use post_reform (= 1 in years >= reform_year, else 0)

bacon_result <- bacon(
  formula = FEI ~ post_reform,
  data    = full_panel %>%
              group_by(country_iso) %>%
              mutate(first_treat = min(year[post_reform == 1], na.rm = TRUE)) %>%
              ungroup() %>%
              mutate(first_treat = replace_na(first_treat, 0)),
  id_var  = "country_iso",
  time_var = "year",
  quietly  = FALSE
)

message("Goodman-Bacon weighted average: ", round(attr(bacon_result, "2x2"), 3))
message("(Paper reports weighted average = 0.141, SE = 0.028)")

# Check all weights positive
n_negative_weights <- sum(bacon_result$weight < 0)
message("Number of 2×2 pairs with negative weights: ", n_negative_weights,
        " (should be 0)")

# --- Plot (Fig. S3) ---
fig_s3 <- ggplot(bacon_result, aes(x = weight, y = estimate,
                                    colour = type, shape = type)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
  scale_colour_brewer(palette = "Set1", name = "Comparison type") +
  scale_shape_manual(values = c(16,17,15,18), name = "Comparison type") +
  labs(x = "Weight", y = "2×2 DiD estimate",
       title = "Fig. S3 | Goodman-Bacon decomposition",
       subtitle = paste("Weighted average =",
                        round(attr(bacon_result, "2x2"), 3))) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "figS3_goodman_bacon.pdf"),
       fig_s3, width = 7, height = 5, device = cairo_pdf)
message("✓ Goodman-Bacon decomposition complete. Fig. S3 saved.")
