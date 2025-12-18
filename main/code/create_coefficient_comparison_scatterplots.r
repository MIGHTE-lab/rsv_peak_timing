## -----------------------------------------------------------------------------
## Script name: create_coefficient_comparison_scatterplots.R
##
## Purpose of script: Compare ridge coefficients with and without intercept (m1 vs m2)
##
## Author: George Dewey
##
## Date Created: 2025-10-01
## -----------------------------------------------------------------------------

library(tidyverse)
library(readxl)
library(ggplot2)
library(scales)

# Read the coefficient data
coef_data = read_excel("data/unified_comparison_single_sheet.xlsx")

# Clean and reshape the data for plotting
coef_data_clean = coef_data %>%
  # Remove any rows with missing data
  filter(!is.na(State), !is.na(Season)) %>%
  # Rename columns for easier handling
  rename(
    state = State,
    season = Season,
    # Model 1 (ridge no intercept)
    m1_intercept = m1_Intercept,
    m1_covid = m1_COVID,
    m1_flu = m1_Flu,
    m1_rsv = m1_RSV,
    # Model 2 (ridge with intercept)
    m2_intercept = m2_Intercept,
    m2_covid = m2_COVID,
    m2_flu = m2_Flu,
    m2_rsv = m2_RSV
  ) %>%
  select(state, season, m1_covid, m1_flu, m1_rsv, m2_covid, m2_flu, m2_rsv)

# Reshape to long format for plotting
coef_long = coef_data_clean %>%
  pivot_longer(
    cols = -c(state, season),
    names_to = c("model", "pathogen"),
    names_pattern = "(m\\d)_(\\w+)",
    values_to = "coefficient"
  ) %>%
  pivot_wider(
    names_from = model,
    values_from = coefficient
  ) %>%
  mutate(
    pathogen = case_when(
      pathogen == "covid" ~ "COVID-19",
      pathogen == "flu" ~ "Influenza",
      pathogen == "rsv" ~ "RSV"
    ),
    # Create season labels for faceting
    season_label = case_when(
      season == "22-23" ~ "2022-23",
      season == "23-24" ~ "2023-24",
      season == "24-25" ~ "2024-25"
    ),
    # Order seasons properly
    season_factor = factor(season_label, levels = c("2022-23", "2023-24", "2024-25"))
  )

# Determine the overall axis limits to make both axes identical
x_range = range(coef_long$m1, na.rm = TRUE)
y_range = range(coef_long$m2, na.rm = TRUE)

# Find the overall range that encompasses both x and y data
overall_min = min(c(x_range[1], y_range[1]))
overall_max = max(c(x_range[2], y_range[2]))

# Add a small buffer around the data
range_buffer = 0.05 * (overall_max - overall_min)
axis_limits = c(overall_min - range_buffer, overall_max + range_buffer)

# Create the scatterplot comparison
coefficient_comparison_plot = coef_long %>%
  ggplot(aes(x = m1, y = m2)) +

  # Add diagonal reference line (y = x)
  geom_abline(intercept = 0, slope = 1,
              linetype = "dashed", color = "gray60", size = 0.8) +

  # Add scatterpoints colored by pathogen
  geom_point(aes(color = pathogen, shape = pathogen),
             size = 2.5, alpha = 0.8) +

  # Customize colors and shapes based on your preferences
  scale_color_manual(
    values = c(
      "COVID-19" = "#EFB75B",  # Orange squares
      "Influenza" = "#C34129", # Red triangles
      "RSV" = "#183A5A"        # Dark blue circles
    ),
    name = ""
  ) +

  scale_shape_manual(
    values = c(
      "COVID-19" = 15,  # Square
      "Influenza" = 17, # Triangle
      "RSV" = 16        # Circle
    ),
    name = ""
  ) +

  # Facet by season with FIXED scales to ensure alignment
  facet_wrap(~ season_factor, ncol = 3, scales = "fixed") +

  # Set identical limits for both axes to make square plots
  scale_x_continuous(
    name = "Beta Coefficient (No Intercept)",
    labels = number_format(accuracy = 0.1),
    limits = axis_limits
  ) +

  scale_y_continuous(
    name = "Beta Coefficient (With Intercept)",
    labels = number_format(accuracy = 0.1),
    limits = axis_limits
  ) +

  # Force square aspect ratio
  coord_fixed(ratio = 1) +

  # Apply clean theme with additional spacing
  theme_minimal() +
  theme(
    # Text sizes
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 14, face = "bold"),

    # Add spacing between axis tick labels and axis titles
    axis.title.x = element_text(margin = margin(t = 15)),  # Add top margin to x-axis title
    axis.title.y = element_text(margin = margin(r = 15)),  # Add right margin to y-axis title

    # Legend
    legend.position = "bottom",
    legend.title = element_blank(),  # Remove legend title
    legend.text = element_text(size = 14),

    # Grid lines
    panel.grid.major = element_line(color = "gray90", size = 0.3),
    panel.grid.minor = element_line(color = "gray95", size = 0.2),

    # Add more spacing between facets
    panel.spacing = unit(2, "lines"),  # Increased from 1 to 2 lines

    # Plot margins
    plot.margin = margin(20, 20, 20, 20),

    # Ensure panels maintain square aspect ratio
    aspect.ratio = 1
  )

# Display the plot
print(coefficient_comparison_plot)

# Save the plot with adjusted dimensions to accommodate square facets
ggsave("figures/coefficient_comparison_scatterplots.png",
       coefficient_comparison_plot,
       width = 15, height = 7,  # Increased height to accommodate square facets
       units = "in",
       dpi = 300,
       bg = "white")

# # Calculate summary statistics for the comparison
# comparison_stats = coef_long %>%
#   group_by(season, pathogen) %>%
#   summarise(
#     n_states = n(),
#
#     # Mean coefficients
#     mean_no_int = round(mean(m1, na.rm = TRUE), 3),
#     mean_with_int = round(mean(m2, na.rm = TRUE), 3),
#
#     # Correlation between the two approaches
#     correlation = round(cor(m1, m2, use = "complete.obs"), 3),
#
#     # How often the with-intercept coefficient is smaller
#     pct_reduced = round(100 * sum(m2 < m1, na.rm = TRUE) / n(), 1),
#
#     # How often the with-intercept coefficient is set to zero
#     pct_zero = round(100 * sum(m2 == 0, na.rm = TRUE) / n(), 1),
#
#     # Mean absolute difference
#     mean_abs_diff = round(mean(abs(m2 - m1), na.rm = TRUE), 3),
#
#     .groups = "drop"
#   ) %>%
#   arrange(season, pathogen)
#
# # Overall summary across all seasons
# overall_stats = coef_long %>%
#   group_by(pathogen) %>%
#   summarise(
#     n_observations = n(),
#
#     # Overall correlations
#     overall_correlation = round(cor(m1, m2, use = "complete.obs"), 3),
#
#     # Overall reduction patterns
#     overall_pct_reduced = round(100 * sum(m2 < m1, na.rm = TRUE) / n(), 1),
#     overall_pct_zero = round(100 * sum(m2 == 0, na.rm = TRUE) / n(), 1),
#
#     # Mean coefficients across all seasons
#     overall_mean_no_int = round(mean(m1, na.rm = TRUE), 3),
#     overall_mean_with_int = round(mean(m2, na.rm = TRUE), 3),
#
#     .groups = "drop"
#   )
#
# # Create a summary table showing the impact of including intercept
# impact_summary = coef_long %>%
#   mutate(
#     coefficient_change = m2 - m1,
#     pct_change = ifelse(m1 != 0,
#                         100 * (m2 - m1) / m1,
#                         NA),
#     was_zeroed = (m1 > 0 & m2 == 0)
#   ) %>%
#   group_by(pathogen) %>%
#   summarise(
#     n_observations = n(),
#
#     # Mean changes
#     mean_absolute_change = round(mean(abs(coefficient_change), na.rm = TRUE), 3),
#     mean_percent_change = round(mean(pct_change, na.rm = TRUE), 1),
#
#     # Frequency of being zeroed out
#     n_zeroed = sum(was_zeroed, na.rm = TRUE),
#     pct_zeroed = round(100 * sum(was_zeroed, na.rm = TRUE) / n(), 1),
#
#     # Frequency of reduction
#     n_reduced = sum(coefficient_change < 0, na.rm = TRUE),
#     pct_reduced = round(100 * sum(coefficient_change < 0, na.rm = TRUE) / n(), 1),
#
#     .groups = "drop"
#   ) %>%
#   arrange(desc(pct_zeroed))
#
# print("\n=== COEFFICIENT COMPARISON ANALYSIS ===")
# print("\nBy Season and Pathogen:")
# print(comparison_stats)
#
# print("\nOverall Summary by Pathogen:")
# print(overall_stats)
#
# print("\nImpact of Including Intercept:")
# print(impact_summary)
#
# # Save summary tables
# write_csv(comparison_stats, "data/coefficient_comparison_by_season.csv")
# write_csv(overall_stats, "data/coefficient_comparison_overall.csv")
# write_csv(impact_summary, "data/intercept_impact_summary.csv")
