## -----------------------------------------------------------------------------
## Script name: create_state_timeline_visualization.R
##
## Purpose of script: Visualize peak timing with states as rows and seasons separated by vertical lines
##
## Author: George Dewey
##
## Date Created: 2025-10-01
## -----------------------------------------------------------------------------

library(tidyverse)
library(lubridate)
library(ggplot2)
library(scales)

# Read the peak dates data
peak_data = read_csv("data/nssp_peaks.csv")

# Convert date columns to proper date format
peak_data = peak_data %>%
  mutate(
    max_rsv_date = as_date(max_rsv_date),
    max_flu_date = as_date(max_flu_date),
    max_cov_date = as_date(max_cov_date),
    max_ili_date = as_date(max_ili_date)
  )

# Reshape data for plotting
peak_data_long = peak_data %>%
  pivot_longer(
    cols = c(max_rsv_date, max_flu_date, max_cov_date, max_ili_date),
    names_to = "pathogen_type",
    values_to = "peak_date"
  ) %>%
  mutate(
    pathogen = case_when(
      pathogen_type == "max_rsv_date" ~ "RSV",
      pathogen_type == "max_flu_date" ~ "Flu",
      pathogen_type == "max_cov_date" ~ "COVID-19",
      pathogen_type == "max_ili_date" ~ "ILI"
    )
  ) %>%
  # Order states alphabetically for consistent y-axis
  arrange(state) %>%
  mutate(
    state_factor = factor(state, levels = sort(unique(state), decreasing = TRUE))
  )

# Define ILI season boundaries (May 31 of each year)
season_boundaries = as_date(c(
  "2022-05-31",  # Start of 22-23 season
  "2023-05-31",  # End of 22-23, start of 23-24
  "2024-05-31",  # End of 23-24, start of 24-25
  "2025-05-31"   # End of 24-25 season
))

# Determine the full date range for x-axis
date_range = range(peak_data_long$peak_date, na.rm = TRUE)
x_start = floor_date(date_range[1], "month") - months(1)
x_end = ceiling_date(date_range[2], "month") + months(1)

# Create the main visualization
state_timeline_plot = peak_data_long %>%
  ggplot(aes(x = peak_date, y = state_factor)) +

  # Add vertical lines for ILI season boundaries
  geom_vline(
    xintercept = season_boundaries,
    linetype = "solid",
    color = "black",
    size = 0.8,
    alpha = 0.7
  ) +

  # Add season labels at the top with more white space
  annotate("text",
           x = c(
             mean(c(season_boundaries[1], season_boundaries[2])),
             mean(c(season_boundaries[2], season_boundaries[3])),
             mean(c(season_boundaries[3], season_boundaries[4]))
           ),
           y = length(unique(peak_data_long$state)) + 3,  # Increased from 1.5 to 3
           label = c("2022-23 Season", "2023-24 Season", "2024-25 Season"),
           size = 4,
           fontface = "bold",
           color = "black") +

  # Add different shapes and colors for each pathogen
  geom_point(aes(color = pathogen, shape = pathogen),
             size = 2.8, alpha = 0.9) +

  # Customize colors and shapes
  scale_color_manual(
    values = c(
      'ILI' = 'grey10',
      'RSV' = '#183A5A',
      'Flu' = '#C34129',
      'COVID-19' = "#EFB75B"
    ),
    name = "Peak Type"
  ) +

  scale_shape_manual(
    values = c(
      "RSV" = 16,      # Circle
      "Flu" = 17,      # Triangle
      "COVID-19" = 15, # Square
      "ILI" = 18       # Diamond
    ),
    name = "Peak Type"
  ) +

  # Customize x-axis to show the full date range
  scale_x_date(
    limits = c(x_start, x_end),
    date_labels = "%b\n%Y",
    date_breaks = "3 months",
    expand = expansion(mult = c(0.01, 0.01))
  ) +

  # Customize y-axis with state names
  scale_y_discrete(
    expand = expansion(mult = c(0.02, 0.08))  # Increased top expansion for season labels
  ) +

  # Themes and labels - removed title, subtitle, x-axis label, and caption
  labs(
    x = "",  # Removed "Date" label
    y = "State"
  ) +

  theme_minimal() +
  theme(
    # Adjust text sizes for readability
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 14, angle = 0),
    axis.title = element_text(size = 20, face = "bold"),

    # Legend positioning and styling
    legend.position = "bottom",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.box = "horizontal",

    # Grid lines - keep horizontal for state separation, remove vertical
    panel.grid.major.y = element_line(color = "gray90", size = 0.3),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_line(color = "gray95", size = 0.2),
    panel.grid.minor.x = element_blank(),

    # Plot margins
    plot.margin = margin(20, 20, 20, 20),

    # Panel background
    panel.background = element_rect(fill = "white", color = NA)
  )

# Display the plot
print(state_timeline_plot)

# Save the plot
ggsave("figures/state_timeline_peak_visualization.png",
       state_timeline_plot,
       width = 16, height = 12,
       units = "in",
       dpi = 300,
       bg = "white")

# Create enhanced summary statistics
peak_summary_enhanced = peak_data %>%
  group_by(state) %>%
  summarise(
    # Count seasons with data
    n_seasons = n(),

    # RSV-Flu timing patterns
    mean_rsv_to_flu_days = round(mean(as.numeric(max_flu_date - max_rsv_date), na.rm = TRUE), 1),

    # Consistency metrics
    rsv_always_first = all(max_rsv_date <= max_flu_date, na.rm = TRUE),
    flu_always_after_rsv = all(max_flu_date >= max_rsv_date, na.rm = TRUE),

    # COVID timing variability (coefficient of variation of day-of-year)
    covid_timing_variability = round(sd(yday(max_cov_date), na.rm = TRUE) / mean(yday(max_cov_date), na.rm = TRUE), 3),

    .groups = "drop"
  ) %>%
  arrange(state)

# Overall summary by season
seasonal_summary = peak_data %>%
  group_by(season) %>%
  summarise(
    n_states = n(),

    # RSV-Flu timing
    mean_rsv_to_flu_days = round(mean(as.numeric(max_flu_date - max_rsv_date), na.rm = TRUE), 1),
    median_rsv_to_flu_days = round(median(as.numeric(max_flu_date - max_rsv_date), na.rm = TRUE), 1),

    # Proportion where RSV peaks before flu
    pct_rsv_before_flu = round(100 * sum(max_rsv_date < max_flu_date, na.rm = TRUE) / n(), 1),

    # COVID timing spread (range in days)
    covid_peak_range_days = as.numeric(max(max_cov_date, na.rm = TRUE) - min(max_cov_date, na.rm = TRUE)),

    .groups = "drop"
  )

print("\n=== ENHANCED PEAK TIMING ANALYSIS ===")
print("\nSummary by Season:")
print(seasonal_summary)

print("\nStates with Most Consistent RSV-before-Flu Pattern:")
consistent_states = peak_summary_enhanced %>%
  filter(rsv_always_first == TRUE) %>%
  select(state, n_seasons, mean_rsv_to_flu_days)
print(consistent_states)

print("\nStates with Highest COVID Timing Variability:")
variable_covid = peak_summary_enhanced %>%
  arrange(desc(covid_timing_variability)) %>%
  head(10) %>%
  select(state, covid_timing_variability, n_seasons)
print(variable_covid)

# Save summary tables
write_csv(peak_summary_enhanced, "data/state_peak_timing_summary.csv")
write_csv(seasonal_summary, "data/seasonal_peak_timing_summary.csv")

cat("\n=== OUTPUT FILES ===\n")
cat("Main visualization: figures/state_timeline_peak_visualization.png\n")
cat("State summary: data/state_peak_timing_summary.csv\n")
cat("Seasonal summary: data/seasonal_peak_timing_summary.csv\n")