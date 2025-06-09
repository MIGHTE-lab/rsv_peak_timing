## -----------------------------------------------------------------------------
## Script name: test_spatial_correlation.R
##
## Purpose of script: To evaluate geographic patterns and determine potential
## spatial correlations in RSV/Flu peak and onset data
##
## Author: George Dewey
##
## Date Created: 2025-06-05
##
## Last Updated: 2025-06-05
## -----------------------------------------------------------------------------

## Load packages
library(ggplot2)
library(maps)
library(dplyr)
library(viridis)
library(lubridate)
library(RColorBrewer)

## Load data
flu_dat_23 = read_csv("data/processed/flu_onset_peak_times_23.csv")
flu_dat_24 = read_csv("data/processed/flu_onset_peak_times_24.csv")
rsv_dat_23 = read_csv("data/processed/rsv_onset_peak_times_23.csv")
rsv_dat_24 = read_csv("data/processed/rsv_onset_peak_times_24.csv")
cov_dat_23 = read_csv("data/processed/cov_onset_peak_times_23.csv")
cov_dat_24 = read_csv("data/processed/cov_onset_peak_times_24.csv")

# Load state map data
states_map = map_data("state")

# Function to determine optimal time groupings based on data distribution
determine_onset_groups = function(flu_data, timing_type = "onset", n_groups = 4) {
  # Choose the timing column
  timing_col = paste0("flu_", timing_type)

  # Find the earliest date to use as week 0
  start_date = min(flu_data[[timing_col]], na.rm = TRUE)

  # Calculate weeks from start for each state
  flu_data = flu_data %>%
    mutate(
      weeks_from_start = as.numeric(difftime(.data[[timing_col]], start_date, units = "weeks")),
      week_number = round(weeks_from_start)
    )

  # Get the range of weeks
  min_week = min(flu_data$week_number, na.rm = TRUE)
  max_week = max(flu_data$week_number, na.rm = TRUE)
  week_range = max_week - min_week

  # Determine groupings based on data distribution
  if (timing_type == "onset") {
    # For onset, create meaningful seasonal groups
    breaks = c(-Inf, 4, 12, 20, Inf)
    labels = c("Very Early\n(Weeks 0-4)", "Early\n(Weeks 5-12)",
               "Mid-Season\n(Weeks 13-20)", "Late\n(Weeks 21+)")
    colors = c("#2166ac", "#5aae61", "#fdae61", "#d73027")  # Blue to red progression
  } else {
    # For peaks, use quartile-based groupings
    quartiles = quantile(flu_data$week_number, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
    breaks = c(-Inf, quartiles[2], quartiles[3], quartiles[4], Inf)
    labels = c(paste0("Earliest\n(Weeks 0-", floor(quartiles[2]), ")"),
               paste0("Early\n(Weeks ", ceiling(quartiles[2]), "-", floor(quartiles[3]), ")"),
               paste0("Late\n(Weeks ", ceiling(quartiles[3]), "-", floor(quartiles[4]), ")"),
               paste0("Latest\n(Weeks ", ceiling(quartiles[4]), "+)"))
    colors = c("#2166ac", "#5aae61", "#fdae61", "#d73027")
  }

  return(list(breaks = breaks, labels = labels, colors = colors, start_date = start_date))
}

# Create single map showing flu spread timing
create_flu_spread_timing_map = function(flu_data, timing_type = "onset") {

  # Get optimal groupings
  grouping_info = determine_onset_groups(flu_data, timing_type)

  # Choose the timing column
  timing_col = paste0("flu_", timing_type)

  # Calculate weeks and create groups
  flu_data_processed = flu_data %>%
    mutate(
      region = tolower(state),
      weeks_from_start = as.numeric(difftime(.data[[timing_col]], grouping_info$start_date, units = "weeks")),
      week_number = round(weeks_from_start),
      timing_group = cut(week_number,
                         breaks = grouping_info$breaks,
                         labels = grouping_info$labels,
                         include.lowest = TRUE),
      timing_date = .data[[timing_col]]
    ) %>%
    arrange(week_number)

  # Join with map data
  map_data = left_join(states_map, flu_data_processed, by = "region", relationship = "many-to-many")

  # Create the map
  p = ggplot(map_data, aes(x = long, y = lat, group = group)) +
    geom_polygon(aes(fill = timing_group), color = "white", linewidth = 0.4) +
    scale_fill_manual(
      values = grouping_info$colors,
      name = paste("Flu", str_to_title(timing_type), "\nTiming"),
      na.value = "lightgray"
    ) +
    coord_fixed(1.3) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(b = 10)),
      plot.subtitle = element_text(hjust = 0.5, size = 12, margin = margin(b = 15)),
      legend.position = "right",
      legend.key.size = unit(0.8, "cm"),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 12, face = "bold"),
      plot.margin = margin(20, 20, 20, 20)
    ) +
    labs(
      title = paste("Geographic Pattern of Flu", str_to_title(timing_type), "Spread: 2023-24 Season"),
      subtitle = paste("Week 0 reference:", format(grouping_info$start_date, "%B %d, %Y"),
                       "| Total spread duration:", max(flu_data_processed$week_number, na.rm = TRUE), "weeks")
    )

  # Add state labels for key early and late states
  state_centers = flu_data_processed %>%
    filter(state %in% c("Florida", "Delaware", "Connecticut")) %>%  # Example key states
    select(state, week_number)

  return(list(map = p, data = flu_data_processed, grouping = grouping_info))
}

# Create summary statistics table
create_spread_summary_table = function(flu_data_processed) {
  summary_table = flu_data_processed %>%
    group_by(timing_group) %>%
    summarise(
      n_states = n(),
      pct_states = round(n() / nrow(flu_data_processed) * 100, 1),
      week_range = paste(min(week_number, na.rm = TRUE), "-", max(week_number, na.rm = TRUE)),
      example_states = paste(head(state, 3), collapse = ", "),
      .groups = "drop"
    ) %>%
    arrange(timing_group)

  return(summary_table)
}

# Alternative version with more granular weekly coloring
create_weekly_gradient_map = function(flu_data, timing_type = "onset") {
  # Choose the timing column
  timing_col = paste0("flu_", timing_type)

  # Find the earliest date
  start_date = min(flu_data[[timing_col]], na.rm = TRUE)

  # Calculate weeks and prepare data
  flu_data_processed = flu_data %>%
    mutate(
      region = tolower(state),
      weeks_from_start = as.numeric(difftime(.data[[timing_col]], start_date, units = "weeks")),
      week_number = round(weeks_from_start)
    )

  # Join with map data
  map_data = left_join(states_map, flu_data_processed, by = "region", relationship = "many-to-many")

  # Create gradient map
  p = ggplot(map_data, aes(x = long, y = lat, group = group)) +
    geom_polygon(aes(fill = week_number), color = "white", linewidth = 0.4) +
    scale_fill_viridis_c(
      name = paste("Week of\n", str_to_title(timing_type)),
      option = "plasma",
      direction = 1,
      na.value = "lightgray",
      guide = guide_colorbar(barwidth = 1, barheight = 8)
    ) +
    coord_fixed(1.3) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(b = 10)),
      plot.subtitle = element_text(hjust = 0.5, size = 12, margin = margin(b = 15)),
      legend.position = "right",
      legend.title = element_text(size = 12, face = "bold"),
      plot.margin = margin(20, 20, 20, 20)
    ) +
    labs(
      title = paste("Flu", str_to_title(timing_type), "Timing Across US States: 2023-24"),
      subtitle = paste("Week 0 =", format(start_date, "%B %d, %Y"),
                       "| Darker colors = later", timing_type)
    )

  return(p)
}

# Usage with your flu_dat_23 data:

# 1. Create the main categorical timing map for ONSET
flu_onset_map_result = create_flu_spread_timing_map(flu_dat_23, timing_type = "onset")
print(flu_onset_map_result$map)

# 2. Create the main categorical timing map for PEAK
flu_peak_map_result = create_flu_spread_timing_map(flu_dat_23, timing_type = "peak")
print(flu_peak_map_result$map)

# 3. Create gradient version for more detailed view
flu_onset_gradient = create_weekly_gradient_map(flu_dat_23, timing_type = "onset")
print(flu_onset_gradient)

# 4. Print summary statistics
cat("ONSET TIMING SUMMARY:\n")
onset_summary = create_spread_summary_table(flu_onset_map_result$data)
print(onset_summary)

cat("\nPEAK TIMING SUMMARY:\n")
peak_summary = create_spread_summary_table(flu_peak_map_result$data)
print(peak_summary)

# 5. Show the earliest and latest states
cat("\nKEY TIMING INSIGHTS:\n")
timing_extremes = flu_onset_map_result$data %>%
  arrange(week_number) %>%
  select(state, timing_date, week_number)

cat("Earliest onset:", timing_extremes$state[1], "- Week", timing_extremes$week_number[1],
    "(", format(timing_extremes$timing_date[1], "%B %d"), ")\n")
cat("Latest onset:", timing_extremes$state[nrow(timing_extremes)], "- Week",
    timing_extremes$week_number[nrow(timing_extremes)],
    "(", format(timing_extremes$timing_date[nrow(timing_extremes)], "%B %d"), ")\n")
cat("Total spread duration:", max(timing_extremes$week_number) - min(timing_extremes$week_number), "weeks\n")