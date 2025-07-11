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

# 1. Setup ---------------------------------------------------------------------
## Load data
flu_dat_23 = read_csv("data/processed/flu_onset_peak_times_23.csv")
flu_dat_24 = read_csv("data/processed/flu_onset_peak_times_24.csv")
rsv_dat_23 = read_csv("data/processed/rsv_onset_peak_times_23.csv")
rsv_dat_24 = read_csv("data/processed/rsv_onset_peak_times_24.csv")
cov_dat_23 = read_csv("data/processed/cov_onset_peak_times_23.csv")
cov_dat_24 = read_csv("data/processed/cov_onset_peak_times_24.csv")

# Load required libraries
library(tidyverse)
library(maps)
library(viridis)
library(lubridate)
library(RColorBrewer)
library(patchwork)

# Load state map data
states_map = map_data("state")

# 2. Main functions ------------------------------------------------------------

# Function to determine optimal time groupings based on data distribution
determine_onset_groups = function(disease_data, timing_type = "onset", disease_type = "flu", n_groups = 4) {
  # Choose the timing column
  timing_col = paste0(disease_type, "_", timing_type)

  # Find the earliest date to use as week 0
  start_date = min(disease_data[[timing_col]], na.rm = TRUE)

  # Calculate weeks from start for each state
  disease_data = disease_data %>%
    mutate(
      weeks_from_start = as.numeric(difftime(.data[[timing_col]], start_date, units = "weeks")),
      week_number = round(weeks_from_start)
    )

  # Get the range of weeks and unique weeks
  min_week = min(disease_data$week_number, na.rm = TRUE)
  max_week = max(disease_data$week_number, na.rm = TRUE)
  unique_weeks = sort(unique(disease_data$week_number))
  week_range = max_week - min_week

  # Create robust groupings that handle clustered data
  if (timing_type == "onset") {
    # For onset, create meaningful seasonal groups
    breaks = c(-Inf, 4, 12, 20, Inf)
    labels = c("Very Early\n(Weeks 0-4)", "Early\n(Weeks 5-12)",
               "Mid-Season\n(Weeks 13-20)", "Late\n(Weeks 21+)")
    colors = c("#2166ac", "#5aae61", "#fdae61", "#d73027")
  } else {
    # For peaks, create adaptive groupings based on actual data distribution
    if (length(unique_weeks) <= 4) {
      # If we have 4 or fewer unique weeks, create one group per week
      breaks = c(-Inf, unique_weeks[-length(unique_weeks)] + 0.5, Inf)
      labels = paste0("Week ", unique_weeks)
      colors = colorRampPalette(c("#2166ac", "#d73027"))(length(unique_weeks))
    } else {
      # Use equal-interval groupings based on the actual data range
      interval_size = ceiling(week_range / 4)

      breaks = c(-Inf,
                 min_week + interval_size,
                 min_week + 2*interval_size,
                 min_week + 3*interval_size,
                 Inf)

      # Ensure breaks are unique by adjusting if necessary
      for (i in 2:(length(breaks)-1)) {
        if (breaks[i] == breaks[i-1]) {
          breaks[i] = breaks[i] + 0.5
        }
      }

      labels = c(paste0("Earliest\n(Weeks ", min_week, "-", floor(breaks[2]), ")"),
                 paste0("Early\n(Weeks ", ceiling(breaks[2]), "-", floor(breaks[3]), ")"),
                 paste0("Late\n(Weeks ", ceiling(breaks[3]), "-", floor(breaks[4]), ")"),
                 paste0("Latest\n(Weeks ", ceiling(breaks[4]), "+)"))
      colors = c("#2166ac", "#5aae61", "#fdae61", "#d73027")
    }
  }

  # Final check to ensure breaks are unique
  if (length(unique(breaks)) != length(breaks)) {
    # Final fallback: create simple equal-width intervals with guaranteed uniqueness
    breaks = c(-Inf,
               min_week + week_range/4 + 0.1,
               min_week + week_range/2 + 0.2,
               min_week + 3*week_range/4 + 0.3,
               Inf)
    labels = c("Earliest", "Early", "Late", "Latest")
    colors = c("#2166ac", "#5aae61", "#fdae61", "#d73027")
  }

  return(list(breaks = breaks, labels = labels, colors = colors, start_date = start_date))
}

determine_onset_groups(flu_dat_24, timing_type = "peak", disease_type = "flu")

# Create single map showing disease spread timing
create_disease_spread_timing_map = function(disease_data, timing_type = "onset", disease_type = "flu", show_title = TRUE) {
  disease_data = rsv_dat_23
  disease_type = "rsv"
  timing_type = "onset"
  # Get optimal groupings
  grouping_info = determine_onset_groups(disease_data, timing_type, disease_type)

  # Choose the timing column
  timing_col = paste0(disease_type, "_", timing_type)
  season = unique(disease_data$season)[1]  # Extract season from data

  # Calculate weeks and create groups
  disease_data_processed = disease_data %>%
    mutate(
      region = tolower(state),
      weeks_from_start = as.numeric(difftime(.data[[timing_col]], grouping_info$start_date, units = "weeks")),
      week_number = round(weeks_from_start)
    ) %>%
    arrange(week_number)

  # Create timing groups with error handling
  tryCatch({
    disease_data_processed = disease_data_processed %>%
      mutate(
        timing_group = cut(week_number,
                           breaks = grouping_info$breaks,
                           labels = grouping_info$labels,
                           include.lowest = TRUE),
        timing_date = .data[[timing_col]]
      )
  }, error = function(e) {
    # Fallback: create simple equal groups if cut() fails
    disease_data_processed <<- disease_data_processed %>%
      mutate(
        timing_group = case_when(
          week_number <= quantile(week_number, 0.25, na.rm = TRUE) ~ "Earliest",
          week_number <= quantile(week_number, 0.5, na.rm = TRUE) ~ "Early",
          week_number <= quantile(week_number, 0.75, na.rm = TRUE) ~ "Late",
          TRUE ~ "Latest"
        ),
        timing_date = .data[[timing_col]]
      )

    # Update grouping info for fallback
    grouping_info$colors <<- c("#2166ac", "#5aae61", "#fdae61", "#d73027")
  })
  # Join with map data
  map_data = left_join(states_map, disease_data_processed, by = "region", relationship = "many-to-many")

  # Create the map
  p = ggplot(map_data, aes(x = long, y = lat, group = group)) +
    geom_polygon(aes(fill = timing_group), color = "white", linewidth = 0.4) +
    scale_fill_manual(
      values = grouping_info$colors,
      name = paste(str_to_upper(disease_type), str_to_title(timing_type), "\nTiming"),
      na.value = "lightgray"
    ) +
    coord_fixed(1.3) +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.key.size = unit(0.8, "cm"),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 12, face = "bold"),
      plot.margin = margin(20, 20, 20, 20)
    )
  p

  # Conditionally add titles
  if (show_title) {
    p = p +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(b = 10)),
        plot.subtitle = element_text(hjust = 0.5, size = 12, margin = margin(b = 15))
      ) +
      labs(
        title = paste("Geographic Pattern of", str_to_upper(disease_type), str_to_title(timing_type), "Spread:", season, "Season"),
        subtitle = paste("Week 0 reference:", format(grouping_info$start_date, "%B %d, %Y"),
                         "| Total spread duration:", max(disease_data_processed$week_number, na.rm = TRUE), "weeks")
      )
  }

  return(list(map = p, data = disease_data_processed, grouping = grouping_info))
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
create_weekly_gradient_map = function(disease_data, timing_type = "onset", disease_type = "flu", show_title = TRUE) {
  # Choose the timing column
  timing_col = paste0(disease_type, "_", timing_type)
  season = unique(disease_data$season)[1]  # Extract season from data

  # Set fixed start date based on season and timing type
  if (timing_type == "onset") {
    # For onsets, use August 1st
    if (season == "23-24") {
      start_date = as.Date("2023-08-01")
    } else if (season == "24-25") {
      start_date = as.Date("2024-08-01")
    } else {
      # Fallback for other seasons - extract year from season and set June 1st
      year = as.numeric(paste0("20", substr(season, 1, 2)))
      start_date = as.Date(paste0(year, "-08-01"))
    }
  } else {
    # For peaks, use December 1st
    if (season == "23-24") {
      start_date = as.Date("2023-11-01")
    } else if (season == "24-25") {
      start_date = as.Date("2024-11-01")
    } else {
      # Fallback for other seasons - extract year from season and set October 1st
      year = as.numeric(paste0("20", substr(season, 1, 2)))
      start_date = as.Date(paste0(year, "-12-01"))
    }
  }

  # Calculate weeks and prepare data
  disease_data_processed = disease_data %>%
    mutate(
      region = tolower(state),
      weeks_from_start = as.numeric(difftime(.data[[timing_col]], start_date, units = "weeks")),
      week_number = round(weeks_from_start)
    )

  # Join with map data
  map_data = left_join(states_map, disease_data_processed, by = "region", relationship = "many-to-many")

  # Create gradient map with fixed scale for comparability
  # Set consistent scale limits across all disease/timing combinations
  # Assuming flu season typically spans from June to May (52 weeks max)

  # Range is smaller for peaks vs onsets - Scale depends on onsets or peaks selection
  if (timing_type == "onset") {
    limits_vec = c(0, 20)
  } else {
    limits_vec = c(2, 10)
  }

  # Different color palette for onsets or peaks
  if (timing_type == "onset") {
    palette_option = "plasma"
  } else {
    palette_option = "viridis"
  }

  # Different label based on season selection


  p = ggplot(map_data, aes(x = long, y = lat, group = group)) +
    geom_polygon(aes(fill = week_number), color = "white", linewidth = 0.4) +
    scale_fill_viridis_c(
      name = paste("Week from\n", ifelse(timing_type == "onset", "August 1st", "October 1st")),
      option = palette_option,
      direction = -1,
      na.value = "lightgray",
      limits = limits_vec,
      oob = scales::squish,
      labels = function(x) {
        # Create custom labels with "+" for maximum values
        max_val = ifelse(timing_type == "onset", 20, 10)
        ifelse(x == max_val, paste0(x, "+"), as.character(x))
      },# This line handles out-of-bounds values
      guide = guide_colorbar(barwidth = 15, barheight = 0.8, direction = "horizontal")
    ) +
    coord_fixed(1.3) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 12, face = "bold", hjust = 0),  # Added hjust = 0
      plot.margin = margin(20, 20, 20, 20)

    )

  # Conditionally add titles
  if (show_title) {
    p = p +
      theme(
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(b = 10)),
        plot.subtitle = element_text(hjust = 0.5, size = 12, margin = margin(b = 15))
      ) +
      labs(
        title = paste(str_to_upper(disease_type), str_to_title(timing_type), "Timing Across US States:", season, "Season"),
        subtitle = paste("Week 0 =", ifelse(timing_type == "onset", "August 1st", "November 1st"),
                         substr(season, 1, 2), "| Darker colors = later", timing_type)
      )
  }

  return(p)
}

# Maps for each disease/season

# Flu - 23-24 season
# # 1. Create the main categorical timing map for ONSET
# flu_onset_23_timing = create_disease_spread_timing_map(flu_dat_23, timing_type = "onset", show_title = FALSE)
# print(flu_onset_23_timing$map)
#
# # 2. Create the main categorical timing map for PEAK
# flu_peak_23_timing = create_disease_spread_timing_map(flu_dat_23, timing_type = "peak", show_title = FALSE)
# print(flu_peak_23_timing$map)

# 3. Create gradient version for ONSET
flu_onset_23_gradient = create_weekly_gradient_map(flu_dat_23, timing_type = "onset", show_title = FALSE)
print(flu_onset_23_gradient)

# 4. Check gradient version for PEAKS
flu_peak_23_gradient = create_weekly_gradient_map(flu_dat_23, timing_type = "peak", show_title = FALSE)
print(flu_peak_23_gradient)

# Flu - 24-25 season
# # 1. Create the main categorical timing map for ONSET
# flu_onset_24_timing = create_disease_spread_timing_map(flu_dat_24, timing_type = "onset", show_title = FALSE)
# print(flu_onset_24_timing$map)
#
# # 2. Create the main categorical timing map for PEAK
# flu_peak_24_timing = create_disease_spread_timing_map(flu_dat_24, timing_type = "peak", show_title = FALSE)
# print(flu_peak_24_timing$map)

# 3. Onsets gradient
flu_onset_24_gradient = create_weekly_gradient_map(flu_dat_24, timing_type = "onset", show_title = FALSE)
print(flu_onset_24_gradient)

# 4. Peaks gradient
flu_peak_24_gradient = create_weekly_gradient_map(flu_dat_24, timing_type = "peak", show_title = FALSE)
print(flu_peak_24_gradient)

# RSV - 23-24 season
# # Onset
# rsv_onset_23_timing = create_disease_spread_timing_map(rsv_dat_23, disease_type = "rsv", timing_type = "onset", show_title = FALSE)
# print(rsv_onset_23_timing$map)
#
# # Peak
# rsv_peak_23_timing = create_disease_spread_timing_map(rsv_dat_23, disease_type = "rsv", timing_type = "peak", show_title = FALSE)
# print(rsv_peak_23_timing$map)

# Onset gradient
rsv_onset_23_gradient = create_weekly_gradient_map(rsv_dat_23, disease_type = "rsv", timing_type = "onset", show_title = FALSE)
print(rsv_onset_23_gradient)

# Peak gradient
rsv_peak_23_gradient = create_weekly_gradient_map(rsv_dat_23, disease_type = "rsv", timing_type = "peak", show_title = FALSE)
print(rsv_peak_23_gradient)

# RSV - 24-25 season
# Onset
# rsv_onset_24_timing = create_disease_spread_timing_map(rsv_dat_24, disease_type = "rsv", timing_type = "onset", show_title = FALSE)
# print(rsv_onset_24_timing$map)
#
# # Peak
# rsv_peak_24_timing = create_disease_spread_timing_map(rsv_dat_24, disease_type = "rsv", timing_type = "peak", show_title = FALSE)
# print(rsv_peak_24_timing$map)

# Onset gradient
rsv_onset_24_gradient = create_weekly_gradient_map(rsv_dat_24, disease_type = "rsv", timing_type = "onset", show_title = FALSE)
print(rsv_onset_24_gradient)

# Peak gradient
rsv_peak_24_gradient = create_weekly_gradient_map(rsv_dat_24, disease_type = "rsv", timing_type = "peak", show_title = FALSE)
print(rsv_peak_24_gradient)

## Next step is to combine all these together into a pub-ready version (and share and discuss)

## For direct comparison, we can put peaks/onsets for flu+rsv in one season on the same figure since they will have the same start point

# All peaks
(flu_peak_23_gradient + rsv_peak_23_gradient) /
  (flu_peak_24_gradient + rsv_peak_24_gradient) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = list(c("Flu (23-24)", "RSV (23-24)", "Flu (24-25)", "RSV (24-25)")),
                  title = "Peaks") &
  theme(plot.tag = element_text(size = 10),
        legend.position = "bottom",
        legend.direction = "horizontal",
        plot.margin = margin(1.5, 1.5, 1.5, 1.5))
ggsave(file = "figures/manuscript_figures/supplementary/peaks_gradient_all.png", width = 10, height = 6, units = "in", bg = "white")

# All onsets
(flu_onset_23_gradient + rsv_onset_23_gradient) /
  (flu_onset_24_gradient + rsv_onset_24_gradient) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = list(c("Flu (23-24)", "RSV (23-24)", "Flu (24-25)", "RSV (24-25)")),
                  title = "Onsets") &
  theme(plot.tag = element_text(size = 10),
        legend.position = "bottom",
        legend.direction = "horizontal",
        plot.margin = margin(1.5, 1.5, 1.5, 1.5))
ggsave(file = "figures/manuscript_figures/supplementary/onsets_gradient_all.png", width = 10, height = 6, units = "in", bg = "white")

# Supplementary subfigures
# Supplementary Figure 11a - 2023-24 season onsets
onsets_gradient_23 = flu_onset_23_gradient + rsv_onset_23_gradient +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = list(c("Flu", "RSV")),
                  title = "23-24 season Onsets") &
  theme(plot.tag = element_text(size = 10),
        legend.position = "bottom",
        legend.direction = "horizontal")
onsets_gradient_23
ggsave(file = "figures/manuscript_figures/supplementary/23_24_onsets_gradient.png", width = 10, height = 4, units = "in", bg = "white")

peak_gradient_23 = flu_peak_23_gradient + rsv_peak_23_gradient +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = list(c("Flu", "RSV")),
                  title = "23-24 season Peaks") &
  theme(plot.tag = element_text(size = 10),
        legend.position = "bottom",
        legend.direction = "horizontal")
peak_gradient_23
ggsave(file = "figures/manuscript_figures/supplementary/23_24_peaks_gradient.png", width = 10, height = 4, units = "in", bg = "white")

# Supplementary Figure 11c - 2024-25 season onsets
onset_gradient_24 = flu_onset_24_gradient + rsv_onset_24_gradient +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = list(c("Flu", "RSV")),
                  title = "24-25 season Onsets") &
  theme(plot.tag = element_text(size = 10),
        legend.position = "bottom",
        legend.direction = "horizontal")
onset_gradient_24
ggsave(file = "figures/manuscript_figures/supplementary/24_25_onsets_gradient.png", width = 10, height = 4, units = "in", bg = "white")

# Supplementary Figure 11b - 2024-25 season peaks
peak_gradient_24 = flu_peak_24_gradient + rsv_peak_24_gradient +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = list(c("Flu", "RSV")),
                  title = "24-25 season Peaks") &
  theme(plot.tag = element_text(size = 10),
        legend.position = "bottom",
        legend.direction = "horizontal")
peak_gradient_24
ggsave(file = "figures/manuscript_figures/supplementary/24_25_peaks_gradient.png", width = 10, height = 4, units = "in", bg = "white")




# ## Combined figures
# ## 23-24 Flu
# (flu_onset_23_timing$map / flu_peak_23_timing$map) +
#   plot_annotation(title = "2023-24 Flu Onsets and Peak Timing")
# ggsave("figures/manuscript_figures/supplementary/flu_23_24_onsets_peak_map.png", height = 10, width = 8, units = "in", bg = "white")
#
# ## 23-24 Flu gradient
# (flu_onset_23_gradient / flu_peak_23_gradient) +
#   plot_annotation(title = "2023-24 Flu Onsets and Peak Timing - Gradient")
# getwd()
# ## 23-24 RSV
# (rsv_onset_23_timing$map / rsv_peak_23_timing$map) +
#   plot_annotation(title = "2023-24 RSV Onsets and Peak Timing")
# ggsave("figures/manuscript_figures/supplementary/rsv_23_24_onsets_peak_map.png", height = 10, width = 8, units = "in", bg = "white")
#
# ## 23-24 RSV gradient
# (rsv_onset_23_gradient / rsv_peak_23_gradient) +
#   plot_annotation(title = "2023-24 RSV Onsets and Peak Timing - Gradient")
#
# ## 24-25 Flu
# (flu_onset_24_timing$map / flu_peak_24_timing$map) +
#   plot_annotation(title = "2024-25 Flu Onsets and Peak Timing")
# ggsave("figures/manuscript_figures/supplementary/flu_24_25_onsets_peak_map.png", height = 10, width = 8, units = "in", bg = "white")
#
# ## 24-25 Flu gradient
# (flu_onset_24_gradient / flu_peaks_24_gradient) +
#   plot_annotation(title = "2024-25 Flu Onsets and Peak Timing - Gradient")
#
# ## 24-25 RSV
# (rsv_onset_24_timing$map / rsv_peak_24_timing$map) +
#   plot_annotation(title = "2024-25 RSV Onsets and Peak Timing")
# ggsave("figures/manuscript_figures/supplementary/rsv_24_25_onsets_peak_map.png", height = 10, width = 8, units = "in", bg = "white")
#
# ## 24-25 RSV gradient
# (rsv_onset_24_gradient / rsv_peak_24_gradient) +
#   plot_annotation(title = "2024-25 RSV Onsets and Peak Timing - Gradient")

# # 4. Print summary statistics
# cat("ONSET TIMING SUMMARY:\n")
# onset_summary = create_spread_summary_table(flu_onset_map_result$data)
# print(onset_summary)
#
# cat("\nPEAK TIMING SUMMARY:\n")
# peak_summary = create_spread_summary_table(flu_peak_map_result$data)
# print(peak_summary)
#
# # 5. Show the earliest and latest states
# cat("\nKEY TIMING INSIGHTS:\n")
# timing_extremes = flu_onset_map_result$data %>%
#   arrange(week_number) %>%
#   select(state, timing_date, week_number)
#
# cat("Earliest onset:", timing_extremes$state[1], "- Week", timing_extremes$week_number[1],
#     "(", format(timing_extremes$timing_date[1], "%B %d"), ")\n")
# cat("Latest onset:", timing_extremes$state[nrow(timing_extremes)], "- Week",
#     timing_extremes$week_number[nrow(timing_extremes)],
#     "(", format(timing_extremes$timing_date[nrow(timing_extremes)], "%B %d"), ")\n")
# cat("Total spread duration:", max(timing_extremes$week_number) - min(timing_extremes$week_number), "weeks\n")