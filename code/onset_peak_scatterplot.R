## -----------------------------------------------------------------------------
## Script name: onset_peak_scatterplot.R
##
## Purpose of script: Creating scatterplot of onset and peak dates
##
## Author: George Dewey
##
## Date Created: 2025-04-23
##
## Last Updated: 2025-04-23
## -----------------------------------------------------------------------------


# 1. Setup ----------------------------------------------------------------

## Load packages
library(tidyverse)
library(lubridate)
library(patchwork)
library(ggrepel)

## Load data
flu_dat = read_csv('data/early_warning/flu_matched_peaks_onsets.csv')
flu_dat2 = read_csv('data/early_warning/flu_matched_peaks_onsets2.csv')
rsv_dat = read_csv('data/early_warning/rsv_matched_peaks_onsets.csv')
cov_dat = read_csv('data/early_warning/covid_matched_peaks_onsets.csv')
cov_dat2 = read_csv('data/early_warning/covid_matched_peaks_onsets_2.csv')

## Rename columns
flu_dat = flu_dat %>% rename(flu_onset = onset_date,
                   flu_end = end_date,
                   flu_peak = date_peak)

flu_dat2 = flu_dat2 %>% rename(flu_onset = start_date,
                               flu_end = end_date,
                               flu_peak = date_peak,
                               state = name_location)

rsv_dat = rsv_dat %>% rename(rsv_onset = onset_date,
                             rsv_end = end_date,
                             rsv_peak = date_peak)

cov_dat = cov_dat %>% rename(covid_onset = onset_date,
                             covid_end = end_date,
                             covid_peak = date_peak)

cov_dat2 = cov_dat2 %>%
  rename(covid_onset = start_date,
         covid_end = end_date,
         covid_peak = date_peak,
         state = name_location) %>%
  select(state, covid_onset, covid_end, covid_peak, season, time_diff)

cov_dat2 = cov_dat2 %>% mutate(time_diff_week = round(time_diff/7, 1))

## Create the season-specific data subsets
flu_dat_23 = flu_dat2 %>% filter(season == "23-24")
flu_dat_24 = flu_dat2 %>% filter(season == "24-25")
rsv_dat_23 = rsv_dat %>% filter(season == "23-24")
rsv_dat_24 = rsv_dat %>% filter(season == "24-25")
cov_dat_23 = cov_dat2 %>% filter(season == "23-24")
cov_dat_24 = cov_dat2 %>% filter(season == "24-25")

## Helper function
state_to_abbrev = function(state_name) {
  state_abbrevs = c(
    "Alabama" = "AL", "Alaska" = "AK", "Arizona" = "AZ", "Arkansas" = "AR",
    "California" = "CA", "Colorado" = "CO", "Connecticut" = "CT", "Delaware" = "DE",
    "District of Columbia" = "DC", "Florida" = "FL", "Georgia" = "GA", "Hawaii" = "HI",
    "Idaho" = "ID", "Illinois" = "IL", "Indiana" = "IN", "Iowa" = "IA",
    "Kansas" = "KS", "Kentucky" = "KY", "Louisiana" = "LA", "Maine" = "ME",
    "Maryland" = "MD", "Massachusetts" = "MA", "Michigan" = "MI", "Minnesota" = "MN",
    "Mississippi" = "MS", "Missouri" = "MO", "Montana" = "MT", "Nebraska" = "NE",
    "Nevada" = "NV", "New Hampshire" = "NH", "New Jersey" = "NJ", "New Mexico" = "NM",
    "New York" = "NY", "North Carolina" = "NC", "North Dakota" = "ND", "Ohio" = "OH",
    "Oklahoma" = "OK", "Oregon" = "OR", "Pennsylvania" = "PA", "Rhode Island" = "RI",
    "South Carolina" = "SC", "South Dakota" = "SD", "Tennessee" = "TN", "Texas" = "TX",
    "Utah" = "UT", "Vermont" = "VT", "Virginia" = "VA", "Washington" = "WA",
    "West Virginia" = "WV", "Wisconsin" = "WI", "Wyoming" = "WY"
  )

  return(state_abbrevs[state_name])
}


# 2. Main plotting function -----------------------------------------------
# Create a function to plot flu dates vs comparison dates (RSV or COVID)
create_comparison_plot = function(flu_data, comparison_data, comparison_type, season, type = "onset", base_color = "darkslateblue") {

  # Validate the type parameter
  if(!type %in% c("onset", "peak")) {
    stop("Type must be either 'onset' or 'peak'")
  }

  # Validate the comparison_type parameter
  if(!comparison_type %in% c("rsv", "covid")) {
    stop("comparison_type must be either 'rsv' or 'covid'")
  }

  # Set column names based on type
  flu_col = paste0("flu_", type)
  comparison_col = paste0(comparison_type, "_", type)

  # Get the proper title capitalization
  type_title = paste0(toupper(substr(type, 1, 1)), substr(type, 2, nchar(type)))
  comparison_title = toupper(comparison_type)

  # Create the selection expressions
  flu_select_expr = c("state", "season", flu_col)
  comparison_select_expr = c("state", "season", comparison_col)

  # Join datasets by state and season

  flu_data %>% select(all_of(flu_select_expr))

  merged_data = inner_join(
    flu_data %>% select(all_of(flu_select_expr)),
    comparison_data %>% select(all_of(comparison_select_expr)),
    by = c("state", "season")
  )

  # Rename columns for easier access
  names(merged_data)[names(merged_data) == flu_col] = "flu_date"
  names(merged_data)[names(merged_data) == comparison_col] = "comparison_date"

  # Calculate point density for overlapping points
  # First, create rounded versions of dates to detect nearby points
  merged_data = merged_data %>%
    mutate(
      comparison_date_rounded = as.Date(format(comparison_date, "%Y-%m-%d")),
      flu_date_rounded = as.Date(format(flu_date, "%Y-%m-%d"))
    )

  # Count the number of points at each rounded position
  point_counts = merged_data %>%
    group_by(comparison_date_rounded, flu_date_rounded) %>%
    summarize(count = n(), .groups = "drop")

  # Join the counts back to the original data
  merged_data = merged_data %>%
    left_join(point_counts, by = c("comparison_date_rounded", "flu_date_rounded"))

  # Find the full date range to ensure equal scaling
  all_dates = c(merged_data$comparison_date, merged_data$flu_date)
  min_date = min(all_dates, na.rm = TRUE)
  max_date = max(all_dates, na.rm = TRUE)
  date_range = max_date - min_date

  # Add some padding to ensure square aspect and 45° line
  min_date = min_date - date_range * 0.1
  max_date = max_date + date_range * 0.1

  # Create the scatterplot
  p = ggplot(merged_data, aes(x = comparison_date, y = flu_date)) +
    # Add reference line (x=y)
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    # Add points with alpha and size based on count
    geom_point(
      aes(alpha = count, size = count),
      color = base_color  # Use the color parameter
    ) +
    # Scale alpha for better visibility of overlaps
    scale_alpha_continuous(
      range = c(0.4, 0.9),
      guide = "none"
    ) +
    # Scale size for better visibility of overlaps
    scale_size_continuous(
      range = c(2, 5),
      guide = "none"
    ) +
    # Set identical limits for both axes to ensure square aspect ratio and 45° line
    scale_x_date(
      date_breaks = "2 week",
      date_labels = "%b %d",
      limits = c(min_date, max_date)
    ) +
    scale_y_date(
      date_breaks = "2 week",
      date_labels = "%b %d",
      limits = c(min_date, max_date)
    ) +
    # Enforce square aspect ratio
    coord_fixed(ratio = 1) +
    # Add labels and title
    labs(
      title = paste0("Comparison of ", comparison_title, " and Flu ", type_title, " Dates (", season, " season)"),
      subtitle = paste0("Points above the dashed line indicate ", comparison_title, " ", type, " occurred before flu ", type),
      x = paste0(comparison_title, " ", type_title, " Date"),
      y = paste0("Flu ", type_title, " Date")
    ) +
    # Use a clean theme
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank(),
      aspect.ratio = 1          # Reinforce square aspect
    )

  return(p)
}



# # 3. Onset and peak plots for each season -------------------------------

## Peaks - Flu vs. RSV
peak_plot_flu_rsv_23_24 = create_comparison_plot(flu_dat_23, rsv_dat_23, comparison_type = "rsv", season = "23-24", type = "peak", base_color = "darkslateblue")
peak_plot_flu_rsv_23_24
ggsave(file = "figures/manuscript_figures/peak_scatterplot_flu_rsv_23_24.png", width = 10, height = 8, units = "in", bg = "white")

peak_plot_flu_rsv_24_25 = create_comparison_plot(flu_dat_24, rsv_dat_24, comparison_type = "rsv", season = "24-25", type = "peak", base_color = "dodgerblue2")
peak_plot_flu_rsv_24_25
ggsave(file = "figures/manuscript_figures/peak_scatterplot_flu_rsv_24_25.png", width = 10, height = 8, units = "in", bg = "white")

## Onsets - Flu vs. RSV
onset_plot_flu_rsv_23_24 = create_comparison_plot(flu_dat_23, rsv_dat_23, comparison_type = "rsv", season = "23-24", type = "onset", base_color = "darkslateblue")
onset_plot_flu_rsv_23_24
ggsave(file = "figures/manuscript_figures/onset_scatterplot_flu_rsv_23_24.png", width = 10, height = 8, units = "in", bg = "white")

onset_plot_flu_rsv_24_25 = create_comparison_plot(flu_dat_24, rsv_dat_24, comparison_type = "rsv", season = "24-25", type = "onset", base_color = "dodgerblue2")
onset_plot_flu_rsv_24_25
ggsave(file = "figures/manuscript_figures/onset_scatterplot_flu_rsv_24_25.png", width = 10, height = 8, units = "in", bg = "white")

## Peaks -- Flu vs. COVID
peak_plot_flu_cov_23_24 = create_comparison_plot(flu_dat_23, cov_dat_23, comparison_type = "covid", season = "23-24", type = "peak", base_color = "darkslategrey")
peak_plot_flu_cov_23_24
ggsave(file = "figures/manuscript_figures/peak_scatterplot_flu_cov_23_24.png", width = 10, height = 8, units = "in", bg = "white")
peak_plot_flu_cov_24_25 = create_comparison_plot(flu_dat_24, cov_dat_24, comparison_type = "covid", season = "24-25", type = "peak", base_color = "mediumseagreen")
peak_plot_flu_cov_24_25
ggsave(file = "figures/manuscript_figures/peak_scatterplot_flu_cov_24_25.png", width = 10, height = 8, units = "in", bg = "white")

## Onsets -- Flu vs. COVID

onset_plot_flu_cov_24_25 = create_comparison_plot(flu_dat_24, cov_dat_24, comparison_type = "covid", season = "24-25", type = "onset", base_color = "mediumseagreen")
onset_plot_flu_cov_24_25

# Onset plots per season
onset_plot_23_24 = create_onset_comparison_plot(flu_dat_23, rsv_dat_23, season = "23-24")
onset_plot_24_25 = create_onset_comparison_plot(flu_dat_24, rsv_dat_24, season = "24-25")
onset_plot_23_24
onset_plot_23_24 + onset_plot_24_25
ggsave(file = "figures/manuscript_figures/onset_scatterplot.png", width = 14, height = 8, units = "in", bg = "white")


