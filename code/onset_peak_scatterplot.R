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

library(tidyverse)
library(lubridate)
library(patchwork)
library(ggrepel)

flu_dat = read_csv('early_warning/flu_matched_peaks_onsets.csv')
rsv_dat = read_csv('early_warning/rsv_matched_peaks_onsets.csv')

flu_dat = flu_dat %>% rename(flu_onset = onset_date,
                   flu_end = end_date,
                   flu_peak = date_peak)

flu_dat
rsv_dat = rsv_dat %>% rename(rsv_onset = onset_date,
                             rsv_end = end_date,
                             rsv_peak = date_peak)
rsv_dat

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

# Create a function to plot RSV vs flu onsets
create_onset_comparison_plot = function(flu_data, rsv_data, season) {
  # Join datasets by state and season
  merged_onsets = inner_join(
    flu_data %>% select(state, season, flu_onset),
    rsv_data %>% select(state, season, rsv_onset),
    by = c("state", "season")
  )

  merged_onsets = merged_onsets %>%
    mutate(state_abbrev = sapply(state, state_to_abbrev))

  # Create labels that include state and season
  # Check for duplicates and modify labels accordingly
  merged_onsets = merged_onsets %>%
    group_by(state, season) %>%
    mutate(
      count = n(),
      row_num = row_number(),
      label = if_else(
        count > 1,
        paste0(state_abbrev, "_", row_num),  # For duplicates, add underscore and number
        state_abbrev                          # For unique entries, just use abbreviation
      )
    ) %>%
    ungroup()

  # Create the scatterplot
  p = ggplot(merged_onsets, aes(x = rsv_onset, y = flu_onset)) +
    # Add reference line (x=y)
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    # Add points
    geom_point(aes(color = state), alpha = 0.7, size = 3) +
    # Add labels
    geom_text_repel(
      aes(label = label),
      size = 3,
      box.padding = 0.5,
      point.padding = 0.3,
      segment.color = "gray50",
      max.overlaps = 30,             # Increase max overlaps
      force = 10,                     # Increase repelling force
      force_pull = 0.5,               # Reduce attraction to points
      min.segment.length = 0.1,       # Allow shorter segments
      direction = "both",             # Allow labels in any direction
      seed = 42                       # Set seed for reproducibility
    ) +
    # Add scale breaks and formatting
    scale_x_date(date_breaks = "2 week", date_labels = "%b %d") +
    scale_y_date(date_breaks = "2 week", date_labels = "%b %d") +
    # Add labels and title
    labs(
      title = paste0("Comparison of RSV and Flu Onset Dates (",season," season)"),
      subtitle = "Points above the dashed line indicate RSV onset occurred before flu onset",
      x = "RSV Onset Date",
      y = "Flu Onset Date",
      color = "State"
    ) +
    # Use a clean theme
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",  # Hide legend as we have labels
      panel.grid.minor = element_blank()
    )


  return(p)
}

# Create the plot
onset_plot = create_onset_comparison_plot(flu_dat, rsv_dat)



# Create a function to convert state names to abbreviations


# Create a function to plot RSV vs flu peaks
create_peak_comparison_plot = function(flu_data, rsv_data, season) {
  # Join datasets by state and season
  merged_peaks = inner_join(
    flu_data %>% select(state, season, flu_peak),
    rsv_data %>% select(state, season, rsv_peak),
    by = c("state", "season")
  )

  # Add state abbreviations
  merged_peaks = merged_peaks %>%
    mutate(state_abbrev = sapply(state, state_to_abbrev))

  # Create labels that include state abbreviation
  # Check for duplicates and modify labels accordingly
  merged_peaks = merged_peaks %>%
    group_by(state, season) %>%
    mutate(
      count = n(),
      row_num = row_number(),
      label = if_else(
        count > 1,
        paste0(state_abbrev, "_", row_num),  # For duplicates, add underscore and number
        state_abbrev                          # For unique entries, just use abbreviation
      )
    ) %>%
    ungroup()

  # Create the scatterplot
  p = ggplot(merged_peaks, aes(x = rsv_peak, y = flu_peak)) +
    # Add reference line (x=y)
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
    # Add points
    geom_point(aes(color = state), alpha = 0.7, size = 3) +
    # Add labels
    geom_text_repel(
      aes(label = label),
      size = 3,
      box.padding = 0.5,
      point.padding = 0.3,
      segment.color = "gray50",
      max.overlaps = 30,             # Increase max overlaps
      force = 10,                     # Increase repelling force
      force_pull = 0.5,               # Reduce attraction to points
      min.segment.length = 0.1,       # Allow shorter segments
      direction = "both",             # Allow labels in any direction
      seed = 42                       # Set seed for reproducibility
    ) +
    # Add scale breaks and formatting
    scale_x_date(date_breaks = "2 week", date_labels = "%b %d") +
    scale_y_date(date_breaks = "2 week", date_labels = "%b %d") +
    # Add labels and title
    labs(
      title = paste0("Comparison of RSV and Flu Peak Dates (",season," season)"),
      subtitle = "Points above the dashed line indicate RSV peak occurred before flu peak",
      x = "RSV Peak Date",
      y = "Flu Peak Date",
      color = "State"
    ) +
    # Use a clean theme
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",  # Hide legend as we have labels
      panel.grid.minor = element_blank()
    )

  return(p)
}

# Create the peak comparison plot
peak_plot = create_peak_comparison_plot(flu_dat, rsv_dat)

# Lots of white space, let's do it by season
flu_dat_23 = flu_dat %>% filter(season == "23-24")
flu_dat_24 = flu_dat %>% filter(season == "24-25")
rsv_dat_23 = rsv_dat %>% filter(season == "23-24")
rsv_dat_24 = rsv_dat %>% filter(season == "24-25")

# Peak plots per season
peak_plot_23_24 = create_peak_comparison_plot(flu_dat_23, rsv_dat_23, season = "23-24")
peak_plot_24_25 = create_peak_comparison_plot(flu_dat_24, rsv_dat_24, season = "24-25")
peak_plot_23_24 + peak_plot_24_25
ggsave(file = "figures/manuscript_figures/peak_scatterplot.png", width = 14, height = 8, units = "in", bg = "white")



# Onset plots per season
onset_plot_23_24 = create_onset_comparison_plot(flu_dat_23, rsv_dat_23, season = "23-24")
onset_plot_24_25 = create_onset_comparison_plot(flu_dat_24, rsv_dat_24, season = "24-25")
onset_plot_23_24 + onset_plot_24_25
ggsave(file = "figures/manuscript_figures/onset_scatterplot.png", width = 14, height = 8, units = "in", bg = "white")
