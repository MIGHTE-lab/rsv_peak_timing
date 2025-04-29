## -----------------------------------------------------------------------------
## Script name: onset_peak_statistics.R
##
## Purpose of script: Calculate statistics for EWS-generated onset and peak
## dates
##
## Author: George Dewey
##
## Date Created: 2025-04-24
##
## Last Updated: 2025-04-24
## -----------------------------------------------------------------------------

# Load packages
library(tidyverse)
library(lubridate)

# Load data
flu_dat = read_csv('early_warning/flu_matched_peaks_onsets.csv')
rsv_dat = read_csv('early_warning/rsv_matched_peaks_onsets.csv')

# Some data wrangling
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

# Onsets
## Mean

