## -----------------------------------------------------------------------------
## Script name: combine_ews_output.R
##
## Purpose of script: Combine EWS output of multiple runs
##
## Author: George Dewey
##
## Date Created: 2025-04-23
##
## Last Updated: 2025-04-23
## -----------------------------------------------------------------------------

library(tidyverse)
library(lubridate)

# Flu Data
# Goal is to make a combined dataset with both peaks and onsets

flu_onsets1 = read_csv("data/early_warning/windows/flu_times_coef.csv") %>% select(-...1)
flu_onsets2 = read_csv("data/early_warning/windows_/flu_times_coef_040725.csv") %>% select(-...1)

# Add the season var to onsets as well (based on start_date)
flu_peaks1 = read_csv("data/early_warning/peaks/flu_times_coef.csv") %>% select(-c(...1, date_detected))
flu_peaks2 = read_csv("data/early_warning/peaks_042425/flu_times_coef_041525.csv") %>% select(-c(...1, date_detected))

flu_peaks2

flu_peaks_all = bind_rows(flu_peaks1, flu_peaks2) %>%
  distinct()

flu_peaks_all = flu_peaks_all %>%
  mutate(
    year = year(date_peak),
    month = month(date_peak),
    season = ifelse(
      month >= 6,
      paste0(year %% 100, "-", (year + 1) %% 100),
      # If June or later: "YY-(YY+1)"
      paste0((year - 1) %% 100, "-", year %% 100)
    )   # If before June: "(YY-1)-YY"
  ) %>%
  select(-c(year, month))

flu_peaks_all

flu_onsets_all = bind_rows(flu_onsets1, flu_onsets2) %>% distinct()

flu_onsets_all = flu_onsets_all %>%
  mutate(
    year = year(start_date),
    month = month(start_date),
    season = ifelse(
      month >= 6,
      paste0(year %% 100, "-", (year + 1) %% 100),
      # If June or later: "YY-(YY+1)"
      paste0((year - 1) %% 100, "-", year %% 100)
    )   # If before June: "(YY-1)-YY"
  ) %>%
  select(-c(year, month))

# Function to match onsets with their nearest peaks
match_onsets_to_peaks <- function(onsets_df, peaks_df) {
  # Initialize an empty dataframe to store the matches
  matches <- tibble()

  # Iterate through each onset
  for (i in 1:nrow(onsets_df)) {
    onset_row <- onsets_df[i, ]

    # Filter peaks that are in the same state/proxy and season
    matching_peaks <- peaks_df %>%
      filter(
        name_proxy == onset_row$name_proxy,
        season == onset_row$season,
        # Ensure peak occurs after onset start date
        date_peak >= onset_row$start_date
      )

    if (nrow(matching_peaks) > 0) {
      # Calculate the time difference between onset start and each peak in weeks
      matching_peaks <- matching_peaks %>%
        mutate(time_diff = as.numeric(difftime(date_peak, onset_row$start_date, units = "weeks")))

      # Find the nearest peak
      nearest_peak <- matching_peaks %>%
        arrange(time_diff) %>%
        slice(1)

      # Create a row for the match
      match_row <- bind_cols(
        onset_row,
        nearest_peak %>% select(peak_name_location = name_location, date_peak, time_diff)
      )

      # Add to matches dataframe
      matches <- bind_rows(matches, match_row)
    }
  }

  matches_out = matches %>% select(name_location, start_date, end_date, season, date_peak, time_diff)

  return(matches_out)
}


# match_peaks_with_onsets = function(peaks_data, onsets_data) {
#   # Ensure consistent column naming
#   peaks_data = peaks_data %>%
#     rename(state = name_proxy) %>%
#     select(state, date_peak, season)
#
#   onsets_data = onsets_data %>%
#     rename(state = name_proxy, onset_date = start_date) %>%
#     select(state, onset_date, end_date, season)
#
#   # Cross join datasets by state and season
#   matched = peaks_data %>%
#     inner_join(onsets_data, by = c("state", "season"), relationship = "many-to-many") %>%
#     # Calculate time difference in days
#     mutate(
#       time_diff = as.numeric(difftime(date_peak, onset_date, units = "days")),
#       # Only consider onsets before peaks (positive time_diff)
#       valid_match = time_diff > 0
#     ) %>%
#     # Filter to keep only valid matches (peaks after onsets)
#     filter(valid_match) %>%
#     # For each peak, find the closest preceding onset
#     group_by(state, season, date_peak) %>%
#     slice_min(time_diff, n = 1) %>%
#     # For each onset, find the closest following peak
#     group_by(state, season, onset_date, end_date) %>%
#     slice_min(time_diff, n = 1) %>%
#     ungroup() %>%
#     # Convert to weeks for analysis
#     mutate(weeks_between = time_diff / 7) %>%
#     select(state, season, onset_date, end_date, date_peak, weeks_between)
#
#   return(matched)
# }

# Apply the function
flu_matched = match_onsets_to_peaks(onsets_df = flu_onsets_all,
                                    peaks_df = flu_peaks_all)
flu_matched %>%
  group_by(name_location, season) %>% tally()

flu_matched

write_csv(flu_matched, "data/early_warning/flu_matched_peaks_onsets2.csv")

# RSV data
rsv_onsets = read_csv("early_warning/windows_/rsv_times_coef_040725.csv") %>%
  select(-c(...1)) %>%
  mutate(
    year = year(start_date),
    month = month(start_date),
    season = ifelse(
      month >= 6,
      paste0(year %% 100, "-", (year + 1) %% 100),
      # If June or later: "YY-(YY+1)"
      paste0((year - 1) %% 100, "-", year %% 100)
    )   # If before June: "(YY-1)-YY"
  ) %>%
  select(-c(year, month))

rsv_peaks1 = read_csv("early_warning/peaks_042425/rsv_times_coef_041525.csv")
rsv_peaks2 = read_csv("early_warning/peaks/rsv_times_coef.csv")

rsv_peaks = bind_rows(rsv_peaks1, rsv_peaks2) %>% distinct() %>%
  select(name_proxy, date_peak) %>%
  mutate(
    year = year(date_peak),
    month = month(date_peak),
    season = ifelse(
      month >= 6,
      paste0(year %% 100, "-", (year + 1) %% 100),
      # If June or later: "YY-(YY+1)"
      paste0((year - 1) %% 100, "-", year %% 100)
    )   # If before June: "(YY-1)-YY"
  ) %>%
  select(-c(year, month))

match_rsv_peaks_with_onsets = function(peaks_data, onsets_data) {
  # Ensure consistent column naming
  peaks_data = peaks_data %>%
    rename(state = name_proxy) %>%
    select(state, date_peak, season)

  onsets_data = onsets_data %>%
    rename(state = name_proxy, onset_date = start_date) %>%
    select(state, onset_date, end_date, season)

  # Cross join datasets by state and season
  matched = peaks_data %>%
    inner_join(onsets_data, by = c("state", "season"), relationship = "many-to-many") %>%
    # Calculate time difference in days
    mutate(
      time_diff = as.numeric(difftime(date_peak, onset_date, units = "days")),
      # Only consider onsets before peaks (positive time_diff)
      valid_match = time_diff > 0
    ) %>%
    # Filter to keep only valid matches (peaks after onsets)
    filter(valid_match) %>%
    # For each peak, find the closest preceding onset
    group_by(state, season, date_peak) %>%
    slice_min(time_diff, n = 1) %>%
    # For each onset, find the closest following peak
    group_by(state, season, onset_date, end_date) %>%
    slice_min(time_diff, n = 1) %>%
    ungroup() %>%
    # Convert to weeks for analysis
    mutate(weeks_between = time_diff / 7) %>%
    select(state, season, onset_date, end_date, date_peak, weeks_between)

  return(matched)
}

# Apply the function
rsv_matched = match_rsv_peaks_with_onsets(rsv_peaks, rsv_onsets)

rsv_matched = rsv_matched %>% distinct()

write_csv(rsv_matched, file = "early_warning/rsv_matched_peaks_onsets.csv")

# COVID data
# Flu Data
# Goal is to make a combined dataset with both peaks and onsets

onsets_c1 = read_csv("data/early_warning/windows/covid_times_coef.csv") %>% select(-...1)
onsets_c2 = read_csv("data/early_warning/windows_/covid_times_coef_040725.csv") %>% select(-...1)
# Add the season var to onsets as well (based on start_date)
peaks_c1 = read_csv("data/early_warning/peaks/covid_times_coef.csv") %>% select(-c(...1, date_detected))
peaks_c2 = read_csv("data/early_warning/peaks_042425/covid_times_coef_041525.csv") %>% select(-c(...1, date_detected))

covid_peaks_all = bind_rows(peaks_c1, peaks_c2) %>% distinct()

covid_onsets_all = bind_rows(onsets_c1, onsets_c2) %>% distinct()

onsets_c = covid_onsets_all %>%
  mutate(
    year = year(start_date),
    month = month(start_date),
    season = ifelse(
      month >= 6,
      paste0(year %% 100, "-", (year + 1) %% 100),
      # If June or later: "YY-(YY+1)"
      paste0((year - 1) %% 100, "-", year %% 100)
    )   # If before June: "(YY-1)-YY"
  ) %>%
  select(-c(year, month))

peaks_c = covid_peaks_all %>%
  mutate(
    year = year(date_peak),
    month = month(date_peak),
    season = ifelse(
      month >= 6,
      paste0(year %% 100, "-", (year + 1) %% 100),
      # If June or later: "YY-(YY+1)"
      paste0((year - 1) %% 100, "-", year %% 100)
    )   # If before June: "(YY-1)-YY"
  ) %>%
  select(-c(year, month))

# How to show things with multiple pairs per season?
match_peaks_with_onsets = function(peaks_data, onsets_data) {
  # Ensure consistent column naming
  peaks_data = peaks_data %>%
    rename(state = name_proxy) %>%
    select(state, date_peak, season)

  onsets_data = onsets_data %>%
    rename(state = name_proxy, onset_date = start_date) %>%
    select(state, onset_date, end_date, season)

  # Cross join datasets by state and season
  matched = peaks_data %>%
    inner_join(onsets_data, by = c("state", "season"), relationship = "many-to-many") %>%
    # Calculate time difference in days
    mutate(
      time_diff = as.numeric(difftime(date_peak, onset_date, units = "days")),
      # Only consider onsets before peaks (positive time_diff)
      valid_match = time_diff > 0
    ) %>%
    # Filter to keep only valid matches (peaks after onsets)
    filter(valid_match) %>%
    # For each peak, find the closest preceding onset
    group_by(state, season, date_peak) %>%
    slice_min(time_diff, n = 1) %>%
    # For each onset, find the closest following peak
    group_by(state, season, onset_date, end_date) %>%
    slice_min(time_diff, n = 1) %>%
    ungroup() %>%
    # Convert to weeks for analysis
    mutate(weeks_between = time_diff / 7) %>%
    select(state, season, onset_date, end_date, date_peak, weeks_between)

  return(matched)
}

# Apply the function
covid_matched = match_peaks_with_onsets(peaks_c, onsets_c)
write_csv(covid_matched, file = "data/early_warning/covid_matched_peaks_onsets.csv")

# Trying the updated matching function
library(tidyverse)
library(lubridate)

# Example data (as shown in your description)
onsets_c <- tribble(
  ~name_location, ~name_proxy, ~start_date, ~end_date, ~season,
  "Alabama", "Alabama", as.Date("2023-07-22"), as.Date("2023-09-09"), "23-24",
  "Alabama", "Alabama", as.Date("2023-12-09"), as.Date("2024-01-27"), "23-24",
  "Alabama", "Alabama", as.Date("2024-07-06"), as.Date("2024-08-31"), "24-25",
  "Alabama", "Alabama", as.Date("2024-12-07"), as.Date("2025-02-15"), "24-25",
  "Alaska", "Alaska", as.Date("2022-12-10"), as.Date("2023-03-04"), "22-23",
  "Alaska", "Alaska", as.Date("2023-06-24"), as.Date("2023-09-30"), "23-24"
)

peaks_c <- tribble(
  ~name_location, ~name_proxy, ~date_peak, ~season,
  "State", "Alabama", as.Date("2023-12-30"), "23-24",
  "State", "Alabama", as.Date("2024-08-17"), "24-25",
  "State", "Alaska", as.Date("2024-01-06"), "23-24",
  "State", "Alaska", as.Date("2024-07-27"), "24-25"
)


covid_matched2 = match_onsets_to_peaks(onsets_c, peaks_c)
write_csv(covid_matched2, "data/early_warning/covid_matched_peaks_onsets_2.csv")

match_onsets_to_peaks
