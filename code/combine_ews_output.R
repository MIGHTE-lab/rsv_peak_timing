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

onsets1 = read_csv("early_warning/windows/flu_times_coef.csv") %>% select(-...1)
onsets1
onsets2 = read_csv("early_warning/windows_/flu_times_coef_040725.csv") %>% select(-...1)
onsets2

# Add the season var to onsets as well (based on start_date)
peaks1 = read_csv("early_warning/peaks/flu_times_coef.csv") %>% select(-c(...1, date_detected))
peaks1
peaks2 = read_csv("early_warning/peaks_/flu_times_coef_040725.csv") %>% select(-c(...1, date_detected))
peaks2

flu_peaks_all = bind_rows(peaks1, peaks2) %>%
  distinct()

onsets2 = onsets2 %>%
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
flu_matched = match_peaks_with_onsets(flu_peaks_all, onsets2)

write_csv(flu_matched, "early_warning/flu_matched_peaks_onsets.csv")

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

rsv_peaks1 = read_csv("early_warning/peaks_/rsv_times_coef_040725.csv")
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
