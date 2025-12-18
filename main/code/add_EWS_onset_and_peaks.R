## -----------------------------------------------------------------------------
## Script name: add_EWS_onset_and_peaks.R
##
## Purpose of script: Add EWS onset and peak dates to RSV/flu timing plots
##
## Author: George Dewey
##
## Date Created: 2025-03-03
##
## Last Updated: 2025-04-22
## using updated EWS data
##
## -----------------------------------------------------------------------------

## Load packages
library(tidyverse)
library(lubridate)
library(ggh4x)
library(patchwork)

## Load NSSP data
nssp = read_csv(file = "data/processed/nssp_all_years_041525.csv",
                show_col_types = FALSE)

## Test to see if we can regenerate the individual season plots from the
## individual years
nssp %>%
  filter(season == "22-23") %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = flu_times_coef, color = 'Flu'), linewidth = 0.65) +
  geom_line(aes(y = rsv_times_coef, color = 'RSV'), linewidth = 0.65) +
  geom_line(aes(y = covid_times_coef, color = "COVID-19"), linewidth = 0.65) +
  geom_line(aes(y = ili_rescaled, color = 'ILI'), linewidth = 0.65) +
  facet_wrap(~ state, scales = 'free', ncol = 10) +
  scale_color_manual(values = c(
    'ILI' = '#E2DBC9',
    'RSV' = '#183A5A',
    'Flu' = '#C34129',
    'COVID-19' = "#EFB75B"
  )) +
  labs(x = '', y = '') +
  scale_x_date(date_labels = "%b\n%y", date_breaks = "2 months") +
  envalysis::theme_publish() +
  theme(
    axis.text.x = element_text(size = 5),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    legend.title = element_blank(),
    legend.position = 'right'
  )
## Import looks fine

## Load the data from EWS and do some processing

## Onsets ----
flu_onsets =
  read_csv(
    file = "early_warning/windows_/flu_times_coef_040725.csv",
    show_col_types = FALSE) %>%
  select(-1)

flu_onsets = flu_onsets %>%
  select(name_proxy, start_date, end_date) %>%
  rename(state = name_proxy, flu_onset = start_date, flu_end = end_date)

rsv_onsets = read_csv(
  file = "early_warning/windows_/rsv_times_coef_040725.csv",
  show_col_types = FALSE) %>%
  select(-1)

rsv_onsets = rsv_onsets %>%
  select(name_proxy, start_date, end_date) %>%
  rename(state = name_proxy, rsv_onset = start_date, rsv_end = end_date)

flu_onsets = flu_onsets %>%
  mutate(
    year = year(flu_onset),
    month = month(flu_onset),
    season = ifelse(
      month >= 6,
      paste0(year %% 100, "-", (year + 1) %% 100),
      # If June or later: "YY-(YY+1)"
      paste0((year - 1) %% 100, "-", year %% 100)
    )   # If before June: "(YY-1)-YY"
  )

rsv_onsets = rsv_onsets %>%
  mutate(
    year = year(rsv_onset),
    month = month(rsv_onset),
    season = ifelse(
      month >= 6,
      paste0(year %% 100, "-", (year + 1) %% 100),
      # If June or later: "YY-(YY+1)"
      paste0((year - 1) %% 100, "-", year %% 100)
    )   # If before June: "(YY-1)-YY"
  )

flu_onsets = flu_onsets %>% select(state, season, flu_onset, flu_end)

rsv_onsets = rsv_onsets %>% select(state, season, rsv_onset, rsv_end)

# Create a function to find the best matching RSV onset for each flu onset
match_onsets = function(flu_data, rsv_data) {
  # Create a crossed dataset with all possible flu-rsv onset combinations
  # within the same state and season
  crossed = flu_data %>%
    inner_join(rsv_data, by = c("state", "season"),
               relationship = "many-to-many") %>%
    # Calculate time difference in days (negative means RSV is before flu)
    mutate(
      time_diff = as.numeric(difftime(flu_onset, rsv_onset, units = "days")),
      # Prioritize RSV onsets that come before flu onsets
      # Small positive penalty for RSV after flu
      priority_diff = ifelse(time_diff < 0,
                             abs(time_diff),
                             time_diff * 1.5)
    )

  # For each flu onset, find the closest RSV onset with priority for those before
  best_matches = crossed %>%
    group_by(state, season, flu_onset, flu_end) %>%
    slice_min(priority_diff, n = 1) %>%
    ungroup() %>%
    # Calculate weeks difference for analysis
    mutate(weeks_diff = time_diff / 7) %>%
    select(state, season, flu_onset, flu_end, rsv_onset, rsv_end, weeks_diff)

  return(best_matches)
}

# Apply the function to your datasets
onsets = match_onsets(flu_onsets, rsv_onsets)

onsets %>%
  filter(season != "22-23") %>%
  ggplot(aes(x = rsv_onset, y = flu_onset)) +
  geom_point()

write_csv(onsets, file = "data/onsets_042225.csv")

## Peaks ----

flu_peaks = read_csv(
  file = "early_warning/peaks_042425/flu_times_coef_041525.csv",
  show_col_types = FALSE) %>%
  select(-1)

flu_peaks = flu_peaks %>%
  select(name_proxy, date_peak) %>%
  rename(state = name_proxy, flu_peak = date_peak)

flu_peaks = flu_peaks %>% mutate(
  year = year(flu_peak),
  month = month(flu_peak),
  season = ifelse(
    month >= 6,
    paste0(year %% 100, "-", (year + 1) %% 100),
    # If June or later: "YY-(YY+1)"
    paste0((year - 1) %% 100, "-", year %% 100)
  )   # If before June: "(YY-1)-YY"
) %>%
  select(-c(year, month))

rsv_peaks = read_csv(
  file = "early_warning/peaks_042425/rsv_times_coef_041525.csv",
  show_col_types = FALSE) %>%
  select(-1)

rsv_peaks = rsv_peaks %>%
  select(name_proxy, date_peak) %>%
  rename(state = name_proxy, rsv_peak = date_peak)

rsv_peaks = rsv_peaks %>%
  mutate(
    year = year(rsv_peak),
    month = month(rsv_peak),
    season = ifelse(
      month >= 6,
      paste0(year %% 100, "-", (year + 1) %% 100),
      # If June or later: "YY-(YY+1)"
      paste0((year - 1) %% 100, "-", year %% 100)
    )   # If before June: "(YY-1)-YY"
  ) %>%
  select(-c(year, month))

# Create a function to find the best matching RSV peak for each flu peak
match_peaks = function(flu_data, rsv_data) {
  # Create a crossed dataset with all possible flu-rsv peak combinations
  # within the same state and season
  crossed = flu_data %>%
    inner_join(rsv_data, by = c("state", "season"),
               relationship = "many-to-many") %>%
    # Calculate time difference in days (negative means RSV peak is before flu peak)
    mutate(
      time_diff = as.numeric(difftime(flu_peak, rsv_peak, units = "days")),
      # Prioritize RSV peaks that come before flu peaks
      # Small positive penalty for RSV after flu
      priority_diff = ifelse(time_diff < 0,
                             abs(time_diff),
                             time_diff * 1.5)
    )

  # For each flu peak, find the closest RSV peak with priority for those before
  best_matches = crossed %>%
    group_by(state, season, flu_peak) %>%
    slice_min(priority_diff, n = 1) %>%
    ungroup() %>%
    # Calculate weeks difference for analysis
    mutate(weeks_diff = time_diff / 7) %>%
    select(state, season, flu_peak, rsv_peak, weeks_diff)

  return(best_matches)
}

# Apply the function to your datasets
peaks = match_peaks(flu_peaks, rsv_peaks)

# Looks reasonable - we see that some diseases have multiple peaks in seasons.

# Calculating the distribution of the differences between onsets and between
# peaks according to EWS

## Using onsets ----
onsets = onsets %>%
  mutate(weeks_diff = ifelse(is.na(rsv_onset), NA,
                             as.numeric(difftime(flu_onset,
                                                 rsv_onset, units = "weeks"))))

### As a histogram ----
onsets %>%
  drop_na() %>% # Includes only 111/150 seasons
  ggplot(aes(x = weeks_diff)) +
  geom_histogram(binwidth = 2) +
  labs(x = "Difference (weeks)", y = "Count") +
  theme_minimal()

### As a boxplot -----
onsets %>%
  drop_na() %>%
  ggplot(aes(x = weeks_diff)) +
  geom_boxplot(notch = TRUE) +
  labs(x = "Difference (weeks)", y = "Count") +
  theme_minimal()

## Using peaks
peaks_shortest = peaks %>%
  mutate(weeks_diff = ifelse(is.na(rsv_peak), NA, as.numeric(difftime(flu_peak, rsv_peak, units = "weeks")))) %>%
  group_by(state, season) %>%  # Group by state and season
  filter(!is.na(weeks_diff)) %>%  # Remove rows where weeks_diff is NA
  slice_min(abs(weeks_diff), n = 1) %>%  # Keep the row with the shortest absolute time difference
  ungroup()

### Histogram
peaks_shortest %>%
  ggplot(aes(x =  weeks_diff)) +
  geom_histogram(binwidth = 2) +
  labs(x = "Difference (weeks)", y = "Count") +
  theme_minimal()

median(onsets$weeks_diff, na.rm = TRUE)
mean(onsets$weeks_diff, na.rm = TRUE)
# 1 week diff in onsets
# Calculate the 95% percentile range
quantile(onsets$weeks_diff, probs = c(0.025, 0.975), na.rm = TRUE)

# Calculate the 95% percentile range of peak difference
median(peaks_shortest$weeks_diff, na.rm = TRUE)
mean(peaks_shortest$weeks_diff, na.rm = TRUE)
quantile(peaks_shortest$weeks_diff, probs = c(0.025, 0.975), na.rm = TRUE)

# 2 week diff in peaks

# Grouped boxplot

onsets$cat = "Onsets"
peaks_shortest$cat = "Peaks"

onsets_plot = onsets %>% select(state, season, cat, weeks_diff)
peaks_plot = peaks_shortest %>% select(state, season, cat, weeks_diff)

plot_onsets_peaks = bind_rows(onsets_plot, peaks_plot)

plot_onsets_peaks %>%
  ggplot(aes(x = weeks_diff, y = cat))+
  geom_boxplot(notch = TRUE) +
  labs(x= "Difference (weeks)", y = "") +
  envalysis::theme_publish()
ggsave(filename = "figures/fig2_v1.png", width = 10,
       height = 8, units = "in", bg = "white")
