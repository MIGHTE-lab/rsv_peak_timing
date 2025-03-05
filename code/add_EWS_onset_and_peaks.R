## -----------------------------------------------------------------------------
## Script name: add_EWS_onset_and_peaks.R
##
## Purpose of script: Add EWS onset and peak dates to RSV/flu timing plots
##
## Author: George Dewey
##
## Date Created: 2025-03-03
##
## Last Updated: 2025-03-03
## -----------------------------------------------------------------------------

## Load packages
library(tidyverse)
library(lubridate)
library(ggh4x)
library(patchwork)

## Load NSSP data
nssp = read_csv(file = "data/nssp_all_years.csv", show_col_types = FALSE)

## Test to see if we can regenerate the indivdual season plots from the individual years
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
flu_onsets = read_csv(file = "early_warning/windows/flu_times_coef.csv", show_col_types = FALSE) %>%
  select(-1)

flu_onsets = flu_onsets %>%
  select(name_proxy, start_date, end_date) %>%
  rename(state = name_proxy, flu_onset = start_date, flu_end = end_date)

rsv_onsets = read_csv(file = "early_warning/windows/rsv_times_coef.csv", show_col_types = FALSE) %>%
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

flu_onsets = flu_onsets %>% slice(-c(8, 133))

rsv_onsets = rsv_onsets %>% select(state, season, rsv_onset, rsv_end)

rsv_onsets = rsv_onsets %>% slice(-98)

onsets = flu_onsets %>% left_join(rsv_onsets, by = c("state", "season"))

onsets = onsets %>% select(state, season, flu_onset, flu_end, rsv_onset, rsv_end)

write_csv(onsets, file = "data/onsets.csv")

## Peaks ----

flu_peaks = read_csv(file = "early_warning/peaks/flu_times_coef.csv", show_col_types = FALSE) %>%
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

flu_peaks

rsv_peaks = read_csv(file = "early_warning/peaks/rsv_times_coef.csv", show_col_types = FALSE) %>%
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

peaks = flu_peaks %>% left_join(rsv_peaks, by = c("state", "season"), relationship = "many-to-many")

# Looks reasonable - we see that some diseases have multiple peaks in seasons. I guess the algorithm could not detect peaks in 22-23 because we started data collection too late in the season

# Calculating the distribution of the differences between onsets and between peaks according to EWS

## Using onsets ----
onsets = onsets %>% mutate(weeks_diff = ifelse(is.na(rsv_onset), NA, as.numeric(difftime(flu_onset, rsv_onset, units = "weeks"))))

### As a histogram ----
onsets %>%
  drop_na() %>% # Includes only 106/150 seasons
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

onsets$cat = "Onsets"
peaks_shortest$cat = "Peaks"

onsets_plot = onsets %>% select(state, season, cat, weeks_diff)
peaks_plot = peaks_shortest %>% select(state, season, cat, weeks_diff)

plot_onsets_peaks = bind_rows(onsets_plot, peaks_plot)

plot_onsets_peaks %>%
  ggplot(aes(x = weeks_diff, y = cat))+
  geom_boxplot(notch = TRUE) +
  labs(x= "Difference (weeks)", y = "Data") +
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

median(peaks_shortest$weeks_diff, na.rm = TRUE)
mean(peaks_shortest$weeks_diff, na.rm = TRUE)
