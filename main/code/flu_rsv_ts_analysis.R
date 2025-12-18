## -----------------------------------------------------------------------------
## Script name: flu_rsv_ts_analysis.R
##
## Purpose of script: Conduct formal time series analysis on seasonal ILI data
##
## Author: George Dewey
##
## Date Created: 2025-01-29
##
## Last Updated: 2025-04-14
##
## Version: 1.3.1
##
## Update notes: Clean up code for manuscript submission
## -----------------------------------------------------------------------------

## 1. Load packages and data ---------------------------------------------------
## Data management
library(tidyverse)
library(zoo)

## Working with dates
library(lubridate)
library(MMWRweek)

## Visualization
library(ggh4x)

## Modeling
library(glmnet)
library(penalized)
library(caret)

## Set working directory
setwd("~/Documents/Projects/rsv_peak_timing")

## Load NSSP data
dat = read_csv(
  "data/NSSP/nssp_full_040525.csv"
)

## Load ILI data and do some data management
## ILI data - need data from 21 - 23
ili_data = read_csv("data/ILInet/ILINet_state_2025_03-29.csv",
                     na = "X",
                     skip = 1)

names(ili_data) = c(
  "region_type",
  "region",
  "year",
  "epiweek",
  "weighted_ili",
  "unweighted_ili",
  "age_0_4",
  "age_25_49",
  "age_25_64",
  "age_5_24",
  "age_50_64",
  "age_65",
  "ilitotal",
  "num_data_providers",
  "total_patients"
)

ili_data_with_dates = ili_data %>%
  filter(year %in% 2022:2025) %>%
  mutate(week_end = as_date(MMWRweek2Date(year, epiweek) + 6)) %>%
  select(region, year, epiweek, week_end, unweighted_ili, ilitotal) %>%
  rename(state = region)

dat = dat %>%
  filter(county == "All") %>%
  select(
    week_end,
    geography,
    percent_visits_covid,
    percent_visits_influenza,
    percent_visits_rsv
  )

# Create season variable
dat = dat %>%
  mutate(
    year = year(week_end),
    month = month(week_end),
    ili_season = ifelse(
      month >= 6,
      paste0(year %% 100, "-", (year + 1) %% 100),
      # If June or later: "YY-(YY+1)"
      paste0((year - 1) %% 100, "-", year %% 100)
    )   # If before June: "(YY-1)-YY"
  )

dat = dat %>% rename(
  state = geography,
  covid = percent_visits_covid,
  flu = percent_visits_influenza,
  rsv = percent_visits_rsv,
  season = ili_season
) %>%
  select(week_end, state, covid, flu, rsv, season)

dat_combined = dat %>%
  left_join(ili_data_with_dates, by = c("state", "week_end")) %>%
  rename(ili = unweighted_ili) %>%
  select(-ilitotal)

minmax = function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

dat_combined_ = dat_combined %>% drop_na()

## 2. Create modeling functions ----
# Build function which takes in a dataset and returns the optimal ridge
# coefficients after using gridsearch for the optimal lambda value

ridge_CV = function(data) {

    # Run the ridge regression using glmnet instead of penalized
    # glmnet is chosen for its efficient implementation and support for cross-validation.
    # Define predictors and outcome

    predictors = c("covid_rescaled", "flu_rescaled", "rsv_rescaled")
    X = as.matrix(data[, predictors])
    y = data$ili_rescaled

    lambda_values = 10^seq(-5, 5, length = 100)

    cv_ridge = cv.glmnet(X, y, alpha = 0, lambda = lambda_values, nfolds = 10,
                         intercept = FALSE, lower.limits = 0)

    # Implement gridsearch for lambda values
    lambda_grid = seq(cv_ridge$lambda.min / 1000, cv_ridge$lambda.min * 1000, length = 500)

    # Fit models for each lambda
    ridge_models = lapply(lambda_grid, function(l) {
      glmnet(X, y, alpha = 0, lambda = l, intercept = FALSE, lower.limits = 0)
    })

    # Extract coefficients for each lambda
    ridge_coefs = sapply(ridge_models, function(model) coef(model)[,1])

    # Identify the first lambda where all coefficients are nonzero
    valid_lambda_idx = which(colSums(ridge_coefs != 0) == length(predictors))[1]
    if (is.null(valid_lambda_idx)) {
      stop("No valid lambda found where all coefficients are nonzero.")
    }
    best_lambda = lambda_grid[valid_lambda_idx]

    # Fit the final ridge regression with this lambda
    ridge_final = glmnet(X, y, alpha = 0, lambda = best_lambda, intercept = FALSE, lower.limits = 0)

    # Extract final coefficients
    coef_values = as.numeric(coef(ridge_final))

    return(coef_values) # 2: covid, 3: flu, 4: rsv
}

# 3. Run models and create state-season plots ----

# 22-23 season
nssp_signals_22_23 = NULL
for (geo in dat_combined_ %>% distinct(state) %>% pull()) {

  # Filter the data to only include data from the selected state and season
  data = dat_combined_ %>% filter(state == geo, season == "22-23") %>%
    mutate(
      ili_rescaled = minmax(ili),
      flu_rescaled = minmax(flu),
      rsv_rescaled = minmax(rsv),
      covid_rescaled = minmax(covid))

  # Use ridge.cv to generate the coefficients
  coefs_tmp = ridge_CV(data)

  signals_tmp = data %>%
    mutate(
      beta_flu = coefs_tmp[3],
      beta_rsv = coefs_tmp[4],
      beta_covid = coefs_tmp[2],
      flu_times_coef = flu_rescaled * beta_flu,
      rsv_times_coef = rsv_rescaled * beta_rsv,
      covid_times_coef = covid_rescaled * beta_covid,
    )

  nssp_signals_22_23 = bind_rows(nssp_signals_22_23, signals_tmp)
  print(paste0("State ",geo," completed."))
}
nssp_signals_22_23 %>%
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
  ),
  breaks = c('ILI', 'Flu', 'RSV', 'COVID-19')  # This sets the legend order
  ) +
  labs(x = '', y = '') +
  scale_x_date(date_labels = "%b\n%y", date_breaks = "2 months") +
  envalysis::theme_publish() +
  theme(
    axis.text.x = element_text(size = 5),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    legend.title = element_blank(),
    legend.position = 'right',
    strip.text = element_text(size = 9)
  )
# ggsave(file = "figures/nssp_22_23_060525.png", width = 16, height = 10, units = "in", bg = "white")

# 23-24 season
nssp_signals_23_24 = NULL
for (geo in dat_combined_ %>% distinct(state) %>% pull()) {

  # Filter the data to only include data from the selected state and season
  data = dat_combined_ %>% filter(state == geo, season == "23-24") %>%
    mutate(
      ili_rescaled = minmax(ili),
      flu_rescaled = minmax(flu),
      rsv_rescaled = minmax(rsv),
      covid_rescaled = minmax(covid))

  # Use ridge.cv to generate the coefficients
  coefs_tmp = ridge_CV(data)

  signals_tmp = data %>%
    mutate(
      beta_flu = coefs_tmp[3],
      beta_rsv = coefs_tmp[4],
      beta_covid = coefs_tmp[2],
      flu_times_coef = flu_rescaled * beta_flu,
      rsv_times_coef = rsv_rescaled * beta_rsv,
      covid_times_coef = covid_rescaled * beta_covid,
    )

  nssp_signals_23_24 = bind_rows(nssp_signals_23_24, signals_tmp)
  print(paste0("State ",geo," completed."))
}
nssp_signals_23_24 %>%
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
  ),
  breaks = c('ILI', 'Flu', 'RSV', 'COVID-19')  # This sets the legend order
  ) +
  labs(x = '', y = '') +
  scale_x_date(date_labels = "%b\n%y", date_breaks = "2 months") +
  envalysis::theme_publish() +
  theme(
    axis.text.x = element_text(size = 5),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    legend.title = element_blank(),
    legend.position = 'right',
    strip.text = element_text(size = 9)
  )
# ggsave(file = "figures/nssp_23_24_060525.png", width = 16, height = 10, units = "in", bg = "white")

# 24-25 season
nssp_signals_24_25 = NULL
for (geo in dat_combined_ %>% filter(state != "Wyoming") %>% distinct(state) %>% pull()) {

  # Filter the data to only include data from the selected state and season
  data = dat_combined_ %>% filter(state == geo, season == "24-25") %>%
    mutate(
      ili_rescaled = minmax(ili),
      flu_rescaled = minmax(flu),
      rsv_rescaled = minmax(rsv),
      covid_rescaled = minmax(covid))

  # Use ridge.cv to generate the coefficients
  coefs_tmp = ridge_CV(data)

  signals_tmp = data %>%
    mutate(
      beta_flu = coefs_tmp[3],
      beta_rsv = coefs_tmp[4],
      beta_covid = coefs_tmp[2],
      flu_times_coef = flu_rescaled * beta_flu,
      rsv_times_coef = rsv_rescaled * beta_rsv,
      covid_times_coef = covid_rescaled * beta_covid,
    )

  nssp_signals_24_25 = bind_rows(nssp_signals_24_25, signals_tmp)
  print(paste0("State ",geo," completed."))
}
nssp_signals_24_25 %>%
  filter(state != "Vermont") %>%
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
  ),
  breaks = c('ILI', 'Flu', 'RSV', 'COVID-19')  # This sets the legend order
  ) +
  labs(x = '', y = '') +
  scale_x_date(date_labels = "%b\n%y", date_breaks = "2 months") +
  envalysis::theme_publish() +
  theme(
    axis.text.x = element_text(size = 5),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    legend.title = element_blank(),
    legend.position = 'right',
    strip.text = element_text(size = 9)
  )
# ggsave(file = "figures/nssp_24_25_060525.png", width = 16, height = 10, units = "in", bg = "white")

# 4. Generate Figure 1 (3x3 state-season example figure) ----

# First combine the three yearly datasets
nssp_all_years = bind_rows(nssp_signals_22_23, nssp_signals_23_24, nssp_signals_24_25)

# Then select only the states of interest (MD, NY, TX)
nssp_all_years %>%
  filter(state %in% c("Maryland", "Texas", "New York")) %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = flu_times_coef, color = 'Flu'), linewidth = 0.85) +
  geom_line(aes(y = rsv_times_coef, color = 'RSV'), linewidth = 0.85) +
  geom_line(aes(y = covid_times_coef, color = "COVID-19"), linewidth = 0.85) +
  geom_line(aes(y = ili_rescaled, color = 'ILI'), linewidth = 0.85) +
  facet_grid(state ~ season, scales = 'free') +
  scale_color_manual(values = c(
    'ILI' = '#E2DBC9',
    'RSV' = '#183A5A',
    'Flu' = '#C34129',
    'COVID-19' = "#EFB75B"
  )) +
  labs(x = '', y = 'Rescaled Volume / Contribution') +
  scale_x_date(date_labels = "%b\n%y", date_breaks = "2 months") +
  envalysis::theme_publish() +
  theme(
    axis.text.x = element_text(size = 5),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    legend.title = element_blank(),
    legend.position = 'bottom'
  )
ggsave(file = "figures/3x3_MD_NY_TX_100325.png", height = 8, width = 8, units = "in", bg = "white")

# 5. Determining temporal order from ridge model results ----
## Create helper functions
find_states_with_early_rsv_peak = function(data) {
  return(
   data %>%
    group_by(state) %>%
    summarize(
      max_rsv_date = week_end[which.max(rsv_rescaled)],
      max_flu_date = week_end[which.max(flu_rescaled)],
      .groups = "drop"
    ) %>%
    filter(max_rsv_date < max_flu_date) %>%
    select(state, max_rsv_date, max_flu_date) %>%
     nrow()
  )
}

find_states_with_concurrent_rsv_peak = function(data) {
  return(
    data %>%
      group_by(state) %>%
      summarize(
        max_rsv_date = week_end[which.max(rsv_rescaled)],
        max_flu_date = week_end[which.max(flu_rescaled)],
        .groups = "drop"
      ) %>%
      filter(max_rsv_date == max_flu_date) %>%
      select(state, max_rsv_date, max_flu_date) %>%
      nrow()
  )
}

find_states_with_early_flu_peak = function(data) {
  return(
    data %>%
      group_by(state) %>%
      summarize(
        max_rsv_date = week_end[which.max(rsv_rescaled)],
        max_flu_date = week_end[which.max(flu_rescaled)],
        .groups = "drop"
      ) %>%
      filter(max_rsv_date > max_flu_date) %>%
      select(state, max_rsv_date, max_flu_date) %>%
      nrow()
  )
}

# Flu before RSV
find_states_with_early_flu_peak(nssp_signals_22_23) #2
find_states_with_early_flu_peak(nssp_signals_23_24) #9
find_states_with_early_flu_peak(nssp_signals_24_25) #10
# 21 - 14.2%

find_states_with_early_rsv_peak(nssp_signals_22_23) #45
find_states_with_early_rsv_peak(nssp_signals_23_24) #35
find_states_with_early_rsv_peak(nssp_signals_24_25) #34
# 114 - 77.0%

find_states_with_concurrent_rsv_peak(nssp_signals_22_23) #3
find_states_with_concurrent_rsv_peak(nssp_signals_23_24) #6
find_states_with_concurrent_rsv_peak(nssp_signals_24_25) #5
# 14 - 9.5%

## Checking start and end dates for peaks
nssp_signals_22_23 %>%
  group_by(state) %>%
  summarize(
    max_rsv_date = week_end[which.max(rsv_rescaled)],
    max_flu_date = week_end[which.max(flu_rescaled)],
    .groups = "drop"
  ) %>%
  filter(max_rsv_date == max_flu_date) %>%
  select(state, max_rsv_date, max_flu_date)

## Creating a table showing the peak timings for each season
nssp_peaks_22_23 = nssp_signals_22_23 %>%
  group_by(state) %>%
  summarize(
    max_rsv_date = week_end[which.max(rsv_rescaled)],
    max_flu_date = week_end[which.max(flu_rescaled)],
    max_cov_date = week_end[which.max(covid_rescaled)],
    max_ili_date = week_end[which.max(ili_rescaled)],
    .groups = "drop"
  ) %>%
  mutate(season = "22-23")

nssp_peaks_23_24 = nssp_signals_23_24 %>%
  group_by(state) %>%
  summarize(
    max_rsv_date = week_end[which.max(rsv_rescaled)],
    max_flu_date = week_end[which.max(flu_rescaled)],
    max_cov_date = week_end[which.max(covid_rescaled)],
    max_ili_date = week_end[which.max(ili_rescaled)],
    .groups = "drop"
  ) %>%
  mutate(season = "23-24")

nssp_peaks_24_25 = nssp_signals_24_25 %>%
  group_by(state) %>%
  summarize(
    max_rsv_date = week_end[which.max(rsv_rescaled)],
    max_flu_date = week_end[which.max(flu_rescaled)],
    max_cov_date = week_end[which.max(covid_rescaled)],
    max_ili_date = week_end[which.max(ili_rescaled)],
    .groups = "drop"
  ) %>%
  mutate(season = "24-25")

nssp_peaks = bind_rows(nssp_peaks_22_23, nssp_peaks_23_24, nssp_peaks_24_25)

write_csv(nssp_peaks, file = "data/nssp_peaks.csv")


# Getting the month ranges for the peak values of RSV and influenza
nssp_signals_22_23 %>%
  group_by(state) %>%
  summarize(
    max_rsv_date = week_end[which.max(rsv_rescaled)],
    max_flu_date = week_end[which.max(flu_rescaled)],
    max_ili_date = week_end[which.max(ili_rescaled)],
    .groups = "drop"
  ) %>%
 summarize(min_rsv = min(max_rsv_date),
           max_rsv = max(max_rsv_date),
           min_flu = min(max_flu_date),
           max_flu = max(max_flu_date),
           min_ili = min(max_ili_date),
           max_ili = max(max_ili_date))
# rsv: min Oct 1, max_rsv Dec 24
# flu: min Nov 5 max Dec 31
# ili: min Nov 5 max Dec 31

nssp_signals_23_24  %>%
  group_by(state) %>%
  summarize(
    max_rsv_date = week_end[which.max(rsv_rescaled)],
    max_flu_date = week_end[which.max(flu_rescaled)],
    max_ili_date = week_end[which.max(ili_rescaled)],
    .groups = "drop"
  ) %>%
  summarize(min_rsv = min(max_rsv_date),
            max_rsv = max(max_rsv_date),
            min_flu = min(max_flu_date),
            max_flu = max(max_flu_date),
            min_ili = min(max_ili_date),
            max_ili = max(max_ili_date))
# rsv: min Nov 4, max Feb 24
# flu: min Nov 4, max March 9
# ili: min Nov 4, Max March 9

nssp_signals_24_25  %>%
  group_by(state) %>%
  summarize(
    max_rsv_date = week_end[which.max(rsv_rescaled)],
    max_flu_date = week_end[which.max(flu_rescaled)],
    max_ili_date = week_end[which.max(ili_rescaled)],
    .groups = "drop"
  ) %>%
  summarize(min_rsv = min(max_rsv_date),
            max_rsv = max(max_rsv_date),
            min_flu = min(max_flu_date),
            max_flu = max(max_flu_date),
            min_ili = min(max_ili_date),
            max_ili = max(max_ili_date))
# RSV: min June 1, max Feb 15
# Flu: min June 15, max Feb 15
# ILI: min Sep 28, max Feb 15

# Exporting the data for early warning system
nssp_all_years = nssp_all_years %>% select(week_end, state, season, covid, flu, rsv, ili,
                          contains("rescaled"), contains("beta"), contains("times"))

write_csv(nssp_all_years, file = "data/nssp_all_years_0605.csv")
