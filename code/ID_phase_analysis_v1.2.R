## -----------------------------------------------------------------------------
## Script name: ID_phase_analysis.R
##
## Purpose of script: Conduct formal time series analysis on seasonal ILI data
##
## Author: George Dewey
##
## Date Created: 2025-01-29
##
## Last Updated: 2025-03-04
##
## Version: 1.1
## Add analysis including Feb 2025 data
##
## -----------------------------------------------------------------------------

## 1. Load packages and data ---------------------------------------------------
## Data management
library(tidyverse)
library(zoo)

## Working with dates
library(lubridate)
library(MMWRweek)

## Visualization
library(firatheme)
library(ggh4x)

## Modeling
library(glmnet)
library(penalized)
library(caret)

## Set working directory
setwd("~/Documents/Projects/rsv_peak_timing")

## Load NSSP data
dat = read_csv(
  "data/NSSP/nssp_full_02282025.csv"
)

## Load ILI data and do some data management
## ILI data - need data from 21 - 23
ili_data = read_csv("data/ILInet/ILINet_state_2025-02-22.csv",
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


dat = dat %>% filter(county == "All") %>% select(
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

dat_combined = dat %>% left_join(ili_data_with_dates, by = c("state", "week_end")) %>%
  rename(ili = unweighted_ili) %>%
  select(-ilitotal)

minmax = function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

dat_combined_ = dat_combined %>% drop_na()

# Let's do some data checks here
dat_MT_2425 = dat_combined_ %>% filter(state == "Montana", season == "24-25")  %>%
  mutate(
    ili_rescaled = minmax(ili),
    flu_rescaled = minmax(flu),
    rsv_rescaled = minmax(rsv),
    covid_rescaled = minmax(covid)
  )

# Use 10-fold CV to identify best lambda
  cv = cv.glmnet(
    x = data.matrix(dat_MT_2425[, c('flu_rescaled', 'rsv_rescaled', 'covid_rescaled')]),
    y = data.matrix(dat_MT_2425[, 'ili_rescaled']),
    intercept = FALSE,
    alpha = 0
  )

  cv$lambda.min
  cv$lambda.1se

# Fit the penalized model using the specified lambda from the CV model
  fit = penalized(
    ili_rescaled ~ flu_rescaled + rsv_rescaled + covid_rescaled ,
    positive = TRUE,
    unpenalized = ~ 0,
    data = dat_MT_2425,
    lambda2 = cv$lambda.1se
  )

coefficients(fit)

# Let's do some data checks here
dat_MI_2324 = dat_combined_ %>% filter(state == "Michigan", season == "23-24")  %>%
  mutate(
    ili_rescaled = minmax(ili),
    flu_rescaled = minmax(flu),
    rsv_rescaled = minmax(rsv),
    covid_rescaled = minmax(covid)
  )

# Use 10-fold CV to identify best lambda
cv = cv.glmnet(
  x = data.matrix(dat_MI_2324[, c('flu_rescaled', 'rsv_rescaled', 'covid_rescaled')]),
  y = data.matrix(dat_MI_2324[, 'ili_rescaled']),
  intercept = FALSE,
  alpha = 0
)

cv$lambda.min
cv$lambda.1se

# Fit the penalized model using the specified lambda from the CV model
fit = penalized(
  ili_rescaled ~ flu_rescaled + rsv_rescaled + covid_rescaled ,
  positive = TRUE,
  unpenalized = ~ 0,
  data = dat_MI_2324,
  lambda2 = cv$lambda.1se
)

library(hdm)
ridge_adaptive <- rlasso(
  x = as.matrix(dat_MT_2425[, c('flu_rescaled', 'rsv_rescaled', 'covid_rescaled')]),
  y = dat_MT_2425$ili_rescaled,
  post = FALSE # Standard adaptive ridge
)


dat_MI_2324

coefficients(fit)

cor(dat_MI_2324$flu_rescaled, dat_MI_2324$rsv_rescaled)


# 2. Running the penalized ridge models ----

# 22-23 season ----
nssp_signals_22_23 = NULL
for (geo in dat_combined_ %>% distinct(state) %>% pull()) {
  # Filter the data to only include data from the selected state and season
  data = dat_combined_ %>% filter(state == geo, season == "22-23") %>%
    mutate(
      ili_rescaled = minmax(ili),
      flu_rescaled = minmax(flu),
      rsv_rescaled = minmax(rsv),
      covid_rescaled = minmax(covid)
    )

  # Use 10-fold CV to identify best lambda
  cv = cv.glmnet(
    x = data.matrix(data[, c('flu_rescaled', 'rsv_rescaled', 'covid_rescaled')]),
    y = data.matrix(data[, 'ili_rescaled']),
    intercept = FALSE,
    alpha = 0
  )

  # Use penalized package to run the ridge w/ no intercept + non-negative
  # and the specified lambda from the CV model
  fit = penalized(
    ili_rescaled ~ flu_rescaled + rsv_rescaled + covid_rescaled ,
    positive = TRUE,
    unpenalized = ~ 0,
    data = data,
    lambda2 = cv$lambda.1se
  )

  coefs_tmp = coefficients(fit)

  signals_tmp = data %>%
    mutate(
      beta_flu = coefs_tmp[1],
      beta_rsv = coefs_tmp[2],
      beta_covid = coefs_tmp[3],
      flu_times_coef = flu_rescaled * beta_flu,
      rsv_times_coef = rsv_rescaled * beta_rsv,
      covid_times_coef = covid_rescaled * beta_covid,
    )

  nssp_signals_22_23 = bind_rows(nssp_signals_22_23, signals_tmp)

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

# ggsave(file = "figures/nssp_22_23.png", width = 16, height = 10, units = "in", bg = "white")

# 23-24 season ----
nssp_signals_23_24 = NULL
for (geo in dat_combined_ %>% distinct(state) %>% pull()) {
  # Filter the data to only include data from the selected state and season
  data = dat_combined_ %>% filter(state == geo, season == "23-24") %>%
    mutate(
      ili_rescaled = minmax(ili),
      flu_rescaled = minmax(flu),
      rsv_rescaled = minmax(rsv),
      covid_rescaled = minmax(covid)
    )

  # Use 10-fold CV to identify best lambda
  cv = cv.glmnet(
    x = data.matrix(data[, c('flu_rescaled', 'rsv_rescaled', 'covid_rescaled')]),
    y = data.matrix(data[, 'ili_rescaled']),
    intercept = FALSE,
    alpha = 0
  )

  # Use penalized package to run the ridge w/ no intercept + non-negative
  # and the specified lambda from the CV model
  fit = penalized(
    ili_rescaled ~ flu_rescaled + rsv_rescaled + covid_rescaled ,
    positive = TRUE,
    unpenalized = ~ 0,
    data = data,
    lambda2 = cv$lambda.1se
  )

  coefs_tmp = coefficients(fit)

  signals_tmp = data %>%
    mutate(
      beta_flu = coefs_tmp[1],
      beta_rsv = coefs_tmp[2],
      beta_covid = coefs_tmp[3],
      flu_times_coef = flu_rescaled * beta_flu,
      rsv_times_coef = rsv_rescaled * beta_rsv,
      covid_times_coef = covid_rescaled * beta_covid,
    )

  nssp_signals_23_24 = bind_rows(nssp_signals_23_24, signals_tmp)

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
# ggsave(file = "figures/nssp_23_24.png", width = 16, height = 10, units = "in", bg = "white")

# 24-25 season ----
nssp_signals_24_25 = NULL
for (geo in dat_combined_ %>% distinct(state) %>% pull()) {
  # Filter the data to only include data from the selected state and season
  data = dat_combined_ %>% filter(state == geo, season == "24-25") %>%
    mutate(
      ili_rescaled = minmax(ili),
      flu_rescaled = minmax(flu),
      rsv_rescaled = minmax(rsv),
      covid_rescaled = minmax(covid)
    )

  # Use 10-fold CV to identify best lambda
  cv = cv.glmnet(
    x = data.matrix(data[, c('flu_rescaled', 'rsv_rescaled', 'covid_rescaled')]),
    y = data.matrix(data[, 'ili_rescaled']),
    intercept = FALSE,
    alpha = 0
  )

  # Use penalized package to run the ridge w/ no intercept + non-negative
  # and the specified lambda from the CV model
  fit = penalized(
    ili_rescaled ~ flu_rescaled + rsv_rescaled + covid_rescaled ,
    positive = TRUE,
    unpenalized = ~ 0,
    data = data,
    lambda2 = cv$lambda.min
  )

  coefs_tmp = coefficients(fit)

  signals_tmp = data %>%
    mutate(
      beta_flu = coefs_tmp[1],
      beta_rsv = coefs_tmp[2],
      beta_covid = coefs_tmp[3],
      flu_times_coef = flu_rescaled * beta_flu,
      rsv_times_coef = rsv_rescaled * beta_rsv,
      covid_times_coef = covid_rescaled * beta_covid,
    )

  nssp_signals_24_25 = bind_rows(nssp_signals_24_25, signals_tmp)

}
nssp_signals_24_25 %>%
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
ggsave(file = "figures/nssp_24_25.png", width = 16, height = 10, units = "in", bg = "white")

# Generating the 3x3 panel Figure 1. Show the seasonal changes for MD, NY, and TX

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
  labs(x = '', y = '') +
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
ggsave(file = "figures/3x3 _MD_NY_TX_v2.png", height = 8, width = 8, units = "in", bg = "white")

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

# Function calls
# Flu before RSV
find_states_with_early_flu_peak(nssp_signals_22_23) #2
find_states_with_early_flu_peak(nssp_signals_23_24) #9
find_states_with_early_flu_peak(nssp_signals_24_25) #8
# 19

find_states_with_early_rsv_peak(nssp_signals_22_23) #45
find_states_with_early_rsv_peak(nssp_signals_23_24) #35
find_states_with_early_rsv_peak(nssp_signals_24_25) #36
# 126

126/150

find_states_with_concurrent_rsv_peak(nssp_signals_22_23) #3
find_states_with_concurrent_rsv_peak(nssp_signals_23_24) #6
find_states_with_concurrent_rsv_peak(nssp_signals_24_25) #6
# 15

all_states <- c(
  "Alabama", "Alaska", "Arizona", "Arkansas", "California", "Colorado",
  "Connecticut", "Delaware", "District of Columbia", "Florida", "Georgia",
  "Hawaii", "Idaho", "Illinois", "Indiana", "Iowa", "Kansas", "Kentucky",
  "Louisiana", "Maine", "Maryland", "Massachusetts", "Michigan", "Minnesota",
  "Mississippi", "Missouri", "Montana", "Nebraska", "Nevada", "New Hampshire",
  "New Jersey", "New Mexico", "New York", "North Carolina", "North Dakota",
  "Ohio", "Oklahoma", "Oregon", "Pennsylvania", "Rhode Island", "South Carolina",
  "South Dakota", "Tennessee", "Texas", "Utah", "Vermont", "Virginia",
  "Washington", "West Virginia", "Wisconsin", "Wyoming"
)

setdiff(all_states, unique(dat_combined_$state))

nssp_22_23_flu_before_rsv = find_states_with_early_flu_peak(nssp_signals_22_23)
nssp_22_23_rsv_before_flu = find_states_with_early_rsv_peak(nssp_signals_22_23)
nssp_22_23_concurrent_peak = find_states_with_concurrent_rsv_peak(nssp_signals_22_23)

nssp_22_23_check = bind_rows(nssp_22_23_flu_before_rsv, nssp_22_23_rsv_before_flu, nssp_22_23_concurrent_peak)
setdiff(all_states, nssp_22_23_check$state)

nssp_signals_22_23

nssp_signals_22_23 %>%
  group_by(state) %>%
  summarize(
    max_rsv_date = week_end[which.max(rsv_rescaled)],
    max_flu_date = week_end[which.max(flu_rescaled)],
    .groups = "drop"
  ) %>%
  filter(max_rsv_date == max_flu_date) %>%
  select(state, max_rsv_date, max_flu_date) %>%
  nrow()

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

# In the current season, this method does not really work. At least for now,
# the data suggests that flu and RSV are currently peaking


# Evaluating temporal relationships between RSV and Flu signals
nssp_signals_22_23 = nssp_signals_22_23 %>% select(-c(year, epiweek))

# write_csv(nssp_signals_22_23, file = "data/nssp_signals_22_23_processed.csv")

# Exporting the data for early warning system
nssp_all_years = nssp_all_years %>% select(week_end, state, season, covid, flu, rsv, ili,
                          contains("rescaled"), contains("beta"), contains("times"))

write_csv(nssp_all_years, file = "data/nssp_all_years.csv")
# Phase Analysis ------

# First bind all the data together and create lags of the signals
nssp_full = bind_rows(nssp_signals_22_23, nssp_signals_23_24, nssp_signals_24_25)

nssp_full = nssp_full %>%
  mutate(
    vir_rhs = covid + flu + rsv,
    vir_rescaled_rhs = covid_rescaled + flu_rescaled + rsv_rescaled
  )



# Find optimal lag for ILI and the regression rhs
ccf(nssp_full$vir_rhs,
    nssp_full$ili_rescaled,
    lag.max = 8,
    plot = TRUE)
# In fact the optimal lag is @ lag 0

nssp_full_lagged = nssp_full %>%
  group_by(state) %>%
  mutate(
    # Create lags 0-4 for each virus
    RSV_0 = rsv,
    RSV_1 = lag(rsv, 1),
    RSV_2 = lag(rsv, 2),
    RSV_3 = lag(rsv, 3),
    RSV_4 = lag(rsv, 4),
    RSV_5 = lag(rsv, 5),
    RSV_6 = lag(rsv, 6),

    FLU_0 = flu,
    FLU_1 = lag(flu, 1),
    FLU_2 = lag(flu, 2),
    FLU_3 = lag(flu, 3),
    FLU_4 = lag(flu, 4),
    FLU_5 = lag(flu, 5),
    FLU_6 = lag(flu, 6),

    # Create seasonal terms
    week_num = week(week_end),
    sin_annual = sin(2 * pi * week_num / 52),
    cos_annual = cos(2 * pi * week_num / 52),
    sin_semiannual = sin(4 * pi * week_num / 52),
    cos_semiannual = cos(4 * pi * week_num / 52)
  )

nssp_full_lagged = nssp_full_lagged %>% select(-c(year, epiweek))
nssp_full_lagged[is.na(nssp_full_lagged) == TRUE] = 0

# Fit the model
model <- nlme::lme(
  ili ~ RSV_0 + RSV_1 + RSV_2 + RSV_3 + RSV_4 +
    FLU_0 + FLU_1 + FLU_2 + FLU_3 + FLU_4 +
    sin_annual + cos_annual +
    sin_semiannual + cos_semiannual,
  random = ~ 1 | state / season,
  correlation = nlme::corAR1(),
  data = nssp_full_lagged
)

summary(model)

# If we want the model to identify whether flu leads RSV or vice-versa, we
# create two models: a lagged mixed-effects model where RSV is the outcome and
# a lagged mixed effects model where Flu is the outcome.

nssp_full_lagged[is.na(nssp_full_lagged) == TRUE] = 0

# Modeling with nlme ----

model_rsv_flu = nlme::lme(
  rsv ~ flu + FLU_1 + FLU_2 + FLU_3 + FLU_4 +
    FLU_5 + FLU_6 + sin_annual + cos_annual +
    sin_semiannual + cos_semiannual,
  random = ~ 1 | state / season,
  correlation = nlme::corAR1(),
  data = nssp_full_lagged
)

model_flu_rsv = nlme::lme(
  flu ~ rsv + RSV_1 + RSV_2 + RSV_3 + RSV_4 +
    RSV_5 + RSV_6 +  sin_annual + cos_annual +
    sin_semiannual + cos_semiannual,
  random = ~ 1 | state / season,
  correlation = nlme::corAR1(),
  data = nssp_full_lagged
)

# Alternative testing using lmer ----
library(lme4)
library(lmerTest)

rsv_model = lmer(
  rsv ~ flu + FLU_1 + FLU_2 + FLU_3 + FLU_4 + FLU_5 + FLU_6 +
    sin_annual + cos_annual + sin_semiannual + cos_semiannual +
    (1 | state) + (1 | season:state),
  data = nssp_full_lagged, REML = TRUE
)
summary(rsv_model)

flu_model = lmer(
  flu ~ rsv + RSV_1 + RSV_2 + RSV_3 + RSV_4 + RSV_5 + RSV_6 +
    sin_annual + cos_annual + sin_semiannual + cos_semiannual +
    (1| state) + (1| season:state),
  data = nssp_full_lagged, REML = TRUE
)
 summary(flu_model)

# If betas of the lagged flu predictors are positive and significant in the model that predicts rsv using flu, this suggests Flu precedes RSV
summary(model_rsv_flu)
ccf(nssp_full_lagged$rsv,
    nssp_full_lagged$flu,
    lag.max = 8,
    plot = TRUE)

# If betas of the lagged rsv predictors are positive and signfiicant in the model that predicts flu using rsv, this suggests RSV precedes flu
summary(model_flu_rsv)

# Based on these coefficients RSV precedes flu at 2 and 4-week intervals

# Cross correlation heatmap

lags <- 0:6
cor_matrix <- sapply(lags, function(l)
  cor(nssp_full_lagged$flu, lag(nssp_full_lagged$rsv, l), use = "complete.obs"))
cor_df <- data.frame(Lag = lags, Correlation = cor_matrix)

# Create a big lagged correlation table for each state per season
nssp_full_lagged %>% select(
  week_end,
  state,
  flu_rescaled,
  ili_rescaled,
  contains("RSV", ignore.case = FALSE),
  contains("FLU", ignore.case = FALSE)
) %>%
  filter(state == "Arizona")

nssp_corrs_flu_lag = nssp_full_lagged %>%
  group_by(state) %>%
  summarise(across(
    starts_with("FLU_", ignore.case = FALSE),
    ~ cor(.x, rsv_rescaled, use = "pairwise.complete.obs"),
    .names = "cor_{.col}"
  )) %>%
  pivot_longer(
    cols = starts_with("cor_"),
    names_prefix = "cor_",
    names_to = "Flu_Lag",
    values_to = "correlation"
  ) %>%
  mutate(state = factor(state, levels = sort(unique(state))))  # Order states alphabetically

nssp_corrs_rsv_lag = nssp_full_lagged %>%
  group_by(state) %>%
  summarise(across(
    starts_with("RSV_", ignore.case = FALSE),
    ~ cor(.x, flu_rescaled, use = "pairwise.complete.obs"),
    .names = "cor_{.col}"
  )) %>%
  pivot_longer(
    cols = starts_with("cor_"),
    names_prefix = "cor_",
    names_to = "RSV_Lag",
    values_to = "correlation"
  ) %>%
  mutate(state = factor(state, levels = sort(unique(state))))  # Order states alphabetically

# Correlation table between RSV and lagged flu data
ggplot(nssp_corrs_flu_lag, aes(
  x = Flu_Lag,
  y = factor(state, levels = sort(unique(state), decreasing = TRUE)),
  fill = correlation
)) +
  geom_tile(linewidth = 10) +
  geom_text(aes(label = round(correlation, 2)), color = "black", size = 4) +  # Adds correlation values
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0
  ) +
  labs(title = "Correlation between RSV and Lagged Flu Data",
       x = "RSV Lag",
       y = "State",
       fill = "Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = "~/Documents/Projects/rsv_peak_timing/figures/nssp_rsv_lagged_flu_corrs.png", height = 8, width = 10, units = "in", bg = "white")

## Correlation table between flu and lagged RSV
ggplot(nssp_corrs_rsv_lag, aes(
  x = RSV_Lag,
  y = factor(state, levels = sort(unique(state), decreasing = TRUE)),
  fill = correlation
)) +
  geom_tile(linewidth = 10) +
  geom_text(aes(label = round(correlation, 2)), color = "black", size = 4) +  # Adds correlation values
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0
  ) +
  labs(title = "Correlation between Flu and Lagged RSV Variables",
       x = "RSV Lag",
       y = "State",
       fill = "Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = "~/Documents/Projects/rsv_peak_timing/figures/nssp_flu_lagged_rsv_corrs.png", height = 8, width = 10, units = "in", bg = "white")

cor(nssp_full_lagged$flu_rescaled, nssp_full_lagged$RSV_1)
cor(nssp_full_lagged$flu_rescaled, nssp_full_lagged$RSV_2)
cor(nssp_full_lagged$flu_rescaled, nssp_full_lagged$RSV_3)
cor(nssp_full_lagged$flu_rescaled, nssp_full_lagged$RSV_4)

cor(nssp_full_lagged$rsv_rescaled, nssp_full_lagged$FLU_1)
cor(nssp_full_lagged$rsv_rescaled, nssp_full_lagged$FLU_2)
cor(nssp_full_lagged$rsv_rescaled, nssp_full_lagged$FLU_3)
cor(nssp_full_lagged$rsv_rescaled, nssp_full_lagged$FLU_4)


# Plot heatmap
ggplot(cor_df, aes(x = Lag, y = 1, fill = Correlation)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue",
                       high = "red",
                       mid = "white") +
  labs(title = "Heatmap of RSV vs. Flu Correlations",
       x = "Lag (Weeks)",
       y = "",
       fill = "Correlation") +
  theme_minimal()

nssp_full_ = nssp_full %>% select(week_end, state, season, ili_rescaled, contains('times_coef'))

# Saving the data for use in EWS
write_csv(nssp_full_, file = "data/processed/nssp_all_years_030425.csv")

