## -----------------------------------------------------------------------------
## Script name: evaluate_NSSP_signal_quality.R
##
## Purpose of script:
##
## Author: George Dewey
##
## Date Created: 2025-01-29
##
## Last Updated:
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
dat = read_csv("data/NSSP/NSSP_Emergency_Department_Visit_Trajectories_by_State_and_Sub_State_Regions-_COVID-19__Flu__RSV__Combined___20250129.csv")

## Load ILI data
## ILI data - need data from 21 - 23
ili_data = read_csv("data/ILInet/ILINet_state_2025-01-18.csv",
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
  filter(year %in% 2016:2025) %>%
  mutate(week_end = as_date(MMWRweek2Date(year, epiweek) + 6)) %>%
  select(region, year, epiweek, week_end, unweighted_ili, ilitotal) %>%
  rename(state = region)

dat = dat %>% filter(county == "All") %>% select(week_end, geography, percent_visits_covid, percent_visits_influenza, percent_visits_rsv)

# Create season variable
dat = dat %>%
  mutate(
    year = year(week_end),
    month = month(week_end),
    ili_season = ifelse(month >= 6,
                        paste0(year %% 100, "-", (year + 1) %% 100),  # If June or later: "YY-(YY+1)"
                        paste0((year - 1) %% 100, "-", year %% 100))   # If before June: "(YY-1)-YY"
  )

dat = dat %>% rename(state = geography,
               covid = percent_visits_covid,
               flu = percent_visits_influenza,
               rsv = percent_visits_rsv,
               season = ili_season) %>%
  select(week_end, state, covid, flu, rsv, season)

dat_combined = dat %>% left_join(ili_data_with_dates, by = c("state", "week_end")) %>%
  rename(ili = unweighted_ili) %>%
  select(-ilitotal)

minmax = function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

dat_combined_ = dat_combined %>% drop_na()

# Now run the penalized model for each season with states as the individual unit of analysis

# 22-23 season
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

  nssp_signals_22_23 = bind_rows(nssp_signals_22_23, signals_tmp)

}
nssp_signals_22_23 %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = flu_times_coef, color = 'Flu'), linewidth = 0.65) +
  geom_line(aes(y = rsv_times_coef, color = 'RSV'), linewidth = 0.65) +
  geom_line(aes(y = covid_times_coef, color = "COVID-19"), linewidth = 0.65) +
  geom_line(aes(y = ili_rescaled, color = 'ILI'), linewidth = 0.65) +
  facet_wrap( ~ state, scales = 'free', ncol = 10) +
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
ggsave(file = "figures/nssp_22_23.png", width = 16, height = 10, units = "in", bg = "white")

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

  nssp_signals_23_24 = bind_rows(nssp_signals_23_24, signals_tmp)

}
nssp_signals_23_24 %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = flu_times_coef, color = 'Flu'), linewidth = 0.65) +
  geom_line(aes(y = rsv_times_coef, color = 'RSV'), linewidth = 0.65) +
  geom_line(aes(y = covid_times_coef, color = "COVID-19"), linewidth = 0.65) +
  geom_line(aes(y = ili_rescaled, color = 'ILI'), linewidth = 0.65) +
  facet_wrap( ~ state, scales = 'free', ncol = 10) +
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
ggsave(file = "figures/nssp_23_24.png", width = 16, height = 10, units = "in", bg = "white")


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
  facet_wrap( ~ state, scales = 'free', ncol = 10) +
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

# Evaluating temporal relationships between RSV and Flu signals
nssp_signals_22_23 = nssp_signals_22_23 %>% select(-c(year, epiweek))

write_csv(nssp_signals_22_23, file = "data/nssp_signals_22_23_processed.csv")

# Testing out phase analysis

# First bind all the data together and create lags of the signals
nssp_full = bind_rows(nssp_signals_22_23, nssp_signals_23_24, nssp_signals_24_25)

nssp_full = nssp_full %>%
  mutate(vir_rhs = covid + flu + rsv,
         vir_rescaled_rhs = covid_rescaled + flu_rescaled + rsv_rescaled)

# Find optimal lag for ILI and the regression rhs
ccf(nssp_full$vir_rhs, nssp_full$ili_rescaled, lag.max = 8, plot = TRUE)
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

    FLU_0 = flu,
    FLU_1 = lag(flu, 1),
    FLU_2 = lag(flu, 2),
    FLU_3 = lag(flu, 3),
    FLU_4 = lag(flu, 4),

    # Create seasonal terms
    week_num = week(week_end),
    sin_annual = sin(2*pi*week_num/52),
    cos_annual = cos(2*pi*week_num/52),
    sin_semiannual = sin(4*pi*week_num/52),
    cos_semiannual = cos(4*pi*week_num/52)
  )

nssp_full_lagged = nssp_full_lagged %>% select(-c(year, epiweek))
nssp_full_lagged[is.na(nssp_full_lagged) == TRUE] = 0

# Fit the model
model <- nlme::lme(
  ili ~ RSV_0 + RSV_1 + RSV_2 + RSV_3 + RSV_4 +
    FLU_0 + FLU_1 + FLU_2 + FLU_3 + FLU_4 +
    sin_annual + cos_annual +
    sin_semiannual + cos_semiannual,
  random = ~ 1 | state/season,
  correlation = nlme::corAR1(),
  data = nssp_full_lagged
)

summary(model)
