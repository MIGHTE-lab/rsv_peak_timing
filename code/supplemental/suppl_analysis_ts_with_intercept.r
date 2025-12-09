## -----------------------------------------------------------------------------
## Script name: suppl_analysis_ts_with_intercept.R
##
## Purpose of script: Conduct secondary time series analysis on seasonal ILI data,
## including the model intercept term.
##
## Author: George Dewey
##
## Date Created: 2025-01-29
##
## Last Updated: 2025-09-29
##
## Version: 1.4.0
##
## Update notes: Modified to include and plot the intercept from ridge models.
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
library(caret)

## Set working directory
# setwd("~/Documents/Projects/rsv_peak_timing") # User should set their own directory

## Load NSSP data
dat = read_csv(
  "data/NSSP/nssp_full_040525.csv"
)

## Load ILI data and do some data management
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
      paste0((year - 1) %% 100, "-", year %% 100)
    )
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
# MODIFIED function to include the intercept in the model
ridge_CV = function(data) {

  predictors = c("covid_rescaled", "flu_rescaled", "rsv_rescaled")
  X = as.matrix(data[, predictors])
  y = data$ili_rescaled

  lambda_values = 10^seq(-5, 5, length = 100)

  # Retain the intercept
  cv_ridge = cv.glmnet(X, y, alpha = 0, lambda = lambda_values, nfolds = 10,
                       intercept = TRUE, lower.limits = 0)

  lambda_grid = seq(cv_ridge$lambda.min / 1000, cv_ridge$lambda.min * 1000, length = 500)

  ridge_models = lapply(lambda_grid, function(l) {
    # CHANGE: Set intercept = TRUE
    glmnet(X, y, alpha = 0, lambda = l, intercept = TRUE, lower.limits = 0)
  })

  ridge_coefs = sapply(ridge_models, function(model) coef(model)[,1])

  # NOTE: The number of predictors for this check is now the original number + intercept
  valid_lambda_idx = which(colSums(ridge_coefs != 0) == (length(predictors) + 1))[1]
  if (is.na(valid_lambda_idx)) {
    warning("No valid lambda found where all coefficients (incl. intercept) are nonzero. Using lambda.min.")
    best_lambda = cv_ridge$lambda.min
  } else {
    best_lambda = lambda_grid[valid_lambda_idx]
  }


  # CHANGE: Set intercept = TRUE for the final model
  ridge_final = glmnet(X, y, alpha = 0, lambda = best_lambda, intercept = TRUE, lower.limits = 0)

  coef_values = as.numeric(coef(ridge_final))

  # CHANGE: The returned vector now includes the intercept as the first element
  return(coef_values) # 1: intercept, 2: covid, 3: flu, 4: rsv
}

# 3. Run models and create state-season plots ----

# 22-23 season
nssp_signals_22_23_intercept = NULL
for (geo in dat_combined_ %>% distinct(state) %>% pull()) {

  data = dat_combined_ %>% filter(state == geo, season == "22-23") %>%
    mutate(
      ili_rescaled = minmax(ili),
      flu_rescaled = minmax(flu),
      rsv_rescaled = minmax(rsv),
      covid_rescaled = minmax(covid))

  coefs_tmp = ridge_CV(data)

  signals_tmp = data %>%
    mutate(
      # CHANGE: Extract intercept and adjust indices for other coefficients
      beta_intercept = coefs_tmp[1],
      beta_covid = coefs_tmp[2],
      beta_flu = coefs_tmp[3],
      beta_rsv = coefs_tmp[4],
      flu_times_coef = flu_rescaled * beta_flu,
      rsv_times_coef = rsv_rescaled * beta_rsv,
      covid_times_coef = covid_rescaled * beta_covid,
    )

  nssp_signals_22_23_intercept = bind_rows(nssp_signals_22_23_intercept, signals_tmp)
  print(paste0("State ",geo," completed."))
}

# Plots including the intercept terms
nssp_signals_22_23_intercept %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = flu_times_coef, color = 'Flu'), linewidth = 0.65) +
  geom_line(aes(y = rsv_times_coef, color = 'RSV'), linewidth = 0.65) +
  geom_line(aes(y = covid_times_coef, color = "COVID-19"), linewidth = 0.65) +
  geom_line(aes(y = ili_rescaled, color = 'ILI'), linewidth = 0.65) +
  geom_line(aes(y = beta_intercept, color = 'Intercept'), linetype = "dashed", linewidth = 0.65) +
  facet_wrap(~ state, scales = 'free', ncol = 10) +
  scale_color_manual(values = c(
    'ILI' = '#E2DBC9',
    'RSV' = '#183A5A',
    'Flu' = '#C34129',
    'COVID-19' = "#EFB75B",
    'Intercept' = 'grey30'
  ),
  breaks = c('ILI', 'Flu', 'RSV', 'COVID-19', 'Intercept')
  ) +
  labs(x = '', y = 'Rescaled Volume / Contribution') +
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
ggsave(file = "figures/manuscript_figures/supplementary/nssp_22_23_intercept_100125.png", width = 16, height = 10, units = "in", bg = "white")

# 23-24 season
nssp_signals_23_24_intercept = NULL
for (geo in dat_combined_ %>% distinct(state) %>% pull()) {

  data = dat_combined_ %>% filter(state == geo, season == "23-24") %>%
    mutate(
      ili_rescaled = minmax(ili),
      flu_rescaled = minmax(flu),
      rsv_rescaled = minmax(rsv),
      covid_rescaled = minmax(covid))

  coefs_tmp = ridge_CV(data)

  signals_tmp = data %>%
    mutate(
      beta_intercept = coefs_tmp[1],
      beta_covid = coefs_tmp[2],
      beta_flu = coefs_tmp[3],
      beta_rsv = coefs_tmp[4],
      flu_times_coef = flu_rescaled * beta_flu,
      rsv_times_coef = rsv_rescaled * beta_rsv,
      covid_times_coef = covid_rescaled * beta_covid,
    )

  nssp_signals_23_24_intercept = bind_rows(nssp_signals_23_24_intercept, signals_tmp)
  print(paste0("State ",geo," completed."))
}

# MODIFIED plot to include the intercept term
nssp_signals_23_24_intercept %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = flu_times_coef, color = 'Flu'), linewidth = 0.65) +
  geom_line(aes(y = rsv_times_coef, color = 'RSV'), linewidth = 0.65) +
  geom_line(aes(y = covid_times_coef, color = "COVID-19"), linewidth = 0.65) +
  geom_line(aes(y = ili_rescaled, color = 'ILI'), linewidth = 0.65) +
  geom_line(aes(y = beta_intercept, color = 'Intercept'), linetype = "dashed", linewidth = 0.65) +
  facet_wrap(~ state, scales = 'free', ncol = 10) +
  scale_color_manual(values = c(
    'ILI' = '#E2DBC9',
    'RSV' = '#183A5A',
    'Flu' = '#C34129',
    'COVID-19' = "#EFB75B",
    'Intercept' = 'grey30'
  ),
  breaks = c('ILI', 'Flu', 'RSV', 'COVID-19', 'Intercept')
  ) +
  labs(x = '', y = 'Rescaled Volume / Contribution') +
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
ggsave(file = "figures/manuscript_figures/supplementary/nssp_23_24_intercept_100125.png", width = 16, height = 10, units = "in", bg = "white")

# 24-25 season
nssp_signals_24_25_intercept = NULL
for (geo in dat_combined_ %>% distinct(state) %>% pull()) {

  data = dat_combined_ %>% filter(state == geo, season == "23-24") %>%
    mutate(
      ili_rescaled = minmax(ili),
      flu_rescaled = minmax(flu),
      rsv_rescaled = minmax(rsv),
      covid_rescaled = minmax(covid))

  coefs_tmp = ridge_CV(data)

  signals_tmp = data %>%
    mutate(
      # CHANGE: Extract intercept and adjust indices for other coefficients
      beta_intercept = coefs_tmp[1],
      beta_covid = coefs_tmp[2],
      beta_flu = coefs_tmp[3],
      beta_rsv = coefs_tmp[4],
      flu_times_coef = flu_rescaled * beta_flu,
      rsv_times_coef = rsv_rescaled * beta_rsv,
      covid_times_coef = covid_rescaled * beta_covid,
    )

  nssp_signals_24_25_intercept = bind_rows(nssp_signals_24_25_intercept, signals_tmp)
  print(paste0("State ",geo," completed."))
}

# MODIFIED plot to include the intercept term
nssp_signals_24_25_intercept %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = flu_times_coef, color = 'Flu'), linewidth = 0.65) +
  geom_line(aes(y = rsv_times_coef, color = 'RSV'), linewidth = 0.65) +
  geom_line(aes(y = covid_times_coef, color = "COVID-19"), linewidth = 0.65) +
  geom_line(aes(y = ili_rescaled, color = 'ILI'), linewidth = 0.65) +
  geom_line(aes(y = beta_intercept, color = 'Intercept'), linetype = "dashed", linewidth = 0.65) +
  facet_wrap(~ state, scales = 'free', ncol = 10) +
  scale_color_manual(values = c(
    'ILI' = '#E2DBC9',
    'RSV' = '#183A5A',
    'Flu' = '#C34129',
    'COVID-19' = "#EFB75B",
    'Intercept' = 'grey30'
  ),
  breaks = c('ILI', 'Flu', 'RSV', 'COVID-19', 'Intercept')
  ) +
  labs(x = '', y = 'Rescaled Volume / Contribution') +
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
ggsave(file = "figures/manuscript_figures/supplementary/nssp_24_25_intercept_100125.png", width = 16, height = 10, units = "in", bg = "white")


