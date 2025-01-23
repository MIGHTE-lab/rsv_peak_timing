## -----------------------------------------------------------------------------
## Script name: compare_ID_peaks.R
##
## Purpose of script: To visualize the ILI, Flu, RSV, and COVID signals to detect
## the sequential pattern of peaks identified through early warning systems
##
## Author: George Dewey
##
## Date Created: 2024-07-01
##
## Last Updated: 2025-01-23
##
## Version: 1.4.1
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

## Load RSVnet data
rsvnet = read_csv("data/RSV/RSV-NET_data/rsv_hosps_010825.csv")

## NSSP (RSV) data -- for this version, unused.
# nssp_data = read_csv("data/ground_truth/RSV/NSSP/NSSP_Combined_07-05-24.csv")

## ILI data - need data from 21 - 23
ili_data = read_csv("data/ILInet/ILINet_state_2024-11-30.csv", na = "X",
                    skip = 1)

## FluSurv Hospitalization data
flu_data = read_csv("data/FluSURV/fluSurv_RSVNET_16_24.csv",
                    na = "null", skip = 2, n_max = 99866)

## COVID data - data limited by available dates
covid_data = read_csv("data/COVID/COVID_aggregated.csv")

## 2. Data cleaning ------------------------------------------------------------

## Fixing variable names
names(rsvnet) = c("state", "season", "week_end", "age_cat", "sex", "race", "rate",
                  "cum_rate", "type")

names(ili_data) = c("region_type", "region", "year", "epiweek", "weighted_ili",
                    "unweighted_ili", "age_0_4", "age_25_49", "age_25_64",
                    "age_5_24", "age_50_64", "age_65", "ilitotal",
                    "num_data_providers",
                    "total_patients")

names(flu_data) = c("state", "network", "season", "year", "epiweek", "age_cat",
                    "sex_cat", "race_cat", "virus", "cum_rate", "weekly_rate",
                    "age_adj_cum_rate", "age_adj_weekly_rate")

## Unify the date format and start from the same week
##
## Based on evaluation of the dataset start, dates are listed in two formats:
## either the week_end date or the MMWR week format (i.e., year + week number,)
## where the week starts on Sunday and ends on Saturday. To proceed with time
## series analysis, we need to unify each of the data inputs with the same date
## format. For easiest interpretation (and consistency with the EWS projects) we
## will unify on the the end of week date (i.e., each data point should be
## associated with the Saturday of the week in which the data was collected -
## the end of the MMWR week).
##
# Since the NSSP data is already in week_end format, we don"t need to do any
# pruning
MMWRweek2Date(2022, 39)+6
# confirms that NSSP data starts on MMWR week 39 of 2022

# NSSP data starts at MMWR week 39 of 2022. We eliminate the first week and
# convert the date column to week_end date. The MMWRweek2date function converts
# MMWR week/epiweek to the ****START OF WEEK/SUNDAY*** date, so we need to add
# 6 days to line the data up with the NSSP data.
# ili_data_short = ili_data %>%
#   select(region, year, epiweek, unweighted_ili, total_patients)

rsvnet_adults = rsvnet %>%
  filter(age_cat == "≥18 years (Adults)", race == "All", sex == "All",
         state != "RSV-NET")

# We need ILI starting 2016-2017 season - for the most part ILI data is not age
# stratified so should be OK to look at the overall peak

# One options - maybe exclude

## ILI data - clean so that we have dates + regional identifier (state)
ili_data_with_dates = ili_data %>%
  filter(year %in% 2016:2024) %>%
  mutate(week_end = as_date(MMWRweek2Date(year, epiweek) + 6)) %>%
  select(region, year, epiweek, week_end, unweighted_ili, ilitotal) %>%
  rename(state = region)

rsvnet_adults %>% filter(week_end == as_date("2022-12-17"))

rsvnet_adults %>% distinct(state)

# Create the figure 1 dataset by merging ILI with RSV
rsv_ili = rsvnet_adults %>%
  left_join(ili_data_with_dates, by = c("state", "week_end"))

# Now fix the COVID data because we need it for 19-22

covid_data = covid_data %>%
  mutate(week_end = MMWRweek2Date(year, epiweek) + 6) %>%
  select(week_end, epiweek, year, cases, geography)

# We also need flu-speciifc data from FluSurv for all data, unify the categories
# do it for all adults (age >= 18, all genders, all races

# NOTE: specific virus confirmations for Flu A and Flu B are not age-restricted
# Better to keep overall status vs. specific flu strains?

flu_data_adults = flu_data %>%
  filter(age_cat == ">= 18", sex_cat == "Overall", race_cat == "Overall",
         !state %in% c("New York - Albany", "New York - Rochester")) %>%
  mutate(week_end = as_date(MMWRweek2Date(year, epiweek) + 6)) %>%
  select(state, week_end, weekly_rate) %>%
  rename(flu = weekly_rate)

flu_data %>% distinct(state)

rsvnet_adults = rsvnet_adults %>%
  select(state, week_end, rate) %>%
  rename(rsv = rate)

ili_data_with_dates = ili_data_with_dates %>%
  select(state, week_end, unweighted_ili) %>%
  rename(ili = unweighted_ili)

ili_flu_rsv = rsvnet_adults %>%
  left_join(ili_data_with_dates, by = c("state", "week_end")) %>%
  left_join(flu_data_adults, by = c("state", "week_end"))

# To do the regression properly, we should consider times when Flu Surveillance
# was not recorded (i.e., flu off-season) as 0
ili_flu_rsv[is.na(ili_flu_rsv) == TRUE] = 0

# Drop NY and NC because of data issues
ili_flu_rsv = ili_flu_rsv %>% filter(!state %in% c("New York", "North Carolina"))

# The date ranges vary for individual 11 retained RSVNET states: let's check each
# The states that got excluded are:  NY, NC, MD


rsvnet_adults %>% distinct(state)

# First, get the list of all unique states in the remaining dataset
ili_flu_rsv %>% distinct(state)
# Utah, Tennessee, Oregon, New Mexico, Minnesota, Michigan, Maryland
# Georgia, Connecticut, Colorado, California

# Next, rescale data for each state as needed

rsvnet_states = c("Utah", "Tennessee", "Oregon", "New Mexico", "Minnesota",
                  "Michigan", "Georgia", "Connecticut",
                  "Colorado", "California")

# We want to replicate the penalized regression process that we used before but
# for each year

# What happens if we do the regression for all years and not just 1?
# It appears yearly modeling looks better

# Make sure all NA values (from off-season times or other data artifacts) are 0
ili_flu_rsv[is.na(ili_flu_rsv) == TRUE] = 0

# Add a season indicator
ili_flu_rsv = ili_flu_rsv %>% mutate(
    season = case_when(
      week_end >= as.Date(paste0(format(week_end, "%Y"), "-06-01")) &
        week_end < as.Date(paste0(format(week_end + months(12), "%Y"), "-05-31")) ~ paste0(format(week_end, "%y"), "-", format(week_end + years(1), "%y")),
      week_end >= as.Date(paste0(format(week_end - months(6), "%Y"), "-06-01")) &
        week_end < as.Date(paste0(format(week_end, "%Y"), "-05-31")) ~ paste0(format(week_end - years(1), "%y"), "-", format(week_end, "%y")),
      TRUE ~ NA_character_
    ))

# Some additional data processing
ili_flu_rsv = ili_flu_rsv %>%
  mutate(season = ifelse(week_end == as_date("2018-03-31") & is.na(season) == TRUE, "17-18", season))

ili_flu_rsv = ili_flu_rsv %>% filter(state %in% rsvnet_states)

ili_flu_rsv = ili_flu_rsv %>% mutate(season = ifelse(season == "NA-20", "19-20", season))

seasons = c("16-17", "17-18", "18-19", "19-20", "20-21", "21-22", "22-23", "23-24")

# 3. Modeling and data aggregation ---------------------------------------------

# 2016-2017 season ----
signals_with_betas_16_17 = NULL
for(geo in ili_flu_rsv %>% filter(season == "16-17") %>% distinct(state) %>% pull()){

  # Filter the data to only include data from the selected state and season
  data = ili_flu_rsv %>% filter(state == geo, season == "16-17") %>%
      mutate(ili_rescaled = ili/max(ili),
             flu_rescaled = flu/max(flu),
             rsv_rescaled = rsv/max(rsv))

  # Use 10-fold CV to identify best lambda
  cv = cv.glmnet(x = data.matrix(data[,c('flu_rescaled', 'rsv_rescaled')]),
                 y = data.matrix(data[,'ili_rescaled']),
                 intercept = FALSE, alpha = 0)

  # Use penalized package to run the ridge w/ no intercept + non-negative
  # and the specified lambda from the CV model
  fit = penalized(ili_rescaled ~ flu_rescaled + rsv_rescaled ,
                  positive = TRUE,
                  unpenalized = ~0,
                  data = data,
                  lambda2 = cv$lambda.min)

  coefs_tmp = coefficients(fit)

  signals_tmp = data %>%
    mutate(beta_flu = coefs_tmp[1],
           beta_rsv = coefs_tmp[2],
           flu_times_coef = flu_rescaled*beta_flu,
           rsv_times_coef = rsv_rescaled*beta_rsv)

  signals_with_betas_16_17 = bind_rows(signals_with_betas_16_17, signals_tmp)

}
signals_with_betas_16_17
# Generate the plot year by year - maybe the ordering is not that consistent

# The criteria we are looking for is: RSV's MAX VALUE occurs before Flu's max value
# The volume of intensity will be lower for RSV since the overall percentages are
# lower across the board.

# That is, in mathematical terms t(max RSV) < t(max flu) or t(max covid)

# 2016-2017 - RSV first in 2/5 state-seasons, flu and RSV synchronous in 3/5 states
signals_with_betas_16_17 %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = flu_times_coef, color = 'Flu'), linewidth = 0.65, na.rm = TRUE) +
  geom_line(aes(y = rsv_times_coef, color = 'RSV'), linewidth = 0.65, na.rm = TRUE) +
  geom_line(aes(y = ili_rescaled, color = 'ILI'), linewidth = 0.65, na.rm = TRUE) +
  facet_grid2(~ state, scales = "free_x", independent = "x") +
  scale_color_manual(values = c('ILI' = 'grey',
                                'RSV' = 'dodgerblue4',
                                'Flu' = 'darkred')) +
  labs(x = '', y = '') +
  envalysis::theme_publish() +
  theme(
    axis.text.x = element_text(size = 5),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    legend.title = element_blank(),
    legend.position = 'right',
    strip.text = element_text(size = 10, face = "bold")  # State labels formatting
  )



the new code nop1 = signals_with_betas_16_17  %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = flu_times_coef, color = 'Flu'), linewidth = 0.65) +
  geom_line(aes(y = rsv_times_coef, color = 'RSV'), linewidth = 0.65) +
  geom_line(aes(y = ili_rescaled, color = 'ILI'), linewidth = 0.65) +
  facet_wrap(~state, scales = 'free', nrow = 1) +
  scale_color_manual(values = c('ILI' = 'grey',
                                'RSV' = 'dodgerblue4',
                                'Flu' = 'darkred')) +
  labs(x = '', y = '') +
  envalysis::theme_publish() +
  theme(axis.text.x = element_text(size = 5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        legend.position = 'right')

ggsave(file = "figures/16_17_peak_timing.png", height = 6, width = 10, units = "in", bg = "white")

# 2017-2018 season ----
signals_with_betas_17_18 = NULL
for(geo in ili_flu_rsv %>% filter(season == "17-18") %>% distinct(state) %>% pull()){

  # Filter the data to only include data from the selected state and season
  data = ili_flu_rsv %>% filter(state == geo, season == "17-18") %>%
    mutate(ili_rescaled = ili/max(ili),
           flu_rescaled = flu/max(flu),
           rsv_rescaled = rsv/max(rsv))

  # Use 10-fold CV to identify best lambda
  cv = cv.glmnet(x = data.matrix(data[,c('flu_rescaled', 'rsv_rescaled')]),
                 y = data.matrix(data[,'ili_rescaled']),
                 intercept = FALSE, alpha = 0)

  # Use penalized package to run the ridge w/ no intercept + non-negative
  # and the specified lambda from the CV model
  fit = penalized(ili_rescaled ~ flu_rescaled + rsv_rescaled ,
                  positive = TRUE,
                  unpenalized = ~0,
                  data = data,
                  lambda2 = cv$lambda.min)

  coefs_tmp = coefficients(fit)

  signals_tmp = data %>%
    mutate(beta_flu = coefs_tmp[1],
           beta_rsv = coefs_tmp[2],
           flu_times_coef = flu_rescaled*beta_flu,
           rsv_times_coef = rsv_rescaled*beta_rsv)

  signals_with_betas_17_18 = bind_rows(signals_with_betas_17_18, signals_tmp)

}
signals_with_betas_17_18  %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = flu_times_coef, color = 'Flu'), linewidth = 0.65) +
  geom_line(aes(y = rsv_times_coef, color = 'RSV'), linewidth = 0.65) +
  geom_line(aes(y = ili_rescaled, color = 'ILI'), linewidth = 0.65) +
  facet_wrap(~state, scales = 'free', nrow = 1) +
  scale_color_manual(values = c('ILI' = 'grey',
                                'RSV' = 'dodgerblue4',
                                'Flu' = 'darkred')) +
  labs(x = '', y = '') +
  envalysis::theme_publish() +
  theme(axis.text.x = element_text(size = 5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        legend.position = 'right')
ggsave(file = "figures/17_18_peak_timing.png", height = 6, width = 10, units = "in", bg = "white")
# RSV first in 1/6
# Flu first in 3/6
# Synchronous in 2/6

# 2018-2019 season ----
signals_with_betas_18_19 = NULL
for(geo in ili_flu_rsv %>% filter(season == "18-19") %>% distinct(state) %>% pull()){

  # Filter the data to only include data from the selected state and season
  data = ili_flu_rsv %>% filter(state == geo, season == "18-19") %>%
    mutate(ili_rescaled = ili/max(ili),
           flu_rescaled = flu/max(flu),
           rsv_rescaled = rsv/max(rsv))

  # Use 10-fold CV to identify best lambda
  cv = cv.glmnet(x = data.matrix(data[,c('flu_rescaled', 'rsv_rescaled')]),
                 y = data.matrix(data[,'ili_rescaled']),
                 intercept = FALSE, alpha = 0)

  # Use penalized package to run the ridge w/ no intercept + non-negative
  # and the specified lambda from the CV model
  fit = penalized(ili_rescaled ~ flu_rescaled + rsv_rescaled ,
                  positive = TRUE,
                  unpenalized = ~0,
                  data = data,
                  lambda2 = cv$lambda.min)

  coefs_tmp = coefficients(fit)

  signals_tmp = data %>%
    mutate(beta_flu = coefs_tmp[1],
           beta_rsv = coefs_tmp[2],
           flu_times_coef = flu_rescaled*beta_flu,
           rsv_times_coef = rsv_rescaled*beta_rsv)

  signals_with_betas_18_19 = bind_rows(signals_with_betas_18_19, signals_tmp)

}
signals_with_betas_18_19  %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = flu_times_coef, color = 'Flu'), linewidth = 0.65) +
  geom_line(aes(y = rsv_times_coef, color = 'RSV'), linewidth = 0.65) +
  geom_line(aes(y = ili_rescaled, color = 'ILI'), linewidth = 0.65) +
  facet_wrap(~state, scales = 'free', ncol = 5) +
  scale_color_manual(values = c('ILI' = 'grey',
                                'RSV' = 'dodgerblue4',
                                'Flu' = 'darkred')) +
  labs(x = '', y = '') +
  envalysis::theme_publish() +
  theme(axis.text.x = element_text(size = 5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        legend.position = 'right')
ggsave(file = "figures/18_19_peak_timing.png", height = 6, width = 10, units = "in", bg = "white")

# Flu first in 4/10
# RSV first in 6/10

# 2019-2020 season ----

# Starting here, we should merge in the COVID data
ili_flu_rsv_covid = ili_flu_rsv %>% left_join(covid_data, by = c("week_end", "state" = "geography")) %>%
  select(state, week_end, rsv, ili, flu, season, cases) %>%
  rename(covid = cases)

# Zero out NAs (meaning either COVID data was not collected or no COVID was present)
ili_flu_rsv_covid[is.na(ili_flu_rsv_covid) == TRUE] = 0

# We have COVID data for 5 seasons (but not for every state)

# I think the plan of action should be to use COVID data whenever possible

signals_with_betas_19_20 = NULL
for(geo in ili_flu_rsv_covid %>% filter(season == "19-20") %>% distinct(state) %>% pull()){

  # Filter the data to only include data from the selected state and season
  data = ili_flu_rsv_covid %>% filter(state == geo, season == "19-20") %>%
    mutate(ili_rescaled = ili/max(ili),
           flu_rescaled = flu/max(flu),
           rsv_rescaled = rsv/max(rsv),
           covid_rescaled = covid/max(covid))

  # Use 10-fold CV to identify best lambda
  cv = cv.glmnet(x = data.matrix(data[,c('flu_rescaled', 'rsv_rescaled', 'covid_rescaled')]),
                 y = data.matrix(data[,'ili_rescaled']),
                 intercept = FALSE, alpha = 0)

  # Use penalized package to run the ridge w/ no intercept + non-negative
  # and the specified lambda from the CV model
  fit = penalized(ili_rescaled ~ flu_rescaled + rsv_rescaled + covid_rescaled,
                  positive = TRUE,
                  unpenalized = ~0,
                  data = data,
                  lambda2 = cv$lambda.min)

  coefs_tmp = coefficients(fit)

  signals_tmp = data %>%
    mutate(beta_flu = coefs_tmp[1],
           beta_rsv = coefs_tmp[2],
           beta_covid = coefs_tmp[3],
           flu_times_coef = flu_rescaled*beta_flu,
           rsv_times_coef = rsv_rescaled*beta_rsv,
           covid_times_coef = covid_rescaled*beta_covid)

  signals_with_betas_19_20 = bind_rows(signals_with_betas_19_20, signals_tmp)

}
signals_with_betas_19_20  %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = flu_times_coef, color = 'Flu'), linewidth = 0.65) +
  geom_line(aes(y = rsv_times_coef, color = 'RSV'), linewidth = 0.65) +
  geom_line(aes(y = covid_times_coef, color = "COVID-19"), linewidth = 0.65) +
  geom_line(aes(y = ili_rescaled, color = 'ILI'), linewidth = 0.65) +
  facet_wrap(~state, scales = 'free', ncol = 5) +
  scale_color_manual(values = c('ILI' = '#E2DBC9',
                                'RSV' = '#183A5A',
                                'Flu' = '#C34129',
                                'COVID-19' = "#EFB75B")) +
  labs(x = '', y = '') +
  envalysis::theme_publish() +
  theme(axis.text.x = element_text(size = 5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        legend.position = 'right')
# RSV in 5/10
ggsave(file = "figures/19_20_peak_timing.png", height = 6, width = 10, units = "in", bg = "white")

ili_flu_rsv_covid %>% filter(season == "20-21") %>% view()

rsvnet_adults %>%
  filter(week_end >= as_date("2020-05-01") & week_end <= as_date("2021-06-01")) %>%
  arrange(state) %>%
  view()

# 2020-2021 season ----
signals_with_betas_20_21 = NULL
for(geo in ili_flu_rsv_covid %>% filter(season == "20-21") %>% distinct(state) %>% pull()){

  # Filter the data to only include data from the selected state and season
  data = ili_flu_rsv_covid %>% filter(state == geo, season == "20-21") %>%
    mutate(ili_rescaled = ili/max(ili, na.rm = TRUE),
           flu_rescaled = flu/max(flu, na.rm = TRUE),
           rsv_rescaled = rsv/max(rsv, na.rm = TRUE),
           covid_rescaled = covid/max(covid, na.rm = TRUE)) %>%
    mutate(across(everything(), ~ replace_na(., 0)))

  # Use 10-fold CV to identify best lambda
  cv = cv.glmnet(x = data.matrix(data[,c('flu_rescaled', 'rsv_rescaled', 'covid_rescaled')]),
                 y = data.matrix(data[,'ili_rescaled']),
                 intercept = FALSE, alpha = 0)

  # Use penalized package to run the ridge w/ no intercept + non-negative
  # and the specified lambda from the CV model
  fit = penalized(ili_rescaled ~ flu_rescaled + rsv_rescaled + covid_rescaled,
                  positive = TRUE,
                  unpenalized = ~0,
                  data = data,
                  lambda2 = cv$lambda.min)

  coefs_tmp = coefficients(fit)

  signals_tmp = data %>%
    mutate(beta_flu = coefs_tmp[1],
           beta_rsv = coefs_tmp[2],
           beta_covid = coefs_tmp[3],
           flu_times_coef = flu_rescaled*beta_flu,
           rsv_times_coef = rsv_rescaled*beta_rsv,
           covid_times_coef = covid_rescaled*beta_covid)

  signals_with_betas_20_21 = bind_rows(signals_with_betas_20_21, signals_tmp)

}
signals_with_betas_20_21 %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = flu_times_coef, color = 'Flu'), linewidth = 0.65) +
  geom_line(aes(y = rsv_times_coef, color = 'RSV'), linewidth = 0.65) +
  geom_line(aes(y = covid_times_coef, color = "COVID-19"), linewidth = 0.65) +
  geom_line(aes(y = ili_rescaled, color = 'ILI'), linewidth = 0.65) +
  facet_wrap(~state, scales = 'free', ncol = 5) +
  scale_color_manual(values = c('ILI' = '#E2DBC9',
                                'RSV' = '#183A5A',
                                'Flu' = '#C34129',
                                'COVID-19' = "#EFB75B")) +
  labs(x = '', y = '') +
  envalysis::theme_publish() +
  theme(axis.text.x = element_text(size = 5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        legend.position = 'right')

# 21-22 season ----
signals_with_betas_21_22 = NULL
for(geo in ili_flu_rsv_covid %>% filter(season == "21-22") %>% distinct(state) %>% pull()){

  # Filter the data to only include data from the selected state and season
  data = ili_flu_rsv_covid %>% filter(state == geo, season == "21-22") %>%
    mutate(ili_rescaled = ili/max(ili, na.rm = TRUE),
           flu_rescaled = flu/max(flu, na.rm = TRUE),
           rsv_rescaled = rsv/max(rsv, na.rm = TRUE),
           covid_rescaled = covid/max(covid, na.rm = TRUE)) %>%
    mutate(across(everything(), ~ replace_na(., 0)))

  # Use 10-fold CV to identify best lambda
  cv = cv.glmnet(x = data.matrix(data[,c('flu_rescaled', 'rsv_rescaled', 'covid_rescaled')]),
                 y = data.matrix(data[,'ili_rescaled']),
                 intercept = FALSE, alpha = 0)
  # cv = cv.glmnet(x = data.matrix(data[,c('flu_rescaled', 'rsv_rescaled', 'covid_rescaled')]),
  #                y = data.matrix(data[,'ili_rescaled']),
  #                intercept = FALSE, alpha = 0)

  # Use penalized package to run the ridge w/ no intercept + non-negative
  # and the specified lambda from the CV model
  fit = penalized(ili_rescaled ~ flu_rescaled + rsv_rescaled + covid_rescaled,
                  positive = TRUE,
                  unpenalized = ~0,
                  data = data,
                  lambda2 = cv$lambda.min)

  coefs_tmp = coefficients(fit)

  signals_tmp = data %>%
    mutate(beta_flu = coefs_tmp[1],
           beta_rsv = coefs_tmp[2],
           beta_covid = coefs_tmp[3],
           flu_times_coef = flu_rescaled*beta_flu,
           rsv_times_coef = rsv_rescaled*beta_rsv,
           covid_times_coef = covid_rescaled*beta_covid)

  signals_with_betas_21_22 = bind_rows(signals_with_betas_21_22, signals_tmp)

}
signals_with_betas_21_22 %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = flu_times_coef, color = 'Flu'), linewidth = 0.65) +
  geom_line(aes(y = rsv_times_coef, color = 'RSV'), linewidth = 0.65) +
  geom_line(aes(y = covid_times_coef, color = "COVID-19"), linewidth = 0.65) +
  geom_line(aes(y = ili_rescaled, color = 'ILI'), linewidth = 0.65) +
  facet_wrap(~state, scales = 'free', ncol = 5) +
  scale_color_manual(values = c('ILI' = '#E2DBC9',
                                'RSV' = '#183A5A',
                                'Flu' = '#C34129',
                                'COVID-19' = "#EFB75B")) +
  labs(x = '', y = '') +
  scale_x_date(date_labels = "%b\n%y",
               date_breaks = "2 months") +
  envalysis::theme_publish() +
  theme(axis.text.x = element_text(size = 5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        legend.position = 'right')
ggsave(file = "figures/21_22_peak_timing.png", height = 6, width = 12, units = "in", bg = "white")

signals_with_betas_22_23 = NULL
for(geo in ili_flu_rsv_covid %>% filter(season == "22-23") %>% distinct(state) %>% pull()){

  # Filter the data to only include data from the selected state and season
  data = ili_flu_rsv_covid %>% filter(state == geo, season == "22-23") %>%
    mutate(ili_rescaled = ili/max(ili, na.rm = TRUE),
           flu_rescaled = flu/max(flu, na.rm = TRUE),
           rsv_rescaled = rsv/max(rsv, na.rm = TRUE),
           covid_rescaled = covid/max(covid, na.rm = TRUE)) %>%
    mutate(across(everything(), ~ replace_na(., 0)))

  # Use 10-fold CV to identify best lambda
  cv = cv.glmnet(x = data.matrix(data[,c('flu_rescaled', 'rsv_rescaled', 'covid_rescaled')]),
                 y = data.matrix(data[,'ili_rescaled']),
                 intercept = FALSE, alpha = 0)
  # cv = cv.glmnet(x = data.matrix(data[,c('flu_rescaled', 'rsv_rescaled', 'covid_rescaled')]),
  #                y = data.matrix(data[,'ili_rescaled']),
  #                intercept = FALSE, alpha = 0)

  # Use penalized package to run the ridge w/ no intercept + non-negative
  # and the specified lambda from the CV model
  fit = penalized(ili_rescaled ~ flu_rescaled + rsv_rescaled + covid_rescaled,
                  positive = TRUE,
                  unpenalized = ~0,
                  data = data,
                  lambda2 = cv$lambda.min)

  coefs_tmp = coefficients(fit)

  signals_tmp = data %>%
    mutate(beta_flu = coefs_tmp[1],
           beta_rsv = coefs_tmp[2],
           beta_covid = coefs_tmp[3],
           flu_times_coef = flu_rescaled*beta_flu,
           rsv_times_coef = rsv_rescaled*beta_rsv,
           covid_times_coef = covid_rescaled*beta_covid)

  signals_with_betas_22_23 = bind_rows(signals_with_betas_22_23, signals_tmp)

}
signals_with_betas_22_23 %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = flu_times_coef, color = 'Flu'), linewidth = 0.65) +
  geom_line(aes(y = rsv_times_coef, color = 'RSV'), linewidth = 0.65) +
  geom_line(aes(y = covid_times_coef, color = "COVID-19"), linewidth = 0.65) +
  geom_line(aes(y = ili_rescaled, color = 'ILI'), linewidth = 0.65) +
  facet_wrap(~state, scales = 'free', ncol = 5) +
  scale_color_manual(values = c('ILI' = '#E2DBC9',
                                'RSV' = '#183A5A',
                                'Flu' = '#C34129',
                                'COVID-19' = "#EFB75B")) +
  labs(x = '', y = '') +
  scale_x_date(date_labels = "%b\n%y",
               date_breaks = "2 months") +
  envalysis::theme_publish() +
  theme(axis.text.x = element_text(size = 5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        legend.position = 'right')
ggsave("figures/22_23_peak_timing.png", height = 6, width = 12, units = "in", bg = "white")

# 23-24 season ----
signals_with_betas_23_24 = NULL
for(geo in ili_flu_rsv_covid %>% filter(season == "23-24") %>% distinct(state) %>% pull()){

  # Filter the data to only include data from the selected state and season
  data = ili_flu_rsv_covid %>% filter(state == geo, season == "23-24") %>%
    mutate(ili_rescaled = ili/max(ili, na.rm = TRUE),
           flu_rescaled = flu/max(flu, na.rm = TRUE),
           rsv_rescaled = rsv/max(rsv, na.rm = TRUE),
           covid_rescaled = covid/max(covid, na.rm = TRUE)) %>%
    mutate(across(everything(), ~ replace_na(., 0)))

  # Use 10-fold CV to identify best lambda
  cv = cv.glmnet(x = data.matrix(data[,c('flu_rescaled', 'rsv_rescaled', 'covid_rescaled')]),
                 y = data.matrix(data[,'ili_rescaled']),
                 intercept = FALSE, alpha = 0)
  # cv = cv.glmnet(x = data.matrix(data[,c('flu_rescaled', 'rsv_rescaled', 'covid_rescaled')]),
  #                y = data.matrix(data[,'ili_rescaled']),
  #                intercept = FALSE, alpha = 0)

  # Use penalized package to run the ridge w/ no intercept + non-negative
  # and the specified lambda from the CV model
  fit = penalized(ili_rescaled ~ flu_rescaled + rsv_rescaled + covid_rescaled,
                  positive = TRUE,
                  unpenalized = ~0,
                  data = data,
                  lambda2 = cv$lambda.min)

  coefs_tmp = coefficients(fit)

  signals_tmp = data %>%
    mutate(beta_flu = coefs_tmp[1],
           beta_rsv = coefs_tmp[2],
           beta_covid = coefs_tmp[3],
           flu_times_coef = flu_rescaled*beta_flu,
           rsv_times_coef = rsv_rescaled*beta_rsv,
           covid_times_coef = covid_rescaled*beta_covid)

  signals_with_betas_23_24 = bind_rows(signals_with_betas_23_24, signals_tmp)

}
signals_with_betas_23_24 %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = flu_times_coef, color = 'Flu'), linewidth = 0.65) +
  geom_line(aes(y = rsv_times_coef, color = 'RSV'), linewidth = 0.65) +
  # geom_line(aes(y = covid_times_coef, color = "COVID-19"), linewidth = 0.65) +
  geom_line(aes(y = ili_rescaled, color = 'ILI'), linewidth = 0.65) +
  facet_wrap(~state, scales = 'free', ncol = 5) +
  scale_color_manual(values = c('ILI' = '#E2DBC9',
                                'RSV' = '#183A5A',
                                'Flu' = '#C34129')) +
                                # 'COVID-19' = "#EFB75B")) +
  labs(x = '', y = '') +
  scale_x_date(date_labels = "%b\n%y",
               date_breaks = "2 months") +
  envalysis::theme_publish() +
  theme(axis.text.x = element_text(size = 5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        legend.position = 'right')
ggsave("figures/23_24_peak_timing.png", height = 6, width = 12, units = "in", bg = "white")

## Determining in how many regions RSV peak occurs before Flu peak
find_states_with_early_rsv_peak = function(data) {
  data %>%
    group_by(state) %>%
    summarize(
      max_rsv_date = week_end[which.max(rsv_rescaled)],
      max_flu_date = week_end[which.max(flu_rescaled)],
      .groups = "drop"
    ) %>%
    filter(max_rsv_date < max_flu_date) %>%
    select(state, max_rsv_date, max_flu_date)
}

# Function calls
find_states_with_early_rsv_peak(signals_with_betas_16_17) #3/6
find_states_with_early_rsv_peak(signals_with_betas_17_18) #1/6
find_states_with_early_rsv_peak(signals_with_betas_18_19) #7/10
find_states_with_early_rsv_peak(signals_with_betas_19_20) #5/10
find_states_with_early_rsv_peak(signals_with_betas_20_21) # COVID year - exclude
find_states_with_early_rsv_peak(signals_with_betas_21_22) #10/10
find_states_with_early_rsv_peak(signals_with_betas_22_23) #4/10
find_states_with_early_rsv_peak(signals_with_betas_23_24) #2/10

# In 32/52 state-seasons RSV was first
# In the Pre-COVID times: 16/32
# In post/during COVID: 16/30

# Make the big combined figure ----
combined_signals = bind_rows(signals_with_betas_16_17,
                             signals_with_betas_17_18,
                             signals_with_betas_18_19,
                             signals_with_betas_19_20,
                             signals_with_betas_21_22,
                             signals_with_betas_22_23,
                             signals_with_betas_23_24)

combined_signals_ = combined_signals %>%
  group_by(season) %>%
  mutate(start_date = min(week_end), end_date = max(week_end)) %>%
  ungroup()

combined_signals_complete = combined_signals_ %>%
  complete(season = unique(season), state = unique(state), fill = list(
    flu_times_coef = NA,
    rsv_times_coef = NA,
    covid_times_coef = NA,
    ili_rescaled = NA
  ))

# Facet wrap version
combined_signals_complete %>% distinct(season)

combined_signals_%>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = flu_times_coef, color = 'Flu'), linewidth = 0.65, na.rm = TRUE) +
  geom_line(aes(y = rsv_times_coef, color = 'RSV'), linewidth = 0.65, na.rm = TRUE) +
  geom_line(aes(y = ili_rescaled, color = 'ILI'), linewidth = 0.65, na.rm = TRUE) +
  facet_grid2(~ state, scales = "free_x", independent = "x") +
  scale_color_manual(values = c('ILI' = 'grey',
                                'RSV' = 'dodgerblue4',
                                'Flu' = 'darkred')) +
  labs(x = '', y = '') +
  envalysis::theme_publish() +
  theme(
    axis.text.x = element_text(size = 5),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    legend.title = element_blank(),
    legend.position = 'right',
    strip.text = element_text(size = 10, face = "bold")  # State labels formatting
  )
# Facet grid version






combined_signals_ %>%
  mutate(season = factor(season, levels = unique(season)),  # Ensure seasons are ordered
         state = factor(state, levels = sort(unique(state)))) %>% # Ensure states are ordered
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = flu_times_coef, color = 'Flu'), linewidth = 0.65) +
  geom_line(aes(y = rsv_times_coef, color = 'RSV'), linewidth = 0.65) +
  geom_line(aes(y = covid_times_coef, color = "COVID-19"), linewidth = 0.65) +
  geom_line(aes(y = ili_rescaled, color = 'ILI'), linewidth = 0.65) +
  facet_grid2(season ~ state, scales = "free_x", independent = "x") +
  scale_x_date(date_labels = "%b",
                date_breaks = "2 months") +
  scale_color_manual(values = c('ILI' = '#E2DBC9',
                                'RSV' = '#183A5A',
                                'Flu' = '#C34129',
                                'COVID-19' = "#EFB75B")) +
  labs(x = '', y = '') +
  guides(color = guide_legend(override.aes = list(size = 5, linewidth = 5))) +
  envalysis::theme_publish() +
  theme(
    strip.text = element_text(size = 20, face = "bold"),  # Format state (column) labels
    axis.text.x = element_text(size = 5),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    legend.title = element_blank(),
    legend.position = 'bottom',
    legend.direction = 'horizontal',
    legend.text = element_text(size = 30),
    legend.spacing.x = unit(5, "in"),
    legend.spacing.y = unit(5, "in")
  )

ggsave(file = "figures/fig1_v1.png", width = 20, height = 15, units = "in", bg = "white")


ranked_signals = combined_signals_ %>%
  pivot_longer(cols = all_of(ends_with("coef")), names_to = "variable", values_to = "value") %>%
  group_by(state, season, variable) %>%
  slice_max(value, n = 1) %>%
  select(state, season, variable, week_end, value) %>%
  drop_na()

ranked_signals %>%
  group_by(state, season) %>%
  slice_max(value, n = 1) %>%
  select(state, season, variable, value) %>%
  view()

ranked_signals %>%
  group_by(state, season) %>%
  slice_min(week_end, n = 1) %>%
  select(state, season, variable, week_end)