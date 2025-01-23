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
## Last Updated: 2025-01-20
##
## Version: 1.4
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

# 2016-2017 season
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
signals_with_betas_16_17  %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = flu_times_coef, color = 'Flu'), linewidth = 0.65) +
  geom_line(aes(y = rsv_times_coef, color = 'RSV'), linewidth = 0.65) +
  geom_line(aes(y = ili_rescaled, color = 'ILI'), linewidth = 0.65) +
  facet_wrap(~state, scales = 'free', ncol = 3) +
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

# 2017-2018 season
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
  facet_wrap(~state, scales = 'free', ncol = 3) +
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

# 2018-2019 season
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

# 2019-2020 season

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

# 2020-2021 season
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

# 21-22 season
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

signals_with_betas_23_24


find_states_with_early_rsv_peak <- function(data) {
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

find_states_with_early_rsv_peak(signals_with_betas_16_17) #3/6
find_states_with_early_rsv_peak(signals_with_betas_17_18) #1/6
find_states_with_early_rsv_peak(signals_with_betas_18_19) #7/10
find_states_with_early_rsv_peak(signals_with_betas_19_20) #5/10
find_states_with_early_rsv_peak(signals_with_betas_20_21) # COVID year
find_states_with_early_rsv_peak(signals_with_betas_21_22) #10/10
find_states_with_early_rsv_peak(signals_with_betas_22_23) #4/10
find_states_with_early_rsv_peak(signals_with_betas_23_24) #2/10

# In 32/52 state-seasons RSV was first

# Example usage
`Ωresult <- find_states_with_early_rsv_peak(data)
print(result)

## Flu Data
## FluSurv
fluSurv_clean = fluSurv %>%
  filter(age_category == "Overall",
         race_category == "Overall",
         sex_category == "Overall") %>%
  mutate(week_end = as_date(MMWRweek2Date(mmwr_year, mmwr_week)) + 6) %>%
  select(week_end, catchment, weekly_rate, week_end) %>%
  drop_na()

# Filter the flu hospitalization data to match the start of the RSV data
# flu_hosps = flu_hosp_data %>% filter(date >= as_date("2022-10-01")) %>%
#   select(date, location_name, true_hosp) %>%
#   rename(week_end = date,
#          geography = location_name)
## Note here.. perhaps the flu hosps data is the cause of the problem? Circular
## logic with ILI data since the imputation relies on it

# Filter the COVID data to only the date ranges we are interested in
covid_data = covid_data %>%
  mutate(week_end = as_date(MMWRweek2Date(year, epiweek)) + 6) %>%
  filter(geography != "United States", week_end >= as_date("2021-10-03"))

## Get the list of states for which we have data - based on FluSurvNet (14) -
## exclude New York because we only have data for 2 counties
fluSurv_states = fluSurv_clean %>% filter(str_detect(catchment, "York") == F) %>%
  select(catchment) %>% unique() %>% pull()

states_data_available = nssp_data %>%
  filter(geography != "United States", county != "All") %>%
  select(geography) %>%
  distinct() %>% pull()

states_data_available =
  states_data_available[states_data_available != "District of Columbia"]

## 3. Merge the 4 datasets together --------------------------------------------
## 1. NSSP (RSV) data - % ED visits
## 2. ILInet data - % healthcare visits attributed to ILI
## 3. Flu hospitalization data - count of hospitalizations attributed to flu
## 4. COVID-19 confirmed cases

nssp_data_pv = nssp_data %>%
  select(week_end, geography, county, percent_visits_rsv,
         percent_visits_influenza, percent_visits_covid) %>%
  filter(county == "All")

## We want to create one dataset for each year:

## 3.1. 21/22 season; only flu, covid, and ILI --------------------------------------
flu_21_22 = fluSurv_clean %>%
  filter(catchment %in% fluSurv_states) %>%
  filter(week_end >= as_date("2021-10-09") & week_end < as_date("2022-10-01")) %>%
  rename(geography = catchment,
         flu_hosp_rate = weekly_rate)
## Note, NC does not have data before the 23/24 season

covid_21_22 = covid_data %>%
  filter(week_end >= as_date("2021-10-09") & week_end < as_date("2022-10-01")) %>%
  select(week_end, geography, cases) %>%
  rename(covid_cases = cases)

range(covid_21_22$week_end)

ili_21_22 = ili_data_with_dates %>%
  filter(region %in% fluSurv_states, week_end >= as_date("2021-10-09") &
           week_end < as_date("2022-10-01")) %>%
  mutate(week_end = week_end - 1) %>%
  select(week_end, region, unweighted_ili) %>%
  rename(geography = region)

range(ili_21_22$week_end)

flu_covid_ili_21_22 = flu_21_22 %>%
  left_join(covid_21_22, by = c("week_end", "geography")) %>%
  left_join(ili_21_22, by = c("week_end", "geography")) %>%
  rename(flu = flu_hosp_rate,
         covid = covid_cases,
         ili = unweighted_ili)

## 21/22 season: only flu, covid, and ILI
flu_21_22 = fluSurv_clean %>%
  filter(catchment %in% fluSurv_states) %>%
  filter(week_end >= as_date("2021-10-09") & week_end < as_date("2022-10-01")) %>%
  rename(geography = catchment,
         flu_hosp_rate = weekly_rate)
## Note, NC does not have data before the 23/24 season
range(flu_21_22$week_end)

covid_21_22 = covid_data %>%
  filter(week_end >= as_date("2021-10-09") & week_end < as_date("2022-10-01")) %>%
  select(week_end, geography, cases) %>%
  rename(covid_cases = cases)

range(covid_21_22$week_end)

ili_21_22 = ili_data_with_dates %>%
  filter(region %in% fluSurv_states, week_end >= as_date("2021-10-09") &
           week_end < as_date("2022-10-01")) %>%
  mutate(week_end = week_end - 1) %>%
  select(week_end, region, unweighted_ili) %>%
  rename(geography = region)

range(ili_21_22$week_end)

flu_covid_ili_21_22 = flu_21_22 %>%
  left_join(covid_21_22, by = c("week_end", "geography")) %>%
  left_join(ili_21_22, by = c("week_end", "geography")) %>%
  rename(flu = flu_hosp_rate,
         covid = covid_cases,
         ili = unweighted_ili)
range(flu_21_22$week_end)

write_csv(flu_covid_ili_21_22, file = "data/ground_truth/signals_21_22.csv")

signals_21_22_with_betas = NULL
for(state in fluSurv_states[-14]){

  tmp_model = lm(ili ~ flu + covid - 1,
                 data = flu_covid_ili_21_22 %>%
                   filter(geography == state))

  # predicted_vals = predict(tmp_model, newdata = signals_22_23 %>%
  #                            filter(geography == state) %>%
  #                            select(rsv, flu, covid) %>%
  #                            as.data.frame)

  coefs_tmp = summary(tmp_model)$coef

  signals_tmp = flu_covid_ili_21_22 %>%
    filter(geography == state) %>%
    mutate(beta_flu = coefs_tmp[1],
           beta_covid = coefs_tmp[2],
           flu_times_coef = flu*beta_flu,
           cov_times_coef = covid*beta_covid)

  signals_21_22_with_betas = bind_rows(signals_21_22_with_betas, signals_tmp)
}

## Figure 3.1: Generate the facet plot for the 21/22 season - ili ~ flu + COVID -----
signals_21_22_with_betas %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = flu_times_coef, color = "Flu")) +
  geom_line(aes(y = cov_times_coef, color = "COVID")) +
  geom_line(aes(y = ili, color = "ILI")) +
  facet_wrap(~geography, scales = "free", ncol = 5) +
  scale_color_manual(values = c("ILI" = "grey",
                                "RSV" = "dodgerblue4",
                                "Flu" = "darkred",
                                "COVID" = "gold2")) +
  theme_minimal() +
  scale_x_date(date_labels = "%b\n%y",
               limits = c(as_date("2021-10-01"), "2022-6-11"),
               date_breaks = "2 months") +
  theme_minimal() +
  labs(x = "", y = "") +
  theme(axis.text.x = element_text(size = 5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "right")
# ggsave("data/figures/ili_flu_covid_21_22.png",
#        width = 10, height = 6, unit = "in", bg = "white")


## 3.2. 22/23 season: flu, covid, RSV, ILI (using NSSP flu data) --------------------
flu_22_23 = nssp_data_pv %>%
  filter(week_end >= as_date("2022-06-01") & week_end < as_date("2023-03-25"),
         geography %in% states_data_available) %>%
  select(week_end, geography, percent_visits_influenza) %>%
  rename(flu = percent_visits_influenza)

range(flu_22_23$week_end)

rsv_22_23 = nssp_data_pv %>%
  filter(week_end >= as_date("2022-10-01") & week_end < as_date("2023-03-25"),
         geography %in% states_data_available) %>%
  select(week_end, geography, percent_visits_rsv) %>%
  rename(rsv = percent_visits_rsv)

range(rsv_22_23$week_end)

covid_22_23 = covid_data %>%
  filter(week_end >= as_date("2022-10-01") & week_end < as_date("2023-03-25")) %>%
  select(week_end, geography, cases) %>%
  rename(covid_cases = cases)

range(covid_22_23$week_end)

ili_22_23 = ili_data_with_dates %>%
  filter(region %in% states_data_available, week_end >= as_date("2022-10-01") &
           week_end < as_date("2023-03-25")) %>%
  mutate(week_end = week_end - 1) %>%
  select(week_end, region, unweighted_ili) %>%
  rename(geography = region)

range(ili_22_23$week_end)

flu_rsv_covid_ili_22_23 = flu_22_23 %>%
  left_join(covid_22_23, by = c("week_end", "geography")) %>%
  left_join(ili_22_23, by = c("week_end", "geography")) %>%
  left_join(rsv_22_23, by = c("week_end", "geography")) %>%
  rename(covid = covid_cases,
         ili = unweighted_ili)

write_csv(flu_rsv_covid_ili_22_23, file = "data/ground_truth/signals_22_23.csv")

signals_22_23_with_betas = NULL
for(state in states_data_available){

  flu_covid_ili_21_22

  tmp_model = lm(ili ~ flu + covid + rsv - 1,
                 data = flu_rsv_covid_ili_22_23 %>%
                   filter(geography == state))

  # predicted_vals = predict(tmp_model, newdata = signals_22_23 %>%
  #                            filter(geography == state) %>%
  #                            select(rsv, flu, covid) %>%
  #                            as.data.frame)

  coefs_tmp = summary(tmp_model)$coef

  signals_tmp = flu_rsv_covid_ili_22_23 %>%
    filter(geography == state) %>%
    mutate(beta_flu = coefs_tmp[1],
           beta_covid = coefs_tmp[2],
           beta_rsv = coefs_tmp[3],
           flu_times_coef = flu*beta_flu,
           cov_times_coef = covid*beta_covid,
           rsv_times_coef = rsv*beta_rsv)

  signals_22_23_with_betas = bind_rows(signals_22_23_with_betas, signals_tmp)
}

## Figure 3.2. Facet plot for 22/23 season - ili ~ flu + COVID + rsv ----------
signals_22_23_with_betas %>%
  filter(geography %in% states_data_available[1:24]) %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = flu_times_coef, color = "Flu")) +
  geom_line(aes(y = cov_times_coef, color = "COVID")) +
  geom_line(aes(y = rsv_times_coef, color = "RSV")) +
  geom_line(aes(y = ili, color = "ILI")) +
  facet_wrap(~geography, scales = "free", ncol = 6) +
  scale_color_manual(values = c("ILI" = "grey",
                                "RSV" = "dodgerblue4",
                                "Flu" = "darkred",
                                "COVID" = "gold2")) +
  theme_minimal() +
  scale_x_date(date_labels = "%b\n%y",
               date_breaks = "2 months") +
  labs(x = "", y = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 5),
        legend.title = element_blank(),
        legend.position = "right")

# ggsave("data/figures/ili_flu_covid_rsv_22_23_A.png",
#        width = 10, height = 6, unit = "in", bg = "white")


signals_22_23_with_betas %>%
  filter(geography %in% states_data_available[25:47]) %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = flu_times_coef, color = "Flu")) +
  geom_line(aes(y = cov_times_coef, color = "COVID")) +
  geom_line(aes(y = rsv_times_coef, color = "RSV")) +
  geom_line(aes(y = ili, color = "ILI")) +
  facet_wrap(~geography, scales = "free", ncol = 6) +
  scale_color_manual(values = c("ILI" = "grey",
                                "RSV" = "dodgerblue4",
                                "Flu" = "darkred",
                                "COVID" = "gold2")) +
  theme_minimal() +
  scale_x_date(date_labels = "%b\n%y",
               date_breaks = "2 months") +
  theme_minimal() +
  labs(x = "", y = "") +
  theme(axis.text.x = element_text(size = 5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 5),
        legend.title = element_blank(),
        legend.position = "right")

# ggsave("data/figures/ili_flu_covid_rsv_22_23_B.png",
#        width = 10, height = 6, unit = "in", bg = "white")

## 3.3. 23/24 season: flu, RSV, ILI (using NSSP flu data) ----------------------
flu_23_24 = nssp_data_pv %>%
  filter(week_end >= as_date("2023-04-01") & week_end < as_date("2024-06-29"),
         geography %in% states_data_available) %>%
  select(week_end, geography, percent_visits_influenza) %>%
  rename(flu = percent_visits_influenza)

range(flu_23_24$week_end)

rsv_23_24 = nssp_data_pv %>%
  filter(week_end >= as_date("2023-04-01") & week_end < as_date("2024-06-29"),
         geography %in% states_data_available) %>%
  select(week_end, geography, percent_visits_rsv) %>%
  rename(rsv = percent_visits_rsv)

range(rsv_23_24$week_end)

ili_23_24 = ili_data_with_dates %>%
  filter(region %in% states_data_available, week_end >= as_date("2023-04-01") &
           week_end < as_date("2024-06-29")) %>%
  mutate(week_end = week_end - 1) %>%
  select(week_end, region, unweighted_ili) %>%
  rename(geography = region)

range(ili_23_24$week_end)

flu_rsv_ili_23_24 = flu_23_24 %>%
  left_join(ili_23_24, by = c("week_end", "geography")) %>%
  left_join(rsv_23_24, by = c("week_end", "geography")) %>%
  rename(ili = unweighted_ili)

write_csv(flu_rsv_ili_23_24, file = "signals_23_24.csv")

signals_23_24_with_betas = NULL
for(state in states_data_available){

  tmp_model = lm(ili ~ flu + rsv - 1,
                 data = flu_rsv_ili_23_24 %>%
                   filter(geography == state))

  coefs_tmp = summary(tmp_model)$coef

  signals_tmp = flu_rsv_ili_23_24 %>%
    filter(geography == state) %>%
    mutate(beta_flu = coefs_tmp[1],
           beta_rsv = coefs_tmp[2],
           flu_times_coef = flu*beta_flu,
           rsv_times_coef = rsv*beta_rsv)

  signals_23_24_with_betas = bind_rows(signals_23_24_with_betas, signals_tmp)
}

## Figure 3.3. Facet plot for 23/24 season - ili ~ flu + COVID + rsv -----------
signals_23_24_with_betas %>%
  filter(geography %in% states_data_available[1:24]) %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = flu_times_coef, color = "Flu")) +
  geom_line(aes(y = rsv_times_coef, color = "RSV")) +
  geom_line(aes(y = ili, color = "ILI")) +
  facet_wrap(~geography, scales = "free", ncol = 6) +
  scale_color_manual(values = c("ILI" = "grey",
                                "RSV" = "dodgerblue4",
                                "Flu" = "darkred",
                                "COVID" = "gold2")) +
  theme_minimal() +
  scale_x_date(date_labels = "%b\n%y",
               date_breaks = "2 months") +
  labs(x = "", y = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 5),
        legend.title = element_blank(),
        legend.position = "right")

# ggsave("data/figures/ili_flu_rsv_23_24_A.png",
#        width = 10, height = 6, unit = "in", bg = "white")


signals_23_24_with_betas %>%
  filter(geography %in% states_data_available[25:47]) %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = flu_times_coef, color = "Flu")) +
  geom_line(aes(y = rsv_times_coef, color = "RSV")) +
  geom_line(aes(y = ili, color = "ILI")) +
  facet_wrap(~geography, scales = "free", ncol = 6) +
  scale_color_manual(values = c("ILI" = "grey",
                                "RSV" = "dodgerblue4",
                                "Flu" = "darkred",
                                "COVID" = "gold2")) +
  theme_minimal() +
  scale_x_date(date_labels = "%b\n%y",
               date_breaks = "2 months") +
  theme_minimal() +
  labs(x = "", y = "") +
  theme(axis.text.x = element_text(size = 5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 5),
        legend.title = element_blank(),
        legend.position = "right")

# ggsave("data/figures/ili_flu_rsv_23_24_B.png",
#        width = 10, height = 6, unit = "in", bg = "white")

## OLD STUFF BELOW -------------------------------------------------------------


















## Generate the facet plot for the 21/22 season - ili ~ flu + COVID
signals_21_22_with_betas %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = flu_times_coef, color = "Flu")) +
  geom_line(aes(y = cov_times_coef, color = "COVID")) +
  geom_line(aes(y = ili, color = "ILI")) +
  facet_wrap(~geography, scales = "free", ncol = 5) +
  scale_color_manual(values = c("ILI" = "grey",
                                "RSV" = "dodgerblue4",
                                "Flu" = "darkred",
                                "COVID" = "gold2")) +
  theme_minimal() +
  scale_x_date(date_labels = "%b\n%y",
               limits = c(as_date("2021-10-01"), "2022-6-11"),
               date_breaks = "2 months") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 5),
        legend.title = element_blank(),
        legend.position = "right")
# ggsave("data/figures/ili_flu_covid_21_22.png",
#        width = 10, height = 6, unit = "in", bg = "white")

## Restrict the data to only include 22/23 season where we have all data
signals_22_23 = nssp_ILI_flu_COVID %>%
  filter(week_end < as_date("2023-04-01")) %>%
  drop_na()

signals_statewise_w_betas = NULL
for(state in high_pop_states){

  tmp_model = lm(ili ~ rsv + flu + covid - 1,
                 data = signals_22_23 %>% filter(geography == state))

  predicted_vals = predict(tmp_model, newdata = signals_22_23 %>%
                             filter(geography == state) %>%
                             select(rsv, flu, covid) %>%
                             as.data.frame)

  coefs_tmp = summary(tmp_model)$coef

  signals_tmp = signals_22_23 %>%
    filter(geography == state) %>%
    mutate(beta_rsv = coefs_tmp[1],
           beta_flu = coefs_tmp[2],
           beta_covid = coefs_tmp[3],
           rsv_times_coef = rsv*beta_rsv,
           flu_times_coef = flu*beta_flu,
           cov_times_coef = covid*beta_covid)

  signals_statewise_w_betas = bind_rows(signals_statewise_w_betas, signals_tmp)
}

signals_statewise_w_betas %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = rsv_times_coef, color = "RSV")) +
  geom_line(aes(y = flu_times_coef, color = "Flu")) +
  geom_line(aes(y = cov_times_coef, color = "COVID")) +
  geom_line(aes(y = ili, color = "ILI")) +
  facet_wrap(~factor(geography, levels = c(c(state_pops$state[1:24])))) +
  scale_color_manual(values = c("ILI" = "grey",
                                "RSV" = "dodgerblue4",
                                "Flu" = "darkred",
                                "COVID" = "gold2"))


























ili_data_RSV_short = ili_data_RSV %>%
  filter(region %in% states_data_available) %>%
  select(week_end, region, unweighted_ili)

nssp_ILI_flu_COVID = nssp_data %>%
  filter(geography != "United States", county == "All") %>%
  left_join(ili_data_RSV_short, by = c("week_end", "geography" = "region")) %>%
  left_join(flu_hosps, by = c("week_end", "geography")) %>%
  left_join(covid_data, by = c("week_end", "geography"))

## Rename variables and fix format
nssp_ILI_flu_COVID = nssp_ILI_flu_COVID %>%
  rename(rsv = "percent_visits_rsv",
         ili = "unweighted_ili",
         flu = true_hosp,
         covid = cases) %>%
  mutate(flu = as.numeric(flu)) %>%
  select(week_end, geography, rsv, ili, flu, covid)

## 4. Rescale each time series -------------------------------------------------
## 1. Use the max value of each trend (so modify it to the %max of each trend)

data_combined_rescaled = nssp_ILI_flu_COVID %>%
  filter(geography != "District of Columbia") %>%
  group_by(geography) %>%
  mutate(rsv_rescaled = rsv/max(rsv),
         ili_rescaled = ili/max(ili),
         flu_rescaled = flu/max(flu, na.rm = T),
         cov_rescaled = covid/max(covid, na.rm = T))

## Since we only want 20 plots per page, we split the data into two subsets

# Get state population data to rank the states and split
`%notin%` <- Negate(`%in%`)

state_pops = read_csv("~/Documents/Projects/rsv_prediction/code/state-population-table.csv")
state_pops = state_pops %>% select(state, population) %>%
  mutate(population_rank = dense_rank(desc(population))) %>%
  filter(state %notin% c("South Dakota",
                         "New Hampshire",
                         "Missouri"))

high_pop_states = state_pops$state[1:24]
low_pop_states = state_pops$state[25:nrow(state_pops)]

data_combined_rescaled_1 =
  data_combined_rescaled %>% filter(geography %in% high_pop_states)
data_combined_rescaled_2 =
  data_combined_rescaled %>% filter(geography %in% low_pop_states)

## 5. Plotting all signals together (rescaled) ---------------------------------

## 5.1 2022-23 flu season (with RSV/COVID superimposed) ------------------------
data_combined_rescaled_1 %>%
  mutate(geography = factor(geography, levels =
                              c(state_pops$state[1:24]))) %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = ili_rescaled, color = "ILI")) +
  geom_line(aes(y = rsv_rescaled, color = "RSV")) +
  geom_line(aes(y = flu_rescaled, color = "Flu")) +
  geom_line(aes(y = cov_rescaled, color = "COVID")) +
  facet_wrap(~geography, scales = "free", ncol = 6) +
  labs(x = "", y = "% of Max Signal") +
  scale_color_manual(values = c("ILI" = "grey",
                                "RSV" = "dodgerblue4",
                                "Flu" = "darkred",
                                "COVID" = "gold2")) +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1), expand = c(0,0)) +
  scale_x_date(date_labels = "%b\n%y",
               limits = c(as_date("2022-10-01"), as_date("2023-03-25"))) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 5),
        legend.title = element_blank(),
        legend.position = "right")

ggsave(file = "~/Google Drive/My Drive/EWS_GW_BP_LC_RG/data/figures/RSV_ILI_Flu_COVID/rescaled_high_pop_states.png", width = 10, height = 6, unit = "in", bg = "white")

data_combined_rescaled_2 %>%
  mutate(geography = factor(geography,
                            levels = c(state_pops$state[25:nrow(state_pops)]))) %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = ili_rescaled, color = "ILI")) +
  geom_line(aes(y = rsv_rescaled, color = "RSV")) +
  geom_line(aes(y = flu_rescaled, color = "Flu")) +
  geom_line(aes(y = cov_rescaled, color = "COVID")) +
  facet_wrap(~geography, scales = "free", ncol = 6) +
  labs(x = "", y = "% of Max Signal") +
  scale_color_manual(values = c("ILI" = "grey",
                                "RSV" = "dodgerblue4",
                                "Flu" = "darkred",
                                "COVID" = "gold2")) +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1), expand = c(0,0)) +
  scale_x_date(date_labels = "%b\n%y",
               limits = c(as_date("2022-10-01"), as_date("2023-03-25"))) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 5),
        legend.title = element_blank(),
        legend.position = "right")
ggsave(file = "~/Google Drive/My Drive/EWS_GW_BP_LC_RG/data/figures/RSV_ILI_Flu_COVID/rescaled_low_pop_states.png", width = 10, height = 6, unit = "in", bg = "white")

## 5.2 Plot both 2022-23 and 23-24 seasons (and exclude COVID) -----------------
data_combined_rescaled_1 %>%
  mutate(geography = factor(geography, levels =
                              c(state_pops$state[1:24]))) %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = ili_rescaled, color = "ILI")) +
  geom_line(aes(y = rsv_rescaled, color = "RSV")) +
  geom_line(aes(y = flu_rescaled, color = "Flu")) +
  facet_wrap(~geography, scales = "free", ncol = 6) +
  labs(x = "", y = "% of Max Signal") +
  scale_color_manual(values = c("ILI" = "grey",
                                "RSV" = "dodgerblue4",
                                "Flu" = "darkred",
                                "COVID" = "gold2")) +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1), expand = c(0,0)) +
  scale_x_date(date_labels = "%b\n%y",
               limits = c(as_date("2022-10-01"), as_date("2024-06-29"))) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 5),
        legend.title = element_blank(),
        legend.position = "right")
ggsave(file = "~/Google Drive/My Drive/EWS_GW_BP_LC_RG/data/figures/RSV_ILI_Flu_COVID/rescaled_high_pop_states_no_covid.png", width = 10, height = 6, unit = "in", bg = "white")

data_combined_rescaled_2 %>%
  mutate(geography = factor(geography, levels =
                              c(state_pops$state[25:nrow(state_pops)]))) %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = ili_rescaled, color = "ILI")) +
  geom_line(aes(y = rsv_rescaled, color = "RSV")) +
  geom_line(aes(y = flu_rescaled, color = "Flu")) +
  facet_wrap(~geography, scales = "free", ncol = 6) +
  labs(x = "", y = "% of Max Signal") +
  scale_color_manual(values = c("ILI" = "grey",
                                "RSV" = "dodgerblue4",
                                "Flu" = "darkred",
                                "COVID" = "gold2")) +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1), expand = c(0,0)) +
  scale_x_date(date_labels = "%b\n%y",
               limits = c(as_date("2022-10-01"), as_date("2024-06-29"))) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 5),
        legend.title = element_blank(),
        legend.position = "right")
ggsave(file = "~/Google Drive/My Drive/EWS_GW_BP_LC_RG/data/figures/RSV_ILI_Flu_COVID/rescaled_low_pop_states_no_covid.png", width = 10, height = 6, unit = "in", bg = "white")

## 6. Normalize the signals using bestNormalize --------------------------------

# States with complete cases only
# states_complete = c("California", "Texas", "Florida", "New York", "Pennsylvania", "Illinois", "Ohio", "Georgia", "North Carolina", "Michigan", "Virginia", "Washington", "Arizona", "Tennessee", "Massachusetts", "Indiana", "Maryland", "Wisconsin", "Colorado", "Minnesota", "South Carolina", "Alabama", "Louisiana", "Kentucky", "Oregon", "Oklahoma", "Connecticut", "Iowa", "Arkansas", "Mississippi", "New Mexico", "Idaho", "Nebraska", "West Virginia", "Hawaii", "Montana", "Vermont", "Wyoming")

bc_transformed_df = tibble()
set.seed(824)

for(state in states_data_available){

  tmp_dates = data_combined_rescaled %>% ungroup() %>%
    filter(geography == state, week_end > as_date("2022-10-01"),
           week_end < as_date("2023-04-01")) %>%
    select(week_end) %>% pull()
  tmp_rsv = data_combined_rescaled %>%  ungroup() %>%
    filter(geography == state, week_end > as_date("2022-10-01"),
           week_end < as_date("2023-04-01"))  %>%
    select(rsv) %>% pull()
  tmp_ili = data_combined_rescaled %>% ungroup() %>%
    filter(geography == state, week_end > as_date("2022-10-01"),
           week_end < as_date("2023-04-01"))  %>%
    select(ili) %>% pull()
  tmp_flu = data_combined_rescaled %>% ungroup() %>%
    filter(geography == state, week_end > as_date("2022-10-01"),
           week_end < as_date("2023-04-01"))  %>%
    select(flu) %>% pull()
  tmp_cov = data_combined_rescaled %>% ungroup() %>%
    filter(geography == state, week_end > as_date("2022-10-01"),
           week_end < as_date("2023-04-01"))  %>%
    select(covid) %>% pull()

  # Use bestNormalize to transform (chooses the best out of 8 methods)
  tmp_rsv_normalized = bestNormalize(tmp_rsv, out_of_sample = F)
  tmp_ili_normalized = bestNormalize(tmp_ili, out_of_sample = F)
  tmp_flu_normalized = bestNormalize(tmp_flu, out_of_sample = F)
  tmp_cov_normalized = bestNormalize(tmp_cov, out_of_sample = F)

  # Unlist and convert to vector
  tmp_rsv_normalized_vec =
    as.numeric(unlist(tmp_rsv_normalized)[1:length(tmp_rsv)])
  tmp_ili_normalized_vec =
    as.numeric(unlist(tmp_ili_normalized)[1:length(tmp_ili)])
  tmp_flu_normalized_vec =
    as.numeric(unlist(tmp_flu_normalized)[1:length(tmp_flu)])
  tmp_cov_normalized_vec =
    as.numeric(unlist(tmp_cov_normalized)[1:length(tmp_cov)])

  # Aggregate the data into a temporary tibble
  tmp_df = tibble(geography = state,
                  date = tmp_dates,
                  rsv_normalized = tmp_rsv_normalized_vec,
                  ili_normalized = tmp_ili_normalized_vec,
                  flu_normalized = tmp_flu_normalized_vec,
                  cov_normalized = tmp_cov_normalized_vec)

  # Add the tibble to the transformed data tibble
  bc_transformed_df = bind_rows(bc_transformed_df, tmp_df)

  print(paste0(state," data completed at ",Sys.Date()))
}
bc_transformed_df

## Get correlations - not needed at this point
# cor_df = NULL
# for(state in states_data_available){
#   tmp_data = nssp_ILI_RSV %>% filter(geography == state)
#   cor_row = c(state,cor(tmp_data$percent_visits_rsv,
#                         tmp_data$unweighted_ili),
#               ccf(tmp_data$percent_visits_rsv,
#                   tmp_data$unweighted_ili)[1]$acf)
#   cor_df = rbind(cor_df, cor_row)
# }
# cor_df = cor_df %>% as_tibble()
# names(cor_df) = c("geography", "correlation", "lag1_correlation")
# cor_df = cor_df %>% mutate(correlation = as.numeric(correlation, 3),
#                            lag1_correlation = as.numeric(lag1_correlation, 3))

## 6.1 Plotting all 4 signals after normalization ----------------------------------

data_combined_normalized_1 =
  bc_transformed_df %>% filter(geography %in% high_pop_states)
data_combined_normalized_2 =
  bc_transformed_df %>% filter(geography %in% low_pop_states)

bc_transformed_df

data_combined_normalized_1 %>%
  mutate(geography = factor(geography, levels =
                              c(state_pops$state[1:24]))) %>%
  ggplot(aes(x = date)) +
  geom_line(aes(y = ili_normalized, color = "ILI")) +
  geom_line(aes(y = rsv_normalized, color = "RSV")) +
  geom_line(aes(y = flu_normalized, color = "Flu")) +
  geom_line(aes(y = cov_normalized, color = "COVID")) +
  facet_wrap(~geography, scales = "free", ncol = 6) +
  labs(x = "", y = "Normalized Signal") +
  scale_color_manual(values = c("ILI" = "grey",
                                "RSV" = "dodgerblue4",
                                "Flu" = "darkred",
                                "COVID" = "gold2")) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_date(date_labels = "%b\n%y",
               limits = c(as_date("2022-10-01"), as_date("2023-03-25"))) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 5),
        legend.title = element_blank(),
        legend.position = "right")
ggsave(file = "~/Google Drive/My Drive/EWS_GW_BP_LC_RG/data/figures/RSV_ILI_Flu_COVID/normalized_high_pop_states.png", width = 10, height = 6, unit = "in", bg = "white")

data_combined_normalized_2 %>%
  mutate(geography = factor(geography, levels =
                              c(state_pops$state[25:nrow(state_pops)]))) %>%
  ggplot(aes(x = date)) +
  geom_line(aes(y = ili_normalized, color = "ILI")) +
  geom_line(aes(y = rsv_normalized, color = "RSV")) +
  geom_line(aes(y = flu_normalized, color = "Flu")) +
  geom_line(aes(y = cov_normalized, color = "COVID")) +
  facet_wrap(~geography, scales = "free", ncol = 6) +
  labs(x = "", y = "Normalized Signal") +
  scale_color_manual(values = c("ILI" = "grey",
                                "RSV" = "dodgerblue4",
                                "Flu" = "darkred",
                                "COVID" = "gold2")) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_date(date_labels = "%b\n%y",
               limits = c(as_date("2022-10-01"), as_date("2023-03-25"))) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 5),
        legend.title = element_blank(),
        legend.position = "right")
ggsave(file = "~/Google Drive/My Drive/EWS_GW_BP_LC_RG/data/figures/RSV_ILI_Flu_COVID/normalized_low_pop_states.png", width = 10, height = 6, unit = "in", bg = "white")

## 6.2 Redo the normalization without COVID and plot for the entire data span --

bc_transformed_df_2 = tibble()
set.seed(824)

for(state in states_data_available){

  tmp_dates = data_combined_rescaled %>% ungroup() %>%
    filter(geography == state, week_end > as_date("2022-10-01")) %>%
    select(week_end) %>% pull()
  tmp_rsv = data_combined_rescaled %>%  ungroup() %>%
    filter(geography == state, week_end > as_date("2022-10-01"))  %>%
    select(rsv) %>% pull()
  tmp_ili = data_combined_rescaled %>% ungroup() %>%
    filter(geography == state, week_end > as_date("2022-10-01"))  %>%
    select(ili) %>% pull()
  tmp_flu = data_combined_rescaled %>% ungroup() %>%
    filter(geography == state, week_end > as_date("2022-10-01"))  %>%
    select(flu) %>% pull()

  # Use bestNormalize to transform (chooses the best out of 8 methods)
  tmp_rsv_normalized = bestNormalize(tmp_rsv, out_of_sample = F)
  tmp_ili_normalized = bestNormalize(tmp_ili, out_of_sample = F)
  tmp_flu_normalized = bestNormalize(tmp_flu, out_of_sample = F)

  # Unlist and convert to vector
  tmp_rsv_normalized_vec = as.numeric(unlist(tmp_rsv_normalized)[1:length(tmp_rsv)])
  tmp_ili_normalized_vec = as.numeric(unlist(tmp_ili_normalized)[1:length(tmp_ili)])
  tmp_flu_normalized_vec = as.numeric(unlist(tmp_flu_normalized)[1:length(tmp_flu)])

  # Aggregate the data into a temporary tibble
  tmp_df = tibble(geography = state,
                  date = tmp_dates,
                  rsv_normalized = tmp_rsv_normalized_vec,
                  ili_normalized = tmp_ili_normalized_vec,
                  flu_normalized = tmp_flu_normalized_vec)

  # Add the tibble to the transformed data tibble
  bc_transformed_df_2 = bind_rows(bc_transformed_df_2, tmp_df)

  print(paste0(state," data completed at ",Sys.Date()))
}
bc_transformed_df_2

data_combined_normalized_nc_1 =
  bc_transformed_df_2 %>% filter(geography %in% high_pop_states)
data_combined_normalized_nc_2 =
  bc_transformed_df_2 %>% filter(geography %in% low_pop_states)

bc_transformed_df

data_combined_normalized_nc_1 %>%
  mutate(geography = factor(geography, levels =
                              c(state_pops$state[1:24]))) %>%
  ggplot(aes(x = date)) +
  geom_line(aes(y = ili_normalized, color = "ILI")) +
  geom_line(aes(y = rsv_normalized, color = "RSV")) +
  geom_line(aes(y = flu_normalized, color = "Flu")) +
  facet_wrap(~geography, scales = "free", ncol = 6) +
  labs(x = "", y = "Normalized Signal") +
  scale_color_manual(values = c("ILI" = "grey",
                                "RSV" = "dodgerblue4",
                                "Flu" = "darkred",
                                "COVID" = "gold2")) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_date(date_labels = "%b\n%y", limits = c(as_date("2022-10-01"),
                                                  as_date("2024-07-01"))) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 5),
        legend.title = element_blank(),
        legend.position = "right")

data_combined_normalized_nc_2 %>%
  mutate(geography = factor(geography, levels =
                              c(state_pops$state[25:nrow(state_pops)]))) %>%
  ggplot(aes(x = date)) +
  geom_line(aes(y = ili_normalized, color = "ILI")) +
  geom_line(aes(y = rsv_normalized, color = "RSV")) +
  geom_line(aes(y = flu_normalized, color = "Flu")) +
  facet_wrap(~geography, scales = "free", ncol = 6) +
  labs(x = "", y = "Normalized Signal") +
  scale_color_manual(values = c("ILI" = "grey",
                                "RSV" = "dodgerblue4",
                                "Flu" = "darkred",
                                "COVID" = "gold2")) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_date(date_labels = "%b\n%y", limits = c(as_date("2022-10-01"),
                                                  as_date("2024-07-01"))) +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 5),
        legend.title = element_blank(),
        legend.position = "right")

## 7. Plotting the true data based on stacked regression curves ----------------

summary(lm(ili ~ rsv + flu + covid - 1, data = nssp_ILI_flu_COVID))


nssp_ILI_flu_COVID %>%
  mutate(rsv_times_coef = rsv*2.8250195,
         flu_times_coef = flu*0.0060596,
         covid_times_coef = covid*0.000022942) %>%
  filter(geography == "California") %>%
  ggplot(aes(x = rsv_times_coef + flu_times_coef + covid_times_coef, y = ili)) +
  geom_point()

## Restrict the data to only include 22/23 season where we have all data
signals_22_23 = nssp_ILI_flu_COVID %>%
  filter(week_end < as_date("2023-04-01")) %>%
  drop_na()

nnls(x = as.matrix(signals_22_23[,c("rsv", "ili")]), y = as.matrix(signals_22_23[,"ili"]))

## Figure out where the issue is.. is it a state-by-state problem?

high_pop_states

signals_22_23

signals_statewise_w_betas = NULL
for(state in high_pop_states){

  tmp_model = lm(ili ~ rsv + flu + covid - 1,
                 data = signals_22_23 %>% filter(geography == state))

  predicted_vals = predict(tmp_model, newdata = signals_22_23 %>%
            filter(geography == state) %>%
            select(rsv, flu, covid) %>%
            as.data.frame)

  coefs_tmp = summary(tmp_model)$coef

  signals_tmp = signals_22_23 %>%
    filter(geography == state) %>%
    mutate(beta_rsv = coefs_tmp[1],
           beta_flu = coefs_tmp[2],
           beta_covid = coefs_tmp[3],
           rsv_times_coef = rsv*beta_rsv,
           flu_times_coef = flu*beta_flu,
           cov_times_coef = covid*beta_covid)

  signals_statewise_w_betas = bind_rows(signals_statewise_w_betas, signals_tmp)
}

signals_statewise_w_betas %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = rsv_times_coef, color = "RSV")) +
  geom_line(aes(y = flu_times_coef, color = "Flu")) +
  geom_line(aes(y = cov_times_coef, color = "COVID")) +
  geom_line(aes(y = ili, color = "ILI")) +
  facet_wrap(~factor(geography, levels = c(c(state_pops$state[1:24])))) +
  scale_color_manual(values = c("ILI" = "grey",
                                "RSV" = "dodgerblue4",
                                "Flu" = "darkred",
                                "COVID" = "gold2"))

signals_tmp = signals_22_23 %>% filter(geography == "Massachusetts")

signals_tmp %>%
  ggplot(aes(x = week_end))+
  geom_line(aes(y = rsv_times_coef, color = "RSV")) +
  geom_line(aes(y = flu_times_coef, color = "Flu")) +
  geom_line(aes(y = cov_times_coef, color = "COVID")) +
  geom_line(aes(y = ili, color = "ILI")) +
  scale_color_manual(values = c("ILI" = "grey",
                                "RSV" = "dodgerblue4",
                                "Flu" = "darkred",
                                "COVID" = "gold2"))

coefs = summary(lm(ili ~ rsv + flu + covid + 0, data = signals_tmp))$coef
coefs

signals_tmp %>%
  mutate(rsv_times_coef = rsv*coefs[1],
         flu_times_coef = flu*coefs[2],
         covid_times_coef = covid*coefs[3]) %>%
  ggplot(aes(x = week_end))+
  geom_line(aes(y= rsv_times_coef, color = "RSV")) +
  geom_line(aes(y= flu_times_coef, color = "Flu")) +
  geom_line(aes(y = covid_times_coef, color = "COVID")) +
  geom_line(aes(y = ili, color = "ILI")) +
  scale_color_manual(values = c("ILI" = "grey",
                                "RSV" = "dodgerblue4",
                                "Flu" = "darkred",
                                "COVID" = "gold2"))

signals_AL %>%
  mutate(rsv_times_coef = rsv*coefs[1],
         flu_times_coef = flu*coefs[2],
         covid_times_coef = covid*coefs[3]) %>%
  ggplot(aes(x = rsv_times_coef + flu_times_coef + covid_times_coef, y = ili)) +
  geom_point() +
  scale_x_continuous(limits = c(0, 10))



signals %>%
  mutate(rsv_times_coef = rsv*3.29136,
         flu_times_coef = flu*-0.0007447,
         covid_times_coef = covid*0.00006518) %>%
  filter(geography == "California") %>%
  ggplot(aes(x = week_end))+
  geom_line(aes(y= rsv_times_coef, color = "RSV")) +
  geom_line(aes(y= flu_times_coef, color = "FLU")) +
  geom_line(aes(y = covid_times_coef, color = "COVID")) +
  geom_line(aes(y = ili, color = "ILI")) +
  scale_color_manual(values = c("ILI" = "grey",
                                "RSV" = "dodgerblue4",
                                "FLU" = "darkred",
                                "COVID" = "gold2"))



# Function to compute the beta coefficient
compute_betas = function(window_data) {
  model = lm(ili ~ rsv + flu + 0, data = window_data)
  return(coef(model)) # Return the beta coefficients
}

nssp_ILI_flu = nssp_ILI_flu_COVID %>%
  select(-covid)

state = "California"

test_df = nssp_ILI_flu %>% filter(geography == state) %>% as.data.frame()

# Perform rolling regression
roll_lm_coefs = roll_lm(x = as.matrix(test_df[,c(3, 5)]),
        y = as.matrix(test_df[,4]),
        width = 3, intercept = F)$coefficients
roll_lm_coefs
test_d
tmp = cbind(test_df$week_end, test_df$ili, roll_lm_coefs) %>% as.data.frame()

# Function to compute the beta coefficients using NNLS
test_df = test_df %>% drop_na()
compute_betas <- function(window_data) {
  X <- as.matrix(window_data[, c("rsv", "flu")])
  Y <- as.matrix(window_data[, "ili"])
  model <- nnls(X, Y)
  # model = glmnet(X, Y, lower.limits = 0, intercept = FALSE)
  return(coef(model)) # Return the beta coefficients
}

window_size = 3

test_df2 = test_df[, c("ili", "rsv", "flu")] %>%
  as.matrix()

test_df2

lm()

# Perform rolling regression
betas <- rollapply(test_df2, width = window_size,
                   FUN = compute_betas,
                   by.column = FALSE,
                   align = "right")

# Create a data frame to hold the results
result <- data.frame(
  Time = seq(window_size, nrow(test_df2)),
  Beta_rsv = sapply(betas, `[`, 1),
  Beta_flu = sapply(betas, `[`, 2)
)

result

date_vec = test_df %>% select(week_end) %>% pull()
date_vec2 = date_vec[-c(1, 2)]

result$date = date_vec2

result %>%
  ggplot(aes(x = date)) +
  geom_line(aes(y = Beta_rsv, color = "RSV")) +
  geom_line(aes(y = Beta_flu, color = "Flu"))

head(tmp)
names(tmp) = c("date", "ili", "beta_rsv", "beta_flu")
head(tmp)

tmp %>% mutate(date = as_date(date)) %>%
  ggplot(aes(x = date)) +
  geom_line(aes(y = ili*10, color = "ILI")) +
  geom_line(aes(y = beta_rsv/5, color = "RSV")) +
  geom_line(aes(y = beta_flu*5, color = "Flu")) +
  scale_color_manual(values = c("ILI" = "grey",
                                "RSV" = "dodgerblue4",
                                "Flu" = "darkred",
                                "COVID" = "gold2"))
  # geom_line(aes(y = beta_flu, color = "Flu"))

tmp
# Create a data frame to hold the results
result <- data.frame(
  Time = seq(window_size, nrow(data)),
  Beta = beta_estimates
)

print(result)



## Original facet plot (unmodified time series) --------------------------------
nssp_ILI_RSV %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = unweighted_ili, color = "ILI")) +
  geom_line(aes(y = percent_visits_rsv, color = "RSV")) +
  # geom_text(data = cor_df,
  #           aes(x = as.Date("2024-04-01"), y = 18,
  #               label = paste("rho", "==", round(correlation, 2))),
  #           parse = T, size = 2) +
  # geom_text(data = cor_df,
  #           aes(x = as.Date("2024-04-01"), y = 14,
  #               label = paste("rho[lag1]", "==", round(lag1_correlation, 2))),
  #           parse = T, size = 2) +
  facet_wrap(~geography) +
  labs(x = "", y = "Percent ED Visits") +
  scale_color_manual(values = c("ILI" = "dodgerblue2", "RSV" = "gold3")) +
  scale_y_continuous(breaks = c(0, 10, 20)) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "none")

## Facet plot using rescaled time series (% of max) ----------------------------
cor_df = cor_df %>% left_join(state_pops, by = c("geography" = "state")) %>%
  arrange(population_rank)

cor_df = cor_df %>%
  mutate(geography_f = factor(geography,
                              levels = c(state_pops$state,
                                         "District of Columbia")))

nssp_ILI_RSV %>%
  mutate(geography_f = factor(geography, levels = c(state_pops$state,
                                                    "District of Columbia"))) %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = ili_rescaled, color = "ILI Visits")) +
  geom_line(aes(y = rsv_rescaled, color = "RSV ED Visits")) +
  geom_line(aes(y = vir_rescaled, color = "ILI Virology")) +
  # geom_text(data = cor_df,
  #           aes(x = as.Date("2024-05-01"), y = 0.95,
  #               label = paste("rho", "==", round(correlation, 2))),
  #           parse = T, size = 1.5) +
  # geom_text(data = cor_df,
  #           aes(x = as.Date("2024-04-01"), y = 0.8,
  #               label = paste("rho[lag1]", "==", round(lag1_correlation, 2))),
  #           parse = T, size = 1.5) +
  facet_wrap(~geography_f) +
  labs(x = "", y = "% of Max Signal") +
  scale_color_manual(values = c("ILI Visits" = "dodgerblue2",
                                "RSV ED Visits" = "gold3",
                                "ILI Virology" = "green4")) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  scale_x_date(date_breaks = "year") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "right")

setdiff(datasets::state.name, nssp_ILI_RSV$geography %>% drop_na() %>% unique())

range(nssp_ILI_RSV$week_end)

nssp_data

ggsave(file = "data/figures/rescaled_RSV_ILI_correlations.png", width = 10, height = 7, units = "in", bg = "white")

## Facet plot using normalized time series -------------------------------------
bc_transformed_df %>%
  mutate(geography_f = factor(geography, levels = c(state_pops$state,
                                                    "District of Columbia"))) %>%
  ggplot(aes(x = date)) +
  geom_line(aes(y = ili_normalized, color = "ILI")) +
  geom_line(aes(y = rsv_normalized, color = "RSV")) +
  geom_line(aes(y = vir_normalized, color = "VIR")) +
  # geom_text(data = cor_df,
  #           aes(x = as.Date("2022-12-22"), y = -2,
  #               label = paste("rho", "==", round(correlation, 2))),
  #           parse = T, size = 1.75) +
  # geom_text(data = cor_df,
  #           aes(x = as.Date("2023-01-15"), y = -3,
  #               label = paste("rho[lag1]", "==", round(lag1_correlation, 2))),
  #           parse = T, size = 1.75) +
  facet_wrap(~geography_f) +
  labs(x = "", y = "Normalized Signal Value") +
  scale_color_manual(values = c("ILI" = "black", "RSV" = "dodgerblue3", "VIR" = "red")) +
  scale_y_continuous(breaks = c(-3, 1.5, 0, -1.5, 3), limits = c(-3, 3)) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "right")
ggsave(file = "data/figures/normalized_RSV_ILI_VIR_v2.png", width = 10, height = 7, units = "in", bg = "white")


## testing zone
ili_data_RSV %>%
  filter(region == "Maine") %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = unweighted_ili, color = "ILI")) +
  # facet_wrap(~region) +
  labs(x = "", y = "Percent ED Visits") +
  theme(axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "none")

nssp_data %>%
  filter(geography == "New Jersey", county == "All") %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = percent_visits_rsv, color = "RSV")) +
  # facet_wrap(~region) +
  labs(x = "", y = "Percent ED Visits") +
  theme(axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "none")

ili_virology_data %>%
  filter(region == "New Jersey") %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = percent_positive, color = "virology")) +
  # facet_wrap(~region) +
  labs(x = "", y = "Percent ED Visits") +
  theme(axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "none")



