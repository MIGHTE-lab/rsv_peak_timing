## -----------------------------------------------------------------------------
## Script name: stacked_regression_non_negative.R
##
## Purpose of script: Use best estimated ridge regression to generate the stacked 
## curves for ID signals
##
## Author: George Dewey
##
## Date Created: 2024-07-15
##
## Last Updated: 2024-07-25
## -----------------------------------------------------------------------------

## Load packages

## Data management
library(tidyverse)
library(zoo)

## Working with dates
library(lubridate)
library(MMWRweek)

## Normalize data
library(bestNormalize)

## Modeling
library(roll)
library(glmnet)
library(penalized)
library(caret)

## Visualization
library(extrafont)

## Set the Google Drive folder as the working directory
setwd('~/Google Drive/My Drive/EWS_GW_BP_LC_RG')

## Load data - read in the combined data from `compare_ID_peaks.R`
signals_21_22 = read_csv('data/ground_truth/signals_21_22.csv', 
                         show_col_types = F)
signals_22_23 = read_csv('data/ground_truth/signals_22_23.csv',
                         show_col_types = F)
signals_23_24 = read_csv('data/ground_truth/signals_23_24.csv',
                         show_col_types = F)

fluSurv_states = signals_21_22 %>% select(geography) %>% unique() %>% pull()
states_data_available = signals_22_23 %>% select(geography) %>% unique() %>% pull()

# 1. Best ridge estimator (penalized ridge with non-negative constraint) -------

### 21-22 season: ili ~ flu + covid ---------------------------------------------
# Do the modeling with penalized ridge instead of LM
signals_21_22_with_betas = NULL
for(state in fluSurv_states[-14]){
  
  # Filter the data to only include data from the selected state
  data = signals_21_22 %>% filter(geography == state) %>% 
    select(ili, flu, covid)

  # Use 10-fold CV to identify best lambda
  cv = cv.glmnet(x = data.matrix(data[,c('flu', 'covid')]), 
                 y = data.matrix(data[,'ili']), 
                 intercept = FALSE, alpha = 0)

  # Use penalized package to run the ridge w/ no intercept + non-negative 
  # and the specified lambda from the CV model
  fit = penalized(ili ~ flu + covid , 
                  positive = TRUE,
                  unpenalized = ~0,
                  data = signals_21_22 %>% filter(geography == state),
                  lambda2 = cv$lambda.min)
  
  coefs_tmp = coefficients(fit)
  
  signals_tmp = signals_21_22 %>% 
    filter(geography == state) %>%
    mutate(beta_flu = coefs_tmp[1],
           beta_covid = coefs_tmp[2],
           flu_times_coef = flu*beta_flu,
           cov_times_coef = covid*beta_covid)
  
  signals_21_22_with_betas = bind_rows(signals_21_22_with_betas, signals_tmp)
}

#### 21/22 Season Figure  -----------------------------------------------------

signals_21_22_with_betas %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = flu_times_coef, color = 'Flu')) +
  geom_line(aes(y = cov_times_coef, color = 'COVID')) +
  geom_line(aes(y = ili, color = 'ILI')) +
  facet_wrap(~geography, scales = 'free', ncol = 5) +
  scale_color_manual(values = c('ILI' = 'grey',
                                'RSV' = 'dodgerblue4',
                                'Flu' = 'darkred',
                                'COVID' = 'gold2')) +
  theme_minimal() +
  scale_x_date(date_labels = '%b\n%y', 
               limits = c(as_date('2021-10-01'), as_date('2022-6-11')),
               date_breaks = '2 months') +
  theme_minimal() +
  labs(x = '', y = '') +
  theme(axis.text.x = element_text(size = 5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        legend.position = 'right')

ggsave('data/figures/peak_decomposition/cv_ridge/signals_21_22_best_ridge_cv.png',
       width = 10, height = 6, unit = 'in', bg = 'white')

### 22-23 season: ili ~ flu + rsv + covid ---------------------------------------
range(signals_22_23$week_end)
signals_22_23_with_betas = NULL
for(state in states_data_available){
  
  # Filter the data to only include data from the selected state
  data = signals_22_23 %>% filter(geography == state) %>% 
    select(ili, flu, covid, rsv)

  # Use 10-fold CV to identify best lambda
  cv = cv.glmnet(x = data.matrix(data[,c('flu', 'covid', 'rsv')]), 
                 y = data.matrix(data[,'ili']), 
                 intercept = FALSE, alpha = 0)

  # Use penalized package to run the ridge w/ no intercept + non-negative 
  # and the specified lambda from the CV model
  fit = penalized(ili ~ flu + covid + rsv, 
                  positive = TRUE,
                  unpenalized = ~0,
                  data = signals_22_23 %>% filter(geography == state),
                  lambda2 = cv$lambda.min)
  coefs_tmp = coefficients(fit)

  if(state == 'California'){
    signals_tmp = signals_22_23 %>% 
      filter(geography == state) %>%
      mutate(beta_covid = coefs_tmp[1],
             beta_rsv = coefs_tmp[2],
             cov_times_coef = covid*beta_covid,
             rsv_times_coef = rsv*beta_rsv)
  }
  else {
    signals_tmp = signals_22_23 %>%
    filter(geography == state) %>%
      mutate(beta_flu = coefs_tmp[1],
             beta_covid = coefs_tmp[2],
             beta_rsv = coefs_tmp[3],
             flu_times_coef = flu*beta_flu,
             cov_times_coef = covid*beta_covid,
             rsv_times_coef = rsv*beta_rsv)
  }
    
  signals_22_23_with_betas = bind_rows(signals_22_23_with_betas, signals_tmp)
}

#### 22/23 Season Figure A -----------------------------------------------------
signals_22_23_with_betas %>%
  filter(geography %in% states_data_available[1:24]) %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = flu_times_coef, color = 'Flu')) +
  geom_line(aes(y = cov_times_coef, color = 'COVID')) +
  geom_line(aes(y = rsv_times_coef, color = 'RSV')) +
  geom_line(aes(y = ili, color = 'ILI')) +
  facet_wrap(~geography, scales = 'free', ncol = 6) +
  scale_color_manual(values = c('ILI' = 'grey',
                                'RSV' = 'dodgerblue4',
                                'Flu' = 'darkred',
                                'COVID' = 'gold2')) +
  theme_minimal() +
  scale_x_date(date_labels = '%b\n%y', 
               date_breaks = '2 months') +
  labs(x = '', y = '') +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 5),
        legend.title = element_blank(),
        legend.position = 'right')

ggsave('data/figures/signals_22_23_best_ridge_cv_A.png',
       width = 10, height = 6, unit = 'in', bg = 'white')

#### 22/23 Season Figure B -----------------------------------------------------

signals_22_23_with_betas %>%
  filter(geography %in% states_data_available[25:47]) %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = flu_times_coef, color = 'Flu')) +
  geom_line(aes(y = cov_times_coef, color = 'COVID')) +
  geom_line(aes(y = rsv_times_coef, color = 'RSV')) +
  geom_line(aes(y = ili, color = 'ILI')) +
  facet_wrap(~geography, scales = 'free', ncol = 6) +
  scale_color_manual(values = c('ILI' = '#919191',
                                'RSV' = '#20AFD1',
                                'Flu' = '#FF5408',
                                'COVID' = '#FFCe00')) +
  theme_minimal() +
  scale_x_date(date_labels = '%b\n%y', 
               date_breaks = '2 months') +
  theme_minimal() +
  labs(x = '', y = '') +
  theme(axis.text.x = element_text(size = 5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 5),
        legend.title = element_blank(),
        legend.position = 'right')

ggsave('data/figures/signals_22_23_best_ridge_cv_B.png',
       width = 10, height = 6, unit = 'in', bg = 'white')

#### 22/23 Statewise Figures ---------------------------------------------------

## For presentation
signals_22_23_with_betas %>%
  filter(geography == 'Maryland') %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = flu_times_coef, color = 'Flu'), linewidth = 1.2, alpha = 0.8) +
  geom_line(aes(y = cov_times_coef, color = 'COVID'), linewidth = 1.2, alpha = 0.8) +
  geom_line(aes(y = rsv_times_coef, color = 'RSV'), linewidth = 1.2, alpha = 0.8) +
  geom_line(aes(y = ili, color = 'ILI'), alpha = 0) +
  geom_area(aes(y = ili), fill = '#d6ede9', alpha = 0.4) +
  scale_color_manual(values = c('ILI' = '#d6ede9',
                                'RSV' = '#20afd1',
                                'Flu' = '#cc3e00',
                                'COVID' = '#ffce00')) +
  scale_x_date(date_labels = '%b\n%Y', 
               date_breaks = '2 months') +
  theme_fira() +
  labs(x = '', y = '') +
  theme(plot.title = element_text(size = 24, family = 'Barlow Semi Condensed'),
        axis.text.x = element_text(size = 16, family = 'Barlow Semi Condensed'),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 16, family = 'Barlow Semi Condensed'),
        legend.position = 'none')+
  guides(col = guide_legend(override.aes = list(alpha = 1)))

ggsave(file = 'data/figures/epistorm_pres_2024_07_30/MD_22_23_season_single_no_legend_6x5.png',
       width = 6, height = 5, units = 'in', bg = 'white')


signals_22_23_with_betas %>%
  filter(geography == 'New York') %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = flu_times_coef, color = 'Flu')) +
  geom_line(aes(y = cov_times_coef, color = 'COVID')) +
  geom_line(aes(y = rsv_times_coef, color = 'RSV')) +
  geom_line(aes(y = ili, color = 'ILI'), alpha = 0) +
  geom_area(aes(y = ili), fill = '#d6ede9', alpha = 0.4) +
  scale_color_manual(values = c('ILI' = '#d6ede9',
                                'RSV' = '#20afd1',
                                'Flu' = '#cc3e00',
                                'COVID' = '#ffce00')) +
  scale_x_date(date_labels = '%b\n%Y', 
               date_breaks = '2 months') +
  theme_fira() +
  labs(x = '', y = '', title = 'New York') +
  theme(plot.title = element_text(size = 24, family = 'Barlow Semi Condensed'),
        axis.text.x = element_text(size = 16, family = 'Barlow Semi Condensed'),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 16, family = 'Barlow Semi Condensed'),
        legend.position = 'bottom')+
  guides(col = guide_legend(override.aes = list(alpha = 1)))

ggsave(file = 'data/figures/epistorm_pres_2024_07_30/NY_22_23_season_single.png',
       width = 8, height = 6, units = 'in', bg = 'white')

signals_22_23_with_betas %>%
  filter(geography == 'Pennsylvania') %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = flu_times_coef, color = 'Flu')) +
  geom_line(aes(y = cov_times_coef, color = 'COVID')) +
  geom_line(aes(y = rsv_times_coef, color = 'RSV')) +
  geom_line(aes(y = ili, color = 'ILI'), alpha = 0) +
  geom_area(aes(y = ili), fill = '#d6ede9', alpha = 0.4) +
  scale_color_manual(values = c('ILI' = '#d6ede9',
                                'RSV' = '#20afd1',
                                'Flu' = '#cc3e00',
                                'COVID' = '#ffce00')) +
  scale_x_date(date_labels = '%b\n%Y', 
               date_breaks = '2 months') +
  theme_fira() +
  labs(x = '', y = '', title = 'Pennsylvania') +
  theme(plot.title = element_text(size = 24, family = 'Barlow Semi Condensed'),
        axis.text.x = element_text(size = 16, family = 'Barlow Semi Condensed'),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 16, family = 'Barlow Semi Condensed'),
        legend.position = 'bottom') + 
  guides(col = guide_legend(override.aes = list(alpha = 1)))

ggsave(file = 'data/figures/epistorm_pres_2024_07_30/PA_22_23_season_single.png',
       width = 8, height = 6, units = 'in', bg = 'white')

signals_22_23_with_betas %>%
  filter(geography == 'Wisconsin') %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = flu_times_coef, color = 'Flu')) +
  geom_line(aes(y = cov_times_coef, color = 'COVID')) +
  geom_line(aes(y = rsv_times_coef, color = 'RSV')) +
  geom_line(aes(y = ili, color = 'ILI'), alpha = 0) +
  geom_area(aes(y = ili), fill = '#d6ede9', alpha = 0.4) +
  scale_color_manual(values = c('ILI' = '#d6ede9',
                                'RSV' = '#20afd1',
                                'Flu' = '#cc3e00',
                                'COVID' = '#ffce00')) +
  scale_x_date(date_labels = '%b\n%Y', 
               date_breaks = '2 months') +
  theme_fira() +
  labs(x = '', y = '', title = 'Wisconsin') +
  theme(plot.title = element_text(size = 24, family = 'Barlow Semi Condensed'),
        axis.text.x = element_text(size = 16, family = 'Barlow Semi Condensed'),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 16, family = 'Barlow Semi Condensed'),
        legend.position = 'bottom') + 
  guides(col = guide_legend(override.aes = list(alpha = 1)))

ggsave(file = 'data/figures/epistorm_pres_2024_07_30/WI_22_23_season_single.png',
       width = 8, height = 6, units = 'in', bg = 'white')

### 23-24 season: ili ~ flu + rsv -----------------------------------------------
signals_23_24_with_betas = NULL
for(state in states_data_available){
  
  # Filter the data to only include data from the selected state
  data = signals_23_24 %>% filter(geography == state) %>% 
    select(ili, flu, rsv)
  
  # Use 10-fold CV to identify best lambda
  cv = cv.glmnet(x = data.matrix(data[,c('flu', 'rsv')]), 
                 y = data.matrix(data[,'ili']), 
                 intercept = FALSE, alpha = 0)
  
  # Use penalized package to run the ridge w/ no intercept + non-negative 
  # and the specified lambda from the CV model
  
  fit = penalized(ili ~ flu + rsv, 
                  positive = TRUE,
                  unpenalized = ~0,
                  data = signals_23_24 %>% filter(geography == state),
                  lambda2 = cv$lambda.min)
  coefs_tmp = coefficients(fit)

  # if(state == 'California'){
  #   signals_tmp = signals_23_23 %>% 
  #     filter(geography == state) %>%
  #     mutate(beta_covid = coefs_tmp[1],
  #            beta_rsv = coefs_tmp[2],
  #            cov_times_coef = covid*beta_covid,
  #            rsv_times_coef = rsv*beta_rsv)
  # }
  # else {
    signals_tmp = signals_23_24 %>%
      filter(geography == state) %>%
      mutate(beta_flu = coefs_tmp[1],
             beta_rsv = coefs_tmp[2],
             flu_times_coef = flu*beta_flu,
             rsv_times_coef = rsv*beta_rsv)
  
  
  signals_23_24_with_betas = bind_rows(signals_23_24_with_betas, signals_tmp)
}

#### 23/24 Season Figure A -----------------------------------------------------
signals_23_24_with_betas %>%
  filter(geography %in% states_data_available[1:24]) %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = flu_times_coef, color = 'Flu')) +
  geom_line(aes(y = rsv_times_coef, color = 'RSV')) +
  geom_line(aes(y = ili, color = 'ILI')) +
  facet_wrap(~geography, scales = 'free', ncol = 6) +
  scale_color_manual(values = c('ILI' = 'grey',
                                'RSV' = 'dodgerblue4',
                                'Flu' = 'darkred',
                                'COVID' = 'gold2')) +
  theme_minimal() +
  scale_x_date(date_labels = '%b\n%y', 
               date_breaks = '2 months') +
  labs(x = '', y = '') +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 5),
        legend.title = element_blank(),
        legend.position = 'right')

ggsave('data/figures/signals_23_24_best_ridge_cv_A.png',
       width = 10, height = 6, unit = 'in', bg = 'white')

#### 23/24 Season Figure B -----------------------------------------------------

signals_23_24_with_betas %>%
  filter(geography %in% states_data_available[25:47]) %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = flu_times_coef, color = 'Flu')) +
  geom_line(aes(y = rsv_times_coef, color = 'RSV')) +
  geom_line(aes(y = ili, color = 'ILI')) +
  facet_wrap(~geography, scales = 'free', ncol = 6) +
  scale_color_manual(values = c('ILI' = 'grey',
                                'RSV' = 'dodgerblue4',
                                'Flu' = 'darkred',
                                'COVID' = 'gold2')) +
  theme_minimal() +
  scale_x_date(date_labels = '%b\n%y', 
               date_breaks = '2 months') +
  theme_minimal() +
  labs(x = '', y = '') +
  theme(axis.text.x = element_text(size = 5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size = 5),
        legend.title = element_blank(),
        legend.position = 'right')

ggsave('data/figures/signals_23_24_best_ridge_cv_B.png',
       width = 10, height = 6, unit = 'in', bg = 'white')


## 2. Stacked regression using linear models as the base estimators -----------

### 21/22 Season --------------------------------------------------------------
signals_21_22_with_stacked_lm_betas = NULL
for(state in fluSurv_states[-14]){
  
  # Filter the data so that it only includes the selected state's data
  data = signals_21_22 %>% filter(geography == state) %>% 
    select(flu, covid, ili)
  
  # Split data into training and testing sets
  set.seed(123)
  trainIndex = createDataPartition(data$ili, p = .8, 
                                   list = FALSE, 
                                   times = 1)
  trainData = data[trainIndex, ]
  testData  = data[-trainIndex, ]
  
  # Train two base linear regression models - each has only one predictor
  model_lm1 = lm(ili ~ flu + 0, data = trainData)
  model_lm2 = lm(ili ~ covid + 0, data = trainData)
  coef_lm1 = coef(model_lm1)
  coef_lm2 = coef(model_lm2)
  
  # Extract predictions from base models
  trainData$pred_lm1 = predict(model_lm1, trainData)
  trainData$pred_lm2 = predict(model_lm2, trainData)
  
  testData$pred_lm1 = predict(model_lm1, testData)
  testData$pred_lm2 = predict(model_lm2, testData)
  
  # Prepare data for the meta-model
  meta_x_train = as.matrix(trainData[, c("pred_lm1", "pred_lm2")])
  meta_y_train = trainData$ili
  meta_x_test = as.matrix(testData[, c("pred_lm1", "pred_lm2")])
  
  # Train the meta-model
  meta_model = lm(meta_y_train ~ meta_x_train + 0)
  
  meta_coefficients = coef(meta_model)
  
  final_coefficients = 
    c(meta_coefficients["meta_x_trainpred_lm1"] * coef_lm1["flu"],
      meta_coefficients["meta_x_trainpred_lm2"] * coef_lm2["covid"])
  
  signals_tmp = signals_21_22 %>%
    filter(geography == state) %>%
    mutate(beta_flu = final_coefficients[1],
           beta_covid = final_coefficients[2],
           flu_times_coef = flu*beta_flu,
           covid_times_coef = covid*beta_covid)
  
  signals_21_22_with_stacked_lm_betas = 
    bind_rows(signals_21_22_with_stacked_lm_betas, signals_tmp)
}

#### 21/22 Season Figure  -----------------------------------------------------
signals_21_22_with_stacked_lm_betas %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = flu_times_coef, color = 'Flu')) +
  geom_line(aes(y = covid_times_coef, color = 'COVID')) +
  geom_line(aes(y = ili, color = 'ILI')) +
  facet_wrap(~geography, scales = 'free', ncol = 5) +
  scale_color_manual(values = c('ILI' = 'grey',
                                'RSV' = 'dodgerblue4',
                                'Flu' = 'darkred',
                                'COVID' = 'gold2')) +
  theme_minimal() +
  scale_x_date(date_labels = '%b\n%y', 
               limits = c(as_date('2021-10-01'), '2022-6-11'),
               date_breaks = '2 months') +
  theme_minimal() +
  labs(x = '', y = '') +
  theme(axis.text.x = element_text(size = 5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        legend.position = 'right')
signals_21_22_with_stacked_lm_betas
signals_21_22_with_betas

ggsave('data/figures/signals_21_22_stacked_lm.png',
       width = 10, height = 6, unit = 'in', bg = 'white')

### 22/23 Season ---------------------------------------------------------------
signals_22_23_with_stacked_lm_betas = NULL
for(state in states_data_available){
  
  # Filter the data so that it only includes the selected state's data
  data = signals_22_23 %>% filter(geography == state) %>% 
    select(flu, covid, rsv, ili)
  
  # Split data into training and testing sets
  set.seed(123)
  trainIndex = createDataPartition(data$ili, p = .8, 
                                   list = FALSE, 
                                   times = 1)
  trainData = data[trainIndex, ]
  testData  = data[-trainIndex, ]
  
  # Train three base linear regression models - each has only one predictor
  model_lm1 = lm(ili ~ flu + 0, data = trainData)
  model_lm2 = lm(ili ~ covid + 0, data = trainData)
  model_lm3 = lm(ili ~ rsv + 0, data = trainData)
  coef_lm1 = coef(model_lm1)
  coef_lm2 = coef(model_lm2)
  coef_lm3 = coef(model_lm3)
  
  # Extract predictions from base models
  trainData$pred_lm1 = predict(model_lm1, trainData)
  trainData$pred_lm2 = predict(model_lm2, trainData)
  trainData$pred_lm3 = predict(model_lm3, trainData)
  
  testData$pred_lm1 = predict(model_lm1, testData)
  testData$pred_lm2 = predict(model_lm2, testData)
  testData$pred_lm3 = predict(model_lm3, testData)
  
  # Prepare data for the meta-model
  meta_x_train = as.matrix(trainData[, c("pred_lm1", "pred_lm2", "pred_lm3")])
  meta_y_train = trainData$ili
  meta_x_test = as.matrix(testData[, c("pred_lm1", "pred_lm2", "pred_lm3")])
  
  # Train the meta-model
  meta_model = lm(meta_y_train ~ meta_x_train + 0)
  
  meta_coefficients = coef(meta_model)
  
  final_coefficients = 
    c(meta_coefficients["meta_x_trainpred_lm1"] * coef_lm1["flu"],
      meta_coefficients["meta_x_trainpred_lm2"] * coef_lm2["covid"],
      meta_coefficients["meta_x_trainpred_lm3"] * coef_lm3["rsv"])
  
  signals_tmp = signals_22_23 %>%
    filter(geography == state) %>%
    mutate(beta_flu = final_coefficients[1],
           beta_covid = final_coefficients[2],
           beta_rsv = final_coefficients[3],
           flu_times_coef = flu*beta_flu,
           covid_times_coef = covid*beta_covid,
           rsv_times_coef = rsv*beta_rsv)
  
  signals_22_23_with_stacked_lm_betas = 
    bind_rows(signals_22_23_with_stacked_lm_betas, signals_tmp)
}

#### 22/23 Season Figure A -----------------------------------------------------
signals_22_23_with_stacked_lm_betas %>%
  filter(geography %in% states_data_available[1:24]) %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = flu_times_coef, color = 'Flu')) +
  geom_line(aes(y = covid_times_coef, color = 'COVID')) +
  geom_line(aes(y = rsv_times_coef, color = 'RSV')) +
  geom_line(aes(y = ili, color = 'ILI')) +
  facet_wrap(~geography, scales = 'free', ncol = 5) +
  scale_color_manual(values = c('ILI' = 'grey',
                                'RSV' = 'dodgerblue4',
                                'Flu' = 'darkred',
                                'COVID' = 'gold2')) +
  theme_minimal() +
  scale_x_date(date_labels = '%b\n%y', 
               date_breaks = '2 months') +
  theme_minimal() +
  labs(x = '', y = '') +
  theme(axis.text.x = element_text(size = 5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        legend.position = 'right')

ggsave('data/figures/peak_decomposition/stacked_lm/signals_22_23_stacked_lm_A.png',
       width = 10, height = 6, unit = 'in', bg = 'white')

#### 22/23 Season Figure B -----------------------------------------------------

signals_22_23_with_stacked_lm_betas %>%
  filter(geography %in% states_data_available[25:47]) %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = flu_times_coef, color = 'Flu')) +
  geom_line(aes(y = covid_times_coef, color = 'COVID')) +
  geom_line(aes(y = rsv_times_coef, color = 'RSV')) +
  geom_line(aes(y = ili, color = 'ILI')) +
  facet_wrap(~geography, scales = 'free', ncol = 5) +
  scale_color_manual(values = c('ILI' = 'grey',
                                'RSV' = 'dodgerblue4',
                                'Flu' = 'darkred',
                                'COVID' = 'gold2')) +
  theme_minimal() +
  scale_x_date(date_labels = '%b\n%y', 
               date_breaks = '2 months') +
  theme_minimal() +
  labs(x = '', y = '') +
  theme(axis.text.x = element_text(size = 5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        legend.position = 'right')

ggsave('data/figures/peak_decomposition/stacked_lm/signals_22_23_stacked_lm_B.png',
       width = 10, height = 6, unit = 'in', bg = 'white')

### 23/24 season ----------------------------------------------------------------
signals_23_24_with_stacked_lm_betas = NULL
for(state in states_data_available){
  
  # Filter the data so that it only includes the selected state's data
  data = signals_23_24 %>% filter(geography == state) %>% 
    select(flu, rsv, ili)
  
  # Split data into training and testing sets
  set.seed(123)
  trainIndex = createDataPartition(data$ili, p = .8, 
                                   list = FALSE, 
                                   times = 1)
  trainData = data[trainIndex, ]
  testData  = data[-trainIndex, ]
  
  # Train three base linear regression models - each has only one predictor
  model_lm1 = lm(ili ~ flu + 0, data = trainData)
  model_lm2 = lm(ili ~ rsv + 0, data = trainData)
  coef_lm1 = coef(model_lm1)
  coef_lm2 = coef(model_lm2)

  # Extract predictions from base models
  trainData$pred_lm1 = predict(model_lm1, trainData)
  trainData$pred_lm2 = predict(model_lm2, trainData)

  testData$pred_lm1 = predict(model_lm1, testData)
  testData$pred_lm2 = predict(model_lm2, testData)

  # Prepare data for the meta-model
  meta_x_train = as.matrix(trainData[, c("pred_lm1", "pred_lm2")])
  meta_y_train = trainData$ili
  meta_x_test = as.matrix(testData[, c("pred_lm1", "pred_lm2")])
  
  # Train the meta-model
  meta_model = lm(meta_y_train ~ meta_x_train + 0)
  
  meta_coefficients = coef(meta_model)
  
  final_coefficients = 
    c(meta_coefficients["meta_x_trainpred_lm1"] * coef_lm1["flu"],
      meta_coefficients["meta_x_trainpred_lm2"] * coef_lm2["rsv"])
  
  signals_tmp = signals_23_24 %>%
    filter(geography == state) %>%
    mutate(beta_flu = final_coefficients[1],
           beta_rsv = final_coefficients[2],
           flu_times_coef = flu*beta_flu,
           rsv_times_coef = rsv*beta_rsv)
  
  signals_23_24_with_stacked_lm_betas = 
    bind_rows(signals_23_24_with_stacked_lm_betas, signals_tmp)
}

#### 22/23 Season Figure A -----------------------------------------------------
signals_23_24_with_stacked_lm_betas %>%
  filter(geography %in% states_data_available[1:24]) %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = flu_times_coef, color = 'Flu')) +
  geom_line(aes(y = rsv_times_coef, color = 'RSV')) +
  geom_line(aes(y = ili, color = 'ILI')) +
  facet_wrap(~geography, scales = 'free', ncol = 6) +
  scale_color_manual(values = c('ILI' = 'grey',
                                'RSV' = 'dodgerblue4',
                                'Flu' = 'darkred',
                                'COVID' = 'gold2')) +
  theme_minimal() +
  scale_x_date(date_labels = '%b\n%y', 
               date_breaks = '2 months') +
  theme_minimal() +
  labs(x = '', y = '') +
  theme(axis.text.x = element_text(size = 5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        legend.position = 'right')

ggsave('data/figures/peak_decomposition/stacked_lm/signals_23_24_stacked_lm_A.png',
       width = 10, height = 6, unit = 'in', bg = 'white')

#### 23/24 Season Figure B -----------------------------------------------------
signals_23_24_with_stacked_lm_betas %>%
  filter(geography %in% states_data_available[25:47]) %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = flu_times_coef, color = 'Flu')) +
  geom_line(aes(y = rsv_times_coef, color = 'RSV')) +
  geom_line(aes(y = ili, color = 'ILI')) +
  facet_wrap(~geography, scales = 'free', ncol = 6) +
  scale_color_manual(values = c('ILI' = 'grey',
                                'RSV' = 'dodgerblue4',
                                'Flu' = 'darkred',
                                'COVID' = 'gold2')) +
  theme_minimal() +
  scale_x_date(date_labels = '%b\n%y', 
               date_breaks = '2 months') +
  theme_minimal() +
  labs(x = '', y = '') +
  theme(axis.text.x = element_text(size = 5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        legend.position = 'right')

ggsave('data/figures/peak_decomposition/stacked_lm/signals_23_24_stacked_lm_B.png',
       width = 10, height = 6, unit = 'in', bg = 'white')

# 3. Stacked regression using ridge as the meta estimator -------------

signals_21_22_with_stacked_ridge_betas = NULL
for(state in fluSurv_states[-14]){
    
    # Filter the data to only include data from the selected state
    data = signals_21_22 %>% filter(geography == state) %>% 
      select(ili, flu, covid)
    
    trainIndex <- createDataPartition(data$ili, p = .8, 
                                      list = FALSE, 
                                      times = 1)
    trainData <- data[trainIndex, ]
    testData  <- data[-trainIndex, ]
    
    x_train = as.matrix(trainData[, c('flu', 'covid')])
    y_train = trainData$ili
    
    # Train first linear regression model
    model_lm1 <- lm(ili ~ flu, data = trainData)
    
    # Train second linear regression model
    model_lm2 <- lm(ili ~ covid, data = trainData)
    
    # Extract predictions from base models
    trainData$pred_lm1 <- predict(model_lm1, trainData)
    trainData$pred_lm2 <- predict(model_lm2, trainData)
    
    testData$pred_lm1 <- predict(model_lm1, testData)
    testData$pred_lm2 <- predict(model_lm2, testData)
    
    # Prepare data for the meta-model
    meta_x_train <- as.matrix(trainData[, c("pred_lm1", "pred_lm2")])
    meta_y_train <- trainData$ili
    
    # Train the meta-model (ridge regression)
    meta_model <- cv.glmnet(meta_x_train, meta_y_train, intercept = FALSE,
                            alpha = 0)
    meta_y_train
    meta_x_train
    meta_data = bind_cols(meta_x_train, meta_y_train)
    names(meta_data) = c('pred_lm1', 'pred_lm2', 'response')
    meta_model2 = penalized(response ~ pred_lm1 + pred_lm2, unpenalized = ~0,
                            positive = TRUE, lambda2 = meta_model$lambda.min,
                            data = meta_data)

    # Make predictions using the meta-model
    meta_x_test <- as.matrix(testData[, c("pred_lm1", "pred_lm2")])

    meta_coefficients = coefficients(meta_model2)
    
    signals_tmp = signals_21_22 %>% 
      filter(geography == state) %>%
      mutate(beta_flu = meta_coefficients[1],
             beta_covid = meta_coefficients[2],
             flu_times_coef = flu*beta_flu,
             cov_times_coef = covid*beta_covid)
    
    signals_21_22_with_stacked_ridge_betas = 
      bind_rows(signals_21_22_with_stacked_ridge_betas, signals_tmp)
}

signals_21_22_with_stacked_ridge_betas %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = flu_times_coef, color = 'Flu')) +
  geom_line(aes(y = cov_times_coef, color = 'COVID')) +
  geom_line(aes(y = ili, color = 'ILI')) +
  facet_wrap(~geography, scales = 'free', ncol = 5) +
  scale_color_manual(values = c('ILI' = 'grey',
                                'RSV' = 'dodgerblue4',
                                'Flu' = 'darkred',
                                'COVID' = 'gold2')) +
  theme_minimal() +
  scale_x_date(date_labels = '%b\n%y', 
               limits = c(as_date('2021-10-01'), '2022-6-11'),
               date_breaks = '2 months') +
  theme_minimal() +
  labs(x = '', y = '') +
  theme(axis.text.x = element_text(size = 5),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        legend.position = 'right')
