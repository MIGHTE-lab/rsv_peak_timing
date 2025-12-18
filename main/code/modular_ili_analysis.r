## -----------------------------------------------------------------------------
## Script name: modular_ili_analysis.R
##
## Purpose of script: Modular analysis of ILI time series with multiple
## modeling approaches
##
## Author: George Dewey
##
## Date Created: 2025-10-01
##
## Last Updated: 2025-10-01
##
## Version: 2.0.0
##
## Update notes: Added modular approach with three modeling options
## -----------------------------------------------------------------------------

## 1. Load packages and setup -------------------------------------------------
library(tidyverse)
library(zoo)
library(lubridate)
library(MMWRweek)
library(ggh4x)
library(glmnet)
library(caret)
library(ggplot2)

## Set working directory (user should modify as needed)
# setwd("~/Documents/Projects/rsv_peak_timing")

## Load and prepare data (assuming this part remains the same)
dat = read_csv("data/NSSP/nssp_full_040525.csv")

ili_data = read_csv("data/ILInet/ILINet_state_2025_03-29.csv",
                    na = "X", skip = 1)

names(ili_data) = c(
  "region_type", "region", "year", "epiweek", "weighted_ili", "unweighted_ili",
  "age_0_4", "age_25_49", "age_25_64", "age_5_24", "age_50_64", "age_65",
  "ilitotal", "num_data_providers", "total_patients"
)

ili_data_with_dates = ili_data %>%
  filter(year %in% 2022:2025) %>%
  mutate(week_end = as_date(MMWRweek2Date(year, epiweek) + 6)) %>%
  select(region, year, epiweek, week_end, unweighted_ili, ilitotal) %>%
  rename(state = region)

dat = dat %>%
  filter(county == "All") %>%
  select(week_end, geography, percent_visits_covid, percent_visits_influenza, percent_visits_rsv)

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

dat = dat %>%
  rename(
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

## 2. Modeling Functions ------------------------------------------------------

# Model 1: Ridge regression without intercept (original approach)
ridge_no_intercept = function(data) {
  predictors = c("covid_rescaled", "flu_rescaled", "rsv_rescaled")
  X = as.matrix(data[, predictors])
  y = data$ili_rescaled

  lambda_values = 10^seq(-5, 5, length = 100)
  cv_ridge = cv.glmnet(X, y, alpha = 0, lambda = lambda_values, nfolds = 10,
                       intercept = FALSE, lower.limits = 0)

  lambda_grid = seq(cv_ridge$lambda.min / 1000, cv_ridge$lambda.min * 1000, length = 500)
  ridge_models = lapply(lambda_grid, function(l) {
    glmnet(X, y, alpha = 0, lambda = l, intercept = FALSE, lower.limits = 0)
  })

  ridge_coefs = sapply(ridge_models, function(model) coef(model)[,1])
  valid_lambda_idx = which(colSums(ridge_coefs != 0) == length(predictors))[1]

  if (is.na(valid_lambda_idx)) {
    warning("No valid lambda found. Using lambda.min.")
    best_lambda = cv_ridge$lambda.min
  } else {
    best_lambda = lambda_grid[valid_lambda_idx]
  }

  ridge_final = glmnet(X, y, alpha = 0, lambda = best_lambda, intercept = FALSE, lower.limits = 0)
  coef_values = as.numeric(coef(ridge_final))

  return(c(0, coef_values[2:4])) # Return [intercept=0, covid, flu, rsv]
}

# Model 2: Ridge regression with intercept
ridge_with_intercept = function(data) {
  predictors = c("covid_rescaled", "flu_rescaled", "rsv_rescaled")
  X = as.matrix(data[, predictors])
  y = data$ili_rescaled

  lambda_values = 10^seq(-5, 5, length = 100)
  cv_ridge = cv.glmnet(X, y, alpha = 0, lambda = lambda_values, nfolds = 10,
                       intercept = TRUE, lower.limits = 0)

  lambda_grid = seq(cv_ridge$lambda.min / 1000, cv_ridge$lambda.min * 1000, length = 500)
  ridge_models = lapply(lambda_grid, function(l) {
    glmnet(X, y, alpha = 0, lambda = l, intercept = TRUE, lower.limits = 0)
  })

  ridge_coefs = sapply(ridge_models, function(model) coef(model)[,1])
  valid_lambda_idx = which(colSums(ridge_coefs != 0) == (length(predictors) + 1))[1]

  if (is.na(valid_lambda_idx)) {
    warning("No valid lambda found. Using lambda.min.")
    best_lambda = cv_ridge$lambda.min
  } else {
    best_lambda = lambda_grid[valid_lambda_idx]
  }

  ridge_final = glmnet(X, y, alpha = 0, lambda = best_lambda, intercept = TRUE, lower.limits = 0)
  coef_values = as.numeric(coef(ridge_final))

  return(coef_values) # Return [intercept, covid, flu, rsv]
}

# Model 3: Linear regression without regularization (lambda = 0)
linear_no_regularization = function(data) {
  predictors = c("covid_rescaled", "flu_rescaled", "rsv_rescaled")
  X = as.matrix(data[, predictors])
  y = data$ili_rescaled

  # Use glmnet with lambda = 0 (equivalent to OLS) and no intercept
  linear_model = glmnet(X, y, alpha = 0, lambda = 0, intercept = FALSE, lower.limits = 0)
  coef_values = as.numeric(coef(linear_model))

  return(c(0, coef_values[2:4])) # Return [intercept=0, covid, flu, rsv]
}

## 3. Main Analysis Function --------------------------------------------------

run_ili_analysis = function(model_type = "ridge_no_intercept",
                            seasons = c("22-23", "23-24", "24-25"),
                            save_outputs = TRUE,
                            output_prefix = NULL) {

  # Validate model type
  valid_models = c("ridge_without_intercept", "ridge_with_intercept", "linear_no_regularization")
  if (!model_type %in% valid_models) {
    stop("model_type must be one of: ", paste(valid_models, collapse = ", "))
  }

  # Set output prefix if not provided
  if (is.null(output_prefix)) {
    output_prefix = model_type
  }

  # Select modeling function
  model_function = switch(model_type,
                          "ridge_without_intercept" = ridge_no_intercept,
                          "ridge_with_intercept" = ridge_with_intercept,
                          "linear_no_regularization" = linear_no_regularization
  )

  # Initialize storage
  all_signals = NULL
  model_results_table = tibble()

  # Helper function for contribution statistics
  calculate_contribution_stats = function(contribution_series) {
    list(
      mean_contribution = mean(contribution_series, na.rm = TRUE),
      max_contribution = max(contribution_series, na.rm = TRUE),
      min_contribution = min(contribution_series, na.rm = TRUE),
      sum_contribution = sum(contribution_series, na.rm = TRUE),
      peak_week = which.max(contribution_series)[1]
    )
  }

  # Run analysis for each season
  for (current_season in seasons) {
    cat("Processing season:", current_season, "\n")

    # Get states for this season (handle Wyoming exclusion for 24-25)
    states_to_process = if (current_season == "24-25") {
      dat_combined_ %>% filter(state != "Wyoming") %>% distinct(state) %>% pull()
    } else {
      dat_combined_ %>% distinct(state) %>% pull()
    }

    for (geo in states_to_process) {

      data = dat_combined_ %>%
        filter(state == geo, season == current_season) %>%
        mutate(
          ili_rescaled = minmax(ili),
          flu_rescaled = minmax(flu),
          rsv_rescaled = minmax(rsv),
          covid_rescaled = minmax(covid)
        )

      if (nrow(data) < 10) {
        warning(paste("Insufficient data for", geo, "in", current_season, "season"))
        next
      }

      # Fit model
      tryCatch({
        coefs_tmp = model_function(data)

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

        # Calculate contribution statistics
        flu_stats = calculate_contribution_stats(signals_tmp$flu_times_coef)
        rsv_stats = calculate_contribution_stats(signals_tmp$rsv_times_coef)
        covid_stats = calculate_contribution_stats(signals_tmp$covid_times_coef)

        # Create summary row
        summary_row = tibble(
          state = geo,
          season = current_season,
          model_type = model_type,
          n_weeks = nrow(data),

          # Model coefficients
          beta_intercept = coefs_tmp[1],
          beta_covid = coefs_tmp[2],
          beta_flu = coefs_tmp[3],
          beta_rsv = coefs_tmp[4],

          # Flu contribution statistics
          flu_mean_contribution = flu_stats$mean_contribution,
          flu_max_contribution = flu_stats$max_contribution,
          flu_sum_contribution = flu_stats$sum_contribution,
          flu_peak_week = flu_stats$peak_week,

          # RSV contribution statistics
          rsv_mean_contribution = rsv_stats$mean_contribution,
          rsv_max_contribution = rsv_stats$max_contribution,
          rsv_sum_contribution = rsv_stats$sum_contribution,
          rsv_peak_week = rsv_stats$peak_week,

          # COVID contribution statistics
          covid_mean_contribution = covid_stats$mean_contribution,
          covid_max_contribution = covid_stats$max_contribution,
          covid_sum_contribution = covid_stats$sum_contribution,
          covid_peak_week = covid_stats$peak_week,

          # Additional diagnostics
          total_explained_variance = flu_stats$sum_contribution +
            rsv_stats$sum_contribution +
            covid_stats$sum_contribution,
          dominant_pathogen = case_when(
            flu_stats$max_contribution == max(flu_stats$max_contribution,
                                              rsv_stats$max_contribution,
                                              covid_stats$max_contribution) ~ "Flu",
            rsv_stats$max_contribution == max(flu_stats$max_contribution,
                                              rsv_stats$max_contribution,
                                              covid_stats$max_contribution) ~ "RSV",
            covid_stats$max_contribution == max(flu_stats$max_contribution,
                                                rsv_stats$max_contribution,
                                                covid_stats$max_contribution) ~ "COVID",
            TRUE ~ "Tie"
          )
        )

        model_results_table = bind_rows(model_results_table, summary_row)
        all_signals = bind_rows(all_signals, signals_tmp)

      }, error = function(e) {
        warning(paste("Error processing", geo, "in", current_season, ":", e$message))
      })

      cat("Completed:", geo, "\n")
    }
  }

  ## 4. Create Visualizations --------------------------------------------------

  # Create the 50-state facet plot
  facet_plot = all_signals %>%
    ggplot(aes(x = week_end)) +
    geom_line(aes(y = flu_times_coef, color = 'Flu'), linewidth = 0.65) +
    geom_line(aes(y = rsv_times_coef, color = 'RSV'), linewidth = 0.65) +
    geom_line(aes(y = covid_times_coef, color = "COVID-19"), linewidth = 0.65) +
    geom_line(aes(y = ili_rescaled, color = 'ILI'), linewidth = 0.65) +
    {if(model_type == "ridge_with_intercept")
      geom_line(aes(y = beta_intercept, color = 'Intercept'), linetype = "dashed", linewidth = 0.65)
    } +
    facet_wrap(~ state, scales = 'free', ncol = 10) +
    scale_color_manual(values = c(
      'ILI' = '#E2DBC9',
      'RSV' = '#183A5A',
      'Flu' = '#C34129',
      'COVID-19' = "#EFB75B",
      'Intercept' = 'grey30'
    ),
    breaks = if(model_type == "ridge_with_intercept") {
      c('ILI', 'Flu', 'RSV', 'COVID-19', 'Intercept')
    } else {
      c('ILI', 'Flu', 'RSV', 'COVID-19')
    }) +
    labs(
      x = '',
      y = ''
      # y = 'Rescaled Volume / Contribution'
      # title = paste("ILI Time Series Analysis:", str_replace_all(str_to_title(model_type), "_", " "))
    ) +
    scale_x_date(date_labels = "%b\n%y", date_breaks = "2 months") +
    envalysis::theme_publish() +
    theme(
      axis.text.x = element_text(size = 5),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.y = element_blank(),
      legend.title = element_blank(),
      legend.position = 'right',
      strip.text = element_text(size = 9),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    )

  # Create coefficient comparison plot
  coef_plot = model_results_table %>%
    select(state, season, beta_intercept, beta_flu, beta_rsv, beta_covid) %>%
    pivot_longer(cols = starts_with("beta_"), names_to = "coefficient", values_to = "value") %>%
    mutate(coefficient = case_when(
      coefficient == "beta_intercept" ~ "Intercept",
      coefficient == "beta_flu" ~ "Flu",
      coefficient == "beta_rsv" ~ "RSV",
      coefficient == "beta_covid" ~ "COVID"
    )) %>%
    filter(!(coefficient == "Intercept" & model_type != "ridge_with_intercept")) %>%
    ggplot(aes(x = coefficient, y = value, fill = coefficient)) +
    geom_boxplot() +
    scale_fill_manual(values =
                      c('ILI' = '#E2DBC9',
                      'RSV' = '#183A5A',
                      'Flu' = '#C34129',
                      'COVID-19' = "#EFB75B",
                      'Intercept' = 'grey30')) +
    facet_wrap(~season) +
    labs(
      title = paste("Distribution of Model Coefficients:", str_replace_all(str_to_title(model_type), "_", " ")),
      x = "Coefficient Type",
      y = "Coefficient Value"
    ) +
    envalysis::theme_publish() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
    )

  ## 5. Save Outputs -----------------------------------------------------------

  if (save_outputs) {
    # Save data
    write_csv(model_results_table, paste0("data/model_results_", output_prefix, ".csv"))
    write_csv(all_signals, paste0("data/time_series_", output_prefix, ".csv"))

    # Save plots
    ggsave(paste0("figures/facet_plot_", output_prefix, ".png"), facet_plot,
           width = 16, height = 10, units = "in", bg = "white")
    ggsave(paste0("figures/coefficient_plot_", output_prefix, ".png"), coef_plot,
           width = 10, height = 6, units = "in", bg = "white")

    cat("Outputs saved with prefix:", output_prefix, "\n")
  }

  ## 6. Print Summary ----------------------------------------------------------

  # cat("\n=== ANALYSIS SUMMARY ===\n")
  # cat("Model Type:", str_replace_all(str_to_title(model_type), "_", " "), "\n")
  # cat("Seasons Analyzed:", paste(seasons, collapse = ", "), "\n")
  # cat("Total state-seasons:", nrow(model_results_table), "\n")
  # cat("Mean intercept:", round(mean(model_results_table$beta_intercept), 3), "\n")
  # cat("Mean flu coefficient:", round(mean(model_results_table$beta_flu), 3), "\n")
  # cat("Mean RSV coefficient:", round(mean(model_results_table$beta_rsv), 3), "\n")
  # cat("Mean COVID coefficient:", round(mean(model_results_table$beta_covid), 3), "\n")
  #
  # cat("\nDominant pathogen distribution:\n")
  # print(table(model_results_table$dominant_pathogen))

  ## 7. Return Results ---------------------------------------------------------

  return(list(
    model_results = model_results_table,
    time_series_data = all_signals,
    facet_plot = facet_plot,
    coefficient_plot = coef_plot,
    model_type = model_type
  ))
}

## 8. Usage ----------------------------------------------------------

# Run for all 3 seasons combined to get the combined coef plots ----
results_all_original = run_ili_analysis(
  model_type = "ridge_without_intercept",  # or your preferred option
  seasons = c("22-23", "23-24", "24-25"),
  save_outputs = FALSE
)
# results_all_original$coefficient_plot
# ggsave(filename = "figures/manuscript_figures/supplementary/coef_plot_ridge_no_intercept.png",
#        width = 12,
#        height = 10,
#        units = "in",
#        bg = "white")

tbl_ridge_no_intercept = results_all_original$model_results

results_with_intercept = run_ili_analysis(
  model_type = "ridge_with_intercept",
  seasons = c("22-23", "23-24", "24-25"),
  save_outputs = FALSE
)
# results_with_intercept$coefficient_plot
# ggsave(filename = "figures/manuscript_figures/supplementary/coef_plot_ridge_with_intercept.png",
#        width = 12,
#        height = 10,
#        units = "in",
#        bg = "white")

results_linear = run_ili_analysis(
  model_type = "linear_no_regularization",
  seasons = c("22-23", "23-24", "24-25"),
  save_outputs = FALSE
)
# results_linear$coefficient_plot
# ggsave(filename = "figures/manuscript_figures/supplementary/coef_plot_linear.png",
#        width = 12,
#        height = 10,
#        units = "in",
#        bg = "white")


# Run original ridge model without intercept to get the season-by-season facet plots ----
results_original_22_23 = run_ili_analysis(
  model_type = "ridge_without_intercept",
  seasons = c("22-23"),
  save_outputs = TRUE,
  output_prefix = "original_ridge"
)
results_original_22_23$facet_plot # Save the results of these separately
ggsave("figures/modular_output/facet_plot_original_22_23.png",
       width = 16,
       height = 10,
       bg = "white")

results_original_23_24 = run_ili_analysis(
  model_type = "ridge_without_intercept",
  seasons = c("23-24"),
  save_outputs = TRUE,
  output_prefix = "original_ridge"
)
results_original_23_24$facet_plot # Save the results of these separately
ggsave("figures/modular_output/facet_plot_original_23_24.png",
       width = 16,
       height = 10,
       bg = "white")

results_original_24_25 = run_ili_analysis(
  model_type = "ridge_without_intercept",
  seasons = c("24-25"),
  save_outputs = TRUE,
  output_prefix = "original_ridge"
)
results_original_24_25$facet_plot
ggsave("figures/modular_output/facet_plot_original_24_25.png",
       width = 16,
       height = 10,
       bg = "white")

# Run ridge model with intercept for seasonal facet plots ----
## 22-23
results_with_intercept_22_23 = run_ili_analysis(
  model_type = "ridge_with_intercept",
  seasons = c("22-23"),
  save_outputs = FALSE,
  output_prefix = "ridge_with_intercept"
)
results_with_intercept_22_23$facet_plot
ggsave("figures/modular_output/facet_plot_intercept_22_23.png",
       width = 16,
       height = 10,
       bg = "white")

# 23-24
results_with_intercept_23_24 = run_ili_analysis(
  model_type = "ridge_with_intercept",
  seasons = c("23-24"),
  save_outputs = FALSE
)
results_with_intercept_23_24$facet_plot
ggsave("figures/modular_output/facet_plot_intercept_23_24.png",
       width = 16,
       height = 10,
       bg = "white")

# 24-25
results_with_intercept_24_25 = run_ili_analysis(
  model_type = "ridge_with_intercept",
  seasons = c("24-25"),
  save_outputs = FALSE
)
results_with_intercept_24_25$facet_plot
ggsave("figures/modular_output/facet_plot_intercept_24_25.png",
       width = 16,
       height = 10,
       bg = "white")

# Example 3: Run linear model without regularization
# 22-23
results_linear_22_23 = run_ili_analysis(
  model_type = "linear_no_regularization",
  seasons = c("22-23"),
  save_outputs = FALSE
)
results_linear_22_23$facet_plot
ggsave("figures/modular_output/facet_plot_linear_22_23.png",
       width = 16,
       height = 10,
       bg = "white")

results_linear_23_24 = run_ili_analysis(
  model_type = "linear_no_regularization",
  seasons = c("23-24"),
  save_outputs = TRUE,
  output_prefix = "linear_no_reg"
)
results_linear_23_24$facet_plot
ggsave("figures/modular_output/facet_plot_linear_23_24.png",
       width = 16,
       height = 10,
       bg = "white")

results_linear_24_25 = run_ili_analysis(
  model_type = "linear_no_regularization",
  seasons = c("24-25"),
  save_outputs = FALSE
)
results_linear_24_25$facet_plot
ggsave("figures/modular_output/facet_plot_linear_24_25.png",
       width = 16,
       height = 10,
       bg = "white")

results_all_ro

# Save the model outputs for building tables ----
save(results_all_original, file = "data/results_ridge_no_intercept.Rdata")
save(results_with_intercept, file = "data/results_ridge_with_intercept.Rdata")
save(results_linear, file = "data/results_linear.Rdata")
