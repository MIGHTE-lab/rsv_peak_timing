## -----------------------------------------------------------------------------
## Script name: lagged_correlation_timing.R
##
## Purpose of script: Check lagged correlations between flu and RSV in NSSP data
##
## Author: George Dewey
##
## Date Created: 2025-03-11
##
## Last Updated: 2025-03-20
##
## Version Update Notes: based on narrative, change correlation comparison
## to using flu and lagged rsv
## -----------------------------------------------------------------------------

# 1. Setup ----
# Load packages
library(tidyverse)
library(lubridate)
library(patchwork)
library(purrr)

# Load data
dat = read_csv("data/nssp_all_years_040725.csv")

dat_cor = dat %>% select(week_end, state, season, contains("rescaled")) %>%
  mutate(rsv_lag1 = lag(rsv_rescaled, 1),
         rsv_lag2 = lag(rsv_rescaled, 2),
         rsv_lag3 = lag(rsv_rescaled, 3),
         rsv_lag4 = lag(rsv_rescaled, 4),
         rsv_lag5 = lag(rsv_rescaled, 5),
         cov_lag1 = lag(covid_rescaled, 1),
         cov_lag2 = lag(covid_rescaled, 2),
         cov_lag3 = lag(covid_rescaled, 3),
         cov_lag4 = lag(covid_rescaled, 4),
         cov_lag5 = lag(covid_rescaled, 5),
         cov_lead1 = lead(covid_rescaled, 1),
         cov_lead2 = lead(covid_rescaled, 2),
         cov_lead3 = lead(covid_rescaled, 3),
         cov_lead4 = lead(covid_rescaled, 4),
         cov_lead5 = lead(covid_rescaled, 5))

dat_cor_fig1 = dat_cor %>%
  filter(state %in% c("Maryland", "Texas", "New York"))

write_csv(dat_cor_fig1, file = "data/data_cor_fig1.csv")

compute_correlations = function(season_data) {
  cor_matrix = season_data %>%
    rename(rsv_lag0 = rsv_rescaled) %>%
    select(state, flu_rescaled, rsv_lag0, rsv_lag1:rsv_lag5) %>%
    group_by(state) %>%
    summarise(across(rsv_lag0:rsv_lag5,
                     ~cor(flu_rescaled, .,
                          use = "pairwise.complete.obs"))) %>%
    pivot_longer(cols = rsv_lag0:rsv_lag5, names_to = "RSV_Variable",
                 values_to = "Correlation")

  return(cor_matrix)
}

nssp_22_23 = dat_cor %>% filter(season == "22-23")
nssp_23_24 = dat_cor %>% filter(season == "23-24")
nssp_24_25 = dat_cor %>% filter(season == "24-25")

cor_list_22_23 = compute_correlations(nssp_22_23)
cor_list_23_24 = compute_correlations(nssp_23_24)
cor_list_24_25 = compute_correlations(nssp_24_25)

lag_levels =
  c("rsv_lag5", "rsv_lag4", "rsv_lag3", "rsv_lag2", "rsv_lag1", "rsv_lag0")
rsv_labels = c("Lag 5", "Lag 4", "Lag 3", "Lag 2", "Lag 1", "Lag 0")

# 2. Correlation heatmaps ----
# Flu vs RSV
## 22-23 season ----
# First, identify the optimal (strongest) lag for each state
optimal_lags_22_23 = cor_list_22_23 %>%
  group_by(state) %>%
  slice_max(order_by = Correlation, n = 1, with_ties = FALSE) %>%
  ungroup()

# Sort the optimal lags from earliest data to latest
optimal_lags_22_23 = optimal_lags_22_23 %>%
  mutate(RSV_Variable =
           factor(RSV_Variable,
                  levels = c("rsv_lag5", "rsv_lag4", "rsv_lag3",
                             "rsv_lag2", "rsv_lag1", "rsv_lag0"))) %>%
  arrange(RSV_Variable, desc(Correlation)) %>%
  mutate(state = factor(state, levels = unique(state)))

optimal_lags_22_23 = optimal_lags_22_23 %>%
  mutate(RSV_Variable = factor(RSV_Variable, levels = lag_levels)) %>%
  arrange(RSV_Variable, desc(Correlation))

# Extract state order from optimal_lags
state_order_22_23 = optimal_lags_22_23 %>%
  pull(state)

state_order_22_23 = rev(state_order_22_23)

# Update correlation lists
cor_list_22_23 = cor_list_22_23 %>%
  mutate(state = factor(state, levels = state_order_22_23),
         RSV_Variable = factor(RSV_Variable, levels = lag_levels))

# Now your ggplot call will use the ordering from optimal_lags
ggplot(cor_list_22_23, aes(x = RSV_Variable, y = state, fill = Correlation)) +
  geom_tile() +
  geom_text(aes(label = round(Correlation, 2)), color = "black", size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_x_discrete(labels = rsv_labels) +
  theme_minimal() +
  labs(title = "Correlations of Flu and Lagged RSV Signals (22-23)",
       x = "RSV Lag", y = "State") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file = "figures/manuscript_figures/correlation_matrices/cor_matrix_flu_rsv_22_23_040725.png",
       height = 9, width = 10, units = "in", bg = "white")

## 23-24 season ----
# First, identify the optimal (strongest) lag for each state
optimal_lags_23_24 = cor_list_23_24 %>%
  group_by(state) %>%
  slice_max(order_by = Correlation, n = 1, with_ties = FALSE) %>%
  ungroup()

# Sort the optimal lags from earliest data to latest
optimal_lags_23_24 = optimal_lags_23_24 %>%
  mutate(RSV_Variable =
           factor(RSV_Variable,
                  levels = c("rsv_lag5", "rsv_lag4", "rsv_lag3", "rsv_lag2",
                             "rsv_lag1", "rsv_lag0"))) %>%
  arrange(RSV_Variable, desc(Correlation)) %>%
  mutate(state = factor(state, levels = unique(state)))
lag_levels = c("rsv_lag5", "rsv_lag4", "rsv_lag3", "rsv_lag2", "rsv_lag1",
               "rsv_lag0")

optimal_lags_23_24 = optimal_lags_23_24 %>%
  mutate(RSV_Variable = factor(RSV_Variable, levels = lag_levels)) %>%
  arrange(RSV_Variable, desc(Correlation))

# Extract state order from optimal_lags
state_order_23_24 = optimal_lags_23_24 %>%
  pull(state)

state_order_23_24 = rev(state_order_23_24)

# Update cor_23_24
cor_list_23_24 = cor_list_23_24 %>%
  mutate(state = factor(state, levels = state_order_23_24),
         RSV_Variable = factor(RSV_Variable, levels = lag_levels))

# Now your ggplot call will use the ordering from optimal_lags
ggplot(cor_list_23_24, aes(x = RSV_Variable, y = state, fill = Correlation)) +
  geom_tile() +
  geom_text(aes(label = round(Correlation, 2)), color = "black", size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_x_discrete(labels = rsv_labels) +
  theme_minimal() +
  labs(title = "Correlations of Flu and Lagged RSV Signals (23-24)",
       x = "RSV Lag", y = "State") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file = "figures/manuscript_figures/correlation_matrices/cor_matrix_flu_rsv_23_24_040725.png",
       height = 9, width = 10, units = "in", bg = "white")

## 24-25 season ----
# First, identify the optimal (strongest) lag for each state
optimal_lags_24_25 = cor_list_24_25 %>%
  group_by(state) %>%
  slice_max(order_by = Correlation, n = 1, with_ties = FALSE) %>%
  ungroup()

# Sort the optimal lags from earliest data to latest
optimal_lags_24_25 = optimal_lags_24_25 %>%
  mutate(RSV_Variable =
           factor(RSV_Variable,
                  levels = c("rsv_lag5", "rsv_lag4", "rsv_lag3",
                             "rsv_lag2", "rsv_lag1", "rsv_lag0"))) %>%
  arrange(RSV_Variable, desc(Correlation)) %>%
  mutate(state = factor(state, levels = unique(state)))

optimal_lags_24_25 = optimal_lags_24_25 %>%
  mutate(RSV_Variable = factor(RSV_Variable, levels = lag_levels)) %>%
  arrange(RSV_Variable, desc(Correlation))

# Extract state order from optimal_lags
state_order_24_25 = optimal_lags_24_25 %>%
  pull(state)

state_order_24_25 = rev(state_order_24_25)

# Update data data to set the factor levels for both state and RSV_Variable
cor_list_24_25 = cor_list_24_25 %>%
  mutate(state = factor(state, levels = state_order_24_25),
         RSV_Variable = factor(RSV_Variable, levels = lag_levels))

# Now your ggplot call will use the ordering from optimal_lags
ggplot(cor_list_24_25, aes(x = RSV_Variable, y = state, fill = Correlation)) +
  geom_tile() +
  geom_text(aes(label = round(Correlation, 2)), color = "black", size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_x_discrete(labels = rsv_labels) +
  theme_minimal() +
  labs(title = "Correlations of Flu and Lagged RSV Signals (24-25)",
       x = "RSV Lag", y = "State") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file = "figures/manuscript_figures/correlation_matrices/cor_matrix_flu_rsv_24_25_040725.png",
       height = 9, width = 10, units = "in", bg = "white")

## 2.1.4. Version 2 for 23-24 season ----
# Initialize state order using rsv_lag0
ordered_states = cor_list_23_24 %>%
  filter(RSV_Variable == "rsv_lag0") %>%
  arrange(Correlation) %>%
  pull(state)

# Iteratively reorder the states for the remaining columns
for (rsv_var in rsv_order[-1]) {
  new_order = cor_list_23_24 %>%
    filter(RSV_Variable == rsv_var) %>%
    arrange(Correlation) %>%
    pull(state)

  # Preserve previous order while inserting newly ranked states
  ordered_states = c(intersect(ordered_states, new_order),
                     setdiff(new_order, ordered_states))
}

# Apply the final order to the dataset
cor_list_23_24 = cor_list_23_24 %>%
  mutate(state = factor(state, levels = ordered_states),
         RSV_Variable = factor(RSV_Variable, levels = rsv_order))
rsv_labels = c("Lag 5", "Lag 4", "Lag 3", "Lag 2", "Lag 1", "Lag 0")

# Generate heatmap with updated x-axis labels
ggplot(cor_list_23_24, aes(x = RSV_Variable, y = state, fill = Correlation)) +
  geom_tile() +
  geom_text(aes(label = round(Correlation, 2)), color = "black", size = 3) +
  scale_fill_gradient2(low = "blue",
                       mid = "white",
                       high = "red",
                       midpoint = 0) +
  scale_x_discrete(labels = rsv_labels) +  # Update x-axis labels
  theme_minimal() +
  labs(title = "Correlation: Flu vs RSV (23-24)",
       x = "RSV Lag", y = "State") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file = "figures/cor_matrix_rsv_lag_23_24.png",
       height = 9, width = 7, units = "in", bg = "white")

# 2.2. Correlations between Flu and COVID ----
# Since we are unsure about the seasonality of COVID, let's check both leads and lags
compute_correlations_flu_cov = function(season_data) {
  cor_matrix = season_data %>%
    rename(cov_lag0 = covid_rescaled) %>%
    select(state, flu_rescaled, cov_lag0, cov_lag1:cov_lag5,
           cov_lead1:cov_lead5) %>%
    group_by(state) %>%
    summarise(across(cov_lag0:cov_lead5,
                     ~cor(flu_rescaled, .,
                          use = "pairwise.complete.obs"))) %>%
    pivot_longer(cols = cov_lag0:cov_lead5, names_to = "COV_Variable",
                 values_to = "Correlation")

  return(cor_matrix)
}

# Compute the correlations between flu and lagged COVID
cor_cov_flu_22_23 = compute_correlations_flu_cov(nssp_22_23)
cor_cov_flu_23_24 = compute_correlations_flu_cov(nssp_23_24)
cor_cov_flu_24_25 = compute_correlations_flu_cov(nssp_24_25)

cov_order = c("cov_lag5", "cov_lag4", "cov_lag3", "cov_lag2", "cov_lag1",
              "cov_lag0", "cov_lead1", "cov_lead2", "cov_lead3", "cov_lead4",
              "cov_lead5")
cov_labels = c("Lag 5", "Lag 4", "Lag 3", "Lag 2", "Lag 1", "Lag 0", "Lead 1",
               "Lead 2", "Lead 3", "Lead 4", "Lead 5")

# Ordered heatmaps for COVID and Flu
# 2.2.1. 22-23 season -----
optimal_lag_flu_COVID_22_23 = cor_cov_flu_22_23 %>%
  group_by(state) %>%
  slice_max(order_by = Correlation, n = 1, with_ties = FALSE) %>%
  ungroup()

# Sort the optimal lags from earliest data to latest
optimal_lag_flu_COVID_22_23 = optimal_lag_flu_COVID_22_23 %>%
  mutate(COV_Variable = factor(COV_Variable,
                               levels = cov_order)) %>%
  arrange(COV_Variable, desc(Correlation)) %>%
  mutate(state = factor(state, levels = unique(state)))

optimal_lag_flu_COVID_22_23

# Extract state order from optimal_lags
state_order_flu_COVID_22_23 = optimal_lag_flu_COVID_22_23 %>%
  pull(state)

state_order_flu_COVID_22_23

state_order_flu_COVID_22_23 = rev(state_order_flu_COVID_22_23)

# Update data data to set the factor levels for both state and RSV_Variable
cor_cov_flu_22_23 = cor_cov_flu_22_23 %>%
  mutate(state = factor(state, levels = state_order_flu_COVID_22_23),
         COV_Variable = factor(COV_Variable, levels = cov_order))


# Generate heatmap with updated x-axis labels
ggplot(cor_cov_flu_22_23, aes(x = COV_Variable, y = state, fill = Correlation)) +
  geom_tile() +
  geom_text(aes(label = round(Correlation, 2)), color = "black", size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_x_discrete(labels = cov_labels) +  # Update x-axis labels
  theme_minimal() +
  labs(title = "Correlation: Flu vs Shifted COVID Values (22-23)",
       x = "", y = "State") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file = "figures/manuscript_figures/correlation_matrices/cor_matrix_covid_flu_22_23_040725.png",
       height = 8, width = 14, units = "in", bg = "white")

## 2.2.2 23-24 season ----

optimal_lag_flu_COVID_23_24 = cor_cov_flu_23_24 %>%
  group_by(state) %>%
  slice_max(order_by = Correlation, n = 1, with_ties = FALSE) %>%
  ungroup()

# Sort the optimal lags from earliest data to latest
optimal_lag_flu_COVID_23_24 = optimal_lag_flu_COVID_23_24 %>%
  mutate(COV_Variable = factor(COV_Variable,
                               levels = cov_order)) %>%
  arrange(COV_Variable, desc(Correlation)) %>%
  mutate(state = factor(state, levels = unique(state)))

# Extract state order from optimal_lags
state_order_flu_COVID_23_24 = optimal_lag_flu_COVID_23_24 %>%
  pull(state)

state_order_flu_COVID_23_24 = rev(state_order_flu_COVID_23_24)

# Update data data to set the factor levels for both state and RSV_Variable
cor_cov_flu_23_24 = cor_cov_flu_23_24 %>%
  mutate(state = factor(state, levels = state_order_flu_COVID_23_24),
         COV_Variable = factor(COV_Variable, levels = cov_order))

# Generate heatmap with updated x-axis labels
ggplot(cor_cov_flu_23_24, aes(x = COV_Variable, y = state, fill = Correlation)) +
  geom_tile() +
  geom_text(aes(label = round(Correlation, 2)), color = "black", size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_x_discrete(labels = cov_labels) +  # Update x-axis labels
  theme_minimal() +
  labs(title = "Correlation: Flu vs Shifted COVID Values (23-24)",
       x = "", y = "State") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file = "figures/manuscript_figures/correlation_matrices/cor_matrix_covid_flu_23_24_040725.png",
       height = 8, width = 14, units = "in", bg = "white")

## 2.2.3 24-25 season ----
optimal_lag_flu_COVID_24_25 = cor_cov_flu_24_25 %>%
  group_by(state) %>%
  slice_max(order_by = Correlation, n = 1, with_ties = FALSE) %>%
  ungroup()

# Sort the optimal lags from earliest data to latest
optimal_lag_flu_COVID_24_25 = optimal_lag_flu_COVID_24_25 %>%
  mutate(COV_Variable = factor(COV_Variable,
                               levels = cov_order)) %>%
  arrange(COV_Variable, desc(Correlation)) %>%
  mutate(state = factor(state, levels = unique(state)))

# Extract state order from optimal_lags
state_order_flu_COVID_24_25 = optimal_lag_flu_COVID_24_25 %>%
  pull(state)

state_order_flu_COVID_24_25 = rev(state_order_flu_COVID_24_25)

# Update data data to set the factor levels for both state and RSV_Variable
cor_cov_flu_24_25 = cor_cov_flu_24_25 %>%
  mutate(state = factor(state, levels = state_order_flu_COVID_24_25),
         COV_Variable = factor(COV_Variable, levels = cov_order))

# Generate heatmap with updated x-axis labels
ggplot(cor_cov_flu_24_25, aes(x = COV_Variable, y = state, fill = Correlation)) +
  geom_tile() +
  geom_text(aes(label = round(Correlation, 2)), color = "black", size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_x_discrete(labels = cov_labels) +  # Update x-axis labels
  theme_minimal() +
  labs(title = "Correlation: Flu vs Shifted COVID Values (24-25)",
       x = "", y = "State") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file = "figures/manuscript_figures/correlation_matrices/cor_matrix_covid_flu_24_25_040725.png",
       height = 8, width = 14, units = "in", bg = "white")

# Sensitivity Analysis ----------------------------------------------------

# Here we perform a sensitivity analysis for the correlation by leaving out the
# double peak of the 24-25 season. Based on inspection of the 24-25 plots,
# the double peak of ILI during the current season occurred during mid-January
# in most states. Let's truncate the data after the first 2 weeks of January
# so we leave in the first peak of ILI but not the second.

nssp_24_25_no_second_peak = nssp_24_25 %>%
  filter(week_end <= as_date("2025-01-07"))

cor_list_24_25_trunc = compute_correlations(nssp_24_25_no_second_peak)

cor_list_24_25_trunc %>%
  mutate(Flu_Variable =
           ifelse(Flu_Variable == "flu_rescaled",
                  "flu_lead0",
                  Flu_Variable)) %>%
  ggplot(aes(x =
               factor(Flu_Variable,
                      levels = c("flu_lead0", "flu_lead1", "flu_lead2",
                                 "flu_lead3", "flu_lead4", "flu_lead5")),
                           y = state, fill = Correlation)) +
  geom_tile() +
  geom_text(aes(label = round(Correlation, 2)), color = "black", size = 3) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(title = "Correlation: RSV vs Flu Variables (24-25), Truncated Data",
        x = "Flu Variables (Lagged & Unlagged)", y = "State") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

cor_list_24_25 = cor_list_24_25 %>%
  mutate(Flu_Variable =
           ifelse(Flu_Variable == "flu_rescaled",
                  "flu_lead0",
                  Flu_Variable)) %>%
  rename(cor_no_trunc = Correlation)

cor_list_24_25_trunc %>%
  mutate(Flu_Variable = ifelse(Flu_Variable ==
                                 "flu_rescaled",
                               "flu_lead0",
                               Flu_Variable)) %>%
  left_join(cor_list_24_25, by = c("state", "Flu_Variable")) %>%
  rename(cor_trunc = Correlation) %>%
  ggplot(aes(x = cor_no_trunc, y = cor_trunc)) +
  geom_point() +
  labs(x = "Correlation (full data)", y = "Correlation (truncated data)") +
  envalysis::theme_publish()

cor_list_24_25_trunc %>%
    mutate(Flu_Variable = ifelse(Flu_Variable == "flu_rescaled",
                                 "flu_lead0", Flu_Variable)) %>%
    left_join(cor_list_24_25, by = c("state", "Flu_Variable")) %>%
    rename(cor_trunc = Correlation) %>%
    mutate(change = abs(cor_no_trunc - cor_trunc)) %>%
  ggplot(aes(x = change)) +
  geom_histogram(binwidth = 0.025) +
  envalysis::theme_publish() +
  labs(x = "Absolute Change in Correlation", y = "Count")

## Exporting data for EWS
nssp_24_25 = nssp_24_25 %>%
  select(week_end, state, season, ili_rescaled, flu_rescaled, rsv_rescaled, covid_rescaled)

dat_24_25 = dat %>%
  filter(season == "24-25")

write_csv(dat_24_25, file = "data/nssp_24_25.csv")
