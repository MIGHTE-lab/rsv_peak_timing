## -----------------------------------------------------------------------------
## Script name: lagged_correlation_timing.R
##
## Purpose of script: Check lagged correlations between flu and RSV in NSSP data
##
## Author: George Dewey
##
## Date Created: 2025-03-11
##
## Last Updated: 2025-03-13
## -----------------------------------------------------------------------------

# 1. Setup ----
# Load packages
library(tidyverse)
library(lubridate)
library(patchwork)
library(purrr)

# Load data
dat = read_csv("data/processed/nssp_all_years.csv")

dat_cor = dat %>% select(week_end, state, season, contains("rescaled")) %>%
  mutate(flu_lead1 = lead(flu_rescaled, 1),
         flu_lead2 = lead(flu_rescaled, 2),
         flu_lead3 = lead(flu_rescaled, 3),
         flu_lead4 = lead(flu_rescaled, 4),
         flu_lead5 = lead(flu_rescaled, 5))

compute_correlations = function(season_data) {
  cor_matrix = season_data %>%
    select(state, rsv_rescaled, flu_rescaled, flu_lead1:flu_lead5) %>%  # Select necessary columns
    group_by(state) %>%
    summarise(across(flu_rescaled:flu_lead5,
                     ~cor(rsv_rescaled, ., use = "pairwise.complete.obs"))) %>%  # Compute correlations with rsv_rescaled
    pivot_longer(cols = flu_rescaled:flu_lead5, names_to = "Flu_Variable", values_to = "Correlation")

  return(cor_matrix)
}

# List of seasons
seasons = c("22-23", "23-24", "24-25")

nssp_22_23 = dat_cor %>% filter(season == "22-23")
nssp_23_24 = dat_cor %>% filter(season == "23-24")
nssp_24_25 = dat_cor %>% filter(season == "24-25")

cor_list_22_23 = compute_correlations(nssp_22_23)
cor_list_23_24 = compute_correlations(nssp_23_24)
cor_list_24_25 = compute_correlations(nssp_24_25)

# 2. Correlation heatmaps ----

## 22-23 season
cor_list_22_23 = cor_list_22_23 %>%
  mutate(state = factor(state, levels = sort(unique(state), decreasing = TRUE)))

ggplot(cor_list_22_23, aes(x = factor(Flu_Variable, levels = c("flu_rescaled", "flu_lead1", "flu_lead2", "flu_lead3", "flu_lead4", "flu_lead5")),
                           y = state, fill = Correlation)) +
  geom_tile() +
  geom_text(aes(label = round(Correlation, 2)), color = "black", size = 3) +  # Add correlation values
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(title = "Correlation: RSV vs Flu Variables (22-23)",
       x = "Flu Variables (Lagged & Unlagged)", y = "State") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file = "figures/cor_matrix_flu_lead_22_23.png", width = 9, height = 10, units = "in", bg = "white")

## 23-24 season
cor_list_23_24 = cor_list_23_24 %>%
  mutate(state = factor(state, levels = sort(unique(state), decreasing = TRUE)))

ggplot(cor_list_23_24, aes(x = factor(Flu_Variable, levels = c("flu_rescaled", "flu_lead1", "flu_lead2", "flu_lead3", "flu_lead4", "flu_lead5")),
                           y = state, fill = Correlation)) +
  geom_tile() +
  geom_text(aes(label = round(Correlation, 2)), color = "black", size = 3) +  # Add correlation values
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(title = "Correlation: RSV vs Flu Variables (23-24)",
       x = "Flu Variables (Lagged & Unlagged)", y = "State") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file = "figures/cor_matrix_flu_lead_23_24.png", width = 9, height = 10, units = "in", bg = "white")

## 24-25 season
cor_list_24_25 = cor_list_24_25 %>%
  mutate(state = factor(state, levels = sort(unique(state), decreasing = TRUE)))
cor_list_24_25

ggplot(cor_list_24_25, aes(x = factor(Flu_Variable, levels = c("flu_lead0", "flu_lead1", "flu_lead2", "flu_lead3", "flu_lead4", "flu_lead5")),
                           y = state, fill = cor_no_trunc)) +
  geom_tile() +
  geom_text(aes(label = round(cor_no_trunc, 2)), color = "black", size = 3) +  # Add correlation values
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(title = "Correlation: RSV vs Flu Variables (24-25), Full Data",
       x = "Flu Variables (Lagged & Unlagged)", y = "State") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file = "figures/cor_matrix_flu_lead_24_25.png", width = 9, height = 10, units = "in", bg = "white")

dat_FL_22_23 = nssp_22_23 %>% filter(state == "Florida")
cor(dat_FL_22_23$rsv_rescaled, dat_FL_22_23$flu_rescaled)

ccf(dat$flu_rescaled, dat$rsv_rescaled, lag.max = 5)


# Sensitivity Analysis ----------------------------------------------------

# Here we perform a sensitivity analysis for the correlation by leaving out the
# double peak of the 24-25 season. Based on inspection of the 24-25 plots,
# the double peak of ILI during the current season occurred during mid-January
# in most states. Let's truncate the data after the first 2 weeks of January
# so we leave in the first peak of ILI but not the second.

nssp_24_25_no_second_peak = nssp_24_25 %>% filter(week_end <= as_date("2025-01-07"))

cor_list_24_25_trunc = compute_correlations(nssp_24_25_no_second_peak)

cor_list_24_25_trunc %>%
  mutate(Flu_Variable = ifelse(Flu_Variable == "flu_rescaled", "flu_lead0", Flu_Variable)) %>%
  ggplot(aes(x = factor(Flu_Variable, levels = c("flu_lead0", "flu_lead1", "flu_lead2", "flu_lead3", "flu_lead4", "flu_lead5")),
                           y = state, fill = Correlation)) +
  geom_tile() +
  geom_text(aes(label = round(Correlation, 2)), color = "black", size = 3) +  # Add correlation values
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(title = "Correlation: RSV vs Flu Variables (24-25), Truncated Data",
       x = "Flu Variables (Lagged & Unlagged)", y = "State") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

cor_list_24_25 = cor_list_24_25 %>%
  mutate(Flu_Variable = ifelse(Flu_Variable == "flu_rescaled", "flu_lead0", Flu_Variable)) %>%
  rename(cor_no_trunc = Correlation)

cor_list_24_25_trunc %>%
  mutate(Flu_Variable = ifelse(Flu_Variable == "flu_rescaled", "flu_lead0", Flu_Variable)) %>%
  left_join(cor_list_24_25, by = c("state", "Flu_Variable")) %>%
  rename(cor_trunc = Correlation) %>%
  ggplot(aes(x = cor_no_trunc, y = cor_trunc)) +
  geom_point() +
  labs(x = "Correlation (full data)", y = "Correlation (truncated data)") +
  envalysis::theme_publish()

cor_list_24_25_trunc %>%
    mutate(Flu_Variable = ifelse(Flu_Variable == "flu_rescaled", "flu_lead0", Flu_Variable)) %>%
    left_join(cor_list_24_25, by = c("state", "Flu_Variable")) %>%
    rename(cor_trunc = Correlation) %>%
    mutate(change = abs(cor_no_trunc - cor_trunc)) %>%
  ggplot(aes(x = change)) +
  geom_histogram(binwidth = 0.025) +
  envalysis::theme_publish() +
  labs(x = "Absolute Change in Correlation", y = "Count")
