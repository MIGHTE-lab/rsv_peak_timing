## -----------------------------------------------------------------------------
## Script name: fig1_plotting.R
##
## Purpose of script: Plotting code for Figure 1 of RSV Timing manuscript
##
## Author: George Dewey
##
## Date Created: 2025-03-26
##
## Last Updated: 2025-04-24
## Update notes: Updated for new data pull and after correction of models
## -----------------------------------------------------------------------------


# Setup -------------------------------------------------------------------

## Load packages
library(tidyverse)
library(lubridate)
library(patchwork)

## Load data
nssp_all_years = read_csv("data/processed/nssp_all_years_041525.csv")
dat_cor_fig1 = read_csv("data/data_cor_fig1.csv")
nssp_fig1 = nssp_all_years %>%
  filter(state %in% c("Maryland", "Texas", "New York"))

# Get the correlations for these states in each of the three seasons
dat_MD = dat_cor_fig1 %>% filter(state == "Maryland")
dat_TX = dat_cor_fig1 %>% filter(state == "Texas")
dat_NY = dat_cor_fig1 %>% filter(state == "New York")

# Function to compute correlations
compute_state_correlations = function(season_data, chosen_season) {
  cor_matrix = season_data %>%
    filter(season == chosen_season) %>%
    rename(rsv_lag0 = rsv_rescaled) %>%
    select(state, flu_rescaled, rsv_lag0, rsv_lag1:rsv_lag5) %>%  # Select necessary columns
    group_by(state) %>%
    summarise(across(rsv_lag0:rsv_lag5,
                     ~cor(flu_rescaled, ., use = "pairwise.complete.obs"))) %>%  # Compute correlations with flu_rescaled
    pivot_longer(cols = rsv_lag0:rsv_lag5, names_to = "RSV_Variable", values_to = "Correlation")

  return(cor_matrix)
}

cor_mat_MD_22_23 = compute_state_correlations(dat_MD, "22-23")
cor_mat_MD_23_24 = compute_state_correlations(dat_MD, "23-24")
cor_mat_MD_24_25 = compute_state_correlations(dat_MD, "24-25")

cor_mat_TX_22_23 = compute_state_correlations(dat_TX, "22-23")
cor_mat_TX_23_24 = compute_state_correlations(dat_TX, "23-24")
cor_mat_TX_24_25 = compute_state_correlations(dat_TX, "24-25")

cor_mat_NY_22_23 = compute_state_correlations(dat_NY, "22-23")
cor_mat_NY_23_24 = compute_state_correlations(dat_NY, "23-24")
cor_mat_NY_24_25 = compute_state_correlations(dat_NY, "24-25")

# Make a correlation matrix for each season to correspond to the columns in the
# trajectory plots

## 22-23
cor_22_23 = bind_rows(cor_mat_MD_22_23, cor_mat_TX_22_23, cor_mat_NY_22_23)

cor_22_23 = cor_22_23 %>%
  mutate(
    state = factor(state, levels = c("Texas", "New York", "Maryland")),
    RSV_Variable = factor(RSV_Variable, levels = c("rsv_lag5", "rsv_lag4", "rsv_lag3", "rsv_lag2", "rsv_lag1", "rsv_lag0"))
  )

# Create the heatmap
h1 = ggplot(cor_22_23, aes(x = RSV_Variable, y = state, fill = Correlation)) +
  geom_tile() +
  geom_text(aes(label = round(Correlation, 3)), color = "black", size = 3.5) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_x_discrete(labels =c ("Lag 5", "Lag 4", "Lag 3", "Lag 2", "Lag 1", "Lag 0")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        axis.title = element_blank())

## 23-24
cor_23_24 = bind_rows(cor_mat_MD_23_24, cor_mat_TX_23_24, cor_mat_NY_23_24)

cor_23_24 = cor_23_24 %>%
  mutate(
    state = factor(state, levels = c("Texas", "New York", "Maryland")),
    RSV_Variable = factor(RSV_Variable, levels = c("rsv_lag5", "rsv_lag4", "rsv_lag3", "rsv_lag2", "rsv_lag1", "rsv_lag0"))
  )

# Create the heatmap
h2 = ggplot(cor_23_24, aes(x = RSV_Variable, y = state, fill = Correlation)) +
  geom_tile() +
  geom_text(aes(label = round(Correlation, 3)), color = "black", size = 3.5) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_x_discrete(labels =c ("Lag 5", "Lag 4", "Lag 3", "Lag 2", "Lag 1", "Lag 0")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        axis.title = element_blank())

## 24-25
cor_24_25 = bind_rows(cor_mat_MD_24_25, cor_mat_TX_24_25, cor_mat_NY_24_25)

cor_24_25 = cor_24_25 %>%
  mutate(
    state = factor(state, levels = c("Texas", "New York", "Maryland")),
    RSV_Variable = factor(RSV_Variable, levels = c("rsv_lag5", "rsv_lag4", "rsv_lag3", "rsv_lag2", "rsv_lag1", "rsv_lag0"))
  )

# Create the heatmap
h3 = ggplot(cor_24_25, aes(x = RSV_Variable, y = state, fill = Correlation)) +
  geom_tile() +
  geom_text(aes(label = round(Correlation, 3)), color = "black", size = 3.5) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_x_discrete(labels =c ("Lag 5", "Lag 4", "Lag 3", "Lag 2", "Lag 1", "Lag 0")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        axis.title = element_blank())


# Then select only the states of interest (MD, NY, TX)
p1 = nssp_all_years %>%
  filter(state %in% c("Maryland", "Texas", "New York")) %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = flu_times_coef, color = 'Flu'), linewidth = 0.85) +
  geom_line(aes(y = rsv_times_coef, color = 'RSV'), linewidth = 0.85) +
  geom_line(aes(y = covid_times_coef, color = "COVID-19"), linewidth = 0.85) +
  geom_line(aes(y = ili_rescaled, color = 'ILI'), linewidth = 0.85) +
  facet_grid(state ~ season, scales = 'free') +
  scale_color_manual(
    limits = c("ILI", "RSV", "Flu", "COVID-19"),
    values = c(
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

p1 / (h1 + h2 + h3) + plot_layout(heights = c(0.8, 0.2))
# ggsave("figures/manuscript_figures/fig1/fig2_040825.png", width = 12, height = 12, units = "in", bg = "white")

## Making individual subfigures for presentation ----

## MA
nssp_all_years %>%
  filter(state %in% c("Texas"), season == "23-24") %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = flu_times_coef, color = 'Flu'), linewidth = 0.85) +
  geom_line(aes(y = rsv_times_coef, color = 'RSV'), linewidth = 0.85) +
  geom_line(aes(y = covid_times_coef, color = "COVID-19"), linewidth = 0.85) +
  geom_line(aes(y = ili_rescaled, color = 'ILI'), linewidth = 0.85) +
  facet_grid(state ~ season, scales = 'free') +
  scale_color_manual(
    limits = c("ILI", "RSV", "Flu", "COVID-19"),
    values = c(
      'ILI' = '#E2DBC9',
      'RSV' = '#183A5A',
      'Flu' = '#C34129',
      'COVID-19' = "#EFB75B"
    )) +
  labs(x = '', y = '') +
  scale_x_date(date_labels = "%b\n%y", date_breaks = "3 months") +
  envalysis::theme_publish() +
  theme(
    axis.text.x = element_text(size = 5),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    legend.title = element_blank(),
    legend.position = 'right'
  )

## MA
nssp_all_years %>%
  filter(state %in% c("Massachusetts"), season == ) %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = flu_times_coef, color = 'Flu'), linewidth = 0.85) +
  geom_line(aes(y = rsv_times_coef, color = 'RSV'), linewidth = 0.85) +
  geom_line(aes(y = covid_times_coef, color = "COVID-19"), linewidth = 0.85) +
  geom_line(aes(y = ili_rescaled, color = 'ILI'), linewidth = 0.85) +
  facet_grid(state ~ season, scales = 'free') +
  scale_color_manual(
    limits = c("ILI", "RSV", "Flu", "COVID-19"),
    values = c(
      'ILI' = '#E2DBC9',
      'RSV' = '#183A5A',
      'Flu' = '#C34129',
      'COVID-19' = "#EFB75B"
    )) +
  labs(x = '', y = '') +
  scale_x_date(date_labels = "%b\n%y", date_breaks = "3 months") +
  envalysis::theme_publish() +
  theme(
    axis.text.x = element_text(size = 5),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    legend.title = element_blank(),
    legend.position = 'right'
  )

## TX
nssp_all_years %>%
  filter(state %in% c("Texas")) %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = flu_times_coef, color = 'Flu'), linewidth = 0.85) +
  geom_line(aes(y = rsv_times_coef, color = 'RSV'), linewidth = 0.85) +
  geom_line(aes(y = covid_times_coef, color = "COVID-19"), linewidth = 0.85) +
  geom_line(aes(y = ili_rescaled, color = 'ILI'), linewidth = 0.85) +
  facet_grid(state ~ season, scales = 'free') +
  scale_color_manual(
    limits = c("ILI", "RSV", "Flu", "COVID-19"),
    values = c(
      'ILI' = '#E2DBC9',
      'RSV' = '#183A5A',
      'Flu' = '#C34129',
      'COVID-19' = "#EFB75B"
    )) +
  labs(x = '', y = '') +
  scale_x_date(date_labels = "%b\n%y", date_breaks = "3 months") +
  envalysis::theme_publish() +
  theme(
    axis.text.x = element_text(size = 5),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    legend.title = element_blank(),
    legend.position = 'right'
  )

## MD
nssp_all_years %>%
  filter(state %in% c("Maryland")) %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = flu_times_coef, color = 'Flu'), linewidth = 0.85) +
  geom_line(aes(y = rsv_times_coef, color = 'RSV'), linewidth = 0.85) +
  geom_line(aes(y = covid_times_coef, color = "COVID-19"), linewidth = 0.85) +
  geom_line(aes(y = ili_rescaled, color = 'ILI'), linewidth = 0.85) +
  facet_grid(state ~ season, scales = 'free') +
  scale_color_manual(
    limits = c("ILI", "RSV", "Flu", "COVID-19"),
    values = c(
      'ILI' = '#E2DBC9',
      'RSV' = '#183A5A',
      'Flu' = '#C34129',
      'COVID-19' = "#EFB75B"
    )) +
  labs(x = '', y = '') +
  scale_x_date(date_labels = "%b\n%y", date_breaks = "3 months") +
  envalysis::theme_publish() +
  theme(
    axis.text.x = element_text(size = 5),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    legend.title = element_blank(),
    legend.position = 'right'
  )

p1_22_23

nssp_all_years %>%
  filter(season == "22-23", state %in% c("Maryland", "Texas", "New York")) %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = flu_times_coef, color = 'Flu'), linewidth = 0.85) +
  geom_line(aes(y = rsv_times_coef, color = 'RSV'), linewidth = 0.85) +
  geom_line(aes(y = covid_times_coef, color = "COVID-19"), linewidth = 0.85) +
  geom_line(aes(y = ili_rescaled, color = 'ILI'), linewidth = 0.85) +
  facet_grid(state ~ season, scales = 'free') +
  scale_color_manual(
    limits = c("ILI", "RSV", "Flu", "COVID-19"),
    values = c(
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

p1

