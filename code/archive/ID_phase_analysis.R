## -----------------------------------------------------------------------------
## Script name: ID_phase_analysis.R
##
## Purpose of script: Formal phase analysis for RSV and influenza timeseries
##
## Author: George Dewey
##
## Date Created: 2025-03-20
##
## Last Updated: 2025-03-20
## -----------------------------------------------------------------------------

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

