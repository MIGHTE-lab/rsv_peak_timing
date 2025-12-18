## -----------------------------------------------------------------------------
## Script name: statewise_aggregate_figures.R
##
## Purpose of script: To build figures for the 2024-07-30 Epistorm meeting
##
## Author: George Dewey
##
## Date Created: 2024-07-16
##
## Last Updated: 2024-07-16
## -----------------------------------------------------------------------------

## Load packages
library(tidyverse)
library(lubridate)
library(firatheme)

signals_aggregate = read_csv(file = 'data/ground_truth/signals_aggregated.csv')

signals_aggregate

## Comparing Figures State-wise (for presentation)
### Ohio -----------------------------------------------------------------------
signals_aggregate %>%
  filter(geography == 'Ohio') %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = flu_times_coef, color = 'Flu'), linewidth = 0.75) +
  geom_line(aes(y = cov_times_coef, color = 'COVID'), linewidth = 0.75) +
  geom_line(aes(y = rsv_times_coef, color = 'RSV'), linewidth = 0.75) +
  geom_line(aes(y = ili, color = 'ILI'), linewidth = 0.75, alpha = 0.8) +
  geom_area(aes(y = ili), fill = '#d6ede9', alpha = 0.4) +
  scale_color_manual(values = c('ILI' = '#d6ede9',
                                'RSV' = '#20afd1',
                                'Flu' = '#cc3e00',
                                'COVID' = '#ffce00')) +
  scale_x_date(date_labels = '%b\n%Y',
               date_breaks = '4 months',
               limits = c(as_date('2021-11-01'), as_date('2024-06-29'))) +
  theme_minimal() +
  labs(x = '', y = '', title = 'Ohio') +
  theme(plot.title = element_text(size = 35, family = 'Barlow Semi Condensed'),
        axis.text.x = element_text(size = 24, family = 'Barlow Semi Condensed'),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 20, family = 'Barlow Semi Condensed'),
        legend.position = 'bottom')

### Tennessee -----------------------------------------------------------------------
signals_aggregate %>%
  filter(geography == 'Tennessee') %>%
  ggplot(aes(x = week_end)) +
  geom_area(aes(y = ili, fill = 'ILI'), alpha = 0.4) +
  geom_line(aes(y = flu_times_coef, color = 'Flu'), linewidth = 0.75, alpha = 1) +
  geom_line(aes(y = cov_times_coef, color = 'COVID'), linewidth = 0.75, alpha = 1) +
  geom_line(aes(y = rsv_times_coef, color = 'RSV'), linewidth = 0.75, alpha = 1) +
  geom_line(aes(y = ili, color = 'ILI'), linewidth = 0.75, alpha = 0) +
  scale_color_manual(values = c('ILI' = '#d6ede9',
                                'RSV' = '#20afd1',
                                'Flu' = '#cc3e00',
                                'COVID' = '#ffce00')) +
  scale_fill_manual(values = c('ILI' = '#d6ede9',
                               'RSV' = '#20afd1',
                               'Flu' = '#cc3e00',
                               'COVID' = '#ffce00'), guide = 'none') +
  scale_x_date(date_labels = '%b\n%Y',
               date_breaks = '4 months',
               limits = c(as_date('2021-11-01'), as_date('2024-06-29'))) +
  theme_fira() +
  labs(x = '', y = '', title = 'Tennessee') +
  theme(plot.title = element_text(size = 35, family = 'Barlow Semi Condensed'),
        axis.text.x = element_text(size = 24, family = 'Barlow Semi Condensed'),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 24, family = 'Barlow Semi Condensed'),
        legend.position = 'bottom') +
  guides(col = guide_legend(override.aes = list(alpha = 1)))

ggsave('data/figures/epistorm_pres_2024_07_30/TN_trend.png', width = 10, height = 7,
       unit = 'in', bg = 'white')

### Maryland -----------------------------------------------------------------------
signals_aggregate %>%
  filter(geography == 'Maryland') %>%
  ggplot(aes(x = week_end)) +
  geom_area(aes(y = ili, fill = 'ILI'), alpha = 0.4) +
  geom_line(aes(y = flu_times_coef, color = 'Flu'), linewidth = 0.75, alpha = 1) +
  geom_line(aes(y = cov_times_coef, color = 'COVID'), linewidth = 0.75, alpha = 1) +
  geom_line(aes(y = rsv_times_coef, color = 'RSV'), linewidth = 0.75, alpha = 1) +
  geom_line(aes(y = ili, color = 'ILI'), linewidth = 0.75, alpha = 0) +
  scale_color_manual(values = c('ILI' = '#d6ede9',
                                'RSV' = '#20afd1',
                                'Flu' = '#cc3e00',
                                'COVID' = '#ffce00')) +
  scale_fill_manual(values = c('ILI' = '#d6ede9',
                               'RSV' = '#20afd1',
                               'Flu' = '#cc3e00',
                               'COVID' = '#ffce00'), guide = 'none') +
  scale_x_date(date_labels = '%b\n%Y',
               date_breaks = '4 months',
               limits = c(as_date('2021-11-01'), as_date('2024-06-29'))) +
  theme_fira() +
  labs(x = '', y = '', title = 'Maryland') +
  theme(plot.title = element_text(size = 35, family = 'Barlow Semi Condensed'),
        axis.text.x = element_text(size = 24, family = 'Barlow Semi Condensed'),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 24, family = 'Barlow Semi Condensed'),
        legend.position = 'bottom') +
  guides(col = guide_legend(override.aes = list(alpha = 1)))

ggsave('data/figures/epistorm_pres_2024_07_30/MD_trend.png', width = 10, height = 7,
       unit = 'in', bg = 'white')

### Ohio -----------------------------------------------------------------------
signals_aggregate %>%
  filter(geography == 'Ohio') %>%
  ggplot(aes(x = week_end)) +
  geom_area(aes(y = ili, fill = 'ILI'), alpha = 0.4) +
  geom_line(aes(y = flu_times_coef, color = 'Flu'), linewidth = 0.75, alpha = 1) +
  geom_line(aes(y = cov_times_coef, color = 'COVID'), linewidth = 0.75, alpha = 1) +
  geom_line(aes(y = rsv_times_coef, color = 'RSV'), linewidth = 0.75, alpha = 1) +
  geom_line(aes(y = ili, color = 'ILI'), linewidth = 0.75, alpha = 0) +
  scale_color_manual(values = c('ILI' = '#d6ede9',
                                'RSV' = '#20afd1',
                                'Flu' = '#cc3e00',
                                'COVID' = '#ffce00')) +
  scale_fill_manual(values = c('ILI' = '#d6ede9',
                               'RSV' = '#20afd1',
                               'Flu' = '#cc3e00',
                               'COVID' = '#ffce00'), guide = 'none') +
  scale_x_date(date_labels = '%b\n%Y',
               date_breaks = '4 months',
               limits = c(as_date('2021-11-01'), as_date('2024-06-29'))) +
  theme_fira() +
  labs(x = '', y = '', title = 'Ohio') +
  theme(plot.title = element_text(size = 35, family = 'Barlow Semi Condensed'),
        axis.text.x = element_text(size = 24, family = 'Barlow Semi Condensed'),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 24, family = 'Barlow Semi Condensed'),
        legend.position = 'bottom') +
  guides(col = guide_legend(override.aes = list(alpha = 1)))

ggsave('data/figures/epistorm_pres_2024_07_30/OH_trend.png', width = 10, height = 7,
       unit = 'in', bg = 'white')

### Oregon -----------------------------------------------------------------------
signals_aggregate %>%
  filter(geography == 'Oregon') %>%
  ggplot(aes(x = week_end)) +
  geom_line(aes(y = flu_times_coef, color = 'Flu'), linewidth = 0.75, alpha = 0.8) +
  geom_line(aes(y = cov_times_coef, color = 'COVID'), linewidth = 0.75, alpha = 0.8) +
  geom_line(aes(y = rsv_times_coef, color = 'RSV'), linewidth = 0.75, alpha = 0.8) +
  geom_line(aes(y = ili, color = 'ILI'), linewidth = 0.75, alpha = 0.8) +
  scale_color_manual(values = c('ILI' = '#d6ede9',
                                'RSV' = '#20afd1',
                                'Flu' = '#cc3e00',
                                'COVID' = '#ffce00')) +
  scale_x_date(date_labels = '%b\n%Y',
               date_breaks = '4 months',
               limits = c(as_date('2021-11-01'), as_date('2024-06-29'))) +
  theme_minimal() +
  labs(x = '', y = '', title = 'Oregon') +
  theme(plot.title = element_text(size = 35, family = 'Barlow Semi Condensed'),
        axis.text.x = element_text(size = 24, family = 'Barlow Semi Condensed'),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(size = 20, family = 'Barlow Semi Condensed'),
        legend.position = 'bottom')
