## -----------------------------------------------------------------------------
## Script name: onset_peak_scatterplot.R
##
## Purpose of script: Creating scatterplot of onset and peak dates
##
## Author: George Dewey
##
## Date Created: 2025-04-23
##
## Last Updated: 2025-07-08
##
## Update note: Code cleanup
## -----------------------------------------------------------------------------

# 1. Setup ----------------------------------------------------------------

## Load packages
library(tidyverse)
library(lubridate)
library(patchwork)
library(ggrepel)

## Load data
flu_dat = read_csv('data/early_warning/flu_matched_peaks_onsets.csv')
rsv_dat = read_csv('data/early_warning/rsv_matched_peaks_onsets.csv')
cov_dat = read_csv('data/early_warning/covid_matched_peaks_onsets_2.csv')

## Rename columns
flu_dat = flu_dat %>% rename(
  flu_onset = start_date,
  flu_end = end_date,
  flu_peak = date_peak,
  state = name_location
)

rsv_dat = rsv_dat %>% rename(rsv_onset = onset_date,
                             rsv_end = end_date,
                             rsv_peak = date_peak)

cov_dat = cov_dat %>%
  rename(
    covid_onset = start_date,
    covid_end = end_date,
    covid_peak = date_peak,
    state = name_location
  ) %>%
  select(state, covid_onset, covid_end, covid_peak, season, time_diff) %>%
  mutate(time_diff_week = round(time_diff / 7, 1))

## Create the season-specific data subsets
flu_dat_23 = flu_dat %>% filter(season == "23-24")
flu_dat_24 = flu_dat %>% filter(season == "24-25")
rsv_dat_23 = rsv_dat %>% filter(season == "23-24")
rsv_dat_24 = rsv_dat %>% filter(season == "24-25")
cov_dat_23 = cov_dat %>% filter(season == "23-24")
cov_dat_24 = cov_dat %>% filter(season == "24-25")

## Helper function
state_to_abbrev = function(state_name) {
  state_abbrevs = c(
    "Alabama" = "AL",
    "Alaska" = "AK",
    "Arizona" = "AZ",
    "Arkansas" = "AR",
    "California" = "CA",
    "Colorado" = "CO",
    "Connecticut" = "CT",
    "Delaware" = "DE",
    "District of Columbia" = "DC",
    "Florida" = "FL",
    "Georgia" = "GA",
    "Hawaii" = "HI",
    "Idaho" = "ID",
    "Illinois" = "IL",
    "Indiana" = "IN",
    "Iowa" = "IA",
    "Kansas" = "KS",
    "Kentucky" = "KY",
    "Louisiana" = "LA",
    "Maine" = "ME",
    "Maryland" = "MD",
    "Massachusetts" = "MA",
    "Michigan" = "MI",
    "Minnesota" = "MN",
    "Mississippi" = "MS",
    "Missouri" = "MO",
    "Montana" = "MT",
    "Nebraska" = "NE",
    "Nevada" = "NV",
    "New Hampshire" = "NH",
    "New Jersey" = "NJ",
    "New Mexico" = "NM",
    "New York" = "NY",
    "North Carolina" = "NC",
    "North Dakota" = "ND",
    "Ohio" = "OH",
    "Oklahoma" = "OK",
    "Oregon" = "OR",
    "Pennsylvania" = "PA",
    "Rhode Island" = "RI",
    "South Carolina" = "SC",
    "South Dakota" = "SD",
    "Tennessee" = "TN",
    "Texas" = "TX",
    "Utah" = "UT",
    "Vermont" = "VT",
    "Virginia" = "VA",
    "Washington" = "WA",
    "West Virginia" = "WV",
    "Wisconsin" = "WI",
    "Wyoming" = "WY"
  )

  return(state_abbrevs[state_name])
}

# Checking to see if the conclusions are correct
flu_dat_23 %>%
  left_join(rsv_dat_23, by = c("state", "season")) %>%
  select(state, flu_onset, flu_peak, rsv_onset, rsv_peak) %>%
  mutate(time_diff =
           as.numeric(difftime(flu_peak, rsv_peak, units = "weeks")))

flu_dat_24 %>%
  left_join(rsv_dat_24, by = c("state", "season")) %>%
  select(state, flu_onset, flu_peak, rsv_onset, rsv_peak) %>%
  mutate(time_diff =
           as.numeric(difftime(flu_peak, rsv_peak, units = "weeks"))) %>%
  ggplot(aes(x = time_diff)) +
  geom_boxplot()

# Export these data to build other visualizations
write_csv(flu_dat_23, "data/processed/flu_onset_peak_times_23.csv")
write_csv(flu_dat_24, "data/processed/flu_onset_peak_times_24.csv")
write_csv(rsv_dat_23, "data/processed/rsv_onset_peak_times_23.csv")
write_csv(rsv_dat_24, "data/processed/rsv_onset_peak_times_24.csv")
write_csv(cov_dat_23, "data/processed/cov_onset_peak_times_23.csv")
write_csv(cov_dat_24, "data/processed/cov_onset_peak_times_24.csv")

# 2. Main plotting function -----------------------------------------------
# Create a function to plot flu dates vs comparison dates (RSV or COVID)
create_comparison_plot = function(flu_data,
                                  comparison_data,
                                  comparison_type,
                                  season,
                                  type = "onset",
                                  base_color = "darkslateblue") {
  # Validate the type parameter
  if (!type %in% c("onset", "peak")) {
    stop("Type must be either 'onset' or 'peak'")
  }
  # Validate the comparison_type parameter
  if (!comparison_type %in% c("rsv", "covid")) {
    stop("comparison_type must be either 'rsv' or 'covid'")
  }
  # Set column names based on type
  flu_col = paste0("flu_", type)
  comparison_col = paste0(comparison_type, "_", type)
  # Get the proper title capitalization
  type_title = paste0(toupper(substr(type, 1, 1)),
                      substr(type, 2, nchar(type)))
  comparison_title = toupper(comparison_type)
  # Create the selection expressions
  flu_select_expr = c("state", "season", flu_col)
  comparison_select_expr = c("state", "season", comparison_col)
  # Join datasets by state and season
  merged_data = inner_join(
    flu_data %>% select(all_of(flu_select_expr)),
    comparison_data %>% select(all_of(comparison_select_expr)),
    by = c("state", "season")
  )
  # Rename columns for easier access
  names(merged_data)[names(merged_data) == flu_col] = "flu_date"
  names(merged_data)[names(merged_data) == comparison_col] = "comparison_date"

  # Calculate point density for overlapping points
  # First, create rounded versions of dates to detect nearby points
  merged_data = merged_data %>%
    mutate(
      comparison_date_rounded = as.Date(format(comparison_date, "%Y-%m-%d")),
      flu_date_rounded = as.Date(format(flu_date, "%Y-%m-%d"))
    )
  # Count the number of points at each rounded position
  point_counts = merged_data %>%
    group_by(comparison_date_rounded, flu_date_rounded) %>%
    summarize(count = n(), .groups = "drop")
  # Join the counts back to the original data
  merged_data = merged_data %>%
    left_join(point_counts,
              by = c("comparison_date_rounded", "flu_date_rounded"))

  # Calculate the ratio of points above the line to points below/on the line
  # Points above the line: flu_date > comparison_date (flu is later)
  # Points on or below the line: flu_date <= comparison_date (flu is earlier or same time)
  points_above =
    sum(merged_data$flu_date > merged_data$comparison_date, na.rm = TRUE)
  points_below_or_on =
    sum(merged_data$flu_date <= merged_data$comparison_date, na.rm = TRUE)

  # Calculate ratio (handle division by zero)
  if (points_below_or_on == 0) {
    ratio_text = "All points above line"
  } else {
    ratio = round(points_above / points_below_or_on, 2)
    ratio_text = paste0("Ratio (above/below line): ", ratio)
  }

  # Find the full date range to ensure equal scaling
  all_dates = c(merged_data$comparison_date, merged_data$flu_date)
  min_date = min(all_dates, na.rm = TRUE)
  max_date = max(all_dates, na.rm = TRUE)
  date_range = max_date - min_date
  # Add some padding to ensure square aspect and 45Â° line
  min_date = min_date - date_range * 0.1
  max_date = max_date + date_range * 0.1

  # Create the scatterplot
  p = ggplot(merged_data, aes(x = comparison_date, y = flu_date)) +
    # Add reference line (x=y)
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dashed",
      color = "gray50"
    ) +
    # Add points with alpha and size based on count
    geom_point(
      aes(alpha = count, size = count),
      color = base_color  # Use the color parameter
    ) +
    # Scale alpha for better visibility of overlaps
    scale_alpha_continuous(range = c(0.4, 0.9), guide = "none") +
    # Scale size for better visibility of overlaps
    scale_size_continuous(range = c(2, 5), guide = "none") +
    # Set identical limits for both axes to ensure square aspect ratio and 45Â° line
    scale_x_date(
      date_breaks = "2 week",
      date_labels = "%b %d",
      limits = c(min_date, max_date)
    ) +
    scale_y_date(
      date_breaks = "2 week",
      date_labels = "%b %d",
      limits = c(min_date, max_date)
    ) +
    # Enforce square aspect ratio
    coord_fixed(ratio = 1) +
    # Add labels and title
    labs(
      x = paste0(comparison_title, " ", type_title, " Date"),
      y = paste0("Flu ", type_title, " Date")
    ) +
    # Use a clean theme
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank(),
      aspect.ratio = 1          # Reinforce square aspect
    ) +
    # Add a season label in a box at the bottom right
    annotate(
      "label",
      x = max_date - date_range * 0.05,
      # Position near right edge
      y = min_date + date_range * 0.05,
      # Position near bottom edge
      label = season,
      hjust = 1,
      vjust = 0,
      size = 3.5,
      fontface = "bold",
      label.size = 0.5,
      # Border thickness
      label.padding = unit(0.35, "lines"),
      # Padding around text
      fill = "white",
      # Background color
      alpha = 0.8
    )

  return(p)
}

# # 3. Onset and peak plots for each season -------------------------------

## Peaks - Flu vs. RSV
peak_plot_flu_rsv_23_24 = create_comparison_plot(
  flu_dat_23,
  rsv_dat_23,
  comparison_type = "rsv",
  season = "23-24",
  type = "peak",
  base_color = "darkslateblue"
)
peak_plot_flu_rsv_23_24
# ggsave(file = "figures/manuscript_figures/fig3/peak_scatterplot_flu_rsv_23_24.png", width = 10, height = 8, units = "in", bg = "white")

peak_plot_flu_rsv_24_25 = create_comparison_plot(
  flu_dat_24,
  rsv_dat_24,
  comparison_type = "rsv",
  season = "24-25",
  type = "peak",
  base_color = "dodgerblue2"
)
peak_plot_flu_rsv_24_25
# ggsave(file = "figures/manuscript_figures/fig3/peak_scatterplot_flu_rsv_24_25.png", width = 10, height = 8, units = "in", bg = "white")

## Onsets - Flu vs. RSV
onset_plot_flu_rsv_23_24 = create_comparison_plot(
  flu_dat_23,
  rsv_dat_23,
  comparison_type = "rsv",
  season = "23-24",
  type = "onset",
  base_color = "darkslateblue"
)
onset_plot_flu_rsv_23_24
# ggsave(file = "figures/manuscript_figures/fig3/onset_scatterplot_flu_rsv_23_24.png", width = 10, height = 8, units = "in", bg = "white")

onset_plot_flu_rsv_24_25 = create_comparison_plot(
  flu_dat_24,
  rsv_dat_24,
  comparison_type = "rsv",
  season = "24-25",
  type = "onset",
  base_color = "dodgerblue2"
)
onset_plot_flu_rsv_24_25
# ggsave(file = "figures/manuscript_figures/fig3/onset_scatterplot_flu_rsv_24_25.png", width = 10, height = 8, units = "in", bg = "white")

## Peaks -- Flu vs. COVID
peak_plot_flu_cov_23_24 = create_comparison_plot(
  flu_dat_23,
  cov_dat_23,
  comparison_type = "covid",
  season = "23-24",
  type = "peak",
  base_color = "darkslategrey"
)
peak_plot_flu_cov_23_24
# ggsave(file = "figures/manuscript_figures/fig3/peak_scatterplot_flu_cov_23_24.png", width = 10, height = 8, units = "in", bg = "white")

peak_plot_flu_cov_24_25 = create_comparison_plot(
  flu_dat_24,
  cov_dat_24,
  comparison_type = "covid",
  season = "24-25",
  type = "peak",
  base_color = "mediumseagreen"
)
peak_plot_flu_cov_24_25
# ggsave(file = "figures/manuscript_figures/fig3/peak_scatterplot_flu_cov_24_25.png", width = 10, height = 8, units = "in", bg = "white")

## Onsets -- Flu vs. COVID
onset_plot_flu_cov_23_24 = create_comparison_plot(
  flu_dat_23,
  cov_dat_23,
  comparison_type = "covid",
  season = "23-24",
  type = "onset",
  base_color = "darkslategrey"
)
onset_plot_flu_cov_23_24
# ggsave(file = "figures/manuscript_figures/fig3/onset_scatterplot_flu_cov_23_24.png", width = 10, height = 8, units = "in", bg = "white")

onset_plot_flu_cov_24_25 = create_comparison_plot(
  flu_dat_24,
  cov_dat_24,
  comparison_type = "covid",
  season = "24-25",
  type = "onset",
  base_color = "mediumseagreen"
)
onset_plot_flu_cov_24_25
# ggsave(file = "figures/manuscript_figures/fig3/onset_scatterplot_flu_cov_24_25.png", width = 10, height = 8, units = "in", bg = "white")


# 4. Grouped violin plot for time difference ------------------------------

# First create a merged dataset for the 23-24 season
result_23_24 = flu_dat_23 %>%
  left_join(rsv_dat_23, by = c("state", "season")) %>%
  left_join(
    cov_dat_23 %>%
      group_by(state, season) %>%
      summarize(
        covid_onset_list = list(covid_onset),
        covid_peak_list = list(covid_peak),
        .groups = "drop"
      ),
    by = c("state", "season")
  )

# Process 23-24 season and calculate time differences
result_23_24 = result_23_24 %>%
  rowwise() %>%
  mutate(
    # Find closest COVID onset to flu onset
    closest_covid_onset = if (length(unlist(covid_onset_list)) > 0) {
      closest_date = unlist(covid_onset_list)[which.min(abs(difftime(
        flu_onset, unlist(covid_onset_list), units = "days"
      )))]
      as.Date(closest_date, origin = "1970-01-01")
    } else {
      NA_Date_
    },
    # Find closest COVID peak to flu peak
    closest_covid_peak = if (length(unlist(covid_peak_list)) > 0) {
      closest_date = unlist(covid_peak_list)[which.min(abs(difftime(
        flu_peak, unlist(covid_peak_list), units = "days"
      )))]
      as.Date(closest_date, origin = "1970-01-01")
    } else {
      NA_Date_
    },

    # Calculate time differences in days
    flu_rsv_peak_diff =
      as.numeric(difftime(flu_peak, rsv_peak, units = "days")),
    flu_covid_peak_diff =
      as.numeric(difftime(flu_peak, closest_covid_peak, units = "days")),
    flu_rsv_onset_diff =
      as.numeric(difftime(flu_onset, rsv_onset, units = "days")),
    flu_covid_onset_diff =
      as.numeric(difftime(flu_onset, closest_covid_onset, units = "days")),

    # Convert to weeks
    flu_rsv_peak_diff_weeks = flu_rsv_peak_diff / 7,
    flu_covid_peak_diff_weeks = flu_covid_peak_diff / 7,
    flu_rsv_onset_diff_weeks = flu_rsv_onset_diff / 7,
    flu_covid_onset_diff_weeks = flu_covid_onset_diff / 7
  ) %>%
  ungroup() %>%
  select(-covid_onset_list, -covid_peak_list)

# Now create a merged dataset for the 24-25 season
result_24_25 = flu_dat_24 %>%
  left_join(rsv_dat_24, by = c("state", "season")) %>%
  left_join(
    cov_dat_24 %>%
      group_by(state, season) %>%
      summarize(
        covid_onset_list = list(covid_onset),
        covid_peak_list = list(covid_peak),
        .groups = "drop"
      ),
    by = c("state", "season")
  )

# Process 24-25 season and calculate time differences
result_24_25 = result_24_25 %>%
  rowwise() %>%
  mutate(
    # Find closest COVID onset to flu onset
    closest_covid_onset = if (length(unlist(covid_onset_list)) > 0) {
      closest_date = unlist(covid_onset_list)[which.min(abs(difftime(
        flu_onset, unlist(covid_onset_list), units = "days"
      )))]
      as.Date(closest_date, origin = "1970-01-01")
    } else {
      NA_Date_
    },
    # Find closest COVID peak to flu peak
    closest_covid_peak = if (length(unlist(covid_peak_list)) > 0) {
      closest_date = unlist(covid_peak_list)[which.min(abs(difftime(
        flu_peak, unlist(covid_peak_list), units = "days"
      )))]
      as.Date(closest_date, origin = "1970-01-01")
    } else {
      NA_Date_
    },

    # Calculate time differences in days
    flu_rsv_peak_diff = as.numeric(difftime(flu_peak, rsv_peak, units = "days")),
    flu_covid_peak_diff = as.numeric(difftime(flu_peak, closest_covid_peak, units = "days")),
    flu_rsv_onset_diff = as.numeric(difftime(flu_onset, rsv_onset, units = "days")),
    flu_covid_onset_diff = as.numeric(difftime(flu_onset, closest_covid_onset, units = "days")),

    # Convert to weeks
    flu_rsv_peak_diff_weeks = flu_rsv_peak_diff / 7,
    flu_covid_peak_diff_weeks = flu_covid_peak_diff / 7,
    flu_rsv_onset_diff_weeks = flu_rsv_onset_diff / 7,
    flu_covid_onset_diff_weeks = flu_covid_onset_diff / 7
  ) %>%
  ungroup() %>%
  select(-covid_onset_list, -covid_peak_list)

# Combine both seasons
result_combined = bind_rows(result_23_24, result_24_25)

# Reshape data to long format for plotting
plot_data = result_combined %>%
  select(
    state,
    season,
    flu_rsv_peak_diff_weeks,
    flu_covid_peak_diff_weeks,
    flu_rsv_onset_diff_weeks,
    flu_covid_onset_diff_weeks
  ) %>%
  # Rename columns to more readable format for the plot
  rename(
    "Flu-RSV Peak" = flu_rsv_peak_diff_weeks,
    "Flu-COVID Peak" = flu_covid_peak_diff_weeks,
    "Flu-RSV Onset" = flu_rsv_onset_diff_weeks,
    "Flu-COVID Onset" = flu_covid_onset_diff_weeks
  ) %>%
  # Convert to long format
  pivot_longer(
    cols = c(
      "Flu-RSV Peak",
      "Flu-COVID Peak",
      "Flu-RSV Onset",
      "Flu-COVID Onset"
    ),
    names_to = "comparison",
    values_to = "weeks_difference"
  )

# Create ordered factors
plot_data = plot_data %>%
  mutate(
    # Create ordered factor for the comparison variable
    comparison = factor(
      comparison,
      levels = c(
        "Flu-RSV Peak",
        "Flu-RSV Onset",
        "Flu-COVID Peak",
        "Flu-COVID Onset"
      )
    ),
    # Make the season a factor for proper ordering
    season = factor(season, levels = c("23-24", "24-25"))
  )

plot_data

# Create a grouped violin plot with season as the grouping factor
violin_peaks = plot_data %>% filter(str_detect(comparison, "Peak")) %>%
  ggplot(aes(x = comparison, y = weeks_difference, fill = season)) +
  # Add violins with position dodge to group by season and lighter transparency
  geom_violin(
    position = position_dodge(width = 0.8),
    width = 1,
    alpha = 0.5,
    trim = FALSE
  ) +

  # Add boxplots inside violins with darker version of the same fill colors
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.15,
    alpha = 0.9,
    outlier.shape = NA,
    color = "black"
  ) +

  # Add reference line at zero
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "darkgray") +

  # Add mean points
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23,
    position = position_dodge(width = 0.8),
    size = 2.5,
    color = "black"
  ) +

  # Custom colors for seasons
  scale_fill_manual(values = c(
    "23-24" = "#1f78b4",
    # Blue
    "24-25" = "#33a02c"   # Green
  )) +

  # Clear labeling
  labs(
    x = "",
    # No label needed for categorical x-axis
    y = "Weeks Difference\n(+ means flu occurred later)",
    fill = "Season"
  ) +

  # Clean theme with good spacing
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(
      size = 11,
      face = "bold",
      angle = 45,
      hjust = 1
    ),
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "none",
  )

violin_peaks_horizontal = plot_data %>%
  filter(str_detect(comparison, "Peak")) %>%
  ggplot(aes(y = comparison, x = weeks_difference, fill = season)) +
  # Add violins with position dodge to group by season and lighter transparency
  geom_violin(
    position = position_dodge(width = 0.8),
    width = 1,
    alpha = 0.5,
    trim = FALSE
  ) +
  # Add boxplots inside violins with darker version of the same fill colors
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.15,
    alpha = 0.9,
    outlier.shape = NA,
    color = "black"
  ) +
  # Add reference line at zero
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "darkgray") +
  # Add mean points
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23,
    position = position_dodge(width = 0.8),
    size = 2.5,
    color = "black"
  ) +
  # Custom colors for seasons
  scale_fill_manual(values = c(
    "23-24" = "#1f78b4",
    # Blue
    "24-25" = "#33a02c"   # Green
  )) +
  # Clear labeling
  labs(
    y = "",
    # No label needed for categorical y-axis
    x = "Weeks Difference\n(+ means flu occurred later)",
    fill = "Season"
  ) +
  # Clean theme with good spacing
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 11, face = "bold"),
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "none"
  )

violin_onset_horizontal = plot_data %>%
  filter(str_detect(comparison, "Onset")) %>%
  ggplot(aes(y = comparison, x = weeks_difference, fill = season)) +
  # Add violins with position dodge to group by season and lighter transparency
  geom_violin(
    position = position_dodge(width = 0.8),
    width = 1,
    alpha = 0.5,
    trim = FALSE
  ) +
  # Add boxplots inside violins with darker version of the same fill colors
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.15,
    alpha = 0.9,
    outlier.shape = NA,
    color = "black"
  ) +
  # Add reference line at zero
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "darkgray") +
  # Add mean points
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 23,
    position = position_dodge(width = 0.8),
    size = 2.5,
    color = "black"
  ) +
  # Custom colors for seasons
  scale_fill_manual(values = c(
    "23-24" = "#1f78b4",
    # Blue
    "24-25" = "#33a02c"   # Green
  )) +
  # Clear labeling
  labs(
    y = "",
    # No label needed for categorical y-axis
    x = "Weeks Difference\n(+ means flu occurred later)",
    fill = "Season"
  ) +
  # Clean theme with good spacing
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 11, face = "bold"),
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "none"
  )

violin_onset_horizontal
# ggsave(file = "figures/manuscript_figures/ews_diff_violin_plot.png",
#         width = 10, height = 8, units = "in", bg = "white")

# 5. Getting point statistics from scatterplots

## Assembling composite figure for manuscript

## RSV
scatter_rsv_2x2 = (peak_plot_flu_rsv_23_24 + onset_plot_flu_rsv_23_24) /
  (peak_plot_flu_rsv_24_25 + onset_plot_flu_rsv_24_25)
scatter_rsv_2x2
ggsave(
  file = "figures/manuscript_figures/fig3/flu_rsv_peaks_onsets.png",
  width = 14,
  height = 10,
  units = "in",
  bg = "white"
)

## COVID
scatter_covid_2x2 = (peak_plot_flu_cov_23_24 + onset_plot_flu_cov_23_24) /
  (peak_plot_flu_cov_24_25 + onset_plot_flu_cov_24_25)
scatter_covid_2x2
ggsave(
  file = "figures/manuscript_figures/fig3/flu_cov_peaks_onsets.png",
  width = 14,
  height = 10,
  units = "in",
  bg = "white"
)

# Function to add a subplot annotation (a, b, c, etc.) with better positioning
add_subplot_label = function(plot, label, position = c(0.02, 0.98)) {
  plot +
    theme(plot.tag = element_text(size = 14, face = "bold"),
          plot.tag.position = position) +
    labs(tag = label)
}

# Add subplot labels to each plot
peak_plot_flu_rsv_23_24_labeled =
  add_subplot_label(peak_plot_flu_rsv_23_24, "A")
peak_plot_flu_cov_23_24_labeled =
  add_subplot_label(peak_plot_flu_cov_23_24, "B")
peak_plot_flu_rsv_24_25_labeled =
  add_subplot_label(peak_plot_flu_rsv_24_25, "C")
peak_plot_flu_cov_24_25_labeled =
  add_subplot_label(peak_plot_flu_cov_24_25, "D")

# Do the same for onsets
onset_plot_flu_rsv_23_24_labeled =
  add_subplot_label(onset_plot_flu_rsv_23_24, "A")
onset_plot_flu_cov_23_24_labeled =
  add_subplot_label(onset_plot_flu_cov_23_24, "B")
onset_plot_flu_rsv_24_25_labeled =
  add_subplot_label(onset_plot_flu_rsv_24_25, "C")
onset_plot_flu_cov_24_25_labeled =
  add_subplot_label(onset_plot_flu_cov_24_25, "D")

# For the violin plot, position the label further to the left to avoid y-axis overlap
violin_labeled =
  add_subplot_label(violin_peaks_horizontal, "E", position = c(-0.08, 0.98))

violin_onset_labeled =
  add_subplot_label(violin_onset_horizontal, "E", position = c(-0.08, 0.98))

# Set consistent dimensions for the scatterplots with reduced margins
set_plot_dimensions = function(plot) {
  plot +
    theme(
      aspect.ratio = 1,
      # Ensure square plots
    )
}

# Apply consistent dimensions to all plots
peak_plot_flu_rsv_23_24_labeled =
  set_plot_dimensions(peak_plot_flu_rsv_23_24_labeled)
peak_plot_flu_cov_23_24_labeled =
  set_plot_dimensions(peak_plot_flu_cov_23_24_labeled)
peak_plot_flu_rsv_24_25_labeled =
  set_plot_dimensions(peak_plot_flu_rsv_24_25_labeled)
peak_plot_flu_cov_24_25_labeled =
  set_plot_dimensions(peak_plot_flu_cov_24_25_labeled)

# Apply consistent dimensions to all plots
onset_plot_flu_cov_23_24_labeled =
  set_plot_dimensions(onset_plot_flu_cov_23_24_labeled)
onset_plot_flu_rsv_23_24_labeled =
  set_plot_dimensions(onset_plot_flu_rsv_23_24_labeled)
onset_plot_flu_rsv_24_25_labeled =
  set_plot_dimensions(onset_plot_flu_rsv_24_25_labeled)
onset_plot_flu_cov_24_25_labeled =
  set_plot_dimensions(onset_plot_flu_cov_24_25_labeled)


# For the violin plot, adjust margins without enforcing aspect ratio
violin_labeled_sq = violin_labeled +
  theme(
    aspect.ratio = 1,
    plot.tag.position = c(0.05, 0.98),
    # Moved further right from the default
    plot.tag = element_text(size = 14, face = "bold")
  )

violin_labeled_onset_sq = violin_onset_labeled +
  theme(
    aspect.ratio = 1,
    plot.tag.position = c(0.05, 0.98),
    # Moved further right from the default
    plot.tag = element_text(size = 14, face = "bold")
  )

# Create empty plot for spacing
empty_plot = ggplot() + theme_void()

# Create the combined plot with the violin plot centered
combined_peak_plot = (
  # First row of scatterplots
  peak_plot_flu_rsv_23_24_labeled + peak_plot_flu_cov_23_24_labeled
) / (
  # Second row of scatterplots
  peak_plot_flu_rsv_24_25_labeled + peak_plot_flu_cov_24_25_labeled
) / (
  # Third row with centered violin plot
  empty_plot + violin_labeled_sq + empty_plot +
    plot_layout(widths = c(0.25, 0.5, 0.25))
) +
  # Overall layout settings
  plot_layout(
    heights = c(1, 1, 1)  # Equal heights for rows
  ) &
  # Reduce spacing between plots
  theme(panel.spacing = unit(1, "lines"))

# Display the combined plot
combined_peak_plot
ggsave(
  "figures/manuscript_figures/fig2/fig2_070125.png",
  height = 10,
  width = 8,
  units = "in",
  bg = "white"
)

# Create the combined plot with the violin plot centered
combined_onset_plot = (
  # First row of scatterplots
  onset_plot_flu_rsv_23_24_labeled + onset_plot_flu_cov_23_24_labeled
) / (
  # Second row of scatterplots
  onset_plot_flu_rsv_24_25_labeled + onset_plot_flu_cov_24_25_labeled
) / (
  # Third row with centered violin plot
  empty_plot + violin_labeled_onset_sq + empty_plot +
    plot_layout(widths = c(0.25, 0.5, 0.25))
) +
  # Overall layout settings
  plot_layout(
    heights = c(1, 1, 1)  # Equal heights for rows
  ) &
  # Reduce spacing between plots
  theme(panel.spacing = unit(1, "lines"))

combined_onset_plot
ggsave(
  "figures/manuscript_figures/supplementary/supp_fig_10_070125.png",
  height = 10,
  width = 8,
  units = "in",
  bg = "white"
)



library(boot)
# Function to calculate median in bootstrap samples
median_boot = function(data, indices) {
  return(median(data[indices], na.rm = TRUE))
}

plot_data %>%
  filter(comparison == "Flu-RSV Peak") %>%
  group_by(comparison) %>%
  summarize(min = min(weeks_difference, na.rm = TRUE),
            max = max(weeks_difference, na.rm = TRUE),
            median = median(weeks_difference, na.rm = TRUE),
            mean = mean(weeks_difference, na.rm = TRUE))

# Expanded summary function with 95% CI for median
onset_peak_stats = plot_data %>%
  group_by(season, comparison) %>%
  summarize(
    median_diff = median(weeks_difference, na.rm = TRUE),
    mean_diff = mean(weeks_difference, na.rm = TRUE),
    # Create bootstrap objects for the median
    boot_obj = list(boot(weeks_difference, median_boot, R = 1000)),
    # Extract confidence intervals from bootstrap
    ci_lower = map_dbl(boot_obj, ~ quantile(.$t, probs = 0.025)),
    ci_upper = map_dbl(boot_obj, ~ quantile(.$t, probs = 0.975)),
    # Calculate the 2.5th and 97.5th percentiles directly (alternate method)
    pct_2.5 = quantile(weeks_difference, probs = 0.025, na.rm = TRUE),
    pct_97.5 = quantile(weeks_difference, probs = 0.975, na.rm = TRUE),
    # Also include the min and max values to calculate the range
    min = min(weeks_difference, na.rm = TRUE),
    max = max(weeks_difference, na.rm = TRUE)
  ) %>%
  # Clean up by removing the boot object column
  select(-boot_obj)

# Save the results
write_csv(onset_peak_stats, file = "data/onset_peak_stats.csv")
