# Calculating effectiveness of ACIP recommendations for RSV prophylaxis

# The general guideline for nirsevimab is to administer it to infants born
# during the RSV season, which typically runs from October 1 to March 31 in the
# Northern Hemisphere.
# The goal is to provide protection against RSV during the peak season when the
# virus is most prevalent.

# We use the NSSP data from the RSV timing project to calculate what percentage
# of RSV cases occur outside the recommended prophylaxis period.

# Load libraries
library(tidyverse)
library(lubridate)

# Load the data
setwd("/Users/gdewey/Documents/Projects/rsv_peak_timing")
nssp_data <- read_csv("data/NSSP/nssp_full_07292025.csv")

# State ED data
ed_data <- read_csv("data/ed_rates_state_per_1000.csv")

# State populations
state_pops <- read_csv("data/us_census_state_pop.csv")

# Do some data wrangling
state_pops_ <- state_pops %>%
  select(NAME, ESTIMATESBASE2020, POPESTIMATE2022,
         POPESTIMATE2023, POPESTIMATE2024) %>%
  filter(NAME %in% c(datasets::state.name, "District of Columbia")) %>%
  rename(state = NAME, pop2020 = ESTIMATESBASE2020, pop2022 = POPESTIMATE2022,
         pop2023 = POPESTIMATE2023, pop2024 = POPESTIMATE2024)

ed_data_ <- ed_data %>%
  select(Location, Total) %>%
  rename(state = Location, ed_visits_per_1000 = Total) %>%
  filter(state != "United States")

total_ed_visits_by_state <- state_pops_ %>%
  left_join(ed_data_, by = "state") %>%
  mutate(total_ed_visits_2023 = (ed_visits_per_1000 * pop2020) / 1000) %>%
  select(state, total_ed_visits_2023)

nssp_data %>%
  select(week_end, geography, county, percent_visits_rsv) %>%
  filter(county == "All") %>%
  left_join(total_ed_visits_by_state, by = c("geography" = "state")) %>%
  mutate(estimated_RSV_visits =
           (percent_visits_rsv / 100) * total_ed_visits_2023) %>%
  select(-county)

rsv_dat <- nssp_data %>%
  select(week_end, geography, county, percent_visits_rsv) %>%
  filter(county == "All", week_end > as_date("2023-06-01")) %>%
  # Create RSV season variable
  mutate(
    rsv_season = case_when(
      week_end >= as.Date("2023-06-01") & week_end <= as.Date("2024-05-31") ~
        "2023-24",
      week_end >= as.Date("2024-06-01") ~ "2024-25",
      TRUE ~ NA_character_
    ),
    # Extract month for seasonal classification
    month = month(week_end),
    # Classify as peak season (Oct-Mar) or off-season
    season_period = case_when(
      month %in% c(5, 6, 7, 8, 9) ~ "Before ACIP range",
      month %in% c(10, 11, 12, 1, 2, 3) ~ "During ACIP range",
      month %in% c(4) ~ "After ACIP range",
      TRUE ~ NA_character_
    )
  ) %>%
  # Filter out any rows without RSV season classification
  filter(!is.na(rsv_season)) %>%
  # Group and summarize
  group_by(geography, rsv_season, season_period) %>%
  summarise(
    sum_percent_visits_rsv = sum(percent_visits_rsv, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Pivot to get columns for each time period
  pivot_wider(
    names_from = season_period,
    values_from = sum_percent_visits_rsv,
    values_fill = 0
  ) %>%
  # Calculate total across all months
  mutate(
  total_all_months = `Before ACIP range` + `During ACIP range` +
    `After ACIP range`
) %>%
# Reorder columns for clarity
select(geography, rsv_season, `Before ACIP range`, `During ACIP range`,
       `After ACIP range`, total_all_months) %>%
  # Reorder columns for clarity
  select(geography, rsv_season, `Before ACIP range`,
         `During ACIP range`, `After ACIP range`, total_all_months) %>%
  left_join(total_ed_visits_by_state, by = c("geography" = "state")) %>%
  mutate(visits_before_ACIP =
           total_ed_visits_2023 * `Before ACIP range` / 1000,
         visits_during_ACIP =
           total_ed_visits_2023 * `During ACIP range` / 1000,
         visits_after_ACIP = total_ed_visits_2023 * `After ACIP range` / 1000,
         visits_total = total_all_months * total_ed_visits_2023 / 1000)

# Calculating population adjusted rates
rsv_dat_pop_adj <- nssp_data %>%
  select(week_end, geography, county, percent_visits_rsv) %>%
  filter(county == "All", week_end > as_date("2023-06-01")) %>%
  # Create RSV season variable
  mutate(
    rsv_season = case_when(
      week_end >= as.Date("2023-06-01") & week_end <= as.Date("2024-05-31") ~
        "2023-24",
      week_end >= as.Date("2024-06-01") ~ "2024-25",
      TRUE ~ NA_character_
    ),
    # Extract month for seasonal classification
    month = month(week_end),
    # Classify as peak season (Oct-Mar) or off-season
    season_period = case_when(
      month %in% c(5, 6, 7, 8, 9) ~ "Before ACIP range",
      month %in% c(10, 11, 12, 1, 2, 3) ~ "During ACIP range",
      month %in% c(4) ~ "After ACIP range",
      TRUE ~ NA_character_
    )
  ) %>%
  # Filter out any rows without RSV season classification
  filter(!is.na(rsv_season)) %>%
  # Group and summarize
  group_by(geography, rsv_season, season_period) %>%
  summarise(
    sum_percent_visits_rsv = sum(percent_visits_rsv, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Pivot to get columns for each time period
  pivot_wider(
    names_from = season_period,
    values_from = sum_percent_visits_rsv,
    values_fill = 0
  ) %>%
  # Calculate total across all months
  mutate(
    total_all_months = `Before ACIP range` + `During ACIP range` +
      `After ACIP range`
  ) %>%
  # Reorder columns for clarity
  select(geography, rsv_season, `Before ACIP range`, `During ACIP range`,
         `After ACIP range`, total_all_months) %>%
  left_join(total_ed_visits_by_state, by = c("geography" = "state")) %>%
  left_join(state_pops_, by = c("geography" = "state")) %>%
  mutate(
    # Calculate actual visit numbers
    visits_before_ACIP = total_ed_visits_2023 * `Before ACIP range` / 1000,
    visits_during_ACIP = total_ed_visits_2023 * `During ACIP range` / 1000,
    visits_after_ACIP = total_ed_visits_2023 * `After ACIP range` / 1000,
    visits_total = total_all_months * total_ed_visits_2023 / 1000,
    # Calculate population-adjusted rates (per 100,000 population)
    rate_before_ACIP_per_100k = (visits_before_ACIP / pop2020) * 100000,
    rate_during_ACIP_per_100k = (visits_during_ACIP / pop2020) * 100000,
    rate_after_ACIP_per_100k = (visits_after_ACIP / pop2020) * 100000,
    rate_total_per_100k = (visits_total / pop2020) * 100000,
    # Calculate what % of all ED visits these RSV visits represent
    pct_ed_visits_before_ACIP = (visits_before_ACIP / total_ed_visits_2023) * 100,
    pct_ed_visits_during_ACIP = (visits_during_ACIP / total_ed_visits_2023) * 100,
    pct_ed_visits_after_ACIP = (visits_after_ACIP / total_ed_visits_2023) * 100,
    pct_ed_visits_total = (visits_total / total_ed_visits_2023) * 100
  )

rsv_dat_pop_adj <- rsv_dat_pop_adj %>%
  select(-c(pop2022, pop2023, pop2024))

rsv_regional_rates = rsv_dat_pop_adj %>%
  left_join(us_regions, by = "geography") %>%
  group_by(rsv_season, region) %>%
  summarize(mean_rate_before_ACIP_per_100k = mean(rate_before_ACIP_per_100k),
            mean_rate_during_ACIP_per_100k = mean(rate_during_ACIP_per_100k),
            mean_rate_after_ACIP_per_100k = mean(rate_after_ACIP_per_100k)) %>%
  filter(region != "National")

rsv_state_level_rates = rsv_dat_pop_adj %>%
  left_join(us_regions, by = "geography") %>%
  group_by(rsv_season, region) %>%
  arrange(region) %>%
  select(geography, rsv_season, starts_with("rate"))

write_csv(rsv_dat_pop_adj, file = "data/rsv_dat_population_adjusted.csv")
write_csv(rsv_regional_rates, file = "data/rsv_regional_mean_rates.csv")
write_csv(rsv_state_level_rates, file = "data/rsv_state_level_rates.csv")


# Define US regions 
us_regions = data.frame(
  geography = c(
    # Northeast
    "Connecticut", "Maine", "Massachusetts", "New Hampshire", "Rhode Island",
    "Vermont", "New Jersey", "New York", "Pennsylvania",
    # Midwest
    "Illinois", "Indiana", "Michigan", "Ohio", "Wisconsin", "Iowa", "Kansas",
    "Minnesota", "Missouri", "Nebraska", "North Dakota", "South Dakota",
    # South
    "Delaware", "Florida", "Georgia", "Maryland", "North Carolina",
    "South Carolina","Virginia", "West Virginia", "Alabama", "Kentucky",
    "Mississippi", "Tennessee","Arkansas", "Louisiana", "Oklahoma", "Texas",
    "District of Columbia",
    # West
    "Arizona", "Colorado", "Idaho", "Montana", "Nevada", "New Mexico", "Utah",
    "Wyoming", "Alaska", "California", "Hawaii", "Oregon", "Washington",
    # National
    "United States"
  ),
  region = c(
    # Northeast (9 states)
    rep("Northeast", 9),
    # Midwest (12 states)
    rep("Midwest", 12),
    # South (17 states/territories)
    rep("South", 17),
    # West (13 states)
    rep("West", 13),
    # National
    "National"
  )
)

#  Prepare data with regional grouping and ordering
plot_data_regional = rsv_dat %>%
  select(geography, rsv_season, `Before ACIP range`, `During ACIP range`, 
         `After ACIP range`, total_all_months) %>%
  mutate(
    perc_before_range = `Before ACIP range`/total_all_months,
    perc_during_range = `During ACIP range`/total_all_months,
    perc_after_range = `After ACIP range`/total_all_months
  ) %>%
  # Add regional classifications
  left_join(us_regions, by = "geography") %>%
  # Calculate % outside ACIP for 2023-24 season only for ordering (before + after)
  filter(rsv_season == "2023-24") %>%
  group_by(geography, region) %>%
  summarise(
    perc_outside_2023_24 = first(perc_before_range + perc_after_range),
    .groups = "drop"
  ) %>%
  # Create ordering: National first, then regions, then states within regions by 2023-24 data
  mutate(
    region_order = case_when(
      region == "National" ~ 1,
      region == "Northeast" ~ 2,
      region == "Midwest" ~ 3,
      region == "South" ~ 4,
      region == "West" ~ 5,
      TRUE ~ 6
    )
  ) %>%
  arrange(region_order, desc(perc_outside_2023_24)) %>%
  mutate(
    geography = factor(geography, levels = geography),
    region = factor(region, levels = c("National", "Northeast", "Midwest", "South", "West"))
  ) %>%
  # Rejoin with original data
  select(geography, region) %>%
  left_join(
    rsv_dat %>%
      select(geography, rsv_season, `Before ACIP range`, `During ACIP range`,
             `After ACIP range`, total_all_months) %>%
      mutate(
        perc_before_range = `Before ACIP range`/total_all_months,
        perc_during_range = `During ACIP range`/total_all_months,
        perc_after_range = `After ACIP range`/total_all_months
      ),
    by = "geography"
  ) %>%
  # Reshape for plotting
  pivot_longer(
    cols = c(perc_before_range, perc_during_range, perc_after_range),
    names_to = "period",
    values_to = "percentage"
  ) %>%
  mutate(
    period = case_when(
      period == "perc_before_range" ~ "Before ACIP Range",
      period == "perc_during_range" ~ "During ACIP Range",
      period == "perc_after_range" ~ "After ACIP Range"
    )
  )


# Define region-specific color pairs (8 colors total)
region_colors = list(
  "National" = c("Outside ACIP Range" = "#2E8B57", "During ACIP Range" = "#98FB98"),      # Sea Green pair
  "Northeast" = c("Outside ACIP Range" = "#B22222", "During ACIP Range" = "#FFB6C1"),    # Fire Brick pair  
  "Midwest" = c("Outside ACIP Range" = "#4169E1", "During ACIP Range" = "#87CEEB"),      # Royal Blue pair
  "South" = c("Outside ACIP Range" = "#FF8C00", "During ACIP Range" = "#FFE4B5"),        # Dark Orange pair
  "West" = c("Outside ACIP Range" = "#8B008B", "During ACIP Range" = "#DDA0DD")          # Dark Magenta pair
)

# Create a mapping dataframe for colors
color_mapping = plot_data_regional %>%
  select(region) %>%
  distinct() %>%
  mutate(
    outside_color = case_when(
      region == "National" ~ "#2E8B57",
      region == "Northeast" ~ "#B22222", 
      region == "Midwest" ~ "#4169E1",
      region == "South" ~ "#FF8C00",
      region == "West" ~ "#8B008B"
    ),
    during_color = case_when(
      region == "National" ~ "#98FB98",
      region == "Northeast" ~ "#FFB6C1",
      region == "Midwest" ~ "#87CEEB", 
      region == "South" ~ "#FFE4B5",
      region == "West" ~ "#DDA0DD"
    )
  )

# Add color information to plot data
plot_data_with_colors <- plot_data_regional %>%
  left_join(color_mapping, by = "region") %>%
  mutate(
    color_value = case_when(
      period == "Outside ACIP Range" ~ outside_color,
      period == "During ACIP Range" ~ during_color
    ),
    # Create unique identifier for color mapping
    region_period = paste(region, period, sep = "_")
  )

# Create all possible color combinations for the scale
all_colors = c(
  "National_Before ACIP Range" = "#8B4513", "National_During ACIP Range" = "#DEB887", "National_After ACIP Range" = "#8B4513",
  "Northeast_Before ACIP Range" = "#DC143C", "Northeast_During ACIP Range" = "#FFC0CB", "Northeast_After ACIP Range" = "#DC143C",
  "Midwest_Before ACIP Range" = "#4682B4", "Midwest_During ACIP Range" = "#B0E0E6", "Midwest_After ACIP Range" = "#4682B4",
  "South_Before ACIP Range" = "#FF6347", "South_During ACIP Range" = "#FFEFD5", "South_After ACIP Range" = "#FF6347",
  "West_Before ACIP Range" = "#9932CC", "West_During ACIP Range" = "#E6E6FA", "West_After ACIP Range" = "#9932CC"
)

# Create the plot
ggplot(plot_data_with_colors, aes(x = percentage, y = geography, fill = region_period)) +
  geom_col(position = "stack", width = 0.7) +
  facet_grid(region ~ rsv_season, scales = "free_y", space = "free_y") +
  scale_fill_manual(
  values = all_colors,
  labels = c(
    "National_Before ACIP Range" = "Outside ACIP Range (Dark)",
    "National_During ACIP Range" = "During ACIP Range (Light)",
    "National_After ACIP Range" = "Outside ACIP Range (Dark)",
    "Northeast_Before ACIP Range" = "Outside ACIP Range (Dark)",
    "Northeast_During ACIP Range" = "During ACIP Range (Light)",
    "Northeast_After ACIP Range" = "Outside ACIP Range (Dark)",
    "Midwest_Before ACIP Range" = "Outside ACIP Range (Dark)",
    "Midwest_During ACIP Range" = "During ACIP Range (Light)",
    "Midwest_After ACIP Range" = "Outside ACIP Range (Dark)",
    "South_Before ACIP Range" = "Outside ACIP Range (Dark)",
    "South_During ACIP Range" = "During ACIP Range (Light)",
    "South_After ACIP Range" = "Outside ACIP Range (Dark)",
    "West_Before ACIP Range" = "Outside ACIP Range (Dark)",
    "West_During ACIP Range" = "During ACIP Range (Light)",
    "West_After ACIP Range" = "Outside ACIP Range (Dark)"
  ),
  name = "Time Period"
) +
  scale_x_continuous(
    labels = scales::percent_format(),
    expand = c(0, 0)
  ) +
  scale_y_discrete(limits = rev) +
  labs(
    title = "RSV Visits by ACIP Range Period and US Region",
    x = "Percentage of Total RSV Visits",
    y = "State"
  ) +
  envalysis::theme_publish() +
  theme(
    axis.text.y = element_text(size = 7),
    axis.text.x = element_text(size = 7),
    plot.title = element_text(size = 14, face = "bold"),
    legend.position = "none",
    strip.text.y = element_text(angle = 0, hjust = 0),
    panel.spacing = unit(1, "lines")
  ) +
  guides(fill = guide_legend(nrow = 2))
ggsave(file = "figures/acip_comparison/acip_comparison_v1.png",
       height = 10, width = 8, units = "in", bg = "white")

  

rsv_dat %>%
  left_join(us_regions, by = "geography") %>%
  group_by(region, rsv_season) %>%
  summarize(avg_before_perc_RSV = 
              mean(`Before ACIP range` / total_all_months, na.rm = TRUE),
            median_before_perc_RSV =
              median(`Before ACIP range` / total_all_months, na.rm = TRUE),
            avg_after_perc_RSV =
              mean(`After ACIP range` / total_all_months, na.rm = TRUE),
            median_after_perc_RSV =
              median(`After ACIP range` / total_all_months, na.rm = TRUE))

rsv_dat %>%
  select(geography, rsv_season, total_ed_visits_2023, visits_before_ACIP,
         visits_during_ACIP, visits_after_ACIP, visits_total) %>%
  mutate(perc_before_ACIP = visits_before_ACIP / visits_total,
         perc_during_ACIP = visits_during_ACIP / visits_total,
         perc_after_ACIP = visits_after_ACIP / visits_total) %>%
  select(geography, rsv_season, starts_with("perc")) %>%
  view()

rsv_dat %>%
  select(geography, rsv_season, `Before ACIP range`, `During ACIP range`, 
         `After ACIP range`, total_all_months) %>%
  mutate(
    perc_before_range = `Before ACIP range`/total_all_months,
    perc_during_range = `During ACIP range`/total_all_months,
    perc_after_range = `After ACIP range`/total_all_months
  ) %>%
  # Add regional classifications
  left_join(us_regions, by = "geography") %>%
  group_by(region, rsv_season) %>%
  summarize(mean_before = mean(perc_before_range),
            mean_after = mean(perc_after_range)) %>%
  arrange(rsv_season) %>%
  filter(region != "National")

rsv_dat %>%
  select(geography, rsv_season, `Before ACIP range`, `During ACIP range`, 
         `After ACIP range`, total_all_months) %>%
  mutate(
    perc_before_range = `Before ACIP range`/total_all_months * 100,
    perc_during_range = `During ACIP range`/total_all_months * 100,
    perc_after_range = `After ACIP range`/total_all_months * 100,
  ) %>%
  # Add regional classifications
  left_join(us_regions, by = "geography") %>%
  filter(region == "South") %>%
  select(geography, rsv_season, contains("perc")) %>%
  view()