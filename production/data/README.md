# Data Directory

This directory should contain the raw and processed data files used in the analysis. **Data files are not included in this repository** due to size and data use agreements.

## Required Directory Structure

```
data/
├── NSSP/                      # NSSP emergency department data
│   └── nssp_full_YYMMDD.csv
├── ILInet/                    # ILINet surveillance data
│   └── ILINet_state_YYYY_MM-DD.csv
├── processed/                 # Intermediate processed files
│   ├── nssp_all_years_YYMMDD.csv
│   └── [other processed files]
├── early_warning/            # EWS outputs
│   ├── flu_matched_peaks_onsets.csv
│   ├── rsv_matched_peaks_onsets.csv
│   └── covid_matched_peaks_onsets.csv
├── results/                  # Model outputs
│   ├── results_ridge_no_intercept.Rdata
│   ├── results_ridge_with_intercept.Rdata
│   └── results_linear.Rdata
└── README.md                 # This file
```

## Data Sources

### 1. ILINet Data

**Source**: CDC FluView Interactive  
**URL**: https://www.cdc.gov/flu/weekly/

**How to Download**:
1. Visit the FluView Interactive website
2. Navigate to "Download Data"
3. Select:
   - Geographic Level: State
   - Time Period: 2022-2025 (all available weeks)
   - Data Type: ILI Activity
4. Download as CSV
5. Save to `data/ILInet/`

**Expected Format**:
```
REGION TYPE, REGION, YEAR, WEEK, %WEIGHTED ILI, %UNWEIGHTED ILI, ...
States, Alabama, 2022, 40, 2.3, 2.1, ...
```

### 2. NSSP Data

**Source**: CDC National Syndromic Surveillance Program  
**URL**: https://data.cdc.gov/

**Datasets Needed**:

**A. Emergency Department Visits - COVID-19**
- Search: "COVID-19 Emergency Department Visits"
- Filter: State-level, Weekly
- Years: 2022-2025
- Save as: `nssp_covid_YYMMDD.csv`

**B. Emergency Department Visits - Influenza**
- Search: "Influenza Emergency Department Visits"
- Filter: State-level, Weekly
- Years: 2022-2025
- Save as: `nssp_flu_YYMMDD.csv`

**C. Emergency Department Visits - RSV**
- Search: "RSV Emergency Department Visits"
- Filter: State-level, Weekly
- Years: 2022-2025
- Save as: `nssp_rsv_YYMMDD.csv`

**Expected Combined Format** (`nssp_full_YYMMDD.csv`):
```
week_end, geography, county, percent_visits_covid, percent_visits_influenza, percent_visits_rsv
2022-10-08, Alabama, All, 5.2, 3.1, 2.8
```

**Note**: You may need to combine these files. See "Combining NSSP Files" section below.

### 3. Early Warning System Data (Optional)

If you have access to CDC's Early Warning System outputs, place them in `data/early_warning/`. Otherwise, the `add_EWS_onset_and_peaks.R` script will generate these files.

## Combining NSSP Files

If you download separate NSSP files, use this R code to combine them:

```r
library(tidyverse)

# Load individual files
covid = read_csv("data/NSSP/nssp_covid_YYMMDD.csv")
flu = read_csv("data/NSSP/nssp_flu_YYMMDD.csv")
rsv = read_csv("data/NSSP/nssp_rsv_YYMMDD.csv")

# Standardize column names if needed
covid = covid %>% 
  rename(
    week_end = [date_column_name],
    geography = [state_column_name],
    percent_visits_covid = [covid_column_name]
  )

flu = flu %>% 
  rename(
    week_end = [date_column_name],
    geography = [state_column_name],
    percent_visits_influenza = [flu_column_name]
  )

rsv = rsv %>% 
  rename(
    week_end = [date_column_name],
    geography = [state_column_name],
    percent_visits_rsv = [rsv_column_name]
  )

# Merge files
nssp_combined = covid %>%
  left_join(flu, by = c("week_end", "geography", "county")) %>%
  left_join(rsv, by = c("week_end", "geography", "county"))

# Save combined file
write_csv(nssp_combined, "data/NSSP/nssp_full_YYMMDD.csv")
```

## Data Quality Checks

Before running analysis, verify your data:

```r
library(tidyverse)

# Load data
nssp = read_csv("data/NSSP/nssp_full_YYMMDD.csv")
ili = read_csv("data/ILInet/ILINet_state_YYYY_MM-DD.csv", skip = 1)

# Check date ranges
cat("NSSP date range:", 
    min(nssp$week_end), "to", max(nssp$week_end), "\n")
cat("ILINet date range:", 
    min(ili$WEEK), "to", max(ili$WEEK), "\n")

# Check for missing values
cat("\nMissing values in NSSP:\n")
print(colSums(is.na(nssp)))

cat("\nMissing values in ILINet:\n")
print(colSums(is.na(ili)))

# Check states
cat("\nNumber of states in NSSP:", 
    length(unique(nssp$geography)), "\n")
cat("Number of states in ILINet:", 
    length(unique(ili$REGION)), "\n")

# View first few rows
cat("\nFirst rows of NSSP:\n")
print(head(nssp))

cat("\nFirst rows of ILINet:\n")
print(head(ili))
```

## Data Use and Citation

**Important**: When using CDC data:

1. **Follow CDC data use policies**
2. **Cite data sources appropriately**:
   ```
   Centers for Disease Control and Prevention. National Syndromic 
   Surveillance Program (NSSP). Accessed: [Date]. Available from: 
   https://data.cdc.gov/
   
   Centers for Disease Control and Prevention. ILINet - Influenza-like 
   Illness Surveillance Network. Accessed: [Date]. Available from: 
   https://www.cdc.gov/flu/weekly/
   ```
3. **Do not redistribute raw data** without permission
4. **Acknowledge data providers** in publications

## File Size Expectations

Typical file sizes (approximate):
- ILINet data: ~5-10 MB per file
- NSSP individual files: ~10-20 MB each
- NSSP combined: ~30-50 MB
- Processed files: ~20-40 MB
- Model results: ~1-5 MB

Total space needed: ~500 MB - 1 GB

## Troubleshooting

### Issue: Cannot find data files

**Check**:
1. Files are in correct subdirectories
2. File names match expected format
3. Update paths in R scripts if using different names

### Issue: Date format errors

**Solution**:
```r
# Convert dates to standard format
nssp = nssp %>%
  mutate(week_end = as.Date(week_end, format = "%Y-%m-%d"))
```

### Issue: Column names don't match

**Solution**:
```r
# View actual column names
names(nssp)

# Rename as needed
nssp = nssp %>%
  rename(
    expected_name = actual_name
  )
```

### Issue: Missing states

**Solution**:
- Verify data download included all states
- Some states may have data quality issues
- Check CDC documentation for known gaps

## Privacy and Security

- **Do not commit data files to GitHub** (they're in `.gitignore`)
- **Do not share personally identifiable information**
- **Aggregate data only** - individual-level data is not used in this analysis
- **Follow institutional IRB/ethics guidelines** if applicable

## Questions?

For data-related questions:
- **General**: Open an issue on GitHub
- **CDC data**: Contact CDC directly or check their documentation
- **Data processing**: See `docs/WORKFLOW.md` or email g.dewey@northeastern.edu

---

**Last Updated**: 2025-12-18
