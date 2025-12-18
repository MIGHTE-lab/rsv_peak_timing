# Repository File Index

Complete listing of all files in this repository with descriptions.

## Root Directory

| File | Purpose |
|------|---------|
| `README.md` | Main repository documentation with overview, installation, and usage |
| `LICENSE` | MIT License for code distribution |
| `CITATION.cff` | Citation information in CFF format for easy referencing |
| `CONTRIBUTING.md` | Guidelines for contributing to the project |
| `.gitignore` | Specifies files to exclude from version control |
| `run_complete_analysis.R` | Master script to run entire analysis pipeline |

## R/ - Analysis Scripts

### 02_modeling/

| Script | Purpose | Inputs | Outputs | Runtime |
|--------|---------|--------|---------|---------|
| `flu_rsv_ts_analysis.R` | Main ridge regression analysis decomposing ILI into pathogen contributions | NSSP data, ILINet data | Model coefficients, fitted values, R² statistics | ~15-20 min |
| `lagged_correlation_timing.R` | Computes cross-correlations between RSV and flu at different time lags | Rescaled time series from modeling | Correlation matrices, optimal lag identification | ~5 min |

**Key Features**:
- Ridge regression with cross-validation
- Positive coefficient constraints
- Season-specific models for each state
- Model diagnostics and validation

### 03_early_warning/

| Script | Purpose | Inputs | Outputs | Runtime |
|--------|---------|--------|---------|---------|
| `add_EWS_onset_and_peaks.R` | Detects epidemic onsets and peaks using anomaly detection | Processed time series | Onset/peak dates by pathogen and state-season | ~5-10 min |
| `combine_ews_output.R` | Consolidates EWS outputs across pathogens | Individual EWS files | Combined timing dataset with time differences | <1 min |

**Key Features**:
- Exponential growth detection algorithm
- Onset identification with 6-week lookback
- Peak detection (global maximum)
- Cross-pathogen timing comparisons

### 04_visualization/

| Script | Purpose | Inputs | Outputs | Runtime |
|--------|---------|--------|---------|---------|
| `fig1_plotting.R` | Creates Figure 1 showing epidemic trajectories | Time series data, correlations | Figure 1 (3×3 grid of example states) | ~2 min |
| `onset_peak_scatterplot.R` | Generates onset/peak timing scatterplots | EWS outputs | Figure 2 (timing comparison plots) | ~3 min |
| `state_timeline_viz.r` | Creates heatmap of peak timing across states | Combined EWS data | Timeline visualization, summary statistics | ~2 min |
| `create_coefficient_comparison_scatterplots.r` | Compares coefficients across modeling approaches | Multiple model results | Comparison scatterplots | ~2 min |

**Key Features**:
- Publication-ready figures
- Consistent color schemes
- State labeling and annotations
- Multiple output formats (PNG, PDF)

### 05_tables/

| Script | Purpose | Inputs | Outputs | Runtime |
|--------|---------|--------|---------|---------|
| `onset_peak_statistics.R` | Calculates comprehensive timing statistics | EWS data | Summary statistics, formatted output | <1 min |
| `create_results_table.r` | Formats model results into Excel tables | Model output files | Excel workbook with multiple sheets | ~1 min |
| `create_unified_comparison_table.r` | Creates comparison across modeling approaches | All model results | Unified comparison table | ~1 min |

**Key Features**:
- Excel output with formatting
- Summary statistics by season and state
- Model comparison metrics
- Data quality indicators

### utils/

| Script | Purpose | Status |
|--------|---------|--------|
| `test_spatial_correlation.R` | Exploratory spatial analysis of timing patterns | Supplementary/Optional |

## docs/ - Documentation

| File | Purpose |
|------|---------|
| `WORKFLOW.md` | Detailed explanation of analysis workflow and methods |
| `CODE_GUIDE.md` | Comprehensive guide to each script with technical details |
| `INSTALLATION.md` | Step-by-step installation and setup instructions |
| `session_info_last_run.txt` | R session information from last complete run (auto-generated) |

## data/ - Data Files (Not Included)

See `data/README.md` for complete data documentation.

### Required Structure:

```
data/
├── NSSP/                          # Raw NSSP data (user provided)
│   └── nssp_full_YYMMDD.csv
├── ILInet/                        # Raw ILINet data (user provided)
│   └── ILINet_state_YYYY_MM-DD.csv
├── processed/                     # Intermediate files (generated)
│   └── nssp_all_years_YYMMDD.csv
├── early_warning/                 # EWS outputs (generated)
│   ├── flu_matched_peaks_onsets.csv
│   ├── rsv_matched_peaks_onsets.csv
│   └── covid_matched_peaks_onsets.csv
└── results/                       # Model results (generated)
    ├── results_ridge_no_intercept.Rdata
    ├── results_ridge_with_intercept.Rdata
    └── results_linear.Rdata
```

## figures/ - Generated Figures (Not Included)

| Figure | Script | Description |
|--------|--------|-------------|
| `figure1_epidemic_trajectories.png` | `fig1_plotting.R` | 3×3 grid showing ILI, RSV, flu, COVID-19 for example states |
| `figure2_onset_peak_scatterplots.png` | `onset_peak_scatterplot.R` | Scatterplots comparing timing between pathogens |
| `state_timeline_peak_visualization.png` | `state_timeline_viz.r` | Heatmap of peak timing across all states and seasons |
| `coefficient_comparisons/*.png` | `create_coefficient_comparison_scatterplots.r` | Multiple comparison plots |

## tables/ - Generated Tables (Not Included)

| Table | Script | Description |
|-------|--------|-------------|
| `model_results_formatted.xlsx` | `create_results_table.r` | Model coefficients and fit statistics |
| `unified_model_comparison.xlsx` | `create_unified_comparison_table.r` | Cross-model comparison |
| `timing_statistics_summary.csv` | `onset_peak_statistics.R` | Summary statistics |
| `state_peak_timing_summary.csv` | `state_timeline_viz.r` | State-level timing patterns |
| `seasonal_peak_timing_summary.csv` | `state_timeline_viz.r` | Season-level aggregates |

## File Dependencies

### Execution Order:

```
1. Data preparation (user)
   ↓
2. flu_rsv_ts_analysis.R
   ↓
   ├─→ 3a. lagged_correlation_timing.R (parallel)
   │
   ├─→ 3b. add_EWS_onset_and_peaks.R
   │        ↓
   │   4. combine_ews_output.R
   │        ↓
   ├─→ 5a. fig1_plotting.R (parallel)
   ├─→ 5b. onset_peak_scatterplot.R (parallel)
   ├─→ 5c. state_timeline_viz.r (parallel)
   └─→ 5d. create_coefficient_comparison_scatterplots.r (parallel)
        ↓
   6a. onset_peak_statistics.R (parallel)
   6b. create_results_table.r (parallel)
   6c. create_unified_comparison_table.r (parallel)
```

### Key Dependencies:

**Core R Packages**:
- `tidyverse` - Data manipulation and visualization
- `glmnet` - Ridge regression
- `lubridate` - Date handling
- `patchwork` - Multi-panel figures

**Data Files**:
- NSSP emergency department data → All scripts
- ILINet surveillance data → Modeling scripts
- Model results → Visualization and table scripts
- EWS outputs → Timing analysis scripts

## File Counts

| Category | Count | Notes |
|----------|-------|-------|
| R Scripts | 11 | Core analysis + 1 utility |
| Documentation | 4 | Markdown files |
| Configuration | 3 | LICENSE, CITATION, .gitignore |
| Root Scripts | 1 | Master pipeline runner |
| **Total Code Files** | **19** | |
| Generated Outputs | Variable | Figures, tables, processed data |

## Size Estimates

| Category | Typical Size |
|----------|-------------|
| R Scripts | 5-30 KB each |
| Documentation | 10-50 KB each |
| Raw Data (user) | 500 MB - 1 GB |
| Processed Data | 50-100 MB |
| Model Results | 5-20 MB |
| Figures | 1-5 MB total |
| Tables | 1-5 MB total |
| **Total Repository** | ~20-50 KB (code only) |
| **With User Data** | ~1-2 GB |

## Quick Reference

### To Run Complete Analysis:
```r
source("run_complete_analysis.R")
```

### To Run Individual Stages:

**Modeling**:
```r
source("R/02_modeling/flu_rsv_ts_analysis.R")
```

**EWS Detection**:
```r
source("R/03_early_warning/add_EWS_onset_and_peaks.R")
source("R/03_early_warning/combine_ews_output.R")
```

**Figures**:
```r
source("R/04_visualization/fig1_plotting.R")
source("R/04_visualization/onset_peak_scatterplot.R")
```

**Tables**:
```r
source("R/05_tables/create_results_table.r")
```

### To Get Help:

1. **Installation**: See `docs/INSTALLATION.md`
2. **Workflow**: See `docs/WORKFLOW.md`
3. **Code Details**: See `docs/CODE_GUIDE.md`
4. **Data Setup**: See `data/README.md`
5. **Contributing**: See `CONTRIBUTING.md`

### Common File Locations:

- **Scripts**: `R/[stage]/script_name.R`
- **Data**: `data/[category]/file.csv`
- **Outputs**: `figures/` or `tables/`
- **Docs**: `docs/GUIDE_NAME.md`
- **Help**: `README.md` or GitHub Issues

## Version Information

**Last Updated**: 2025-12-18  
**Repository Version**: 1.0.0  
**Manuscript Version**: 5.6 clean  

## Maintainers

- George Dewey (g.dewey@northeastern.edu) - Lead Developer
- Mauricio Santillana (m.santill@g.harvard.edu) - Principal Investigator

## Additional Resources

- **GitHub Repository**: https://github.com/MIGHTE-lab/rsv_peak_timing
- **Paper Preprint**: [To be added]
- **Published Article**: [To be added upon acceptance]
- **CDC Data Sources**: https://data.cdc.gov/

---

*This index is automatically maintained. Last generated: 2025-12-18*
