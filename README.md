# Post-Pandemic Timing of Influenza, RSV, and COVID-19 in the United States

[![DOI](https://img.shields.io/badge/DOI-pending-blue)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

This repository contains the code and analysis pipeline for the manuscript:

**"Uncovering the Post-Pandemic Timing of Influenza, RSV, and COVID-19 Driving Seasonal Influenza-like Illness in the United States: A Retrospective Ecological Study"**

Authors: George Dewey, Austin G. Meyer, Raul Garrido Garcia, Mauricio Santillana

## Overview

This project analyzes the temporal relationships between RSV, influenza, and COVID-19 epidemics in the post-pandemic United States (2022-2025). Using ridge regression models and anomaly detection algorithms applied to CDC surveillance data, we characterize when these epidemics peak relative to each other and their contribution to seasonal influenza-like illness (ILI).

### Key Findings

- In 77% of analyzed state-seasons (114/148), RSV peaks occurred before influenza peaks
- Median time difference between RSV and influenza peaks: +3.0 weeks (IQR: 5.0 weeks)
- RSV and influenza onset timing showed increasing synchrony in recent seasons
- COVID-19 outbreak timing showed no consistent seasonal pattern

## Repository Structure

```
rsv_flu_timing_analysis/
├── R/                          # Analysis scripts
│   ├── 01_data_processing/    # Data loading and preprocessing
│   ├── 02_modeling/           # Ridge regression and time series analysis
│   ├── 03_early_warning/      # Onset and peak detection algorithms
│   ├── 04_visualization/      # Figure generation scripts
│   └── 05_tables/             # Results table creation
├── data/                       # Data directory (not included in repo)
│   ├── raw/                   # Raw CDC data
│   ├── processed/             # Processed datasets
│   └── results/               # Model outputs
├── figures/                    # Generated figures
├── docs/                       # Additional documentation
└── README.md                   # This file
```

## Data Sources

All data are publicly available from the U.S. Centers for Disease Control and Prevention:

1. **ILINet (Influenza-like Illness Surveillance Network)**
   - Source: https://www.cdc.gov/flu/weekly/
   - Weekly state-level ILI surveillance data

2. **NSSP (National Syndromic Surveillance Program)**
   - Source: https://data.cdc.gov/
   - Emergency department visit data for COVID-19, influenza, and RSV

**Note:** Due to data usage agreements and file size, raw data files are not included in this repository. Users must download data directly from CDC sources using the links above.

## Installation

### System Requirements

- R version ≥ 4.0.0
- RStudio (recommended)
- Operating System: Windows, macOS, or Linux

### Required R Packages

```r
# Data manipulation
install.packages(c("tidyverse", "zoo", "lubridate", "MMWRweek"))

# Modeling
install.packages(c("glmnet", "penalized", "caret"))

# Visualization
install.packages(c("ggplot2", "patchwork", "ggh4x", "ggrepel", 
                   "viridis", "RColorBrewer", "envalysis"))

# Spatial analysis
install.packages(c("maps", "sf"))

# Tables and output
install.packages(c("openxlsx", "knitr", "kableExtra"))
```

## Usage

### Quick Start

The analysis pipeline consists of five main stages:

```r
# Set your working directory
setwd("path/to/rsv_flu_timing_analysis")

# 1. Process and combine data sources
source("R/01_data_processing/data_loading_preprocessing.R")

# 2. Run ridge regression models
source("R/02_modeling/flu_rsv_ts_analysis.R")
source("R/02_modeling/lagged_correlation_timing.R")

# 3. Detect epidemic onsets and peaks
source("R/03_early_warning/add_EWS_onset_and_peaks.R")
source("R/03_early_warning/combine_ews_output.R")

# 4. Generate figures
source("R/04_visualization/fig1_plotting.R")
source("R/04_visualization/onset_peak_scatterplot.R")
source("R/04_visualization/state_timeline_viz.R")

# 5. Create results tables
source("R/05_tables/create_results_table.R")
source("R/05_tables/create_unified_comparison_table.R")
```

### Detailed Workflow

See [docs/WORKFLOW.md](docs/WORKFLOW.md) for a detailed description of each analysis step.

## Script Descriptions

### 01_data_processing/
- **Data loading and preprocessing scripts** (to be organized from existing code)

### 02_modeling/
- **`flu_rsv_ts_analysis.R`**: Main ridge regression analysis modeling ILI as a function of RSV, influenza, and COVID-19
- **`lagged_correlation_timing.R`**: Compute lagged cross-correlations between RSV and influenza time series

### 03_early_warning/
- **`add_EWS_onset_and_peaks.R`**: Apply anomaly detection algorithms to identify epidemic onsets and peaks
- **`combine_ews_output.R`**: Consolidate early warning system outputs across seasons

### 04_visualization/
- **`fig1_plotting.R`**: Generate Figure 1 (epidemic trajectories for example states)
- **`onset_peak_scatterplot.R`**: Create scatterplots showing timing relationships between epidemics
- **`state_timeline_viz.R`**: Visualize state-by-state peak timing patterns
- **`create_coefficient_comparison_scatterplots.R`**: Compare model coefficients across approaches

### 05_tables/
- **`create_results_table.R`**: Format model results into publication-ready tables
- **`create_unified_comparison_table.R`**: Compare results across different modeling approaches
- **`onset_peak_statistics.R`**: Calculate summary statistics for onset and peak timing

### Additional Scripts/
- **`test_spatial_correlation.R`**: Exploratory spatial analysis of epidemic timing patterns

## Reproducibility

### Software Versions

This analysis was conducted using:
- R version 4.3.x
- Key package versions documented in `docs/session_info.txt`

### Random Seeds

Where applicable, random seeds are set within scripts for reproducibility of model fitting procedures.

### Expected Runtime

- Data processing: ~5 minutes
- Model fitting: ~15-20 minutes (all state-seasons)
- Figure generation: ~5 minutes
- Total: ~30 minutes on standard laptop

## Citation

If you use this code or data in your research, please cite:

```
Dewey G, Meyer AG, Garcia RG, Santillana M. Uncovering the Post-Pandemic 
Timing of Influenza, RSV, and COVID-19 Driving Seasonal Influenza-like 
Illness in the United States: A Retrospective Ecological Study. 
[Journal Name]. [Year];[Volume]:[Pages]. DOI: [pending]
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

For questions or issues, please contact:
- George Dewey: g.dewey@northeastern.edu
- Mauricio Santillana: m.santill@g.harvard.edu

Or open an issue on this repository.

## Funding

GD, RG, and MS were supported by cooperative agreement CDC-RFA-FT-23-0069 from the CDCs Center for Forecasting and Outbreak Analytics. AM was supported by the National Institutes of Health LRP #1L70AI194328-01. The study contents are solely the responsibility of the authors and do not necessarily represent the official views of the Centers for Disease Control and Prevention.

## Acknowledgments

We thank Marc Lipsitch for helpful comments. 

## Repository Status

✅ Code reviewed and documented
✅ Ready for public release upon manuscript acceptance
🔄 Data: Users must download from CDC sources (see Data Sources section)
