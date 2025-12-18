# Quick Start Guide

Get up and running with the RSV/Flu timing analysis in 15 minutes.

## Prerequisites

- ✅ R version 4.0+ installed
- ✅ RStudio (recommended)
- ✅ ~2 GB free disk space
- ✅ CDC data downloaded (see below)

## 5-Step Setup

### Step 1: Clone Repository (2 minutes)

```bash
git clone https://github.com/MIGHTE-lab/rsv_peak_timing.git
cd rsv_peak_timing
```

Or download ZIP and extract.

### Step 2: Install Packages (5 minutes)

Open R/RStudio and run:

```r
# Install all required packages
required_packages = c(
  "tidyverse", "zoo", "lubridate", "MMWRweek", "ggh4x",
  "glmnet", "penalized", "caret", "patchwork", "ggrepel",
  "viridis", "RColorBrewer", "maps", "openxlsx", "knitr", "kableExtra"
)

new_packages = required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages, dependencies = TRUE)

# Verify
sapply(required_packages, require, character.only = TRUE)
```

### Step 3: Download Data (5 minutes)

**ILINet Data**:
1. Visit: https://www.cdc.gov/flu/weekly/
2. Download → State-level ILI data (2022-2025)
3. Save to `data/ILInet/ILINet_state_YYYY_MM-DD.csv`

**NSSP Data**:
1. Visit: https://data.cdc.gov/
2. Search "Emergency Department Visits" (COVID, Flu, RSV)
3. Download state-level weekly data
4. Save to `data/NSSP/nssp_full_YYMMDD.csv`

Or use example data (if provided):
```r
# Copy example files to data directories
# (Example files to be provided separately)
```

### Step 4: Update File Paths (2 minutes)

Edit `R/02_modeling/flu_rsv_ts_analysis.R`:

```r
# Lines ~37-38, update these paths:
dat = read_csv("data/NSSP/nssp_full_[YOUR_DATE].csv")
ili_data = read_csv("data/ILInet/ILINet_state_[YOUR_DATE].csv", skip = 1)
```

### Step 5: Run Analysis (1 minute to start, ~30 min to complete)

```r
# Set working directory
setwd("path/to/rsv_peak_timing")

# Run complete pipeline
source("run_complete_analysis.R")
```

Or run stages individually:
```r
# Just modeling
source("R/02_modeling/flu_rsv_ts_analysis.R")

# Just figures
source("R/04_visualization/fig1_plotting.R")
```

## What You'll Get

After running the complete pipeline:

**📊 Figures** (`figures/` directory):
- Figure 1: Epidemic trajectories for example states
- Figure 2: Onset/peak timing scatterplots
- State timeline visualization
- Coefficient comparison plots

**📋 Tables** (`tables/` directory):
- Model results (Excel format)
- Timing statistics
- Cross-model comparisons

**💾 Data** (`data/` directory):
- Processed time series
- Model coefficients
- Onset/peak dates
- Statistical summaries

## Quick Commands

### Check Installation
```r
source("run_complete_analysis.R")  # Will check packages first
```

### Run Specific Analysis
```r
# Modeling only
source("R/02_modeling/flu_rsv_ts_analysis.R")

# Early warning detection
source("R/03_early_warning/add_EWS_onset_and_peaks.R")

# Generate all figures
source("R/04_visualization/fig1_plotting.R")
source("R/04_visualization/onset_peak_scatterplot.R")
source("R/04_visualization/state_timeline_viz.r")

# Create tables
source("R/05_tables/onset_peak_statistics.R")
```

### View Outputs
```r
# List generated figures
list.files("figures", pattern = "\\.png$")

# List generated tables
list.files("tables", pattern = "\\.xlsx$")

# Load model results
load("data/results/results_ridge_no_intercept.Rdata")
```

### Get Help
```r
# View script documentation
?tidyverse  # Package help
file.edit("docs/WORKFLOW.md")  # Detailed workflow
file.edit("docs/CODE_GUIDE.md")  # Script details
```

## Common Issues

### "Cannot find file"
```r
# Check working directory
getwd()

# Should end with "rsv_peak_timing"
# If not:
setwd("path/to/rsv_peak_timing")
```

### "Package not found"
```r
# Install missing package
install.packages("package_name")

# Or reinstall all
source("docs/INSTALLATION.md")  # Follow package installation section
```

### "Error in ridge regression"
```r
# Check data format
dat = read_csv("data/NSSP/nssp_full_YYMMDD.csv")
head(dat)  # Should have columns: week_end, geography, percent_visits_*

# Verify date range
range(dat$week_end)  # Should cover 2022-2025
```

### "Figure not generated"
```r
# Check dependencies
# Figures require model outputs to exist
file.exists("data/results/results_ridge_no_intercept.Rdata")

# If FALSE, run modeling first
source("R/02_modeling/flu_rsv_ts_analysis.R")
```

## Next Steps

Once you have results:

1. **Explore outputs**: Check `figures/` and `tables/` directories
2. **Customize analysis**: Modify scripts for your needs
3. **Read documentation**: See `docs/WORKFLOW.md` for details
4. **Contribute**: See `CONTRIBUTING.md` to improve the code

## Time Estimates

| Task | Duration |
|------|----------|
| Setup (Steps 1-4) | ~15 minutes |
| Data processing | ~5 minutes |
| Modeling (all states) | ~15-20 minutes |
| EWS detection | ~5-10 minutes |
| Visualization | ~5 minutes |
| Tables | ~2 minutes |
| **Total** | **~45 minutes** |

*Times are estimates on standard laptop (8GB RAM, modern processor)*

## System Requirements

**Minimum**:
- R 4.0+
- 4 GB RAM
- 2 GB disk space

**Recommended**:
- R 4.3+
- 8+ GB RAM
- 5 GB disk space
- RStudio

## Getting Help

**Documentation**:
- 📖 Full workflow: `docs/WORKFLOW.md`
- 💻 Code details: `docs/CODE_GUIDE.md`
- ⚙️ Installation: `docs/INSTALLATION.md`
- 📂 File index: `INDEX.md`

**Support**:
- 🐛 Report issues: GitHub Issues
- 📧 Email: g.dewey@northeastern.edu
- 💬 Discussions: GitHub Discussions (if enabled)

## Resources

- **Paper**: [To be added upon publication]
- **CDC Data**: https://data.cdc.gov/
- **ILINet**: https://www.cdc.gov/flu/weekly/
- **R Documentation**: https://www.rdocumentation.org/

---

**Ready to start?** Run this in R:

```r
source("run_complete_analysis.R")
```

Good luck! 🚀
