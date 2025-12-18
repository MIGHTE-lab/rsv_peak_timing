## -----------------------------------------------------------------------------
## Script name: run_complete_analysis.R
##
## Purpose: Master script to run the complete RSV/Flu timing analysis pipeline
##
## Author: Analysis Team
##
## Date Created: 2025-12-18
##
## Usage: source("run_complete_analysis.R")
##
## Note: Ensure all data files are in place before running (see data/README.md)
## -----------------------------------------------------------------------------

## Setup -----------------------------------------------------------------------

cat("\n╔═══════════════════════════════════════════════════════════════════╗\n")
cat("║   RSV/Flu Timing Analysis - Complete Pipeline                    ║\n")
cat("╚═══════════════════════════════════════════════════════════════════╝\n\n")

# Record start time
start_time = Sys.time()

# Check if working directory is correct
if (!file.exists("R/02_modeling/flu_rsv_ts_analysis.R")) {
  stop("ERROR: Please set working directory to project root.\n",
       "Use: setwd('path/to/rsv_flu_timing_analysis')")
}

# Create output directories if they don't exist
dir.create("figures", showWarnings = FALSE, recursive = TRUE)
dir.create("tables", showWarnings = FALSE, recursive = TRUE)
dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)
dir.create("data/early_warning", showWarnings = FALSE, recursive = TRUE)
dir.create("data/results", showWarnings = FALSE, recursive = TRUE)

cat("✓ Output directories verified\n\n")

## Check Required Packages -----------------------------------------------------

cat("Checking required packages...\n")

required_packages = c(
  "tidyverse", "zoo", "lubridate", "MMWRweek", "ggh4x",
  "glmnet", "penalized", "caret", "patchwork", "ggrepel",
  "viridis", "RColorBrewer", "maps", "openxlsx", "knitr", "kableExtra"
)

missing_packages = required_packages[!(required_packages %in% installed.packages()[,"Package"])]

if (length(missing_packages) > 0) {
  cat("ERROR: Missing required packages:", paste(missing_packages, collapse = ", "), "\n")
  cat("Install with: install.packages(c('", paste(missing_packages, collapse = "', '"), "'))\n", sep = "")
  stop("Please install missing packages before continuing.")
}

cat("✓ All required packages installed\n\n")

## Stage 1: Data Processing ----------------------------------------------------

cat("═══════════════════════════════════════════════════════════════════\n")
cat("STAGE 1: Data Processing\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

# Note: User must implement data loading with their local data files
# The flu_rsv_ts_analysis.R script includes data loading code

cat("NOTE: Data loading is included in the modeling script.\n")
cat("      Ensure data files are in place (see data/README.md)\n\n")

## Stage 2: Statistical Modeling -----------------------------------------------

cat("═══════════════════════════════════════════════════════════════════\n")
cat("STAGE 2: Statistical Modeling\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

cat("Running ridge regression analysis...\n")
stage2_start = Sys.time()

tryCatch({
  source("R/02_modeling/flu_rsv_ts_analysis.R")
  cat("✓ Ridge regression completed\n")
}, error = function(e) {
  cat("✗ Error in ridge regression:\n")
  cat(conditionMessage(e), "\n")
  stop("Pipeline stopped due to error in modeling")
})

cat("\nRunning lagged correlation analysis...\n")
tryCatch({
  source("R/02_modeling/lagged_correlation_timing.R")
  cat("✓ Lagged correlation analysis completed\n")
}, error = function(e) {
  cat("✗ Error in lagged correlation:\n")
  cat(conditionMessage(e), "\n")
  cat("Warning: Continuing with pipeline...\n")
})

stage2_time = difftime(Sys.time(), stage2_start, units = "mins")
cat(sprintf("\n   Stage 2 completed in %.1f minutes\n\n", stage2_time))

## Stage 3: Early Warning Detection --------------------------------------------

cat("═══════════════════════════════════════════════════════════════════\n")
cat("STAGE 3: Early Warning Detection\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

stage3_start = Sys.time()

cat("Detecting epidemic onsets and peaks...\n")
tryCatch({
  source("R/03_early_warning/add_EWS_onset_and_peaks.R")
  cat("✓ Onset and peak detection completed\n")
}, error = function(e) {
  cat("✗ Error in EWS detection:\n")
  cat(conditionMessage(e), "\n")
  stop("Pipeline stopped due to error in EWS detection")
})

cat("\nCombining EWS outputs...\n")
tryCatch({
  source("R/03_early_warning/combine_ews_output.R")
  cat("✓ EWS outputs combined\n")
}, error = function(e) {
  cat("✗ Error combining EWS outputs:\n")
  cat(conditionMessage(e), "\n")
  cat("Warning: Continuing with pipeline...\n")
})

stage3_time = difftime(Sys.time(), stage3_start, units = "mins")
cat(sprintf("\n   Stage 3 completed in %.1f minutes\n\n", stage3_time))

## Stage 4: Visualization ------------------------------------------------------

cat("═══════════════════════════════════════════════════════════════════\n")
cat("STAGE 4: Visualization\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

stage4_start = Sys.time()

cat("Generating Figure 1 (epidemic trajectories)...\n")
tryCatch({
  source("R/04_visualization/fig1_plotting.R")
  cat("✓ Figure 1 created\n")
}, error = function(e) {
  cat("✗ Error creating Figure 1:\n")
  cat(conditionMessage(e), "\n")
  cat("Warning: Continuing with pipeline...\n")
})

cat("\nGenerating onset/peak scatterplots...\n")
tryCatch({
  source("R/04_visualization/onset_peak_scatterplot.R")
  cat("✓ Onset/peak scatterplots created\n")
}, error = function(e) {
  cat("✗ Error creating scatterplots:\n")
  cat(conditionMessage(e), "\n")
  cat("Warning: Continuing with pipeline...\n")
})

cat("\nGenerating state timeline visualization...\n")
tryCatch({
  source("R/04_visualization/state_timeline_viz.r")
  cat("✓ State timeline visualization created\n")
}, error = function(e) {
  cat("✗ Error creating timeline:\n")
  cat(conditionMessage(e), "\n")
  cat("Warning: Continuing with pipeline...\n")
})

cat("\nGenerating coefficient comparison plots...\n")
tryCatch({
  source("R/04_visualization/create_coefficient_comparison_scatterplots.r")
  cat("✓ Coefficient comparison plots created\n")
}, error = function(e) {
  cat("✗ Error creating coefficient plots:\n")
  cat(conditionMessage(e), "\n")
  cat("Warning: Continuing with pipeline...\n")
})

stage4_time = difftime(Sys.time(), stage4_start, units = "mins")
cat(sprintf("\n   Stage 4 completed in %.1f minutes\n\n", stage4_time))

## Stage 5: Tables and Statistics ----------------------------------------------

cat("═══════════════════════════════════════════════════════════════════\n")
cat("STAGE 5: Tables and Statistics\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

stage5_start = Sys.time()

cat("Calculating onset and peak statistics...\n")
tryCatch({
  source("R/05_tables/onset_peak_statistics.R")
  cat("✓ Statistics calculated\n")
}, error = function(e) {
  cat("✗ Error calculating statistics:\n")
  cat(conditionMessage(e), "\n")
  cat("Warning: Continuing with pipeline...\n")
})

cat("\nCreating results tables...\n")
tryCatch({
  source("R/05_tables/create_results_table.r")
  cat("✓ Results tables created\n")
}, error = function(e) {
  cat("✗ Error creating results tables:\n")
  cat(conditionMessage(e), "\n")
  cat("Warning: Continuing with pipeline...\n")
})

cat("\nCreating unified comparison table...\n")
tryCatch({
  source("R/05_tables/create_unified_comparison_table.r")
  cat("✓ Unified comparison table created\n")
}, error = function(e) {
  cat("✗ Error creating comparison table:\n")
  cat(conditionMessage(e), "\n")
  cat("Warning: Continuing with pipeline...\n")
})

stage5_time = difftime(Sys.time(), stage5_start, units = "mins")
cat(sprintf("\n   Stage 5 completed in %.1f minutes\n\n", stage5_time))

## Summary ---------------------------------------------------------------------

cat("═══════════════════════════════════════════════════════════════════\n")
cat("ANALYSIS COMPLETE\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

total_time = difftime(Sys.time(), start_time, units = "mins")

cat("Time Summary:\n")
cat(sprintf("  Stage 2 (Modeling):        %.1f minutes\n", stage2_time))
cat(sprintf("  Stage 3 (EWS Detection):   %.1f minutes\n", stage3_time))
cat(sprintf("  Stage 4 (Visualization):   %.1f minutes\n", stage4_time))
cat(sprintf("  Stage 5 (Tables):          %.1f minutes\n", stage5_time))
cat(sprintf("  ────────────────────────────────────────\n"))
cat(sprintf("  Total Time:                %.1f minutes\n\n", total_time))

cat("Output Locations:\n")
cat("  Figures:        figures/\n")
cat("  Tables:         tables/\n")
cat("  Processed Data: data/processed/\n")
cat("  Model Results:  data/results/\n\n")

# Save session info
cat("Saving session information...\n")
writeLines(capture.output(sessionInfo()), "docs/session_info_last_run.txt")
cat("✓ Session info saved to docs/session_info_last_run.txt\n\n")

cat("✓ Analysis pipeline completed successfully!\n\n")

# Optional: Open results directory
if (interactive()) {
  response = readline(prompt = "Open results directory? (y/n): ")
  if (tolower(response) == "y") {
    if (.Platform$OS.type == "windows") {
      shell.exec("figures")
    } else {
      system("open figures")
    }
  }
}

## -----------------------------------------------------------------------------
## End of script
## -----------------------------------------------------------------------------
