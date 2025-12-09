## -----------------------------------------------------------------------------
## Script name: create_unified_comparison_table.R
##
## Purpose of script: Create unified comparison table across three modeling approaches
##
## Author: George Dewey
##
## Date Created: 2025-10-01
##
## Last Updated: 2025-10-01
## -----------------------------------------------------------------------------

## Load libraries
library(tidyverse)
library(openxlsx)

## Load the model outputs
load("data/results_ridge_no_intercept.Rdata")
load("data/results_ridge_with_intercept.Rdata")
load("data/results_linear.Rdata")

## Function to create unified comparison table
create_unified_comparison_table = function(results_ridge_no_int,
                                           results_ridge_with_int,
                                           results_linear,
                                           output_file = "unified_model_comparison.xlsx",
                                           version = "full") {

  # Validate version parameter
  if (!version %in% c("full", "simplified")) {
    stop("version must be either 'full' or 'simplified'")
  }

  # Extract model results from each approach
  ridge_no_int = results_ridge_no_int$model_results %>%
    mutate(model_approach = "Ridge No Intercept")

  ridge_with_int = results_ridge_with_int$model_results %>%
    mutate(model_approach = "Ridge With Intercept")

  linear_model = results_linear$model_results %>%
    mutate(model_approach = "Linear No Regularization")

  # Create state-season identifier for joining
  add_state_season_id = function(df) {
    df %>% mutate(state_season = paste(state, season, sep = "_"))
  }

  ridge_no_int = add_state_season_id(ridge_no_int)
  ridge_with_int = add_state_season_id(ridge_with_int)
  linear_model = add_state_season_id(linear_model)

  # Define column sets for each version
  if (version == "simplified") {
    # Simplified version: only basic info and beta coefficients
    select_cols = function(df, suffix) {
      df %>% select(
        state_season,
        state,
        season,
        n_weeks,
        beta_intercept,
        beta_covid,
        beta_flu,
        beta_rsv
      ) %>%
        rename_with(~ paste0(.x, "_", suffix), .cols = -c(state_season, state, season))
    }

  } else {
    # Full version: all columns
    select_cols = function(df, suffix) {
      df %>% select(
        state_season,
        state,
        season,
        n_weeks,
        beta_intercept,
        beta_covid,
        beta_flu,
        beta_rsv,
        flu_mean_contribution,
        flu_max_contribution,
        flu_sum_contribution,
        flu_peak_week,
        rsv_mean_contribution,
        rsv_max_contribution,
        rsv_sum_contribution,
        rsv_peak_week,
        covid_mean_contribution,
        covid_max_contribution,
        covid_sum_contribution,
        covid_peak_week,
        total_explained_variance,
        dominant_pathogen
      ) %>%
        rename_with(~ paste0(.x, "_", suffix), .cols = -c(state_season, state, season))
    }
  }

  # Apply column selection with appropriate suffixes
  ridge_no_int_selected = select_cols(ridge_no_int, "ridge_no_int")
  ridge_with_int_selected = select_cols(ridge_with_int, "ridge_with_int")
  linear_selected = select_cols(linear_model, "linear")

  # Join all three datasets
  unified_table = ridge_no_int_selected %>%
    full_join(ridge_with_int_selected, by = c("state_season", "state", "season")) %>%
    full_join(linear_selected, by = c("state_season", "state", "season")) %>%
    select(-state_season) %>%
    arrange(season, state)

  # Round numeric columns for better presentation
  unified_table = unified_table %>%
    mutate(across(where(is.numeric), ~ round(.x, 3)))

  # Create better column names for display
  if (version == "simplified") {
    display_names = c(
      "State", "Season",
      # Ridge No Intercept block
      "N Weeks", "β Intercept", "β COVID", "β Flu", "β RSV",
      # Ridge With Intercept block
      "N Weeks", "β Intercept", "β COVID", "β Flu", "β RSV",
      # Linear block
      "N Weeks", "β Intercept", "β COVID", "β Flu", "β RSV"
    )
  } else {
    display_names = c(
      "State", "Season",
      # Ridge No Intercept block (21 columns)
      "N Weeks", "β Intercept", "β COVID", "β Flu", "β RSV",
      "Flu Mean", "Flu Max", "Flu Sum", "Flu Peak Wk",
      "RSV Mean", "RSV Max", "RSV Sum", "RSV Peak Wk",
      "COVID Mean", "COVID Max", "COVID Sum", "COVID Peak Wk",
      "Total Explained", "Dominant Path",
      # Ridge With Intercept block (21 columns)
      "N Weeks", "β Intercept", "β COVID", "β Flu", "β RSV",
      "Flu Mean", "Flu Max", "Flu Sum", "Flu Peak Wk",
      "RSV Mean", "RSV Max", "RSV Sum", "RSV Peak Wk",
      "COVID Mean", "COVID Max", "COVID Sum", "COVID Peak Wk",
      "Total Explained", "Dominant Path",
      # Linear block (21 columns)
      "N Weeks", "β Intercept", "β COVID", "β Flu", "β RSV",
      "Flu Mean", "Flu Max", "Flu Sum", "Flu Peak Wk",
      "RSV Mean", "RSV Max", "RSV Sum", "RSV Peak Wk",
      "COVID Mean", "COVID Max", "COVID Sum", "COVID Peak Wk",
      "Total Explained", "Dominant Path"
    )
  }

  # Apply display names
  names(unified_table) = display_names

  # Create Excel workbook
  wb = createWorkbook()
  addWorksheet(wb, "Unified Comparison")
  writeData(wb, "Unified Comparison", unified_table, startRow = 2)

  # Add model headers
  if (version == "simplified") {
    model_headers = c("", "", rep("Ridge No Intercept", 5),
                      rep("Ridge With Intercept", 5),
                      rep("Linear No Regularization", 5))
    block_positions = list(
      basic = 1:2,
      ridge_no_int = 3:7,
      ridge_with_int = 8:12,
      linear = 13:17
    )
  } else {
    model_headers = c("", "", rep("Ridge No Intercept", 19),
                      rep("Ridge With Intercept", 19),
                      rep("Linear No Regularization", 19))
    block_positions = list(
      basic = 1:2,
      ridge_no_int = 3:21,
      ridge_with_int = 22:40,
      linear = 41:59
    )
  }

  # Write model headers
  writeData(wb, "Unified Comparison", t(model_headers), startRow = 1)

  # Formatting
  # Model header formatting
  addStyle(wb, "Unified Comparison",
           style = createStyle(
             fontSize = 12,
             fontColour = "white",
             fgFill = "#2F5597",
             halign = "center",
             valign = "center",
             textDecoration = "bold",
             border = "TopBottomLeftRight",
             borderColour = "white"
           ),
           rows = 1, cols = block_positions$ridge_no_int)

  addStyle(wb, "Unified Comparison",
           style = createStyle(
             fontSize = 12,
             fontColour = "white",
             fgFill = "#C5504B",
             halign = "center",
             valign = "center",
             textDecoration = "bold",
             border = "TopBottomLeftRight",
             borderColour = "white"
           ),
           rows = 1, cols = block_positions$ridge_with_int)

  addStyle(wb, "Unified Comparison",
           style = createStyle(
             fontSize = 12,
             fontColour = "white",
             fgFill = "#70AD47",
             halign = "center",
             valign = "center",
             textDecoration = "bold",
             border = "TopBottomLeftRight",
             borderColour = "white"
           ),
           rows = 1, cols = block_positions$linear)

  # Column header formatting
  addStyle(wb, "Unified Comparison",
           style = createStyle(
             fontSize = 10,
             fontColour = "black",
             fgFill = "#D9E2F3",
             halign = "center",
             valign = "center",
             textDecoration = "bold",
             border = "TopBottomLeftRight",
             borderColour = "#8EA9DB"
           ),
           rows = 2, cols = 1:ncol(unified_table))

  # Data formatting
  addStyle(wb, "Unified Comparison",
           style = createStyle(
             fontSize = 9,
             halign = "center",
             valign = "center",
             border = "TopBottomLeftRight",
             borderColour = "#D9D9D9"
           ),
           rows = 3:(nrow(unified_table) + 2),
           cols = 1:ncol(unified_table),
           gridExpand = TRUE)

  # Alternate row coloring
  addStyle(wb, "Unified Comparison",
           style = createStyle(fgFill = "#F8F9FA"),
           rows = seq(3, nrow(unified_table) + 2, by = 2),
           cols = 1:ncol(unified_table),
           gridExpand = TRUE)

  # Add subtle background colors for each model block
  addStyle(wb, "Unified Comparison",
           style = createStyle(fgFill = "#EBF1FF"),
           rows = 3:(nrow(unified_table) + 2),
           cols = block_positions$ridge_no_int,
           gridExpand = TRUE)

  addStyle(wb, "Unified Comparison",
           style = createStyle(fgFill = "#FCE4EC"),
           rows = 3:(nrow(unified_table) + 2),
           cols = block_positions$ridge_with_int,
           gridExpand = TRUE)

  addStyle(wb, "Unified Comparison",
           style = createStyle(fgFill = "#E8F5E8"),
           rows = 3:(nrow(unified_table) + 2),
           cols = block_positions$linear,
           gridExpand = TRUE)

  # Auto-size columns
  setColWidths(wb, "Unified Comparison", cols = 1:ncol(unified_table), widths = "auto")

  # Freeze panes (first two columns)
  freezePane(wb, "Unified Comparison", firstActiveRow = 3, firstActiveCol = 3)

  # Save workbook
  saveWorkbook(wb, output_file, overwrite = TRUE)

  # Print summary
  cat("\n=== UNIFIED COMPARISON TABLE CREATED ===\n")
  cat("Version:", version, "\n")
  cat("Output file:", output_file, "\n")
  cat("Dimensions:", nrow(unified_table), "rows ×", ncol(unified_table), "columns\n")
  cat("State-seasons included:", nrow(unified_table), "\n")

  if (version == "simplified") {
    cat("Columns per model: 5 (N, β_intercept, β_covid, β_flu, β_rsv)\n")
  } else {
    cat("Columns per model: 19 (full statistics)\n")
  }

  return(unified_table)
}

## Usage Examples =============================================================

# Example 1: Create simplified version
unified_simple = create_unified_comparison_table(
  results_ridge_no_int = results_all_original,
  results_ridge_with_int = results_with_intercept,
  results_linear = results_linear,
  output_file = "figures/unified_comparison_simplified_raw.xlsx",
  version = "simplified"
)

# Example 2: Create full version
# unified_full = create_unified_comparison_table(
#   results_ridge_no_int = results_original,
#   results_ridge_with_int = results_with_intercept,
#   results_linear = results_linear,
#   output_file = "unified_comparison_full.xlsx",
#   version = "full"
# )

cat("Unified comparison table function loaded successfully!\n")
cat("Use create_unified_comparison_table() with version = 'simplified' or 'full'\n")