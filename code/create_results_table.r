## -----------------------------------------------------------------------------
## Script name: create_excel_results_table.R
##
## Purpose of script: Create formatted Excel tables from model results
##
## Author: George Dewey
##
## Date Created: 2025-10-01
##
## Last Updated: 2025-10-01
## -----------------------------------------------------------------------------

# Load required libraries
library(tidyverse)
library(openxlsx)
library(knitr)
library(kableExtra)

## Function to create formatted Excel tables
create_results_excel = function(model_results,
                                output_file = "model_results_formatted.xlsx",
                                include_summary = TRUE) {

  # Create a new workbook
  wb = createWorkbook()

  ## 1. COMPREHENSIVE RESULTS TABLE ==========================================

  # Prepare comprehensive data with better column names
  comprehensive_data = model_results %>%
    arrange(season, state) %>%
    select(
      State = state,
      Season = season,
      `Model Type` = model_type,
      `N Weeks` = n_weeks,
      `β Intercept` = beta_intercept,
      `β COVID` = beta_covid,
      `β Flu` = beta_flu,
      `β RSV` = beta_rsv,
      `Flu Max Contrib` = flu_max_contribution,
      `RSV Max Contrib` = rsv_max_contribution,
      `COVID Max Contrib` = covid_max_contribution,
      `Flu Peak Week` = flu_peak_week,
      `RSV Peak Week` = rsv_peak_week,
      `COVID Peak Week` = covid_peak_week,
      `Total Explained Var` = total_explained_variance,
      `Dominant Pathogen` = dominant_pathogen
    ) %>%
    mutate(
      # Round numeric columns for better presentation
      across(starts_with("β"), ~ round(.x, 3)),
      across(contains("Contrib"), ~ round(.x, 3)),
      `Total Explained Var` = round(`Total Explained Var`, 3)
    )

  # Add comprehensive sheet
  addWorksheet(wb, "Comprehensive Results")
  writeData(wb, "Comprehensive Results", comprehensive_data, startRow = 1)

  # Format comprehensive sheet
  # Header formatting
  addStyle(wb, "Comprehensive Results",
           style = createStyle(
             fontSize = 11,
             fontColour = "white",
             fgFill = "#4472C4",
             halign = "center",
             valign = "center",
             textDecoration = "bold",
             border = "TopBottomLeftRight",
             borderColour = "white"
           ),
           rows = 1, cols = 1:ncol(comprehensive_data))

  # Data formatting
  addStyle(wb, "Comprehensive Results",
           style = createStyle(
             fontSize = 10,
             halign = "center",
             valign = "center",
             border = "TopBottomLeftRight",
             borderColour = "#D9D9D9"
           ),
           rows = 2:(nrow(comprehensive_data) + 1),
           cols = 1:ncol(comprehensive_data),
           gridExpand = TRUE)

  # Alternate row coloring
  addStyle(wb, "Comprehensive Results",
           style = createStyle(fgFill = "#F2F2F2"),
           rows = seq(2, nrow(comprehensive_data) + 1, by = 2),
           cols = 1:ncol(comprehensive_data),
           gridExpand = TRUE)

  # Auto-size columns
  setColWidths(wb, "Comprehensive Results", cols = 1:ncol(comprehensive_data), widths = "auto")

  ## 2. SUMMARY STATISTICS TABLE ==============================================

  if (include_summary) {
    # Create summary statistics
    summary_stats = model_results %>%
      group_by(season, model_type) %>%
      summarise(
        `N State-Seasons` = n(),
        `Mean β Intercept` = round(mean(beta_intercept, na.rm = TRUE), 3),
        `Mean β COVID` = round(mean(beta_covid, na.rm = TRUE), 3),
        `Mean β Flu` = round(mean(beta_flu, na.rm = TRUE), 3),
        `Mean β RSV` = round(mean(beta_rsv, na.rm = TRUE), 3),
        `Mean Flu Max Contrib` = round(mean(flu_max_contribution, na.rm = TRUE), 3),
        `Mean RSV Max Contrib` = round(mean(rsv_max_contribution, na.rm = TRUE), 3),
        `Mean COVID Max Contrib` = round(mean(covid_max_contribution, na.rm = TRUE), 3),
        `% Flu Dominant` = round(100 * sum(dominant_pathogen == "Flu") / n(), 1),
        `% RSV Dominant` = round(100 * sum(dominant_pathogen == "RSV") / n(), 1),
        `% COVID Dominant` = round(100 * sum(dominant_pathogen == "COVID") / n(), 1),
        .groups = "drop"
      ) %>%
      arrange(season)

    # Add summary sheet
    addWorksheet(wb, "Summary Statistics")
    writeData(wb, "Summary Statistics", summary_stats, startRow = 1)

    # Format summary sheet (similar to comprehensive)
    addStyle(wb, "Summary Statistics",
             style = createStyle(
               fontSize = 11,
               fontColour = "white",
               fgFill = "#70AD47",
               halign = "center",
               valign = "center",
               textDecoration = "bold",
               border = "TopBottomLeftRight",
               borderColour = "white"
             ),
             rows = 1, cols = 1:ncol(summary_stats))

    addStyle(wb, "Summary Statistics",
             style = createStyle(
               fontSize = 10,
               halign = "center",
               valign = "center",
               border = "TopBottomLeftRight",
               borderColour = "#D9D9D9"
             ),
             rows = 2:(nrow(summary_stats) + 1),
             cols = 1:ncol(summary_stats),
             gridExpand = TRUE)

    # Alternate row coloring
    addStyle(wb, "Summary Statistics",
             style = createStyle(fgFill = "#F2F2F2"),
             rows = seq(2, nrow(summary_stats) + 1, by = 2),
             cols = 1:ncol(summary_stats),
             gridExpand = TRUE)

    setColWidths(wb, "Summary Statistics", cols = 1:ncol(summary_stats), widths = "auto")
  }

  ## 3. COEFFICIENT COMPARISON TABLE =========================================

  # Create coefficient comparison table
  coef_comparison = model_results %>%
    select(state, season, model_type, beta_intercept, beta_covid, beta_flu, beta_rsv) %>%
    pivot_longer(cols = starts_with("beta_"),
                 names_to = "coefficient",
                 values_to = "value") %>%
    mutate(
      coefficient = case_when(
        coefficient == "beta_intercept" ~ "Intercept",
        coefficient == "beta_covid" ~ "COVID",
        coefficient == "beta_flu" ~ "Flu",
        coefficient == "beta_rsv" ~ "RSV"
      ),
      value = round(value, 3)
    ) %>%
    group_by(season, model_type, coefficient) %>%
    summarise(
      Mean = round(mean(value, na.rm = TRUE), 3),
      Median = round(median(value, na.rm = TRUE), 3),
      `Std Dev` = round(sd(value, na.rm = TRUE), 3),
      Min = round(min(value, na.rm = TRUE), 3),
      Max = round(max(value, na.rm = TRUE), 3),
      .groups = "drop"
    ) %>%
    arrange(season, coefficient)

  # Add coefficient comparison sheet
  addWorksheet(wb, "Coefficient Comparison")
  writeData(wb, "Coefficient Comparison", coef_comparison, startRow = 1)

  # Format coefficient comparison sheet
  addStyle(wb, "Coefficient Comparison",
           style = createStyle(
             fontSize = 11,
             fontColour = "white",
             fgFill = "#E7E6E6",
             halign = "center",
             valign = "center",
             textDecoration = "bold",
             border = "TopBottomLeftRight",
             borderColour = "black"
           ),
           rows = 1, cols = 1:ncol(coef_comparison))

  addStyle(wb, "Coefficient Comparison",
           style = createStyle(
             fontSize = 10,
             halign = "center",
             valign = "center",
             border = "TopBottomLeftRight",
             borderColour = "#D9D9D9"
           ),
           rows = 2:(nrow(coef_comparison) + 1),
           cols = 1:ncol(coef_comparison),
           gridExpand = TRUE)

  setColWidths(wb, "Coefficient Comparison", cols = 1:ncol(coef_comparison), widths = "auto")

  ## 4. SAVE WORKBOOK ========================================================

  saveWorkbook(wb, output_file, overwrite = TRUE)
  cat("Excel file saved as:", output_file, "\n")

  # Return the data frames for potential further use
  return(list(
    comprehensive = comprehensive_data,
    summary = if(include_summary) summary_stats else NULL,
    coefficients = coef_comparison
  ))
}

## Function to create Word-ready HTML table
create_word_ready_table = function(model_results,
                                   table_type = "summary",
                                   caption = "Model Results Summary") {

  if (table_type == "summary") {
    # Create summary table
    table_data = model_results %>%
      group_by(season) %>%
      summarise(
        `N States` = n(),
        `Mean β Flu` = round(mean(beta_flu, na.rm = TRUE), 3),
        `Mean β RSV` = round(mean(beta_rsv, na.rm = TRUE), 3),
        `Mean β COVID` = round(mean(beta_covid, na.rm = TRUE), 3),
        `% RSV Dominant` = round(100 * sum(dominant_pathogen == "RSV") / n(), 1),
        `% Flu Dominant` = round(100 * sum(dominant_pathogen == "Flu") / n(), 1),
        `% COVID Dominant` = round(100 * sum(dominant_pathogen == "COVID") / n(), 1),
        .groups = "drop"
      )

  } else if (table_type == "coefficients") {
    # Create coefficient summary table
    table_data = model_results %>%
      select(season, beta_flu, beta_rsv, beta_covid) %>%
      pivot_longer(cols = starts_with("beta_"),
                   names_to = "Pathogen",
                   values_to = "Coefficient") %>%
      mutate(
        Pathogen = case_when(
          Pathogen == "beta_flu" ~ "Influenza",
          Pathogen == "beta_rsv" ~ "RSV",
          Pathogen == "beta_covid" ~ "COVID-19"
        )
      ) %>%
      group_by(season, Pathogen) %>%
      summarise(
        Mean = round(mean(Coefficient, na.rm = TRUE), 3),
        `Std Dev` = round(sd(Coefficient, na.rm = TRUE), 3),
        Median = round(median(Coefficient, na.rm = TRUE), 3),
        .groups = "drop"
      ) %>%
      arrange(season, Pathogen)
  }

  # Create formatted HTML table
  html_table = table_data %>%
    kable("html",
          caption = caption,
          align = "c") %>%
    kable_styling(
      bootstrap_options = c("striped", "hover", "condensed", "responsive"),
      full_width = FALSE,
      position = "center"
    ) %>%
    row_spec(0, bold = TRUE, background = "#4472C4", color = "white") %>%
    column_spec(1, bold = TRUE, background = "#F8F9FA")

  return(html_table)
}

## Usage Examples =============================================================

# Example 1: Create comprehensive Excel file
excel_tables = create_results_excel(
  model_results = tbl_ridge_no_intercept,
  output_file = "comprehensive_model_results_ridge_no_intercept.xlsx",
  include_summary = TRUE
)

# Example 2: Create Word-ready summary table
word_summary_ridge_no_intercept = create_word_ready_table(
#   model_results = tbl_ridge_no_intercept,
#   table_type = "summary",
#   caption = "Table 1. Summary of Ridge Regression Model Results by Season"
# )
#
# # # Save as HTML file that can be opened in Word
# writeLines(as.character(word_summary_ridge_no_intercept), "summary_table_for_word_ridge_no_intercept.html")

# Example 3: Create Word-ready coefficient table
# word_coefficients = create_word_ready_table(
#   model_results = results_all_original$model_results,
#   table_type = "coefficients",
#   caption = "Table 2. Model Coefficient Statistics by Season and Pathogen"
# )
#
# writeLines(as.character(word_coefficients), "coefficient_table_for_word.html")

cat("Excel and Word table creation functions loaded successfully!\n")
cat("Use create_results_excel() for Excel output\n")
cat("Use create_word_ready_table() for Word-compatible HTML tables\n")