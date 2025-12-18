## -----------------------------------------------------------------------------
## Script name: covid_data_aggregation.R
##
## Purpose of script: To combine state-level COVID case data into a single
## aggregated file
##
## Author: George Dewey
##
## Date Created: 2024-07-09
##
## Last Updated: 2024-07-10
## -----------------------------------------------------------------------------
library(tidyverse)
library(readr)

all_covid_data =  list.files(path="data/ground_truth/COVID/", full.names = TRUE) %>% 
  lapply(read_csv) %>% 
  bind_rows 

write_csv(all_covid_data, file = 'data/ground_truth/COVID/COVID_aggregated.csv')
