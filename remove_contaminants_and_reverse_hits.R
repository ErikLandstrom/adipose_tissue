### remove_contaminants_and_reverse_hits
### Date: 181210
### Author: Erik Ländström


# Description -------------------------------------------------------------

# Removes contaminants and reverse hits from the maxquant output file
# proteinGroups.txt.


# Arguments ---------------------------------------------------------------

# tb = output tibble from maxquant


# Function ----------------------------------------------------------------


remove_contaminants_and_reverse_hits <- function(tb) {
  # Library
  library(tidyverse)
  
  data_raw <<- tb %>%
    filter(is.na(`Potential contaminant`)) %>%
    filter(is.na(Reverse))
}