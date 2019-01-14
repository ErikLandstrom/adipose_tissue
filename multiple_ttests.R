### multiple_ttests
### Author: Erik Ländström
### Date: 181217


# Description -------------------------------------------------------------

# Performs a t-test for each feature of the dataset (e.g. protein),
# using a tidy tibble in a longformat as input. Outputs a tibble with the 
# the statistics for ttest object using the tidy function from the broom package.

# Save output with t_test_+ <- multiple_ttests()

# Arguments ---------------------------------------------------------------

# tb = tidy tibble in long format
# col_name = column name of the variable of interest, unquoted

# Function ----------------------------------------------------------------

multiple_ttests <- function(tb, col_name) {
  
  # Libraries
  library(tidyverse)
  library(broom)
  
  # Quote
  col_name <- enquo(col_name)
  
  # ttest
  t_test <-  tb %>%
    dplyr::select(name, !!col_name, LFQ) %>%
    group_by(name, !!col_name) %>%
    nest() %>%
    spread(key = !!col_name, value = data) %>%
    group_by(name) %>%
    mutate(p_value = t.test(unlist(midy), unlist(wt))$p.value,
           mean_midy = mean(unlist(midy)),
           mean_wt = mean(unlist(wt)),
           log2_difference = mean_midy - mean_wt,
           average_mean = (mean_midy + mean_wt) / 2)
  
  return(t_test)
}

