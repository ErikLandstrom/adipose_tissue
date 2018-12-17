### multiple_ttests
### Author: Erik Ländström
### Date: 181217


# Description -------------------------------------------------------------

# Performs a t-test for each feature of the dataset (e.g. protein),
# using a tidy tibble in a longformat as input. Outputs a tibble with the 
# the statistics for ttest object using the tidy function from the broom package.


# Arguments ---------------------------------------------------------------

# tb = tidy tibble in long format
# col_name = character string with column name of the variable of interest

# Function ----------------------------------------------------------------

multiple_ttests <- function(tb, col_name) {
  
  # Libraries
  library(tidyverse)
  library(broom)
  
  # t-test
  sym(col_name)
  
  t_test <<-  tb %>%
    select(name, !!sym(col_name), LFQ) %>%
    group_by(name, !!sym(col_name)) %>%
    nest() %>%
    spread(key = !!sym(col_name), value = data) %>%
    group_by(name) %>%
    mutate(test = t.test(unlist(midy), unlist(wt))$p.value,
           mean_midy = mean(unlist(midy)),
           mean_wt = mean(unlist(wt)),
           log2_difference = mean_midy - mean_wt)
}

