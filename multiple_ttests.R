### multiple_ttests
### Author: Erik Ländström
### Date: 181217


# Description -------------------------------------------------------------

# Performs a t-test for each feature of the dataset (e.g. protein),
# using a tidy tibble in a longformat as input. Outputs a tibble with the 
# the statistics for ttest object using the tidy function from the broom package.

# Save output with t_test_+ <- multiple_ttests()

# Arguments ---------------------------------------------------------------

# tb        = tidy tibble in long format
# group_var = column name of the variable of interest, unquoted
# group_1   = name of group 1, unquoted
# group_2   = name of group 2, unquoted
# mean_1    = name of mean column for group_1
# mean_2    = name of mean column for group_2
# values    = column name of observations
# equal_var = logical, if equal variance should be assumed

# Function ----------------------------------------------------------------

multiple_ttests <- function(tb, group_var, 
                            group_1,
                            group_2,
                            mean_1,
                            mean_2,
                            values = LFQ, 
                            equal_var = FALSE) {
  
  # Libraries
  library(tidyverse)
  library(broom)
  
  # Quote
  group_var <- enquo(group_var)
  values    <- enquo(values)
  group_1   <- enquo(group_1)
  group_2   <- enquo(group_2)
  mean_1_name <- enquo(mean_1)
  mean_2_name <- enquo(mean_2)
  
  # Calculate the mean for each protein
  means <- tb %>%
    group_by(name) %>%
    summarise(mean_total = mean(!!values))
  
  
  # ttest
  t_test <-  tb %>%
    dplyr::select(name, !!group_var, !!values) %>%
    group_by(name, !!group_var) %>%
    nest() %>%
    spread(key = !!group_var, value = data) %>%
    group_by(name) %>%
    mutate(p_value = t.test(unlist(!!group_1), unlist(!!group_2), var.equal = equal_var)$p.value,
           !!mean_1_name := mean(unlist(!!group_1)),
           !!mean_2_name := mean(unlist(!!group_2)),
           log2_difference = !!mean_1_name - !!mean_2_name)
  
  # Join data sets
  t_test <- left_join(t_test, means, by = "name")
  
  return(t_test)
}

