### separate_condition_from_DEP_output
### Author: Erik Ländström
### Date: 190110


# Description -------------------------------------------------------------

# Separates the condition column from a DEP output into separate factors.


# Arguments ---------------------------------------------------------------

# tb = tibble, DEP output long fromat
# conditions = character vector, column names for the new factors

# Function ----------------------------------------------------------------

separate_condition_from_DEP_output <- function (tb, conditions) {
  temp <- tb %>% 
    separate(condition, into = conditions, sep = "_") %>% 
    rename(LFQ = intensity)
  return(temp)
}