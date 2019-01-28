### join_qvalues_with_results
### Author: Erik Ländström
### Date: 190110


# Description -------------------------------------------------------------

# Joins qvalues to ANOVA output data.

# Save data from function with <-.

# Arguments ---------------------------------------------------------------

# tb = output from multiple_ttests.R
# tb2 = output from fdr_correction_with_qvalue.R


# Function ----------------------------------------------------------------

join_qvalues_with_results <- function(tb, tb2, p_value_col = "p_value") {
  
  # Make a tibble with p and qvalues
  qvalues <- tibble(p_value = tb2[["pvalues"]], qvalue = tb2[["qvalues"]])
  
  # Join with ANOVA results
  data <- full_join(tb, qvalues, by = p_value_col)
  
  return(data)
}