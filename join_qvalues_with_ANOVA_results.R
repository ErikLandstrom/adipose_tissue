### join_qvalues_with_ANOVA_results
### Author: Erik Ländström
### Date: 190110


# Description -------------------------------------------------------------

# Joins qvalues to ANOVA output data.

# Save data from function with <-.

# Arguments ---------------------------------------------------------------

# tb = output from fdr_correction_with_qvalue.R


# Function ----------------------------------------------------------------

join_qvalues_with_ANOVA_results <- function(tb) {
  
  # Make a tibble with p and qvalues
  qvalues <- tibble(p.value = qvalues[["pvalues"]], qvalue = qvalues[["qvalues"]])
  
  # Join with ANOVA results
  anova_results <- full_join(tb, qvalues, by = "p.value")
  
  return(anova_results)
}
