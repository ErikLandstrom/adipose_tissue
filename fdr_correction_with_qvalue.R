### fdr_correction_with qvalue
### Author: Erik Ländström
### Date: 190110


# Description -------------------------------------------------------------

# Uses the package qvalue to do fdr correction from a pvalue vector from a 
# statistical test.

# Save data from function with <-.

# Arguments ---------------------------------------------------------------

# tb = tidy tibble containing pvalues
# p_value_col = column containing p-values

# Libraries ---------------------------------------------------------------

# qvalue
# tidyverse


# Function ----------------------------------------------------------------

fdr_correction_with_qvalue <- function(tb, p_value_col) {
  
  # Quote
  p_value_col <- enquo(p_value_col)
  
  # fdr correction
  qvalues <- qvalue(as_vector(tb %>%
                                ungroup() %>% # To make sure to only get one column
                                dplyr::select(!!p_value_col)))
  
  return(qvalues)
}


