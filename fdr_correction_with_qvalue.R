### fdr_correction_with qvalue
### Author: Erik Ländström
### Date: 190110


# Description -------------------------------------------------------------

# Uses the package qvalue to do fdr correction from a pvalue vector from a 
# statistical test.

# Save data from function with <-.

# Arguments ---------------------------------------------------------------

# tb = tidy tibble containing pvalues


# Libraries ---------------------------------------------------------------

# qvalue
# tidyverse


# Function ----------------------------------------------------------------

fdr_correction_with_qvalue <- function(tb) {
  qvalues <- qvalue(as_vector(tb %>%
                                dplyr::select(p.value)))
  
  return(qvalues)
}