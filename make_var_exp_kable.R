### make_var_exp_kable
### Author: Erik Ländström
### Date: 190124


# Description -------------------------------------------------------------

# Prints a variance explained table using kable.


# Arguments ---------------------------------------------------------------

# tb = pca results
# print_format = character, specifying which format to use. Uses the default
# formats from the kable function


# Function ----------------------------------------------------------------

make_var_exp_kable <- function(tb, print_format) {
  
  library(knitr)
  library(kableExtra)
  
  kable(
    tb, 
    format    = print_format,
    digits    = 4,
    col.names = c("PC", "Variance",
                  "Proportion explained", 
                  "Cumulative proportion explained"),
    align     = "c",
    caption   = "Proportion of the variance that each principal component explains"
  ) %>% 
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = F)
}
