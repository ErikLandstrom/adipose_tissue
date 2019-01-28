### linear_model_1_ind_var
### Author: Erik Ländström
### Date: 190128


# Description -------------------------------------------------------------

# Performs a linear model using one independent variable and extracts 
# residuals from the statistical object for further analysis.


# Arguments ---------------------------------------------------------------

# tb = tidy tibble
# cases = column containing the name of the cases
# observations = string, column containing dependent variable
# independent_var = string, column containing independent factor 

# Function ----------------------------------------------------------------

linear_model_1_ind_var <- function(tb, cases, observations, ind_variable) {
  
  
  cases <- enquo(cases)
  
  # Paste formula
  lm_formula <- as.formula(paste(observations, "~", ind_variable))
  
  # Linear model
  residuals <- tb %>% 
    group_by(!!cases) %>% 
    nest() %>% 
    mutate(model = map(data, ~ lm(lm_formula, data = .x)),
           resids = map2(data, model, add_residuals)) %>% 
    unnest(resids)
  
  # Return dataframe
  return(residuals)
}