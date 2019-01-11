### multiple_2way_anova_in_tidyverse_with_tukeyHSD
### Author: Erik Ländström
### Date: 181217


# Description -------------------------------------------------------------

# Performs a two way ANOVA for each feature of the dataset (e.g. protein),
# using a tidy tibble in a longformat as input. Outputs a tibble with the 
# the ANOVA statistics using the tidy function from the broom package.
# Added tukeyHSD_test.

# Save output like "name" <- function() {}

# Use unnest(tidy_.+) to extract p-values


# Arguments ---------------------------------------------------------------

# tb = tidy tibble
# cases = column containing the name of the cases
# observations = string, column containing dependent variable
# independent_var1 = string, column containing independent factor 1
# independent_var2 = string, column containing independent factor 2

# Function ----------------------------------------------------------------

multiple_2way_anova_in_tidyverse_with_tukeyHSD <- function(tb, cases, observations,
                                                           independent_var1,
                                                           independent_var2) {
  
  # Libraries
  library(tidyverse)
  library(broom)
  
  # Quoting
  cases <- enquo(cases)
  
  # Create formula object
  lm_formula <- as.formula(paste(observations, "~", independent_var1, "*", independent_var2))
  
  # Create string with the name of the interaction effect
  interaction_effect <- paste(independent_var1, ":", independent_var2, sep = "")
  
  # 2-way ANOVA
  anova <- tb %>%
    group_by(!!cases) %>%
    nest() %>%
    mutate(anova = map(data, ~ aov(lm_formula, data = .x)),
           tukey = map(data, ~ TukeyHSD(aov(lm_formula, data = .x),
                                        which = interaction_effect)),
           tidy_anova = map(anova, tidy),
           tidy_tukey = map(tukey, tidy))
  
  return(anova)
}
