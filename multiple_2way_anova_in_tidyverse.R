### multiple_2way_anova_in_tidyverse
### Author: Erik Ländström
### Date: 181217


# Description -------------------------------------------------------------

# Performs a two way ANOVA for each feature of the dataset (e.g. protein),
# using a tidy tibble in a longformat as input. Outputs a tibble with the 
# the ANOVA statistics using the tidy function from the broom package.


# Arguments ---------------------------------------------------------------

# tb = tidy tibble
# factor1 = first independent variable as character
# factor2 = second independent variable as character
# 

# Function ----------------------------------------------------------------

multiple_2way_anova_in_tidyverse <- function(tb, factor1, factor2) {
  
  # Libraries
  library(tidyverse)
  library(broom)
  
  # 2-way ANOVA
  anova <- tb %>%
    group_by(name) %>%
    nest() %>%
    mutate(anova = map(data, ~ aov(LFQ ~ eval(as.name(factor1)) * eval(as.name(factor2)), data = .x)),
           tidied = map(anova, tidy)) %>%
    unnest(tidied)
  return(anova)
}
