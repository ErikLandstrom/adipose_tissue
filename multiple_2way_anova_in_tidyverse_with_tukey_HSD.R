### multiple_2way_anova_in_tidyverse
### Author: Erik Ländström
### Date: 181217


# Description -------------------------------------------------------------

# Performs a two way ANOVA for each feature of the dataset (e.g. protein),
# using a tidy tibble in a longformat as input. Outputs a tibble with the 
# the ANOVA statistics using the tidy function from the broom package.

# Save output like "name" <- function() {}

# Use unnest(tidied) to extract p-values


# Arguments ---------------------------------------------------------------

# tb = tidy tibble


# Function ----------------------------------------------------------------

multiple_2way_anova_in_tidyverse <- function(tb) {
  
  # Libraries
  library(tidyverse)
  library(broom)
  
  # 2-way ANOVA
  anova <- tb %>%
    group_by(name) %>%
    nest() %>%
    mutate(anova = map(data, ~ aov(LFQ ~ genotype * tissue, data = .x)),
           tukey = map(data, ~ TukeyHSD(aov(LFQ ~ genotype * tissue, data = .x))),
           tidy_anova = map(anova, tidy),
           tidy_tukey = map(tukey, tidy))
  
  return(anova)
}
