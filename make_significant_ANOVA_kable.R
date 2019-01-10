### make_significant_ANOVA_kable
### Author: Erik Ländström
### Date: 190110


# Description -------------------------------------------------------------

# Calculates and prints a table of significant differences
# detected after performing a statistical test.


# Arguments ---------------------------------------------------------------

# tb = tibble with tidy ANOVA results


# Function ----------------------------------------------------------------


make_significant_ANOVA_kable <- function(tb) {
    
  kable(tb %>%
          dplyr::filter(term != "Residuals") %>%
          group_by(term) %>%
          summarise(n = sum(p.value < 0.05)),
        caption = "Number of significant differences per factor.", format = "pandoc",
        align = "c")
}
