### plot_pvalue_histogram_ANOVAe
### Author: Erik Ländström
### Date: 190110


# Description -------------------------------------------------------------

# Plots a pvalue histogram for each factor in a 2-way ANOVA.


# Arguments ---------------------------------------------------------------

# tb = tibble with tidy ANOVA results
# limit_max = integer, max y axis limit

# Function ----------------------------------------------------------------

plot_pvalue_histogram_ANOVA <- function(tb, limit_max) {
  tb %>% 
    dplyr::filter(term != "Residuals") %>%
    ggplot(aes(x = p.value)) +
    geom_histogram(breaks = seq(0, 1, 0.05), fill = "black", color = "white") +
    xlab("p-value") +
    ylab("Number of proteins") +
    facet_wrap(~ term) +
    scale_y_continuous(limits = c(0, limit_max), expand = c(0, 0)) + 
    labs(
      title = "p-value histograms for each ANOVA factor"
    ) +
    theme(panel.background = element_blank(),
          axis.line = element_line(color = "black"))
}










