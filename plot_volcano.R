### plot_volcano
### Author: Erik Ländström
### Date: 181217


# Description -------------------------------------------------------------

# Plots a volcano from data generated with multiple_ttests.R.


# Arguments ---------------------------------------------------------------

# tb = output from multiple_ttests.R
# plot_title = plot title


# Function ----------------------------------------------------------------

plot_volcano <- function(tb, plot_title = "Volcano plot") {
  
  # Volcano plot
  tb %>% 
    ggplot(aes(log2_difference, -log10(p_value))) +
    geom_point() +
    xlab("Log2 fold difference") +
    ylab("-log10(p-value") +
    labs(
      title = plot_title
    ) +
    theme(
      panel.background = element_blank(),
      axis.line = element_line("black")) 
}