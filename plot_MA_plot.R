### plot_MA_plot
### Author: Erik Ländström
### Date: 181217


# Description -------------------------------------------------------------

# Plots an MA-plot from data generated with multiple_ttests.R.


# Arguments ---------------------------------------------------------------

# tb = output from multiple_ttests.R
# plot_title = plot title


# Function ----------------------------------------------------------------

plot_MA_plot <- function(tb, plot_title = "MA plot") {
  tb %>% 
    ggplot(aes(mean_total, log2_difference)) +
    geom_point() +
    geom_hline(yintercept = 0, color = "red") +
    xlab("Mean intensity") +
    ylab("Log2 fold difference") +
    labs(
      title = plot_title
    ) +
    theme(
      panel.background = element_blank(),
      axis.line = element_line("black")) 
}