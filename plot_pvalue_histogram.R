### plot_pvalue_histogram
### Author: Erik Ländström
### Date: 190110


# Description -------------------------------------------------------------

# Plots a pvalue histogram for tidy results.


# Arguments ---------------------------------------------------------------

# tb = tibble with tidy results
# col_name = name of column name with p-values
# limit_max = integer, max y axis limit
# plot_title = plot title, default "p-value histogram"

# Function ----------------------------------------------------------------

plot_pvalue_histogram <- function(tb, col_name, limit_max,
                                  plot_title = "p-value histogram") {
  
  # Quote
  col_name <- enquo(col_name)
  
  # pvalue plot
  tb %>% 
    ggplot(aes(x = !!col_name)) +
    geom_histogram(breaks = seq(0, 1, 0.05), fill = "black", color = "white") +
    xlab("p-value") +
    ylab("Number of proteins") +
    scale_y_continuous(limits = c(0, limit_max), expand = c(0, 0)) + 
    labs(
      title = plot_title
    ) +
    theme(panel.background = element_blank(),
          axis.line = element_line(color = "black"))
}










