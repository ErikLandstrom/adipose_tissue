### plot_pca
### Author: Erik Ländström
### Date: 190125


# Description -------------------------------------------------------------

# Plots the first two principal components of a pca.


# Arguments ---------------------------------------------------------------

# tb = nested pca object
# plot_title = character, plot title


# Function ----------------------------------------------------------------

plot_pca <- function(tb, plot_title) {
  
  library(tidyverse)
  library(purrr)
  
  # Plot the first two principal components
  p <- tb %>% 
    mutate(
      pca_graph = map2(
        .x = pca,
        .y = data,
        ~ autoplot(
          .x,
          data         = .y,
          colour       = "group",
          size         = 3,
          loadings     = FALSE,
          label        = TRUE,
          label.label  = "sample",
          label.repel  = TRUE,
          label.colour = "black"
        ) +
          theme(panel.background = element_blank(),
                axis.line        = element_line(color = "black"),
                plot.title       = element_text(hjust = 0.5),
                legend.key       = element_blank()) +
          labs(x     = "Principal Component 1",
               y     = "Principal Component 2",
               title = plot_title)
      )
    ) %>%
    pull(pca_graph)
  
  # Add horizontal and vertical lines to the plot
  p1 <- p[[1]] + # p[[1]] is used to get the "pure" plot
    geom_hline(yintercept = 0,
               color = "grey",
               linetype = "dashed") +
    geom_vline(xintercept = 0,
               color = "grey",
               linetype = "dashed")
  
  # Reverse the layers to plot lines beneath the scatterplot
  p1$layers <- rev(p1$layers)
  
  # Plot
  p1
}