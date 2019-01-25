### plot_multiple_PCAs
### Author: Erik Ländström
### Date: 190125


# Description -------------------------------------------------------------

# Plots the specified variations of principal components.


# Arguments ---------------------------------------------------------------

# tb       = nested pca object
# pc_start = double, indicates the first principal component to be plotted.
#            Used for the x axis.
# pc_end   = double, indicates the final pc to be plotted.


# Function ----------------------------------------------------------------


plot_multiple_PCAs <- function(tb, pc_start, pc_end) {
  # Create empty list
  pca_plot_list <- list()
  
  # Save plots in the list
  for (i in (pc_start + 1):pc_end) {
    
    # Plot
    temp_plot <-  autoplot(tb[["pca"]][[1]], data = tb[["data"]][[1]],
                                x = pc_start,
                                y = i,
                                colour       = "group",
                                size         = 3,
                                loadings     = FALSE,
                                label        = TRUE,
                                label.label  = "sample",
                                label.repel  = TRUE,
                                label.colour = "black") +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
      theme(plot.title = element_text(hjust = 0.5),
            panel.background = element_blank(),
            panel.grid = element_blank(),
            axis.line = element_line("black")) 
    
    # Reverse the plot layers
    temp_plot$layers <- rev(temp_plot$layers)
    
    # Save plot in the list
    pca_plot_list[[i - pc_start]] <- temp_plot
  }
  
  # Return the final list
  return(pca_plot_list)
}
