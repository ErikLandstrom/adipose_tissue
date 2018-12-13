### plot_venn_diagram
### Author: Erik Ländström
### Date: 181212


# Description -------------------------------------------------------------


# Plots venn diagram of quantified proteins, needs output from 
# create_long_observed_matrix_tibbles.R, The original_data_named_matrix tibble.


# Arguments ---------------------------------------------------------------

# tb = named observed matrix tibble
# filename = output filename (.tiff)
# title_text = title


# Function ----------------------------------------------------------------

plot_venn_diagram <- function(tb, file, title_text) {
  # Libraries
  library(VennDiagram)
  library(tidyverse)
  
  venn.diagram(list(
    m_midy  = as_vector(tb$m_midy),
    m_wt    = as_vector(tb$m_wt),
    sc_midy = as_vector(tb$sc_midy),
    sc_wt   = as_vector(tb$sc_wt)
  ),
  filename = file, na = "remove",
  fill = c("red", "green", "blue", "yellow"),
  width = 4100,
  main = title_text,
  category.names = c("Mesenterial MIDY", "Mesenterial WT", "Subcutaneous MIDY", "Subcutaneous WT"))
}