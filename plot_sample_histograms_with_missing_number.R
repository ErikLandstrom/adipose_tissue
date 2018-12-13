### plot_histogram_with_missing_values
### Author: Erik Ländström
### Date: 181213


# Description -------------------------------------------------------------

# plot_sample_histograms_with_missing_values takes a nabular tibble in long 
# format and plots a histogram for each sample with missing values imputed at 
# 10% of minimum LFQ intensity with optimal number of bins according to the 
# Freedman-Diakonis rule.


# Arguments ---------------------------------------------------------------

# tb = nabular tibble in long format. Created with nabular()
# sample_number = sample ID number
# tissue_type = string with tissue type 


# Libraries ---------------------------------------------------------------

# tidyverse
# naniar


# Function ----------------------------------------------------------------

plot_sample_histograms_with_missing_values <- function(tb, sample_number, tissue_type) {
  # Libraries
  library(tidyverse)
  library(naniar)
  
   vec <- as_vector(tb %>%
     filter(sample == sample_number,
            tissue == tissue_type) %>%
     select(LFQ))
  
  # Calculates the optimal number of bins for one sample according to the Freedman-Diaconis rule
  b <- diff(range(vec, na.rm = TRUE)) / (2 * IQR(vec, na.rm = TRUE) / length(vec)^(1 / 3))
  
  # Plot histogram
  tb %>%
    mutate(LFQ = impute_below(LFQ)) %>%
    ggplot(aes(LFQ, fill = LFQ_NA)) +
    geom_histogram(bins = b) +
    facet_wrap(tissue ~ sample)
}
