### create_long_observed_matrix_tibbles
### Author: Erik Ländström
### Date: 181213


# Description -------------------------------------------------------------

# Takes an output file (wide tibble) from DEP (Maxquant) and converts it to long
# format. Additionally calculates the number of observed proteins per sample
# group and creates a new tibble, and then replaces observed values (not NA) 
# with the gene name and saves output in a new tibble.


# Arguments ---------------------------------------------------------------

# tb = wide tibble, output from DEP

create_long_observed_matrix_tibbles <- function(tb) {
  # Library
  library(tidyverse)
  
  # Generate a long tibble
  original_data_long <<- tb %>%
    gather(m_wt_736:sc_wt_745, key = "Sample", value = "Intensity") %>%
    separate(col = "Sample", into = c("tissue", "genotype", "sample"), sep = "_")
  
  # Calculate number of times there is a LFQ intensity per group
  original_data_observed_matrix <<- original_data_long %>%
    group_by(tissue, genotype, name) %>%
    summarise(n = sum(!is.na(Intensity))) %>%
    replace_with_na(replace = list(n = 0)) %>%
    unite(group, tissue, genotype) %>%
    spread(key = group, value = n)
  
  # Replaces observed values with gene name
  original_data_name_matrix <<- original_data_observed_matrix %>%
    group_by(name) %>%
    mutate(m_midy = ifelse(m_midy > 0, name),
           m_wt = ifelse(m_wt > 0, name),
           sc_midy = ifelse(sc_midy > 0, name),
           sc_wt = ifelse(sc_wt > 0, name))
}
