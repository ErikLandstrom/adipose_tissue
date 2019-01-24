### g_test
### Author: Erik Ländström
### Date: 181214


# Description -------------------------------------------------------------

# Performs a G-test to analyse if there are any significant differences in the
# distribution of missing and observed values between the different sample groups.
# Script is adapted from Webb-Robertson, 2010, "Combined Statistical Analyses of Peptide Intensities and Peptide
# Occurrences Improves Identification of Significant Peptides from
# MS-Based Proteomics Data"


# Arguments ---------------------------------------------------------------

# tb = tibble, preprocessed with the DEP package


# Function ----------------------------------------------------------------

g_test <- function (tb) {
  # Library
  library(tidyverse)
  library(naniar)
  
  # Make long format tibble
  original_data_long <- tb %>%
    gather(m_wt_736:sc_wt_745, key = "Sample", value = "LFQ") %>%
    separate(col = "Sample", into = c("tissue", "genotype", "sample"), sep = "_") %>%
    unite(group, tissue, genotype) %>%
    dplyr::select(name, group, sample, LFQ, everything())
  
  # Make nabular tibble
  original_data_long_nab <- nabular(original_data_long)
  
  # Count number of samples
  n <- n_distinct(original_data_long_nab$sample) * 2 # Should change this in the future
  
  # Count number of groups
  n_groups <- as.integer(original_data_long_nab %>%
    summarise(n_groups = n_distinct(group)))
  
  # Count number of samples per group
  samples_per_group <- original_data_long_nab %>%
                                   group_by(group) %>%
                                   summarise(n_distinct(sample))
  
  # Calculates expected values for each peptide and sample group, and
  # each term for the G-test
  count_matrix <- original_data_long_nab %>%
    group_by(name, group) %>%
    summarise(missing = sum(is.na(LFQ)), # Number of missing values per group and protein
              observed = sum(!is.na(LFQ))) %>% # Number of observed values per group and protein
    ungroup() %>%
    group_by(name) %>%
    mutate(total_missing = sum(missing), # Calculates the number of missing values per protein (total, regardless of sample group)
           total_observed = sum(observed)) %>% # Calculates the number of missing values per protein (total, regardless of sampel group)
    group_by(group) %>%
    mutate(samples_per_group = missing + observed, # Total number of samples per group
           expected_missing = (total_missing * samples_per_group) / n, # Calculates the expected number of missing values per protein and group by random chance
           expected_observed = (total_observed * samples_per_group) / n, # Calculates the expected number of observed values per protein and group by random chance
           missing_term = missing * log(missing / expected_missing), # Calculates the G-test missing term for each protein
           observed_term = observed * log(observed / expected_observed)) # Calculates the G-test observed term for each protein
    
  # Perform G-test
  g_test <- count_matrix %>%
    group_by(name, group) %>%
    summarise(gtest_sample_term = sum(missing_term, observed_term, na.rm = TRUE),
              missing = missing,
              observed = observed) %>% # Calculates the "g-value" for each sample group
    ungroup() %>%
    group_by(name) %>%
    summarise(g_value = 2 * sum(gtest_sample_term)) %>% # Calculates the g-value for each protein
    group_by(name) %>%
    mutate(gtest_pvalue = 1 - pchisq(g_value, n_groups - 1))
  
  # Join g_test with number of observed values per group
  g_test <<- full_join(g_test, 
            count_matrix %>%
              dplyr::select(name, group, observed) %>%
              spread(key = "group", value = observed),
            by = "name")
}

