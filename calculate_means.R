### calculate_means
### Author: Erik Ländström
### Date: 190306


# Description -------------------------------------------------------------

# Calculate means from output separate_condition_from_DEP_output.R function

# Arguments ---------------------------------------------------------------

# tb = output from separate_condition_from_DEP_output

# Function ----------------------------------------------------------------

calculate_means <- function (tb) {
  
  # Total mean
  total_mean <- tb %>% 
    group_by(name) %>% 
    summarise(mean_all = mean(LFQ))
  
  # Means for each group
  group_means <- tb %>% 
    group_by(name, tissue, genotype) %>% 
    summarise(mean = mean(LFQ)) %>% 
    unite(group, tissue, genotype, sep = "_") %>%
    spread(key = "group", value = "mean")
  
  # Means for each tissue group
  tissue_means <- tb %>% 
    group_by(name, tissue) %>% 
    summarise(
      mean = mean(LFQ)
    ) %>% 
    spread(key = "tissue", value = "mean") %>% 
    mutate(tissue_diff = m - sc)
  
  # Means for each genotype group
  genotype_means <- tb %>% 
    group_by(name, genotype) %>% 
    summarise(
      mean = mean(LFQ)
    ) %>% 
    spread(key = "genotype", value = "mean") %>% 
    mutate(geno_diff = midy - wt)
  
  # join mean tables
  adipose_means <- full_join(total_mean[, 1:2], group_means, by = "name")
  adipose_means <- full_join(adipose_means, tissue_means, by = "name")
  adipose_means <- full_join(adipose_means, genotype_means, by = "name")
  
  # Return
  return(adipose_means)
}