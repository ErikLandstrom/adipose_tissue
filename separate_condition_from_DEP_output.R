adipose_imputed <- adipose_imputed %>%
  separate(condition, into = c("tissue", "genotype"), sep = "_") %>%
  rename(LFQ = intensity)

separate_condition_from_DEP_output <- function (tb) {
  temp <- tb %>% 
    separate(condition, into = c("tissue", "genotype"), sep = "_") %>% 
    rename(LFQ = intensity)
  return(temp)
}