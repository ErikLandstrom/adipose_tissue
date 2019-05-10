
impute_NAs_with_10_percent_of_mean <- function(tb) {
  temp <- read_tsv(tb$datapath)
  
  temp <- temp %>% 
    group_by(name) %>% 
    mutate(missing = if_else(is.na(intensity), "Missing", "aObserved"),
           intensity = replace_na(intensity, log2(mean(2 ^ intensity, na.rm = TRUE) * 0.1)))
  
  return(temp)
}