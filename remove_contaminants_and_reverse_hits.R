### remove_contaminants_and_reverse_hits
### Date: 181210
### Author: Erik Ländström




remove_contaminants_and_reverse_hits <- function(tb) {
  data_raw <<- tb %>%
    dplyr::filter(is.na(`Potential contaminant`)) %>%
    dplyr::filter(is.na(Reverse))
}