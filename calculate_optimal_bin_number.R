### calculate_optimal_bin_number
### Author: Erik Ländström
### Date: 181212

# Calculates the optimal number of bins for one sample accroding to the 
# Freedman-Diaconis rule

# Arguments:
# tb = a vector of values, e.g.: LFQ intensities for one sample

calculate_optimal_bin_number <- function(tb) {
  b <<- diff(range(original_data$m_wt_736, na.rm = TRUE)) / 
    (2 * IQR(original_data$m_wt_736, na.rm = TRUE) / 
       length(original_data$m_wt_736)^(1/3))
}  
  
  diff(range(original_data$m_wt_736, na.rm = TRUE)) / (2 * IQR(original_data$m_wt_736, na.rm = TRUE) / length(original_data$m_wt_736)^(1/3))
