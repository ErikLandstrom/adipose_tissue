




prepare_data_for_volcano_plot <- function (tb, group_neg, group_pos) {
  data_for_volcano <- tb %>% 
    mutate(significant = case_when(
      log2_difference < 0 & qvalue < 0.05 ~ group_neg,
      log2_difference > 0 & qvalue < 0.05 ~ group_pos,
      qvalue >= 0.05 ~ "Not significant"
    ), significant = as_factor(significant))
  
  return(data_for_volcano)
}

