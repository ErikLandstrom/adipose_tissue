### prepare_anova_without_int_for_volcano_plot
### Author: Erik Ländström
### Date: 181217


# Description -------------------------------------------------------------

# Prepares ANOVA without interaction for volcano plot. Using the output from
# calculate means and multiple_2way_anova_in_tidyverse_with_no_interaction.

# Arguments ---------------------------------------------------------------

# tb       = tibble with means
# tb_anova = ANOVA tibble

# Function ----------------------------------------------------------------

prepare_anova_without_int_for_volcano_plot <- function (tb, tb_anova) {
  diffs <- tb %>% 
    rename(tissue = tissue_diff,
           genotype = geno_diff) %>% 
    dplyr::select(name, tissue, genotype) %>%
    gather(tissue_diff, geno_diff, key = "term", value = "lf2c") # term is used for easy joining with anova results
  
  anova_results_without_interaction <- full_join(tb_anova, adipose_lfc2,
                                                 by = c("name", "term")) %>% 
    rename(group = term) %>% 
    mutate(enriched_1 = as.numeric(lfc2 > 0 & qvalue < 0.05),
           enriched_2 = as.numeric(lfc2 < 0 & qvalue < 0.05),
           significant_temp = (enriched_1 + enriched_2 * 2) + 1,
           significant = case_when(group == "tissue" & significant_temp > 1 ~ significant_temp + 2,
                                   group == "tissue" & significant_temp == 1 ~ significant_temp,
                                   group == "genotype" ~ significant_temp),
           significant = as.factor(significant))
}