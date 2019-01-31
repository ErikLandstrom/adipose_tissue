### make_top20_proteins_kable
### Author: Erik Ländström
### Date: 190130


# Description -------------------------------------------------------------

# Creates html table of the top 20 most significant/abundant proteins for a
# group after t-test or ANOVA.

# Can be saved in html with save_kable("filename.html") and to image with
# webshot::webshot("filename.html", file = "filename.png")

# Arguments ---------------------------------------------------------------

# tb = statistical test results, filtered by q-value and group
# top_row = numeric, number of rows to include. Negative for bottom.
# variable = which variable to use
# caption_title = title
# background_color = character, background color for header row


# Function ----------------------------------------------------------------

make_top20_proteins_kable <- function (tb, top_rows, variable, caption_title, background_color = NULL) {
  
  variable = enquo(variable)
  
  kable(tb %>% 
    ungroup() %>% 
    top_n(top_rows, !!variable) %>% 
    dplyr::select(name, log2_difference, p_value, qvalue) %>%
    left_join(protein_names, by = "name") %>%
    arrange(!!variable) %>%
    mutate_if(is.numeric, funs(as.character(signif(., 3)))) %>% 
    rename(`Gene name` = name, `Log2 fold change` = log2_difference, `p-value` = p_value, `q-value` = qvalue, `Protein name` = protein_name),
  caption = caption_title,
  align = "c") %>% 
  kable_styling(bootstrap_options = c("striped", "hover"), full_width = F) %>%
  row_spec(0, background = background_color) %>% 
  row_spec(1:20, color = "black")
}
