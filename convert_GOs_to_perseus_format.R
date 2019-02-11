### convert_GOs_to_perseus_format
### Author: Erik Ländström
### Date: 190207


# Description -------------------------------------------------------------

# Converts the clean and tidy GO categories data to clean perseus format for
# further joining with maxquant output.

# Save output with object <- function...


# Arguments ---------------------------------------------------------------

# tb = tibble, output fromtidy_gene_ontology_output_from_biodb


# Function ----------------------------------------------------------------

tidy_gene_ontology_output_from_biodb


convert_GOs_to_perseus_format <- function(tb) {
  
  # Libraries
  library(dplyr)
  library(tidyr)
  
  GOs_clean_untidy <- tb %>% 
    gather(annotation_id, annotation, key = "var", value = "val") %>%
    unite(type_var, type, var, sep = "_") %>% # Unites GO category and IDs
    group_by(`Majority protein IDs`, type_var) %>% 
    mutate(val = paste(val, collapse = ";")) %>% # Paste all GO categories per protein group and category type
    distinct() %>% # Removes all duplicates
    group_by(`Majority protein IDs`) %>% 
    spread(key = "type_var", value = "val")
  
  return(GOs_clean_untidy)
}



