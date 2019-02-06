### spread_gene_names
### Author: Erik Ländström
### Date: 190205


# Description -------------------------------------------------------------

# Converts the format of the output file from biodb, to create a tibble with
# protein group in column 1, and with separate columns for gene symbols (5 
# columns). Fills observations with NA if a protein group has fewer than
# 5 protein IDs.

# The gene names from biodb should be copied to the output of extract_five_
# protein_ids_for_biodb.R and then read with read_tsv or read_xlsx depending
# on format.

# Save object with object_name <- function()


# Arguments ---------------------------------------------------------------

# tb = tibble, data from biodb
# value_col = character, column name for gene name column. Default = Gene Symbol


# Function ----------------------------------------------------------------

spread_gene_names <- function(tb, value_col = "Gene Symbol") {
  # Libraries
  library(dplyr)
  library(tidyr)
  
  gene_names <- tb %>% 
    mutate(`Gene Symbol` = ifelse(`Gene Symbol` == "-", NA, `Gene Symbol`)) %>% # Replaces - with NAs
    group_by(`Majority protein IDs`) %>% 
    mutate(index = 1:n()) %>%  # Creates an index for use with spread
    select(-protein_ID) %>% 
    spread(key = "index", value = value_col) %>% 
    rename(gene_1 = `1`, gene_2 = `2`, gene_3 = `3`, gene_4 = `4`, gene_5 = `5`)
  
  # Return object
  return(gene_names)
}
