### replace_LOC_and_NA_with_gene_name
### Author: Erik Ländström
### Date: 190205


# Description -------------------------------------------------------------

# Replaces gene_1 with gene_2 if gene_1 is a locus and gene_2 has a gene symbol
# and creates a new column, gene_name. Or with third gene name, or 4th, or 5th.

# Save object with gene_names <- function ...


# Arguments ---------------------------------------------------------------

# tb = tibble, data processed with spread_gene_names.R


# Function ----------------------------------------------------------------


replace_LOC_and_NA_with_gene_name <- function(tb) {
  # Libraries
  library(dplyr)
  
  gene_names <- tb %>% 
    mutate(gene_name = if_else(str_detect(gene_1, "LOC.+") & !str_detect(gene_2,"LOC.+"), gene_2, gene_1), # Replace LOC with 2nd name
           gene_name = if_else(str_detect(gene_name, "LOC.+") & !str_detect(gene_3,"LOC.+"), gene_3, gene_name), # 3rd
           gene_name = if_else(str_detect(gene_name, "LOC.+") & !str_detect(gene_4,"LOC.+"), gene_4, gene_name), # 4th
           gene_name = if_else(str_detect(gene_name, "LOC.+") & !str_detect(gene_3,"LOC.+"), gene_5, gene_name), # 5th
           # Replace NAs with gene name or locus (for contaminants)
           gene_name = if_else(is.na(gene_name) & !is.na(gene_2), gene_2, gene_name), # Replace NA with 2nd name
           gene_name = if_else(is.na(gene_name) & !is.na(gene_3), gene_3, gene_name),
           gene_name = if_else(is.na(gene_name) & !is.na(gene_4), gene_4, gene_name),
           gene_name = if_else(is.na(gene_name) & !is.na(gene_5), gene_5, gene_name),
           gene_name = replace_na(gene_name, gene_1)) %>% # Replaces artifact NAs with 1st name
    dplyr::select(`Majority protein IDs`, gene_name, everything())
  
  return(gene_names)
}
