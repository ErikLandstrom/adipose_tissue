### extract_five_protein_ids_for_biodb
### Author: Erik Ländström
### Date: 190205


# Description -------------------------------------------------------------

# Extracts five (if possible) protein IDs (ref_seq) for every protein group from
# an output file from MaxQuant (proteinGroups.txt). The proteinGroups.txt file 
# should first be processed with make_tibble_for_DEP.R and reverse hits can be 
# filtered out. 

# Save object with object_name <- extract_five_protein_ids_for_biodb ...

# The output can be saved as a txt file (or xlsx) using write_tsv from
# the readr package (or write_xlsx from readxl).

# The protein IDs can then be copied into the bioDB database to retrieve gene
# symbols.


# Arguments ---------------------------------------------------------------

# tb = tibble, output from maked_tibble_for_DEP.R


# Function ----------------------------------------------------------------

extract_five_protein_ids_for_biodb <- function(tb) {
  
  # Libraries
  library(dplyr)
  library(tidyr)
  
  
  protein_ids <- tb %>% 
    dplyr::select(`Majority protein IDs`) %>% 
    mutate(Protein_IDs = `Majority protein IDs`,
           protein_ID = str_split(Protein_IDs, ";")) %>% # Split protein group 
    dplyr::select(`Majority protein IDs`, protein_ID) %>% 
    unnest(protein_ID) %>% 
    group_by(`Majority protein IDs`) %>% 
    top_n(5) # Extracts IDs for first five ocurrences per protein group
  
  # Return tibble
  return(protein_ids)
}
