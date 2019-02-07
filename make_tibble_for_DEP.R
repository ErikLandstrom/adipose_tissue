### make_tibble_for_DEP
### Date: 181210
### Author: Erik Ländström


# Description -------------------------------------------------------------

# Reads and converts the proteinGroups.txt file from maxquant into file that is
# ready to use with the DEP package.


# Arguments ---------------------------------------------------------------

# file_name = a text file, output from maxquant called "proteinGroups.txt"


# Function ----------------------------------------------------------------

make_tibble_for_DEP <- function(file_name = "proteinGroups.txt") {
  # Library 
  library(tidyverse)
  
  # Read file
  data_raw <- read_tsv(file_name, col_types = list("Reverse" = col_character())
  
  # Extract LFQ column names
  lfq_col_names <- str_extract(colnames(data_raw), "LFQ.+")
  lfq_col_names <- lfq_col_names[which(!is.na(lfq_col_names))]
  
  data_raw <- data_raw %>%
    dplyr::select(`Protein IDs`, `Majority protein IDs`,
                  `Fasta headers`, Peptides, 
                  `Razor + unique peptides`, `Unique peptides`, 
                  lfq_col_names,
                  `Only identified by site`, `Reverse`,
                  `Potential contaminant`)
  
  # Return object
  return(data_raw)
}