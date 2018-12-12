### make_tibble_for_DEP
### Date: 181210
### Author: Erik Ländström

# Reads and converts the proteinGroups.txt file from maxquant into file that is
# ready to use with the DEP package.

make_tibble_for_DEP <- function(tb = "ProteinGroups.txt") {
  data_raw <- read_tsv(tb)
  
  # Extract LFQ column names
  lfq_col_names <- str_extract(colnames(data_raw), "LFQ.+")
  lfq_col_names <- lfq_col_names[which(!is.na(lfq_col_names))]
  
  data_raw <<- data_raw %>%
    dplyr::select(`Protein IDs`, `Majority protein IDs`,
                  `Fasta headers`, Peptides, 
                  `Razor + unique peptides`, `Unique peptides`, 
                  lfq_col_names,
                  `Only identified by site`, `Reverse`,
                  `Potential contaminant`)
}