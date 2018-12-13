### convert_NCBIIDs_to_gene_names
### Date: 181210
### Author: Erik Ländström


# Description -------------------------------------------------------------

# Takes refseq protein IDs from NCBI and converts them to gene names using the
# bioMart package. Saves gene names and IDs in a tibble called gene names that 
# can be joined with the original dataset 


# Arguments ---------------------------------------------------------------

# vec = a vector of protein refseq IDs


# Function ----------------------------------------------------------------


convert_NCBIIDs_to_gene_names <- function(vec) {
  # Library
  library(biomaRt)
  library(tidyverse)
  
  # Load sus scrofa dataset from biomaRt and change column name to Protein_ID
  ensembl <- useMart("ensembl", dataset = "sscrofa_gene_ensembl")  
  
  # Get gene names for peptide IDs
  NP <- as_tibble(getBM(attributes = c("refseq_peptide", "external_gene_name"),
                        filters = "refseq_peptide",
                        values = vec,
                        mart = ensembl
  ))
  
  NP <- NP %>%
    rename(Protein_ID = refseq_peptide)
  
  # Load sus scrofa dataset from biomaRt and change column name to Protein_ID
  XP <- as.tibble(getBM(attributes = c("refseq_peptide_predicted", "external_gene_name"),
                        filters = "refseq_peptide_predicted",
                        values = vec,
                        mart = ensembl
  ))
  
  XP <- XP %>%
    rename(Protein_ID = refseq_peptide_predicted)
  
  # Use full_join to join the two tables
  gene_names <- full_join(NP, XP, by = "Protein_ID")
  
  # Tidy data, gene names in one column and remove NAs
  gene_names <- gene_names %>%
    gather(external_gene_name.x, external_gene_name.y, key = "Source", value = "Gene_name") %>%
    dplyr::select(Protein_ID, Gene_name) %>%
    dplyr::filter(!is.na(Gene_name)) 
  
  # If there are any "" artifacts from the biomaRt search, replace them with NA
  gene_names[gene_names == ""] <- NA
  
  # Save output in a tibble
  gene_names <<- gene_names
}