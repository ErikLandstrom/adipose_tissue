### convert_NCBIIDs_to_gene_names
### Date: 181210
### Author: Erik Ländström




convert_NCBIIDs_to_gene_names <- function(tb) {
  
  # Load sus scrofa dataset from biomaRt and change column name to Protein_ID
  ensembl <- useMart("ensembl", dataset = "sscrofa_gene_ensembl")  
  
  # Get gene names for peptide IDs
  NP <- as_tibble(getBM(attributes = c("refseq_peptide", "external_gene_name"),
                        filters = "refseq_peptide",
                        values = tb,
                        mart = ensembl
  ))
  
  NP <- NP %>%
    rename(Protein_ID = refseq_peptide)
  
  # Load sus scrofa dataset from biomaRt and change column name to Protein_ID
  XP <- as.tibble(getBM(attributes = c("refseq_peptide_predicted", "external_gene_name"),
                        filters = "refseq_peptide_predicted",
                        values = tb,
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