### join_genes_GOs_mqdata
### Author: Erik Ländström
### Date: 190207


# Description -------------------------------------------------------------

# Joins together maxquant data, gene names and GOs.


# Arguments ---------------------------------------------------------------

# genes_data, output from replace_LOC_and_NA_with_gene_name.R
# GO_data, output from convert_GOs_to_perseus_format
# mq_data, original data imported from maxquant


# Function ----------------------------------------------------------------

join_genes_GOs_mqdata <- function(genes_data, GO_data, mq_data) {
  
  # Libraries
  library(tidyr)
  library(dplyr)
  
  # Unite all gene name columns wiht ";" as separator
  gene_names_untidy <- genes_data %>% 
    unite(alt_gene_names, gene_1:gene_5, sep = ";")
  
  # Join gene names with GOs ready in perseus format
  genes_GOs_untidy <- full_join(gene_names_untidy, GO_data, by = "Majority protein IDs")
  
  # Join gene names and GOs with maxquant data
  mq_data <- full_join(mq_data, genes_GOs_untidy, by = "Majority protein IDs")
  
  return(mq_data)
}