### tidy_gene_ontology_output_from_biodb
### Author: Erik Ländström
### Date: 190206


# Description -------------------------------------------------------------

# Tidies and cleans up annotation data from biodb. Input data should be in the
# format:

# gene_name, with each gene ontology type in a separated column.

# The scripts works with the three gene ontology categories and KEGG pathways
# (KEGG pathway info in biodb, maybe othe KEGG options work too). But the script
# can easily be modified to handle more or fewer categories (although there could
# be a need to update regex expressions then)

# Save output with: object_name <- function()

# Messages means that duplicated categories are removed, I don't know why they
# are created though...


# Arguments ---------------------------------------------------------------

# tb = tibble, gene names with annotation data from biodb
# gene_col, unqouted name of gene name column

# Function ----------------------------------------------------------------

tidy_gene_ontology_output_from_biodb <- function(tb, gene_col = gene_name) {
  
  # Libraries
  library(dplyr)
  library(stringr)
  library(tidyr)
  
  # Quote
  gene_col = enquo(gene_col)
  
  # Tidy data
  ontologies <- tb %>% 
    rename(gene_name = !!gene_col, 
           bp = `GO - Biological Process`,
           cc =`GO - Cellular Component`,
           mf = `GO - Molecular Function`,
           kegg = `KEGG Pathway Info`) %>%
    dplyr::select(gene_name, bp, cc, mf, kegg) %>% 
    gather(bp, cc, mf, kegg, key = "type", value = "annotation") %>% # gathers all annotations into one column
    group_by(gene_name, type) %>% 
    nest() %>% 
    mutate(data = str_split(data, ";")) %>% # Separates all annotations per gene
    unnest() %>% 
    separate(data, into = c("annotation_id", "annotation"), sep = " \\[") %>% 
    ungroup() %>%  
    group_by(gene_name, type) %>% 
    distinct(annotation_id, .keep_all = TRUE) %>% # Removes duplicated annotations per gene that for some reason gets created
    dplyr::filter(!is.na(gene_name)) # Removes NA artifacts in gene_name
  
  # Clean up strings in annotaion_id column (identifiers)
  ontologies$annotation_id <- str_replace_all(ontologies$annotation_id, 'list\\(annotation = \\"', "")
  ontologies$annotation_id <- str_replace_all(ontologies$annotation_id, 'list\\(annotation = c\\(\\"', "")
  ontologies$annotation_id <- str_replace_all(ontologies$annotation_id, '.*\\"\\)', "")
  
  # Clean up strings in annotation column (descriptions)
  ontologies$annotation <- str_replace_all(ontologies$annotation, "Name: ", "")
  ontologies$annotation <- str_replace_all(ontologies$annotation, "\\].*", "")
  ontologies$annotation <- str_replace_all(ontologies$annotation, "Pathway Title:", "")
  
  # Replace empty strings and ")" with NA in annotation_id column
  ontologies <- ontologies %>% 
    mutate(annotation_id = ifelse(annotation_id == "", NA, annotation_id),
           annotation_id = ifelse(annotation_id == ")", NA, annotation_id)) 
  
  # Return object
  return(ontologies)
}