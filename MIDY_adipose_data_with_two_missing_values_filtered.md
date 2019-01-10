---
title: "MIDY_adipose_data_with_two_missing_values_filtered"
author: "Erik Ländström"
date: "14 December 2018"
output: 
  html_document:
    keep_md: TRUE
---



# Libraries


```r
library(tidyverse)
library(DEP)
```

```
## Warning in fun(libname, pkgname): mzR has been built against a different Rcpp version (0.12.16)
## than is installed on your system (1.0.0). This might lead to errors
## when loading mzR. If you encounter such issues, please send a report,
## including the output of sessionInfo() to the Bioc support forum at 
## https://support.bioconductor.org/. For details see also
## https://github.com/sneumann/mzR/wiki/mzR-Rcpp-compiler-linker-issue.
```

```r
library(conflicted)
library(biomaRt)
library(naniar)
library(knitr)
library(readxl)
library(gridExtra)
library(visdat)
library(qvalue)
```


# 2. Preprocessing

## 2.1 Read rawdata from Maxquant


```r
make_tibble_for_DEP <- function(tb = "ProteinGroups.txt") {
  data_raw <- read_tsv(tb)
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
```


```r
make_tibble_for_DEP("proteinGroups.txt")
```

## 2.2 Remove contaminant and reverse hits

Calculate number of contaminants and reverse hits.


```r
data_raw %>%
  summarise(n_reverse = length(Reverse) - sum(is.na(Reverse)),
            n_contaminant = length(Reverse) - sum(is.na(`Potential contaminant`)))
```

```
## # A tibble: 1 x 2
##   n_reverse n_contaminant
##       <int>         <int>
## 1        33            49
```


```r
remove_contaminants_and_reverse_hits <- function(tb) {
  data_raw <<- tb %>%
    dplyr::filter(is.na(`Potential contaminant`)) %>%
    dplyr::filter(is.na(Reverse))
}
```


```r
remove_contaminants_and_reverse_hits(data_raw)
dim(data_raw)
```

```
## [1] 2779   29
```

## 2.3 Add gene names from BiomaRt


```r
# Remove isoform from refseq identifier
data_raw <- data_raw %>%
  mutate(Protein_ID = str_extract(`Majority protein IDs`, "NP_[0-9]+|XP_[0-9]+"))
```

Extract gene names from biomart.


```r
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

convert_NCBIIDs_to_gene_names(data_raw$Protein_ID)
```

Join gene names with raw data by using protein IDs as identifiers.


```r
data <- full_join(data_raw, gene_names, by = "Protein_ID")
```

Check if there are any duplicated gene names and identifiers.


```r
# gene names
data$Gene_name %>% duplicated() %>% any()
```

```
## [1] TRUE
```

```r
# identifiers
data$Protein_ID %>% duplicated() %>% any()
```

```
## [1] FALSE
```

Duplicated gene names were found in the data.


```r
# Number of NA gene names
data %>%
  summarise(sum(is.na(Gene_name)))
```

```
## # A tibble: 1 x 1
##   `sum(is.na(Gene_name))`
##                     <int>
## 1                     442
```

```r
# Count number of instances a gene name appears
data %>%
  group_by(Gene_name) %>%
  summarise(n = n()) %>%
  dplyr::filter(n > 1) %>%
  arrange()
```

```
## # A tibble: 36 x 2
##    Gene_name     n
##    <chr>     <int>
##  1 ATP5F1C       2
##  2 C4BPA         2
##  3 CALU          2
##  4 CAP1          2
##  5 CKM           2
##  6 CMBL          2
##  7 COL14A1       2
##  8 CRK           2
##  9 EPB41L2       2
## 10 EPB41L3       2
## # ... with 26 more rows
```

```r
# Make unique "row names", either gene name or identifier
data_unique <- make_unique(data, "Gene_name", "Protein_ID")

data_unique %>%
  dplyr::select(name, ID) %>% 
  duplicated() %>%
  any()
```

```
## [1] FALSE
```

## 2.4 Generate a SummarizedExperiment


```r
# Grep for lfq columns
lfq_columns <- grep("LFQ.", colnames(data_unique))

# Read experimental design
experimental_design <- read_xlsx("MIDY_adipose_experimental_design.xlsx")

# Generate summarized experiment
data_se <- make_se(data_unique, lfq_columns, experimental_design)

data_se
```

```
## class: SummarizedExperiment 
## dim: 2779 20 
## metadata(0):
## assays(1): ''
## rownames(2779): NP_001001534 DES ... TIMM9 L3HYPDH
## rowData names(13): Protein IDs Majority protein IDs ... name ID
## colnames(20): m_wt_736 sc_wt_736 ... m_wt_745 sc_wt_745
## colData names(4): label ID condition replicate
```

# 3. Analysis

## 3.1 Missing value filtering


```r
plot_frequency(data_se)
```

![](MIDY_adipose_data_with_two_missing_values_filtered_files/figure-html/frequency-1.png)<!-- -->



```r
# Filter for proteins present in all samples of at least one condition
data_2_valid_values <- filter_missval(data_se, thr = 3)
```


```r
# Imputation
imputed_data <- impute(data_2_valid_values, fun = "man", shift = 1.8, scale = 0.3)
```


```r
# Save imputed data in a tibble
adipose_imputed <- get_df_long(imputed_data)

# Tidy up the condition column
adipose_imputed <- adipose_imputed %>%
  separate(condition, into = c("tissue", "genotype"), sep = "_") %>%
  rename(LFQ = intensity)
```

# 2-way ANOVA

Perform a 2-way ANOVA on the imputed dataset.




```r
adipose_anova <- multiple_2way_anova_in_tidyverse_with_tukeyHSD(adipose_imputed)
```

Unnest anova list-column.




```r
# Extract p-values and more
anova_results <- adipose_anova %>%
  unnest(tidy_anova) %>% 
  dplyr::filter(term != "Residuals")

# Count number of significant p-values
make_significant_ANOVA_kable(anova_results)
```



Table: Number of significant differences per factor.

      term           n  
-----------------  -----
    genotype        180 
 genotype:tissue    105 
     tissue         684 

```r
# Plot pvalue histograms for each factor
plot_pvalue_histogram_ANOVA(anova_results, 750)
```

![](MIDY_adipose_data_with_two_missing_values_filtered_files/figure-html/anova_results-1.png)<!-- -->


```r
# Function
fdr_correction_with_qvalue <- function(tb) {
  qvalues <- qvalue(as_vector(tb %>%
                                dplyr::select(p.value)))

  return(qvalues)
}

# fdr correction with the qvalue package
qvalues <- fdr_correction_with_qvalue(anova_results)

# Plot fdr statistics
plot(qvalues)
```

![](MIDY_adipose_data_with_two_missing_values_filtered_files/figure-html/unnamed-chunk-11-1.png)<!-- -->

```r
hist(qvalues)
```

![](MIDY_adipose_data_with_two_missing_values_filtered_files/figure-html/unnamed-chunk-11-2.png)<!-- -->

```r
# Function
join_qvalues_with_ANOVA_results <- function(tb) {
  qvalues <- tibble(p.value = qvalues[["pvalues"]], qvalue = qvalues[["qvalues"]])
  
  anova_results <- full_join(anova_results, qvalues, by = "p.value")
  return(anova_results)
}

# Join qvalues with anova results
anova_results <- join_qvalues_with_ANOVA_results(qvalues)

kable(anova_results %>%
        dplyr::filter(term != "Residuals") %>%
        group_by(term) %>%
        summarise(n = sum(qvalue < 0.05)),
      caption = "Number of significant qvalues.", format = "pandoc")
```



Table: Number of significant qvalues.

term                 n
----------------  ----
genotype            31
genotype:tissue      9
tissue             342

```r
kable(anova_results %>%
        dplyr::filter(term != "Residuals") %>%
        group_by(term) %>%
        summarise(n = sum(qvalue < 0.10)),
      caption = "Number of significant q-values.", format = "pandoc")
```



Table: Number of significant q-values.

term                 n
----------------  ----
genotype            72
genotype:tissue     29
tissue             463

## Tukey results


```r
tukey_results <- adipose_anova %>%
  unnest(tidy_tukey)
```

# 2-way ANOVA without interaction effects

Since there doesn't seem to be many interaction effect, I want to see what 
will happen if you do an ANOVA without the interaction effect between genotype
and tissue.




```r
# Perform ANOVA without interaction effect
adipose_anova_without_interaction <- multiple_2way_anova_in_tidyverse_with_no_interaction(adipose_imputed)

# Extract p-values and more
anova_results <- adipose_anova_without_interaction %>%
  unnest(tidy_anova) %>% 
  dplyr::filter(term != "Residuals")

# Count number of significant p-values and print a kable
make_significant_ANOVA_kable(anova_results)
```



Table: Number of significant differences per factor.

   term       n  
----------  -----
 genotype    180 
  tissue     684 

```r
# Plot pvalue histograms for each factor
plot_pvalue_histogram_ANOVA(anova_results, 750)
```

![](MIDY_adipose_data_with_two_missing_values_filtered_files/figure-html/no_interaction1-1.png)<!-- -->

The p-value histograms doesn't look any different, but with fewer p-values to
correct maybe more q-values will be significant.


```r
# fdr correction with the qvalue package
qvalues <- fdr_correction_with_qvalue(anova_results)

# Plot fdr statistics
plot(qvalues)
```

![](MIDY_adipose_data_with_two_missing_values_filtered_files/figure-html/no_interaction2-1.png)<!-- -->

```r
hist(qvalues)
```

![](MIDY_adipose_data_with_two_missing_values_filtered_files/figure-html/no_interaction2-2.png)<!-- -->

```r
# Join qvalues with anova results
anova_results <- join_qvalues_with_ANOVA_results(qvalues)

kable(anova_results %>%
        dplyr::filter(term != "Residuals") %>%
        group_by(term) %>%
        summarise(n = sum(qvalue < 0.05)),
      caption = "Number of significant qvalues.", format = "pandoc")
```



Table: Number of significant qvalues.

term          n
---------  ----
genotype     59
tissue      429

```r
kable(anova_results %>%
        dplyr::filter(term != "Residuals") %>%
        group_by(term) %>%
        summarise(n = sum(qvalue < 0.10)),
      caption = "Number of significant q-values.", format = "pandoc")
```



Table: Number of significant q-values.

term          n
---------  ----
genotype    124
tissue      581

When using no interaction more q-values are deemed significant than when an
interaction factor is used.

## FDR correction using `fdrtool`

I want to see what the difference will be when using the `fdrtool` package 
instead of the `qvalue` package for generating FDR.


```r
library(fdrtool)

# Calculate qvalues uisng the fdrtool package
fdr_qvalues <- fdrtool(as_vector(anova_results %>%
                    dplyr::select(p.value)),
        statistic = "pvalue")
```

```
## Step 1... determine cutoff point
## Step 2... estimate parameters of null distribution and eta0
## Step 3... compute p-values and estimate empirical PDF/CDF
## Step 4... compute q-values and local fdr
## Step 5... prepare for plotting
```

![](MIDY_adipose_data_with_two_missing_values_filtered_files/figure-html/unnamed-chunk-14-1.png)<!-- -->

```r
# Create a tibble with pvalues and qvalues
fdr_qvalues <- tibble(
  p.value = fdr_qvalues[["pval"]],
  fdrtool_qvalue = fdr_qvalues[["qval"]],
  local_FDR = fdr_qvalues[["lfdr"]]
)

# Join fdrtool results with previous anova results (qvalues generated with qvalue)
anova_results <- full_join(anova_results, fdr_qvalues, by = "p.value")

# Print table with significant differences
kable(anova_results %>%
        dplyr::filter(term != "Residuals") %>%
        group_by(term) %>%
        summarise(n_fdrtool = sum(fdrtool_qvalue < 0.05),
                  n_qvalue = sum(qvalue < 0.05)),
      caption = "Number of significant qvalues at < 0.05.", format = "pandoc",
      align = "c")
```



Table: Number of significant qvalues at < 0.05.

   term      n_fdrtool    n_qvalue 
----------  -----------  ----------
 genotype       55           59    
  tissue        418         429    

```r
# Plot qvalues from the two fdr packages to see correlation
anova_results %>% 
  ggplot(aes(qvalue, fdrtool_qvalue)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ term)
```

![](MIDY_adipose_data_with_two_missing_values_filtered_files/figure-html/unnamed-chunk-14-2.png)<!-- -->

There doesn't seem to be much difference between the results between the two
packages, I will stick to using the `qvalue` package right now.

# t-test



## Mesenterial


```r
# t_test for mesenterial tissue
t_test_m <- multiple_ttests(adipose_imputed, "genotype")
```


```r
t_test_m %>%
  ggplot(aes(x = p_value)) +
  geom_histogram(breaks = seq(0, 1, 0.05), color = "black", fill = "white")
```

![](MIDY_adipose_data_with_two_missing_values_filtered_files/figure-html/unnamed-chunk-17-1.png)<!-- -->


```r
q_values_ttest_m <- qvalue(t_test_m$p_value)
q_values_ttest_m <- tibble(p_value = q_values_ttest_m[["pvalues"]],
                           q_value = q_values_ttest_m[["qvalues"]])
```

