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

## Calculate mean


```r
adipose_means <- adipose_imputed %>% 
  group_by(name, tissue, genotype) %>% 
  summarise(
    mean = mean(LFQ)
  ) %>% 
  unite(group, tissue, genotype, sep = "_") %>%
  spread(key = "group", value = "mean") %>% 
  mutate(sc_diff = sc_midy - sc_wt)

adipose_means %>% 
  ggplot(aes(m_wt, m_midy, color = sc_diff)) +
  geom_point() +
  scale_color_viridis_c()
```

![](MIDY_adipose_data_with_two_missing_values_filtered_files/figure-html/unnamed-chunk-9-1.png)<!-- -->

```r
adipose_imputed %>% 
  group_by(name, tissue) %>% 
  summarise(
    mean = mean(LFQ)
  ) %>% 
  spread(key = "tissue", value = "mean")
```

```
## # A tibble: 1,959 x 3
## # Groups:   name [1,959]
##    name        m    sc
##    <chr>   <dbl> <dbl>
##  1 A1BG     30.3  30.4
##  2 A2M      32.6  31.8
##  3 AACS     25.9  25.0
##  4 ABCA9    23.4  24.1
##  5 ABHD14B  28.4  28.3
##  6 ABHD5    24.2  24.6
##  7 ABI3BP   25.4  28.7
##  8 ABLIM1   23.1  23.5
##  9 ACAA1    28.9  28.1
## 10 ACAA2    27.7  27.5
## # ... with 1,949 more rows
```

```r
adipose_imputed %>% 
  group_by(name, genotype) %>% 
  summarise(
    mean = mean(LFQ)
  ) %>% 
  spread(key = "genotype", value = "mean")
```

```
## # A tibble: 1,959 x 3
## # Groups:   name [1,959]
##    name     midy    wt
##    <chr>   <dbl> <dbl>
##  1 A1BG     30.2  30.5
##  2 A2M      32.2  32.3
##  3 AACS     26.1  24.8
##  4 ABCA9    23.9  23.6
##  5 ABHD14B  28.4  28.3
##  6 ABHD5    24.7  24.2
##  7 ABI3BP   27.0  27.1
##  8 ABLIM1   23.1  23.5
##  9 ACAA1    28.7  28.3
## 10 ACAA2    27.8  27.4
## # ... with 1,949 more rows
```



# 2-way ANOVA

Perform a 2-way ANOVA on the imputed dataset.




```r
adipose_anova <- multiple_2way_anova_in_tidyverse_with_tukeyHSD(adipose_imputed, name, "LFQ", "genotype", "tissue")
```

Unnest anova list-column.




```r
# Extract p-values and more
anova_results <- adipose_anova %>%
  unnest(tidy_anova) %>% 
  dplyr::filter(term != "Residuals")

anova_results
```

```
## # A tibble: 5,877 x 7
##    name  term               df     sumsq    meansq statistic  p.value
##    <chr> <chr>           <dbl>     <dbl>     <dbl>     <dbl>    <dbl>
##  1 A1BG  genotype            1 0.297     0.297      2.24     0.154   
##  2 A1BG  tissue              1 0.124     0.124      0.931    0.349   
##  3 A1BG  genotype:tissue     1 0.0000599 0.0000599  0.000451 0.983   
##  4 A2M   genotype            1 0.113     0.113      0.921    0.352   
##  5 A2M   tissue              1 3.08      3.08      25.1      0.000129
##  6 A2M   genotype:tissue     1 0.138     0.138      1.12     0.305   
##  7 AACS  genotype            1 7.93      7.93       3.76     0.0704  
##  8 AACS  tissue              1 4.05      4.05       1.92     0.185   
##  9 AACS  genotype:tissue     1 2.01      2.01       0.952    0.344   
## 10 ABCA9 genotype            1 0.475     0.475      1.12     0.305   
## # ... with 5,867 more rows
```

```r
# Count number of significant p-values
make_significant_ANOVA_kable(anova_results)
```



Table: Number of significant differences per factor.

      term           n  
-----------------  -----
    genotype        186 
 genotype:tissue    101 
     tissue         689 

```r
# Plot pvalue histograms for each factor
plot_pvalue_histogram_ANOVA(anova_results, p.value, 750)
```

![](MIDY_adipose_data_with_two_missing_values_filtered_files/figure-html/anova_results-1.png)<!-- -->


```r
# Function
fdr_correction_with_qvalue <- function(tb, p_value_col) {
  
  # Quote
  p_value_col <- enquo(p_value_col)
  
  # fdr correction
  qvalues <- qvalue(as_vector(tb %>%
                                ungroup() %>% # To make sure to only get one column
                                dplyr::select(!!p_value_col)))
  
  return(qvalues)
}

# fdr correction with the qvalue package
qvalues <- fdr_correction_with_qvalue(anova_results, p.value)

# Plot fdr statistics
plot(qvalues)
```

![](MIDY_adipose_data_with_two_missing_values_filtered_files/figure-html/unnamed-chunk-12-1.png)<!-- -->

```r
hist(qvalues)
```

![](MIDY_adipose_data_with_two_missing_values_filtered_files/figure-html/unnamed-chunk-12-2.png)<!-- -->

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
genotype            34
genotype:tissue      5
tissue             338

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
genotype            67
genotype:tissue     24
tissue             448

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
 genotype    187 
  tissue     690 

```r
# Plot pvalue histograms for each factor
plot_pvalue_histogram_ANOVA(anova_results, p.value, 750)
```

![](MIDY_adipose_data_with_two_missing_values_filtered_files/figure-html/no_interaction1-1.png)<!-- -->

The p-value histograms doesn't look any different, but with fewer p-values to
correct maybe more q-values will be significant.


```r
# fdr correction with the qvalue package
qvalues <- fdr_correction_with_qvalue(anova_results, p.value)

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
genotype     52
tissue      423

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
genotype    120
tissue      578

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

![](MIDY_adipose_data_with_two_missing_values_filtered_files/figure-html/unnamed-chunk-15-1.png)<!-- -->

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
 genotype       49           52    
  tissue        422         423    

```r
# Plot qvalues from the two fdr packages to see correlation
anova_results %>% 
  ggplot(aes(qvalue, fdrtool_qvalue)) +
  geom_point() +
  geom_smooth(method = "lm") +
  facet_wrap(~ term)
```

![](MIDY_adipose_data_with_two_missing_values_filtered_files/figure-html/unnamed-chunk-15-2.png)<!-- -->

There doesn't seem to be much difference between the results between the two
packages, I will stick to using the `qvalue` package right now.

# t-test

Function for performing multiple t-tests.



## Mesenterial

First select mesenterial tissue.


```r
mesenterial_data <- adipose_imputed %>%
  dplyr::filter(tissue == "m")
```

### t-test, mesenterial, unequal variance


```r
# t_test for mesenterial tissue
t_test_m <- multiple_ttests(mesenterial_data, genotype, equal_var = FALSE)
```


```r
plot_pvalue_histogram <- function(tb, col_name, limit_max, plot_title = "p-value histogram") {
  
  # Quote
  col_name <- enquo(col_name)
  
  # pvalue plot
  tb %>% 
    ggplot(aes(x = !!col_name)) +
    geom_histogram(breaks = seq(0, 1, 0.05), fill = "black", color = "white") +
    xlab("p-value") +
    ylab("Number of proteins") +
    scale_y_continuous(limits = c(0, limit_max), expand = c(0, 0)) + 
    labs(
      title = plot_title
    ) +
    theme(panel.background = element_blank(),
          axis.line = element_line(color = "black"))
}
```

Plot p-value histogram.


```r
plot_pvalue_histogram(t_test_m, p_value, 150)
```

![](MIDY_adipose_data_with_two_missing_values_filtered_files/figure-html/unnamed-chunk-20-1.png)<!-- -->

Use `qvalue` to calculate q-values.


```r
fdr_correction_with_qvalue <- function(tb, p_value_col) {
  
  # Quote
  p_value_col <- enquo(p_value_col)
  
  # fdr correction
  qvalues <- qvalue(as_vector(tb %>%
                                ungroup() %>% # To make sure to only get one column
                                dplyr::select(!!p_value_col)))
  
  return(qvalues)
}
```


```r
q_values <- fdr_correction_with_qvalue(t_test_m, p_value)

plot(q_values)
```

![](MIDY_adipose_data_with_two_missing_values_filtered_files/figure-html/unnamed-chunk-22-1.png)<!-- -->

```r
hist(q_values)
```

![](MIDY_adipose_data_with_two_missing_values_filtered_files/figure-html/unnamed-chunk-22-2.png)<!-- -->


```r
join_qvalues_with_results <- function(tb, p_value_col = "p_value") {
  
  # Make a tibble with p and qvalues
  q_values <- tibble(p_value = q_values[["pvalues"]], q_value = q_values[["qvalues"]])
  
  # Join with ANOVA results
  tb <- full_join(tb, q_values, by = p_value_col)
  
  return(tb)
}
```


```r
t_test_m_unequal <- join_qvalues_with_results(t_test_m)
```

__MA plot__


```r
plot_MA_plot <- function(tb, plot_title = "MA plot") {
  tb %>% 
    ggplot(aes(mean_total, log2_difference)) +
    geom_point() +
    geom_hline(yintercept = 0, color = "red") +
    xlab("Mean intensity") +
    ylab("Log2 fold difference") +
    labs(
      title = plot_title
    ) +
    theme(
      panel.background = element_blank(),
      axis.line = element_line("black")) 
}
```


```r
plot_MA_plot(t_test_m_unequal)
```

![](MIDY_adipose_data_with_two_missing_values_filtered_files/figure-html/unnamed-chunk-26-1.png)<!-- -->

__Volcano plot__


```r
plot_volcano <- function(tb, plot_title = "Volcano plot") {
  tb %>% 
    ggplot(aes(log2_difference, -log10(p_value))) +
    geom_point() +
    xlab("Log2 fold difference") +
    ylab("-log10(p-value") +
    labs(
      title = plot_title
    ) +
    theme(
      panel.background = element_blank(),
      axis.line = element_line("black")) 
}
```


```r
plot_volcano(t_test_m_unequal)
```

![](MIDY_adipose_data_with_two_missing_values_filtered_files/figure-html/unnamed-chunk-28-1.png)<!-- -->


### t-test, mesenterial, equal variance


```r
# t_test for mesenterial tissue
t_test_m_equal <- multiple_ttests(mesenterial_data, genotype, equal_var = TRUE)
```


```r
plot_pvalue_histogram(t_test_m_equal, p_value, 150)
```

![](MIDY_adipose_data_with_two_missing_values_filtered_files/figure-html/unnamed-chunk-30-1.png)<!-- -->

Use `qvalue` to calculate q-values.



```r
q_values <- fdr_correction_with_qvalue(t_test_m_equal, p_value)

plot(q_values)
```

```
## geom_path: Each group consists of only one observation. Do you need to
## adjust the group aesthetic?
## geom_path: Each group consists of only one observation. Do you need to
## adjust the group aesthetic?
## geom_path: Each group consists of only one observation. Do you need to
## adjust the group aesthetic?
```

![](MIDY_adipose_data_with_two_missing_values_filtered_files/figure-html/unnamed-chunk-31-1.png)<!-- -->

```r
hist(q_values)
```

![](MIDY_adipose_data_with_two_missing_values_filtered_files/figure-html/unnamed-chunk-31-2.png)<!-- -->


```r
t_test_m_equal <- join_qvalues_with_results(t_test_m_equal)
```

### Subcutaneous

First select subcutaneous tissue.


```r
sc_data <- adipose_imputed %>%
  dplyr::filter(tissue == "sc")
```

### t-test, sc, unequal variance


```r
t_test_sc_unequal <- multiple_ttests(sc_data, genotype, equal_var = FALSE)
```

