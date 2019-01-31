---
title: "MIDY adipose tissue - 742 removed"
author: "Erik Ländström"
date: "21 January 2019"
output: 
  html_document:
    toc: TRUE
    toc_float: TRUE
    keep_md: TRUE
---





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
library(kableExtra)
library(readxl)
library(gridExtra)
library(visdat)
library(qvalue)
library(ggfortify)
library(broom)
library(pheatmap)
library(modelr)
```

# 1. Preprocessing

## 1.1 Read rawdata from Maxquant




```r
make_tibble_for_DEP("proteinGroups.txt")
```

Remove sample 742 since it has a insulinoma tumor.


```r
# The column names containing 742
str_extract(colnames(data_raw), ".+742.+") 
```

```
##  [1] NA                          NA                         
##  [3] NA                          NA                         
##  [5] NA                          NA                         
##  [7] NA                          NA                         
##  [9] NA                          NA                         
## [11] NA                          NA                         
## [13] NA                          NA                         
## [15] NA                          NA                         
## [17] NA                          NA                         
## [19] "LFQ intensity 742_m_midy"  "LFQ intensity 742_sc_midy"
## [21] NA                          NA                         
## [23] NA                          NA                         
## [25] NA                          NA                         
## [27] NA                          NA                         
## [29] NA
```

```r
# Remove columns with data from sample 742
data_raw <- data_raw %>% 
  dplyr::select(-`LFQ intensity 742_m_midy`, -`LFQ intensity 742_sc_midy`)
```

## 1.2 Remove contaminant and reverse hits

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
## [1] 2779   27
```

## 1.3 Add gene names from BiomaRt


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

# Check if there are any duplicated IDs
data_unique %>%
  dplyr::select(name, ID) %>% 
  duplicated() %>%
  any()
```

```
## [1] FALSE
```


```r
# Grep for lfq columns
lfq_columns <- grep("LFQ.", colnames(data_unique))

# Read experimental design
experimental_design <- read_xlsx("MIDY_adipose_experimental_design.xlsx")
experimental_design
```

```
## # A tibble: 20 x 3
##    label       condition replicate
##    <chr>       <chr>         <dbl>
##  1 736_m_wt    m_wt            736
##  2 736_sc_wt   sc_wt           736
##  3 737_m_midy  m_midy          737
##  4 737_sc_midy sc_midy         737
##  5 738_m_wt    m_wt            738
##  6 738_sc_wt   sc_wt           738
##  7 739_m_midy  m_midy          739
##  8 739_sc_midy sc_midy         739
##  9 740_m_midy  m_midy          740
## 10 740_sc_midy sc_midy         740
## 11 741_m_wt    m_wt            741
## 12 741_sc_wt   sc_wt           741
## 13 742_m_midy  m_midy          742
## 14 742_sc_midy sc_midy         742
## 15 743_m_wt    m_wt            743
## 16 743_sc_wt   sc_wt           743
## 17 744_m_midy  m_midy          744
## 18 744_sc_midy sc_midy         744
## 19 745_m_wt    m_wt            745
## 20 745_sc_wt   sc_wt           745
```

```r
# Remove 742 samples from experimental design
experimental_design <- experimental_design %>% 
  dplyr::filter(replicate != 742)

# Generate summarized experiment
data_se <- make_se(data_unique, lfq_columns, experimental_design)

data_se
```

```
## class: SummarizedExperiment 
## dim: 2779 18 
## metadata(0):
## assays(1): ''
## rownames(2779): NP_001001534 DES ... TIMM9 L3HYPDH
## rowData names(13): Protein IDs Majority protein IDs ... name ID
## colnames(18): m_wt_736 sc_wt_736 ... m_wt_745 sc_wt_745
## colData names(4): label ID condition replicate
```

## Extract protein name from fasta headers


```r
protein_names <- data_unique %>% 
  dplyr::select(name, `Fasta headers`) %>% 
  rename(protein_name = `Fasta headers`)

# Extract the protein names from fasta headers
protein_names <- protein_names %>% 
  mutate(protein_name = str_extract(protein_name, "^[^\\[]+"))
```


# 2. Analysis

## 2.1 Missing value filtering

Plot protein identification overlap plot.


```r
plot_frequency(data_se)
```

![](MIDY_adipose_tissue_742_removed_files/figure-html/frequency-1.png)<!-- -->

I think there are quite a few protein that have been identified in 0 samples, 
since I removed two samples (742).




```r
data_unfiltered <- get_df_wide(data_se)

g_test(data_unfiltered)
```

I filtered for three missing values before and from the g-test, I will try with
bought 2 and 3 missing values.

Since it's not possible to filter on identified values in `DEP`, I will use the
values from the g-test.


```r
# Filter for 2 valid values
filtered_names_2_vv <- g_test %>%
  dplyr::filter(m_midy >= 2 | m_wt >= 2 | sc_midy >= 2 | sc_wt >= 2) %>% 
  dplyr::select(name)

# Join with data_unique
data_filtered <- left_join(filtered_names_2_vv, data_unique, by = "name")

# Make a summarized experiment
lfq_columns <- grep("LFQ.", colnames(data_filtered))
data_se <- make_se(data_filtered, lfq_columns, experimental_design)
```



```r
plot_missval(data_se)
```

![](MIDY_adipose_tissue_742_removed_files/figure-html/unnamed-chunk-13-1.png)<!-- -->


```r
plot_detect(data_se)
```

![](MIDY_adipose_tissue_742_removed_files/figure-html/unnamed-chunk-14-1.png)<!-- -->


```r
plot_numbers(data_se)
```

![](MIDY_adipose_tissue_742_removed_files/figure-html/unnamed-chunk-15-1.png)<!-- -->


```r
plot_coverage(data_se)
```

![](MIDY_adipose_tissue_742_removed_files/figure-html/unnamed-chunk-16-1.png)<!-- -->

## 2.2 Normalization


```r
data_norm <- normalize_vsn(data_se)
meanSdPlot(data_se)$gg + ggtitle("2 valid values")
```

![](MIDY_adipose_tissue_742_removed_files/figure-html/unnamed-chunk-17-1.png)<!-- -->![](MIDY_adipose_tissue_742_removed_files/figure-html/unnamed-chunk-17-2.png)<!-- -->

```r
meanSdPlot(data_norm)$gg + ggtitle("2 valid values normalized")
```

![](MIDY_adipose_tissue_742_removed_files/figure-html/unnamed-chunk-17-3.png)<!-- -->![](MIDY_adipose_tissue_742_removed_files/figure-html/unnamed-chunk-17-4.png)<!-- -->


```r
p1 <- plot_normalization(data_se)
p2 <- plot_normalization(data_norm)
grid.arrange(p1, p2)
```

![](MIDY_adipose_tissue_742_removed_files/figure-html/unnamed-chunk-18-1.png)<!-- -->

```r
plot_normalization(data_se)
```

![](MIDY_adipose_tissue_742_removed_files/figure-html/unnamed-chunk-18-2.png)<!-- -->

```r
plot_normalization(data_norm)
```

![](MIDY_adipose_tissue_742_removed_files/figure-html/unnamed-chunk-18-3.png)<!-- -->


```r
# Imputation
imputed_data <- impute(data_norm, fun = "man", shift = 1.8, scale = 0.3)
```


```r
plot_imputation(data_norm, imputed_data)
```

![](MIDY_adipose_tissue_742_removed_files/figure-html/unnamed-chunk-19-1.png)<!-- -->

## 2.3 Save normalized data
Save normalized data as tibble for further analysis with `naniar` and more.


```r
# Save unfiltered data that has been run through DEP
dep_data_unmodified_wide <- get_df_wide(data_se)

write_tsv(dep_data_unmodified_wide, "MIDY_adipose_no742_filtered_prots_190121.txt")
```



```r
# Save imputed data in a tibble
adipose_imputed <- get_df_long(imputed_data)

# Tidy up the condition column
adipose_imputed <- adipose_imputed %>%
  separate(condition, into = c("tissue", "genotype"), sep = "_") %>%
  rename(LFQ = intensity)
```

## 2.4 Venn diagram

Create matrix for plotting venn diagram.




```r
create_long_observed_matrix_tibbles(dep_data_unmodified_wide)
plot_venn_diagram(original_data_name_matrix, "190122_MIDY_adipose_no742_venn_diagram.tiff", "MIDY adipose tissue (without 742) with 2 valid values")
```

```
## Loading required package: grid
```

```
## Loading required package: futile.logger
```

```
## [1] 1
```
![](190122_MIDY_adipose_no742_venn_diagram.tiff)

## 2.5 PCA


```r
# Get wide imputed tibble
wide_imputed_tibble <- get_df_wide(imputed_data)
```




```r
# Perform PCA
pca_imputed_data <- format_data_for_PCA(wide_imputed_tibble, 18, c("tissue", "genotype", "sample"))
```


```r
# Calculate variance explained
var_explained <- calculate_variance_explained(pca_imputed_data)

# Calculate variance explained
var_explained <- pca_imputed_data %>%
  unnest(pca_aug) %>%
  summarise_at(.vars = vars(contains("fittedPC")), .funs = funs(var)) %>%
  gather(key = pc, value = variance) %>%
  mutate(
    var_exp     = variance / sum(variance),
    cum_var_exp = cumsum(var_exp),
    pc          = str_replace(pc, ".fittedPC", "")
  )


# Print variance explained table
make_var_exp_kable(var_explained, "html")
```

<table class="table table-striped table-hover" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>Proportion of the variance that each principal component explains</caption>
 <thead>
  <tr>
   <th style="text-align:center;"> PC </th>
   <th style="text-align:center;"> Variance </th>
   <th style="text-align:center;"> Proportion explained </th>
   <th style="text-align:center;"> Cumulative proportion explained </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;"> 1 </td>
   <td style="text-align:center;"> 443.4420 </td>
   <td style="text-align:center;"> 0.2285 </td>
   <td style="text-align:center;"> 0.2285 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 2 </td>
   <td style="text-align:center;"> 236.1344 </td>
   <td style="text-align:center;"> 0.1217 </td>
   <td style="text-align:center;"> 0.3501 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 3 </td>
   <td style="text-align:center;"> 170.8470 </td>
   <td style="text-align:center;"> 0.0880 </td>
   <td style="text-align:center;"> 0.4381 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 4 </td>
   <td style="text-align:center;"> 130.7735 </td>
   <td style="text-align:center;"> 0.0674 </td>
   <td style="text-align:center;"> 0.5055 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 5 </td>
   <td style="text-align:center;"> 123.6906 </td>
   <td style="text-align:center;"> 0.0637 </td>
   <td style="text-align:center;"> 0.5692 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 6 </td>
   <td style="text-align:center;"> 103.1267 </td>
   <td style="text-align:center;"> 0.0531 </td>
   <td style="text-align:center;"> 0.6224 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 7 </td>
   <td style="text-align:center;"> 93.0788 </td>
   <td style="text-align:center;"> 0.0480 </td>
   <td style="text-align:center;"> 0.6703 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 8 </td>
   <td style="text-align:center;"> 84.7829 </td>
   <td style="text-align:center;"> 0.0437 </td>
   <td style="text-align:center;"> 0.7140 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 9 </td>
   <td style="text-align:center;"> 81.0859 </td>
   <td style="text-align:center;"> 0.0418 </td>
   <td style="text-align:center;"> 0.7558 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 10 </td>
   <td style="text-align:center;"> 71.8647 </td>
   <td style="text-align:center;"> 0.0370 </td>
   <td style="text-align:center;"> 0.7928 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 11 </td>
   <td style="text-align:center;"> 64.9454 </td>
   <td style="text-align:center;"> 0.0335 </td>
   <td style="text-align:center;"> 0.8263 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 12 </td>
   <td style="text-align:center;"> 63.7053 </td>
   <td style="text-align:center;"> 0.0328 </td>
   <td style="text-align:center;"> 0.8591 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 13 </td>
   <td style="text-align:center;"> 60.2867 </td>
   <td style="text-align:center;"> 0.0311 </td>
   <td style="text-align:center;"> 0.8901 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 14 </td>
   <td style="text-align:center;"> 56.7686 </td>
   <td style="text-align:center;"> 0.0292 </td>
   <td style="text-align:center;"> 0.9194 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 15 </td>
   <td style="text-align:center;"> 56.1059 </td>
   <td style="text-align:center;"> 0.0289 </td>
   <td style="text-align:center;"> 0.9483 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 16 </td>
   <td style="text-align:center;"> 54.0657 </td>
   <td style="text-align:center;"> 0.0279 </td>
   <td style="text-align:center;"> 0.9761 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 17 </td>
   <td style="text-align:center;"> 46.2958 </td>
   <td style="text-align:center;"> 0.0239 </td>
   <td style="text-align:center;"> 1.0000 </td>
  </tr>
  <tr>
   <td style="text-align:center;"> 18 </td>
   <td style="text-align:center;"> 0.0000 </td>
   <td style="text-align:center;"> 0.0000 </td>
   <td style="text-align:center;"> 1.0000 </td>
  </tr>
</tbody>
</table>

```r
# Plot variance explained graphs
plot_var_exp(var_explained)
```

![](MIDY_adipose_tissue_742_removed_files/figure-html/unnamed-chunk-24-1.png)<!-- -->

From these plots I would say that the "elbow" is somewhere between the third and 
fifth principal component.




```r
# Plot first two PCs
plot_pca(pca_imputed_data, "First two principal components of PCA performed on MIDY adipose tissue")
```

![](MIDY_adipose_tissue_742_removed_files/figure-html/unnamed-chunk-25-1.png)<!-- -->

The first two PCs separates both tissues pretty well, and to a lesser extent
the genotype groups.

I'd like to do plot of more of the PCs.




```r
pca_1 <- plot_multiple_PCAs(pca_imputed_data, 1, 5)
pca_1
```

```
## [[1]]
```

![](MIDY_adipose_tissue_742_removed_files/figure-html/echo-1.png)<!-- -->

```
## 
## [[2]]
```

![](MIDY_adipose_tissue_742_removed_files/figure-html/echo-2.png)<!-- -->

```
## 
## [[3]]
```

![](MIDY_adipose_tissue_742_removed_files/figure-html/echo-3.png)<!-- -->

```
## 
## [[4]]
```

![](MIDY_adipose_tissue_742_removed_files/figure-html/echo-4.png)<!-- -->


```r
pca_2 <- plot_multiple_PCAs(pca_imputed_data, 2, 5)
pca_2
```

```
## [[1]]
```

![](MIDY_adipose_tissue_742_removed_files/figure-html/unnamed-chunk-27-1.png)<!-- -->

```
## 
## [[2]]
```

![](MIDY_adipose_tissue_742_removed_files/figure-html/unnamed-chunk-27-2.png)<!-- -->

```
## 
## [[3]]
```

![](MIDY_adipose_tissue_742_removed_files/figure-html/unnamed-chunk-27-3.png)<!-- -->



```r
pca_3 <- plot_multiple_PCAs(pca_imputed_data, 3, 5)
pca_3
```

```
## [[1]]
```

![](MIDY_adipose_tissue_742_removed_files/figure-html/unnamed-chunk-28-1.png)<!-- -->

```
## 
## [[2]]
```

![](MIDY_adipose_tissue_742_removed_files/figure-html/unnamed-chunk-28-2.png)<!-- -->

These plots suggest that there might be two potential outliers __m_wt_743__
and __sc_wt_741__.

## 2.6 Heatmap


```r
pheatmap(imputed_tibble_for_pca[, 2:18], method = "ward.D2", 
         show_rownames = FALSE)

pheatmap(imputed_tibble_for_pca[, 2:18], method = "ward.D2", 
         show_rownames = FALSE,
         color = colorRampPalette(c("red", "white", "blue"))(100))
```


# 3 2-way ANOVA

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
## # A tibble: 5,823 x 7
##    name  term               df   sumsq  meansq statistic  p.value
##    <chr> <chr>           <dbl>   <dbl>   <dbl>     <dbl>    <dbl>
##  1 A1BG  genotype            1 0.685   0.685      6.62   0.0221  
##  2 A1BG  tissue              1 0.0999  0.0999     0.965  0.343   
##  3 A1BG  genotype:tissue     1 0.00116 0.00116    0.0112 0.917   
##  4 A2M   genotype            1 0.289   0.289      2.24   0.157   
##  5 A2M   tissue              1 3.01    3.01      23.2    0.000272
##  6 A2M   genotype:tissue     1 0.144   0.144      1.11   0.310   
##  7 AACS  genotype            1 5.64    5.64       2.86   0.113   
##  8 AACS  tissue              1 1.70    1.70       0.861  0.369   
##  9 AACS  genotype:tissue     1 0.393   0.393      0.200  0.662   
## 10 ABCA9 genotype            1 0.178   0.178      0.440  0.518   
## # ... with 5,813 more rows
```

```r
# Count number of significant p-values
make_significant_ANOVA_kable(anova_results)
```



Table: Number of significant differences per factor.

      term           n  
-----------------  -----
    genotype        178 
 genotype:tissue    78  
     tissue         624 

```r
# Plot pvalue histograms for each factor
plot_pvalue_histogram_ANOVA(anova_results, p.value, 750)
```

![](MIDY_adipose_tissue_742_removed_files/figure-html/anova_results-1.png)<!-- -->




```r
# fdr correction with the qvalue package
qvalues <- fdr_correction_with_qvalue(anova_results, p.value)

# Plot fdr statistics
plot(qvalues)
```

![](MIDY_adipose_tissue_742_removed_files/figure-html/unnamed-chunk-33-1.png)<!-- -->

```r
hist(qvalues)
```

![](MIDY_adipose_tissue_742_removed_files/figure-html/unnamed-chunk-33-2.png)<!-- -->




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

term                 n
----------------  ----
genotype            18
genotype:tissue      3
tissue             261

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
genotype            53
genotype:tissue     15
tissue             390

## 3.1 Tukey results


```r
tukey_results <- adipose_anova %>%
  unnest(tidy_tukey)
```

# 4 2-way ANOVA without interaction effects

Since there doesn't seem to be many interaction effect, I want to see what 
will happen if you do an ANOVA without the interaction effect between genotype
and tissue.




```r
# Perform ANOVA without interaction effect
adipose_anova_without_interaction <- multiple_2way_anova_in_tidyverse_with_no_interaction(adipose_imputed)

# Extract p-values and more
anova_results_without_interaction <- adipose_anova_without_interaction %>%
  unnest(tidy_anova) %>% 
  dplyr::filter(term != "Residuals")

# Count number of significant p-values and print a kable
make_significant_ANOVA_kable(anova_results_without_interaction)
```



Table: Number of significant differences per factor.

   term       n  
----------  -----
 genotype    178 
  tissue     634 

```r
# Plot pvalue histograms for each factor
plot_pvalue_histogram_ANOVA(anova_results_without_interaction, p.value, 750)
```

![](MIDY_adipose_tissue_742_removed_files/figure-html/no_interaction1-1.png)<!-- -->


```r
# fdr correction with the qvalue package
qvalues <- fdr_correction_with_qvalue(anova_results_without_interaction, p.value)

# Plot fdr statistics
plot(qvalues)
```

![](MIDY_adipose_tissue_742_removed_files/figure-html/no_interaction2-1.png)<!-- -->

```r
hist(qvalues)
```

![](MIDY_adipose_tissue_742_removed_files/figure-html/no_interaction2-2.png)<!-- -->

```r
# Join qvalues with anova results
qvalues <- tibble(p.value = qvalues[["pvalues"]], qvalue = qvalues[["qvalues"]])
anova_results_without_interaction <- full_join(anova_results_without_interaction, qvalues, by = "p.value")

kable(anova_results_without_interaction %>%
        dplyr::filter(term != "Residuals") %>%
        group_by(term) %>%
        summarise(n = sum(qvalue < 0.05)),
      caption = "Number of significant qvalues.", format = "pandoc")
```



Table: Number of significant qvalues.

term          n
---------  ----
genotype     47
tissue      385

```r
kable(anova_results_without_interaction %>%
        dplyr::filter(term != "Residuals") %>%
        group_by(term) %>%
        summarise(n = sum(qvalue < 0.10)),
      caption = "Number of significant q-values.", format = "pandoc")
```



Table: Number of significant q-values.

term          n
---------  ----
genotype    107
tissue      510
## 4.1 Calculate means


```r
# Means for each group
adipose_means <- adipose_imputed %>% 
  group_by(name, tissue, genotype) %>% 
  summarise(
    mean = mean(LFQ)
  ) %>% 
  unite(group, tissue, genotype, sep = "_") %>%
  spread(key = "group", value = "mean") %>% 
  mutate(sc_diff = sc_midy - sc_wt,
         m_diff  = m_midy  - m_wt)

# Means for each tissue group
tissue_means <- adipose_imputed %>% 
  group_by(name, tissue) %>% 
  summarise(
    mean = mean(LFQ)
  ) %>% 
  spread(key = "tissue", value = "mean")

# Means for each genotype group
genotype_means <- adipose_imputed %>% 
  group_by(name, genotype) %>% 
  summarise(
    mean = mean(LFQ)
  ) %>% 
  spread(key = "genotype", value = "mean")

# join mean tables
adipose_means <- full_join(adipose_means[, 1:5], tissue_means, by = "name")
adipose_means <- full_join(adipose_means, genotype_means, by = "name")
```

Join means with qvalue data.


```r
# Function for spreading two columns
spread_nt <- function(data,key_col,...,fill = NA,
                      convert = TRUE,drop = TRUE,sep = "_"){
  key_quo <- rlang::enquo(key_col)
  val_quos <- rlang::quos(...)
  value_cols <- unname(tidyselect::vars_select(names(data),!!!val_quos))
  key_col <- unname(tidyselect::vars_select(names(data),!!key_quo))

  data %>%
    gather(key = ..var..,value = ..val..,!!!val_quos) %>%
    unite(col = ..grp..,c(key_col,"..var.."),sep = sep) %>%
    spread(key = ..grp..,value = ..val..,fill = fill,
           convert = convert,drop = drop,sep = NULL)
}

# Spread pvalue, qvalues and join with means
anova_results_means_withhout_interaction <- anova_results_without_interaction %>% 
  dplyr::select(name, term, p.value, qvalue) %>% 
  spread_nt(key_col = term, p.value, qvalue) %>% 
  full_join(adipose_means, by = "name")
```

## 4.2 Filter for significant q-values


```r
genotype_anova_results_without_interaction_qval_0.05 <- anova_results_means_withhout_interaction %>% 
  dplyr::filter(genotype_qvalue < 0.05) %>% 
  dplyr::select(name, genotype_p.value, genotype_qvalue, midy, wt) %>% 
  mutate(l2fc = midy - wt) %>% 
  arrange(desc(l2fc))

genotype_anova_results_without_interaction_qval_0.05 <- genotype_anova_results_without_interaction_qval_0.05 %>% 
  left_join(protein_names, by = "name") %>% 
  rename(`Gene name` = name, 
         `Genotype p-value` = genotype_p.value,
         `Genotype q-value` = genotype_qvalue, 
         `MIDY mean` = midy, 
         `WT mean` = wt, 
         `log2 fold change` = l2fc,
         `Protein name` = protein_name)

# genotype_anova_results_without_interaction_qval_0.05 %>% 
#   write_tsv("genotype_anova_results_without_interaction_qval_0.05.tsv")

tissue_anova_results_without_interaction_qval_0.05 <- anova_results_means_withhout_interaction %>% 
  dplyr::filter(tissue_qvalue < 0.05) %>% 
  dplyr::select(name, tissue_p.value, tissue_qvalue, sc, m) %>% 
  mutate(l2fc = sc - m) %>% 
  arrange(desc(l2fc))

tissue_anova_results_without_interaction_qval_0.05 <- tissue_anova_results_without_interaction_qval_0.05 %>% 
  left_join(protein_names, by = "name") %>% 
  rename(`Gene name` = name, 
         `Tissue p-value` = tissue_p.value,
         `Tissue q-value` = tissue_qvalue, 
         `SC mean` = sc, 
         `M mean` = m, 
         `log2 fold change` = l2fc,
         `Protein name` = protein_name)
# tissue_anova_results_without_interaction_qval_0.05 %>% 
#   write_tsv("tissue_anova_results_without_interaction_qval_0.05.tsv")

kable(tissue_anova_results_without_interaction_qval_0.05 %>% 
        summarise(`Subcutaneous adipose tissue` = sum(`log2 fold change` > 0),
                  `Mesenterial adipose tissue` = sum(`log2 fold change` < 0)),
      caption = "Number of proteins signifcantly more abundant according to 2-way ANOVA",
      format = "pandoc",
      align = "c"
      )
```



Table: Number of proteins signifcantly more abundant according to 2-way ANOVA

 Subcutaneous adipose tissue    Mesenterial adipose tissue 
-----------------------------  ----------------------------
             205                           180             

Are there proteins found to be significantly differentially expressed in both
groupings.


```r
left_join(genotype_anova_results_without_interaction_qval_0.05,
          tissue_anova_results_without_interaction_qval_0.05, by = "Gene name") %>%
  dplyr::filter(!is.na(`Tissue q-value`))
```

```
## # A tibble: 15 x 13
##    `Gene name` `Genotype p-val~ `Genotype q-val~ `MIDY mean` `WT mean`
##    <chr>                  <dbl>            <dbl>       <dbl>     <dbl>
##  1 SLC25A22          0.00316            0.0259          24.8      23.4
##  2 PRDX3             0.00394            0.0298          29.0      28.4
##  3 XP_0209515~       0.00144            0.0157          26.8      26.2
##  4 HADHA             0.000692           0.00919         29.5      29.1
##  5 HADHB             0.000814           0.0105          28.0      27.6
##  6 ACOT4             0.00171            0.0176          28.3      27.9
##  7 AIFM1             0.00727            0.0435          26.8      26.4
##  8 LNPEP             0.00659            0.0407          27.0      26.6
##  9 ACADVL            0.00402            0.0298          26.8      26.4
## 10 PSMD5             0.00118            0.0134          27.5      27.2
## 11 RAP1A             0.00378            0.0291          29.5      29.3
## 12 HSP90AA1          0.00188            0.0189          31.3      31.5
## 13 XP_0209436~       0.00348            0.0277          25.5      25.8
## 14 AGT               0.00000233         0.000231        28.4      29.1
## 15 PCOLCE            0.00584            0.0373          25.7      26.8
## # ... with 8 more variables: `log2 fold change.x` <dbl>, `Protein
## #   name.x` <chr>, `Tissue p-value` <dbl>, `Tissue q-value` <dbl>, `SC
## #   mean` <dbl>, `M mean` <dbl>, `log2 fold change.y` <dbl>, `Protein
## #   name.y` <chr>
```

15 proteins are found to be significantly differentially expressed in both 
tissue and genotype groupings.

## 4.3 Volcano plots


```r
# Calculate lfc2 for each group
adipose_lfc2 <- adipose_means %>% 
  group_by(name) %>% 
  mutate(tissue = m - sc,
         genotype = midy - wt) %>%
  dplyr::select(name, tissue, genotype) %>% 
  gather(tissue, genotype, key = "term", value = "lfc2") # term is used for easy joining with anova results

# Join with anova results
anova_results_without_interaction <- full_join(anova_results_without_interaction, adipose_lfc2,
                                               by = c("name", "term")) %>% 
  rename(group = term)

anova_results_without_interaction_for_volcano <- anova_results_without_interaction %>% 
  mutate(enriched_1 = as.numeric(lfc2 > 0 & qvalue < 0.05),
         enriched_2 = as.numeric(lfc2 < 0 & qvalue < 0.05),
         significant_temp = (enriched_1 + enriched_2 * 2) + 1,
         significant = case_when(group == "tissue" & significant_temp > 1 ~ significant_temp + 2,
                                 group == "tissue" & significant_temp == 1 ~ significant_temp,
                                 group == "genotype" ~ significant_temp),
         significant = as.factor(significant))


anova_results_without_interaction_for_volcano %>%
  ggplot(aes(lfc2, -log10(p.value), color = significant)) +
  geom_point() +
  scale_color_manual("More abundant in:",
                     values = c("black", "#00BFC4", "#F8666D", "#7CAE00", "#C77CFF"),
                     breaks = c(2, 3, 4, 5),
                     labels = c("MIDY", "WT", "Adipose, mesenterial", "Adipose, subcutaneous")) +
  facet_grid(~ group, scales = "free_x") +
  theme(panel.background = element_blank(),
        axis.line = element_line("black"))
```

![](MIDY_adipose_tissue_742_removed_files/figure-html/unnamed-chunk-42-1.png)<!-- -->

# 5 t-test

Function for performing multiple t-tests.



## 5.1 Mesenterial

First select mesenterial tissue.


```r
mesenterial_data <- adipose_imputed %>%
  dplyr::filter(tissue == "m")
```

### 5.1.1 t-test, mesenterial, unequal variance


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

![](MIDY_adipose_tissue_742_removed_files/figure-html/unnamed-chunk-47-1.png)<!-- -->

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

![](MIDY_adipose_tissue_742_removed_files/figure-html/unnamed-chunk-49-1.png)<!-- -->

```r
hist(q_values)
```

![](MIDY_adipose_tissue_742_removed_files/figure-html/unnamed-chunk-49-2.png)<!-- -->


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

## 5.2 Subcutaneous

First select subcutaneous tissue.


```r
sc_data <- adipose_imputed %>%
  dplyr::filter(tissue == "sc")
```

### 5.2.1 t-test, sc, unequal variance


```r
t_test_sc_unequal <- multiple_ttests(sc_data, genotype, equal_var = FALSE)
```

# 6 Modelling

## 6.1 Separate out the effect of tissue



Remove tissue influence on variation.


```r
lm_tissue_residuals <- linear_model_1_ind_var(adipose_imputed, name, "LFQ", "tissue")
```

### 6.1.1 t-test




```r
lm_genotype_resids_ttest <- multiple_ttests(lm_tissue_residuals, genotype, values = resid, equal_var = FALSE)

# Count number of significant proteins
kable(lm_genotype_resids_ttest %>% 
  ungroup() %>% 
  summarise(n = sum(p_value < 0.05)),
  caption = "Number of significant proteins between genotype groups (p-value < 0.05") %>% 
  kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE)
```

<table class="table table-striped table-hover" style="margin-left: auto; margin-right: auto;">
<caption>Number of significant proteins between genotype groups (p-value &lt; 0.05</caption>
 <thead>
  <tr>
   <th style="text-align:right;"> n </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> 218 </td>
  </tr>
</tbody>
</table>


```r
plot_pvalue_histogram(lm_genotype_resids_ttest, col_name = p_value, limit_max = 250)
```

![](MIDY_adipose_tissue_742_removed_files/figure-html/unnamed-chunk-58-1.png)<!-- -->


### 6.1.2 Correlation between t-test on residuals and 2-way ANOVA

I want to compare the results of the 2-way ANOVA and linear modelling with 
t-test.


```r
# Protein list with ranks according to significance (t-test)
lm_genotype_rank_list <- lm_genotype_resids_ttest %>% 
  ungroup() %>% 
  arrange(p_value) %>% 
  dplyr::select(name, p_value) %>% 
  mutate(rank = rank(p_value))

# Protein list with ranks according to significance (ANOVA, genotype)
anova_genotype_no_inter_rank_list <- anova_results_without_interaction %>% 
  dplyr::filter(group == "genotype") %>% 
  dplyr::select(name, p.value) %>% 
  arrange(p.value) %>% 
  ungroup() %>% 
  mutate(rank_anova = rank(p.value))

# Join data
genotype_pvals_ranks <- full_join(lm_genotype_rank_list, anova_genotype_no_inter_rank_list, by = "name")
```

Compare ranks between the methods.


```r
# Test correlation with a spearman rank test
cor.test(genotype_pvals_ranks$p_value, genotype_pvals_ranks$p.value,
         method ="spearman")
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  genotype_pvals_ranks$p_value and genotype_pvals_ranks$p.value
## S = 2736300, p-value < 2.2e-16
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##       rho 
## 0.9977549
```

```r
cor.test(genotype_pvals_ranks$rank, genotype_pvals_ranks$rank_anova,
         method ="spearman")
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  genotype_pvals_ranks$rank and genotype_pvals_ranks$rank_anova
## S = 2736300, p-value < 2.2e-16
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##       rho 
## 0.9977549
```

The two methods are higly correlated with a rho value of 0.9978.


```r
genotype_pvals_ranks %>% 
  ggplot(aes(p_value, p.value)) +
  geom_point() +
  geom_smooth(method = "lm")
```

![](MIDY_adipose_tissue_742_removed_files/figure-html/unnamed-chunk-61-1.png)<!-- -->

As can be seen the p-values correlate well.


```r
# Extract only signicant proteins
lm_genotype_sig_prot <- genotype_pvals_ranks %>% 
  dplyr::filter(p_value < 0.05)
anova_genotype_sig_prot <- genotype_pvals_ranks %>% 
  dplyr::filter(p.value < 0.05)

left_join(anova_genotype_sig_prot %>% 
            dplyr::select(name), lm_genotype_sig_prot, by = "name") %>% 
  summarise(n = sum(is.na(p_value)))
```

```
## # A tibble: 1 x 1
##       n
##   <int>
## 1     1
```

```r
left_join(anova_genotype_sig_prot %>% 
            dplyr::select(name), lm_genotype_sig_prot, by = "name") %>% 
  dplyr::filter(is.na(p_value))
```

```
## # A tibble: 1 x 5
##   name  p_value  rank p.value rank_anova
##   <chr>   <dbl> <dbl>   <dbl>      <dbl>
## 1 BST1       NA    NA      NA         NA
```

```r
genotype_pvals_ranks %>% 
  dplyr::filter(name == "SCD" | name == "XP_020955891")
```

```
## # A tibble: 2 x 5
##   name         p_value  rank p.value rank_anova
##   <chr>          <dbl> <dbl>   <dbl>      <dbl>
## 1 SCD           0.0248   132  0.0195         87
## 2 XP_020955891  0.300    862  0.309         825
```

There are only two proteins that are not considered to be significant with 
linear model and t-test but that were significant in the 2-way ANOVA. These two
proteins are probably not going to be significant after multiple hypotheses 
correction.

## 6.2 Separate out the effect of tissue


```r
lm_genotype_resids <- linear_model_1_ind_var(adipose_imputed, cases = name, observations = "LFQ", ind_variable = "genotype")
```

### 6.2.1 t-test




```r
lm_tissue_resids_ttest <- multiple_ttests(lm_genotype_resids, 
                                          group_var = tissue,
                                          group_1 = sc,
                                          group_2 = m,
                                          mean_1 = mean_sc,
                                          mean_2 = mean_m,
                                          values = resid)

# Count number of significant proteins
kable(lm_tissue_resids_ttest %>% 
  ungroup() %>% 
  summarise(n = sum(p_value < 0.05)), 
  caption = "Number of significant proteins found between tissue groups (p-value y 0.05)") %>% 
  kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE)
```

<table class="table table-striped table-hover" style="margin-left: auto; margin-right: auto;">
<caption>Number of significant proteins found between tissue groups (p-value y 0.05)</caption>
 <thead>
  <tr>
   <th style="text-align:right;"> n </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> 638 </td>
  </tr>
</tbody>
</table>


```r
plot_pvalue_histogram(lm_tissue_resids_ttest, col_name = p_value, limit_max = 700)
```

![](MIDY_adipose_tissue_742_removed_files/figure-html/unnamed-chunk-66-1.png)<!-- -->


### 6.2.2 Correlation between t-test on residuals and 2-way ANOVA

I want to compare the results of the 2-way ANOVA and linear modelling with 
t-test.


```r
# Protein list with ranks according to significance (t-test)
lm_tissue_rank_list <- lm_tissue_resids_ttest %>% 
  ungroup() %>% 
  arrange(p_value) %>% 
  dplyr::select(name, p_value) %>% 
  mutate(rank = rank(p_value))

# Protein list with ranks according to significance (ANOVA, genotype)
anova_tissue_no_inter_rank_list <- anova_results_without_interaction %>% 
  dplyr::filter(group == "tissue") %>% 
  dplyr::select(name, p.value) %>% 
  arrange(p.value) %>% 
  ungroup() %>% 
  mutate(rank_anova = rank(p.value))

# Join data
tissue_pvals_ranks <- full_join(lm_tissue_rank_list, anova_tissue_no_inter_rank_list, by = "name")
```

Compare ranks between the methods.


```r
# Test correlation with a spearman rank test
cor.test(tissue_pvals_ranks$p_value, tissue_pvals_ranks$p.value,
         method ="spearman")
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  tissue_pvals_ranks$p_value and tissue_pvals_ranks$p.value
## S = 437390, p-value < 2.2e-16
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##       rho 
## 0.9996411
```

```r
cor.test(tissue_pvals_ranks$rank, tissue_pvals_ranks$rank_anova,
         method ="spearman")
```

```
## 
## 	Spearman's rank correlation rho
## 
## data:  tissue_pvals_ranks$rank and tissue_pvals_ranks$rank_anova
## S = 437390, p-value < 2.2e-16
## alternative hypothesis: true rho is not equal to 0
## sample estimates:
##       rho 
## 0.9996411
```

The two methods are higly correlated with a rho value of 0.9996.


```r
tissue_pvals_ranks %>% 
  ggplot(aes(p_value, p.value)) +
  geom_point() +
  geom_smooth(method = "lm")
```

![](MIDY_adipose_tissue_742_removed_files/figure-html/unnamed-chunk-69-1.png)<!-- -->

As can be seen the p-values correlate well, even better than for genotype.


```r
# Extract only signicant proteins
lm_tissue_sig_prot <- tissue_pvals_ranks %>% 
  dplyr::filter(p_value < 0.05)
anova_tissue_sig_prot <- tissue_pvals_ranks %>% 
  dplyr::filter(p.value < 0.05)

left_join(anova_tissue_sig_prot %>% 
            dplyr::select(name), lm_tissue_sig_prot, by = "name") %>% 
  summarise(n = sum(is.na(p_value)))
```

```
## # A tibble: 1 x 1
##       n
##   <int>
## 1     4
```

```r
left_join(anova_tissue_sig_prot %>% 
            dplyr::select(name), lm_tissue_sig_prot, by = "name") %>% 
  dplyr::filter(is.na(p_value))
```

```
## # A tibble: 4 x 5
##   name     p_value  rank p.value rank_anova
##   <chr>      <dbl> <dbl>   <dbl>      <dbl>
## 1 C7            NA    NA      NA         NA
## 2 CKAP4         NA    NA      NA         NA
## 3 DLD           NA    NA      NA         NA
## 4 SERPINA7      NA    NA      NA         NA
```

```r
tissue_pvals_ranks %>% 
  dplyr::filter(name == "C7" | name == "DLD")
```

```
## # A tibble: 2 x 5
##   name  p_value  rank p.value rank_anova
##   <chr>   <dbl> <dbl>   <dbl>      <dbl>
## 1 C7     0.0515   645  0.0438        608
## 2 DLD    0.0522   649  0.0482        626
```

There are only two proteins that are not considered to be significant with 
linear model and t-test but that were significant in the 2-way ANOVA. As for 
genotype, these two proteins are probably not going to be significant after 
multiple hypotheses correction.

## 6.3 Multiple hypotheses correction

Since ther results are so similar I will pool the p-values from the two linear
models and perform multiple hypotheses testing.


```r
# Join results from the t-tests on the different variables
p_vals_lm <- full_join(lm_genotype_resids_ttest %>% 
            ungroup() %>% 
            dplyr::select(name, p_value), 
          lm_tissue_resids_ttest %>% 
            ungroup() %>% 
            dplyr::select(name, p_value),
          by = "name") %>% 
  gather(p_value.x, p_value.y, key = "type", value = p_value)

# fdr correction
qvals_lm <- fdr_correction_with_qvalue(p_vals_lm, p_value)
```



### 6.3.1 Genotype FDRs


```r
# Join qvalues with pvalues: genotype
lm_geno_resids_fdr <- join_qvalues_with_results(lm_genotype_resids_ttest, qvals_lm) %>% 
  dplyr::filter(!is.na(name))

# Print number of proteins significant at q-value < 0.05
kable(
  lm_geno_resids_fdr %>% 
    ungroup() %>% 
    summarise(n = sum(qvalue < 0.05)),
  caption = "Number of significant proteins with q-value < 0.05"
) %>% 
  kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE)
```

<table class="table table-striped table-hover" style="margin-left: auto; margin-right: auto;">
<caption>Number of significant proteins with q-value &lt; 0.05</caption>
 <thead>
  <tr>
   <th style="text-align:right;"> n </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> 63 </td>
  </tr>
</tbody>
</table>

Check whether the proteins are significant both in the ANOVA and the t-test 
following a linear model.


```r
# Extract significant proteins from ANOVA after FDR correction
anova_geno_sig_qval <- anova_results_means_withhout_interaction %>% 
  dplyr::filter(genotype_qvalue < 0.05) %>% 
  dplyr::select(name)

# Extract significant proteins from t-test after FDR correction
lm_geno_resids_fdr_sig <- lm_geno_resids_fdr %>% 
  dplyr::filter(qvalue < 0.05)

# Check whether all proteins from ANOVA are significant in the t-test as well
left_join(anova_geno_sig_qval, 
          lm_geno_resids_fdr_sig,
          by = "name") %>% 
  summarise(n = sum(is.na(p_value)))
```

```
## # A tibble: 1 x 1
##       n
##   <int>
## 1     2
```

All proteins found to be signicant in the ANOVA are found to be significant in
the t-test as well.

Which proteins are significant in the t-test but not the ANOVA?



#### 6.3.1.1 Volcano plot




```r
lm_geno_resid_fdr_volcano <- prepare_data_for_volcano_plot(lm_geno_resids_fdr, 
                                                           group_neg = "WT", group_pos = "MIDY")
```




```r
plot_fdr_volcano(lm_geno_resid_fdr_volcano, 
                 plot_title = "Volcano plot showing significant changes in protein expression\nbetween MIDY and WT pigs in adipose tissue",
                 break_vector = c("MIDY", "WT"),
                 color_vector = c("#00BFC4", "black", "#F8666D"))
```

![](MIDY_adipose_tissue_742_removed_files/figure-html/unnamed-chunk-79-1.png)<!-- -->

### 6.3.2 Tissue FDRs


```r
# Join qvalues with pvalues: tissue
lm_tissue_resids_fdr <- join_qvalues_with_results(lm_tissue_resids_ttest, qvals_lm) %>% 
  dplyr::filter(!is.na(name))

# Print number of proteins significant at q-value < 0.05

kable(
  lm_tissue_resids_fdr %>% 
    ungroup() %>% 
    summarise(n = sum(qvalue < 0.05)),
  caption = "Number of significant proteins with q-value < 0.05"
) %>% 
  kable_styling(bootstrap_options = c("striped", "hover"), full_width = TRUE)
```

<table class="table table-striped table-hover" style="margin-left: auto; margin-right: auto;">
<caption>Number of significant proteins with q-value &lt; 0.05</caption>
 <thead>
  <tr>
   <th style="text-align:right;"> n </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:right;"> 394 </td>
  </tr>
</tbody>
</table>

Check whether the proteins are significant both in the ANOVA and the t-test 
following a linear model.


```r
# Extract significant proteins from ANOVA after FDR correction
anova_tissue_sig_qval <- anova_results_means_withhout_interaction %>% 
  dplyr::filter(tissue_qvalue < 0.05) %>% 
  dplyr::select(name)

# Extract significant proteins from t-test after FDR correction
lm_tissue_resids_fdr_sig <- lm_tissue_resids_fdr %>% 
  dplyr::filter(qvalue < 0.05)

# Check whether all proteins from ANOVA are significant in the t-test as well
left_join(anova_tissue_sig_qval, 
          lm_tissue_resids_fdr_sig,
          by = "name") %>% 
  summarise(n = sum(is.na(p_value)))
```

```
## # A tibble: 1 x 1
##       n
##   <int>
## 1     4
```

```r
left_join(anova_results_means_withhout_interaction %>%
            dplyr::filter(tissue_qvalue < 0.05) %>% 
            dplyr::select(name), 
          lm_tissue_resids_fdr_sig, 
          by = "name") %>% 
  dplyr::filter(is.na(p_value))
```

```
## # A tibble: 4 x 9
##   name  m     sc    p_value mean_sc mean_m log2_difference mean_total
##   <chr> <lis> <lis>   <dbl>   <dbl>  <dbl>           <dbl>      <dbl>
## 1 CAPG  <NUL~ <NUL~      NA      NA     NA              NA         NA
## 2 CLU   <NUL~ <NUL~      NA      NA     NA              NA         NA
## 3 XP_0~ <NUL~ <NUL~      NA      NA     NA              NA         NA
## 4 XP_0~ <NUL~ <NUL~      NA      NA     NA              NA         NA
## # ... with 1 more variable: qvalue <dbl>
```

```r
anova_results_means_withhout_interaction %>% 
  dplyr::filter(name == "FKBP5" | name == "SLC5A10")
```

```
## # A tibble: 2 x 13
##   name  genotype_p.value genotype_qvalue tissue_p.value tissue_qvalue
##   <chr>            <dbl>           <dbl>          <dbl>         <dbl>
## 1 FKBP5            0.425           0.435        0.00496        0.0338
## 2 SLC5~            0.713           0.549        0.00479        0.0333
## # ... with 8 more variables: m_midy <dbl>, m_wt <dbl>, sc_midy <dbl>,
## #   sc_wt <dbl>, m <dbl>, sc <dbl>, midy <dbl>, wt <dbl>
```

```r
lm_tissue_resids_fdr %>% 
  dplyr::filter(name == "FKBP5" | name == "SLC5A10")
```

```
## # A tibble: 2 x 9
## # Groups:   name [2]
##   name  m     sc    p_value mean_sc mean_m log2_difference mean_total
##   <chr> <lis> <lis>   <dbl>   <dbl>  <dbl>           <dbl>      <dbl>
## 1 FKBP5 <tib~ <tib~ 0.00733  -0.692  0.692           -1.38  -5.13e-15
## 2 SLC5~ <tib~ <tib~ 0.00846   0.785 -0.785            1.57  -9.08e-15
## # ... with 1 more variable: qvalue <dbl>
```

All but 2 proteins were found to be significant in the ANOVA are significant in
the t-test as well!

#### 6.3.2.1 Volcano plot


```r
lm_tissue_resids_fdr_volcano <- prepare_data_for_volcano_plot(lm_tissue_resids_fdr, group_neg = "Mesenterial",
                              group_pos = "Subcutaneous")

plot_fdr_volcano(lm_tissue_resids_fdr_volcano,
                 plot_title = "Volcano plot showing significant changes in protein expression\nbetween mesenterial and subcutaneous adipose tissue in pigs",
                 break_vector = c("Mesenterial", "Subcutaneous"),
                 color_vector = c("#7CAE00", "black", "#C77CFF"))
```

![](MIDY_adipose_tissue_742_removed_files/figure-html/unnamed-chunk-82-1.png)<!-- -->

#### 6.3.2.2 Proteins significantly more abundant in subcutaneous


```r
# Filter for proteins significanly upregulated in subcutaneous
lm_t_sig_sc <- lm_tissue_resids_fdr %>% 
  dplyr::filter(log2_difference > 0 & qvalue < 0.05)

# Get 20 most significant proteins
kable(lm_t_sig_sc %>%
        ungroup() %>%
        top_n(-20, qvalue) %>%
        dplyr::select(name, log2_difference, p_value, qvalue) %>%
        left_join(protein_names, by = "name") %>%
        arrange(qvalue) %>%
        mutate_if(is.numeric, funs(as.character(signif(., 3)))) %>% 
        rename(`Gene name` = name, `Log2 fold change` = log2_difference, `p-value` = p_value, `q-value` = qvalue, `Protein name` = protein_name),
      caption = "Top 20 most significant proteins more abundant in subcutaneous adipose tissue",
      align = "c"
) %>%
  kable_styling(bootstrap_options = c("striped", "hover"), full_width = F) %>%
  row_spec(0, background = "#C77CFF") %>% 
  row_spec(1:20, color = "black")
```

<table class="table table-striped table-hover" style="width: auto !important; margin-left: auto; margin-right: auto;">
<caption>Top 20 most significant proteins more abundant in subcutaneous adipose tissue</caption>
 <thead>
  <tr>
   <th style="text-align:center;background-color: #C77CFF;"> Gene name </th>
   <th style="text-align:center;background-color: #C77CFF;"> Log2 fold change </th>
   <th style="text-align:center;background-color: #C77CFF;"> p-value </th>
   <th style="text-align:center;background-color: #C77CFF;"> q-value </th>
   <th style="text-align:center;background-color: #C77CFF;"> Protein name </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:center;color: black;"> CBR2 </td>
   <td style="text-align:center;color: black;"> 5.08 </td>
   <td style="text-align:center;color: black;"> 1.26e-10 </td>
   <td style="text-align:center;color: black;"> 2.8e-07 </td>
   <td style="text-align:center;color: black;"> carbonyl reductase </td>
  </tr>
  <tr>
   <td style="text-align:center;color: black;"> OGN </td>
   <td style="text-align:center;color: black;"> 1.48 </td>
   <td style="text-align:center;color: black;"> 7.06e-10 </td>
   <td style="text-align:center;color: black;"> 5.53e-07 </td>
   <td style="text-align:center;color: black;"> mimecan isoform X1 </td>
  </tr>
  <tr>
   <td style="text-align:center;color: black;"> FERMT2 </td>
   <td style="text-align:center;color: black;"> 0.72 </td>
   <td style="text-align:center;color: black;"> 3.23e-08 </td>
   <td style="text-align:center;color: black;"> 1.27e-05 </td>
   <td style="text-align:center;color: black;"> fermitin family homolog 2 isoform X6 </td>
  </tr>
  <tr>
   <td style="text-align:center;color: black;"> CFD </td>
   <td style="text-align:center;color: black;"> 3.75 </td>
   <td style="text-align:center;color: black;"> 4.13e-08 </td>
   <td style="text-align:center;color: black;"> 1.39e-05 </td>
   <td style="text-align:center;color: black;"> complement factor D </td>
  </tr>
  <tr>
   <td style="text-align:center;color: black;"> COL3A1 </td>
   <td style="text-align:center;color: black;"> 3.95 </td>
   <td style="text-align:center;color: black;"> 6.43e-08 </td>
   <td style="text-align:center;color: black;"> 1.68e-05 </td>
   <td style="text-align:center;color: black;"> collagen alpha-1(III) chain precursor </td>
  </tr>
  <tr>
   <td style="text-align:center;color: black;"> SDR16C5 </td>
   <td style="text-align:center;color: black;"> 1.81 </td>
   <td style="text-align:center;color: black;"> 6.39e-08 </td>
   <td style="text-align:center;color: black;"> 1.68e-05 </td>
   <td style="text-align:center;color: black;"> epidermal retinol dehydrogenase 2 isoform X1 </td>
  </tr>
  <tr>
   <td style="text-align:center;color: black;"> DCN </td>
   <td style="text-align:center;color: black;"> 1.55 </td>
   <td style="text-align:center;color: black;"> 5.13e-07 </td>
   <td style="text-align:center;color: black;"> 8.6e-05 </td>
   <td style="text-align:center;color: black;"> decorin precursor </td>
  </tr>
  <tr>
   <td style="text-align:center;color: black;"> LUM </td>
   <td style="text-align:center;color: black;"> 1.55 </td>
   <td style="text-align:center;color: black;"> 6.41e-07 </td>
   <td style="text-align:center;color: black;"> 9.4e-05 </td>
   <td style="text-align:center;color: black;"> lumican precursor </td>
  </tr>
  <tr>
   <td style="text-align:center;color: black;"> MMP2 </td>
   <td style="text-align:center;color: black;"> 2.18 </td>
   <td style="text-align:center;color: black;"> 7.05e-07 </td>
   <td style="text-align:center;color: black;"> 9.4e-05 </td>
   <td style="text-align:center;color: black;"> 72 kDa type IV collagenase precursor </td>
  </tr>
  <tr>
   <td style="text-align:center;color: black;"> PCOLCE </td>
   <td style="text-align:center;color: black;"> 3.48 </td>
   <td style="text-align:center;color: black;"> 6.32e-07 </td>
   <td style="text-align:center;color: black;"> 9.4e-05 </td>
   <td style="text-align:center;color: black;"> procollagen C-endopeptidase enhancer 1 </td>
  </tr>
  <tr>
   <td style="text-align:center;color: black;"> PHGDH </td>
   <td style="text-align:center;color: black;"> 4.15 </td>
   <td style="text-align:center;color: black;"> 7.6e-07 </td>
   <td style="text-align:center;color: black;"> 9.4e-05 </td>
   <td style="text-align:center;color: black;"> D-3-phosphoglycerate dehydrogenase </td>
  </tr>
  <tr>
   <td style="text-align:center;color: black;"> SNX3 </td>
   <td style="text-align:center;color: black;"> 0.62 </td>
   <td style="text-align:center;color: black;"> 8.55e-07 </td>
   <td style="text-align:center;color: black;"> 1e-04 </td>
   <td style="text-align:center;color: black;"> sorting nexin-3 isoform X1 </td>
  </tr>
  <tr>
   <td style="text-align:center;color: black;"> SORBS1 </td>
   <td style="text-align:center;color: black;"> 0.972 </td>
   <td style="text-align:center;color: black;"> 1.33e-06 </td>
   <td style="text-align:center;color: black;"> 0.000149 </td>
   <td style="text-align:center;color: black;"> sorbin and SH3 domain-containing protein 1 isoform X23 </td>
  </tr>
  <tr>
   <td style="text-align:center;color: black;"> PLCD1 </td>
   <td style="text-align:center;color: black;"> 0.501 </td>
   <td style="text-align:center;color: black;"> 2.02e-06 </td>
   <td style="text-align:center;color: black;"> 0.000198 </td>
   <td style="text-align:center;color: black;"> 1-phosphatidylinositol 4,5-bisphosphate phosphodiesterase delta-1 </td>
  </tr>
  <tr>
   <td style="text-align:center;color: black;"> XP_020930866 </td>
   <td style="text-align:center;color: black;"> 0.789 </td>
   <td style="text-align:center;color: black;"> 2.73e-06 </td>
   <td style="text-align:center;color: black;"> 0.000229 </td>
   <td style="text-align:center;color: black;"> LOW QUALITY PROTEIN: tensin-1 </td>
  </tr>
  <tr>
   <td style="text-align:center;color: black;"> TUBB </td>
   <td style="text-align:center;color: black;"> 0.508 </td>
   <td style="text-align:center;color: black;"> 3.91e-06 </td>
   <td style="text-align:center;color: black;"> 0.000296 </td>
   <td style="text-align:center;color: black;"> tubulin beta chain </td>
  </tr>
  <tr>
   <td style="text-align:center;color: black;"> F13A1 </td>
   <td style="text-align:center;color: black;"> 1.14 </td>
   <td style="text-align:center;color: black;"> 4.98e-06 </td>
   <td style="text-align:center;color: black;"> 0.000365 </td>
   <td style="text-align:center;color: black;"> coagulation factor XIII A chain </td>
  </tr>
  <tr>
   <td style="text-align:center;color: black;"> UGDH </td>
   <td style="text-align:center;color: black;"> 1.03 </td>
   <td style="text-align:center;color: black;"> 5.57e-06 </td>
   <td style="text-align:center;color: black;"> 0.000396 </td>
   <td style="text-align:center;color: black;"> UDP-glucose 6-dehydrogenase </td>
  </tr>
  <tr>
   <td style="text-align:center;color: black;"> AGT </td>
   <td style="text-align:center;color: black;"> 0.601 </td>
   <td style="text-align:center;color: black;"> 6.09e-06 </td>
   <td style="text-align:center;color: black;"> 0.000421 </td>
   <td style="text-align:center;color: black;"> angiotensinogen </td>
  </tr>
  <tr>
   <td style="text-align:center;color: black;"> TMSB4X </td>
   <td style="text-align:center;color: black;"> 0.76 </td>
   <td style="text-align:center;color: black;"> 7.4e-06 </td>
   <td style="text-align:center;color: black;"> 0.000483 </td>
   <td style="text-align:center;color: black;"> thymosin beta-4 isoform X1 </td>
  </tr>
</tbody>
</table>

```r
# Get 20 proteins with highest l2fc
lm_t_sig_sc %>% 
  ungroup() %>% 
  top_n(20, log2_difference)
```

```
## # A tibble: 20 x 9
##    name  m     sc     p_value mean_sc mean_m log2_difference mean_total
##    <chr> <lis> <lis>    <dbl>   <dbl>  <dbl>           <dbl>      <dbl>
##  1 ABI3~ <tib~ <tib~ 2.50e- 4    1.54  -1.54            3.07  -7.11e-15
##  2 CBR2  <tib~ <tib~ 1.26e-10    2.54  -2.54            5.08  -3.55e-15
##  3 CFD   <tib~ <tib~ 4.13e- 8    1.87  -1.87            3.75   5.92e-16
##  4 COL1~ <tib~ <tib~ 8.26e- 5    1.59  -1.59            3.19  -8.09e-15
##  5 COL3~ <tib~ <tib~ 6.43e- 8    1.98  -1.98            3.95  -6.71e-15
##  6 DNAJ~ <tib~ <tib~ 8.11e- 4    1.35  -1.35            2.71  -5.72e-15
##  7 ECM1  <tib~ <tib~ 8.22e- 4    1.47  -1.47            2.93  -1.97e-15
##  8 FKBP9 <tib~ <tib~ 1.26e- 3    1.31  -1.31            2.61   4.14e-15
##  9 IKBIP <tib~ <tib~ 3.60e- 5    1.40  -1.40            2.80  -1.18e-15
## 10 MFAP2 <tib~ <tib~ 1.43e- 5    2.05  -2.05            4.11  -9.87e-16
## 11 MRC2  <tib~ <tib~ 3.13e- 4    1.68  -1.68            3.36  -3.36e-15
## 12 MYOC  <tib~ <tib~ 3.30e- 4    1.38  -1.38            2.76  -6.12e-15
## 13 NXN   <tib~ <tib~ 6.27e- 4    2.43  -2.43            4.86   9.87e-16
## 14 PCOL~ <tib~ <tib~ 6.32e- 7    1.74  -1.74            3.48  -5.92e-16
## 15 PHGDH <tib~ <tib~ 7.60e- 7    2.08  -2.08            4.15  -8.29e-15
## 16 SCP2~ <tib~ <tib~ 4.49e- 4    1.69  -1.69            3.37  -9.08e-15
## 17 VCAN  <tib~ <tib~ 2.73e- 3    1.47  -1.47            2.95  -6.12e-15
## 18 VPS2~ <tib~ <tib~ 2.02e- 4    1.39  -1.39            2.78  -1.13e-14
## 19 XP_0~ <tib~ <tib~ 5.52e- 4    1.46  -1.46            2.93  -5.92e-15
## 20 XP_0~ <tib~ <tib~ 5.96e- 5    1.35  -1.35            2.69  -5.92e-15
## # ... with 1 more variable: qvalue <dbl>
```



```r
lm_sig_m <- lm_tissue_resids_fdr %>% 
  dplyr::filter(log2_difference < 0 & qvalue < 0.05)
```

