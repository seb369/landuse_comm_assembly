Calculating βNTI
================
Samuel Barnett
06 December, 2019

## Introduction

In this notebook we are calculating the βNTI across samples to examine
soil bacterial community assembly in a finger lakes region. The data
generated here will be used in other notebooks for analysis.

This analysis uses the same (-2, 2) range for significance testing as in
Stegen et al. 2012. This means that the following conclusison can be
drawn from this data: βNTI \> 2: Community assembly driven by variable
selection βNTI \< -2: Community assembly driven by homgenizing selection
|βNTI| \< 2: Community assembly is stochastic

Please note that some of these cell have eval=FALSE in the rmarkdown
version so that they wont run when knitted. This was done because they
take a really long time and are better run in terminal.

### Initiate libraries

``` r
# Packages needed for analysis
library(dplyr)
library(tidyr)
library(tibble)
library(phyloseq)
library(picante)
library(parallel)

# Packages needed for plotting
library(ggplot2)
```

### Import data

``` r
# Import bulk soil phyloseq data
bulk.physeq = readRDS("/home/sam/data/fullCyc2_data/bulk_soil_physeq.RDS")

## Check how many reads you have in each of the samples. This will tell you if you need to re-do anything
# Get read counts and make a new dataframe with this data
read_count = data.frame("count" = colSums(otu_table(bulk.physeq))) %>%
  rownames_to_column(var="X.Sample") %>%
  inner_join(data.frame(sample_data(bulk.physeq)), by="X.Sample") %>%
  arrange(-count) %>%
  mutate(X.Sample=factor(X.Sample, levels=X.Sample))

# Now plot read count for each sample. The horizontal line represents a 2000 read threshold
ggplot(data=read_count, aes(x=X.Sample, y=log10(count), fill=ecosystem)) +
  geom_bar(stat="identity") +
  labs(x="Sample", y="Log10(Read count)") +
  geom_hline(yintercept=log10(10000)) +
  theme(text = element_text(size=16),
        axis.text.x = element_blank())
```

![](bNTI_calculation_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
# Everything seems to be at or above 10000 total reads

bulk.physeq
```

    ## phyloseq-class experiment-level object
    ## otu_table()   OTU Table:         [ 15112 taxa and 30 samples ]
    ## sample_data() Sample Data:       [ 30 samples by 30 sample variables ]
    ## tax_table()   Taxonomy Table:    [ 15112 taxa by 7 taxonomic ranks ]
    ## phy_tree()    Phylogenetic Tree: [ 15112 tips and 15111 internal nodes ]

Now we need to rarefy the data to normalize the sequencing depth. We
should also get a normalized dataset which gives relative abundance
rather than readcounts.

``` r
# Rarefy to an even depth
set.seed(72)  # setting seed for reproducibility
bulk.physeq.rare = rarefy_even_depth(bulk.physeq)

# Normalize read counts (this gives relative abundance)
bulk.physeq.norm = transform_sample_counts(bulk.physeq.rare, function(x) x/sum(x))
```

## Calculating βNTI

### The function

This will the the function used to calculate the βNTI for samples from a
phyloseq.

``` r
# Function for calculating the βMNTD for each random null community
bMNTD_null_func <- function(i, OTU.table, tree){
  tree$tip.label = sample(tree$tip.label)
  bMNTD_s = comdistnt(OTU.table, cophenetic(tree), abundance.weighted = TRUE)
  A <- attr(bMNTD_s, "Size")
  B <- if (is.null(attr(bMNTD_s, "Labels"))) sequence(A) else attr(bMNTD_s, "Labels")
  if (isTRUE(attr(bMNTD_s, "Diag"))) attr(bMNTD_s, "Diag") <- FALSE
  if (isTRUE(attr(bMNTD_s, "Upper"))) attr(bMNTD_s, "Upper") <- FALSE
  bMNTD_s.df = data.frame(Sample_1 = B[unlist(lapply(sequence(A)[-1], function(x) x:A))],
                          Sample_2 = rep(B[-length(B)], (length(B)-1):1),
                          bMNTD = as.vector(bMNTD_s),
                          rep=i)
  return(bMNTD_s.df)
}

# The main function for calculating βNTI
Phylo_turnover <- function(physeq, reps, nproc){
  # Extract OTU table
  OTU.table = t(otu_table(physeq))
  # Extract phylogenetic tree
  tree = phy_tree(physeq)

  # Get βMNTD between all communities
  bMNTD_o = comdistnt(OTU.table, cophenetic(tree), abundance.weighted = TRUE)
  A <- attr(bMNTD_o, "Size")
  B <- if (is.null(attr(bMNTD_o, "Labels"))) sequence(A) else attr(bMNTD_o, "Labels")
  if (isTRUE(attr(bMNTD_o, "Diag"))) attr(bMNTD_o, "Diag") <- FALSE
  if (isTRUE(attr(bMNTD_o, "Upper"))) attr(bMNTD_o, "Upper") <- FALSE
  bMNTD_o.df = data.frame(Sample_1 = B[unlist(lapply(sequence(A)[-1], function(x) x:A))],
                          Sample_2 = rep(B[-length(B)], (length(B)-1):1),
                          bMNTD = as.vector(bMNTD_o))
  
  # Get βMNTD for randomized null communities
  rep.list = seq(1, reps)
  bMNTD_s.df.list = mclapply(rep.list, bMNTD_null_func, OTU.table=OTU.table, tree=tree, mc.cores=nproc)
  
  # Combine all data together and calculate βNTI for each sample pair
  bMNTD_s.df <- do.call("rbind", bMNTD_s.df.list)
  bMNTD_s.means.df = bMNTD_s.df %>%
    group_by(Sample_1, Sample_2) %>%
    dplyr::summarize(mean_bMNTD = mean(bMNTD),
                     sd_bMNTD = sd(bMNTD))
  
  bMNTD_o.df = inner_join(bMNTD_o.df, bMNTD_s.means.df, by=c("Sample_1", "Sample_2")) %>%
    mutate(bNTI = (bMNTD - mean_bMNTD)/sd_bMNTD)

  return(bMNTD_o.df)
}
```

### Running function to get βNTI

Only run this cell once (eval=FALSE) because it take a really long time.
Also best to do this in terminal not in RStudio

``` r
# Get permutation beta null deviations
full.bNTI.df = Phylo_turnover(bulk.physeq.rare, 1000, 10)

# Save bNTI table
write.table(full.bNTI.df, "/home/sam/data/fullcyc2_data/Final_data/community_assembly/full_bNTI.txt")
```

What does this table look
like?

``` r
head(read.table("/home/sam/data/fullCyc2_data/Final_data/community_assembly/full_bNTI.txt"))
```

    ##      Sample_1    Sample_2      bMNTD mean_bMNTD     sd_bMNTD      bNTI
    ## 1 EL.A.151026 CC.A.151027 0.02016516 0.02375343 0.0006834651 -5.250116
    ## 2 CF.A.151027 CC.A.151027 0.01910554 0.02356232 0.0006401618 -6.961971
    ## 3 EL.M.151026 CC.A.151027 0.03881704 0.04066696 0.0013836678 -1.336973
    ## 4 MF.F.151026 CC.A.151027 0.03210584 0.03561841 0.0015856144 -2.215275
    ## 5 MP.F.151026 CC.A.151027 0.01923891 0.02469787 0.0007593102 -7.189362
    ## 6 MW.M.151027 CC.A.151027 0.02239621 0.02642288 0.0007702880 -5.227487

# Session info

``` r
sessionInfo()
```

    ## R version 3.4.4 (2018-03-15)
    ## Platform: x86_64-pc-linux-gnu (64-bit)
    ## Running under: Ubuntu 18.04.3 LTS
    ## 
    ## Matrix products: default
    ## BLAS: /usr/lib/x86_64-linux-gnu/blas/libblas.so.3.7.1
    ## LAPACK: /usr/lib/x86_64-linux-gnu/lapack/liblapack.so.3.7.1
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## attached base packages:
    ## [1] parallel  stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] ggplot2_3.2.1   picante_1.8     nlme_3.1-142    vegan_2.5-6    
    ##  [5] lattice_0.20-38 permute_0.9-5   ape_5.3         phyloseq_1.22.3
    ##  [9] tibble_2.1.3    tidyr_1.0.0     dplyr_0.8.3    
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] tidyselect_0.2.5    xfun_0.10           reshape2_1.4.3     
    ##  [4] purrr_0.3.2         splines_3.4.4       rhdf5_2.22.0       
    ##  [7] colorspace_1.4-1    vctrs_0.2.0         htmltools_0.4.0    
    ## [10] stats4_3.4.4        mgcv_1.8-31         yaml_2.2.0         
    ## [13] survival_3.1-7      rlang_0.4.0         pillar_1.4.2       
    ## [16] withr_2.1.2         glue_1.3.1          BiocGenerics_0.24.0
    ## [19] foreach_1.4.7       lifecycle_0.1.0     plyr_1.8.4         
    ## [22] stringr_1.4.0       zlibbioc_1.24.0     Biostrings_2.46.0  
    ## [25] munsell_0.5.0       gtable_0.3.0        codetools_0.2-16   
    ## [28] evaluate_0.14       labeling_0.3        Biobase_2.38.0     
    ## [31] knitr_1.25          IRanges_2.12.0      biomformat_1.6.0   
    ## [34] Rcpp_1.0.2          backports_1.1.5     scales_1.0.0       
    ## [37] S4Vectors_0.16.0    jsonlite_1.6        XVector_0.18.0     
    ## [40] digest_0.6.21       stringi_1.4.3       ade4_1.7-13        
    ## [43] grid_3.4.4          tools_3.4.4         magrittr_1.5       
    ## [46] lazyeval_0.2.2      cluster_2.1.0       crayon_1.3.4       
    ## [49] pkgconfig_2.0.3     zeallot_0.1.0       MASS_7.3-51.4      
    ## [52] Matrix_1.2-17       data.table_1.12.4   assertthat_0.2.1   
    ## [55] rmarkdown_1.16      iterators_1.0.12    R6_2.4.0           
    ## [58] multtest_2.34.0     igraph_1.2.4.1      compiler_3.4.4
