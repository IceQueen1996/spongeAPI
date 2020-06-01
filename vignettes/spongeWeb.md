---
title: "Sparse Partial correlation ON Gene Expression with SPONGE - A Webresource for ceRNA Interaction Networks"
author: "Markus List, Pachl Elisabeth, Markus Hoffmann"
date: "2020-06-01"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{spongeWeb}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---




```r
library(spongeWeb)
```

## Purpose

<div style="float:right;"><img src="https://raw.githubusercontent.com/biomedbigdata/SPONGE-web-R/master/vignettes/spongeLogo.png" alt="SPONGE logo" style="width: 200px; display: block; margin-left: auto; margin-right: auto;"/></div>

With SPONGE being an outstanding approach regarding calculation speed and accuracy, the goal, making the data available in an easy way for as many researchers as possible, is the next logical step. 
Furthermore, the data should become visualized in an interactive network, for uncomplicated research within a small part of interest of the networks. Available ceRNA interaction networks are based on paired gene and miRNA expression data taken from "The Cancer Genome Atlas" (TCGA). 

## Introduction

MicroRNAs (miRNAs) are important non-coding, post-transcriptional regulators that are involved in many biological processes and human diseases. miRNAs regulate their target mRNAs, so called competing endogenous RNAs (ceRNAs), by either degrading them or by preventing their translation. ceRNAs share similar miRNA recognition elements (MREs), sequences with a specific pattern where the belonging miRNA binds to.

miRNAs act as rheostats that regulate gene expression and maintain the functional balance of various gene networks. Furthermore, the miRNA-ceRNA-interactions follow a many-to-many relationship where one miRNA can affect multiple ceRNA targets and one ceRNA can contain multiple MREs for various miRNAs, leading to complex cross-talk. Since failures in these systems may lead to cancer, it is crucial to determine the networks of interactions and to assess their structure.

[SPONGE](https://academic.oup.com/bioinformatics/article/35/14/i596/5529172) is a method for the fast construction of a ceRNA network using ’multiple sensitivity correlation’. SPONGE was applied to paired miRNA and gene expression data from “The Cancer Genome Atlas” (TCGA) for studying global effects of miRNA-mediated cross-talk. The outcome highlights already established and novel protein-coding and non-coding ceRNAs which could serve as biomarkers in cancer. Further information about the SPONGE R package are available under https://bioconductor.org/packages/release/bioc/html/SPONGE.html.

With SPONGE being an outstanding approach regarding calculation speed and accuracy, the goal was to make its results easily accessible to researchers studying ceRNAs in cancer for further analysis. Therefore, an application programming interface (API) was developed, enabling other developers to query a database consisting of the SPONGE results on the TCGA dataset for their ceRNAs, miRNAs or cancer-types of interest. Containing a large amount of statistically related ceRNAs and their associated miRNAs from over 20 cancer types, the database allows to detect the importance of a single gene on a large scale, thus assessing its relevancy in a cancer background.

Additionally, an interactive web interface was set up to provide the possibility to browse and to search the database via a graphical user interface. The website also facilitates processing the data returned from the database and visualizing the ceRNA interactions as a network. The website is available under https://exbio.wzw.tum.de/sponge/home.
By help of these tools, third party developers like data scientists and biomedical researchers become able to carry out in depth cancer analyses and detect correlations between different cancer-causing factors on a new level while benefiting from an easy to use interface, which may lead to an uncomplicated and better understanding of cancer.

A complete definition of all API endpoints can be found under https://exbio.wzw.tum.de/sponge-api/ui.


## General Workflow

We have built a webresource to query SPONGE results easily via an API which requests data from the database. The API can be accessed through an graphical user interface of the API (Flask). If the data needs to be loaded inside a programming environment a R package and a python package is available for easy implementation. For visualisation of the data and for medical researches we provide a web application where graphs, networks and data are shown. The general build of the project is shown below.

<img src="https://raw.githubusercontent.com/biomedbigdata/SPONGE-web-R/master/vignettes/proj_structure.PNG" alt="Structure SPONGE-Web project" style="width: 600px; display: block; margin-left: auto; margin-right: auto;"/>

## How to start requests?
To start with further analysis with SPONGE data, it is important to get an overview about the available disease_types and the number of ceRNA interactions. This can be retrieved with:


```r
get_overallCounts()
```

To get an overview plot very easy, consider:

```r
library(ggplot2)
library(reshape2)
counts <- get_overallCounts()[c("count_interactions_sign", "count_shared_miRNAs", "disease_name")]
counts.m <- melt(counts, id.vars='disease_name')
ggplot(counts.m, aes(disease_name, value)) +   
     geom_bar(aes(fill = variable), position = "dodge", stat="identity") +
     theme(axis.text.x = element_text(angle=90)) +
     labs(y= "counts", x = "disease name")
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-1.png)

To get just an overview about the available datasets without any numbers use:

```{ r eval=FALSE}
get_datasetInformation()
```

To retrieve all used parameters of the SPONGE method to re-create published results for the cancer type of interest, use the following function:


```r
get_runInformation(disease_name = "kidney clear cell carcinoma")
```

Another way to get an overview of the results is to search for a specific gene and get an idea in which ceRNA interaction network the gene of interest contributes most to.


```r
get_geneCount(gene_symbol = c("HOXA1"))
```

\tiny

| count_all| count_sign|gene.ensg_number |gene.gene_symbol |run.dataset.data_origin | run.dataset.dataset_ID|run.dataset.disease_name              | run.run_ID|
|---------:|----------:|:----------------|:----------------|:-----------------------|----------------------:|:-------------------------------------|----------:|
|      3919|         58|ENSG00000105991  |HOXA1            |TCGA                    |                      1|head & neck squamous cell carcinoma   |          7|
|      5708|        108|ENSG00000105991  |HOXA1            |TCGA                    |                      2|breast invasive carcinoma             |          3|
|      3970|          8|ENSG00000105991  |HOXA1            |TCGA                    |                      3|kidney clear cell carcinoma           |          8|
|      2856|        162|ENSG00000105991  |HOXA1            |TCGA                    |                      4|sarcoma                               |         17|
|       827|          0|ENSG00000105991  |HOXA1            |TCGA                    |                      5|ovarian serous cystadenocarcinoma     |         13|
|      5262|        737|ENSG00000105991  |HOXA1            |TCGA                    |                      6|thymoma                               |         20|
|      4302|         15|ENSG00000105991  |HOXA1            |TCGA                    |                      7|brain lower grade glioma              |          2|
|      3242|          0|ENSG00000105991  |HOXA1            |TCGA                    |                      8|pheochromocytoma & paraganglioma      |         15|
|      5754|          4|ENSG00000105991  |HOXA1            |TCGA                    |                      9|liver hepatocellular carcinoma        |         10|
|      5093|        419|ENSG00000105991  |HOXA1            |TCGA                    |                     10|pancreatic adenocarcinoma             |         14|
|      3280|        126|ENSG00000105991  |HOXA1            |TCGA                    |                     12|colon adenocarcinoma                  |          5|
|      6201|        213|ENSG00000105991  |HOXA1            |TCGA                    |                     13|kidney papillary cell carcinoma       |          9|
|      4420|          0|ENSG00000105991  |HOXA1            |TCGA                    |                     14|lung adenocarcinoma                   |         11|
|      2390|          0|ENSG00000105991  |HOXA1            |TCGA                    |                     15|uterine corpus endometrioid carcinoma |         22|
|      3523|        166|ENSG00000105991  |HOXA1            |TCGA                    |                     16|thyroid carcinoma                     |         21|
|      1016|          0|ENSG00000105991  |HOXA1            |TCGA                    |                     17|bladder urothelial carcinoma          |          1|
|      8269|         56|ENSG00000105991  |HOXA1            |TCGA                    |                     18|prostate adenocarcinoma               |         16|
|      3247|        120|ENSG00000105991  |HOXA1            |TCGA                    |                     19|testicular germ cell tumor            |         19|
|      4906|          0|ENSG00000105991  |HOXA1            |TCGA                    |                     20|stomach adenocarcinoma                |         18|
|       811|          0|ENSG00000105991  |HOXA1            |TCGA                    |                     21|lung squamous cell carcinoma          |         12|
|      2169|          0|ENSG00000105991  |HOXA1            |TCGA                    |                     22|esophageal carcinoma                  |          6|
|      2739|         11|ENSG00000105991  |HOXA1            |TCGA                    |                     23|cervical & endocervical cancer        |          4|
|      8251|        522|ENSG00000105991  |HOXA1            |UCSC Xena               |                     11|pancancer                             |         23|
\normalsize


```r
library(ggplot2)
library(reshape2)
geneCounts <- get_geneCount(gene_symbol = c("HOXA1"))[c("count_all", "count_sign", "run.dataset.disease_name")]
geneCounts.m <- melt(geneCounts, id.vars='run.dataset.disease_name')
ggplot(geneCounts.m, aes(run.dataset.disease_name, value)) +   
     geom_bar(aes(fill = variable), position = "dodge", stat="identity") +
     theme(axis.text.x = element_text(angle=90)) +
     labs(y= "counts", x = "disease name")
```

![plot of chunk unnamed-chunk-7](figure/unnamed-chunk-7-1.png)

## How to find a subnetwork?
To find a sub network of nodes of interest use the functions: Get all ceRNA interactions by given identifications (ensg_number, gene_symbol or gene_type), specific cancer type or different filter possibilities according different statistical values (e.g. FDR adjusted p-value).

```r
# Retrieve all possible ceRNAs for gene, identified by ensg_number,
# and threshold for pValue and mscor.
get_all_ceRNAInteractions(disease_name = "pancancer", ensg_number=c("ENSG00000259090","ENSG00000217289"),
                          pValue=0.05, pValueDirection="<",
                          limit=15, information=FALSE)
```

Get all ceRNAs in a disease of interest (search not for a specific ceRNA, but search for all ceRNAs satisfying filter functions).


```r
get_ceRNA(disease_name = "kidney clear cell carcinoma",
          gene_type = "lincRNA", minBetweenness = 0.8)
```

Get all interactions between the given identifiers (ensg_number or gene_symbol).


```r
get_specific_ceRNAInteractions(disease_name = "pancancer",
                                gene_symbol = c("PTENP1","VCAN","FN1")
```

## How to find sponged miRNA?
Find sponged miRNAs (the reason for a edge between two ceRNAs) with

```r
get_sponged_miRNA(disease_name="kidney", gene_symbol = c("TCF7L1", "SEMA4B"), between=True)
```

or find a miRNA induced ceRNA interaction with


```r
# Retrieve all possible ceRNA interactions where miRNA(s) of interest contribute to
get_specific_miRNAInteraction(disease_name = "kidney clear cell carcinoma",
                              mimat_number = c("MIMAT0000076", "MIMAT0000261"),
                              limit = 15, information = FALSE)
```

## Further Information and Analysis
The database also contains information about the raw expression values and survival analyis data, which can be used to for Kaplan-Meyer-Plots (KMPs) for example. These information can be adressed with package functions.
To retrieve expression data use 


```r
# Retrieve gene expression values for specific genes by ensg_numbers
get_geneExprValues(disease_name = "kidney clear cell carcinoma",
                   ensg_number = c("ENSG00000259090","ENSG00000217289"))
# Retrieve gene expression values for specific miRNAs by mimat_numbers
get_mirnaExprValues(disease_name = "kidney clear cell carcinoma",
                    mimat_number = c("MIMAT0000076", "MIMAT0000261"))
```

To get survival analysis data use the function 


```r
get_survAna_rates(disease_name="kidney clear cell carcinoma",
                 ensg_number=c("ENSG00000259090", "ENSG00000217289"),
                 sample_ID = c("TCGA-BP-4968","TCGA-B8-A54F"))
```

It returns a data_frame with gene and patient/sample information and the "group information" encoded by column "overexpressed". Information about expression value of the gene (FALSE = underexpression, gene expression <= mean gene expression over all samples, TRUE = overexpression, gene expression >= mean gene expression over all samples)

For further patient/sample information:

```r
get_survAna_sampleInformation(disease_name = "kidney clear cell carcinoma",
                              sample_ID = c("TCGA-BP-4968","TCGA-B8-A54F"))
```

Further analysis and more complex information like associated cancer hallmarks, GO terms or wikipathway keys about the genes contributing to a network can be received by using this three functions:


```r
get_geneOntology(gene_symbol=c("PTEN","TIGAR"))
get_hallmark(gene_symbol=c("PTEN"))
get_WikiPathwayKey(gene_symbol=c("PTEN","GJA1"))
```

## Session Info

```r
sessionInfo()
#> R version 4.0.0 (2020-04-24)
#> Platform: x86_64-w64-mingw32/x64 (64-bit)
#> Running under: Windows 10 x64 (build 18363)
#> 
#> Matrix products: default
#> 
#> locale:
#> [1] LC_COLLATE=C                    LC_CTYPE=German_Germany.1252   
#> [3] LC_MONETARY=German_Germany.1252 LC_NUMERIC=C                   
#> [5] LC_TIME=German_Germany.1252    
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] reshape2_1.4.4  ggplot2_3.3.1   spongeWeb_0.1.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] Rcpp_1.0.4.6      highr_0.8         plyr_1.8.6        compiler_4.0.0   
#>  [5] pillar_1.4.4      bitops_1.0-6      prettyunits_1.1.1 tools_4.0.0      
#>  [9] digest_0.6.25     testthat_2.3.2    pkgbuild_1.0.8    pkgload_1.0.2    
#> [13] evaluate_0.14     gtable_0.3.0      jsonlite_1.6.1    tibble_3.0.1     
#> [17] lifecycle_0.2.0   pkgconfig_2.0.3   rlang_0.4.6       cli_2.0.2        
#> [21] rstudioapi_0.11   curl_4.3          xfun_0.14         httr_1.4.1       
#> [25] withr_2.2.0       stringr_1.4.0     dplyr_0.8.5       roxygen2_7.1.0   
#> [29] xml2_1.3.2        knitr_1.28        desc_1.2.0        vctrs_0.3.0      
#> [33] grid_4.0.0        tidyselect_1.1.0  rprojroot_1.3-2   glue_1.4.1       
#> [37] R6_2.4.1          processx_3.4.2    fansi_0.4.1       farver_2.0.3     
#> [41] callr_3.4.3       purrr_0.3.4       magrittr_1.5      scales_1.1.1     
#> [45] backports_1.1.6   ps_1.3.3          ellipsis_0.3.1    assertthat_0.2.1 
#> [49] colorspace_1.4-1  labeling_0.3      stringi_1.4.6     munsell_0.5.0    
#> [53] RCurl_1.98-1.2    crayon_1.3.4
```

## Citation
If you use any results from spongeWeb, please cite as follow:

```r
Pachl E., Hoffmann M., Hartung M., Stiegler V., Baumbach J., Schulz M., List M., SPONGE-web: A pan-cancer resource for competing endogenous RNA interactions (manuscript in preparation)
```
