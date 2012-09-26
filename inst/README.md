# GOglm
=============================================

The R package **GOglm** implements the `GOglm` approach discussed in Mi *et al.* [1]. It includes a summarized RNA-Seq data example (the prostate cancer dataset [2]) for methodological illustrations. This is a beta version under development. In this README file, we briefly outline the setup of the logistic regression model for length bias corrections, followed by discussions on access to the RNA-Seq datasets, Gene Ontology (GO) annotations, variable transformations, potential computational issues, choice of DE testing procedures and R/Bioconductor software information.

## GOglm and logistic regression

### Logistic regression and 2-by-2 contingency table

In the generalized linear model (GLM) framework, we used continuous measures of DE as predictors and interpret the results in terms of odds. Sartor *et al.* [3] proposed the LRpath method in the microarray context with the following setup

$$
logit\left[\pi(x)\right]=\beta_{0}+\beta_{1}x
$$

where `x` is the significance statistic defined as -log($p$-value). In the traditional 2-by-2 contingency table framework, a $p$-value cut-off for declaring DE genes is pre-specified so that each gene is associated with a indicator of 0 or 1. All genes are then cross-classified into a 2-by-2 table ready for Fisher's exact test or any other contingency-table-based approaches. If we use this binary indicator as the predictor in the equation above, the $p$-values for testing $\beta_{1}=0$  are roughly equivalent to those obtained using contingency tables, though with minor scale differences.

### GOglm using logistic regression for length bias adjustment

Equation (2) in paper:

$$
logit\left[\pi(x)\right]=\beta_{0}+\beta_{1}x+\beta_{2}L
$$

is what we proposed to use for correcting length bias in GO enrichment analysis. The fundamental cause of length bias is that the transcript length becomes a confounding factor when it correlates with both the GO membership and the DE test significance. When we include transcript length `L` as a covariate, the coefficient $\beta_{1}$ now captures the correlation between the log odds of being in the specified GO category and the DE test significance -- conditional on gene lengths. A significant result from the hypothesis test $H0: \beta_{1} = 0$ indicates that the GO membership is correlated with the DE test significance even after adjusting for length bias.

## Access to the dataset

The first prostate cancer dataset first analyzed in [4] can be obtained from the `goseq` Bioconductor package. It is also included in our `GOglm` package as a data example. The following R codes will load `goseq` into R session and display the first five ENSEMBL genes with 1 indicating differential expression (DE) and 0 otherwise.


```r
library(goseq)
data(genes)
head(genes, 5)
```

```
## ENSG00000230758 ENSG00000182463 ENSG00000124208 ENSG00000230753 
##               0               0               0               0 
## ENSG00000224628 
##               0
```


Testing DE genes was implemented using the `edgeR` package. We adopted the common dispersion parameter approach to fit the negative binomial (NB) model using the quantile-adjusted conditional maximum likelihood (qCML) method, but the tagwise dispersion approach should also work. The following table shows part of the original raw counts table.


```
##                 lane1 lane2 lane3 lane4 lane5 lane6 lane8
## ENSG00000089057   174    98   164   208   252   173    96
## ENSG00000125520     6     1    17    17    23    11     5
## ENSG00000207427     0     0     0     0     0     0     0
## ENSG00000101152    41    33    32    46    35    44    13
## ENSG00000089199     0     0     0     0     0     2     0
```


This `n`-by-`r` matrix of RNA-seq read counts is also required for the NBP exact test `nbp.test()` implemented in `NBPSeq`, and for the nonparametric modeling of the variance `nbinomTest()` implemented in `DESeq`.

For ease of comparison with published results in our paper, in analyzing the prostate cancer data, we used `edgeR` with a common dispersion estimate to obtain DE test $p$-values. For the Arabidopsis dataset (to be included in package), we used `NBPSeq` to obtain DE test $p$-values. EdgeR and NBPSeq are both based on NB models for RNA-Seq read frequencies. The NB model captures potential extra-Poisson variation in RNA-Seq read frequencies between independent biological samples using a dispersion parameter. Other methods based on NB model include the tagwise or trend options in `edgeR`, or the `DESeq` approach. All of these methods use the same exact NB test for assessing DE, but differ in how they estimate the dispersion parameter as a function of the mean frequency.

The table below shows the partial result of the NB exact test using `edgeR`. The genes are ordered by DE testing $p$-values, and significance statistics will be obtained by proper transformations of these (un-adjusted) $p$-values.


```
##                  logFC logCPM    PValue       FDR
## ENSG00000127954 12.373  6.663 2.252e-80 1.115e-75
## ENSG00000151503  5.403  8.495 1.532e-65 3.793e-61
## ENSG00000096060  4.888  9.444 6.908e-60 1.140e-55
## ENSG00000091879  5.669  6.259 1.094e-54 1.354e-50
## ENSG00000132437 -5.931  7.945 2.638e-52 2.612e-48
```


In `GOglm` we load the data and obtain a data frame ready for use by the `goglm` function. Essentially we need two columns: `$p$-value` and `length`. The subsetted data frame `psc.subset` will be passed to the first argument of `goglm`:


```r
library(GOglm)
data(ProsCancer)
rownames(ProsCancer) <- as.character(ProsCancer$X)
head(ProsCancer)
```

```
##                               X  logFC logCPM    PValue       FDR length
## ENSG00000127954 ENSG00000127954 12.373  6.663 2.252e-80 1.115e-75   3923
## ENSG00000151503 ENSG00000151503  5.403  8.495 1.532e-65 3.793e-61   5640
## ENSG00000096060 ENSG00000096060  4.888  9.444 6.908e-60 1.140e-55   4106
## ENSG00000091879 ENSG00000091879  5.669  6.259 1.094e-54 1.354e-50   2536
## ENSG00000132437 ENSG00000132437 -5.931  7.945 2.638e-52 2.612e-48   1968
## ENSG00000166451 ENSG00000166451  4.560  8.453 6.313e-52 5.209e-48   1567
##                 log.gl sig.stat de.ind
## ENSG00000127954  8.275    5.217      1
## ENSG00000151503  8.638    5.012      1
## ENSG00000096060  8.320    4.922      1
## ENSG00000091879  7.838    4.830      1
## ENSG00000132437  7.585    4.786      1
## ENSG00000166451  7.357    4.778      1
```

```r
## subset to two columns required by goglm:
psc.subset <- ProsCancer[, c("PValue", "length")]
head(psc.subset)
```

```
##                    PValue length
## ENSG00000127954 2.252e-80   3923
## ENSG00000151503 1.532e-65   5640
## ENSG00000096060 6.908e-60   4106
## ENSG00000091879 1.094e-54   2536
## ENSG00000132437 2.638e-52   1968
## ENSG00000166451 6.313e-52   1567
```

```r
## result <- goglm(data = psc.subset, trans.p = 'd.log', genome = 'hg18',
## id = 'ensGene')
```


## GO annotations

### Access to the GO database
As shown above, each row corresponds to one gene with "accessible" gene names, such as `ENSG00000127954` in the prostate cancer dataset. What makes "accessible" so important is that we need these identifiers to obtain corresponding GO annotations. We will then know which genes are annotated to any particular category so that the response variable (indicator of whether a gene is annotated to this category) can be constructed. 

The `getgo` function in `goseq` makes use of Bioconductor organism packages to obtain mappings between gene and GO IDs. These packages all have similar names as `org.<Genome>.<GeneID>.db` so that it's convenient to find the database for a particular organism. We loaded `org.Hs.eg.db` for the prostate cancer dataset. This database can be installed using

```r
source("http://www.bioconductor.org/biocLite.R")
biocLite(c("org.Hs.eg.db"))
```

Besides gene names, we still need to pass the "genome"" and  "gene identifier" used by the first argument `genes` to `getgo` and `getlength` functions in `goseq`. For organisms available in the UCSC genome browser, it is relatively easy to find such information using `supportedGenomes` and `supportedGeneIDs`, respectively. In this case, we specify `genome=hg18` and `id=ensGene`. So far, all information required is collected for the prostate cancer dataset. If the organism under study does not have a Bioconductor annotation database for GO mapping, then the user has to provide GO annotations in other ways.

The Bioconductor package `GO.db` is the only set of annotation maps describing the entire GO. It will be used after we identified top-ranked categories and would like to get more details about categories to see if they make biological sense.

### GO availability and gene filtering

In testing our GLM approach we find that the `goseq` function simplifies a key aspect in constructing the contingency table. The total number of genes $N$ is defined as `nrow(pwf)`. This works well if all genes under study have GO annotations available. However this is not true when we obtain a complete list of

```r
getgo(rownames(dataset),"hg18","ensGene")
```

Those `$<NA>` values for genes are associated with `NULL` GO information. No matter which category is tested for enrichment, those genes without GO annotations will automatically be considered as "not in this category" but in fact we don't know which categories it belongs to. Among 22743 genes, only 12641 genes (55.58%) have annotations, so we believe it is more sensible to discard genes without annotations. This will drastically change the final list obtained by GOseq Wallenius.

## Variable transformations, computational issues and choice of DE testing procedures

Predictors in the logistic regression model, namely significance statistics and gene lengths, should be properly transformed in order to avoid computational issues. As a result of many DE testing $p$-values exactly equal to 1, log(-log($p$-value)) will produce `NaN` values not usable in R. Therefore, we use log(1-log($p$-value)) so that those genes' information can be retained.

Non-convergence of Newton-Raphson algorithm is more likely to occur when samples are small. This can be reflected in two ways from our data example: either more required iterations are required to achieve convergence or over-dispersion parameter not close to 1. Non-convergence is mostly suffered by specialized categories with very few annotated genes (say, less than 5).

We didn't explore how different methods for identifying DE genes will influence the downstream enrichment analyses. Some researchers believe that methods for DE testing are in some sense irrelevant to subsequent analysis (see `GOstats` package vignette, `GOvis`), but adopting the currently widely used approaches such as `edgeR`, `DESeq` and `NBPSeq` will make the results more trustworthy. We note in paper that these NB-based approaches are superior when biological replicates present, but potential differences in the final list may result from gene filtering criteria on the population genes, such as deciding fold-changes for (non-)expressed genes and/or discarding genes with unavailable $p$-values.

Bioconductor annotation databases are updated regularly as the state of biological knowledge changes, so that results might be slightly different as releases of packages change. Our analyses in the paper were based on R version 2.15.1 (2012-06-22) and Bioconductor release 2.9 (November 1, 2011).

## References

[1] Mi G, Di Y, Emerson S, Cumbie JS and Chang JH (2012) "Length bias correction in Gene Ontology enrichment analysis using logistic regression", PLOS ONE, in press.

[2] Li H, Lovci M, Kwon Y, Rosenfeld M, Fu X, et al. (2008) "Determination of tag density required for digital transcriptome analysis: application to an androgen-sensitive prostate cancer model", Proc Natl Acad Sci U S A 105: 20179-20184.

[3] Sartor M, Leikauf G, Medvedovic M (2009) "LRpath: a logistic regression approach for identifying enriched biological groups in gene expression data", Bioinformatics 25: 211-217.

[4] Young M, Wakefield M, Smyth G, Oshlack A (2010) "Gene ontology analysis for RNA-seq: accounting for selection bias", Genome Biol 11: R14.
