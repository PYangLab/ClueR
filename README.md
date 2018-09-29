### Cluster Evaluation R package (ClueR) for detecting key kinases or pathways from time-series phosphoproteomics or gene expression (microarray, RNA-seq, or proteomics) datasets

#### Description
CLUster Evaluation (CLUE), a computational method, and its implementation in R (ClueR) is for detecting key kinases or pathways from a given time-series phosphoproteomics or gene expression dataset clustered by cmeans or kmeans algorithms. It firstly identifies the optimal number of clusters in the time-servies dataset; Then, it partition the dataset based on the optimal number of clusters determined in the first step; It finally detects kinases or pathways enriched in each cluster from optimally partitioned dataset.

The above three steps rely extensively on Fisher's exact test, Fisher's combined statistics, cluster regularisations, and they are performed against a user-specified reference annotation database such phosphoSitePlus in the case of phosphoproteomics data or KEGG in the case of gene expression data. There is a large selection of built-in annotation databases for both phosphoproteomics data and gene expression data but users can supply their own annotation database.

CLUE was initially designed for analysing time-course phosphoproteomics dataset using kinase-substrate annotation as reference (e.g. PhosphoSitePlus). It is now extended to identify key pathways from time-series microarray, RNA-seq or proteomics datasets by searching and testing against gene set annotation databases such as KEGG, GO, or Reactome etc.

Previously published phosphoproteomics dataset and gene expression dataset are included in the package to demonstrate how to use CLUE package.

#### Reference
Yang P, Zheng X, Jayaswal V, Hu G, Yang JYH, Jothi R (2015) Knowledge-Based Analysis for Detecting Key Signaling Events from Time-Series Phosphoproteomics Data. PLoS Comput Biol 11(8): e1004403. [fulltext](http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004403)

#### Download and install
The release version can be downloaded from CRAN [link](https://CRAN.R-project.org/package=ClueR);

whereas the latest release of development version can be downloaded from [here](https://github.com/PengyiYang/ClueR/releases)

1. Install the release version from CRAN with `install.packages("ClueR")`

2. You can also install the latest development version from github with:
```r
devtools::install_github("PengyiYang/ClueR")
```
Make sure that you have [Rtools](https://cran.r-project.org/bin/windows/Rtools/) install in your system for building the package from the source.

#### Examples
(1) This example demonstrate how CLUE can be applied to discover the optimal number of clusters from a simulated data.

``` r
## install the latest release of ClueR package

# load the package into R session
library(ClueR) 

# simulate a time-series data with 6 clusters and each cluster with a size of 500 entries
simuData <- temporalSimu(seed=1, groupSize=500, sdd=1, numGroups=6)

# create an artificial annotation database. Generate 100 kinase-substrate groups each
# comprising 50 substrates assigned to a kinase.
# among them, create 5 groups each contains phosphorylation sites defined to have the
# same temporal profile.
kinaseAnno <- list()
groupSize <- 500
for (i in 1:5) {
 kinaseAnno[[i]] <- paste("p", (groupSize*(i-1)+1):(groupSize*(i-1)+50), sep="_")
}

for (i in 6:100) {
 set.seed(i)
 kinaseAnno[[i]] <- paste("p", sample.int(nrow(simuData), size = 50), sep="_")
}
names(kinaseAnno) <- paste("KS", 1:100, sep="_")

# run CLUE with a repeat of 5 times and a range from 2 to 20
set.seed(2)
clueObj <- runClue(Tc=simuData, annotation=kinaseAnno, rep=5, kRange=20)

# visualize the evaluation outcome
xl <- "Number of clusters"
yl <- "Enrichment score"
boxplot(clueObj$evlMat, col=rainbow(ncol(clueObj$evlMat)), las=2, xlab=xl, ylab=yl, main="CLUE")
abline(v=(clueObj$maxK-1), col=rgb(1,0,0,.3))

# generate optimal clustering results using the optimal k determined by CLUE
best <- clustOptimal(clueObj, rep=5, mfrow=c(2, 3))

# list enriched clusters
best$enrichList

# obtain the optimal clustering object
best$clustObj
```

(2) This example shows the application of CLUE on a hES phosphoproteomics data set (Rigbolt et al. Sci Signal. 4(164):rs3, 2011) and uses kinase-substrate annotation compiled from [PhosphoSitePlus](http://www.phosphosite.org).

``` r
## install the latest release of ClueR package

# load the package into R session
library(ClueR) 

# load the human ES phosphoprotoemics data 
data(hES) 

# load the PhosphoSitePlus annotations (Hornbeck et al. Nucleic Acids Res. 40:D261-70, 2012). Note that one can instead use PhosphoELM database by typing "data(PhosphoELM)"
data(PhosphoSite)

# run CLUE with a repeat of 5 times and a range from 2 to 20
set.seed(2)
clueObj <- runClue(Tc=hES, annotation=PhosphoSite.human, rep=5, kRange=20)

# visualize the evaluation outcome
xl <- "Number of clusters"
yl <- "Enrichment score"
boxplot(clueObj$evlMat, col=rainbow(ncol(clueObj$evlMat)), las=2, xlab=xl, ylab=yl, main="CLUE")
abline(v=(clueObj$maxK-1), col=rgb(1,0,0,.3))

# generate the optimal clustering results
best <- clustOptimal(clueObj, rep=5, mfrow=c(3, 4))

# list enriched clusters
best$enrichList

# obtain the optimal clustering object
best$clustObj
```

(3) This example shows the application of CLUE to a gene expression dataset, discover optimal number of clusters, clustering data accordingly, and identify key pathway involved in each cluster.

``` r
# load mouse adipocyte gene expression data (Ma et al. Molecular and Cellular Biology. 2014, 34(19):3607-17)
data(adipocyte)

# load the KEGG annotations. note that one can instead use reactome, GOBP, biocarta database
data(Pathways)

# run CLUE with a repeat of 3 times and a range from 2 to 13
set.seed(3)
clueObj <- runClue(Tc=adipocyte, annotation=Pathways.KEGG, rep=3, kRange=13)

# visualize the evaluation outcome
xl <- "Number of clusters"
yl <- "Enrichment score"
boxplot(clueObj$evlMat, col=rainbow(ncol(clueObj$evlMat)), las=2, xlab=xl, ylab=yl, main="CLUE")
abline(v=(clueObj$maxK-1), col=rgb(1,0,0,.3))

# generate the optimal clustering results
best <- clustOptimal(clueObj, rep=5, mfrow=c(3, 3))

# list enriched clusters
best$enrichList

# obtain the optimal clustering object
best$clustObj
```
