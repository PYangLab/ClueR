\name{ClueR-package}
\alias{ClueR-package}
\alias{CLUE}
\alias{ClueR}
\docType{package}
\title{
CLUster Evaluation R package
}
\description{
CLUster Evaluation (or "CLUE") is an R package for detecting kinases or pathways from a given time-series phosphoproteomics or gene expression dataset clustered by cmeans or kmeans algorithms. It firstly identifies the optimal number of clusters in the time-servies dataset; Then, it partition the dataset based on the optimal number of clusters determined in the first step; It finally detects kinases or pathways enriched in each cluster from optimally partitioned dataset.

The above three steps rely extensively on Fisher's exact test, Fisher's combined statistics, cluster regularisations, and they are performed against a user-specified reference annotation database such phosphoSitePlus in the case of phosphoproteomics data or KEGG in the case of gene expression data. There is a selection of built-in annotation databases for both phosphoproteomics data and gene expression data but users can supply their own annotation database.

CLUE was initially designed for analysing time-course phosphoproteomics dataset using kinase-substrate annotation as reference (e.g. PhosphoSitePlus). It is now extended to identify key pathways from time-series microarray, RNA-seq or proteomics datasets by searching and testing against gene set annotation databases such as KEGG, GO, or Reactome etc.

Previously published phosphoproteomics dataset and gene expression dataset are included in the package to demonstrate how to use CLUE package.

See help from the main function by typing '?runClue' for more details and examples on how to use CLUE.

You can also install the latest development version from github with:

devtools::install_github("PengyiYang/ClueR")

Make sure that you have Rtools install in your system for building the package from the source.
}
\details{
\tabular{ll}{
Package: \tab CLUE\cr
Type: \tab Package\cr
Version: \tab 1.2\cr
Date: \tab 2017-04-28\cr
License: \tab GPL-3\cr
}

}
\author{
Pengyi Yang
}
\references{
Yang P, Zheng X, Jayaswal V, Hu G, Yang JYH, Jothi R (2015) Knowledge-Based Analysis for Detecting Key Signaling Events from Time-Series Phosphoproteomics Data. PLoS Comput Biol 11(8): e1004403.
}
