\name{ClueR-package}
\alias{ClueR-package}
\alias{CLUE}
\docType{package}
\title{
CLUster Evaluation R package
}
\description{
CLUster Evaluation (or "CLUE") is an R package for identifying optimal number of clusters in a given time-course dataset clustered by cmeans or kmeans algorithms. It relies on a reference annotation set to test for enrichment in
each cluster using Fisher's Exact Test and then test for overall enrichment of the entire clusters using Fisher's
combined probability test.

CLUE is designed for analyzing time-course phosphoproteomics dataset using kinase-substrate annotation as reference. However, it can be applied to time-course microarray dataset as well by replacing the kinase-substrate annotation with gene sets annotation.
}
\details{
\tabular{ll}{
Package: \tab CLUE\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-2-1\cr
License: \tab GPL-3\cr
}

}
\author{
Pengyi Yang
}
\references{
Pengyi Yang, Xiaofeng Zheng, Vivek Jayaswal, Guang Hu, Yee Hwa Yang, Raja Jothi, Detecting key signaling events from time-series phosphoproteomics data, submitted.
}
