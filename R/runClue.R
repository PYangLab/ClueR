#' Run CLUster Evaluation
#' 
#' Takes in a time-course matrix and test for enrichment of the clustering using cmeans or kmeans clustering algorithm with a reference annotation.
#' @param Tc a numeric matrix to be clustered. The columns correspond to the time-course and the rows correspond to phosphorylation sites.
#' @param annotation a list with names correspond to kinases and elements correspond to substrates belong to each kinase.
#' @param rep number of times the clustering is to be applied. This is to account for variability in the clustering algorithm. Default is 5.
#' @param kRange the range of k to be tested for clustering. Default is 2:10
#' @param clustAlg the clustering algorithm to be used. The default is cmeans clustering.
#' @param effectiveSize the size of annotation groups to be considered for calculating enrichment. Groups that are too small
#' or too large will be removed from calculating overall enrichment of the clustering.
#' @param pvalueCutoff a pvalue cutoff for determining which kinase-substrate groups to be included in calculating overall enrichment of the clustering.
#' @param alpha a regularisation factor for penalizing large number of clusters.
#' @return a clue output that contains the input parameters used for evaluation and the evaluation results. Use ls(x) to see details of output. 'x' be the output here.
#' @export
#' @examples
#' ## Example 1. Running CLUE with a simulated phosphoproteomics data
#' 
#' ## simulate a time-series phosphoproteomics data with 4 clusters and
#' ## each cluster with a size of 100 phosphosites
#' simuData <- temporalSimu(seed=1, groupSize=100, sdd=1, numGroups=4)
#' 
#' ## create an artificial annotation database. Specifically, Generate 50
#' ## kinase-substrate groups each comprising 20 substrates assigned to a kinase. 
#' ## Among them, create 5 groups each contains phosphosites defined 
#' ## to have the same temporal profile.
#' 
#' kinaseAnno <- list()
#' groupSize <- 100
#' for (i in 1:5) {
#'   kinaseAnno[[i]] <- paste("p", (groupSize*(i-1)+1):(groupSize*(i-1)+20), sep="_")
#' }
#' 
#' for (i in 6:50) {
#'   set.seed(i)
#'   kinaseAnno[[i]] <- paste("p", sample.int(nrow(simuData), size = 20), sep="_")
#' }
#' names(kinaseAnno) <- paste("KS", 1:50, sep="_")
#' 
#' ## run CLUE with a repeat of 3 times and a range from 2 to 8
#' set.seed(1)
#' cl <- runClue(Tc=simuData, annotation=kinaseAnno, rep=3, kRange=2:8)
#' 
#' ## visualize the evaluation outcome
#' boxplot(cl$evlMat, col=rainbow(8), las=2, xlab="# cluster", ylab="Enrichment", main="CLUE")
#' 
#' ## generate optimal clustering results using the optimal k determined by CLUE
#' best <- clustOptimal(cl, rep=3, mfrow=c(2, 3))
#' 
#' ## list enriched clusters
#' best$enrichList
#' 
#' ## obtain the optimal clustering object
#' \donttest{best$clustObj}
#' 
#' ## Example 2. Running CLUE with a phosphoproteomics dataset, discover optimal number of clusters, 
#' ## clustering data accordingly, and identify key kinases involved in each cluster.
#' 
#' ## load the human ES phosphoprotoemics data (Rigbolt et al. Sci Signal. 4(164):rs3, 2011)
#' data(hES)
#' # load the PhosphoSitePlus annotations (Hornbeck et al. Nucleic Acids Res. 40:D261-70, 2012)
#' # note that one can instead use PhosphoELM database by typing "data(PhosphoELM)".
#' data(PhosphoSite)
#' 
#' ## run CLUE with a repeat of 5 times and a range from 2 to 15
#' \donttest{set.seed(1)
#' cl <- runClue(Tc=hES, annotation=PhosphoSite.human, rep=5, kRange=2:15)
#' 
#' boxplot(cl$evlMat, col=rainbow(15), las=2, xlab="# cluster", ylab="Enrichment", main="CLUE")
#' 
#' best <- clustOptimal(cl, rep=3, mfrow=c(4, 4))
#' 
#' best$enrichList}
#' 
#' ## Example 3. Running CLUE with a gene expression dataset, discover optimal number of clusters, 
#' ## clustering data accordingly, and identify key pathway involved in each cluster.
#' 
#' ## load mouse adipocyte gene expression data 
#' # (Ma et al. Molecular and Cellular Biology. 2014, 34(19):3607-17)
#' data(adipocyte)
#' 
#' ## load the KEGG annotations
#' ## note that one can instead use reactome, GOBP, biocarta database
#' data(Pathways)
#' 
#' ## select genes that are differentially expressed during adipocyte differentiation
#' adipocyte.selected <- adipocyte[adipocyte[,"DE"] == 1,]
#' 
#' ## run CLUE with a repeat of 5 times and a range from 10 to 22
#' \donttest{
#' set.seed(3)
#' cl <- runClue(Tc=adipocyte.selected, annotation=Pathways.KEGG, rep=3, kRange=10:20)
#' 
#' xl <- "Number of clusters"
#' yl <- "Enrichment score"
#' boxplot(cl$evlMat, col=rainbow(ncol(cl$evlMat)), las=2, xlab=xl, ylab=yl, main="CLUE")}
#' 
runClue <- function(Tc, annotation, rep=5, kRange=2:10, clustAlg="cmeans", effectiveSize=c(5, 100), pvalueCutoff=0.05, alpha=0.5) {
  
  # standardize the matrix by row
  means <- apply(Tc, 1, mean)
  stds <- apply(Tc, 1, sd)
  tmp <- sweep(Tc, 1, means, FUN="-")
  Tc <- sweep(tmp, 1, stds, FUN="/")
  
  ## filter the annotation groups that has no entry from the Tc
  annotation.intersect <- lapply(annotation, intersect, rownames(Tc))
  annotation.filtered <- annotation.intersect[lapply(annotation.intersect, length) > 0]
  
  # apply CLUE
  repeat.list <- mclapply(1:rep, function(rp){
    cat("repeat", rp, "\n");
    enrichment <- c()
    for (k in kRange) {
      clustered <- c()
      if (clustAlg == "cmeans") {
        clustered <- cmeans(Tc, centers=k, iter.max=50, m=1.25)
      } else if (clustAlg == "kmeans"){
        clustered <- kmeans(Tc, centers=k, iter.max=50)
      } else {
        print("Unknown clustering algorithm specified. Using cmeans clustering instead")
        clustered <- cmeans(Tc, centers=k, iter.max=50, m=1.25)
      }
      
      # compute fisher's p-value for each cluster
      evaluate <- clustEnrichment(clustered, annotation.filtered, effectiveSize, pvalueCutoff)
      fisher.pvalue <- evaluate$fisher.pvalue
      
      # compute clustering enrichment (this is regularised by the number of clusters) 
      escore <- -log10(fisher.pvalue) - alpha * nrow(clustered$centers)
      
      enrichment <- c(enrichment, escore)
    }
    enrichment
  })
  
  # combine the multiple testing results
  x <- do.call(rbind, repeat.list)
  # transform the pvalue 
  #x.transform <- -log10(x)
  # scale the values into [0, 1]
  x.normalize <- (x - min(x)) / (max(x) - min(x))
  rownames(x.normalize) <- paste("repeat", 1:rep, sep="")
  colnames(x.normalize) <- paste("k", kRange, sep="=")
  
  # identify the k that maximize the enrichment
  maxK <- which.max(apply(x.normalize, 2, median)) + (kRange[1] - 1)
  
  # return the evaluation results
  result <- list()
  # input parameters
  result$Tc <- Tc
  result$annotation <- annotation.filtered
  result$clustAlg <- clustAlg
  result$effectiveSize <- effectiveSize
  result$pvalueCutoff <- pvalueCutoff
  # clue parameters
  result$evlMat <- x.normalize
  result$maxK <- maxK
  
  return(result)
}
