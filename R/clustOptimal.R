#' Generate optimal clustering
#' 
#' Takes a clue output and generate the optimal clustering of the time-course data.
#' @param clueObj the output from runClue.
#' @param rep number of times (default is 5) the clustering is to be repeated to find the best clustering result.
#' @param user.maxK user defined optimal k value for generating optimal clustering. If not provided, the optimal k that is identified by clue will be used.
#' @param visualize a boolean parameter indicating whether to visualize the clustering results.
#' @param effectiveSize the size of kinase-substrate groups to be considered for calculating enrichment. Groups that are too small
#' or too large will be removed from calculating overall enrichment of the clustering.
#' @param pvalueCutoff a pvalue cutoff for determining which kinase-substrate groups to be included in calculating overall enrichment of the clustering.
#' @param universe the universe of genes/proteins/phosphosites etc. that the enrichment is calculated against. The default are the row names of the dataset.
#' @param mfrow control the subplots in graphic window.
#' @return return a list containing optimal clustering object and enriched kinases or gene sets. 
#'
#' @import stats
#'
#' @export
#' @examples
#' # simulate a time-series data with 4 distinctive profile groups and each group with
#' # a size of 50 phosphorylation sites.
#' simuData <- temporalSimu(seed=1, groupSize=50, sdd=1, numGroups=4)
#' 
#' # create an artificial annotation database. Generate 20 kinase-substrate groups each
#' # comprising 10 substrates assigned to a kinase.
#' # among them, create 4 groups each contains phosphorylation sites defined to have the
#' # same temporal profile.
#' kinaseAnno <- list()
#' groupSize <- 50
#' for (i in 1:4) {
#'  kinaseAnno[[i]] <- paste("p", (groupSize*(i-1)+1):(groupSize*(i-1)+10), sep="_")
#' }
#'
#' for (i in 5:20) {
#'  set.seed(i)
#'  kinaseAnno[[i]] <- paste("p", sample.int(nrow(simuData), size = 10), sep="_")
#' }
#' names(kinaseAnno) <- paste("KS", 1:20, sep="_")
#'
#' # run CLUE with a repeat of 2 times and a range from 2 to 7
#' set.seed(1)
#' clueObj <- runClue(Tc=simuData, annotation=kinaseAnno, rep=5, kRange=2:7)
#' 
#' # visualize the evaluation outcome
#' xl <- "Number of clusters"
#' yl <- "Enrichment score"
#' boxplot(clueObj$evlMat, col=rainbow(ncol(clueObj$evlMat)), las=2, xlab=xl, ylab=yl, main="CLUE")
#' abline(v=(clueObj$maxK-1), col=rgb(1,0,0,.3))
#' 
#' # generate optimal clustering results using the optimal k determined by CLUE
#' best <- clustOptimal(clueObj, rep=3, mfrow=c(2, 3))
#' 
#' # list enriched clusters
#' best$enrichList
#' 
#' # obtain the optimal clustering object
#' \donttest{
#'   best$clustObj
#' }
#' 
#' 
clustOptimal <- function(clueObj, rep=5, user.maxK=NULL, effectiveSize=NULL, pvalueCutoff=0.05, visualize=TRUE, universe=NULL, mfrow = c(1, 1)) {
  
  bstPvalue <- 1
  bst.cobj <- c()
  bst.evaluation <- c()
  for(i in 1:rep) {
    cobj <- c()
    if (clueObj$clustAlg == "cmeans") {
      if (is.null(user.maxK)) {
        # use clue determined max k value
        cobj <- e1071::cmeans(clueObj$Tc, centers=clueObj$maxK, iter.max=50, m=1.25)
      } else {
        # use user specific k value
        cobj <- e1071::cmeans(clueObj$Tc, centers=user.maxK, iter.max=50, m=1.25)
      }
    } else {
      if (is.null(user.maxK)) {
        # use clue determined max k value
        cobj <- stats::kmeans(clueObj$Tc, centers=clueObj$maxK, iter.max=50)
        # set the membership as the normalised distance to the mean of the cluster
        cobj$membership <- matrix(NA, nrow=nrow(clueObj$Tc), ncol=nrow(cobj$centers))
        
        for(i in 1:ncol(cobj$membership)) {
          cobj$membership[,i] <- (as.numeric(stats::cor(t(clueObj$Tc), cobj$centers[i,])) + 1)/2
        }
        rownames(cobj$membership) <- names(cobj$cluster)
        
      } else {
        # use user specific k value
        cobj <- stats::kmeans(clueObj$Tc, centers=user.maxK, iter.max=50)
        # set the membership as the normalised distance to the mean of the cluster
        cobj$membership <- matrix(NA, nrow=nrow(clueObj$Tc), ncol=nrow(cobj$centers))
        
        for(i in 1:ncol(cobj$membership)) {
          cobj$membership[,i] <- (as.numeric(cor(t(clueObj$Tc), cobj$centers[i,])) + 1)/2
        }
        rownames(cobj$membership) <- names(cobj$cluster)
        
      }
    }
    
    evaluation <- c()
    currentPvalue <- c()
    if (!is.null(clueObj$effectiveSizeK) & !is.null(clueObj$pvalueCutoff)) {
      evaluation <- clustEnrichment(cobj, clueObj$annotation, effectiveSize, pvalueCutoff, universe)
      currentPvalue <- evaluation$fisher.pvalue
    } else if (!is.null(clueObj$effectiveSizeK) & is.null(clueObj$pvalueCutoff)) {
      evaluation <- clustEnrichment(cobj, clueObj$annotation, effectiveSize, clueObj$pvalueCutoff, universe)
      currentPvalue <- evaluation$fisher.pvalue
    } else if (is.null(clueObj$effectiveSizeK) & !is.null(clueObj$pvalueCutoff)) {
      evaluation <- clustEnrichment(cobj, clueObj$annotation, clueObj$effectiveSize, pvalueCutoff, universe)
      currentPvalue <- evaluation$fisher.pvalue
    } else {
      evaluation <- clustEnrichment(cobj, clueObj$annotation, clueObj$effectiveSize, clueObj$pvalueCutoff, universe)
      currentPvalue <- evaluation$fisher.pvalue
    }
      
    if (currentPvalue < bstPvalue) {
      bstPvalue <- currentPvalue
      bst.cobj <- cobj
      bst.enrichList <- evaluation$enrich.list
    }
  }
  
  if(visualize) {
    fuzzPlot(clueObj$Tc, clustObj = bst.cobj, mfrow)
  }
  
  results <- list()
  results$clustObj <- bst.cobj
  results$enrichList <- bst.enrichList
  return(results)
}
