### Cluster Evaluation R package (ClueR) for detecting key signaling events from time-series phosphoproteomics data

##### Description
CLUster Evaluation (or "CLUE") is an R package for identifying optimal number of clusters in a given time-course dataset clustered by cmeans or kmeans algorithms. It relies on a reference annotation set to test for enrichment in
each cluster using Fisher's Exact Test and then test for overall enrichment of the entire clusters using Fisher's
combined probability test.

CLUE is designed for analyzing time-course phosphoproteomics dataset using kinase-substrate annotation as reference. However, it can be applied to time-course microarray dataset as well by replacing the kinase-substrate annotation with gene sets annotation.

Download the latest releases [here](https://github.com/PengyiYang/ClueR/releases)

###### Examples

    ``` r
# load the human ES phosphoprotoemics data (Rigbolt et al. Sci Signal. 4(164):rs3, 2011)
data(hES)
# load the PhosphoSitePlus annotations (Hornbeck et al. Nucleic Acids Res. 40:D261-70, 2012)
data(PhosphoSite)

    ```
    
