#'Calculate the variant burden per cell using big matrices.
#'@import Matrix SummarizedExperiment bigmemory
#'@param se SummarizedExperiment object
#'@export
VariantBurden_big <- function(se){
  burden <- colSums(assays(se)[["fraction"]][,])
  colData(se)[,"Burden"] <- burden
  return(se)
}
