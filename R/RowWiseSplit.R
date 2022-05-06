#'We get the filtered results, we have from the 01_FilteringSampleWise. We now split the results into the
#'individual columns and save them in a list object.
#'@import Matrix SummarizedExperiment parallel
#'@param se SummarizedExperiment object.
#'@param n_cores Number of cores to use.
#'@export
RowWiseSplit <- function(se, n_cores = 1){
  consensus <- assays(se)$consensus
  consensus_list <- mclapply(rownames(se), SeparatingMatrixToList, total_matrix = consensus, mc.cores = n_cores)
  names(consensus_list) <- rownames(se)
  return(consensus_list)
}
