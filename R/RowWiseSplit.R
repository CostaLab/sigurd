#'RowWiseSplit
#'@description
#'Performing the correlation or Fisher test association for a SummarizedExperiment object requires extreme amounts of memory.
#'To reduce the amount of memory necessary, we instead get the individual rows from the consensus assay.
#'We can then remove the NoCalls (no reads) from the individual vectors, further reducing the amount of memory needed.
#'
#'When the NoCalls are removed, the repective cells are also removed if a different assay is selected.
#'@importFrom parallel mclapply
#'@importFrom SummarizedExperiment assays
#'@param se SummarizedExperiment object.
#'@param n_cores Number of cores to use.
#'@param remove_nocalls Do you want to remove NoCall cells?
#'@param assay_to_split Which assay to you want to split?
#'@export
RowWiseSplit <- function(se, n_cores = 1, assay_to_split = "consensus", remove_nocalls = TRUE){
  selected_assay <- SummarizedExperiment::assays(se)[[assay_to_split]]
  consensus <- SummarizedExperiment::assays(se)[["consensus"]]
  consensus_list <- parallel::mclapply(rownames(se), SeparatingMatrixToList, total_matrix = selected_assay, assay_to_split = assay_to_split, consensus = consensus, remove_nocalls = remove_nocalls, mc.cores = n_cores)
  names(consensus_list) <- rownames(se)
  return(consensus_list)
}
