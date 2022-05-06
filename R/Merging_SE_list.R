#'Merging list of SummarizedExperiment objects.
#'@import BiocGenerics
#'@param se SummarizedExperiment object
#'@export
Merging_SE_list <- function(se){
    result <- do.call("cbind", se)
}
