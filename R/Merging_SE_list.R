#'Merging list of SummarizedExperiment objects.
#'@description
#'This function is a wrapper for do.all("cbind", se).
#'@import BiocGenerics
#'@param se SummarizedExperiment object
#'@export
Merging_SE_list <- function(se){
    result <- do.call("cbind", se)
}
