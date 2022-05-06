#'Merging list of SummarizedExperiment objects.
#'@import BiocGenerics
#'@param x
Merging_SE_list <- function(x){
    result <- do.call("cbind", x)
}
