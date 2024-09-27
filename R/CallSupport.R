#'CallSupport
#'@description
#'Check if a cell is supported by a set of variants.
#'@importFrom SummarizedExperiment colData
#'@param SE SummarizedExperiment object.
#'@param VOI_group1 The variants supporting the first group.
#'@param VOI_group2 The variants supporting the second group.
#'@param group1_name The name used for the first group.
#'@param group2_name The name used for the second group.
#'@param min_mutated_reads The minimum number of mutated reads in a cell supporting a group.
#'@param min_reads Minimum number of reads per cell for a classification.
#'@param group_of_interest The column data that divides the cells.
#'@param group_factor How much higher has the mean allele frequency to be in group 1 when compared to group 2 and vice versa? Can be a vector of length 2.
#'@param return_nonsupport Should nonsupporting cells be return too? Default FALSE.
#'@param verbose Should the function be verbose? Default = TRUE
#'@export
CallSupport <- function(SE, VOI_group1, VOI_group2, group1_name = "group1", group2_name = "group2", min_mutated_reads = 3, min_reads = 30, group_factor = NULL, verbose = TRUE, return_nonsupport = FALSE){
  if(!is.numeric(min_reads) & length(min_reads) != 1){
    stop("Your minimum number of reads is not a single number.")
  }
  if(!all(VOI_group1 %in% rownames(SE))){
    stop("No all your VOIs for group 1 are in the data.")
  }
  if(!all(VOI_group2 %in% rownames(SE))){
    stop("No all your VOIs for group 2 are in the data.")
  }
  if(!is.numeric(group_factor)){
    stop("Your group factor is not a number.")
  }
  if(!length(group_factor) %in% 1:2){
    stop("Your group factor is not 1 or to factors.")
  }
  if(!is.numeric(min_mutated_reads) & length(min_mutated_reads) != 1){
    stop("Your minimum number of mutated reads is not a single number.")
  }
  if(!is.character(group1_name) & length(group1_name) != 1){
    stop("Your name for group 1 is not a single string name.")
  }
  if(!is.character(group2_name) & length(group1_name) != 1){
    stop("Your name for group 2 is not a single string name.")
  }
  group_factor1 <- group_factor[1]
  group_factor2 <- group_factor[1]
  if(length(group_factor) == 2){
    group_factor2 <- group_factor[2]
  }
  supporting_reads_group1 <- Matrix::colSums(SummarizedExperiment::assays(SE)[["alts"]][VOI_group1, , drop = FALSE]) + Matrix::colSums(SummarizedExperiment::assays(SE)[["refs"]][VOI_group2, , drop = FALSE])
  supporting_reads_group2 <- Matrix::colSums(SummarizedExperiment::assays(SE)[["alts"]][VOI_group2, , drop = FALSE]) + Matrix::colSums(SummarizedExperiment::assays(SE)[["refs"]][VOI_group1, , drop = FALSE])
  cells_df <- data.frame(Cell = names(supporting_reads_group1), Group1_Support = supporting_reads_group1, Group2_Support = supporting_reads_group2)
  cells_df$Support <- apply(cells_df[,2:3], 1, function(x){
    if(sum(x) < min_reads){
      return("NoCoverage")
    } else if(x[1] > min_mutated_reads & x[1] > (group_factor1 * x[2])){
      return(group1_name)
    } else if(x[2] > min_mutated_reads & x[2] > (group_factor2 * x[1])){
      return(group2_name)
    } else{
      return("NoSupport")
    }
  })
  if(!return_nonsupport){
    cells_df <- subset(cells_df, Support %in% c(group1_name, group2_name))
  }
  return(cells_df)
}
