#'AllelFrequencyFoldChange
#'@description
#'Check if a cell is supported by a set of variants.
#'@importFrom SummarizedExperiment assays colData rowData
#'@importFrom Matrix rowMeans
#'@param SE SummarizedExperiment object.
#'@param VOI The variants variants to be analyzed. If NULL all are used.
#'@param group_of_interest The column data that divides the cells.
#'@param group1 The first group.
#'@param group2 The second group.
#'@param minimum_foldchange Minimum fold change.
#'@param maximum_foldchange Maximum fold change.
#'@param minimum_coverage Minimum coverage for a variant.
#'@param minimum_allele_freq Minimum allele frequency in both groups.
#'@param maximum_allele_freq Maximum allele frequency in both groups.
#'@param verbose Should the function be verbose? Default = TRUE
#'@export
AllelFrequencyFoldChange <- function(SE, VOI = NULL, group_of_interest, group1 = "group1", group2 = "group2", maximum_foldchange = NULL, minimum_foldchange = NULL, minimum_coverage = NULL, minimum_allele_freq = NULL, maximum_allele_freq = NULL, verbose = FALSE){
  if(is.null(VOI)){
    VOI <- row.names(SE)
  } else{
    if(!all(VOI %in% rownames(SE))){
      stop("Not all your variants are in the data.")
    } else{
      SE <- SE[VOI,]
    }
  }
  if(!group_of_interest %in% colnames(SummarizedExperiment::colData(SE))){
    stop("Your group of interest is not in the data.")
  }
  if(!group1 %in% SummarizedExperiment::colData(SE)[,group_of_interest]){
    stop("Your group 1 is not in the data.")
  }
  if(!group2 %in% SummarizedExperiment::colData(SE)[,group_of_interest]){
    stop("Your group 2 is not in the data.")
  }
  cell_data <- SummarizedExperiment::colData(SE)
  cells_group1 <- subset(cell_data, cell_data[, group_of_interest] == group1)
  mean_allele_frequency_group1 <- SummarizedExperiment::assays(SE)[["fraction"]][, cells_group1$Cell, drop = FALSE]
  mean_allele_frequency_group1 <- Matrix::rowMeans(mean_allele_frequency_group1)
  cells_group2 <- subset(cell_data, cell_data[, group_of_interest] == group2)
  mean_allele_frequency_group2 <- SummarizedExperiment::assays(SE)[["fraction"]][, cells_group2$Cell, drop = FALSE]
  mean_allele_frequency_group2 <- Matrix::rowMeans(mean_allele_frequency_group2)

  result <- data.frame(Variant = names(mean_allele_frequency_group1), Coverage = SummarizedExperiment::rowData(SE)[,"Depth"], Group1 = mean_allele_frequency_group1, Group2 = mean_allele_frequency_group2)
  result$FoldChange <- apply(result[,3:4], 1, function(x){
    vaf_1 <- min(x[1], abs(1 - x[1]))
    vaf_2 <- min(x[2], abs(1 - x[2]))
    change_in_frequency <- abs(x[1] - x[2])
    base_value <- min(vaf_1, vaf_2)
    fc <- (base_value + change_in_frequency) / base_value
    return(fc)
  })
  if(!is.null(minimum_foldchange)){
    result <- subset(result, FoldChange > minimum_foldchange)
  }
  if(!is.null(maximum_foldchange)){
    result <- subset(result, FoldChange < maximum_foldchange)
  }
  if(!is.null(minimum_coverage)){
    result <- subset(result, Coverage > minimum_coverage)
  }
  if(!is.null(minimum_allele_freq)){
    result <- subset(result, Group1 > minimum_allele_freq & Group2 > minimum_allele_freq)
  }
  if(!is.null(maximum_allele_freq)){
    result <- subset(result, Group1 < maximum_allele_freq & Group2 < maximum_allele_freq)
  }
  return(result)
}
