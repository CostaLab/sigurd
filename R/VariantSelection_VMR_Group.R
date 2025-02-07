#'VariantSelection_VMR_Group
#'@description
#'This is an adaption of the VMR based selection of variants from Miller et al. 
#'This selection process is designed for ATAC data and not scRNAseq data. 
#'The variants are first selected as described. Then the average allele frequency between the the two groups of interest are compared.
#'If the allele frequency in group1 is higher than the required threshold, the variant is retained.
#'@importFrom SummarizedExperiment rowRanges assays
#'@importFrom data.table data.table
#'@importFrom dplyr sym
#'@importFrom MatrixGenerics rowVars
#'@param SE SummarizedExperiment object.
#'@param stabilize_variance Should the variance be stabilized by using the mean allele frequency for cells with a low coverage? Coverage threshold is set by low_coverage_threshold.
#'@param low_coverage_threshold Cells below this threshold are set to the mean. 
#'@param minimum_fw_rev_reads How many forward and reverse reads should a cell have? Default = 2
#'@param group_of_interest The group of interest in the column data.
#'@param group1 The group of interest.
#'@param group2 The second group.
#'@param group_factor How much higher should the VAF be in group 1 comapred to group 2? Default = 5
#'@param verbose Should the function be verbose?
#'@export
VariantSelection_VMR_Group <- function(SE, stabilize_variance = TRUE, low_coverage_threshold = 10, minimum_fw_rev_reads = 2, 
                                       group_of_interest = NULL, group1 = NULL, group2 = NULL, group_factor = 5, verbose = TRUE){
  
  if(!group_of_interest %in% colnames(SummarizedExperiment::colData(SE))) stop("Error: Your group_of_interest is not in the colData.")
  if(!group1 %in% SummarizedExperiment::colData(SE)[, group_of_interest]) stop("Error: Your group1 is not in the group_of_interest.")
  if(!group2 %in% SummarizedExperiment::colData(SE)[, group_of_interest]) stop("Error: Your group2 is not in the group_of_interest.")
  
  # Processing group1 and group2
  cells_group1 <- SummarizedExperiment::colData(SE)[, group_of_interest, drop = FALSE]
  cells_group1 <- cells_group1[cells_group1[, group_of_interest] == group1, , drop = FALSE]
  SE_group1    <- SE[, rownames(cells_group1)]
  cells_group2 <- SummarizedExperiment::colData(SE)[, group_of_interest, drop = FALSE]
  cells_group2 <- cells_group2[cells_group2[, group_of_interest] == group2, , drop = FALSE]
  SE_group2    <- SE[, rownames(cells_group2)]
  
  if(verbose) print("We get the average allele frequency for group 1.")
  fraction_group1 <- SummarizedExperiment::assays(SE_group1)[["alts"]] / SummarizedExperiment::assays(SE_group1)[["coverage"]]
  fraction_group1[is.na(fraction_group1)] <- 0
  mean_af_group1 <- Matrix::rowSums(SummarizedExperiment::assays(SE_group1)[["alts"]], na.rm = TRUE) / Matrix::rowSums(SummarizedExperiment::assays(SE_group1)[["coverage"]], na.rm = TRUE)
  mean_af_group1[is.na(mean_af_group1)] <- 0
  mean_af_group2 <- Matrix::rowSums(SummarizedExperiment::assays(SE_group2)[["alts"]], na.rm = TRUE) / Matrix::rowSums(SummarizedExperiment::assays(SE_group2)[["coverage"]], na.rm = TRUE)
  mean_af_group2[is.na(mean_af_group2)] <- 0
  mean_cov_group1 <- Matrix::rowMeans(SummarizedExperiment::assays(SE_group1)[["coverage"]])
  mean_cov_group2 <- Matrix::rowMeans(SummarizedExperiment::assays(SE_group2)[["coverage"]])
  
  if(verbose) print("We calculate the concordance.")
  reads_forward_group1 <- SummarizedExperiment::assays(SE_group1)[["forward"]]
  reads_reverse_group1 <- SummarizedExperiment::assays(SE_group1)[["reverse"]]
  dt <- merge(data.table::data.table(summary(reads_forward_group1)), 
              data.table::data.table(summary(reads_reverse_group1)), 
              by.x = c("i", "j"), by.y = c("i", "j"), 
              all = TRUE)
  dt <- dt[dt[,x.x] > 0 | dt[,x.y] > 0,]
  dt <- data.table::data.table(variant = rownames(SE_group1)[dt[[1]]],
                               cell_id = dt[[2]],
                               fw = dt[[3]], rev = dt[[4]])
  concordance <- dt[, .(cor = suppressWarnings(stats::cor(c(fw), c(rev), method = "pearson", use = "pairwise.complete"))), by = list(variant)]
  
  if(stabilize_variance){
    if(verbose) print("We stabilize the variance.")
    coverage <- SummarizedExperiment::assays(SE_group1)[["coverage"]]
    coverage_test <- which(data.matrix(coverage < low_coverage_threshold), arr.ind = TRUE)
    means_for_stabilization <- mean_af_group1[coverage_test[,1]]
    # This matrix contains all the positions that have a low coverage.
    low_coverage_positions <- 1 - Matrix::sparseMatrix(
      i = c(coverage_test[,1], dim(fraction_group1)[1]),
      j = c(coverage_test[,2], dim(fraction_group1)[2]),
      x = 1
    )
    # This matrix has the mean values for the low coverage positions.
    means_for_stabilization_matrix <- sparseMatrix(
      i = c(coverage_test[,1], dim(fraction_group1)[1]),
      j = c(coverage_test[,2], dim(fraction_group1)[2]),
      x = c(means_for_stabilization, 0)
    )
    fraction_group1 <- fraction_group1 * low_coverage_positions + means_for_stabilization_matrix
  }
  variants_variance <- MatrixGenerics::rowVars(fraction_group1)
  
  if(verbose) print(paste0("We check if a cells has more than ", minimum_fw_rev_reads, " reads."))
  detected <- (SummarizedExperiment::assays(SE_group1)[["forward"]] >= minimum_fw_rev_reads) + (assays(SE_group1)[["reverse"]] >= minimum_fw_rev_reads)
  
  # We get the variants in all objects.
  variants_kept <- intersect(rownames(detected), names(variants_variance))
  variants_kept <- intersect(variants_kept, concordance$variant)
  variants_kept <- intersect(variants_kept, rownames(fraction_group1))
  variants_kept <- intersect(variants_kept, names(mean_af_group1))
  variants_kept <- intersect(variants_kept, names(mean_af_group2))
  variants_kept <- intersect(variants_kept, names(mean_cov_group1))
  detected <- detected[variants_kept, ]
  variants_variance <- variants_variance[variants_kept]
  concordance <- subset(concordance, variant %in% variants_kept)
  fraction_group1 <- fraction_group1[variants_kept,]
  mean_af_group1 <- mean_af_group1[variants_kept]
  mean_af_group2 <- mean_af_group2[variants_kept]
  mean_cov_group1 <- mean_cov_group1[variants_kept]
  mean_cov_group2 <- mean_cov_group2[variants_kept]
  
  voi <- data.frame(
    variant = rownames(detected),
    vmr = variants_variance / (mean_af_group1 + 0.00000000001),
    vmr_log10 = log10(variants_variance / (mean_af_group1 + 0.00000000001)),
    mean_group1 = round(mean_af_group1, 7),
    mean_group2 = round(mean_af_group2, 7),
    mean_af_check = mean_af_group1 > (group_factor * mean_af_group2),
    variance = round(variants_variance, 7),
    cells_detected = Matrix::rowSums(detected == 2), # We only select cells that have enough reads for both directions.
    strand_correlation = concordance$cor, 
    mean_coverage_group1 = mean_cov_group1,
    mean_coverage_group2 = mean_cov_group2,
    stringsAsFactors = FALSE, row.names = rownames(detected)
  )
  return(voi)
}
