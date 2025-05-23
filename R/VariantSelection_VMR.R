#'VariantSelection_VMR
#'@description
#'This is an adaption of the VMR based selection of variants from Miller et al. 
#'This selection process is designed for ATAC data and not scRNAseq data. 
#'@importFrom SummarizedExperiment rowRanges assays
#'@importFrom data.table data.table
#'@importFrom dplyr sym
#'@importFrom MatrixGenerics rowVars
#'@param SE SummarizedExperiment object.
#'@param stabilize_variance Should the variance be stabilized by using the mean allele frequency for cells with a low coverage? Coverage threshold is set by low_coverage_threshold.
#'@param low_coverage_threshold Cells below this threshold are set to the mean. 
#'@param minimum_fw_rev_reads How many forward and reverse reads should a cell have? Default = 2
#'@param verbose Should the function be verbose?
#'@export
VariantSelection_VMR <- function(SE, stabilize_variance = TRUE, low_coverage_threshold = 10, minimum_fw_rev_reads = 2, verbose = TRUE){
  if(verbose) print("We get the average allele frequency for the entire data set.")
  fraction <- SummarizedExperiment::assays(SE)[["alts"]] / SummarizedExperiment::assays(SE)[["coverage"]]
  fraction[is.na(fraction)] <- 0
  mean_af <- Matrix::rowSums(SummarizedExperiment::assays(SE)[["alts"]]) / Matrix::rowSums(SummarizedExperiment::assays(SE)[["coverage"]])
  mean_af[is.na(mean_af)] <- 0
  mean_cov <- Matrix::rowMeans(SummarizedExperiment::assays(SE)[["coverage"]])

  if(verbose) print("We calculate the concordance.")
  reads_forward <- SummarizedExperiment::assays(SE)[["forward"]]
  reads_reverse <- SummarizedExperiment::assays(SE)[["reverse"]]
  dt <- merge(data.table::data.table(summary(reads_forward)), 
              data.table::data.table(summary(reads_reverse)), 
              by.x = c("i", "j"), by.y = c("i", "j"), 
              all = TRUE)
  dt <- dt[dt[,x.x] > 0 | dt[,x.y] > 0,]
  dt <- data.table::data.table(variant = rownames(SE)[dt[[1]]],
                               cell_id = dt[[2]],
                               fw = dt[[3]], rev = dt[[4]])
  concordance <- dt[, .(cor = suppressWarnings(stats::cor(c(fw), c(rev), method = "pearson", use = "pairwise.complete"))), by = list(variant)]
  
  if(stabilize_variance){
    if(verbose) print("We stabilize the variance.")
    coverage <- SummarizedExperiment::assays(SE)[["coverage"]]
    coverage_test <- which(data.matrix(coverage < low_coverage_threshold), arr.ind = TRUE)
    means_for_stabilization <- mean_af[coverage_test[,1]]
    # This matrix contains all the positions that have a low coverage.
    low_coverage_positions <- 1 - Matrix::sparseMatrix(
      i = c(coverage_test[,1], dim(fraction)[1]),
      j = c(coverage_test[,2], dim(fraction)[2]),
      x = 1
    )
    # This matrix has the mean values for the low coverage positions.
    means_for_stabilization_matrix <- sparseMatrix(
      i = c(coverage_test[,1], dim(fraction)[1]),
      j = c(coverage_test[,2], dim(fraction)[2]),
      x = c(means_for_stabilization, 0)
    )
    fraction <- fraction * low_coverage_positions + means_for_stabilization_matrix
  }
  variants_variance <- MatrixGenerics::rowVars(fraction)
  names(variants_variance) <- rownames(fraction)

  if(verbose) print(paste0("We check if a cells has more than ", minimum_fw_rev_reads, " reads."))
  detected <- (SummarizedExperiment::assays(SE)[["forward"]] >= minimum_fw_rev_reads) + (assays(SE)[["reverse"]] >= minimum_fw_rev_reads)

  # We get the variants in all objects.
  variants_kept <- intersect(rownames(detected), names(variants_variance))
  variants_kept <- intersect(variants_kept, concordance$variant)
  variants_kept <- intersect(variants_kept, rownames(fraction))
  variants_kept <- intersect(variants_kept, names(mean_af))
  variants_kept <- intersect(variants_kept, names(mean_cov))
  detected <- detected[variants_kept, ]
  variants_variance <- variants_variance[variants_kept]
  concordance <- subset(concordance, variant %in% variants_kept)
  fraction <- fraction[variants_kept,]
  mean_af <- mean_af[variants_kept]
  mean_cov <- mean_cov[variants_kept]

  voi <- data.frame(
    variant = rownames(detected),
    vmr = variants_variance / (mean_af + 0.00000000001),
    vmr_log10 = log10(variants_variance / (mean_af + 0.00000000001)),
    mean = round(mean_af, 7),
    variance = round(variants_variance, 7),
    cells_detected = Matrix::rowSums(detected == 2), # We only select cells that have enough reads for both directions.
    strand_correlation = concordance$cor, 
    mean_coverage = mean_cov,
    stringsAsFactors = FALSE, row.names = rownames(detected)
  )
  return(voi)
}
