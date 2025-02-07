#'Filtering_VAF_Threshold
#'@description
#'We determine the number of mutated cells for a given variant using a vector of VAF thresholds.
#'@importFrom SummarizedExperiment assays
#'@param SE SummarizedExperiment object.
#'@param variant Which variant should be evaluated? Cannot be more than 1 variant at a time.
#'@param thresholds The thresholds to be evaluated. Default between 0.01 and 0.99 in steps of 0.01.
#'@param verbose Should the function be verbose?
#'@export
Filtering_VAF_Threshold <- function(SE, variant = NULL, thresholds = NULL, verbose = TRUE){
  if(!is.null(thresholds)){
    if(min(thresholds) < 0) stop("Error: Minimum threshold is smaller than 0.")
    if(max(thresholds) > 1) stop("Error: Maximum threshold is larger than 1.")
  }
  if(is.null(variant)) stop("Error: You must provide a valid variant name.")
  if(length(variant) > 1) stop("Error: You cannot provide more than 1 variant.")
  if(!variant %in% rownames(SE)) stop("Error: Your variant is not in the SE object.")
  if(!"fraction" %in% names(SummarizedExperiment::assays(SE))) stop("Error: Your SE object does not have a fraction assay.")

  # We generate the thresholds if they are not provided.
  if(is.null(thresholds)){
    if(verbose) print("Use the default thresholds of 0.01 to 0.99.")
    thresholds <- seq(from = 0, to = 0.99, by = 0.01)
  }

  # We retrieve the VAF value for the variant of interest.
  vaf_values <- SummarizedExperiment::assays(SE)[["fraction"]][variant ,]
  vaf_filtered <- c()
  for(i in 1:length(thresholds)){
    print(paste0("Threshold ", thresholds[i]))
    vaf_filtered[paste0("Threshold_", thresholds[i])] <- sum(vaf_values > thresholds[i])
  }
  return(vaf_filtered)
}
