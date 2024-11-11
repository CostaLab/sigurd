#'VariantWiseCorrelation
#'@description
#'We correlate the variants with each other using the Pearson correlation.
#'This function calls CalculateCorrelationPValue to perform the actual correlation.
#'@importFrom parallel mclapply
#'@importFrom stats p.adjust
#'@param variants_list List of fraction values.
#'@param n_cores Number of cores you want to use. Numeric.
#'@param p_value_adjustment Method for P value adjustment. See p.adjust for details.
#'@param verbose Should the function be verbose? Default = TRUE
#'@param value_type Are we using consensus or other information?
#'@export
VariantWiseCorrelation <- function(variants_list, n_cores = 1, p_value_adjustment = "fdr", value_type = "consensus", verbose = TRUE){
  # We correlate the somatic variants with each other and the MT variants.
  # Since we have tens of thousands of MT variants, we do not correlate them with each other.
  variants <- names(variants_list)

  results_total <- c()
  number_of_variants <- length(variants)
  for(i in 1:length(variants)){
    variant_use <- variants[i]
    if(verbose) print(paste0("Correlating Variant: ", variant_use, ", ", i, " out of ", number_of_variants))
    variants_values_use <- variants_list[[variant_use]]
    variants_list_use <- variants_list[names(variants_list) != variant_use]
    all_variants <- names(variants_list_use)
    results <- parallel::mclapply(X = all_variants, CalculateCorrelationPValue, variant_values = variants_values_use, all_variants_list = variants_list_use, mc.cores = n_cores, value_type = value_type)
    results <- do.call("rbind", results)
    results <- data.frame(Variant1 = variant_use, Variant2 = all_variants, P = results[,1], Corr = results[,2],
                          Cells_1_Alt = results[,3], Cells_1_Ref = results[,4], Cells_2_Alt = results[,5], Cells_2_Ref = results[,6])
    results_total <- rbind(results_total, results)
  }

  if(verbose) print("We remove the NA P values.")
  results_total <- results_total[!is.na(results_total$P),]

  if(verbose) print("We remove the negative corrlated SNPs.")
  results_total <- subset(results_total, Corr > 0)

  if(verbose) print("Remove duplicative results.")
  results_check <- apply(results_total[, c("Variant1", "Variant2")], 1, function(x) paste(sort(x), collapse = "_"))
  results_total <- results_total[!duplicated(results_check), , drop = FALSE]

  if(verbose) print(paste0("Adjusting P values using ", p_value_adjustment, "."))
  results_total$P_adj <- stats::p.adjust(results_total$P, method = p_value_adjustment)
  rownames(results_total) <- NULL
  return(results_total)
}
