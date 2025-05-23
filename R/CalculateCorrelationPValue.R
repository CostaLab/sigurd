#'Correlating the SNVs
#'
#'@description
#'We perform the correlation of SNVs and calculate the P values.
#'@importFrom stats cor.test
#'@param variant_values The fraction values you are analyzing. A vector.
#'@param other_mutation All other variants you have. A vector of variant names.
#'@param all_variants_list List of fraction values for all the variants you want to compare your variant with.
#'@param min_intersecting_cells Minimum number of intersecting cells. Correlations with less than this will not be performed.
#'@param value_type Are we using consensus or other information?
#'@export
CalculateCorrelationPValue <- function(variant_values, other_mutation, all_variants_list, value_type = "consensus", min_intersecting_cells = 5){
  other_variant_values <- all_variants_list[[other_mutation]]
  if(sum(names(variant_values) %in% names(other_variant_values)) == 0){
    #print("No intersection")
    return(c(NA,NA,NA,NA,NA,NA))
  } else if(sum(names(variant_values) %in% names(other_variant_values)) < min_intersecting_cells){
    #print("Less than 5 cells intersect.")
    return(c(NA,NA,NA,NA,NA,NA))
  } else{
    variant_values <- variant_values[names(variant_values) %in% names(other_variant_values)]
    other_variant_values <- other_variant_values[names(other_variant_values) %in% names(variant_values)]
    other_variant_values <- other_variant_values[names(variant_values)]
    # We need to have some mutations in both sets.
    # If this is not the case, we are returning NA.
    if(any(all(diff(variant_values) == 0), all(diff(other_variant_values) == 0))){
      #print("There aren't mutated cells in both sets.")
      result <- c(NA,NA,NA,NA,NA,NA)
      return(result)
    } else if(length(variant_values) > 2){
      result <- stats::cor.test(variant_values, other_variant_values)
      if(value_type == "consensus"){
        cells_som_alt <- sum(variant_values == 1)
        cells_som_ref <- sum(variant_values == 0)
        cells_MT_alt  <- sum(other_variant_values == 1)
        cells_MT_ref  <- sum(other_variant_values == 0)
      } else{
        cells_som_alt <- sum(variant_values > 0)
        cells_som_ref <- sum(variant_values == 0)
        cells_MT_alt  <- sum(other_variant_values > 0)
        cells_MT_ref  <- sum(other_variant_values == 0)
      }
      result <- c(result$p.value, result$estimate, cells_som_alt, cells_som_ref, cells_MT_alt, cells_MT_ref)
    } else{
      #print("We do not have more than 2 cells for the somatic variant.")
      result <- c(NA,NA,NA,NA,NA,NA)
    }
  }
  return(result)
}
