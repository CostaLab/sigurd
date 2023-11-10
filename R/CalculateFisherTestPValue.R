#'CalculateFisherTestPValue
#'
#'@description
#'We perform the Fisher Test of SNVs and calculate the P values.
#'@importFrom stats fisher.test
#'@param variant_values The fraction values you are analysing. A vector.
#'@param other_mutation All other variants you have. A vector of variant names.
#'@param all_variants_list List of fraction values for all the variants you want to compare your variant with.
#'@param min_intersecting_cells Minimum number of intersecting cells. Correlations with less than this will not be performed.
#'@export
CalculateFisherTestPValue <- function(variant_values, other_mutation, all_variants_list, min_intersecting_cells = 5){
  other_variant_values <- all_variants_list[[other_mutation]]
  if(sum(names(variant_values) %in% names(other_variant_values)) == 0){
    #print("No intersection")
    return(c(NA,NA,NA,NA,NA,NA))
  } else if(sum(names(variant_values) %in% names(other_variant_values)) < min_intersecting_cells){
    #print("Less than min_intersecting_cells cells intersect.")
    return(c(NA,NA,NA,NA,NA,NA))
  } else{
    variant_values <- variant_values[names(variant_values) %in% names(other_variant_values)]
    other_variant_values <- other_variant_values[names(other_variant_values) %in% names(variant_values)]
    other_variant_values <- other_variant_values[names(variant_values)]
    # We need to have some mutations in both sets.
    # If this is not the case, we are returning NA.
    if(any(!any(variant_values == 1), !any(other_variant_values == 1))){
      #print("There aren't mutated cells in both sets.")
      result <- c(NA,NA,NA,NA,NA,NA)
      return(result)
    } else if(length(variant_values) > 2){
      # We get the following count matrix.
      # count_matrix[1,1] <- All cells that have Variant 1 AND Variant 2.
      # count_matrix[2,1] <- All cells that only have Variant 1 and not Variant 2.
      # count_matrix[1,2] <- All cells that have only Variant 2, but not Variant 1.
      # count_matrix[2,2] <- All cells that have neither Variant 1 or Variant 2.
      count_matrix <- matrix(0, ncol = 2, nrow = 2)
      colnames(count_matrix) <- c("Variant1_Alt", "Variant1_Ref")
      rownames(count_matrix) <- c("Variant2_Alt", "Variant2_Ref")
      count_matrix[1,1] <- sum(variant_values == 1 & other_variant_values == 1)
      count_matrix[2,1] <- sum(variant_values == 1 & other_variant_values == 0)
      count_matrix[1,2] <- sum(variant_values == 0 & other_variant_values == 1)
      count_matrix[2,2] <- sum(variant_values == 0 & other_variant_values == 0)
      result <- stats::fisher.test(x = count_matrix)
      result <- c(result$p.value, result$estimate, count_matrix[1,1], count_matrix[2,1], count_matrix[1,2], count_matrix[2,2])
    } else{
      #print("We do not have more than 2 cells for the somatic variant.")
      result <- c(NA,NA,NA,NA,NA,NA)
    }
  }
  return(result)
}
