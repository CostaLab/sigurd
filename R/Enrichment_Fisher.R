#'Enrichment_FisherTest
#'@description
#'This function performs the Fisher test. The specified variant is 
#'@importFrom SummarizedExperiment assays colData
#'@importFrom stats fisher.test p.adjust
#'@param se SummarizedExperiment object.
#'@param variant The variant to be analysed. The values for this variant are retrieved from the consensus assay.
#'@param use_nocall Should the NoCall values be used? Then the non genotyped cells would be a separate group in addition to Alt and Ref. If FALSE, the NoCalls are removed. Default = FALSE
#'@param both_as_alt Should the consensus value for heterozyguous (3 or both) be treated as homozyguous (Alt)? Default = TRUE
#'@param trait The trait the variant should be associated to. This can be any column from the colData of the object.
#'@param trait_levels Should only trait levels be used for the comparison? For example, if Sample 1-3 are present but only 1 and 2 should be used in the comparison use trait_levels = c("Sample_1", "Sample_2"). NULL uses all levels. Default = NULL
#'@export
Enrichment_FisherTest <- function(se = NULL, variant = NULL, use_nocall = FALSE, both_as_alt = TRUE, trait = NULL, trait_levels = NULL){
  # Testing if all variables have been set.
  if(is.null(se)) stop("You did not set an input object.")
  if(is.null(variant)) stop("You did not set a variant.")
  if(is.null(trait)) stop("You did not set a trait.")
  if(!is.logical(use_nocall)) stop("use_nocall has to be TRUE or FALSE.")

  # Testing if the consensus assay is in the SE object.
  if(!"consensus" %in% names(SummarizedExperiment::assays(se))) stop("The consensus assay is not in the SE object.")

  # Testing if the traits are in the SE.
  if(!variant %in% rownames(se)) stop(paste0(variant, "is not a variant in the SE."))
  if(!trait %in% colnames(SummarizedExperiment::colData(se))) stop(paste0(trait, "is not a column in the SE colData."))

  print(paste0("We get a count matrix of the ", variant, " per ", trait, "."))
  # It can happen that there simply are no observations for a specific level in the data.
  # Then we can supply the levels by and and construct the count matrix by hand.
  if(both_as_alt){
    consensus_values <- c("0" = "NoCall", "1" = "Ref", "2" = "Alt", "3" = "Alt")
  } else {
    consensus_values <- c("0" = "NoCall", "1" = "Ref", "2" = "Alt", "3" = "Both")
  }
  variant_values <- as.character(SummarizedExperiment::assays(se)[["consensus"]][variant,])
  variant_values <- consensus_values[variant_values]
  names(variant_values) <- colnames(se)
  InputDataFrame <- data.frame(Variant = variant_values, Trait = SummarizedExperiment::colData(se)[,trait])
  if(!use_nocall) InputDataFrame <- subset(InputDataFrame, Variant != "NoCall")

  # We check if the trait levels are all in the data.
  if(!is.null(trait_levels) & !all(trait_levels %in% InputDataFrame[,"Trait"])) stop("Not all your trait levels are in the data.")
  if(!is.null(trait_levels)){
    variant_levels <- unique(InputDataFrame[, "Variant"])
    # We generate a matrix with the levels by hand.
    count_matrix <- matrix(0, nrow = length(variant_levels), ncol = length(trait_levels))
    rownames(count_matrix) <- variant_levels
    colnames(count_matrix) <- trait_levels
    for(i in 1:length(variant_levels)){
      row_use <- variant_levels[i]
      for(j in 1:length(trait_levels)){
        col_use <- trait_levels[j]
        count_matrix[i,j] <- sum(InputDataFrame[, "Variant"] == row_use & InputDataFrame[, "Trait"] == col_use)
      }
    }
  } else{
    count_matrix <- as.matrix(table(InputDataFrame[, "Variant"], InputDataFrame[, "Trait"]))
  }

  # We check if the count matrix has more than 1 column and row.
  number_rows <- nrow(count_matrix)
  row_check <- number_rows > 1
  if(!row_check) stop(paste0("You only have ", number_rows, " row (variant). You need at least two (two levels)."))
  number_cols <- ncol(count_matrix)
  col_check <- number_cols > 1
  if(!col_check) stop(paste0("You only have ", number_cols, " column (trait). You need at least two (two levels)."))

  results_matrix_p <- matrix(NA, ncol = nrow(count_matrix), nrow = ncol(count_matrix))
  colnames(results_matrix_p) <- rownames(count_matrix)
  rownames(results_matrix_p) <- colnames(count_matrix)
  results_matrix_odds <- results_matrix_p

  for(i in 1:ncol(count_matrix)){
    for(j in 1:nrow(count_matrix)){
      cells <- matrix(0, ncol = 2, nrow = 2)
      rownames(cells) <- c(rownames(count_matrix)[j], "Rest")
      colnames(cells) <- c(colnames(results_matrix_p)[i], "Rest")
      cells[1,1] <- count_matrix[j,i]
      cells[2,1] <- sum(count_matrix[rownames(count_matrix)[-j], i])
      cells[1,2] <- sum(count_matrix[j,colnames(count_matrix)[-i]])
      cells[2,2] <- sum(count_matrix[rownames(count_matrix)[-j], colnames(count_matrix)[-i]])
      ft <- stats::fisher.test(cells)
      results_matrix_p[i,j] <- ft$p.value
      results_matrix_odds[i,j] <- ft$estimate
    }
  }

  # We check if one of the Traits has only two conditions.
  # Then we only need to perform one comparison. This is important for the P value correction.
  # Otherwise, we would over correct the P values.
  variant_conditions <- nrow(count_matrix)
  trait_conditions <- ncol(count_matrix)
  if(all(c(variant_conditions, trait_conditions) == 2)){
    # If both Traits have only two conditions, we do not need to correct the P values.
    # We have only conducted one test.
    results_matrix_p_adj <- results_matrix_p
  } else if(variant_conditions == 2){
    # If we only have 2 variant_conditions, the columns of results_matrix_p are exactly the same.
    # The values are not exactly the same, but this is due to the float point gods.
    # But the test are simply from the opposite direction.
    # So, we only have to correct the P values per column.
    results_matrix_p_adj <- apply(results_matrix_p, 2, p.adjust, method = "fdr")
  } else if(trait_conditions == 2){
    # The same as above, but now we correct along the rows.
    results_matrix_p_adj <- t(apply(results_matrix_p, 1, p.adjust, method = "fdr"))
  } else{
    # If we have more than 2 conditions in both Traits, we have to correct all values.
    results_matrix_p_adj <- matrix(p.adjust(results_matrix_p, method = "fdr"), nrow = nrow(results_matrix_p), ncol = ncol(results_matrix_p))
  }
  colnames(results_matrix_p_adj) <- colnames(results_matrix_p)
  rownames(results_matrix_p_adj) <- rownames(results_matrix_p)

  # We want to have a heat map with the trait_to_compare2 on the y axis and the trait_to_compare1 on the x axis.
  # Then we can easily indicate, which traits are significantly different in which trait.
  # We generate a matrix with the -log10(P) values for a heatmap.
  heatmap_matrix <- matrix(NA, nrow = nrow(results_matrix_p_adj), ncol = ncol(results_matrix_p_adj))
  rownames(heatmap_matrix) <- rownames(results_matrix_p_adj)
  colnames(heatmap_matrix) <- colnames(results_matrix_p_adj)
  for(i in 1:nrow(heatmap_matrix)){
    for(j in 1:ncol(heatmap_matrix)){
      logged_pvalue <- -log10(results_matrix_p_adj[i,j])
      if(logged_pvalue <= -log10(0.05)) logged_pvalue <- 0
      logged_pvalue <- logged_pvalue * ifelse(results_matrix_odds[i,j] > 1, 1, 0)
      logged_pvalue <- ifelse(logged_pvalue < -log10(0.05), NA, logged_pvalue)
      heatmap_matrix[i,j] <- logged_pvalue
    }
  }

  counts <- data.frame(Variant = variant, as.data.frame(count_matrix))
  results <- list(P_Val = as.data.frame(results_matrix_p), P_Val_adj = as.data.frame(results_matrix_p_adj), Heatmap_Values = as.data.frame(heatmap_matrix), Counts = counts)
  return(results)
}
