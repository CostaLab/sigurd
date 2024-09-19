#'Defining clones.
#'
#'@description
#'This function defines clones using a combination of variants supplied. Each set of variants can be and should be paired with a identity, like sample or cell type. The clones are then only determined in respect to this group. If this information is not provided, the clones are determined for all cells present.
#'
#'@importFrom Matrix summary
#'@importFrom SummarizedExperiment assays
#'@importFrom utils combn read.table
#'@importFrom Matrix rowSums colSums
#'@importFrom S4Vectors metadata
#'@param se SummarizedExperiment object.
#'@param variants_ls List of variants for clonal definition
#'@param grouping The meta data	column used to split the cells into groups. Default = NULL
#'@param identities Vector of groups, like samples.
#'@param verbose Should the function be verbose? Default = TRUE
#'@export
ClonalDefinition <- function(se, variants_ls, grouping = NULL, identities = NULL, verbose = TRUE){
  # Checking if the group variable is in the column data.
  if(!is.null(grouping)){
    # We check if the grouping variable is in the column data.
    if(!grouping %in% colnames(SummarizedExperiment::colData(se))){
      stop(paste0("Your grouping variable (", grouping, ") is not in the column data."))
    }
    # We check if the identities are not NULL.
    if(is.null(identities)){
      stop("No identities supplied.")
    }
    # We check if all the identities of interest are in the column data.
    check_identities <- all(identities %in% SummarizedExperiment::colData(se)[,grouping])
    if(!check_identities){
      stop("Not all identities are in the column data. Please check for missing groups.")
    }
  }
  # Checking if the supplied variants are in a list.
  if(!is.list(variants_ls)){
    stop("You have to supply the variants as a list.")
  }
  # Checking if the length of the variants list is equal to the length of the identities.
  if(!is.null(identities)){
    if(length(identities) != length(variants_ls)){
      stop(paste0("You have ", length(identities), " identities and ", length(variants_ls), " sets of variants. Please check."))
    }
  }
  # Checking if any of the variants is not in the SE object.
  for(i in 1:length(variants_ls)){
    check_variants <- all(variants_ls[[i]] %in% rownames(se))
    if(!check_variants){
      stop(paste0("The set of variants number ", i, " has variants not in the object. Please check."))
    }
  }

  # For each set of variants, we get all possible combinations.
  combinations_ls <- list()
  for(i in 1:length(variants_ls)){
    combinations_ls[i] <- list(lapply(1:length(variants_ls[[i]]), function(x) combn(variants_ls[[i]], x)))
  }

  # We prepare a new meta data column.
  new_meta_data <- rep(NA, ncol(se))
  names(new_meta_data) <- colnames(se)

  # We get the lineages.
  combinations_meta_data <- list()
  for(i in 1:length(combinations_ls)){
    combinations <- combinations_ls[[i]]
    # We get the variants used to create the combinations.
    # We need them for the negative lineage.
    variants_uncombined <- variants_ls[[i]]

    # If the grouping variables is not NULL, we get the relevant identity.
    if(!is.null(grouping)){
      identity_use <- identities[i]
      cells_use <- SummarizedExperiment::colData(se)
      cells_use <- cells_use[cells_use[,grouping] == identity_use, ]
      cells_use <- rownames(cells_use)
      se_use <- se[,cells_use]
    } else{
      se_use <- se
    }

    # We prepare the new meta data for the subset of the cells.
    new_meta_data_subset <- rep(NA, ncol(se_use))
    names(new_meta_data_subset) <- colnames(se_use)

    # We get the negative clone.
    cells <- as.matrix(SummarizedExperiment::assays(se_use)[["consensus"]][variants_uncombined, , drop = FALSE])
    cells <- matrix(cells %in% 2:3, ncol = ncol(cells), nrow = nrow(cells), dimnames = list(rownames(cells), colnames(cells)))
    cells <- colSums(cells) == 0
    cells <- cells[cells]
    new_meta_data_subset[names(cells)] <- "Negative"
    combinations_meta_data[[i]] <- data.frame(Clone = "Negative", Variants = "Negative")

    combination_number <- 0
    for(j in 1:length(combinations)){
      for(k in 1:ncol(combinations[[j]])){
        combination_number <- combination_number + 1
        if(verbose) print(paste0("Combination: ", combination_number))
        combination_use <- combinations[[j]][,k]
        cells <- as.matrix(SummarizedExperiment::assays(se_use)[["consensus"]][combination_use, , drop = FALSE])
        cells <- matrix(cells %in% 2:3, ncol = ncol(cells), nrow = nrow(cells), dimnames = list(rownames(cells), colnames(cells)))
        cells <- colSums(cells) == length(combination_use)
        # The cells should only be positive for the current combination.
        variants_rest <- variants_uncombined[!variants_uncombined %in% combination_use]
        if(length(variants_rest) > 0){
          cells_rest <- as.matrix(SummarizedExperiment::assays(se_use)[["consensus"]][variants_rest, , drop = FALSE])
          cells_rest <- matrix(cells_rest %in% 2:3, ncol = ncol(cells_rest), nrow = nrow(cells_rest), dimnames = list(rownames(cells_rest), colnames(cells_rest)))
          cells_rest <- colSums(cells_rest) == 0
        } else{
          cells_rest <- rep(TRUE, ncol(se_use))
        }
        cells <- colSums(rbind(cells, cells_rest)) == 2
        coverage <- mean(colMeans(SummarizedExperiment::assays(se_use)[["coverage"]][combination_use, cells, drop = FALSE]))
        cells <- cells[cells]
        new_meta_data_subset[names(cells)] <- paste0("C_", combination_number)
        combinations_meta_data[[i]] <- rbind(combinations_meta_data[[i]], 
                                             data.frame(Clone = paste0("C_", combination_number), Variants = paste0(combinations[[j]][,k], collapse = ",")))
      }
    }
    rownames(combinations_meta_data[[i]]) <- combinations_meta_data[[i]][, "Clone"]

    # We rename the clones according to their size.
    variant_order <- sort(table(new_meta_data_subset), decreasing = TRUE)
    new_names <- data.frame(OldName = names(variant_order), NewName = NA, row.names = names(variant_order))
    new_rank <- 1
    for(j in 1:nrow(new_names)){
      if(new_names[j, "OldName"] == "Negative"){
        new_names[j, "NewName"] <- "Negative"
      } else{
        new_names[j, "NewName"] <- paste0("C_", new_rank)
        new_rank <- new_rank + 1
      }
    }
    new_names <- rbind(new_names, c("NoLineage", "NoLineage"))
    rownames(new_names)[nrow(new_names)] <- "NoLineage"
    new_meta_data_subset <- new_names[new_meta_data_subset, "NewName"]
    names(new_meta_data_subset) <- colnames(se_use)

    # We remove unused clones.
    combinations_meta_data[[i]] <- combinations_meta_data[[i]][rownames(combinations_meta_data[[i]]) %in% rownames(new_names),]
    combinations_meta_data[[i]][,"Clone"] <- new_names[rownames(combinations_meta_data[[i]]),"NewName"]
    rownames(combinations_meta_data[[i]]) <- combinations_meta_data[[i]][,"Clone"]

    # We add the new meta data to the combined meta data.
    new_meta_data[names(new_meta_data_subset)] <- new_meta_data_subset
  }
  # If we used identities, we set the names of the combinations_meta_data list to them.
  if(!is.null(identities)){
    names(combinations_meta_data) <- identities
  }

  # We add the new meta data to the column data.
  new_column_data <- SummarizedExperiment::colData(se)
  new_column_data$Clones <- new_meta_data
  SummarizedExperiment::colData(se) <- new_column_data
  S4Vectors::metadata(se)[["ClonalLineages"]] <- combinations_meta_data

  # We return the SE object as a result.
  return(se)
}
