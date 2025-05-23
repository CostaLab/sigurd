#'Defining a single clone using multiple variants.
#'
#'@description
#'This function defines a single clone using a set of variants. A cell is considered to be part of the clone, if it has all of the provided variants. A maximum number of missing variants can be set.
#'
#'@importFrom SummarizedExperiment assays
#'@importFrom S4Vectors metadata
#'@param se SummarizedExperiment object.
#'@param variants Set of variants for clonal definition.
#'@param dropout The number of missing variants allowed.
#'@param variants_essential Which variants are considered essential? These will not be dropped if dropout is not NULL. The values can either be the variants names as character or indices for variants. A list of the same length as variants.
#'@param variants_excluding Which variants are not allowed in the clone? Variants names.  A list of the same length as variants.
#'@param min_vaf Minimum VAF for a cell to be considered positive.
#'@param clone_name What is the name of the clone? Default: SingularClone
#'@param verbose Should the function be verbose? Default = TRUE
#'@export
ClonalDefinition_Singular <- function(se, variants, dropout = FALSE, variants_essential = NULL, variants_excluding = NULL, min_vaf = 0.01, clone_name = "SingularClone", verbose = TRUE){
  # Ensuring variants are characters.
  variants <- lapply(variants, as.character)

  if(length(variants) > 1){
    # Checking that the sets of variants are not intersecting.
    variant_set_combinations <- combn(seq_along(variants), 2, simplify = FALSE)
    intersection_check <- any(
      vapply(variant_set_combinations, function(x) {
        length(intersect(variants[[x[1]]], variants[[x[2]]])) > 0
      }, logical(1))
    )
    if(intersection_check){
      stop(paste0("Your variants are intersecting. Please check the provided list of variants."))
    }
  }

  # Checking if all variants are in the object.
  variants_missing <- unlist(variants)
  variants_missing <- variants_missing[!variants_missing %in% rownames(se)]
  if(length(variants_missing) > 0){
    stop(paste0("The following variants are not in the object: ", paste0(variants_missing, collapse = ", ")))
  }
  
  # We generate a clone per set of variants. 
  new_clones <- setNames(rep("", ncol(se)), colnames(se))

  for(i in seq_along(variants)){
    variants_use <- variants[[i]]
    # Checking if the essential variants are characters or indices.
    variants_essential_usage <- FALSE
    if(!is.null(variants_essential)){
      variants_essential_usage <- TRUE
      # Checking if the length of variants_essential is the same as variants.
      if(length(variants_essential) != length(variants)){
        stop("The length of variants_essential is not equal to the length of variants.")
      }
      variants_essential_use <- variants_essential[[i]]
      variants_essential_character_check <- all(is.character(variants_essential_use))
      variants_essential_index_check <- all(is.numeric(variants_essential_use))
      if(variants_essential_character_check & variants_essential_index_check){
        stop("Your essential variants have to be either all characters (names) or numeric (indices). Please check you essential variants.")
      }
      # Checking if all essential variants are present.
      # This only happens when the essential variants are characters.
      if(variants_essential_character_check){
        variants_essential_missing <- variants_essential_use[!variants_essential_use %in% rownames(se)]
        if(length(variants_essential_missing) > 0){
          stop(paste0("The following essential variants are not in the object: ", paste0(variants_essential_missing, collapse = ", ")))
        }
        variants_essential_missing <- variants_essential_use[!variants_essential_use %in% variants_use]
        if(length(variants_essential_missing) > 0){
          stop(paste0("The following essential variants are not in the set of variants: ", paste0(variants_essential_missing, collapse = ", ")))
        }
      }
      # Replacing the indices with the names.
      if(variants_essential_index_check){
        variants_essential_use <- variants_use[variants_essential_use]
      }
    }

    # Checking for variants_excluding.
    variants_excluding_usage <- FALSE
    if(!is.null(variants_excluding)){
      variants_excluding_usage <- TRUE

      # Checking if the length of variants_essential is the same as variants.
      if(length(variants_excluding) != length(variants)){
        stop("The length of variants_excluding is not equal to the length of variants.")
      }

      variants_excluding_use <- variants_excluding[[i]]
      if(any(!is.character(variants_excluding_use))) stop("Not all your excluding variants are characters.")
      variants_excluding_in_variants <- variants_excluding_use[variants_excluding_use %in% variants_use]
      if(length(variants_excluding_in_variants) > 0) stop(paste0("The following variants are in variants and variants_excluding: ", paste0(variants_excluding_in_variants, collpase = ", ")))
      variants_excluding_presence <- variants_excluding_use[!variants_excluding_use %in% rownames(se)]
      if(length(variants_excluding_presence) > 0){
        if(verbose){
          print(paste0("The following excluding variants are not in the object: ", paste0(variants_excluding, collapse = ", ")))
          print("These are removed and the analysis proceeds.")
          variants_excluding_use <- variants_excluding_presence
        }
      }
    }

    # We check if a cell is mutated for the set of variants.
    check_use <- as.matrix(SummarizedExperiment::assays(se)[["consensus"]][variants_use, , drop = FALSE])	
    if(min_vaf > 0){
      frac_check <- as.matrix(SummarizedExperiment::assays(se)[["fraction"]][variants_use, , drop = FALSE])
      check_use[frac_check < min_vaf] <- 1
    }
    check_use <- matrix(check_use %in% 2:3, ncol = ncol(check_use), nrow = nrow(check_use), dimnames = list(rownames(check_use), colnames(check_use)))

    if(verbose) print("We checking how many variants a cell has.")
    check_all <- colSums(check_use)
  
    if(variants_essential_usage){
      if(verbose) print("Checking if the cells have the essential variants.")
      check_essential <- check_use[variants_essential_use, , drop = FALSE]
      check_essential <- colSums(check_essential) != length(variants_essential_use)
      check_all[check_essential] <- 0
    }
  
    if(dropout > 0){
      minimum_number_of_variants <- length(variants_use) - dropout
      if(verbose) print(paste0("Checking if a variant has at ", minimum_number_of_variants, " of ", length(variants_use), " variants."))
      check_all[check_all < minimum_number_of_variants] <- 0
    }
  
    if(variants_excluding_usage){
      if(verbose) print("Checking if a cell has an excluding variant.")
      check_excluding <- as.matrix(SummarizedExperiment::assays(se)[["consensus"]][variants_excluding_use, , drop = FALSE])	
      if(min_vaf > 0){
        frac_check <- as.matrix(SummarizedExperiment::assays(se)[["fraction"]][variants_excluding_use, , drop = FALSE])
        check_excluding[frac_check < min_vaf] <- 1
      }
      check_excluding <- matrix(check_excluding %in% 2:3, ncol = ncol(check_excluding), nrow = nrow(check_excluding), dimnames = list(rownames(check_excluding), colnames(check_excluding)))
      check_excluding <- colSums(check_excluding) > 0
      check_all[check_excluding] <- 0
    }

    new_meta_data <- setNames(ifelse(check_all > 0, paste0("C", i), ""), names(check_all))
    new_clones[names(new_meta_data)] <- paste0(new_clones[names(new_meta_data)], new_meta_data)
  }
  new_clones[new_clones == ""] <- "Other"
  SummarizedExperiment::colData(se)[, clone_name] <- new_clones
  singular_clone_meta_data <- list(variants = variants, dropout = dropout, variants_essential = variants_essential, variants_excluding = variants_excluding, min_vaf = min_vaf)
  S4Vectors::metadata(se)[[clone_name]] <- singular_clone_meta_data
  # We return the SE object as a result.
  return(se)
}
