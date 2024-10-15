#'Determining the clonal diversity.
#'
#'@description
#'This function determines the clonal diversity. The clones are defined using the function ClonalDefinition. It is not required to use clones, but any column in the column data of the object can be used.
#'
#'@importFrom SummarizedExperiment assays colData
#'@importFrom vegan diversity
#'@param se SummarizedExperiment object.
#'@param cells Which cells should be used? NULL uses all cells.
#'@param grouping The meta data	column used.
#'@param diversity_measure What diversity measure should be calculated? Can be EffectiveSpecies, shannon, simpson, invsimpson
#'@param base The base for the diversity 
#'@param verbose Should the function be verbose? Default = TRUE
#'@export
ClonalDiversity <- function(se, grouping = "Clones", cells = NULL, diversity_measure = "EffectiveSpecies", base = exp(1), verbose = TRUE){
  # We check if the requested diversity index is supported.
  if(!diversity_measure %in% c("EffectiveSpecies", "shannon", "simpson", "invsimpson")){
    stop(paste0("Your diversity measure variable ", diversity_measure, "is not supported."))
  }

  # We check if all the cells to be used are in the SE object.
  if(!is.null(cells)){
    cell_check <- all(cells %in% rownames(se))
    if(!cells){
      stop(paste0("Not all your cells are present in the SE object."))
    } else{
      if(verbose) print("We subset the SE object to only include the requested cells.")
      se <- se[, cells]
    }
  }

  # We check if the grouping variable is present in the SE object.
  if(!is.null(grouping)){
    if(!grouping %in% colnames(SummarizedExperiment::colData(se))){
      stop(paste0("Your grouping variable ", grouping, "is not in the SE object."))
    }
  }

  # We get the diversity for the grouping variable.
  clones <- SummarizedExperiment::colData(se)[, grouping]
  clones <- data.frame(table(clones))
  clones <- clones[,"Freq"]

  # We get the diversity.
  if(diversity_measure == "EffectiveSpecies"){
    diversity <- vegan::diversity(clones, index = "shannon", equalize.groups = FALSE, base = base)
    diversity <- base ^ diversity
  }
  if(diversity_measure == "shannon"){
    diversity <- vegan::diversity(clones, index = "shannon", equalize.groups = FALSE, base = base)
  }
  if(diversity_measure == "simpson"){
    diversity <- vegan::diversity(clones, index = "simpson", equalize.groups = FALSE)
  }
  if(diversity_measure == "invsimpson"){
    diversity <- vegan::diversity(clones, index = "invsimpson", equalize.groups = FALSE)
  }
  return(diversity)
}
