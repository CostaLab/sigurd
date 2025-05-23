#'Filtering the loaded genotyping data.
#'
#'@description
#'We filter a SummarizedExperiment object to exclude variants and cells.
#'
#'@details
#'We do this for one sample at a time.
#'We want to remove:
#' \enumerate{
#'   \item all cells not in an allow list,
#'   \item all cells in an exclusion list,
#'   \item all cells that do not have at least one variant with >1 (Reference),
#'   \item all variants that are for alternative transcripts,
#'   \item all variants that are always NoCall,
#'   \item set variants with a VAF below a threshold to NoCall or Reference.
#' }
#'@importFrom Matrix summary
#'@importFrom SummarizedExperiment assays
#'@importFrom utils read.table
#'@importFrom Matrix rowSums colSums
#'@param SE SummarizedExperiment object.
#'@param cells_include A vector of cell barcodes. Only these cells will be retained. 
#'@param cells_exclude A vector of cell barcodes. These cells will be removed from the output.
#'@param fraction_threshold Variants with an VAF below this threshold are set to 0. Numeric. Default = NULL.
#'@param keep_refs If a cell is below the fraction_threshold, should the number of reference reads be kept? If FALSE, the reference reads will be removed. Default = TRUE
#'@param alts_threshold Variants with a number of alt reads less than this threshold are set to 0. Numeric. Default = NULL.
#'@param min_cells_per_variant In how many cells should a variant be present to be included? Numeric. Default = 2.
#'@param min_variants_per_cell How many variants should be covered in a cell have to be included? Default = 1.
#'@param reject_value Should cells that fall below a threshold (fraction_threshold or alts_threshold) be treated as Reference or NoCall? Default = NoCall.
#'@param verbose Should the function be verbose? Default = TRUE
#'@examples
#' \dontrun{
#'   # Removing all variants that are not detected in at least 2 cells.
#'   # Before we remove the variants, we set fraction value of variants below 0.05 to 0.
#'   se_geno <- Filtering(se_geno, min_cells_per_variant = 2, fraction_threshold = 0.05)
#' }
#'@export
Filtering <- function(SE, cells_include = NULL, cells_exclude = NULL, fraction_threshold = NULL, keep_refs = TRUE, alts_threshold = NULL, min_cells_per_variant = NULL, min_variants_per_cell = NULL, reject_value = "NoCall", verbose = TRUE){
  # Checking if the reject_value variable is correct.
  if(!reject_value %in% c("Reference", "NoCall")){
    stop(paste0("Your reject_value is ", reject_value, ".\nIt should be Reference or NoCall."))
  } else{
    if(reject_value == "NoCall") reject_value_numeric <- 0
    if(reject_value == "Reference") reject_value_numeric <- 1
  }

  # Checking if the keep_refs value is boolean.
  if(!is.logical(keep_refs)){
    stop("Your value for keep_refs is not either TRUE or FALSE.")
  }

  if(!is.null(cells_include)){
    if(verbose) print("We remove all cells not in the allow list.")
    # Check if any cells would be left.
    check_leftover <- any(colnames(SE) %in% cells_include)
    if(!check_leftover) stop("No cells are in your supplied list.")
    SE <- SE[, colnames(SE) %in% cells_include]
    cells_include <- cells_include[cells_include %in% colnames(SE)]
    SE <- SE[, cells_include]
  }


  if(!is.null(cells_exclude)){
    if(verbose) print("We remove the unwanted cell barcodes.")
    # Check if any cells would be left.
    check_leftover <- all(colnames(SE) %in% cells_exclude)
    if(check_leftover) stop("No cells are in your supplied list.")
    SE <- SE[, !colnames(SE) %in% cells_exclude]
  }


  # If the fraction_threshold is 0, we skip the thresholding. 0 and NULL would be the same.
  # We might want to use a fraction_threshold of 0 to use the same variable to create file paths later.
  if(!is.null(fraction_threshold)){
    if(fraction_threshold == 0){
      fraction_threshold <- NULL
    }
  }


  if(!is.null(fraction_threshold)){
    if(any(fraction_threshold >= 1, fraction_threshold <= 0)){
      stop("Your fraction threshold is not 0 < x < 1.")
    }
    if(verbose) print(paste0("We assume that cells with a fraction smaller than our threshold are actually ", reject_value, "."))
    if(verbose) print(paste0("We set consensus values to ", reject_value_numeric, " (", reject_value, ") and fraction values to 0."))
    if(verbose) print(paste0("We do not set fractions between ", fraction_threshold, " and 1 to 1."))
    if(verbose) print("This way, we retain the heterozygous information.")
    # Filtering using sparse matrices.
    consensus_matrix <- SummarizedExperiment::assays(SE)$consensus
    fraction_matrix <- SummarizedExperiment::assays(SE)$fraction
    coverage_matrix <- SummarizedExperiment::assays(SE)$coverage
    alts_matrix <- SummarizedExperiment::assays(SE)$alts
    position_matrix <- Matrix::summary(fraction_matrix)
    position_matrix <- subset(position_matrix, x > 0 & x < fraction_threshold)
    # If no elements fall between 0 and the fraction_threshold, we do not have to change the matrices.
    if(nrow(position_matrix) > 0){
      ij <- as.matrix(position_matrix[, 1:2])
      consensus_matrix[ij] <- reject_value_numeric
      fraction_matrix[ij] <- 0
      coverage_matrix[ij] <- coverage_matrix[ij] - alts_matrix[ij]
      alts_matrix[ij] <- 0
      if(!keep_refs){
        refs_matrix <- SummarizedExperiment::assays(SE)$refs
        coverage_matrix[ij] <- coverage_matrix[ij] - refs_matrix[ij]
        refs_matrix[ij] <- 0
        SummarizedExperiment::assays(SE)$refs <- refs_matrix
      }
      SummarizedExperiment::assays(SE)$consensus <- consensus_matrix
      SummarizedExperiment::assays(SE)$fraction <- fraction_matrix
      SummarizedExperiment::assays(SE)$coverage <- coverage_matrix
      SummarizedExperiment::assays(SE)$alts <- alts_matrix
    }
  }


  if(!is.null(alts_threshold)){
    if(any(!is.numeric(alts_threshold), is.infinite(alts_threshold))){
      stop("Your alts_threshold should be a numeric value.")
    }
    if(verbose) print(paste0("We assume that cells with a number of alternative reads smaller than our threshold are actually ", reject_value, "."))
    if(verbose) print(paste0("We set consensus values to ", reject_value_numeric, " (", reject_value, "), fraction values to 0."))
    if(reject_value == "NoCall") if(verbose) print("We set Alts, Refs and Coverage to 0.")
    if(reject_value == "Reference") if(verbose) print("We set Alts to 0 and adjust the Coverage.")
    # Filtering using sparse matrices.
    consensus_matrix <- SummarizedExperiment::assays(SE)$consensus
    fraction_matrix <- SummarizedExperiment::assays(SE)$fraction
    coverage_matrix <- SummarizedExperiment::assays(SE)$coverage
    alts_matrix <- SummarizedExperiment::assays(SE)$alts
    refs_matrix <- SummarizedExperiment::assays(SE)$refs
    position_matrix <- summary(alts_matrix)
    position_matrix <- subset(position_matrix, x < alts_threshold)
    # If no elements fall between 0 and the alts_threshold, we do not have to change the matrices.
    if(nrow(position_matrix) > 0){
      ij <- as.matrix(position_matrix[, 1:2])
      consensus_matrix[ij] <- reject_value_numeric
      fraction_matrix[ij] <- 0
      # We check if the positions to be changed have ref reads. If not, we set them to NoCall.
      refs_check <- refs_matrix[ij] == 0
      consensus_matrix[ij][refs_check] <- 0
      # Since we remove the Alt reads in either case (NoCall or Reference), we remove the Alts from the Coverage.
      coverage_matrix[ij] <- coverage_matrix[ij] - alts_matrix[ij]
      # We now set the alts matrix to 0.
      alts_matrix[ij] <- 0
      # If we set the cells to NoCall, we also set the Refs to 0.
      # We then also set the coverage matrix to 0.
      if(reject_value == "NoCall"){
        refs_matrix[ij] <- 0
        coverage_matrix[ij] <- 0
      }
      SummarizedExperiment::assays(SE)$consensus <- consensus_matrix
      SummarizedExperiment::assays(SE)$fraction <- fraction_matrix
      SummarizedExperiment::assays(SE)$coverage <- coverage_matrix
      SummarizedExperiment::assays(SE)$alts <- alts_matrix
      SummarizedExperiment::assays(SE)$refs <- refs_matrix
    }
  }


  if(verbose) print("We remove all the variants that are always NoCall.")
  consensus_test <- SummarizedExperiment::assays(SE)$consensus > 0
  keep_variants <- Matrix::rowSums(consensus_test) > 0
  SE <- SE[keep_variants,]

  if(!is.null(min_cells_per_variant)){
    if(verbose) print(paste0("We remove variants, that are not at least detected in ", min_cells_per_variant, " cells."))
    keep_variants <- Matrix::rowSums(SummarizedExperiment::assays(SE)$consensus >= 2)
    keep_variants <- keep_variants >= min_cells_per_variant
    SE <- SE[keep_variants,]
  }

  if(!is.null(min_variants_per_cell)){
    if(verbose) print(paste0("We remove all cells that are not >= 1 (Ref) for at least ", min_variants_per_cell, " variant."))
    consensus_test <- SummarizedExperiment::assays(SE)$consensus >= 1
    keep_cells <- Matrix::colSums(consensus_test) >= min_variants_per_cell
    SE <- SE[,keep_cells]
  }
  return(SE)
}
