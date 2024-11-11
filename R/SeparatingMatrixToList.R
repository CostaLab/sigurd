#'SeparatingMatrixToList
#'@description
#'We separate a matrix of variant information to a list.
#'Each variant is an entry in the list.
#'NoCalls (cells with no reads covering a variant) can be removed.
#'This function gets called by RowWiseSplit in return.
#'@importFrom stats na.omit
#'@param row_use The row the separate.
#'@param total_matrix The matrix to be split.
#'@param assay_to_split Which assay are we splitting?
#'@param remove_nocalls Do you want to remove NoCall cells?
#'@export
SeparatingMatrixToList <- function(row_use, total_matrix, consensus, assay_to_split = "consensus", remove_nocalls = TRUE){
  selected_row <- total_matrix[row_use,]
  selected_row <- stats::na.omit(selected_row)

  if(remove_nocalls == TRUE){
    # We remove the NoCall cells.
    selected_consensus <- consensus[row_use,]
    selected_consensus <- stats::na.omit(selected_consensus)
    selected_row <- selected_row[selected_consensus != 0]
    if(assay_to_split == "consensus"){
      selected_row[selected_row == 1] <- 0
      selected_row[selected_row >= 2] <- 1
    }
  } else if(remove_nocalls == FALSE){
    if(assay_to_split == "consensus"){
      selected_row[selected_row <= 1] <- 0
      selected_row[selected_row >= 2] <- 1
    }
  }
  return(selected_row)
}
