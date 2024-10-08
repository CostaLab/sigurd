#'SeparatingMatrixToList
#'@description
#'We separate a matrix of variant information to a list.
#'Each variant is an entry in the list.
#'NoCalls (cells with no reads covering a variant) can be removed.
#'This function gets called by RowWiseSplit in return.
#'@importFrom stats na.omit
#'@param row_use The row the separate.
#'@param total_matrix The matrix to be split.
#'@param remove_nocalls Do you want to remove NoCall cells?
#'@export
SeparatingMatrixToList <- function(row_use, total_matrix, remove_nocalls = TRUE){
  selected_row <- total_matrix[row_use,]
  selected_row <- stats::na.omit(selected_row)

  if(remove_nocalls == TRUE){
    # We remove the NoCall cells.
    selected_row <- selected_row[selected_row != 0]
    selected_row[selected_row == 1] <- 0
    selected_row[selected_row >= 2] <- 1
  } else if(remove_nocalls == FALSE){
    selected_row[selected_row <= 1] <- 0
    selected_row[selected_row >= 2] <- 1
  }
  return(selected_row)
}
