#'We separate a matrix of variant information to a list.
#'Each variant is an entry in the list.
#'@param row_use The row the separate.
#'@param total_matrix The matrix to be split.
#'@export
SeparatingMatrixToList <- function(row_use, total_matrix){
  selected_row <- total_matrix[row_use,]
  selected_row <- na.omit(selected_row)
  selected_row <- selected_row[selected_row != 0]
  selected_row[selected_row == 1] <- 0
  selected_row[selected_row >= 2] <- 1
  return(selected_row)
}
