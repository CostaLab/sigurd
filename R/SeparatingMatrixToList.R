#'We separate a matrix of variant information to a list.
#'Each variant is an entry in the list.
#'@param row_to_save total_matrix
SeparatingMatrixToList <- function(row_to_save, total_matrix){
  selected_row <- total_matrix[row_to_save,]
  selected_row <- na.omit(selected_row)
  selected_row <- selected_row[selected_row != 0]
  selected_row[selected_row == 1] <- 0
  selected_row[selected_row >= 2] <- 1
  return(selected_row)
}
