SeparatingMatrixToList <- function(row_to_save, total_matrix){
  selected_row <- total_matrix[row_to_save,]
  selected_row <- selected_row[selected_row != 0]
  selected_row[selected_row == 1] <- 0
  selected_row[selected_row >= 2] <- 1
  return(selected_row)
}
